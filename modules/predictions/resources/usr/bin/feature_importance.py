#!/usr/bin/env python
"""
Feature Importance Analysis by Relationship Type

Trains Random Forest classifiers on knowledge graph features for each
KG × embedding × indication combination and extracts Gini feature importances,
aggregated by relationship (entity) type.

Usage:
    feature_importance.py <scores_date>

Arguments:
    scores_date: Clinical scores date (YYYY-MM-DD format)

Output:
    - data/feature_importance.csv: Per-type aggregated importances
    - plots/feature_importance_by_type.png: Boxplot of importance by type across KGs/indications
    - plots/feature_importance_by_type_per_kg.png: Faceted by KG
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier

from alethiotx.artemis.cv.pipeline import prepare as prepare_model
from alethiotx.artemis.clinical.scores import load as load_clinical_scores
from alethiotx.artemis.clinical.scores import unique as unique_clinical_scores
from alethiotx.artemis.clinical.scores import all_targets


# ─── Configuration ───────────────────────────────────────────────────────────

INDICATIONS = ['breast', 'lung', 'bowel', 'prostate', 'melanoma', 'diabetes', 'cardiovascular']
KGS = ['hetionet', 'biokg', 'openbiolink', 'primekg']
EMBEDDINGS = ['RotatE']
N_ITERATIONS = 3  # Use 3 seeds for stability without excessive computation
RANDOM_SEEDS = [1, 2, 3]

KG_DATA_PREFIX = 's3://alethiotx-artemis/data/kgs-no-data-leakage/associations'


# ─── Helper Functions ────────────────────────────────────────────────────────

def load_feature_type_mapping(kg: str, embedding: str) -> dict:
    """
    Load the mapping from feature column name to entity/relationship type.

    For each KG, loads terms.csv (ordered entity IDs used as features) and
    terms_hash_table.csv (entity metadata) to build a positional mapping
    from feature columns in predictions.parquet to their entity types.

    Parameters
    ----------
    kg : str
        Knowledge graph name
    embedding : str
        Embedding method name

    Returns
    -------
    dict
        Mapping from feature column name (str) → entity type (str)
    """
    base = f'{KG_DATA_PREFIX}/{kg}/{embedding}/prepare'

    # Load terms.csv: ordered entity IDs (header is first entity, rows are rest)
    terms_df = pd.read_csv(f'{base}/terms.csv')
    entity_ids = [terms_df.columns[0]] + terms_df.iloc[:, 0].tolist()

    # Load terms_hash_table.csv
    ht = pd.read_csv(f'{base}/terms_hash_table.csv', header=None)

    if ht.shape[1] == 2:
        # Hetionet format: col0 = "Type::ID", col1 = human_readable_name
        # Entity type is prefix of col0
        id_to_type = {row[0]: row[0].split('::')[0] for _, row in ht.iterrows()}
    else:
        # BioKG / OpenBioLink / PrimeKG format: col0 = entity_id, col1 = name, col2 = type
        id_to_type = dict(zip(ht[0], ht[2]))

    # Load predictions.parquet - columns correspond to entity_ids in terms.csv order
    kg_path = f'{KG_DATA_PREFIX}/{kg}/{embedding}/summarize/predictions.parquet'
    kg_features = pd.read_parquet(kg_path)

    # Build positional type mapping: position i → entity_id i → type
    feature_types = []
    for eid in entity_ids:
        etype = id_to_type.get(eid, 'Unknown')
        feature_types.append(etype)

    col_to_type = dict(zip(kg_features.columns, feature_types))

    return col_to_type, kg_features


def compute_feature_importance(
    kg_features: pd.DataFrame,
    clinical_scores: pd.DataFrame,
    known_targets: list,
    rand_seed: int
) -> tuple:
    """
    Train RF model and return feature importances with column names.

    Parameters
    ----------
    kg_features : pd.DataFrame
        Knowledge graph feature matrix (genes × features)
    clinical_scores : pd.DataFrame
        Clinical scores for one indication
    known_targets : list
        All known targets across indications
    rand_seed : int
        Random seed for reproducibility

    Returns
    -------
    tuple
        (feature importance array, column names used for training)
    """
    training_data = prepare_model(
        kg_features,
        clinical_scores,
        known_targets=known_targets,
        rand_seed=rand_seed
    )

    clf = RandomForestClassifier(random_state=rand_seed)
    clf.fit(training_data['X'], training_data['y_binary'])

    return clf.feature_importances_, training_data['X'].columns


def aggregate_importances_by_type(
    importances: np.ndarray,
    col_to_type: dict,
    columns: pd.Index
) -> pd.DataFrame:
    """
    Aggregate feature importances by entity/relationship type.

    Computes both sum and mean (per-feature) importance for each type.
    Mean importance normalizes for the number of features per type,
    revealing which types carry more signal per feature than expected.

    Parameters
    ----------
    importances : np.ndarray
        Per-feature importance values
    col_to_type : dict
        Mapping from column name to entity type
    columns : pd.Index
        Feature column names (ordered)

    Returns
    -------
    pd.DataFrame
        Per-type importance with columns: sum, mean, count
    """
    imp_df = pd.DataFrame({
        'importance': importances,
        'type': [col_to_type.get(c, 'Unknown') for c in columns]
    })
    agg = imp_df.groupby('type')['importance'].agg(['sum', 'mean', 'count'])
    return agg


# ─── Main Execution ──────────────────────────────────────────────────────────

def main():
    """Execute feature importance analysis."""
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)

    scores_date = sys.argv[1]

    print(f"Feature Importance Analysis")
    print(f"  Scores date: {scores_date}")
    print(f"  KGs: {KGS}")
    print(f"  Embeddings: {EMBEDDINGS}")
    print(f"  Indications: {INDICATIONS}")
    print()

    # Load clinical data (using Unique filter as per optimal config)
    clinical_data_raw = load_clinical_scores(date=scores_date)
    clinical_data_raw = unique_clinical_scores(list(clinical_data_raw))
    known_targets = all_targets(list(clinical_data_raw))

    results = []

    for kg in KGS:
        for embedding in EMBEDDINGS:
            print(f"Processing {kg}/{embedding}...")

            # Load feature type mapping and KG features
            try:
                col_to_type, kg_features = load_feature_type_mapping(kg, embedding)
            except Exception as e:
                print(f"  Error loading {kg}/{embedding}: {e}")
                continue

            # Filter clinical data to KG features
            clinical_data = {}
            for indication, scores in zip(INDICATIONS, clinical_data_raw):
                filtered = scores[scores['Target Gene'].isin(kg_features.index)]
                clinical_data[indication] = filtered

            for indication in INDICATIONS:
                for seed in RANDOM_SEEDS:
                    print(f"  {indication} (seed={seed})...")

                    try:
                        importances, train_columns = compute_feature_importance(
                            kg_features,
                            clinical_data[indication],
                            known_targets,
                            seed
                        )

                        # Aggregate by type
                        type_importances = aggregate_importances_by_type(
                            importances, col_to_type, train_columns
                        )

                        for etype, row in type_importances.iterrows():
                            results.append({
                                'kg': kg,
                                'embedding': embedding,
                                'indication': indication,
                                'seed': seed,
                                'entity_type': etype,
                                'importance_sum': row['sum'],
                                'importance_mean': row['mean'],
                                'feature_count': int(row['count'])
                            })
                    except Exception as e:
                        print(f"    Error: {e}")
                        continue

    # Compile results
    results_df = pd.DataFrame(results)

    # Save raw data
    data_dir = Path('data')
    data_dir.mkdir(exist_ok=True)
    results_df.to_csv(data_dir / 'feature_importance.csv', index=False)
    print(f"\n✓ Saved data/feature_importance.csv ({len(results_df)} rows)")

    # ─── Visualization ───────────────────────────────────────────────────────

    plots_dir = Path('plots')
    plots_dir.mkdir(exist_ok=True)

    # Compute feature fraction (%) per type within each model run
    results_df['feature_fraction_pct'] = results_df.groupby(
        ['kg', 'embedding', 'indication', 'seed']
    )['feature_count'].transform(lambda x: x / x.sum() * 100)

    # Compute sum-based importance percentage
    results_df['importance_sum_pct'] = results_df.groupby(
        ['kg', 'embedding', 'indication', 'seed']
    )['importance_sum'].transform(lambda x: x / x.sum() * 100)

    # Compute normalized mean importance (relative to uniform baseline)
    results_df['mean_importance_normalized'] = results_df.groupby(
        ['kg', 'embedding', 'indication', 'seed']
    )['importance_mean'].transform(lambda x: x / x.mean())

    # Determine type order by normalized mean importance
    type_order = (
        results_df.groupby('entity_type')['mean_importance_normalized']
        .median()
        .sort_values(ascending=False)
        .index.tolist()
    )

    # ─── Plot 1: Feature fraction vs Sum importance (side-by-side bars) ──────

    # Aggregate medians for comparison
    comparison = results_df.groupby('entity_type').agg(
        feature_fraction=('feature_fraction_pct', 'median'),
        importance_sum=('importance_sum_pct', 'median')
    ).loc[type_order]

    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(type_order))
    width = 0.35
    bars1 = ax.bar(x - width/2, comparison['feature_fraction'], width,
                   label='Feature fraction (%)', color='#7fbfff', edgecolor='black', linewidth=0.5)
    bars2 = ax.bar(x + width/2, comparison['importance_sum'], width,
                   label='Gini importance (%)', color='#ff7f7f', edgecolor='black', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(type_order, rotation=45, ha='right')
    ax.set_xlabel('Relationship Type')
    ax.set_ylabel('Percentage (%)')
    ax.set_title('Feature Fraction vs Gini Importance by Relationship Type')
    ax.legend()
    plt.tight_layout()
    fig.savefig(plots_dir / 'feature_fraction_vs_importance.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("✓ Saved plots/feature_fraction_vs_importance.png")

    # ─── Plot 2: Normalized mean importance (fold-enrichment over uniform) ───

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(
        data=results_df,
        x='entity_type',
        y='mean_importance_normalized',
        order=type_order,
        ax=ax
    )
    ax.axhline(y=1.0, color='red', linestyle='--', alpha=0.7, label='Uniform baseline')
    ax.set_xlabel('Relationship Type')
    ax.set_ylabel('Normalized Mean Importance\n(fold over uniform)')
    ax.set_title('Per-Feature Importance by Relationship Type\n(normalized for feature count)')
    ax.legend()
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    fig.savefig(plots_dir / 'feature_importance_by_type.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("✓ Saved plots/feature_importance_by_type.png")

    # ─── Plot 3: Normalized mean importance, faceted by KG ───────────────────

    g = sns.catplot(
        data=results_df,
        x='entity_type',
        y='mean_importance_normalized',
        col='kg',
        kind='box',
        order=type_order,
        col_wrap=2,
        height=4,
        aspect=1.5,
        sharey=True
    )
    for ax in g.axes.flat:
        ax.axhline(y=1.0, color='red', linestyle='--', alpha=0.7)
    g.set_xticklabels(rotation=45, ha='right')
    g.set_axis_labels('Relationship Type', 'Normalized Mean Importance')
    g.figure.suptitle('Per-Feature Importance by Relationship Type per KG', y=1.02)
    g.tight_layout()
    g.savefig(plots_dir / 'feature_importance_by_type_per_kg.png', dpi=150, bbox_inches='tight')
    plt.close(g.figure)
    print("✓ Saved plots/feature_importance_by_type_per_kg.png")

    # ─── Plot 4: Normalized mean importance, faceted by indication ───────────

    g2 = sns.catplot(
        data=results_df,
        x='entity_type',
        y='mean_importance_normalized',
        col='indication',
        kind='box',
        order=type_order,
        col_wrap=4,
        height=3.5,
        aspect=1.3,
        sharey=True
    )
    for ax in g2.axes.flat:
        ax.axhline(y=1.0, color='red', linestyle='--', alpha=0.7)
    g2.set_xticklabels(rotation=45, ha='right')
    g2.set_axis_labels('Relationship Type', 'Normalized Mean Importance')
    g2.figure.suptitle('Per-Feature Importance by Relationship Type per Indication', y=1.02)
    g2.tight_layout()
    g2.savefig(plots_dir / 'feature_importance_by_type_per_indication.png', dpi=150, bbox_inches='tight')
    plt.close(g2.figure)
    print("✓ Saved plots/feature_importance_by_type_per_indication.png")

    # Print summary table
    print("\n=== Summary: Median Normalized Mean Importance (fold over uniform) ===")
    summary = results_df.pivot_table(
        values='mean_importance_normalized',
        index='entity_type',
        columns='kg',
        aggfunc='median'
    ).round(2)
    summary = summary.loc[type_order]
    print(summary.to_string())


if __name__ == '__main__':
    main()
