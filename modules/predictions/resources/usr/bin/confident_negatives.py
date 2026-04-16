#!/usr/bin/env python
"""
Confident Negative Selection Experiment

Two-step negative sampling robustness check:
  Step 1: Train RF with standard balanced random negatives (existing approach)
  Step 2: Score all unlabeled genes; keep only low-probability ones as "reliable" negatives
  Step 3: Retrain RF with positives + reliable negatives; compare AUROC & predictions

Runs across the full parameter grid (all KGs × embeddings × CT filters × RF thresholds
× pathway gene counts × iterations).

Usage:
    confident_negatives.py <scores_date>

Arguments:
    scores_date: Clinical scores date (YYYY-MM-DD format)

Output:
    - data/confident_neg_auroc.csv: AUROC comparison (standard vs confident negatives)
    - data/confident_neg_overlap.csv: Target overlap comparison
    - plots/auroc_comparison.png: Side-by-side AUROC
    - plots/overlap_comparison.png: Prediction overlap change
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score

from alethiotx.artemis.cv.pipeline import prepare as prepare_model
from alethiotx.artemis.clinical.scores import load as load_clinical_scores
from alethiotx.artemis.clinical.scores import approved as approved_clinical_scores
from alethiotx.artemis.clinical.scores import unique as unique_clinical_scores
from alethiotx.artemis.clinical.scores import all_targets
from alethiotx.artemis.pathway.genes import load as load_pathway_genes
from alethiotx.artemis.pathway.genes import unique as unique_pathway_genes


# ─── Configuration ───────────────────────────────────────────────────────────

INDICATIONS = ['breast', 'lung', 'bowel', 'prostate', 'melanoma', 'diabetes', 'cardiovascular']
KGS = ['hetionet', 'biokg', 'openbiolink', 'primekg']
EMBEDDINGS = ['ComplEx', 'DistMult', 'RotatE', 'TransE']
CT_FILTERS = ['All', 'Unique', 'Approved']
RF_THRESHOLDS = [0.5, 0.6, 0.7, 0.8, 0.9]
PG_NUMBERS = [0, 100, 300]
N_ITERATIONS = 10
MAX_CV_FOLDS = 5
CONFIDENT_NEG_QUANTILE = 0.25  # Bottom 25% of predicted probabilities


# ─── Helper Functions ────────────────────────────────────────────────────────

def load_kg_features(kg: str, embedding: str) -> pd.DataFrame:
    """Load knowledge graph features from S3."""
    kg_path = f's3://alethiotx-artemis/data/kgs-no-data-leakage/associations/{kg}/{embedding}/summarize/predictions.parquet'
    return pd.read_parquet(kg_path)


def load_clinical_data(scores_date: str, ct_filter: str, kg_features: pd.DataFrame):
    """Load and filter clinical scores for all indications."""
    clinical_data_raw = load_clinical_scores(date=scores_date)

    if ct_filter == 'Approved':
        clinical_data_raw = approved_clinical_scores(list(clinical_data_raw))
    elif ct_filter == 'Unique':
        clinical_data_raw = unique_clinical_scores(list(clinical_data_raw))

    clinical_data = {}
    for indication, scores in zip(INDICATIONS, clinical_data_raw):
        clinical_data[indication] = scores[scores['Target Gene'].isin(kg_features.index)]

    known_targets = all_targets(list(clinical_data_raw))
    return clinical_data, known_targets


def load_pathway_genes_map(scores_date: str, pg_number: int):
    """Load pathway genes for all indications."""
    if pg_number == 0:
        return {ind: [] for ind in INDICATIONS}
    pathway_genes = load_pathway_genes(date=scores_date, n=pg_number)
    pathway_genes = unique_pathway_genes(list(pathway_genes))
    return dict(zip(INDICATIONS, pathway_genes))


def prepare_confident_negatives(kg_features, clinical_scores, known_targets,
                                 pathway_genes, rand_seed):
    """
    Two-step confident negative selection.

    Step 1: Standard prepare_model → train initial RF
    Step 2: Score unlabeled genes → select reliable negatives → re-prepare
    """
    # Step 1: Standard training
    prepare_kwargs = {'known_targets': known_targets, 'rand_seed': rand_seed}
    if pathway_genes:
        prepare_kwargs['pathway_genes'] = pathway_genes

    initial_data = prepare_model(kg_features, clinical_scores, **prepare_kwargs)

    initial_clf = RandomForestClassifier(random_state=rand_seed)
    initial_clf.fit(initial_data['X'], initial_data['y_binary'])

    # Step 2: Score all KG genes
    all_probs = initial_clf.predict_proba(kg_features)[:, 1]
    prob_series = pd.Series(all_probs, index=kg_features.index)

    # Identify unlabeled genes (not positive, not known target, not pathway)
    positive_genes = clinical_scores['Target Gene'].tolist()
    pg_list = pathway_genes if pathway_genes else []
    unlabeled_mask = (
        ~kg_features.index.isin(positive_genes) &
        ~kg_features.index.isin(known_targets) &
        ~kg_features.index.isin(pg_list)
    )
    unlabeled_probs = prob_series[unlabeled_mask]

    # Select reliable negatives: bottom quantile of predicted probability
    threshold = unlabeled_probs.quantile(CONFIDENT_NEG_QUANTILE)
    reliable_negatives = unlabeled_probs[unlabeled_probs <= threshold].index.tolist()

    # Step 3: Re-prepare with restricted negative pool
    # We replicate the logic of prepare_model but restrict the negative pool
    from numpy import log2, random as np_random
    random_state = np.random.RandomState(rand_seed)

    y = clinical_scores[['Target Gene', 'Clinical Score']].copy()
    y['Clinical Score'] = log2(y['Clinical Score'] + 1)
    y.index = y['Target Gene']
    y = y.drop(columns=['Target Gene'])
    y = y.join(kg_features, how='right')

    # Positives
    y_pos = y[~y['Clinical Score'].isna()]
    try:
        y_pos = y_pos.copy()
        y_pos['Clinical Score Binned'] = pd.qcut(
            y_pos['Clinical Score'], q=3,
            labels=['target_0', 'target_1', 'target_2']
        )
    except Exception:
        y_pos = y_pos.copy()
        y_pos['Clinical Score Binned'] = 'target'

    # Pathway genes
    y_pg = pd.DataFrame()
    if pathway_genes:
        pg_filtered = [g for g in pathway_genes if g not in y_pos.index.tolist() and g not in known_targets]
        y_pg = y[y.index.isin(pg_filtered)].copy()
        y_pg['Clinical Score'] = 1
        y_pg['Clinical Score Binned'] = 'pathway_gene'

    n_neg_needed = y_pos.shape[0] + y_pg.shape[0]

    # Sample from reliable negatives only
    reliable_neg_pool = y[y.index.isin(reliable_negatives)]
    if len(reliable_neg_pool) < n_neg_needed:
        # Fall back to all available if not enough reliable negatives
        y_neg = reliable_neg_pool.copy()
    else:
        y_neg = reliable_neg_pool.sample(n_neg_needed, random_state=rand_seed)

    y_neg = y_neg.copy()
    y_neg['Clinical Score'] = -1
    y_neg['Clinical Score Binned'] = 'not_target'

    y_combined = pd.concat([y_pos, y_neg, y_pg])
    y_combined = y_combined.sample(frac=1, random_state=rand_seed)

    X_out = y_combined.iloc[:, 1:-1]
    y_out = y_combined['Clinical Score']
    y_binary = (y_out > 0).astype(int)

    return {
        'X': X_out,
        'y': y_out,
        'y_binary': y_binary,
        'n_reliable_negatives': len(reliable_negatives),
        'confident_threshold': threshold,
    }


# ─── Main Analysis ───────────────────────────────────────────────────────────

def run_comparison(kg_features, clinical_data, known_targets, pathway_genes_map,
                   kg, embedding, ct_filter, pg_number):
    """Run standard vs confident negative comparison for all indications."""
    results = []

    for indication in INDICATIONS:
        standard_aurocs = []
        confident_aurocs = []

        for iteration in range(1, N_ITERATIONS + 1):
            # Standard approach
            prepare_kwargs = {'known_targets': known_targets, 'rand_seed': iteration}
            pg = pathway_genes_map[indication]
            if pg:
                prepare_kwargs['pathway_genes'] = pg

            standard_data = prepare_model(
                kg_features, clinical_data[indication], **prepare_kwargs
            )
            clf = RandomForestClassifier(random_state=iteration)
            min_class_std = int(standard_data['y_binary'].value_counts().min())
            n_splits = min(MAX_CV_FOLDS, min_class_std)
            if n_splits < 2:
                print(f"    Skipping {indication} iter {iteration}: too few samples (min_class={min_class_std})")
                continue
            cv = StratifiedKFold(n_splits=n_splits)
            std_auroc = np.mean(cross_val_score(
                clf, standard_data['X'], standard_data['y_binary'],
                scoring='roc_auc', cv=cv
            ))
            standard_aurocs.append(std_auroc)

            # Confident negative approach
            confident_data = prepare_confident_negatives(
                kg_features, clinical_data[indication], known_targets,
                pg, iteration
            )
            clf2 = RandomForestClassifier(random_state=iteration)
            min_class_conf = int(confident_data['y_binary'].value_counts().min())
            n_splits_conf = min(MAX_CV_FOLDS, min_class_conf)
            if n_splits_conf < 2:
                print(f"    Skipping confident {indication} iter {iteration}: too few samples")
                continue
            cv2 = StratifiedKFold(n_splits=n_splits_conf)
            conf_auroc = np.mean(cross_val_score(
                clf2, confident_data['X'], confident_data['y_binary'],
                scoring='roc_auc', cv=cv2
            ))
            confident_aurocs.append(conf_auroc)

        if not standard_aurocs or not confident_aurocs:
            print(f"    No valid iterations for {indication}, skipping")
            continue

        results.append({
            'kg': kg,
            'embedding': embedding,
            'ct_filter': ct_filter,
            'pg_number': pg_number,
            'indication': indication,
            'standard_auroc_mean': np.mean(standard_aurocs),
            'standard_auroc_std': np.std(standard_aurocs),
            'confident_auroc_mean': np.mean(confident_aurocs),
            'confident_auroc_std': np.std(confident_aurocs),
            'auroc_diff': np.mean(confident_aurocs) - np.mean(standard_aurocs),
        })

    return results


def run_overlap_comparison(kg_features, clinical_data, known_targets, pathway_genes_map,
                           kg, embedding, ct_filter, rf_threshold, pg_number):
    """Compare predicted target sets between standard and confident approaches."""
    results = []

    for indication in INDICATIONS:
        standard_targets_all = set()
        confident_targets_all = set()

        for iteration in range(1, N_ITERATIONS + 1):
            prepare_kwargs = {'known_targets': known_targets, 'rand_seed': iteration}
            pg = pathway_genes_map[indication]
            if pg:
                prepare_kwargs['pathway_genes'] = pg

            # Standard
            standard_data = prepare_model(
                kg_features, clinical_data[indication], **prepare_kwargs
            )
            clf = RandomForestClassifier(random_state=iteration)
            clf.fit(standard_data['X'], standard_data['y_binary'])
            std_probs = clf.predict_proba(kg_features)[:, 1]
            std_targets = set(kg_features.index[std_probs >= rf_threshold])
            standard_targets_all |= std_targets

            # Confident
            confident_data = prepare_confident_negatives(
                kg_features, clinical_data[indication], known_targets,
                pg, iteration
            )
            clf2 = RandomForestClassifier(random_state=iteration)
            clf2.fit(confident_data['X'], confident_data['y_binary'])
            conf_probs = clf2.predict_proba(kg_features)[:, 1]
            conf_targets = set(kg_features.index[conf_probs >= rf_threshold])
            confident_targets_all |= conf_targets

        # Compute overlap
        intersection = standard_targets_all & confident_targets_all
        union = standard_targets_all | confident_targets_all
        jaccard = len(intersection) / len(union) if union else 1.0

        results.append({
            'kg': kg,
            'embedding': embedding,
            'ct_filter': ct_filter,
            'rf_threshold': rf_threshold,
            'pg_number': pg_number,
            'indication': indication,
            'n_standard_targets': len(standard_targets_all),
            'n_confident_targets': len(confident_targets_all),
            'n_shared_targets': len(intersection),
            'jaccard': jaccard,
        })

    return results


# ─── Plotting ────────────────────────────────────────────────────────────────

def plot_auroc_comparison(df, output_dir):
    """Plot AUROC comparison between standard and confident negative approaches."""
    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    axes = axes.flatten()

    for i, indication in enumerate(INDICATIONS):
        ax = axes[i]
        subset = df[df['indication'] == indication].copy()
        subset['label'] = subset['kg'] + '\n' + subset['embedding']

        x = np.arange(len(subset))
        width = 0.35

        ax.bar(x - width/2, subset['standard_auroc_mean'],
               width, yerr=subset['standard_auroc_std'],
               label='Standard', alpha=0.8, capsize=2)
        ax.bar(x + width/2, subset['confident_auroc_mean'],
               width, yerr=subset['confident_auroc_std'],
               label='Confident Neg.', alpha=0.8, capsize=2)

        ax.set_title(indication)
        ax.set_ylabel('AUROC')
        ax.set_xticks(x)
        ax.set_xticklabels(subset['label'], rotation=45, ha='right', fontsize=6)
        if i == 0:
            ax.legend(fontsize=8)

    if len(INDICATIONS) < len(axes):
        for j in range(len(INDICATIONS), len(axes)):
            axes[j].set_visible(False)

    fig.suptitle('AUROC: Standard vs Confident Negative Sampling', fontsize=14)
    plt.tight_layout()
    plt.savefig(output_dir / 'auroc_comparison.png', dpi=300)
    plt.close()


def plot_overlap_comparison(df, output_dir):
    """Plot Jaccard overlap between standard and confident approaches."""
    fig, ax = plt.subplots(figsize=(12, 6))

    pivot = df.pivot_table(
        values='jaccard',
        index='indication',
        columns=['kg', 'embedding'],
        aggfunc='mean'
    )

    sns.heatmap(pivot, annot=True, fmt='.2f', cmap='YlOrRd',
                vmin=0, vmax=1, ax=ax)
    ax.set_title('Jaccard Overlap: Standard vs Confident Negative Target Sets')
    plt.tight_layout()
    plt.savefig(output_dir / 'overlap_comparison.png', dpi=300)
    plt.close()


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)

    scores_date = sys.argv[1]

    data_dir = Path('data')
    plots_dir = Path('plots')
    data_dir.mkdir(exist_ok=True)
    plots_dir.mkdir(exist_ok=True)

    # ── AUROC comparison across full grid ──
    print("=== AUROC Comparison: Standard vs Confident Negatives ===")
    all_auroc_results = []

    for kg in KGS:
        for embedding in EMBEDDINGS:
            print(f"Loading {kg}/{embedding}...")
            kg_features = load_kg_features(kg, embedding)

            for ct_filter in CT_FILTERS:
                clinical_data, known_targets = load_clinical_data(
                    scores_date, ct_filter, kg_features
                )

                for pg_number in PG_NUMBERS:
                    pathway_genes_map = load_pathway_genes_map(scores_date, pg_number)
                    print(f"  Comparing: CT={ct_filter} PG={pg_number}")

                    results = run_comparison(
                        kg_features, clinical_data, known_targets, pathway_genes_map,
                        kg, embedding, ct_filter, pg_number
                    )
                    all_auroc_results.extend(results)

    auroc_df = pd.DataFrame(all_auroc_results)
    auroc_df.to_csv(data_dir / 'confident_neg_auroc.csv', index=False)
    plot_auroc_comparison(auroc_df, plots_dir)

    print(f"\n  Mean AUROC diff (confident - standard): {auroc_df['auroc_diff'].mean():.4f}")
    print(f"  Max |AUROC diff|: {auroc_df['auroc_diff'].abs().max():.4f}")

    # ── Overlap comparison across full grid ──
    print("\n=== Target Overlap Comparison ===")
    all_overlap_results = []

    for kg in KGS:
        for embedding in EMBEDDINGS:
            print(f"Loading {kg}/{embedding}...")
            kg_features = load_kg_features(kg, embedding)

            for ct_filter in CT_FILTERS:
                clinical_data, known_targets = load_clinical_data(
                    scores_date, ct_filter, kg_features
                )

                for pg_number in PG_NUMBERS:
                    pathway_genes_map = load_pathway_genes_map(scores_date, pg_number)

                    for rf_threshold in RF_THRESHOLDS:
                        print(f"  Overlap: CT={ct_filter} PG={pg_number} RF={rf_threshold}")

                        results = run_overlap_comparison(
                            kg_features, clinical_data, known_targets, pathway_genes_map,
                            kg, embedding, ct_filter, rf_threshold, pg_number
                        )
                        all_overlap_results.extend(results)

    overlap_df = pd.DataFrame(all_overlap_results)
    overlap_df.to_csv(data_dir / 'confident_neg_overlap.csv', index=False)
    plot_overlap_comparison(overlap_df, plots_dir)

    print(f"\n  Mean Jaccard overlap: {overlap_df['jaccard'].mean():.3f}")
    print(f"  Min Jaccard overlap: {overlap_df['jaccard'].min():.3f}")

    print("\n✓ Confident negative analysis complete.")


if __name__ == '__main__':
    main()
