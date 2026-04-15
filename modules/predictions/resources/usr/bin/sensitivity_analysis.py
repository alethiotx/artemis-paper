#!/usr/bin/env python
"""
Negative Sampling Sensitivity Analysis

Quantifies how sensitive model predictions are to the random draw of negative
samples by measuring variance across the 10 training iterations (each using
a different random seed and thus different negatives).

Analyses:
1. AUROC variance across iterations per KG × embedding × indication
2. Per-gene prediction stability: fraction of iterations predicting each gene as target
3. Core target set consistency: Jaccard similarity across iteration pairs

Usage:
    sensitivity_analysis.py <scores_date>

Arguments:
    scores_date: Clinical scores date (YYYY-MM-DD format)

Output:
    - data/auroc_variance.csv: AUROC mean ± SD across iterations
    - data/gene_stability.csv: Per-gene prediction frequency across iterations
    - data/jaccard_similarity.csv: Pairwise Jaccard between iteration target sets
    - plots/auroc_variance.png: AUROC variance by KG × embedding × indication
    - plots/gene_stability.png: Distribution of gene prediction stability
    - plots/jaccard_heatmap.png: Pairwise Jaccard similarity heatmap
"""

import sys
from pathlib import Path
from itertools import combinations

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
RANDOM_STATE = 42


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


# ─── Analysis 1: AUROC Variance ─────────────────────────────────────────────

def compute_auroc_variance(kg_features, clinical_data, known_targets, pathway_genes_map,
                           kg, embedding, ct_filter, pg_number):
    """Compute AUROC across iterations for each indication."""
    results = []

    for indication in INDICATIONS:
        scores_per_iter = []

        for iteration in range(1, N_ITERATIONS + 1):
            prepare_kwargs = {
                'known_targets': known_targets,
                'rand_seed': iteration
            }
            pg = pathway_genes_map[indication]
            if pg:
                prepare_kwargs['pathway_genes'] = pg

            training_data = prepare_model(
                kg_features, clinical_data[indication], **prepare_kwargs
            )

            clf = RandomForestClassifier(random_state=iteration)
            cv = StratifiedKFold()
            auroc = np.mean(cross_val_score(
                clf, training_data['X'], training_data['y_binary'],
                scoring='roc_auc', cv=cv
            ))
            scores_per_iter.append(auroc)

        results.append({
            'kg': kg,
            'embedding': embedding,
            'ct_filter': ct_filter,
            'pg_number': pg_number,
            'indication': indication,
            'auroc_mean': np.mean(scores_per_iter),
            'auroc_std': np.std(scores_per_iter),
            'auroc_min': np.min(scores_per_iter),
            'auroc_max': np.max(scores_per_iter),
            'auroc_range': np.max(scores_per_iter) - np.min(scores_per_iter),
        })

    return results


# ─── Analysis 2: Gene Prediction Stability ──────────────────────────────────

def compute_gene_stability(kg_features, clinical_data, known_targets, pathway_genes_map,
                           kg, embedding, ct_filter, rf_threshold, pg_number):
    """For each gene, compute fraction of iterations predicting it as target."""
    prediction_counts = {}

    for indication in INDICATIONS:
        gene_hits = pd.Series(0, index=kg_features.index, dtype=int)

        for iteration in range(1, N_ITERATIONS + 1):
            prepare_kwargs = {
                'known_targets': known_targets,
                'rand_seed': iteration
            }
            pg = pathway_genes_map[indication]
            if pg:
                prepare_kwargs['pathway_genes'] = pg

            training_data = prepare_model(
                kg_features, clinical_data[indication], **prepare_kwargs
            )

            clf = RandomForestClassifier(random_state=iteration)
            clf.fit(training_data['X'], training_data['y_binary'])

            probs = clf.predict_proba(kg_features)[:, 1]
            predicted = probs >= rf_threshold
            gene_hits += predicted.astype(int)

        stability = gene_hits / N_ITERATIONS
        prediction_counts[indication] = stability

    return prediction_counts


# ─── Analysis 3: Jaccard Similarity ──────────────────────────────────────────

def compute_jaccard_similarity(kg_features, clinical_data, known_targets, pathway_genes_map,
                               kg, embedding, ct_filter, rf_threshold, pg_number, indication):
    """Compute pairwise Jaccard similarity of predicted target sets across iterations."""
    target_sets = []

    for iteration in range(1, N_ITERATIONS + 1):
        prepare_kwargs = {
            'known_targets': known_targets,
            'rand_seed': iteration
        }
        pg = pathway_genes_map[indication]
        if pg:
            prepare_kwargs['pathway_genes'] = pg

        training_data = prepare_model(
            kg_features, clinical_data[indication], **prepare_kwargs
        )

        clf = RandomForestClassifier(random_state=iteration)
        clf.fit(training_data['X'], training_data['y_binary'])

        probs = clf.predict_proba(kg_features)[:, 1]
        predicted = set(kg_features.index[probs >= rf_threshold])
        target_sets.append(predicted)

    # Pairwise Jaccard
    n = len(target_sets)
    jaccard_matrix = np.ones((n, n))
    for i, j in combinations(range(n), 2):
        intersection = len(target_sets[i] & target_sets[j])
        union = len(target_sets[i] | target_sets[j])
        jac = intersection / union if union > 0 else 1.0
        jaccard_matrix[i, j] = jac
        jaccard_matrix[j, i] = jac

    return jaccard_matrix, target_sets


# ─── Plotting ────────────────────────────────────────────────────────────────

def plot_auroc_variance(df, output_dir):
    """Plot AUROC mean ± SD across iterations."""
    fig, ax = plt.subplots(figsize=(14, 8))
    df['label'] = df['kg'] + ' / ' + df['embedding']

    for i, indication in enumerate(INDICATIONS):
        subset = df[df['indication'] == indication]
        ax.errorbar(
            subset['label'], subset['auroc_mean'],
            yerr=subset['auroc_std'],
            label=indication, marker='o', capsize=3, alpha=0.7
        )

    ax.set_xlabel('KG / Embedding')
    ax.set_ylabel('AUROC (mean ± SD across 10 iterations)')
    ax.set_title('Sensitivity of AUROC to Negative Sampling (across iterations)')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_dir / 'auroc_variance.png', dpi=300)
    plt.close()


def plot_gene_stability(stability_dict, output_dir, kg, embedding, ct_filter, rf_threshold, pg_number):
    """Plot distribution of gene prediction stability."""
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()

    for i, indication in enumerate(INDICATIONS):
        ax = axes[i]
        stability = stability_dict[indication]
        # Only show genes predicted at least once
        nonzero = stability[stability > 0]
        ax.hist(nonzero, bins=10, edgecolor='black', alpha=0.7)
        ax.set_title(f'{indication}\n(n={len(nonzero)} genes predicted ≥1 iter)')
        ax.set_xlabel('Fraction of iterations predicted as target')
        ax.set_ylabel('Number of genes')
        ax.axvline(x=0.5, color='red', linestyle='--', alpha=0.5)

    # Hide unused subplot
    if len(INDICATIONS) < len(axes):
        for j in range(len(INDICATIONS), len(axes)):
            axes[j].set_visible(False)

    fig.suptitle(
        f'Gene Prediction Stability ({kg}/{embedding}, CT={ct_filter}, '
        f'RF≥{rf_threshold}, PG={pg_number})',
        fontsize=14
    )
    plt.tight_layout()
    plt.savefig(output_dir / 'gene_stability.png', dpi=300)
    plt.close()


def plot_jaccard_heatmap(jaccard_matrix, output_dir, indication, kg, embedding):
    """Plot pairwise Jaccard similarity heatmap."""
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(
        jaccard_matrix, annot=True, fmt='.3f',
        xticklabels=[f'Iter {i+1}' for i in range(jaccard_matrix.shape[0])],
        yticklabels=[f'Iter {i+1}' for i in range(jaccard_matrix.shape[0])],
        vmin=0, vmax=1, cmap='YlOrRd', ax=ax
    )
    ax.set_title(f'Jaccard Similarity of Predicted Targets\n({indication}, {kg}/{embedding})')
    plt.tight_layout()
    plt.savefig(output_dir / f'jaccard_{indication}_{kg}_{embedding}.png', dpi=300)
    plt.close()


# ─── Main Execution ──────────────────────────────────────────────────────────

def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)

    scores_date = sys.argv[1]

    # Create output directories
    data_dir = Path('data')
    plots_dir = Path('plots')
    data_dir.mkdir(exist_ok=True)
    plots_dir.mkdir(exist_ok=True)

    # ── Analysis 1: AUROC variance across full grid ──
    print("=== Analysis 1: AUROC Variance ===")
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
                    print(f"  AUROC: {kg}/{embedding} CT={ct_filter} PG={pg_number}")

                    results = compute_auroc_variance(
                        kg_features, clinical_data, known_targets, pathway_genes_map,
                        kg, embedding, ct_filter, pg_number
                    )
                    all_auroc_results.extend(results)

    auroc_df = pd.DataFrame(all_auroc_results)
    auroc_df.to_csv(data_dir / 'auroc_variance.csv', index=False)
    plot_auroc_variance(auroc_df, plots_dir)
    print(f"  Mean AUROC SD across all combos: {auroc_df['auroc_std'].mean():.4f}")
    print(f"  Max AUROC SD: {auroc_df['auroc_std'].max():.4f}")

    # ── Analysis 2 & 3: Gene stability + Jaccard (representative config) ──
    # Run on all KG × embedding combos with default params
    print("\n=== Analysis 2 & 3: Gene Stability + Jaccard ===")
    all_stability_rows = []
    all_jaccard_rows = []
    ct_filter = 'All'
    rf_threshold = 0.5
    pg_number = 0

    for kg in KGS:
        for embedding in EMBEDDINGS:
            print(f"  Stability: {kg}/{embedding}...")
            kg_features = load_kg_features(kg, embedding)
            clinical_data, known_targets = load_clinical_data(
                scores_date, ct_filter, kg_features
            )
            pathway_genes_map = load_pathway_genes_map(scores_date, pg_number)

            # Gene stability
            stability_dict = compute_gene_stability(
                kg_features, clinical_data, known_targets, pathway_genes_map,
                kg, embedding, ct_filter, rf_threshold, pg_number
            )

            for indication in INDICATIONS:
                stability = stability_dict[indication]
                nonzero = stability[stability > 0]
                all_stability_rows.append({
                    'kg': kg,
                    'embedding': embedding,
                    'indication': indication,
                    'n_genes_predicted': len(nonzero),
                    'n_stable_genes': int((stability >= 0.8).sum()),
                    'n_unstable_genes': int(((stability > 0) & (stability < 0.5)).sum()),
                    'mean_stability': float(nonzero.mean()) if len(nonzero) > 0 else 0,
                })

            # Jaccard (for all indications)
            for indication in INDICATIONS:
                jaccard_matrix, _ = compute_jaccard_similarity(
                    kg_features, clinical_data, known_targets, pathway_genes_map,
                    kg, embedding, ct_filter, rf_threshold, pg_number, indication
                )
                # Store mean off-diagonal Jaccard
                mask = ~np.eye(jaccard_matrix.shape[0], dtype=bool)
                all_jaccard_rows.append({
                    'kg': kg,
                    'embedding': embedding,
                    'indication': indication,
                    'mean_jaccard': float(jaccard_matrix[mask].mean()),
                    'min_jaccard': float(jaccard_matrix[mask].min()),
                })
                plot_jaccard_heatmap(jaccard_matrix, plots_dir, indication, kg, embedding)

            plot_gene_stability(
                stability_dict, plots_dir, kg, embedding, ct_filter, rf_threshold, pg_number
            )

    stability_df = pd.DataFrame(all_stability_rows)
    stability_df.to_csv(data_dir / 'gene_stability.csv', index=False)
    print(f"\n  Mean gene prediction stability: {stability_df['mean_stability'].mean():.3f}")
    print(f"  Mean stable genes (≥80% iters): {stability_df['n_stable_genes'].mean():.0f}")

    jaccard_df = pd.DataFrame(all_jaccard_rows)
    jaccard_df.to_csv(data_dir / 'jaccard_similarity.csv', index=False)
    print(f"  Mean Jaccard similarity: {jaccard_df['mean_jaccard'].mean():.3f}")
    print(f"  Min Jaccard similarity: {jaccard_df['min_jaccard'].min():.3f}")

    print("\n✓ Sensitivity analysis complete.")


if __name__ == '__main__':
    main()
