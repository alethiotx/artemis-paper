#!/usr/bin/env python
"""
Negative Sampling Sensitivity Analysis

Quantifies how sensitive model predictions are to the random draw of negative
samples by measuring variance across the 10 training iterations (each using
a different random seed and thus different negatives).

Analyses:
1. AUROC variance across iterations per KG × embedding × indication
   (reads existing CV results from S3 — no retraining needed)
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

import gc

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier

from alethiotx.artemis.cv.pipeline import prepare as prepare_model
from alethiotx.artemis.clinical.scores import load as load_clinical_scores
from alethiotx.artemis.clinical.scores import all_targets


# ─── Configuration ───────────────────────────────────────────────────────────

INDICATIONS = ['breast', 'lung', 'bowel', 'prostate', 'melanoma', 'diabetes', 'cardiovascular']
KGS = ['hetionet', 'biokg', 'openbiolink', 'primekg']
EMBEDDINGS = ['ComplEx', 'DistMult', 'RotatE', 'TransE']
N_ITERATIONS = 10
RF_THRESHOLD = 0.5

# S3 paths
CV_DATA_PREFIX = 's3://alethiotx-artemis/figs_review/cv/data'
KG_DATA_PREFIX = 's3://alethiotx-artemis/data/kgs-no-data-leakage/associations'


# ─── Helper Functions ────────────────────────────────────────────────────────

def load_kg_features(kg: str, embedding: str) -> pd.DataFrame:
    """Load knowledge graph features from S3."""
    return pd.read_parquet(f'{KG_DATA_PREFIX}/{kg}/{embedding}/summarize/predictions.parquet')


def load_clinical_data(scores_date: str, kg_features: pd.DataFrame):
    """Load clinical scores for all indications (All filter only)."""
    clinical_data_raw = load_clinical_scores(date=scores_date)
    clinical_data = {}
    for indication, scores in zip(INDICATIONS, clinical_data_raw):
        clinical_data[indication] = scores[scores['Target Gene'].isin(kg_features.index)]
    known_targets = all_targets(list(clinical_data_raw))
    return clinical_data, known_targets


# ─── Analysis 1: AUROC Variance (from existing CV results) ──────────────────

def load_auroc_variance():
    """Read existing CV binary results from S3 and compute per-iteration variance."""
    results = []

    for kg in KGS:
        for embedding in EMBEDDINGS:
            for indication in INDICATIONS:
                csv_path = f'{CV_DATA_PREFIX}/{kg}_{embedding}_{indication}_binary.csv'
                try:
                    df = pd.read_csv(csv_path)
                except FileNotFoundError:
                    print(f"    Missing: {kg}/{embedding}/{indication}, skipping")
                    continue

                # Filter to Real targets, Random Forest, roc_auc
                real_rf_auroc = df[
                    (df['targets'] == 'Real') &
                    (df['classifier'] == 'Random Forest') &
                    (df['scoring'] == 'roc_auc')
                ]['score'].values

                if len(real_rf_auroc) == 0:
                    continue

                results.append({
                    'kg': kg,
                    'embedding': embedding,
                    'indication': indication,
                    'auroc_mean': np.mean(real_rf_auroc),
                    'auroc_std': np.std(real_rf_auroc),
                    'auroc_min': np.min(real_rf_auroc),
                    'auroc_max': np.max(real_rf_auroc),
                    'auroc_range': np.max(real_rf_auroc) - np.min(real_rf_auroc),
                    'n_iterations': len(real_rf_auroc),
                })

    return pd.DataFrame(results)


# ─── Analysis 2 & 3: Gene Stability + Jaccard (combined) ────────────────────

def compute_stability_and_jaccard(kg_features, clinical_data, known_targets,
                                   kg, embedding):
    """
    Train models once per iteration, extract both gene stability and Jaccard
    similarity from the same predictions.
    """
    stability_rows = []
    jaccard_rows = []

    for indication in INDICATIONS:
        gene_hits = pd.Series(0, index=kg_features.index, dtype=int)
        target_sets = []

        for iteration in range(1, N_ITERATIONS + 1):
            training_data = prepare_model(
                kg_features, clinical_data[indication],
                known_targets=known_targets,
                rand_seed=iteration
            )

            clf = RandomForestClassifier(random_state=iteration)
            clf.fit(training_data['X'], training_data['y_binary'])

            probs = clf.predict_proba(kg_features)[:, 1]
            predicted_mask = probs >= RF_THRESHOLD
            gene_hits += predicted_mask.astype(int)
            target_sets.append(set(kg_features.index[predicted_mask]))

        # Gene stability
        stability = gene_hits / N_ITERATIONS
        nonzero = stability[stability > 0]
        stability_rows.append({
            'kg': kg,
            'embedding': embedding,
            'indication': indication,
            'n_genes_predicted': len(nonzero),
            'n_stable_genes': int((stability >= 0.8).sum()),
            'n_unstable_genes': int(((stability > 0) & (stability < 0.5)).sum()),
            'mean_stability': float(nonzero.mean()) if len(nonzero) > 0 else 0,
        })

        # Pairwise Jaccard
        n = len(target_sets)
        jaccard_matrix = np.ones((n, n))
        for i, j in combinations(range(n), 2):
            intersection = len(target_sets[i] & target_sets[j])
            union = len(target_sets[i] | target_sets[j])
            jac = intersection / union if union > 0 else 1.0
            jaccard_matrix[i, j] = jac
            jaccard_matrix[j, i] = jac

        mask = ~np.eye(n, dtype=bool)
        jaccard_rows.append({
            'kg': kg,
            'embedding': embedding,
            'indication': indication,
            'mean_jaccard': float(jaccard_matrix[mask].mean()),
            'min_jaccard': float(jaccard_matrix[mask].min()),
        })

    return stability_rows, jaccard_rows


# ─── Plotting ────────────────────────────────────────────────────────────────

def plot_auroc_variance(df, output_dir):
    """Plot AUROC mean ± SD across iterations."""
    fig, ax = plt.subplots(figsize=(14, 8))
    df['label'] = df['kg'] + ' / ' + df['embedding']

    for indication in INDICATIONS:
        subset = df[df['indication'] == indication]
        if subset.empty:
            continue
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


def plot_gene_stability(stability_df, output_dir):
    """Plot summary of gene prediction stability across all KG/embedding combos."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: mean stability by indication
    pivot = stability_df.pivot_table(values='mean_stability', index='indication',
                                      columns=['kg', 'embedding'], aggfunc='mean')
    sns.heatmap(pivot, annot=True, fmt='.2f', cmap='YlOrRd', vmin=0, vmax=1, ax=axes[0])
    axes[0].set_title('Mean Gene Prediction Stability')

    # Right: stable genes count
    pivot2 = stability_df.pivot_table(values='n_stable_genes', index='indication',
                                       columns=['kg', 'embedding'], aggfunc='mean')
    sns.heatmap(pivot2, annot=True, fmt='.0f', cmap='Blues', ax=axes[1])
    axes[1].set_title('Genes Stable in ≥80% Iterations')

    plt.tight_layout()
    plt.savefig(output_dir / 'gene_stability.png', dpi=300)
    plt.close()


def plot_jaccard_summary(jaccard_df, output_dir):
    """Plot Jaccard similarity summary across all combos."""
    fig, ax = plt.subplots(figsize=(10, 6))
    pivot = jaccard_df.pivot_table(values='mean_jaccard', index='indication',
                                    columns=['kg', 'embedding'], aggfunc='mean')
    sns.heatmap(pivot, annot=True, fmt='.3f', cmap='YlOrRd', vmin=0, vmax=1, ax=ax)
    ax.set_title('Mean Pairwise Jaccard Similarity of Predicted Targets\n(across 10 iterations)')
    plt.tight_layout()
    plt.savefig(output_dir / 'jaccard_summary.png', dpi=300)
    plt.close()


# ─── Main Execution ──────────────────────────────────────────────────────────

def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)

    scores_date = sys.argv[1]

    data_dir = Path('data')
    plots_dir = Path('plots')
    data_dir.mkdir(exist_ok=True)
    plots_dir.mkdir(exist_ok=True)

    # ── Analysis 1: AUROC variance from existing CV results (fast) ──
    print("=== Analysis 1: AUROC Variance (from existing CV results) ===")
    auroc_df = load_auroc_variance()
    auroc_df.to_csv(data_dir / 'auroc_variance.csv', index=False)
    plot_auroc_variance(auroc_df, plots_dir)
    print(f"  Mean AUROC SD: {auroc_df['auroc_std'].mean():.4f}")
    print(f"  Max AUROC SD: {auroc_df['auroc_std'].max():.4f}")

    # ── Analysis 2 & 3: Gene stability + Jaccard (combined, single training pass) ──
    print("\n=== Analysis 2 & 3: Gene Stability + Jaccard ===")
    all_stability_rows = []
    all_jaccard_rows = []

    for kg in KGS:
        for embedding in EMBEDDINGS:
            print(f"  {kg}/{embedding}...")
            kg_features = load_kg_features(kg, embedding)
            clinical_data, known_targets = load_clinical_data(scores_date, kg_features)

            stability_rows, jaccard_rows = compute_stability_and_jaccard(
                kg_features, clinical_data, known_targets, kg, embedding
            )
            all_stability_rows.extend(stability_rows)
            all_jaccard_rows.extend(jaccard_rows)

            del kg_features, clinical_data, known_targets
            gc.collect()

    stability_df = pd.DataFrame(all_stability_rows)
    stability_df.to_csv(data_dir / 'gene_stability.csv', index=False)
    print(f"\n  Mean stability: {stability_df['mean_stability'].mean():.3f}")
    print(f"  Mean stable genes (≥80% iters): {stability_df['n_stable_genes'].mean():.0f}")

    jaccard_df = pd.DataFrame(all_jaccard_rows)
    jaccard_df.to_csv(data_dir / 'jaccard_similarity.csv', index=False)
    print(f"  Mean Jaccard: {jaccard_df['mean_jaccard'].mean():.3f}")
    print(f"  Min Jaccard: {jaccard_df['min_jaccard'].min():.3f}")

    plot_gene_stability(stability_df, plots_dir)
    plot_jaccard_summary(jaccard_df, plots_dir)

    print("\n✓ Sensitivity analysis complete.")


if __name__ == '__main__':
    main()
