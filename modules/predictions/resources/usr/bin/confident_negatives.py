#!/usr/bin/env python
"""
Confident Negative Selection Experiment

Two-step negative sampling robustness check:
  Step 1: Train RF with standard balanced random negatives (existing approach)
  Step 2: Score all unlabeled genes; keep only low-probability ones as "reliable" negatives
  Step 3: Retrain RF with positives + reliable negatives; compare predictions

For each KG × embedding × indication × iteration, trains standard and confident
models side-by-side and compares predicted target sets via Jaccard similarity.

Usage:
    confident_negatives.py <scores_date>

Arguments:
    scores_date: Clinical scores date (YYYY-MM-DD format)

Output:
    - data/confident_neg_comparison.csv: Per-iteration prediction comparison
    - data/confident_neg_summary.csv: Summary statistics per KG × embedding × indication
    - plots/jaccard_comparison.png: Jaccard similarity heatmap
    - plots/target_count_comparison.png: Target count scatter
"""

import sys
from pathlib import Path

import gc

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from numpy import log2
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
CONFIDENT_NEG_QUANTILE = 0.25

KG_DATA_PREFIX = 's3://alethiotx-artemis/data/kgs-no-data-leakage/associations'


# ─── Helper Functions ────────────────────────────────────────────────────────

def load_kg_features(kg: str, embedding: str) -> pd.DataFrame:
    return pd.read_parquet(f'{KG_DATA_PREFIX}/{kg}/{embedding}/summarize/predictions.parquet')


def load_clinical_data(scores_date: str, kg_features: pd.DataFrame):
    clinical_data_raw = load_clinical_scores(date=scores_date)
    clinical_data = {}
    for indication, scores in zip(INDICATIONS, clinical_data_raw):
        clinical_data[indication] = scores[scores['Target Gene'].isin(kg_features.index)]
    known_targets = all_targets(list(clinical_data_raw))
    return clinical_data, known_targets


def prepare_confident_negatives(kg_features, clinical_scores, known_targets, rand_seed):
    """
    Two-step confident negative selection.

    1. Train initial RF with standard random negatives
    2. Score unlabeled genes, keep only bottom-quantile as reliable negatives
    3. Re-prepare training data from reliable negatives only
    """
    # Step 1: Standard training
    initial_data = prepare_model(
        kg_features, clinical_scores,
        known_targets=known_targets, rand_seed=rand_seed
    )
    initial_clf = RandomForestClassifier(random_state=rand_seed, n_jobs=1)
    initial_clf.fit(initial_data['X'], initial_data['y_binary'])
    del initial_data

    # Step 2: Score all KG genes
    all_probs = initial_clf.predict_proba(kg_features)[:, 1]
    del initial_clf
    prob_series = pd.Series(all_probs, index=kg_features.index)
    del all_probs

    positive_genes = clinical_scores['Target Gene'].tolist()
    unlabeled_mask = (
        ~kg_features.index.isin(positive_genes) &
        ~kg_features.index.isin(known_targets)
    )
    unlabeled_probs = prob_series[unlabeled_mask]
    del prob_series

    threshold = unlabeled_probs.quantile(CONFIDENT_NEG_QUANTILE)
    reliable_negatives = unlabeled_probs[unlabeled_probs <= threshold].index.tolist()
    del unlabeled_probs
    gc.collect()

    # Step 3: Re-prepare with restricted negative pool
    y = clinical_scores[['Target Gene', 'Clinical Score']].copy()
    y['Clinical Score'] = log2(y['Clinical Score'] + 1)
    y.index = y['Target Gene']
    y = y.drop(columns=['Target Gene'])
    y = y.join(kg_features, how='right')

    y_pos = y[~y['Clinical Score'].isna()].copy()
    n_neg_needed = y_pos.shape[0]

    reliable_neg_pool = y[y.index.isin(reliable_negatives)]
    if len(reliable_neg_pool) < n_neg_needed:
        y_neg = reliable_neg_pool.copy()
    else:
        y_neg = reliable_neg_pool.sample(n_neg_needed, random_state=rand_seed)

    y_neg = y_neg.copy()
    y_neg['Clinical Score'] = -1

    y_combined = pd.concat([y_pos, y_neg])
    y_combined = y_combined.sample(frac=1, random_state=rand_seed)

    X_out = y_combined.drop(columns=['Clinical Score'])
    y_out = y_combined['Clinical Score']
    y_binary = (y_out > 0).astype(int)

    return {'X': X_out, 'y_binary': y_binary}


# ─── Main Analysis ───────────────────────────────────────────────────────────

def run_comparison(kg_features, clinical_data, known_targets, kg, embedding):
    """
    For each indication × iteration, train standard + confident models once,
    compare predicted target sets.
    """
    rows = []

    for indication in INDICATIONS:
        for iteration in range(1, N_ITERATIONS + 1):
            # Standard model
            std_data = prepare_model(
                kg_features, clinical_data[indication],
                known_targets=known_targets, rand_seed=iteration
            )
            std_clf = RandomForestClassifier(random_state=iteration, n_jobs=1)
            std_clf.fit(std_data['X'], std_data['y_binary'])
            del std_data
            std_probs = std_clf.predict_proba(kg_features)[:, 1]
            del std_clf
            std_targets = set(kg_features.index[std_probs >= RF_THRESHOLD])

            # Confident negative model
            conf_data = prepare_confident_negatives(
                kg_features, clinical_data[indication], known_targets, iteration
            )
            conf_clf = RandomForestClassifier(random_state=iteration, n_jobs=1)
            conf_clf.fit(conf_data['X'], conf_data['y_binary'])
            del conf_data
            conf_probs = conf_clf.predict_proba(kg_features)[:, 1]
            del conf_clf
            conf_targets = set(kg_features.index[conf_probs >= RF_THRESHOLD])

            # Compare
            intersection = std_targets & conf_targets
            union = std_targets | conf_targets
            jaccard = len(intersection) / len(union) if union else 1.0

            # Probability correlation
            prob_corr = np.corrcoef(std_probs, conf_probs)[0, 1]
            del std_probs, conf_probs
            gc.collect()

            rows.append({
                'kg': kg,
                'embedding': embedding,
                'indication': indication,
                'iteration': iteration,
                'n_std_targets': len(std_targets),
                'n_conf_targets': len(conf_targets),
                'n_shared': len(intersection),
                'jaccard': jaccard,
                'prob_correlation': prob_corr,
            })

    return rows


# ─── Plotting ────────────────────────────────────────────────────────────────

def plot_jaccard(df, output_dir):
    """Heatmap of mean Jaccard similarity and recovery rate per KG × embedding × indication."""
    df = df.copy()
    df['recovery'] = df['n_shared'] / df['n_std_targets']

    fig, axes = plt.subplots(2, 1, figsize=(18, 14))

    # Top: Jaccard similarity
    summary_jac = df.groupby(['kg', 'embedding', 'indication'])['jaccard'].mean().reset_index()
    pivot_jac = summary_jac.pivot_table(
        values='jaccard', index='indication',
        columns=['kg', 'embedding'], aggfunc='mean'
    )
    pivot_jac.columns = [f'{kg}-{emb}' for kg, emb in pivot_jac.columns]
    sns.heatmap(pivot_jac, annot=True, fmt='.2f', cmap='YlOrRd', vmin=0, vmax=1,
                ax=axes[0], annot_kws={'size': 11})
    axes[0].set_title('Mean Jaccard Similarity: Standard vs Confident Negative Targets', fontsize=13)
    axes[0].set_xlabel('KG-Embedding', fontsize=11)
    axes[0].tick_params(axis='x', labelsize=9, rotation=45)
    axes[0].tick_params(axis='y', labelsize=10)

    # Bottom: Recovery rate (fraction of standard targets also in confident)
    summary_rec = df.groupby(['kg', 'embedding', 'indication'])['recovery'].mean().reset_index()
    pivot_rec = summary_rec.pivot_table(
        values='recovery', index='indication',
        columns=['kg', 'embedding'], aggfunc='mean'
    )
    pivot_rec.columns = [f'{kg}-{emb}' for kg, emb in pivot_rec.columns]
    sns.heatmap(pivot_rec, annot=True, fmt='.3f', cmap='YlGn', vmin=0.95, vmax=1,
                ax=axes[1], annot_kws={'size': 11})
    axes[1].set_title('Recovery Rate: Fraction of Standard Targets Also Predicted by Confident Model', fontsize=13)
    axes[1].set_xlabel('KG-Embedding', fontsize=11)
    axes[1].tick_params(axis='x', labelsize=9, rotation=45)
    axes[1].tick_params(axis='y', labelsize=10)

    plt.tight_layout()
    plt.savefig(output_dir / 'jaccard_comparison.png', dpi=300)
    plt.close()


def plot_target_counts(df, output_dir):
    """Scatter: standard vs confident target counts."""
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(df['n_std_targets'], df['n_conf_targets'], alpha=0.3, s=10)
    max_val = max(df['n_std_targets'].max(), df['n_conf_targets'].max())
    ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.5)
    ax.set_xlabel('Standard Negative Targets')
    ax.set_ylabel('Confident Negative Targets')
    ax.set_title('Number of Predicted Targets: Standard vs Confident')
    plt.tight_layout()
    plt.savefig(output_dir / 'target_count_comparison.png', dpi=300)
    plt.close()


def plot_prob_correlation(df, output_dir):
    """Heatmap of probability correlation."""
    summary = df.groupby(['kg', 'embedding', 'indication'])['prob_correlation'].mean().reset_index()
    pivot = summary.pivot_table(
        values='prob_correlation', index='indication',
        columns=['kg', 'embedding'], aggfunc='mean'
    )
    fig, ax = plt.subplots(figsize=(14, 6))
    sns.heatmap(pivot, annot=True, fmt='.3f', cmap='YlOrRd', vmin=0, vmax=1, ax=ax)
    ax.set_title('Mean Probability Correlation: Standard vs Confident Models')
    plt.tight_layout()
    plt.savefig(output_dir / 'prob_correlation.png', dpi=300)
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

    print("=== Confident Negative Selection Comparison ===")
    all_rows = []

    for kg in KGS:
        for embedding in EMBEDDINGS:
            print(f"  {kg}/{embedding}...")
            kg_features = load_kg_features(kg, embedding)
            clinical_data, known_targets = load_clinical_data(scores_date, kg_features)

            rows = run_comparison(kg_features, clinical_data, known_targets, kg, embedding)
            all_rows.extend(rows)

            del kg_features, clinical_data, known_targets
            gc.collect()

    df = pd.DataFrame(all_rows)
    df.to_csv(data_dir / 'confident_neg_comparison.csv', index=False)

    # Summary
    summary = df.groupby(['kg', 'embedding', 'indication']).agg({
        'jaccard': ['mean', 'std'],
        'prob_correlation': ['mean'],
        'n_std_targets': ['mean'],
        'n_conf_targets': ['mean'],
    }).reset_index()
    summary.columns = ['_'.join(c).rstrip('_') for c in summary.columns]
    summary.to_csv(data_dir / 'confident_neg_summary.csv', index=False)

    print(f"\n  Mean Jaccard: {df['jaccard'].mean():.3f}")
    print(f"  Min Jaccard: {df['jaccard'].min():.3f}")
    print(f"  Mean prob correlation: {df['prob_correlation'].mean():.3f}")

    plot_jaccard(df, plots_dir)
    plot_target_counts(df, plots_dir)
    plot_prob_correlation(df, plots_dir)

    print("\n✓ Confident negative analysis complete.")


if __name__ == '__main__':
    main()
