#!/usr/bin/env python
"""
Feature Range Analysis Cross-Validation

Evaluates Random Forest binary classification performance across different
feature subsamples to analyze feature importance and model stability.

Usage:
    range.py <subsample> <date>

Arguments:
    subsample: Number of features to randomly sample
    date: Clinical scores date (YYYY-MM-DD format)

Output:
    CSV file with cross-validation scores for each subsample iteration
"""

import sys
from pathlib import Path
from typing import List, Dict, Any

import pandas as pd
from sklearn.ensemble import RandomForestClassifier

from alethiotx.artemis.cv.pipeline import run


# ─── Configuration ───────────────────────────────────────────────────────────

# Fixed random seed for reproducibility
RANDOM_STATE = 42

# Fixed analysis parameters
KNOWLEDGE_GRAPH = 'hetionet_full'
INDICATION = 'breast'

# Number of subsampling iterations
N_ITERATIONS = 20

# Classifier definition
CLASSIFIER = RandomForestClassifier(random_state=RANDOM_STATE)

# Target configurations: (shuffle_labels, target_label)
TARGET_CONFIGS = [
    (False, 'Real'),    # Real targets
    (True, 'Random')    # Random labels baseline
]


# ─── Helper Functions ────────────────────────────────────────────────────────

def load_data(kg: str, indication: str, date: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load knowledge graph features and clinical scores.
    
    Parameters
    ----------
    kg : str
        Knowledge graph name
    indication : str
        Disease indication
    date : str
        Clinical scores date (YYYY-MM-DD)
    
    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        KG feature matrix and clinical scores DataFrame
    """
    kg_path = f's3://alethiotx-artemis/data/kgs/associations/{kg}/summarize/predictions.csv'
    scores_path = f's3://alethiotx-artemis/data/clinical_scores/{date}/{indication}.csv'
    
    kg_features = pd.read_csv(kg_path, index_col=0)
    clinical_scores = pd.read_csv(scores_path)
    
    return kg_features, clinical_scores


def evaluate_subsample(
    kg_features_subsample: pd.DataFrame,
    clinical_scores: pd.DataFrame,
    n_features: int,
    shuffle: bool,
    targets: str
) -> List[Dict[str, Any]]:
    """
    Run cross-validation on a feature subsample.
    
    Parameters
    ----------
    kg_features_subsample : pd.DataFrame
        Subsampled feature matrix from knowledge graph
    clinical_scores : pd.DataFrame
        Clinical trial scores
    n_features : int
        Number of features in subsample
    shuffle : bool
        Whether to shuffle target labels (for baseline)
    targets : str
        Target type label ('Real' or 'Random')
    
    Returns
    -------
    List[Dict[str, Any]]
        List of dictionaries, one per score value
    """
    scores = run(
        kg_features_subsample,
        clinical_scores,
        y_slot='y_binary',
        scoring='roc_auc',
        classifier=CLASSIFIER,
        shuffle_scores=shuffle,
        n_iterations=1
    )
    
    # Create one result row per score value
    return [
        {
            'score': float(score),
            'scoring': 'roc_auc',
            'classifier': 'Random Forest',
            'targets': targets,
            'num_features': n_features,
            'bins': 'binary'
        }
        for score in scores
    ]


# ─── Main Execution ──────────────────────────────────────────────────────────

def main():
    """Execute feature range analysis cross-validation."""
    # Parse command-line arguments
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    
    n_features = int(sys.argv[1])
    date = sys.argv[2]
    
    # Load data
    print(f"Loading data for {KNOWLEDGE_GRAPH} / {INDICATION} ({date})...")
    kg_features, clinical_scores = load_data(KNOWLEDGE_GRAPH, INDICATION, date)
    
    print(f"Running {N_ITERATIONS} iterations with {n_features} features...")
    
    # Run all subsample iterations
    results: List[Dict[str, Any]] = []
    
    for iteration in range(N_ITERATIONS):
        # Subsample features for this iteration
        kg_subsample = kg_features.sample(
            n=n_features,
            axis=1,
            random_state=iteration
        )
        
        print(f"  Iteration {iteration + 1}/{N_ITERATIONS}...")
        
        for shuffle, targets in TARGET_CONFIGS:
            result_list = evaluate_subsample(
                kg_subsample,
                clinical_scores,
                n_features,
                shuffle,
                targets
            )
            results.extend(result_list)
    
    # Compile results
    df_results = pd.DataFrame(results)
    
    # Save output
    output_dir = Path('data')
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / f'{n_features}.csv'
    df_results.to_csv(output_path, index=False)
    print(f"✓ Results saved to {output_path}")


if __name__ == '__main__':
    main()
