#!/usr/bin/env python
"""
Binary Classification Cross-Validation

Evaluates binary classification performance of Random Forest and SVM classifiers
on knowledge graph features for drug target prediction.

Usage:
    binary.py <kg> <indication> <date>

Arguments:
    kg: Knowledge graph name (e.g., 'hetionet', 'biokg')
    indication: Disease indication (e.g., 'breast', 'diabetes')
    date: Clinical scores date (YYYY-MM-DD format)

Output:
    CSV file with cross-validation scores for AUROC and accuracy metrics
"""

import sys
from pathlib import Path
from typing import List, Dict, Any

import pandas as pd
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier

from alethiotx.artemis.cv.pipeline import run


# ─── Configuration ───────────────────────────────────────────────────────────

# Fixed random seed for reproducibility
RANDOM_STATE = 42

# Classifier definitions
CLASSIFIERS = {
    'Random Forest': RandomForestClassifier(random_state=RANDOM_STATE),
    'SVM': SVC(kernel='linear', C=1, random_state=RANDOM_STATE)
}

# Evaluation configurations: (metric, shuffle_labels)
EVALUATIONS = [
    ('roc_auc', False, 'Real'),    # Real targets, AUROC
    ('roc_auc', True, 'Random'),   # Random labels, AUROC baseline
    ('accuracy', False, 'Real'),   # Real targets, Accuracy
    ('accuracy', True, 'Random')   # Random labels, Accuracy baseline
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


def evaluate_classifier(
    kg_features: pd.DataFrame,
    clinical_scores: pd.DataFrame,
    classifier_name: str,
    classifier_obj: Any,
    scoring: str,
    shuffle: bool,
    targets: str
) -> Dict[str, Any]:
    """
    Run cross-validation for a single classifier configuration.
    
    Parameters
    ----------
    kg_features : pd.DataFrame
        Feature matrix from knowledge graph
    clinical_scores : pd.DataFrame
        Clinical trial scores
    classifier_name : str
        Human-readable classifier name
    classifier_obj : sklearn estimator or None
        Classifier instance (None for default RF)
    scoring : str
        Metric to evaluate ('roc_auc' or 'accuracy')
    shuffle : bool
        Whether to shuffle target labels (for baseline)
    targets : str
        Target type label ('Real' or 'Random')
    
    Returns
    -------
    Dict[str, Any]
        Dictionary with scores and metadata
    """
    scores = run(
        kg_features,
        clinical_scores,
        y_slot='y_binary',
        scoring=scoring,
        classifier=classifier_obj,
        shuffle_scores=shuffle
    )

    # Create one result row per score value
    return [
        {
            'score': float(score),
            'scoring': scoring,
            'classifier': classifier_name,
            'targets': targets,
            'bins': ''
        }
        for score in scores
    ]


# ─── Main Execution ──────────────────────────────────────────────────────────

def main():
    """Execute binary classification cross-validation."""
    # Parse command-line arguments
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    
    kg, indication, date = sys.argv[1:4]
    
    # Load data
    print(f"Loading data for {kg} / {indication} ({date})...")
    kg_features, clinical_scores = load_data(kg, indication, date)
    
    # Run all evaluation combinations
    results: List[Dict[str, Any]] = []
    
    for scoring, shuffle, targets in EVALUATIONS:
        for classifier_name, classifier_obj in CLASSIFIERS.items():
            print(f"  {targets} {classifier_name} ({scoring})...")
            
            result_list = evaluate_classifier(
                kg_features,
                clinical_scores,
                classifier_name,
                classifier_obj,
                scoring,
                shuffle,
                targets
            )
            results.extend(result_list)
    
    # Compile results
    df_results = pd.DataFrame(results)
    df_results['bins'] = 'binary'
    
    # Save output
    output_dir = Path('data')
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / f'{kg}_{indication}_binary.csv'
    df_results.to_csv(output_path, index=False)
    print(f"✓ Results saved to {output_path}")


if __name__ == '__main__':
    main()
