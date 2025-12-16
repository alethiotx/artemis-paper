#!/usr/bin/env python
"""
Regression Cross-Validation

Evaluates regression performance of Linear Regression and KNN Regressor
on knowledge graph features for continuous drug target prediction.

Usage:
    regression.py <kg> <indication> <date>

Arguments:
    kg: Knowledge graph name (e.g., 'hetionet', 'biokg')
    indication: Disease indication (e.g., 'breast', 'diabetes')
    date: Clinical scores date (YYYY-MM-DD format)

Output:
    CSV file with cross-validation R² scores
"""

import sys
from pathlib import Path
from typing import List, Dict, Any

import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import KNeighborsRegressor
from sklearn.model_selection import KFold

from alethiotx.artemis.cv.pipeline import run


# ─── Configuration ───────────────────────────────────────────────────────────

# Fixed random seed for reproducibility
RANDOM_STATE = 42

# Cross-validation configuration
CV_SPLITS = 5

# Regressor definitions
REGRESSORS = {
    'Linear Regression': LinearRegression(),
    'KNeighbors Regressor': KNeighborsRegressor(n_neighbors=5, weights='uniform')
}

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


def evaluate_regressor(
    kg_features: pd.DataFrame,
    clinical_scores: pd.DataFrame,
    regressor_name: str,
    regressor_obj: Any,
    shuffle: bool,
    targets: str
) -> List[Dict[str, Any]]:
    """
    Run cross-validation for a single regressor configuration.
    
    Parameters
    ----------
    kg_features : pd.DataFrame
        Feature matrix from knowledge graph
    clinical_scores : pd.DataFrame
        Clinical trial scores
    regressor_name : str
        Human-readable regressor name
    regressor_obj : sklearn estimator
        Regressor instance
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
        kg_features,
        clinical_scores,
        y_slot='y',
        scoring='r2',
        classifier=regressor_obj,
        cv=KFold(n_splits=CV_SPLITS, shuffle=True, random_state=RANDOM_STATE),
        shuffle_scores=shuffle
    )
    
    # Create one result row per score value
    return [
        {
            'score': float(score),
            'scoring': 'r2',
            'classifier': regressor_name,
            'targets': targets,
            'bins': ''
        }
        for score in scores
    ]


# ─── Main Execution ──────────────────────────────────────────────────────────

def main():
    """Execute regression cross-validation."""
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
    
    for shuffle, targets in TARGET_CONFIGS:
        for regressor_name, regressor_obj in REGRESSORS.items():
            print(f"  {targets} {regressor_name}...")
            
            result_list = evaluate_regressor(
                kg_features,
                clinical_scores,
                regressor_name,
                regressor_obj,
                shuffle,
                targets
            )
            results.extend(result_list)
    
    # Compile results
    df_results = pd.DataFrame(results)
    
    # Save output
    output_dir = Path('data')
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / f'{kg}_{indication}_regression.csv'
    df_results.to_csv(output_path, index=False)
    print(f"✓ Results saved to {output_path}")


if __name__ == '__main__':
    main()