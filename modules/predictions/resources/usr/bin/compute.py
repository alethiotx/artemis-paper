#!/usr/bin/env python
"""
Drug Target Prediction

Trains Random Forest classifiers on knowledge graph features to predict
drug targets across multiple disease indications, validates against clinical
trials, and compares predictions with SABCS breast cancer targets.

Usage:
    compute.py <kg> <ct_unique> <rf_threshold> <pg_number> <iteration> <date>

Arguments:
    kg: Knowledge graph name (e.g., 'hetionet', 'biokg')
    ct_unique: Clinical trial filter ('All', 'Approved', or 'Unique')
    rf_threshold: Probability threshold for target classification (0.0-1.0)
    pg_number: Number of pathway genes to include (0 = none)
    iteration: Random seed for reproducibility
    date: Clinical scores date (YYYY-MM-DD format)

Output:
    - training/<indication>_<params>.csv: Training labels
    - targets/<indication>_<params>.csv: Predicted targets with probabilities
    - predictions/<params>.csv: Cross-indication overlap matrix
    - sabcs/<params>.csv: SABCS target overlap by indication
"""

import sys
from pathlib import Path
from typing import List, Dict, Tuple, Optional

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier

from alethiotx.artemis.cv.pipeline import prepare as prepare_model
from alethiotx.artemis.pathway.genes import load as load_pathway_genes
from alethiotx.artemis.pathway.genes import unique as unique_pathway_genes
from alethiotx.artemis.clinical.scores import load as load_clinical_scores
from alethiotx.artemis.clinical.scores import approved as approved_clinical_scores
from alethiotx.artemis.clinical.scores import unique as unique_clinical_scores
from alethiotx.artemis.clinical.scores import all_targets


# ─── Configuration ───────────────────────────────────────────────────────────

# Disease indications in analysis
INDICATIONS = ['breast', 'lung', 'bowel', 'prostate', 'melanoma', 'diabetes', 'cardiovascular']

# SABCS data source
SABCS_DATA_PATH = 's3://alethiotx-artemis/data/sabcs24/breast_cancer_targets.csv'


# ─── Helper Functions ────────────────────────────────────────────────────────

def parse_arguments() -> Dict[str, any]:
    """
    Parse and validate command-line arguments.
    
    Returns
    -------
    Dict[str, any]
        Configuration dictionary with parsed arguments
    """
    if len(sys.argv) != 7:
        print(__doc__)
        sys.exit(1)
    
    return {
        'kg': sys.argv[1],
        'ct_unique': sys.argv[2],
        'rf_threshold': float(sys.argv[3]),
        'pg_number': int(sys.argv[4]),
        'iteration': int(sys.argv[5]),
        'scores_date': sys.argv[6]
    }

def build_filename(config: Dict[str, any]) -> str:
    """
    Build output filename from configuration parameters.
    
    Parameters
    ----------
    config : Dict[str, any]
        Configuration dictionary
    
    Returns
    -------
    str
        Formatted filename base
    """
    return (
        f"{config['kg']}_{config['ct_unique']}_"
        f"{config['rf_threshold']}_{config['pg_number']}_"
        f"{config['iteration']}.csv"
    )


def load_pathway_genes_for_indications(
    scores_date: str,
    pg_number: int
) -> Dict[str, List[str]]:
    """
    Load pathway genes for all indications.
    
    Parameters
    ----------
    scores_date : str
        Clinical scores date
    pg_number : int
        Number of pathway genes (0 = none)
    
    Returns
    -------
    Dict[str, List[str]]
        Pathway genes by indication
    """
    if pg_number == 0:
        return {indication: [] for indication in INDICATIONS}
    
    # Load and make unique
    pathway_genes = load_pathway_genes(date=scores_date, n=pg_number)
    pathway_genes = unique_pathway_genes(list(pathway_genes))
    
    # Map to indications
    return dict(zip(INDICATIONS, pathway_genes))

def compute_overlap(
    known_targets: pd.Series,
    predicted_targets: pd.Index
) -> float:
    """
    Calculate overlap ratio between known and predicted targets.
    
    Parameters
    ----------
    known_targets : pd.Series
        Known clinical targets
    predicted_targets : pd.Index
        Predicted targets
    
    Returns
    -------
    float
        Overlap ratio (0.0-1.0)
    """
    if len(known_targets) == 0:
        return 0.0
    
    intersection = np.intersect1d(known_targets, predicted_targets)
    return len(intersection) / len(known_targets)


def predict_targets(
    kg_features: pd.DataFrame,
    clinical_scores: pd.DataFrame,
    indication: str,
    known_targets_all: List[str],
    pathway_genes: List[str],
    sabcs_targets: List[str],
    clinical_data: Dict[str, pd.DataFrame],
    config: Dict[str, any],
    filename_base: str
) -> Tuple[pd.DataFrame, pd.Series, float]:
    """
    Train model and predict drug targets for an indication.
    
    Parameters
    ----------
    kg_features : pd.DataFrame
        Knowledge graph feature matrix
    clinical_scores : pd.DataFrame
        Clinical trial scores for this indication
    indication : str
        Disease indication name
    known_targets_all : List[str]
        All known targets across indications
    pathway_genes : List[str]
        Pathway genes for this indication
    sabcs_targets : List[str]
        SABCS breast cancer targets
    clinical_data : Dict[str, pd.DataFrame]
        Clinical scores for all indications
    config : Dict[str, any]
        Configuration parameters
    filename_base : str
        Base filename for outputs
    
    Returns
    -------
    Tuple[pd.DataFrame, pd.Series, float]
        Predictions, cross-indication overlaps, SABCS overlap
    """
    # Prepare training data
    prepare_kwargs = {
        'known_targets': known_targets_all,
        'rand_seed': config['iteration']
    }
    
    if pathway_genes:
        prepare_kwargs['pathway_genes'] = pathway_genes
    
    training_data = prepare_model(kg_features, clinical_scores, **prepare_kwargs)
    
    # Save training labels
    training_dir = Path('training')
    training_dir.mkdir(exist_ok=True)
    training_data['y_binary'].to_csv(training_dir / f"{indication}_{filename_base}")
    
    # Train classifier
    classifier = RandomForestClassifier(random_state=config['iteration'])
    classifier.fit(training_data['X'], training_data['y_binary'])
    
    # Generate predictions
    probabilities = classifier.predict_proba(kg_features)
    binary_predictions = classifier.predict(kg_features)
    
    # Format predictions
    predictions = pd.DataFrame({
        'prob_neg': probabilities[:, 0],
        'prob_target': probabilities[:, 1],
        'is_target': binary_predictions
    }, index=kg_features.index)
    
    # Save predictions
    targets_dir = Path('targets')
    targets_dir.mkdir(exist_ok=True)
    predictions.to_csv(targets_dir / f"{indication}_{filename_base}")
    
    # Get predicted targets above threshold
    predicted_targets = predictions[
        predictions['prob_target'] >= config['rf_threshold']
    ].index
    
    # Compute overlaps with all indications
    overlaps = {}
    for ind_name in INDICATIONS:
        overlaps[ind_name] = compute_overlap(
            clinical_data[ind_name]['Target Gene'],
            predicted_targets
        )
    
    # Compute SABCS overlap
    sabcs_overlap = compute_overlap(sabcs_targets, predicted_targets)
    
    return predictions, pd.Series(overlaps), sabcs_overlap

def load_and_filter_clinical_data(
    kg_features: pd.DataFrame,
    scores_date: str,
    ct_filter: str
) -> Tuple[Dict[str, pd.DataFrame], List[str]]:
    """
    Load and filter clinical trial scores.
    
    Parameters
    ----------
    kg_features : pd.DataFrame
        Knowledge graph features (for filtering)
    scores_date : str
        Clinical scores date
    ct_filter : str
        Filter type ('All', 'Approved', or 'Unique')
    
    Returns
    -------
    Tuple[Dict[str, pd.DataFrame], List[str]]
        Clinical data by indication and all known targets
    """
    # Load all clinical scores
    clinical_data_raw = load_clinical_scores(date=scores_date)
    
    # Apply filter
    if ct_filter == 'Approved':
        clinical_data_raw = approved_clinical_scores(list(clinical_data_raw))
    elif ct_filter == 'Unique':
        clinical_data_raw = unique_clinical_scores(list(clinical_data_raw))
    
    # Map to indication names and filter by KG features
    clinical_data = {}
    for indication, scores in zip(INDICATIONS, clinical_data_raw):
        filtered_scores = scores[scores['Target Gene'].isin(kg_features.index)]
        clinical_data[indication] = filtered_scores
    
    # Get all known targets across indications
    known_targets = all_targets(list(clinical_data_raw))
    
    return clinical_data, known_targets


def load_sabcs_targets(
    kg_features: pd.DataFrame,
    breast_known_targets: pd.Series
) -> List[str]:
    """
    Load SABCS breast cancer targets excluding known clinical targets.
    
    Parameters
    ----------
    kg_features : pd.DataFrame
        Knowledge graph features (for filtering)
    breast_known_targets : pd.Series
        Known breast cancer clinical targets
    
    Returns
    -------
    List[str]
        SABCS target gene list
    """
    sabcs_df = pd.read_csv(SABCS_DATA_PATH, index_col=0)
    
    # Remove metadata column if present
    if 'Unique count(Title)' in sabcs_df.columns:
        sabcs_df.drop('Unique count(Title)', axis=1, inplace=True)
    
    # Filter to KG features and exclude known breast targets
    sabcs_df = sabcs_df[sabcs_df.index.isin(kg_features.index)]
    sabcs_df = sabcs_df[~sabcs_df.index.isin(breast_known_targets)]
    
    return sabcs_df.index.tolist()


# ─── Main Execution ──────────────────────────────────────────────────────────

def main():
    """Execute drug target prediction pipeline."""
    # Parse configuration
    config = parse_arguments()
    
    print(f"Configuration:")
    print(f"  KG: {config['kg']}")
    print(f"  Clinical trial filter: {config['ct_unique']}")
    print(f"  RF threshold: {config['rf_threshold']}")
    print(f"  Pathway genes: {config['pg_number']}")
    print(f"  Iteration: {config['iteration']}")
    print(f"  Date: {config['scores_date']}")
    
    # Build filename
    filename_base = build_filename(config)
    
    # Load knowledge graph features
    print("Loading knowledge graph features...")
    kg_path = f"s3://alethiotx-artemis/data/kgs/associations/{config['kg']}/summarize/predictions.csv"
    kg_features = pd.read_csv(kg_path, index_col=0)
    
    # Load clinical data
    print("Loading clinical trial data...")
    clinical_data, known_targets_all = load_and_filter_clinical_data(
        kg_features,
        config['scores_date'],
        config['ct_unique']
    )
    
    # Load pathway genes
    print("Loading pathway genes...")
    pathway_genes_map = load_pathway_genes_for_indications(
        config['scores_date'],
        config['pg_number']
    )
    
    # Load SABCS targets
    print("Loading SABCS targets...")
    sabcs_targets = load_sabcs_targets(
        kg_features,
        clinical_data['breast']['Target Gene']
    )
    
    # Run predictions for all indications
    print("Running predictions...")
    overlap_results = []
    sabcs_results = []
    
    for indication in INDICATIONS:
        print(f"  {indication}...")
        
        _, overlaps, sabcs_overlap = predict_targets(
            kg_features,
            clinical_data[indication],
            indication,
            known_targets_all,
            pathway_genes_map[indication],
            sabcs_targets,
            clinical_data,
            config,
            filename_base
        )
        
        overlap_results.append(overlaps)
        sabcs_results.append(sabcs_overlap)
    
    # Compile overlap matrices
    overlap_matrix = pd.DataFrame(overlap_results, index=INDICATIONS)
    sabcs_matrix = pd.DataFrame(sabcs_results, index=INDICATIONS)
    
    # Save results
    predictions_dir = Path('predictions')
    predictions_dir.mkdir(exist_ok=True)
    overlap_matrix.to_csv(predictions_dir / filename_base)
    
    sabcs_dir = Path('sabcs')
    sabcs_dir.mkdir(exist_ok=True)
    sabcs_matrix.to_csv(sabcs_dir / filename_base)
    
    print(f"✓ Overlap matrix saved to predictions/{filename_base}")
    print(f"✓ SABCS overlaps saved to sabcs/{filename_base}")


if __name__ == '__main__':
    main()