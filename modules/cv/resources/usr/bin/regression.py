#!/usr/bin/env python
"""
Regression Cross-Validation

Evaluates regression performance of Linear Regression and KNN Regressor
on knowledge graph features for continuous drug target prediction.

Usage:
    regression.py <kg> <embedding> <indication> <date>

Arguments:
    kg: Knowledge graph name (e.g., 'hetionet', 'biokg')
    embedding: Embedding method (e.g., 'ComplEx', 'DistMult', 'RotatE', 'TransE')
    indication: Disease indication (e.g., 'breast', 'diabetes')
    date: Clinical scores date (YYYY-MM-DD format)

Output:
    CSV file with cross-validation R² scores
"""

import sys
import tempfile
from pathlib import Path
from typing import List, Dict, Any

import fsspec
import pandas as pd
import torch
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

def load_data(kg: str, embedding: str, indication: str, date: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load knowledge graph features and clinical scores.
    
    Parameters
    ----------
    kg : str
        Knowledge graph name
    embedding : str
        Embedding method name
    indication : str
        Disease indication
    date : str
        Clinical scores date (YYYY-MM-DD)
    
    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        KG feature matrix and clinical scores DataFrame
    """
    kg_path = f's3://alethiotx-artemis/data/kgs-no-data-leakage/associations/{kg}/{embedding}/summarize/predictions.parquet'
    scores_path = f's3://alethiotx-artemis/data/clinical_scores/{date}/{indication}.csv'
    
    kg_features = pd.read_parquet(kg_path)
    clinical_scores = pd.read_csv(scores_path)
    
    return kg_features, clinical_scores


def load_raw_embeddings(kg: str, embedding: str) -> pd.DataFrame:
    """
    Load raw entity embeddings from trained PyKEEN model and filter to genes.

    Downloads the trained model from S3, extracts entity embedding weights,
    and maps entity indices to gene symbols using the genes hash table.

    Parameters
    ----------
    kg : str
        Knowledge graph name
    embedding : str
        Embedding method name

    Returns
    -------
    pd.DataFrame
        Gene embedding matrix with gene symbols as index
    """
    # Load genes hash table (entity_label -> gene_symbol)
    genes_ht = pd.read_csv(
        f's3://alethiotx-artemis/data/kgs-no-data-leakage/associations/{kg}/{embedding}/prepare/genes_hash_table.csv',
        header=None, names=['entity_label', 'gene_symbol']
    )

    # Load entity-to-id mapping
    entity_to_id = pd.read_csv(
        f's3://alethiotx-artemis/data/kgs-no-data-leakage/embeddings/{kg}/{embedding}/training_triples/entity_to_id.tsv.gz',
        sep='\t', compression='gzip'
    )

    # Download trained model to temp file (torch.load doesn't support S3)
    model_path = f's3://alethiotx-artemis/data/kgs-no-data-leakage/embeddings/{kg}/{embedding}/trained_model.pkl'
    with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as f:
        tmp = f.name
    try:
        fs = fsspec.filesystem('s3')
        fs.get(model_path, tmp)
        model = torch.load(tmp, map_location='cpu', weights_only=False)
        embeddings = model.entity_representations[0]._embeddings.weight.detach().numpy()
    finally:
        Path(tmp).unlink(missing_ok=True)

    # Map entity labels to embedding indices
    label_to_id = dict(zip(entity_to_id['label'], entity_to_id['id']))

    # Build gene embedding DataFrame
    gene_indices = []
    gene_symbols = []
    for _, row in genes_ht.iterrows():
        entity_label = row['entity_label']
        if entity_label in label_to_id:
            gene_indices.append(label_to_id[entity_label])
            gene_symbols.append(row['gene_symbol'])

    gene_embeddings = embeddings[gene_indices]
    return pd.DataFrame(gene_embeddings, index=gene_symbols)


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
    if len(sys.argv) != 5:
        print(__doc__)
        sys.exit(1)
    
    kg, embedding, indication, date = sys.argv[1:5]
    
    # Load data
    print(f"Loading data for {kg} / {embedding} / {indication} ({date})...")
    kg_features, clinical_scores = load_data(kg, embedding, indication, date)

    print(f"Loading raw embeddings for {kg} / {embedding}...")
    raw_embeddings = load_raw_embeddings(kg, embedding)

    # Run all evaluation combinations for both feature types
    feature_sets = [
        ('LP Scores', kg_features),
        ('Raw Embeddings', raw_embeddings)
    ]

    results: List[Dict[str, Any]] = []

    for feature_type, features in feature_sets:
        for shuffle, targets in TARGET_CONFIGS:
            for regressor_name, regressor_obj in REGRESSORS.items():
                print(f"  [{feature_type}] {targets} {regressor_name}...")

                result_list = evaluate_regressor(
                    features,
                    clinical_scores,
                    regressor_name,
                    regressor_obj,
                    shuffle,
                    targets
                )
                for r in result_list:
                    r['feature_type'] = feature_type
                results.extend(result_list)
    
    # Compile results
    df_results = pd.DataFrame(results)
    
    # Save output
    output_dir = Path('data')
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / f'{kg}_{embedding}_{indication}_regression.csv'
    df_results.to_csv(output_path, index=False)
    print(f"✓ Results saved to {output_path}")


if __name__ == '__main__':
    main()