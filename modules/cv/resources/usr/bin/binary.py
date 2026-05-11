#!/usr/bin/env python
"""
Binary Classification Cross-Validation

Evaluates binary classification performance of Random Forest and SVM classifiers
on knowledge graph features for drug target prediction.

Usage:
    binary.py <kg> <embedding> <indication> <date>

Arguments:
    kg: Knowledge graph name (e.g., 'hetionet', 'biokg')
    embedding: Embedding method (e.g., 'ComplEx', 'DistMult', 'RotatE', 'TransE')
    indication: Disease indication (e.g., 'breast', 'diabetes')
    date: Clinical scores date (YYYY-MM-DD format)

Output:
    CSV file with cross-validation scores for AUROC and accuracy metrics
"""

import sys
import tempfile
from pathlib import Path
from typing import List, Dict, Any

import fsspec
import pandas as pd
import torch
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
    ('accuracy', True, 'Random'),  # Random labels, Accuracy baseline
    ('average_precision', False, 'Real'),   # Real targets, AUPRC
    ('average_precision', True, 'Random')   # Random labels, AUPRC baseline
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
        for scoring, shuffle, targets in EVALUATIONS:
            for classifier_name, classifier_obj in CLASSIFIERS.items():
                print(f"  [{feature_type}] {targets} {classifier_name} ({scoring})...")

                result_list = evaluate_classifier(
                    features,
                    clinical_scores,
                    classifier_name,
                    classifier_obj,
                    scoring,
                    shuffle,
                    targets
                )
                for r in result_list:
                    r['feature_type'] = feature_type
                results.extend(result_list)

    # Compile results
    df_results = pd.DataFrame(results)
    df_results['bins'] = 'binary'
    
    # Save output
    output_dir = Path('data')
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / f'{kg}_{embedding}_{indication}_binary.csv'
    df_results.to_csv(output_path, index=False)
    print(f"✓ Results saved to {output_path}")


if __name__ == '__main__':
    main()
