#!/usr/bin/env python
"""
Clinical Scores Computation for Disease Indications

Computes clinical trial phase scores and clinical scores for drug targets
across multiple disease indications. Queries ChEMBL for drug-target-disease
associations and clinical trial data, then generates scored target lists.

Clinical scores combine:
- Phase scores: Based on highest clinical trial phase reached
- Trial recency: Prioritizes trials from recent years
- Target-disease specificity: Weights by indication relevance

Usage:
    clinical_scores.py

Output:
    - clinical_scores/<YYYY-MM-DD>/breast.csv
    - clinical_scores/<YYYY-MM-DD>/lung.csv
    - clinical_scores/<YYYY-MM-DD>/prostate.csv
    - clinical_scores/<YYYY-MM-DD>/melanoma.csv
    - clinical_scores/<YYYY-MM-DD>/bowel.csv
    - clinical_scores/<YYYY-MM-DD>/diabetes.csv
    - clinical_scores/<YYYY-MM-DD>/cardiovascular.csv
"""

from datetime import datetime
from pathlib import Path
from typing import Dict

import pandas as pd

from alethiotx.artemis.clinical.scores import compute


# ─── Configuration ───────────────────────────────────────────────────────────

# ChEMBL database version
CHEMBL_VERSION = '36'

# Only include trials from last N years
TRIALS_RECENT_YEARS = 6

# Disease MeSH terms to process
DISEASES = [
    'Breast Neoplasms',
    'Lung Neoplasms',
    'Prostatic Neoplasms',
    'Skin Neoplasms',
    'Intestinal Neoplasms',
    'Diabetes Mellitus, Type 2',
    'Cardiovascular Diseases'
]

# Mapping from MeSH disease terms to output filenames
DISEASE_TO_FILENAME = {
    'Breast Neoplasms': 'breast.csv',
    'Lung Neoplasms': 'lung.csv',
    'Prostatic Neoplasms': 'prostate.csv',
    'Skin Neoplasms': 'melanoma.csv',
    'Intestinal Neoplasms': 'bowel.csv',
    'Diabetes Mellitus, Type 2': 'diabetes.csv',
    'Cardiovascular Diseases': 'cardiovascular.csv'
}

# Column name standardization
COLUMN_MAPPING = {
    'target_gene_name': 'Target Gene',
    'phase_scores': 'Phase Score',
    'clinical_scores': 'Clinical Score'
}


# ─── Helper Functions ────────────────────────────────────────────────────────

def get_current_date() -> str:
    """
    Get current date in YYYY-MM-DD format.
    
    Returns
    -------
    str
        Current date string
    """
    now = datetime.now()
    return now.strftime('%Y-%m-%d')


def save_disease_results(
    results: Dict[str, pd.DataFrame],
    output_dir: Path
) -> None:
    """
    Save disease results to individual CSV files.
    
    Parameters
    ----------
    results : dict
        Dictionary mapping disease names to result DataFrames
    output_dir : Path
        Output directory for CSV files
    """
    print("\n" + "=" * 80)
    print("Saving results...")
    print("=" * 80)
    
    for disease, filename in DISEASE_TO_FILENAME.items():
        if disease in results:
            df = results[disease].copy()
            
            # Standardize column names
            df.rename(columns=COLUMN_MAPPING, inplace=True)
            
            # Save to CSV
            output_path = output_dir / filename
            df.to_csv(output_path, index=False)
            
            print(f"  ✓ {filename}: {len(df)} targets")
        else:
            print(f"  ✗ {filename}: No results found")
    
    print("=" * 80)


# ─── Main Execution ──────────────────────────────────────────────────────────

def main():
    """Compute and save clinical scores for all disease indications."""
    # Get current date for output directory
    date = get_current_date()
    print(f"Clinical Scores Computation - {date}")
    print("=" * 80)
    
    # Create output directory
    output_dir = Path('clinical_scores') / date
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_dir}")
    
    # Compute clinical scores for all diseases
    print(f"\nComputing clinical scores for {len(DISEASES)} diseases...")
    print(f"ChEMBL version: {CHEMBL_VERSION}")
    print(f"Trial window: Last {TRIALS_RECENT_YEARS} years")
    print()
    
    results = compute(
        DISEASES,
        chembl_version=CHEMBL_VERSION,
        trials_only_last_n_years=TRIALS_RECENT_YEARS
    )
    
    print(f"\n✓ Computed scores for {len(results)} diseases")
    
    # Save results to CSV files
    save_disease_results(results, output_dir)
    
    print("\n✓ All processing complete!")


if __name__ == '__main__':
    main()