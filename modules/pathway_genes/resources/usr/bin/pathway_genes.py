#!/usr/bin/env python
"""
Pathway Genes Extraction for Disease Indications

Queries pathway databases to extract disease-associated genes for multiple
indications. Uses the Reactome pathway database to identify genes involved
in biological pathways relevant to each disease.

Usage:
    pathway_genes.py

Output:
    - pathway_genes/<YYYY-MM-DD>/breast.csv
    - pathway_genes/<YYYY-MM-DD>/lung.csv
    - pathway_genes/<YYYY-MM-DD>/bowel.csv
    - pathway_genes/<YYYY-MM-DD>/prostate.csv
    - pathway_genes/<YYYY-MM-DD>/melanoma.csv
    - pathway_genes/<YYYY-MM-DD>/diabetes.csv
    - pathway_genes/<YYYY-MM-DD>/cardiovascular.csv
"""

from datetime import datetime
from pathlib import Path
from typing import Dict

import pandas as pd

from alethiotx.artemis.pathway.genes import get


# ─── Configuration ───────────────────────────────────────────────────────────

# Map disease search terms to output filenames
SEARCH_TO_FILENAME = {
    'Breast Cancer': 'breast.csv',
    'Lung Cancer': 'lung.csv',
    'Bowel Cancer': 'bowel.csv',
    'Prostate Cancer': 'prostate.csv',
    'Melanoma': 'melanoma.csv',
    'Diabetes Mellitus Type 2': 'diabetes.csv',
    'Cardiovascular Disease': 'cardiovascular.csv'
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


def extract_pathway_genes(
    search_term: str,
    output_path: Path
) -> int:
    """
    Extract pathway genes for a disease and save to CSV.
    
    Parameters
    ----------
    search_term : str
        Disease search term for pathway query
    output_path : Path
        Path to save CSV file
    
    Returns
    -------
    int
        Number of genes extracted
    """
    genes_df = get(search=search_term)
    genes_df.to_csv(output_path)
    return len(genes_df)


# ─── Main Execution ──────────────────────────────────────────────────────────

def main():
    """Extract pathway genes for all disease indications."""
    # Get current date for output directory
    date = get_current_date()
    print(f"Pathway Genes Extraction - {date}")
    print("=" * 80)
    
    # Create output directory
    output_dir = Path('pathway_genes') / date
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_dir}\n")
    
    # Extract pathway genes for each disease
    print(f"Extracting pathway genes for {len(SEARCH_TO_FILENAME)} diseases...")
    print()
    
    total_genes = 0
    
    for search_term, filename in SEARCH_TO_FILENAME.items():
        output_path = output_dir / filename
        
        try:
            num_genes = extract_pathway_genes(search_term, output_path)
            total_genes += num_genes
            print(f"  ✓ {filename}: {num_genes} genes")
        except Exception as e:
            print(f"  ✗ {filename}: {e}")
    
    print()
    print("=" * 80)
    print(f"✓ Extracted {total_genes} total pathway genes")
    print("=" * 80)


if __name__ == '__main__':
    main()
