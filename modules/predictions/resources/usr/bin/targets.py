#!/usr/bin/env python
"""
Target Prediction Probability Aggregation

Aggregates target prediction probabilities across all indications, knowledge graphs,
clinical trial filters, RF thresholds, pathway genes, and iterations.
Creates a single pickled nested dictionary for downstream analysis.

Usage:
    targets.py <file1.csv> <file2.csv> ... <fileN.csv>

Arguments:
    files: CSV files with target predictions (format: <ind>_<kg>_<ct>_<threshold>_<pg>_<iter>.csv)

Output:
    - all_targets.pickle: Nested dictionary of target probabilities by indication
"""

import sys
import pickle
from pathlib import Path
from typing import Dict, List, Tuple
from collections import defaultdict

import pandas as pd


# ─── Helper Functions ────────────────────────────────────────────────────────

def parse_filename(filepath: str) -> Tuple[str, str, str, str, str, str]:
    """
    Extract parameters from target prediction filename.
    
    Parameters
    ----------
    filepath : str
        Path to target CSV (format: <ind>_<kg>_<ct>_<threshold>_<pg>_<iter>.csv)
    
    Returns
    -------
    Tuple[str, str, str, str, str, str]
        Indication, KG, clinical trial filter, RF threshold, pathway genes, iteration
    """
    filename = Path(filepath).name
    parts = filename.split('_')
    
    indication = parts[0]
    kg = parts[1]
    ct = parts[2]
    rf_threshold = parts[3][:3]  # Extract first 3 chars (e.g., '0.5')
    pg_number = parts[4]
    iteration = parts[5].split('.')[0]
    
    return indication, kg, ct, rf_threshold, pg_number, iteration


def load_target_probabilities(files: List[str]) -> Dict:
    """
    Load all target prediction files and organize by indication and parameters.
    
    Parameters
    ----------
    files : List[str]
        List of CSV file paths
    
    Returns
    -------
    Dict
        Nested dictionary: [indication][ct][rf_threshold][pg_number][kg][iteration] -> Series
    """
    results = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))))
    
    print(f"Loading {len(files)} target prediction files...")
    
    for filepath in files:
        try:
            # Parse filename
            indication, kg, ct, rf_threshold, pg_number, iteration = parse_filename(filepath)
            
            # Load target probabilities only
            data = pd.read_csv(filepath, encoding='utf8', index_col=0)['prob_target']
            
            # Store in nested structure with indication as top-level key
            results[indication][ct][rf_threshold][pg_number][kg][iteration] = data
            
            print(f"  ✓ {Path(filepath).name}")
        except Exception as e:
            print(f"  ✗ {Path(filepath).name}: {e}")
    
    # Convert to regular dict recursively for pickling (removes lambdas)
    def convert_to_dict(d):
        """Recursively convert defaultdict to regular dict."""
        if isinstance(d, defaultdict):
            d = {k: convert_to_dict(v) for k, v in d.items()}
        return d
    
    return convert_to_dict(results)


# ─── Main Execution ──────────────────────────────────────────────────────────

def main():
    """Aggregate target prediction probabilities for all indications."""
    # Parse command-line arguments
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    
    files = sys.argv[1:]
    
    print(f"Processing {len(files)} target prediction files...")
    
    # Load and organize predictions
    results = load_target_probabilities(files)
    
    # Print summary
    indications = sorted(results.keys())
    print(f"\nProcessed {len(indications)} indications: {', '.join(indications)}")
    
    # Save to single pickle file
    output_path = Path('all_targets.pickle')
    
    with open(output_path, 'wb') as f:
        pickle.dump(results, f)
    
    print(f"✓ Saved aggregated predictions to {output_path}")


if __name__ == '__main__':
    main()
