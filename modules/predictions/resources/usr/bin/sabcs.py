#!/usr/bin/env python
"""
SABCS Target Overlap Analysis and Visualization

Aggregates SABCS breast cancer target overlap results across knowledge graphs,
clinical trial filters, RF thresholds, and pathway gene settings. Generates
heatmaps showing overlap percentages by indication and parameters.

Usage:
    sabcs.py <file1.csv> <file2.csv> ... <fileN.csv>

Arguments:
    files: CSV files with SABCS overlap results (format: <kg>_<ct>_<threshold>_<pg>_<iter>.csv)

Output:
    - data/all.pickle: Raw nested dictionary
    - data/all_mean.pickle: Iteration-averaged results
    - plots/<threshold>.png: Heatmap grids for each threshold
"""

import sys
import pickle
from functools import reduce
from pathlib import Path
from typing import Dict, List, Tuple
from collections import defaultdict

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# ─── Configuration ───────────────────────────────────────────────────────────

# Number of iterations to average over
N_ITERATIONS = 10

# RF thresholds to generate heatmaps for
RF_THRESHOLDS = ['0.5', '0.6', '0.7', '0.8', '0.9']

# Pathway gene numbers to include
PG_NUMBERS = ['0', '100', '300']

# Knowledge graphs in grid
KGS = ['biokg', 'hetionet', 'openbiolink', 'primekg']

# Clinical trial filters in grid
CTS = ['All', 'Approved', 'Unique']

# Heatmap color settings
HEATMAP_VMIN = 0
HEATMAP_VMAX = 100
HEATMAP_CENTER = 50
HEATMAP_CMAP = 'vlag'

# Figure dimensions
HEATMAP_FIGSIZE = (9, 8)


# ─── Helper Functions ────────────────────────────────────────────────────────

def parse_filename(filepath: str) -> Tuple[str, str, str, str, str]:
    """
    Extract parameters from SABCS result filename.
    
    Parameters
    ----------
    filepath : str
        Path to SABCS CSV (format: <kg>_<ct>_<threshold>_<pg>_<iter>.csv)
    
    Returns
    -------
    Tuple[str, str, str, str, str]
        Knowledge graph, clinical trial filter, RF threshold, pathway genes, iteration
    """
    filename = Path(filepath).name
    parts = filename.split('_')
    
    kg = parts[0]
    ct = parts[1]
    rf_threshold = parts[2][:3]  # Extract first 3 chars (e.g., '0.5')
    pg_number = parts[3]
    iteration = parts[4].split('.')[0]
    
    return kg, ct, rf_threshold, pg_number, iteration


def load_results(files: List[str]) -> Dict:
    """
    Load all SABCS overlap result files into nested dictionary.
    
    Parameters
    ----------
    files : List[str]
        List of CSV file paths
    
    Returns
    -------
    Dict
        Nested dictionary: [kg][ct][rf_threshold][pg_number][iteration] -> DataFrame
    """
    results = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(dict))))
    
    print(f"Loading {len(files)} SABCS result files...")
    
    for filepath in files:
        try:
            # Parse filename
            kg, ct, rf_threshold, pg_number, iteration = parse_filename(filepath)
            
            # Load data (convert to percentage)
            data = pd.read_csv(filepath, encoding='utf8', index_col=0) * 100
            
            # Store in nested structure
            results[kg][ct][rf_threshold][pg_number][iteration] = data
            
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


def average_iterations(results: Dict) -> Dict:
    """
    Average results over iterations.
    
    Parameters
    ----------
    results : Dict
        Nested results dictionary
    
    Returns
    -------
    Dict
        Averaged results dictionary
    """
    print("Averaging over iterations...")
    
    for kg in sorted(results.keys()):
        for ct in sorted(results[kg].keys()):
            for rf_threshold in sorted(results[kg][ct].keys()):
                for pg_number in sorted(results[kg][ct][rf_threshold].keys()):
                    # Average over iterations (sort iteration keys for determinism)
                    iteration_dfs = [
                        results[kg][ct][rf_threshold][pg_number][iter_key]
                        for iter_key in sorted(results[kg][ct][rf_threshold][pg_number].keys())
                    ]
                    averaged = reduce(
                        lambda a, b: a.add(b, fill_value=0),
                        iteration_dfs
                    ) / N_ITERATIONS
                    
                    results[kg][ct][rf_threshold][pg_number] = averaged
    
    return results


def create_heatmap_grid(results: Dict, rf_threshold: str) -> None:
    """
    Create heatmap grid for a specific RF threshold.
    
    Parameters
    ----------
    results : Dict
        Averaged results dictionary
    rf_threshold : str
        RF threshold to plot
    """
    fig, axes = plt.subplots(
        ncols=len(CTS),
        nrows=len(KGS),
        sharey=True,
        sharex=True,
        figsize=HEATMAP_FIGSIZE
    )
    
    for row_idx, kg in enumerate(KGS):
        for col_idx, ct in enumerate(CTS):
            ax = axes[row_idx, col_idx]
            
            # Combine pathway gene results
            pg_data = pd.concat(
                [results[kg][ct][rf_threshold][pg] for pg in PG_NUMBERS],
                axis=1,
                keys=PG_NUMBERS
            )
            
            # Clean up multi-level columns if present
            if isinstance(pg_data.columns, pd.MultiIndex):
                pg_data.columns = pg_data.columns.get_level_values(0)
            
            pg_data = pg_data[PG_NUMBERS].astype(int)
            
            # Transpose for heatmap
            heatmap_data = pg_data.T
            
            # Show colorbar only on top row
            show_cbar = (row_idx == 0)
            cbar_kws = dict(location="bottom", fraction=0.15) if show_cbar else None
            
            # Create heatmap
            sns.heatmap(
                heatmap_data,
                annot=True,
                fmt='d',
                cbar=show_cbar,
                cbar_kws=cbar_kws,
                ax=ax,
                vmin=HEATMAP_VMIN,
                vmax=HEATMAP_VMAX,
                center=HEATMAP_CENTER,
                cmap=HEATMAP_CMAP
            )
            
            ax.set_xticklabels(ax.get_xticklabels(), rotation=60)
            
            # Add KG label on right
            if col_idx == len(CTS) - 1:
                ax.set_ylabel(kg, rotation=90, size='large')
                ax.yaxis.set_label_position('right')
            
            # Add CT title on top
            if row_idx == 0:
                ax.set_title(ct)
    
    plt.tight_layout()
    return fig


def create_horizontal_unique_heatmap(results: Dict, rf_threshold: str) -> None:
    """
    Create horizontal heatmap showing only Unique CT filter across all KGs.
    
    Parameters
    ----------
    results : Dict
        Averaged results dictionary
    rf_threshold : str
        RF threshold to plot
    """
    fig, axes = plt.subplots(
        ncols=len(KGS),
        nrows=1,
        sharey=True,
        sharex=True,
        figsize=(12, 3)
    )
    
    ct = 'Unique'
    
    for col_idx, kg in enumerate(KGS):
        ax = axes[col_idx]
        
        # Combine pathway gene results
        pg_data = pd.concat(
            [results[kg][ct][rf_threshold][pg] for pg in PG_NUMBERS],
            axis=1,
            keys=PG_NUMBERS
        )
        
        # Clean up multi-level columns if present
        if isinstance(pg_data.columns, pd.MultiIndex):
            pg_data.columns = pg_data.columns.get_level_values(0)
        
        pg_data = pg_data[PG_NUMBERS].astype(int)
        
        # Transpose for heatmap
        heatmap_data = pg_data.T
        
        # Show colorbar only on last column
        show_cbar = (col_idx == len(KGS) - 1)
        cbar_kws = dict(location="right", fraction=0.15) if show_cbar else None
        
        # Create heatmap
        sns.heatmap(
            heatmap_data,
            annot=True,
            fmt='d',
            cbar=show_cbar,
            cbar_kws=cbar_kws,
            ax=ax,
            vmin=HEATMAP_VMIN,
            vmax=HEATMAP_VMAX,
            center=HEATMAP_CENTER,
            cmap=HEATMAP_CMAP
        )
        
        ax.set_xticklabels(ax.get_xticklabels(), rotation=60)
        ax.set_title(kg, size='large')
        
        # Only show y-axis label on first subplot
        if col_idx > 0:
            ax.set_ylabel('')
    
    plt.tight_layout()
    return fig


# ─── Main Execution ──────────────────────────────────────────────────────────

def main():
    """Combine SABCS overlap results and generate visualizations."""
    # Parse command-line arguments
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    
    files = sys.argv[1:]
    
    # Load results
    results = load_results(files)
    
    # Save raw results
    data_dir = Path('data')
    data_dir.mkdir(exist_ok=True)
    
    with open(data_dir / 'all.pickle', 'wb') as f:
        pickle.dump(results, f)
    print("✓ Saved raw results to data/all.pickle")
    
    # Average over iterations
    results_averaged = average_iterations(results)
    
    # Save averaged results
    with open(data_dir / 'all_mean.pickle', 'wb') as f:
        pickle.dump(results_averaged, f)
    print("✓ Saved averaged results to data/all_mean.pickle")
    
    # Create output directory for plots
    plots_dir = Path('plots')
    plots_dir.mkdir(exist_ok=True)
    
    # Generate heatmaps for each RF threshold
    print("Generating heatmaps...")
    for rf_threshold in RF_THRESHOLDS:
        fig = create_heatmap_grid(results_averaged, rf_threshold)
        
        output_path = plots_dir / f'{rf_threshold}.png'
        fig.savefig(output_path)
        plt.close(fig)
        
        print(f"  ✓ {rf_threshold}.png")
    
    # Generate special horizontal plot for RF threshold 0.7 with Unique only
    print("Generating horizontal Unique plot for RF 0.7...")
    fig_horizontal = create_horizontal_unique_heatmap(results_averaged, '0.7')
    output_path_horizontal = plots_dir / '0.7_unique_horizontal.png'
    fig_horizontal.savefig(output_path_horizontal, dpi=300, bbox_inches='tight')
    plt.close(fig_horizontal)
    print(f"  ✓ 0.7_unique_horizontal.png")
    
    print("✓ All visualizations generated")


if __name__ == '__main__':
    main()
