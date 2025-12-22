#!/usr/bin/env python
"""
Drug Target Prediction Results Aggregation and Visualization

Combines target prediction results across knowledge graphs, clinical trial filters,
RF thresholds, and pathway gene settings. Computes sensitivity and specificity,
generates comparative plots, and creates heatmaps for cross-indication analysis.

Usage:
    combine.py <file1.csv> <file2.csv> ... <fileN.csv>

Arguments:
    files: CSV files with prediction results (format: <kg>_<ct>_<threshold>_<pg>_<iter>.csv)

Output:
    - data/all.pickle: Raw nested dictionary
    - data/all_mean.pickle: Iteration-averaged results
    - data/all.csv: Melted metrics for plotting
    - plots/indications/: Specificity and sensitivity by indication
    - plots/kgs/all.png: All KG comparison
    - plots/heatmaps/: Cross-indication heatmaps
"""

import sys
import pickle
from functools import reduce
from pathlib import Path
from typing import Dict, List, Tuple
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from plotnine import (
    ggplot, aes, geom_bar, geom_boxplot, facet_grid,
    theme_seaborn, xlab, ylab, ylim
)


# ─── Configuration ───────────────────────────────────────────────────────────

# Number of iterations to average over
N_ITERATIONS = 10

# Indication order for categorical plotting
INDICATION_ORDER = ['breast', 'lung', 'bowel', 'prostate', 'melanoma', 'diabetes', 'cardiovascular']

# Main text heatmap configuration
MAIN_TEXT_CONFIG = {
    'kg': 'hetionet',
    'ct': 'All',
    'rf_threshold': '0.5',
    'pg_number': '0'
}

# Second subfigure configuration
MAIN_TEXT_CONFIG_2 = {
    'kg': 'hetionet',
    'ct': 'Unique',
    'rf_threshold': '0.7',
    'pg_number': '0'
}

# Full heatmap grid parameters
HEATMAP_RF_THRESHOLDS = ['0.5', '0.6', '0.7', '0.8', '0.9']
HEATMAP_PG_NUMBERS = ['0', '100', '300']
HEATMAP_KGS = ['biokg', 'hetionet', 'openbiolink', 'primekg']
HEATMAP_CTS = ['All', 'Approved', 'Unique']

# Heatmap color settings
HEATMAP_VMIN = 0
HEATMAP_VMAX = 100
HEATMAP_CENTER = 50
HEATMAP_CMAP = 'vlag'


# ─── Helper Functions ────────────────────────────────────────────────────────

def parse_filename(filepath: str) -> Tuple[str, str, str, str, str]:
    """
    Extract parameters from prediction result filename.
    
    Parameters
    ----------
    filepath : str
        Path to prediction CSV (format: <kg>_<ct>_<threshold>_<pg>_<iter>.csv)
    
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
    Load all prediction result files into nested dictionary.
    
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
    
    print(f"Loading {len(files)} result files...")
    
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

def compute_metrics(results: Dict) -> Tuple[Dict, pd.DataFrame]:
    """
    Average over iterations and compute sensitivity/specificity metrics.
    
    Parameters
    ----------
    results : Dict
        Nested results dictionary
    
    Returns
    -------
    Tuple[Dict, pd.DataFrame]
        Averaged results and metrics dataframe
    """
    print("Computing metrics...")
    metrics_list = []
    
    for kg in sorted(results.keys()):
        for ct in sorted(results[kg].keys()):
            for rf_threshold in sorted(results[kg][ct].keys()):
                for pg_number in sorted(results[kg][ct][rf_threshold].keys()):
                    # Average over iterations (sort for determinism)
                    iteration_dfs = [
                        results[kg][ct][rf_threshold][pg_number][iter_key]
                        for iter_key in sorted(results[kg][ct][rf_threshold][pg_number].keys())
                    ]
                    averaged = reduce(
                        lambda a, b: a.add(b, fill_value=0),
                        iteration_dfs
                    ) / N_ITERATIONS
                    
                    results[kg][ct][rf_threshold][pg_number] = averaged
                    
                    # Compute metrics
                    transposed = averaged.T
                    row_sums = transposed.sum()
                    diagonal = pd.Series(np.diag(transposed), index=row_sums.index)
                    
                    # Specificity: 1 - (off-diagonal sum / possible off-diagonal comparisons)
                    specificity = 1 - (row_sums - diagonal) / (len(diagonal) - 1) / 100
                    
                    # Sensitivity: diagonal / 100 (already percentage)
                    sensitivity = diagonal / 100
                    
                    # Add metadata
                    metadata = pd.Series({
                        'kg': kg,
                        'ct': ct,
                        'rf_prob': float(rf_threshold),
                        'pathway_genes': int(pg_number)
                    })
                    
                    # Combine and label
                    specificity_row = pd.concat([specificity, metadata])
                    specificity_row['name'] = 'specificity'
                    
                    sensitivity_row = pd.concat([sensitivity, metadata])
                    sensitivity_row['name'] = 'sensitivity'
                    
                    metrics_list.append(specificity_row.to_frame().T)
                    metrics_list.append(sensitivity_row.to_frame().T)
    
    metrics_df = pd.concat(metrics_list, ignore_index=True)
    return results, metrics_df

def prepare_plotting_data(metrics_df: pd.DataFrame) -> pd.DataFrame:
    """
    Melt metrics DataFrame for plotting and set categorical order.
    
    Parameters
    ----------
    metrics_df : pd.DataFrame
        Metrics with sensitivity/specificity
    
    Returns
    -------
    pd.DataFrame
        Melted DataFrame ready for plotting
    """
    # Melt to long format
    melted = metrics_df.melt(id_vars=['kg', 'ct', 'rf_prob', 'pathway_genes', 'name'])
    
    # Fix type issue with pandas melt creating object column
    melted['value'] = melted['value'].astype(np.float64)
    
    # Set categorical order for indications
    melted['variable'] = pd.Categorical(
        melted['variable'],
        categories=INDICATION_ORDER
    )
    
    return melted

def create_indication_plots(data: pd.DataFrame, results: Dict) -> None:
    """
    Generate specificity and sensitivity plots by indication for each KG.
    
    Parameters
    ----------
    data : pd.DataFrame
        Melted plotting data
    results : Dict
        Results dictionary with KGs
    """
    indications_dir = Path('plots/indications')
    indications_dir.mkdir(parents=True, exist_ok=True)
    
    print("Generating indication plots...")
    
    for kg in sorted(results.keys()):
        # Specificity plot
        specificity_data = data[(data['kg'] == kg) & (data['name'] == 'specificity')]
        
        plot = (
            ggplot(
                specificity_data,
                aes('rf_prob', 'value', color='variable', fill='variable')
            )
            + ylim(0, 1)
            + geom_bar(stat='identity', position='dodge')
            + facet_grid('pathway_genes ~ ct', scales='free')
            + theme_seaborn()
            + xlab('Random Forest probability threshold')
            + ylab('Specificity')
        )
        
        plot.save(str(indications_dir / f'specificity_{kg}.png'), dpi=300, width=10, height=5)
        print(f"  ✓ specificity_{kg}.png")
        
        # Sensitivity plot
        sensitivity_data = data[(data['kg'] == kg) & (data['name'] == 'sensitivity')]
        
        plot = (
            ggplot(
                sensitivity_data,
                aes('rf_prob', 'value', color='variable', fill='variable')
            )
            + ylim(0, 1)
            + geom_bar(stat='identity', position='dodge')
            + facet_grid('pathway_genes ~ ct', scales='free')
            + theme_seaborn()
            + xlab('Random Forest probability threshold')
            + ylab('Sensitivity')
        )
        
        plot.save(str(indications_dir / f'sensitivity_{kg}.png'), dpi=300, width=10, height=5)
        print(f"  ✓ sensitivity_{kg}.png")


def create_kg_comparison_plot(data: pd.DataFrame) -> None:
    """
    Generate overall KG comparison plot.
    
    Parameters
    ----------
    data : pd.DataFrame
        Melted plotting data
    """
    kgs_dir = Path('plots/kgs')
    kgs_dir.mkdir(parents=True, exist_ok=True)
    
    print("Generating KG comparison plot...")
    
    plot = (
        ggplot(
            data,
            aes('factor(rf_prob)', 'value', color='name', fill='name')
        )
        + geom_boxplot()
        + ylim(0, 1)
        + facet_grid('pathway_genes ~ ct', scales='free')
        + theme_seaborn()
        + xlab('Random Forest probability threshold')
        + ylab('Metrics')
    )
    
    plot.save(str(kgs_dir / 'all.png'), dpi=300, width=10, height=6)
    print("  ✓ all.png")

def create_main_text_heatmap(results: Dict) -> None:
    """
    Generate simplified heatmap for main text.
    
    Parameters
    ----------
    results : Dict
        Averaged results dictionary
    """
    heatmaps_dir = Path('plots/heatmaps')
    heatmaps_dir.mkdir(parents=True, exist_ok=True)
    
    print("Generating main text heatmap...")
    
    # Create figure with 2 subfigures
    fig, axes = plt.subplots(ncols=2, nrows=1, sharey=True, sharex=True, figsize=(8, 4))
    
    # First subfigure (All CT, RF 0.5)
    config1 = MAIN_TEXT_CONFIG
    data1 = results[config1['kg']][config1['ct']][config1['rf_threshold']][config1['pg_number']]
    heatmap_data1 = data1.T.round(0).astype(int)
    
    sns.heatmap(
        heatmap_data1,
        annot=True,
        fmt='d',
        cbar=True,
        ax=axes[0],
        vmin=HEATMAP_VMIN,
        vmax=HEATMAP_VMAX,
        center=HEATMAP_CENTER,
        cmap=HEATMAP_CMAP
    )
    
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=60)
    axes[0].set_title(f"{config1['ct']} Clinical Targets, RF {config1['rf_threshold']}", fontsize=10)
    
    # Second subfigure (Unique CT, RF 0.7)
    config2 = MAIN_TEXT_CONFIG_2
    data2 = results[config2['kg']][config2['ct']][config2['rf_threshold']][config2['pg_number']]
    heatmap_data2 = data2.T.round(0).astype(int)
    
    sns.heatmap(
        heatmap_data2,
        annot=True,
        fmt='d',
        cbar=True,
        ax=axes[1],
        vmin=HEATMAP_VMIN,
        vmax=HEATMAP_VMAX,
        center=HEATMAP_CENTER,
        cmap=HEATMAP_CMAP
    )
    
    axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=60)
    axes[1].set_title(f"{config2['ct']} Clinical Targets, RF {config2['rf_threshold']}", fontsize=10)
    
    plt.tight_layout()
    plt.savefig(heatmaps_dir / 'for_main_text.png')
    plt.close()
    
    print("  ✓ for_main_text.png")


def create_full_heatmaps(results: Dict) -> None:
    """
    Generate comprehensive heatmap grid across all parameters.
    
    Parameters
    ----------
    results : Dict
        Averaged results dictionary
    """
    heatmaps_dir = Path('plots/heatmaps')
    heatmaps_dir.mkdir(parents=True, exist_ok=True)
    
    print("Generating comprehensive heatmaps...")
    
    for rf_threshold in HEATMAP_RF_THRESHOLDS:
        for pg_number in HEATMAP_PG_NUMBERS:
            fig, axes = plt.subplots(
                ncols=len(HEATMAP_CTS),
                nrows=len(HEATMAP_KGS),
                sharey=True,
                sharex=True,
                figsize=(9, 11)
            )
            
            for row_idx, kg in enumerate(HEATMAP_KGS):
                for col_idx, ct in enumerate(HEATMAP_CTS):
                    ax = axes[row_idx, col_idx]
                    
                    # Get data and format for heatmap
                    data = results[kg][ct][rf_threshold][pg_number]
                    heatmap_data = data.T.round(0).astype(int)
                    
                    # Show colorbar only on top row
                    show_cbar = (row_idx == 0)
                    cbar_kws = dict(location="bottom", pad=0.05) if show_cbar else None
                    
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
                    if col_idx == len(HEATMAP_CTS) - 1:
                        ax.set_ylabel(kg, rotation=90, size='large')
                        ax.yaxis.set_label_position('right')
                    
                    # Add CT title on top
                    if row_idx == 0:
                        ax.set_title(ct)
            
            plt.tight_layout()
            output_path = heatmaps_dir / f'{rf_threshold}_{pg_number}.png'
            plt.savefig(output_path)
            plt.close()
            
            print(f"  ✓ {rf_threshold}_{pg_number}.png")


# ─── Main Execution ──────────────────────────────────────────────────────────

def main():
    """Combine prediction results and generate visualizations."""
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
    
    # Compute metrics and average over iterations
    results_averaged, metrics_df = compute_metrics(results)
    
    # Save averaged results
    with open(data_dir / 'all_mean.pickle', 'wb') as f:
        pickle.dump(results_averaged, f)
    print("✓ Saved averaged results to data/all_mean.pickle")
    
    # Prepare data for plotting
    plotting_data = prepare_plotting_data(metrics_df)
    plotting_data.to_csv(data_dir / 'all.csv', index=False)
    print("✓ Saved plotting data to data/all.csv")
    
    # Generate all plots
    create_indication_plots(plotting_data, results_averaged)
    create_kg_comparison_plot(plotting_data)
    create_main_text_heatmap(results_averaged)
    create_full_heatmaps(results_averaged)
    
    print("✓ All visualizations generated")


if __name__ == '__main__':
    main()