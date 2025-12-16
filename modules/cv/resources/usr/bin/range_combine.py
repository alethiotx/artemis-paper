#!/usr/bin/env python
"""
Feature Range Analysis Results Aggregation and Visualization

Combines feature subsampling CV results, computes median trends,
and generates visualizations showing performance vs feature count.

Usage:
    range_combine.py <file1.csv> <file2.csv> ... <fileN.csv>

Arguments:
    files: CSV files with feature range CV results

Output:
    - data/all.csv: Combined results
    - plots/all.png: Feature count vs AUROC (real targets)
    - plots/all_baseline.png: Feature count vs AUROC (random baseline)
"""

import sys
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
from plotnine import (
    ggplot, aes, geom_point, geom_smooth, geom_vline, theme_seaborn,
    xlab, ylab, ylim, scale_x_log10, ggsave
)


# ─── Configuration ───────────────────────────────────────────────────────────

# Polynomial degree for smoothing
POLY_DEGREE = 7

# Reference line for full feature set
FULL_FEATURE_COUNT = 10000

# Log scale breaks for x-axis
LOG_BREAKS = [100, 1000, 10000]
LOG_LABELS = ['100', '1000', '10000']

# Plot configurations
PLOT_CONFIGS = [
    {
        'name': 'all',
        'target': 'Real',
        'ylim': (0.7, 1),
        'vline': FULL_FEATURE_COUNT,
        'title': 'Real Targets'
    },
    {
        'name': 'all_baseline',
        'target': 'Random',
        'ylim': (0, 1),
        'vline': None,
        'title': 'Random Baseline'
    }
]


# ─── Helper Functions ────────────────────────────────────────────────────────

def poly(x: np.ndarray, degree: int = 1) -> pd.DataFrame:
    """
    Generate polynomial features for fitting.
    
    Note: These are non-orthogonal polynomial features, suitable
    for smoothing but not for extrapolated predictions.
    
    Parameters
    ----------
    x : np.ndarray
        Input feature vector
    degree : int
        Polynomial degree (default: 1)
    
    Returns
    -------
    pd.DataFrame
        DataFrame with polynomial feature columns
    """
    poly_features = {}
    
    for i in range(degree + 1):
        if i == 1:
            poly_features['x'] = x
        else:
            poly_features[f'x**{i}'] = np.power(x, i)
    
    return pd.DataFrame(poly_features)


def compute_median_trends(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute median scores by feature count and target type.
    
    Parameters
    ----------
    df : pd.DataFrame
        Combined CV results
    
    Returns
    -------
    pd.DataFrame
        Median scores grouped by num_features and targets
    """
    # Ensure score column is numeric and drop any rows with non-numeric scores
    df_clean = df[['num_features', 'targets', 'score']].copy()
    df_clean['score'] = pd.to_numeric(df_clean['score'], errors='coerce')
    
    # Drop rows with NaN scores
    df_clean = df_clean.dropna(subset=['score'])
    
    grouped = df_clean.groupby(
        ['num_features', 'targets'],
        as_index=False
    )['score'].median()
    
    grouped.columns = ['num_features', 'targets', 'score']
    
    return grouped


def create_feature_plot(
    data: pd.DataFrame,
    trends: pd.DataFrame,
    config: dict
) -> object:
    """
    Create a feature count vs performance plot.
    
    Parameters
    ----------
    data : pd.DataFrame
        Raw CV results
    trends : pd.DataFrame
        Median trend data
    config : dict
        Plot configuration with target, ylim, and vline settings
    
    Returns
    -------
    plotnine.ggplot
        Configured plot object
    """
    target = config['target']
    
    # Filter data for target type
    data_filtered = data[data['targets'] == target]
    trends_filtered = trends[trends['targets'] == target]
    
    # Base plot
    plot = (
        ggplot(aes(x='num_features', y='score'), data_filtered)
        + geom_point()
        + geom_smooth(
            aes(group='targets'),
            trends_filtered,
            method='lm',
            formula=f'y ~ poly(x, {POLY_DEGREE})',
            se=False
        )
        + theme_seaborn()
        + xlab('Number of features')
        + ylab('AUROC')
        + ylim(*config['ylim'])
        + scale_x_log10(breaks=LOG_BREAKS, labels=LOG_LABELS)
    )
    
    # Add vertical line if specified
    if config['vline'] is not None:
        plot = plot + geom_vline(xintercept=config['vline'])
    
    return plot


# ─── Main Execution ──────────────────────────────────────────────────────────

def test_with_existing_csv(csv_path: str = 'all.csv') -> None:
    """
    Test visualization generation with an existing all.csv file.
    
    This bypasses the file loading and combination steps, directly
    loading the combined CSV and generating plots.
    
    Parameters
    ----------
    csv_path : str
        Path to the combined all.csv file (default: 'all.csv')
    """
    csv_file = Path(csv_path)
    
    if not csv_file.exists():
        print(f"Error: File '{csv_path}' not found")
        sys.exit(1)
    
    print(f"Loading existing data from {csv_path}...")
    combined = pd.read_csv(csv_file)
    print(f"Loaded {len(combined)} results")
    
    # Ensure score column is numeric
    print("Cleaning data...")
    combined['score'] = pd.to_numeric(combined['score'], errors='coerce')
    combined = combined.dropna(subset=['score'])
    print(f"After cleaning: {len(combined)} valid results")
    
    # Compute median trends
    print("Computing median trends...")
    trends = compute_median_trends(combined)
    
    # Create output directory for plots
    plot_dir = Path('plots')
    plot_dir.mkdir(exist_ok=True)
    
    # Generate all plots
    print("Generating visualizations...")
    for config in PLOT_CONFIGS:
        plot_path = plot_dir / f"{config['name']}.png"
        
        # Check if there's data for this target type
        target_data = combined[combined['targets'] == config['target']]
        
        if len(target_data) == 0:
            print(f"  ⚠ {plot_path.name}: No data for target type '{config['target']}', skipping")
            continue
        
        if len(target_data) < 3:
            print(f"  ⚠ {plot_path.name}: Insufficient data ({len(target_data)} points), skipping")
            continue
        
        try:
            plot = create_feature_plot(combined, trends, config)
            ggsave(
                plot,
                filename=str(plot_path),
                width=8,
                height=4,
                dpi=300
            )
            print(f"  ✓ {plot_path}")
        except Exception as e:
            print(f"  ✗ {plot_path.name}: {e}")
    
    print("✓ Visualization generation complete")


def main():
    """Combine feature range results and generate visualizations."""
    # Parse command-line arguments
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    
    files = sys.argv[1:]
    
    # Load and combine all files
    print(f"Loading {len(files)} feature range result files...")
    dfs: List[pd.DataFrame] = []
    
    for filepath in files:
        try:
            df = pd.read_csv(filepath)
            dfs.append(df)
            print(f"  ✓ {Path(filepath).name}")
        except Exception as e:
            print(f"  ✗ {Path(filepath).name}: {e}")
    
    combined = pd.concat(dfs, ignore_index=True)
    print(f"Combined {len(combined)} results")
    
    # Ensure score column is numeric
    print("Cleaning data...")
    combined['score'] = pd.to_numeric(combined['score'], errors='coerce')
    combined = combined.dropna(subset=['score'])
    print(f"After cleaning: {len(combined)} valid results")
    
    # Compute median trends
    print("Computing median trends...")
    trends = compute_median_trends(combined)
    
    # Save combined data
    output_dir = Path('data')
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / 'all.csv'
    combined.to_csv(output_path, index=False)
    print(f"✓ Saved combined data to {output_path}")
    
    # Create output directory for plots
    plot_dir = Path('plots')
    plot_dir.mkdir(exist_ok=True)
    
    # Generate all plots
    print("Generating visualizations...")
    for config in PLOT_CONFIGS:
        plot_path = plot_dir / f"{config['name']}.png"
        
        # Check if there's data for this target type
        target_data = combined[combined['targets'] == config['target']]
        
        if len(target_data) == 0:
            print(f"  ⚠ {plot_path.name}: No data for target type '{config['target']}', skipping")
            continue
        
        if len(target_data) < 3:
            print(f"  ⚠ {plot_path.name}: Insufficient data ({len(target_data)} points), skipping")
            continue
        
        try:
            plot = create_feature_plot(combined, trends, config)
            ggsave(
                plot,
                filename=str(plot_path),
                width=8,
                height=4,
                dpi=300
            )
            print(f"  ✓ {plot_path}")
        except Exception as e:
            print(f"  ✗ {plot_path.name}: {e}")
    
    print("✓ Visualization generation complete")


if __name__ == '__main__':
    main()
    # test_with_existing_csv()