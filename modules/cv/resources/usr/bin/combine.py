#!/usr/bin/env python
"""
Cross-Validation Results Aggregation and Visualization

Combines CV results from multiple knowledge graphs and indications,
standardizes labels, and generates comparative visualizations.

Usage:
    combine.py <file1.csv> <file2.csv> ... <fileN.csv>

Arguments:
    files: CSV files with CV results (format: <kg>_<indication>_<type>.csv)

Output:
    - data/cv_all.csv: Combined results
    - plots/all.png: KG comparison across metrics
    - plots/all_indications.png: Full comparison by indication
    - plots/all_baseline.png: Random baseline comparisons
"""

import sys
from pathlib import Path
from typing import List, Tuple

import pandas as pd
from plotnine import (
    ggplot, aes, geom_boxplot, facet_grid, theme_seaborn,
    theme, element_text, xlab, ylab, ylim, ggsave
)


# ─── Configuration ───────────────────────────────────────────────────────────

# Label mappings for standardization
BIN_LABELS = {
    '2': '3 classes',
    '3': '4 classes',
    '5': '6 classes',
    'nan': 'No binning'
}

SCORING_LABELS = {
    'accuracy': 'Accuracy',
    'r2': 'R^2 (regression)',
    'roc_auc': 'AUROC'
}

CLASSIFIER_LABELS = {
    'KNeighbors Regressor': 'KN',
    'Linear Regression': 'LR',
    'Random Forest': 'RF',
    'SVM': 'SVM'
}

# Scoring metric order for plots
SCORING_ORDER = ["R^2 (regression)", "Accuracy", "AUROC"]

# Plot configurations
PLOT_CONFIGS = [
    {
        'name': 'all',
        'filter': lambda df: df['targets'] == 'Real',
        'fill': 'kg',
        'color': 'kg',
        'facet': '~scoring+bins',
        'ylim': (0, 1),
        'width': 16,
        'height': 5
    },
    {
        'name': 'all_indications',
        'filter': lambda df: df['targets'] == 'Real',
        'fill': 'indication',
        'color': 'indication',
        'facet': 'kg~scoring+bins',
        'ylim': (0, 1),
        'width': 16,
        'height': 9
    },
    {
        'name': 'all_baseline',
        'filter': lambda df: df['targets'] == 'Random',
        'fill': 'indication',
        'color': 'indication',
        'facet': 'kg~scoring+bins',
        'ylim': (-1, 1),
        'width': 16,
        'height': 9
    }
]


# ─── Helper Functions ────────────────────────────────────────────────────────

def parse_filename(filepath: str) -> Tuple[str, str]:
    """
    Extract knowledge graph and indication from filename.
    
    Parameters
    ----------
    filepath : str
        Path to CV results file (format: <kg>_<indication>_<type>.csv)
    
    Returns
    -------
    Tuple[str, str]
        Knowledge graph name and indication
    """
    filename = Path(filepath).name
    parts = filename.split('_')
    kg = parts[0]
    indication = parts[1]
    return kg, indication


def load_and_annotate_file(filepath: str) -> pd.DataFrame:
    """
    Load a CV results file and add metadata columns.
    
    Parameters
    ----------
    filepath : str
        Path to CV results CSV file
    
    Returns
    -------
    pd.DataFrame
        DataFrame with added 'kg' and 'indication' columns
    """
    kg, indication = parse_filename(filepath)
    
    df = pd.read_csv(filepath)
    df = df.reindex(sorted(df.columns), axis=1)
    df['kg'] = kg
    df['indication'] = indication
    
    return df


def standardize_labels(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize bin, scoring, and classifier labels.
    
    Parameters
    ----------
    df : pd.DataFrame
        Combined CV results
    
    Returns
    -------
    pd.DataFrame
        DataFrame with standardized labels
    """
    # Standardize bins
    df['bins'] = df['bins'].astype('str')
    df = df[~(df['bins'] == '1')]  # Remove single bin results
    
    for old_label, new_label in BIN_LABELS.items():
        df.loc[df['bins'] == old_label, 'bins'] = new_label
    
    # Standardize scoring metrics
    for old_label, new_label in SCORING_LABELS.items():
        df.loc[df['scoring'] == old_label, 'scoring'] = new_label
    
    # Standardize classifier names
    for old_label, new_label in CLASSIFIER_LABELS.items():
        df.loc[df['classifier'] == old_label, 'classifier'] = new_label
    
    # Set categorical order for scoring (only use categories that exist in the data)
    df['scoring'] = df['scoring'].astype('category')
    existing_scores = df['scoring'].cat.categories
    ordered_scores = [s for s in SCORING_ORDER if s in existing_scores]
    if ordered_scores:  # Only reorder if we have matching categories
        df['scoring'] = df['scoring'].cat.reorder_categories(ordered_scores)
    
    return df


def create_plot(df: pd.DataFrame, config: dict) -> object:
    """
    Create a boxplot visualization based on configuration.
    
    Parameters
    ----------
    df : pd.DataFrame
        Combined CV results
    config : dict
        Plot configuration with filter, aesthetics, and layout
    
    Returns
    -------
    plotnine.ggplot
        Configured plot object
    """
    filtered_df = df[config['filter'](df)]
    
    plot = (
        ggplot(
            aes(x='classifier', y='score', fill=config['fill'], color=config['color']),
            filtered_df
        )
        + geom_boxplot()
        + facet_grid(config['facet'], scales="free_x")
        + theme_seaborn()
        + theme(text=element_text(size=14))
        + xlab('')
        + ylab('Score')
        + ylim(*config['ylim'])
    )
    
    return plot


# ─── Main Execution ──────────────────────────────────────────────────────────

def main():
    """Combine CV results and generate visualizations."""
    # Parse command-line arguments
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    
    files = sys.argv[1:]
    
    # Load and combine all files
    print(f"Loading {len(files)} CV result files...")
    dfs: List[pd.DataFrame] = []
    
    for filepath in files:
        try:
            df = load_and_annotate_file(filepath)
            dfs.append(df)
            print(f"  ✓ {Path(filepath).name}")
        except Exception as e:
            print(f"  ✗ {Path(filepath).name}: {e}")
    
    combined = pd.concat(dfs, ignore_index=True)
    print(f"Combined {len(combined)} results")
    
    # Standardize labels
    print("Standardizing labels...")
    combined = standardize_labels(combined)
    
    # Save combined data
    output_dir = Path('data')
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / 'cv_all.csv'
    combined.to_csv(output_path, index=False)
    print(f"✓ Saved combined data to {output_path}")
    
    # Create output directory for plots
    plot_dir = Path('plots')
    plot_dir.mkdir(exist_ok=True)
    
    # Generate all plots
    print("Generating visualizations...")
    for config in PLOT_CONFIGS:
        plot_path = plot_dir / f"{config['name']}.png"
        
        plot = create_plot(combined, config)
        ggsave(
            plot,
            filename=str(plot_path),
            width=config['width'],
            height=config['height'],
            dpi=300
        )
        
        print(f"  ✓ {plot_path}")
    
    print("✓ All visualizations generated")


if __name__ == '__main__':
    main()
