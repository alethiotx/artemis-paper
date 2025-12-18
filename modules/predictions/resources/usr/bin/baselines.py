#!/usr/bin/env python
"""
Baseline Prediction Analysis

Computes baseline statistics for target predictions by calculating the percentage
of targets predicted above threshold across iterations. Generates heatmaps showing
baseline percentages by KG, CT filter, and pathway genes, plus mean baseline trends
for each indication.

Usage:
    baselines.py <pickle_file>

Arguments:
    pickle_file: Path to pickled prediction results dictionary (all indications)

Output:
    - data/<indication>.pickle: Baseline statistics per indication
    - plots/<indication>.png: Heatmap grid of baseline percentages per indication
    - plots/<indication>_mean.png: Mean baseline trends plot per indication
"""

import sys
import pickle
import copy
from pathlib import Path
from typing import Dict, List

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from plotnine import (
    ggplot, aes, theme_seaborn, theme, element_text,
    xlab, ylab, geom_boxplot
)


# ─── Configuration ───────────────────────────────────────────────────────────

# Knowledge graphs to analyze
KNOWLEDGE_GRAPHS = ['biokg', 'hetionet', 'openbiolink', 'primekg']

# RF thresholds for analysis
RF_THRESHOLDS = ['0.5', '0.6', '0.7', '0.8', '0.9']

# Clinical trial filters for heatmaps
CT_FILTERS = ['Approved', 'Unique']

# Pathway gene numbers
PG_NUMBERS = ['0', '100', '300']

# Heatmap color settings
HEATMAP_VMIN = 0
HEATMAP_VMAX = 100
HEATMAP_CMAP = 'vlag'

# Figure dimensions
HEATMAP_FIGSIZE = (6, 10)


# ─── Helper Functions ────────────────────────────────────────────────────────

def compute_baseline_statistics(results: Dict) -> Dict:
    """
    Compute baseline statistics for each parameter combination.
    
    Calculates the percentage of targets predicted above threshold across
    iterations, then averages across targets.
    
    Parameters
    ----------
    results : Dict
        Nested results dictionary with predictions
    
    Returns
    -------
    Dict
        Baseline statistics dictionary
    """
    print("Computing baseline statistics...")
    baseline = copy.deepcopy(results)
    
    for ct in results.keys():
        for rf_threshold in results[ct].keys():
            for pg_number in results[ct][rf_threshold].keys():
                kg_baselines = {}
                
                for kg in KNOWLEDGE_GRAPHS:
                    # Concatenate all iterations
                    iterations_df = pd.concat(
                        results[ct][rf_threshold][pg_number][kg].values(),
                        axis=1
                    )
                    
                    # Calculate percentage of iterations above threshold per target
                    threshold_value = float(rf_threshold)
                    above_threshold = (iterations_df >= threshold_value).sum(axis=1)
                    percentage_per_target = (above_threshold / iterations_df.shape[1]) * 100
                    
                    # Mean across all targets
                    kg_baselines[kg] = percentage_per_target.mean()
                
                # Store as Series
                baseline[ct][rf_threshold][pg_number] = pd.Series(kg_baselines)
    
    return baseline

def create_heatmap_grid(baseline: Dict, indication: str) -> None:
    """
    Create heatmap grid showing baseline percentages.
    
    Parameters
    ----------
    baseline : Dict
        Baseline statistics dictionary
    indication : str
        Disease indication name
    """
    plots_dir = Path('plots')
    plots_dir.mkdir(exist_ok=True)
    
    print("Generating heatmap grid...")
    
    fig, axes = plt.subplots(
        ncols=len(CT_FILTERS),
        nrows=len(RF_THRESHOLDS),
        sharey=True,
        sharex=True,
        figsize=HEATMAP_FIGSIZE
    )
    
    for row_idx, rf_threshold in enumerate(RF_THRESHOLDS):
        for col_idx, ct in enumerate(CT_FILTERS):
            ax = axes[row_idx, col_idx]
            
            # Combine pathway gene results
            pg_data = pd.concat(
                [baseline[ct][rf_threshold][pg] for pg in PG_NUMBERS],
                axis=1
            )
            pg_data.columns = PG_NUMBERS
            heatmap_data = pg_data[PG_NUMBERS].T.astype(int)
            
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
                cmap=HEATMAP_CMAP
            )
            
            ax.set_xticklabels(ax.get_xticklabels(), rotation=60)
            
            # Add RF threshold label on right
            if col_idx == len(CT_FILTERS) - 1:
                ax.set_ylabel(rf_threshold, size='large')
                ax.yaxis.set_label_position('right')
            
            # Add CT title on top
            if row_idx == 0:
                ax.set_title(ct)
    
    plt.tight_layout()
    output_path = plots_dir / f'{indication}.png'
    plt.savefig(output_path)
    plt.close()
    
    print(f"  ✓ {indication}.png")

def create_combined_baseline_plot(all_baselines: Dict[str, Dict]) -> None:
    """
    Create comprehensive barplot showing baseline percentages across all indications.
    
    Averages across knowledge graphs and pathway genes for cleaner visualization.
    
    Parameters
    ----------
    all_baselines : Dict[str, Dict]
        Dictionary mapping indication names to their baseline statistics
    Returns
    -------
    pd.DataFrame
        DataFrame used for plotting
    """
    plots_dir = Path('plots')
    plots_dir.mkdir(exist_ok=True)
    
    print("\nGenerating combined baseline plot across all indications...")
    
    # Collect all data points
    data_records = []
    
    for indication, baseline in all_baselines.items():
        for ct in baseline.keys():
            for rf_threshold in baseline[ct].keys():
                # Average across pathway genes and knowledge graphs
                all_values = []
                for pg_number in baseline[ct][rf_threshold].keys():
                    # baseline[ct][rf_threshold][pg_number] is a Series with kg as index
                    all_values.extend(baseline[ct][rf_threshold][pg_number].values)
                
                avg_baseline = sum(all_values) / len(all_values)
                
                data_records.append({
                    'Indication': indication,
                    'Training_Targets': ct,
                    'RF_Threshold': float(rf_threshold),
                    'Baseline_Pct': avg_baseline
                })
    
    df = pd.DataFrame(data_records)
    
    # Create single boxplot showing distribution across all indications
    # X-axis is RF threshold, boxplots show distribution of 7 indications
    # Color by training target type
    plot = (
        ggplot(df, aes(x='factor(RF_Threshold)', y='Baseline_Pct', fill='Training_Targets', color='Training_Targets'))
        + geom_boxplot()
        + theme_seaborn()
        + xlab('Random Forest Probability Threshold')
        + ylab('Baseline Percentage (%) - Averaged across KGs and PGs')
        + theme(
            text=element_text(size=14),
            axis_text=element_text(size=12),
            axis_title=element_text(size=14),
            legend_text=element_text(size=12),
            legend_title=element_text(size=14),
            figure_size=(10, 6)
        )
    )
    
    output_path = plots_dir / 'all_indications_baseline.png'
    plot.save(str(output_path), dpi=300)
    
    print(f"  ✓ all_indications_baseline.png")

    return(df)


# ─── Main Execution ──────────────────────────────────────────────────────────

def main():
    """Compute baseline statistics and generate visualizations."""
    # Parse command-line arguments
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    
    pickle_file = sys.argv[1]
    
    # Load prediction results (all indications)
    print("Loading prediction results...")
    with open(pickle_file, 'rb') as f:
        all_results = pickle.load(f)
    
    # Process each indication
    data_dir = Path('data')
    data_dir.mkdir(exist_ok=True)
    
    all_baselines = {}
    
    for indication, results in all_results.items():
        print(f"\nProcessing indication: {indication}")
        
        # Compute baseline statistics
        baseline = compute_baseline_statistics(results)
        all_baselines[indication] = baseline
        
        # Save baseline results
        output_path = data_dir / f'{indication}.pickle'
        with open(output_path, 'wb') as f:
            pickle.dump(baseline, f)
        print(f"✓ Saved baseline statistics to {output_path}")
        
        # Generate visualizations
        create_heatmap_grid(baseline, indication)
    
    # Generate combined plot across all indications
    plot_df = create_combined_baseline_plot(all_baselines)

    plot_df.to_csv(data_dir / 'all_indications_baseline_data.csv', index=False)
    
    print("\n✓ All baselines computed and visualized")


if __name__ == '__main__':
    main()