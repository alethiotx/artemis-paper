#!/usr/bin/env python
"""
SABCS Consensus Target Analysis

Computes consensus predictions across knowledge graphs for SABCS breast cancer
targets by averaging probabilities across iterations. Generates heatmaps and
clustermaps showing prediction consistency across different models.

Usage:
    consensus_sabcs.py <date> <pickle_file>

Arguments:
    date: Clinical scores date (YYYY-MM-DD format)
    pickle_file: Path to pickled prediction results dictionary

Output:
    - data/all.pickle: Processed consensus results
    - plots/heatmap_<pg>.png: Full heatmaps by pathway genes
    - plots/clustermap_<pg>.png: Clustered heatmaps (targets with all KG predictions)
"""

import sys
import pickle
from pathlib import Path
from typing import Dict, List

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from alethiotx.artemis.clinical.scores import load as load_clinical_scores


# ─── Configuration ───────────────────────────────────────────────────────────

# SABCS data source
SABCS_DATA_PATH = 's3://alethiotx-artemis/data/sabcs24/breast_cancer_targets.csv'

# Knowledge graphs to analyze
KNOWLEDGE_GRAPHS = ['biokg', 'hetionet', 'openbiolink', 'primekg']

# Pathway gene numbers for heatmaps
PG_NUMBERS = ['0', '100', '300']

# Fixed parameters for heatmap generation
# (RF threshold doesn't matter since we average all iterations before thresholding)
HEATMAP_RF_THRESHOLD = '0.7'
HEATMAP_CT_FILTER = 'Unique'

# Heatmap color settings
HEATMAP_VMIN = 0
HEATMAP_VMAX = 1
HEATMAP_CENTER = 0.5
HEATMAP_CMAP = 'vlag'

# Font size for heatmap labels
HEATMAP_FONTSIZE = 12

# Substrings to strip from x-axis labels
XLABEL_STRIP = ['_RotatE', '_rotate']

# Figure dimensions
CLUSTER_HEATMAP_FIGSIZE = (7, 20)
DENDROGRAM_RATIO = (0.3, 0.05)

# Manually defined cluster boundaries (boundary genes in dendrogram order)
# These correspond to the three clusters described in the manuscript
CLUSTER_BOUNDARIES = {
    'high': {'last_gene': 'TF', 'color': '#d73027'},
    'mid': {'last_gene': 'PTK6', 'color': '#ffffff'},
    'low': {'last_gene': None, 'color': '#4575b4'},  # None = rest of genes
}

# Genes explicitly discussed in the manuscript text
MANUSCRIPT_GENES = {
    'BRCA1', 'PTEN', 'MED1', 'BRAF', 'MDM2', 'GATA3', 'YAP1', 'RAD51',
    'FAT1', 'CDK9', 'RUNX2', 'RUNX3',
    'CDK12', 'RAD52', 'TEAD4', 'CD276',
    'MUC1', 'CCR7', 'CCL21', 'CXCR3', 'KAT6A', 'KAT6B',
}

# ─── Helper Functions ────────────────────────────────────────────────────────

def load_sabcs_targets(scores_date: str) -> List[str]:
    """
    Load SABCS breast cancer targets excluding known clinical targets.
    
    Parameters
    ----------
    scores_date : str
        Clinical scores date (YYYY-MM-DD)
    
    Returns
    -------
    List[str]
        SABCS target gene list
    """
    # Load SABCS data
    sabcs_df = pd.read_csv(SABCS_DATA_PATH, index_col=0)
    
    # Remove metadata column if present
    if 'Unique count(Title)' in sabcs_df.columns:
        sabcs_df.drop('Unique count(Title)', axis=1, inplace=True)
    
    # Load known breast cancer targets
    clinical_data = load_clinical_scores(date=scores_date)
    breast_targets = clinical_data[0]['Target Gene']  # First is breast
    
    # Exclude known targets
    sabcs_df = sabcs_df[~sabcs_df.index.isin(breast_targets)]
    
    return sabcs_df.index.tolist()

def compute_consensus(results: Dict, sabcs_targets: List[str]) -> Dict:
    """
    Compute consensus predictions by averaging across iterations.
    
    Parameters
    ----------
    results : Dict
        Nested results dictionary with predictions
    sabcs_targets : List[str]
        SABCS target genes to filter
    
    Returns
    -------
    Dict
        Consensus results by CT, RF threshold, and PG number
    """
    print("Computing consensus predictions...")
    
    for ct in sorted(results.keys()):
        for rf_threshold in sorted(results[ct].keys()):
            for pg_number in sorted(results[ct][rf_threshold].keys()):
                kg_results = []
                
                for kg in KNOWLEDGE_GRAPHS:
                    for embedding in sorted(results[ct][rf_threshold][pg_number][kg].keys()):
                        # Concatenate all iterations for this KG+embedding (sort iteration keys for determinism)
                        iterations_df = pd.concat(
                            [results[ct][rf_threshold][pg_number][kg][embedding][iter_key] 
                             for iter_key in sorted(results[ct][rf_threshold][pg_number][kg][embedding].keys())],
                            axis=1
                        )
                        
                        # Filter to SABCS targets and average across iterations
                        consensus = iterations_df[iterations_df.index.isin(sabcs_targets)].mean(axis=1)
                        kg_results.append(consensus)
                
                # Combine all KG+embedding combos into single DataFrame
                kg_embedding_labels = [
                    f"{kg}_{embedding}"
                    for kg in KNOWLEDGE_GRAPHS
                    for embedding in sorted(results[ct][rf_threshold][pg_number][kg].keys())
                ]
                combined = pd.concat(kg_results, axis=1)
                combined.columns = kg_embedding_labels
                
                results[ct][rf_threshold][pg_number] = combined
    
    return results

def create_clustered_heatmap(data: pd.DataFrame, pg_number: str, output_dir: Path, suffix: str = '') -> None:
    """
    Create clustered heatmap for targets with predictions from all KGs.
    
    Parameters
    ----------
    data : pd.DataFrame
        Consensus prediction data
    pg_number : str
        Pathway gene number for filename
    output_dir : Path
        Output directory for plots
    suffix : str
        Optional suffix for filename (default: '')
    """
    # Strip substrings from column names for cleaner x-axis labels
    rename_map = {}
    for col in data.columns:
        new_col = col
        for s in XLABEL_STRIP:
            new_col = new_col.replace(s, '')
        rename_map[col] = new_col
    data = data.rename(columns=rename_map)

    # Filter to targets with predictions in all KGs (no NaN)
    complete_data = data[~data.isna().any(axis=1)]
    
    if len(complete_data) == 0:
        print(f"  ⚠ No complete targets for pg={pg_number}, skipping clustermap")
        return
    
    # Sort data deterministically to ensure reproducible clustering
    # This ensures the same input order produces the same dendrogram
    complete_data = complete_data.sort_index()
    
    # Adjust figure size for no-biokg clustermaps
    figsize = CLUSTER_HEATMAP_FIGSIZE
    if suffix == '_no_biokg':
        figsize = (CLUSTER_HEATMAP_FIGSIZE[0], CLUSTER_HEATMAP_FIGSIZE[1] * 1.5)
    
    # Only add cluster color bar for the main manuscript figure (pg=0, no suffix)
    use_cluster_colors = (pg_number == '0' and suffix == '')
    
    if use_cluster_colors:
        # First pass: cluster to determine dendrogram order
        g = sns.clustermap(
            complete_data,
            vmin=HEATMAP_VMIN,
            vmax=HEATMAP_VMAX,
            center=HEATMAP_CENTER,
            cmap=HEATMAP_CMAP,
            xticklabels=1,
            yticklabels=1,
            figsize=figsize,
            dendrogram_ratio=DENDROGRAM_RATIO,
            method='average',
            metric='euclidean'
        )
        
        # Assign cluster colors based on dendrogram order
        reordered_genes = list(complete_data.index[g.dendrogram_row.reordered_ind])
        boundaries = list(CLUSTER_BOUNDARIES.values())
        cluster_idx = 0
        gene_colors = {}
        for gene in reordered_genes:
            gene_colors[gene] = boundaries[cluster_idx]['color']
            if boundaries[cluster_idx]['last_gene'] and gene == boundaries[cluster_idx]['last_gene']:
                cluster_idx += 1
        
        # Redraw with row_colors
        plt.close()
        row_colors = complete_data.index.map(gene_colors)
        g = sns.clustermap(
            complete_data,
            vmin=HEATMAP_VMIN,
            vmax=HEATMAP_VMAX,
            center=HEATMAP_CENTER,
            cmap=HEATMAP_CMAP,
            xticklabels=1,
            yticklabels=1,
            figsize=figsize,
            dendrogram_ratio=DENDROGRAM_RATIO,
            method='average',
            metric='euclidean',
            row_colors=row_colors,
        )
        # Bold the gene labels mentioned in the manuscript and set font size
        for label in g.ax_heatmap.get_yticklabels():
            label.set_fontsize(HEATMAP_FONTSIZE)
            if label.get_text() in MANUSCRIPT_GENES:
                label.set_fontweight('bold')
        for label in g.ax_heatmap.get_xticklabels():
            label.set_fontsize(HEATMAP_FONTSIZE)
    else:
        g = sns.clustermap(
            complete_data,
            vmin=HEATMAP_VMIN,
            vmax=HEATMAP_VMAX,
            center=HEATMAP_CENTER,
            cmap=HEATMAP_CMAP,
            xticklabels=1,
            yticklabels=1,
            figsize=figsize,
            dendrogram_ratio=DENDROGRAM_RATIO,
            method='average',
            metric='euclidean'
        )
        for label in g.ax_heatmap.get_yticklabels():
            label.set_fontsize(HEATMAP_FONTSIZE)
        for label in g.ax_heatmap.get_xticklabels():
            label.set_fontsize(HEATMAP_FONTSIZE)
    
    filename = f'clustermap_{pg_number}{suffix}.png' if suffix else f'clustermap_{pg_number}.png'
    plt.savefig(output_dir / filename, bbox_inches='tight')
    plt.close()


def generate_heatmaps(results: Dict) -> None:
    """
    Generate all clustermaps for consensus predictions.
    
    Creates two sets of clustermaps:
    1. Standard: Targets with predictions in all KGs
    2. No-BioKG: Targets present in Hetionet, OpenBioLink, and PrimeKG but not BioKG
    
    Parameters
    ----------
    results : Dict
        Consensus results dictionary
    """
    plots_dir = Path('plots')
    plots_dir.mkdir(exist_ok=True)
    
    print("Generating clustermaps...")
    
    # RF threshold doesn't matter here since we averaged all iterations
    # All runs use the same random seeds (1-10) and we don't threshold
    ct = HEATMAP_CT_FILTER
    rf_threshold = HEATMAP_RF_THRESHOLD
    
    for pg_number in PG_NUMBERS:
        print(f"  pg={pg_number}...")
        
        data = results[ct][rf_threshold][pg_number]
        
        # Standard clustermap (all KG+embedding combos)
        create_clustered_heatmap(data, pg_number, plots_dir)
        print(f"    ✓ clustermap_{pg_number}.png")
        
        # Filtered clustermap (no BioKG)
        # Select columns not starting with 'biokg'
        non_biokg_cols = [c for c in data.columns if not c.startswith('biokg')]
        filtered_data_no_biokg = data[non_biokg_cols]
        filtered_data_no_biokg = filtered_data_no_biokg[~filtered_data_no_biokg.isna().any(axis=1)]
        
        if len(filtered_data_no_biokg) > 0:
            create_clustered_heatmap(filtered_data_no_biokg, pg_number, plots_dir, suffix='_no_biokg')
            print(f"    ✓ clustermap_{pg_number}_no_biokg.png ({len(filtered_data_no_biokg)} targets)")
        else:
            print(f"    ⚠ No targets with all predictions for pg={pg_number}_no_biokg")


# ─── Main Execution ──────────────────────────────────────────────────────────

def main():
    """Compute consensus SABCS predictions and generate visualizations."""
    # Parse command-line arguments
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    
    scores_date = sys.argv[1]
    pickle_file = sys.argv[2]
    
    print(f"Configuration:")
    print(f"  Clinical scores date: {scores_date}")
    print(f"  Pickle file: {pickle_file}")
    
    # Load prediction results (all indications)
    print("Loading prediction results...")
    with open(pickle_file, 'rb') as f:
        all_results = pickle.load(f)
    
    # Extract breast indication data
    if 'breast' not in all_results:
        raise ValueError("'breast' indication not found in pickle file")
    
    results = all_results['breast']
    
    # Load SABCS targets
    print("Loading SABCS targets...")
    sabcs_targets = load_sabcs_targets(scores_date)
    print(f"  Found {len(sabcs_targets)} SABCS targets")
    
    # Compute consensus predictions
    consensus_results = compute_consensus(results, sabcs_targets)
    
    # Save consensus results
    data_dir = Path('data')
    data_dir.mkdir(exist_ok=True)
    
    with open(data_dir / 'all.pickle', 'wb') as f:
        pickle.dump(consensus_results, f)
    print("✓ Saved consensus results to data/all.pickle")
    
    # Generate heatmaps
    generate_heatmaps(consensus_results)
    
    print("✓ All visualizations generated")


def test_with_existing_pickle(pickle_path: str = 'all.pickle') -> None:
    """
    Test visualization generation using existing all.pickle file.
    
    This function loads pre-computed consensus results and regenerates
    visualizations, useful for testing changes to plot generation without
    recomputing consensus predictions.
    
    Parameters
    ----------
    pickle_path : str
        Path to existing all.pickle file (default: 'all.pickle' in current directory)
    
    Usage
    -----
    >>> test_with_existing_pickle('all.pickle')
    >>> test_with_existing_pickle('data/all.pickle')
    """
    print(f"Testing with existing pickle: {pickle_path}")
    
    # Load pre-computed consensus results
    print("Loading consensus results...")
    with open(pickle_path, 'rb') as f:
        consensus_results = pickle.load(f)
    
    # Generate heatmaps
    generate_heatmaps(consensus_results)
    
    print("✓ Test visualizations generated")


if __name__ == '__main__':
    main()
    # test_with_existing_pickle()
