#!/usr/bin/env python
"""
UpSet Plot Generation for Clinical Targets and Pathway Genes

Generates UpSet plots showing overlaps and unique sets across disease indications
for both clinical targets and pathway genes. Creates visualizations for:
- All clinical targets
- Approved clinical targets only
- Unique clinical targets per indication
- Top 100 pathway genes
- Top 300 pathway genes

Usage:
    upset.py <scores_date>

Arguments:
    scores_date: Date string for clinical scores dataset

Output:
    - ct_all.png: All clinical targets overlap
    - ct_approved.png: Approved clinical targets overlap
    - ct_unique.png: Unique clinical targets per indication
    - pg_100.png: Top 100 pathway genes overlap
    - pg_300.png: Top 300 pathway genes overlap
"""

import sys
from typing import Tuple

from matplotlib import pyplot as plt

from alethiotx.artemis.clinical.scores import (
    load as load_clinical,
    approved,
    unique as unique_clinical
)
from alethiotx.artemis.pathway.genes import (
    load as load_pathway,
    unique as unique_pathway
)
from alethiotx.artemis.upset.plot import prepare, create


# ─── Configuration ───────────────────────────────────────────────────────────

# Disease indications in order
INDICATIONS = ['breast', 'lung', 'prostate', 'melanoma', 'bowel', 'diabetes', 'cardiovascular']

# Clinical target plot configurations
CT_CONFIGS = [
    {
        'name': 'ct_all',
        'title': 'All',
        'min_subset_size': 10,
        'filter_fn': None
    },
    {
        'name': 'ct_approved',
        'title': 'Approved',
        'min_subset_size': 10,
        'filter_fn': approved
    },
    {
        'name': 'ct_unique',
        'title': 'Unique',
        'min_subset_size': 1,
        'filter_fn': unique_clinical
    }
]

# Pathway gene configurations
PG_CONFIGS = [
    {'n': 100, 'min_subset_size': 1},
    {'n': 300, 'min_subset_size': 1}
]


# ─── Helper Functions ────────────────────────────────────────────────────────

def create_clinical_upset(
    data: Tuple,
    config: dict,
    mode: str = 'ct'
) -> None:
    """
    Create and save an UpSet plot for clinical targets.
    
    Parameters
    ----------
    data : tuple
        Tuple of indication dataframes
    config : dict
        Plot configuration with name, title, min_subset_size, and filter_fn
    mode : str
        Data mode ('ct' for clinical targets)
    """
    indications = prepare(*data, mode=mode)
    upset_plot = create(indications, min_subset_size=config['min_subset_size'])
    upset_plot.plot()
    plt.suptitle(config['title'])
    plt.savefig(f"{config['name']}.png")
    plt.close()
    print(f"  ✓ {config['name']}.png")


def create_pathway_upset(
    scores_date: str,
    n: int,
    min_subset_size: int
) -> None:
    """
    Create and save an UpSet plot for pathway genes.
    
    Parameters
    ----------
    scores_date : str
        Date string for clinical scores dataset
    n : int
        Number of top pathway genes to include
    min_subset_size : int
        Minimum size for subsets to display
    """
    # Load and filter pathway genes
    data = load_pathway(date=scores_date, n=n)
    data_unique = unique_pathway(list(data))
    
    # Create and save plot
    indications = prepare(*data_unique, mode='pg')
    upset_plot = create(indications, min_subset_size=min_subset_size)
    upset_plot.plot()
    plt.suptitle(f'{n} pathway genes')
    plt.savefig(f"pg_{n}.png")
    plt.close()
    print(f"  ✓ pg_{n}.png")


# ─── Main Execution ──────────────────────────────────────────────────────────

def main():
    """Generate all UpSet plots for clinical targets and pathway genes."""
    # Parse command-line arguments
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    
    scores_date = sys.argv[1]
    
    # Load clinical scores
    print(f"Loading clinical scores (date: {scores_date})...")
    clinical_data = load_clinical(date=scores_date)
    print(f"✓ Loaded data for {len(clinical_data)} indications")
    
    # Generate clinical target UpSet plots
    print("\nGenerating clinical target UpSet plots...")
    for config in CT_CONFIGS:
        # Apply filter if specified
        if config['filter_fn'] is not None:
            data = config['filter_fn'](list(clinical_data))
        else:
            data = clinical_data
        
        create_clinical_upset(data, config, mode='ct')
    
    # Generate pathway gene UpSet plots
    print("\nGenerating pathway gene UpSet plots...")
    for pg_config in PG_CONFIGS:
        create_pathway_upset(
            scores_date,
            n=pg_config['n'],
            min_subset_size=pg_config['min_subset_size']
        )
    
    print("\n✓ All UpSet plots generated successfully")


if __name__ == '__main__':
    main()