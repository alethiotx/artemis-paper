#!/usr/bin/env python
"""
MeSH Tree Structure Extraction

Downloads and processes MeSH (Medical Subject Headings) tree structure data
from the NLM (National Library of Medicine) and serializes it to a pickle file.
The tree structure represents the hierarchical organization of medical terms.

Usage:
    mesh.py <mesh_file_base>

Arguments:
    mesh_file_base: Base filename for MeSH data (e.g., 'd2024' for 2024 MeSH)

Output:
    - <mesh_file_base>.pkl: Serialized MeSH tree structure
"""

import pickle
import sys
from pathlib import Path

from alethiotx.artemis.mesh.get import tree


# ─── Configuration ───────────────────────────────────────────────────────────

# NLM MeSH data repository URL
MESH_URL_BASE = 'https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/asciimesh/'


# ─── Main Execution ──────────────────────────────────────────────────────────

def main():
    """Download MeSH tree structure and serialize to pickle file."""
    # Parse command-line arguments
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    
    mesh_file_base = sys.argv[1]
    
    # Download and parse MeSH tree
    print(f"Downloading MeSH tree structure (file: {mesh_file_base})...")
    print(f"Source: {MESH_URL_BASE}")
    
    mesh_tree = tree(
        url_base=MESH_URL_BASE,
        file_base=mesh_file_base
    )
    
    print(f"✓ MeSH tree structure retrieved")
    
    # Serialize to pickle file
    output_path = Path(f'{mesh_file_base}.pkl')
    
    with open(output_path, 'wb') as f:
        pickle.dump(mesh_tree, f)
    
    print(f"✓ Saved MeSH tree to {output_path}")


if __name__ == '__main__':
    main()
