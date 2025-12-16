#!/usr/bin/env python
"""
ChEMBL Molecule Data Extraction

Queries the ChEMBL database for molecule information and exports to CSV.
Removes duplicate entries before saving.

Usage:
    chembl.py <chembl_version>

Arguments:
    chembl_version: ChEMBL database version identifier (e.g., 'chembl_34')

Output:
    - <chembl_version>/molecules.csv: Deduplicated molecule data
"""

import sys
from pathlib import Path

from alethiotx.artemis.chembl.query import molecules


# ─── Main Execution ──────────────────────────────────────────────────────────

def main():
    """Query ChEMBL database and extract molecule data."""
    # Parse command-line arguments
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    
    chembl_version = sys.argv[1]
    
    # Query ChEMBL database
    print(f"Querying ChEMBL database (version: {chembl_version})...")
    df = molecules(version=chembl_version)
    print(f"✓ Retrieved {len(df)} molecule records")
    
    # Remove duplicates
    initial_count = len(df)
    df_unique = df.drop_duplicates()
    duplicates_removed = initial_count - len(df_unique)
    
    if duplicates_removed > 0:
        print(f"✓ Removed {duplicates_removed} duplicate records")
    
    # Create output directory
    output_dir = Path(chembl_version)
    output_dir.mkdir(exist_ok=True)
    
    # Save to CSV
    output_path = output_dir / 'molecules.csv'
    df_unique.to_csv(output_path, index=False)
    print(f"✓ Saved {len(df_unique)} unique molecules to {output_path}")


if __name__ == '__main__':
    main()