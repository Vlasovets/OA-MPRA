#!/usr/bin/env python3
"""
MPRA Analysis - Step 1: Load and Prepare Count Data with MPRAlib
Author: Generated for MPRA analysis
Date: 2026-01-20
"""

import pandas as pd
import numpy as np
import gzip
import os
from pathlib import Path

# Check for analysis suffix from environment (for comparison runs)
ANALYSIS_SUFFIX = os.environ.get('ANALYSIS_SUFFIX', '')
FILTER_ARTIFACTS = os.environ.get('FILTER_ARTIFACTS', 'True').lower() == 'true'

# Set paths
DATA_DIR = Path("../data")
if ANALYSIS_SUFFIX:
    OUTPUT_DIR = Path(f"../results/{ANALYSIS_SUFFIX}/tables")
else:
    OUTPUT_DIR = Path("../results/tables")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def load_count_files():
    """Load merged count files for all replicates"""
    print("Loading count data from merged files...")
    
    counts_dir = DATA_DIR / "counts"
    replicates = []
    
    for rep in [1, 2, 3, 4]:
        file_path = counts_dir / f"OA_{rep}.merged.config.default.tsv.gz"
        
        if not file_path.exists():
            print(f"  WARNING: File not found: {file_path}")
            continue
        
        # Read gzipped TSV
        with gzip.open(file_path, 'rt') as f:
            df = pd.read_csv(f, sep='\t', names=['barcode', 'DNA', 'RNA'])
        
        df['replicate'] = f'OA_{rep}'
        replicates.append(df)
        print(f"  Loaded OA_{rep}: {len(df):,} barcodes")
    
    if not replicates:
        raise FileNotFoundError("No count files found!")
    
    # Combine all replicates
    combined = pd.concat(replicates, ignore_index=True)
    print(f"\nTotal observations: {len(combined):,}")
    
    return combined

def load_assignment():
    """Load barcode-to-variant assignment"""
    print("\nLoading assignment data...")
    
    assignment_file = DATA_DIR / "assignments" / "fromFile.tsv.gz"
    
    if not assignment_file.exists():
        raise FileNotFoundError(f"Assignment file not found: {assignment_file}")
    
    # Read gzipped TSV
    with gzip.open(assignment_file, 'rt') as f:
        assignments = pd.read_csv(f, sep='\t')
    
    print(f"  Loaded {len(assignments):,} barcode assignments")
    print(f"  Columns: {', '.join(assignments.columns.tolist())}")
    
    # Determine barcode and variant columns (typically first two columns)
    barcode_col = assignments.columns[0]
    variant_col = assignments.columns[1] if len(assignments.columns) > 1 else assignments.columns[0]
    
    print(f"  Barcode column: {barcode_col}")
    print(f"  Variant column: {variant_col}")
    print(f"  Unique variants: {assignments[variant_col].nunique():,}")
    
    return assignments, barcode_col, variant_col

def merge_and_filter(counts, assignments, barcode_col, variant_col, min_dna=10, min_rna=3, 
                     filter_artifacts=True, artifact_dna_threshold=10, artifact_rna_threshold=50):
    """Merge counts with assignments and apply filters"""
    print("\nMerging counts with assignments...")
    
    # Merge
    data = counts.merge(assignments, left_on='barcode', right_on=barcode_col, how='inner')
    print(f"  After merging: {len(data):,} observations")
    print(f"  Matched barcodes: {data['barcode'].nunique():,}")
    
    # Check for artifacts BEFORE any filtering
    if filter_artifacts:
        print(f"\nIdentifying artifacts (DNA<{artifact_dna_threshold} & RNA≥{artifact_rna_threshold}) in RAW data...")
        artifacts = (data['DNA'] < artifact_dna_threshold) & (data['RNA'] >= artifact_rna_threshold)
        n_artifacts = artifacts.sum()
        
        if n_artifacts > 0:
            print(f"  Found {n_artifacts:,} artifact observations ({n_artifacts/len(data)*100:.2f}%)")
            # Remove artifacts first
            data = data[~artifacts].copy()
            print(f"  After removing artifacts: {len(data):,} observations")
        else:
            print(f"  No artifacts found (0.00%)")
    
    # Apply basic filters
    print(f"\nApplying basic filters (DNA≥{min_dna}, RNA≥{min_rna})...")
    n_before = len(data)
    data_filtered = data[(data['DNA'] >= min_dna) & (data['RNA'] >= min_rna)].copy()
    n_after = len(data_filtered)
    n_removed = n_before - n_after
    
    print(f"  Removed {n_removed:,} low-count observations ({n_removed/n_before*100:.2f}%)")
    print(f"  Remaining: {n_after:,} observations")
    print(f"  Final barcodes: {data_filtered['barcode'].nunique():,}")
    print(f"  Final variants: {data_filtered[variant_col].nunique():,}")
    
    return data_filtered

def summarize_by_replicate(data, variant_col):
    """Generate summary statistics by replicate"""
    print("\nSummary by replicate:")
    
    summary = data.groupby('replicate').agg({
        'barcode': 'nunique',
        'DNA': ['sum', 'mean', 'median'],
        'RNA': ['sum', 'mean', 'median'],
        variant_col: 'nunique'
    }).round(2)
    
    summary.columns = ['_'.join(col).strip() for col in summary.columns.values]
    print(summary)
    
    return summary

def main():
    """Main function"""
    print("="*60)
    print("MPRA Analysis - Load and Prepare Data")
    print("="*60)
    
    # Load data
    counts = load_count_files()
    assignments, barcode_col, variant_col = load_assignment()
    
    # Merge and filter (with optional artifact removal)
    data_filtered = merge_and_filter(counts, assignments, barcode_col, variant_col,
                                     filter_artifacts=FILTER_ARTIFACTS,
                                     artifact_dna_threshold=10,
                                     artifact_rna_threshold=50)
    
    # Summary statistics
    summary = summarize_by_replicate(data_filtered, variant_col)
    
    # Save processed data
    print("\nSaving processed data...")
    data_filtered.to_pickle(OUTPUT_DIR / "data_filtered.pkl")
    data_filtered.to_csv(OUTPUT_DIR / "data_filtered.csv.gz", index=False, compression='gzip')
    summary.to_csv(OUTPUT_DIR / "summary_by_replicate.csv")
    
    # Save metadata
    metadata = {
        'total_observations': len(data_filtered),
        'unique_barcodes': int(data_filtered['barcode'].nunique()),
        'unique_variants': int(data_filtered[variant_col].nunique()),
        'replicates': sorted(data_filtered['replicate'].unique().tolist()),
        'barcode_column': barcode_col,
        'variant_column': variant_col
    }
    
    import json
    with open(OUTPUT_DIR / "metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"\n✓ Data preparation complete!")
    print(f"  Files saved to: {OUTPUT_DIR}")
    print(f"\nNext step: Run 02_mpra_analysis.py for statistical analysis")

if __name__ == "__main__":
    main()
