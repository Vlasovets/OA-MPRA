#!/usr/bin/env python3
"""
MPRA Analysis - Step 3: Integration with External Data
Author: Generated for MPRA analysis
Date: 2026-01-20

This script provides a template for integrating MPRA results with:
- Genomic annotations
- Transcription factor binding sites
- ChIP-seq data
- Other functional genomics data
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Paths
DATA_DIR = Path("../results/tables")
OUTPUT_DIR = Path("../results/tables")

def load_mpra_results():
    """Load MPRA results"""
    print("Loading MPRA results...")
    df = pd.read_pickle(DATA_DIR / "bcalm_results_python.pkl")
    print(f"  Loaded {len(df)} variants")
    return df

def example_annotation_integration(df):
    """
    Example: Integrate with genomic annotations
    Replace this with your actual annotation data
    """
    print("\nExample: Integrating with annotations...")
    
    # Example: Add chromosome position if available in variant names
    # Adjust based on your actual variant naming scheme
    
    # This is just a template - modify based on your data
    df_annotated = df.copy()
    
    # Example: Parse variant names (if they contain genomic coordinates)
    # Format: chr1:12345-67890 or similar
    # df_annotated['chr'] = df_annotated.iloc[:, 0].str.extract(r'(chr\d+)')
    # df_annotated['start'] = df_annotated.iloc[:, 0].str.extract(r':(\d+)').astype(float)
    
    print("  Template ready for annotation integration")
    return df_annotated

def analyze_by_category(df, category_col=None):
    """
    Analyze activity by variant category
    Useful if you have different types of variants (e.g., promoters, enhancers)
    """
    print("\nAnalyzing by category...")
    
    if category_col and category_col in df.columns:
        summary = df.groupby(category_col).agg({
            'log2FC': ['mean', 'median', 'std', 'count'],
            'classification': lambda x: (x == 'Enhancer').sum()
        })
        
        print(summary)
        return summary
    else:
        print("  No category column specified")
        return None

def compare_to_controls(df, control_ids=None):
    """
    Compare test variants to negative controls
    """
    print("\nComparing to controls...")
    
    if control_ids is None:
        print("  No control IDs provided (pass as list)")
        return None
    
    # Identify controls
    df_copy = df.copy()
    df_copy['is_control'] = df_copy.iloc[:, 0].isin(control_ids)
    
    # Compare distributions
    controls = df_copy[df_copy['is_control']]
    test = df_copy[~df_copy['is_control']]
    
    print(f"  Controls: n={len(controls)}, mean log2FC={controls['log2FC'].mean():.3f}")
    print(f"  Test: n={len(test)}, mean log2FC={test['log2FC'].mean():.3f}")
    
    # Statistical test
    from scipy import stats
    if len(controls) > 0 and len(test) > 0:
        stat, pval = stats.mannwhitneyu(
            test['log2FC'].dropna(), 
            controls['log2FC'].dropna(),
            alternative='two-sided'
        )
        print(f"  Mann-Whitney U test p-value: {pval:.2e}")
    
    return df_copy

def export_for_igv(df, output_file="variants_for_igv.bed"):
    """
    Export significant variants in BED format for IGV visualization
    Requires variants to have genomic coordinates
    """
    print("\nExporting for IGV...")
    
    # This is a template - adjust based on your actual coordinate format
    sig_variants = df[df['classification'].isin(['Enhancer', 'Silencer'])]
    
    print(f"  Template created for {len(sig_variants)} significant variants")
    print("  Modify this function to parse your variant coordinates")
    
    # Example BED format (uncomment and modify):
    # bed_df = pd.DataFrame({
    #     'chr': sig_variants['chr'],
    #     'start': sig_variants['start'],
    #     'end': sig_variants['end'],
    #     'name': sig_variants.iloc[:, 0],
    #     'score': (sig_variants['log2FC'] * 100).astype(int),
    #     'strand': '.'
    # })
    # bed_df.to_csv(OUTPUT_DIR / output_file, sep='\t', header=False, index=False)

def main():
    """Main integration analysis"""
    print("="*60)
    print("MPRA Analysis - Integration with External Data")
    print("="*60)
    
    # Load MPRA results
    df = load_mpra_results()
    
    # Run analyses
    df_annotated = example_annotation_integration(df)
    analyze_by_category(df, category_col=None)  # Specify your category column
    
    # Compare to controls (provide your control IDs)
    # control_ids = ['control_1', 'control_2', 'control_3']
    # df_with_controls = compare_to_controls(df, control_ids)
    
    # Export for genome browser
    export_for_igv(df)
    
    print("\n" + "="*60)
    print("INTEGRATION ANALYSIS TEMPLATE COMPLETE")
    print("="*60)
    print("\nThis is a template script.")
    print("Modify the functions above based on your specific:")
    print("  - Variant naming/coordinate scheme")
    print("  - External annotation files")
    print("  - Control variant IDs")
    print("  - Analysis goals")

if __name__ == "__main__":
    main()
