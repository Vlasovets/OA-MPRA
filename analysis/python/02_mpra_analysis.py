#!/usr/bin/env python3
"""
MPRA Analysis - Step 2: Statistical Analysis with MPRAlib
Author: Generated for MPRA analysis
Date: 2026-01-20
"""

import pandas as pd
import numpy as np
import json
import os
from pathlib import Path
from scipy import stats

# Check for analysis suffix from environment
ANALYSIS_SUFFIX = os.environ.get('ANALYSIS_SUFFIX', '')

# Set paths
if ANALYSIS_SUFFIX:
    DATA_DIR = Path(f"../results/{ANALYSIS_SUFFIX}/tables")
    OUTPUT_DIR = Path(f"../results/{ANALYSIS_SUFFIX}/mpra_analysis")
else:
    DATA_DIR = Path("../results/tables")
    OUTPUT_DIR = Path("../results/mpra_analysis")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def load_data():
    """Load filtered data"""
    print("Loading filtered data...")
    
    data = pd.read_pickle(DATA_DIR / "data_filtered.pkl")
    
    with open(DATA_DIR / "metadata.json", 'r') as f:
        metadata = json.load(f)
    
    variant_col = metadata['variant_column']
    
    print(f"  Loaded {len(data):,} observations")
    print(f"  Variants: {metadata['unique_variants']:,}")
    print(f"  Replicates: {', '.join(metadata['replicates'])}")
    
    return data, variant_col

def calculate_activities(data, variant_col):
    """Calculate RNA/DNA ratios (activity scores)"""
    print("\nCalculating activity scores...")
    
    # Calculate log2(RNA/DNA) with pseudocount
    data['log2FC'] = np.log2((data['RNA'] + 1) / (data['DNA'] + 1))
    
    print(f"  log2FC range: {data['log2FC'].min():.3f} to {data['log2FC'].max():.3f}")
    print(f"  log2FC mean: {data['log2FC'].mean():.3f}")
    
    return data

def aggregate_by_variant(data, variant_col):
    """Aggregate barcodes by variant"""
    print("\nAggregating by variant...")
    
    # Calculate variant-level statistics
    variant_stats = data.groupby(variant_col).agg({
        'log2FC': ['mean', 'median', 'std', 'count'],
        'DNA': ['sum', 'mean'],
        'RNA': ['sum', 'mean'],
        'barcode': 'nunique',
        'replicate': 'nunique'
    })
    
    # Flatten column names
    variant_stats.columns = ['_'.join(col).strip() for col in variant_stats.columns.values]
    variant_stats = variant_stats.reset_index()
    
    # Rename for clarity
    variant_stats.rename(columns={
        'log2FC_mean': 'log2FC',
        'log2FC_std': 'log2FC_sd',
        'log2FC_count': 'n_observations',
        'barcode_nunique': 'n_barcodes',
        'replicate_nunique': 'n_replicates'
    }, inplace=True)
    
    print(f"  Aggregated to {len(variant_stats):,} variants")
    print(f"  Mean barcodes per variant: {variant_stats['n_barcodes'].mean():.1f}")
    
    return variant_stats

def statistical_testing(variant_stats, control_mean=0):
    """Perform statistical tests"""
    print("\nPerforming statistical tests...")
    
    # One-sample t-test against control mean (typically 0 for log2FC)
    # For each variant, test if mean log2FC is significantly different from 0
    
    def calc_pvalue(row):
        if row['n_observations'] < 2 or pd.isna(row['log2FC_sd']) or row['log2FC_sd'] == 0:
            return np.nan
        # Calculate t-statistic
        t_stat = (row['log2FC'] - control_mean) / (row['log2FC_sd'] / np.sqrt(row['n_observations']))
        # Two-tailed p-value
        p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df=row['n_observations']-1))
        return p_value
    
    variant_stats['p_value'] = variant_stats.apply(calc_pvalue, axis=1)
    
    # FDR correction (Benjamini-Hochberg)
    valid_pvals = variant_stats['p_value'].dropna()
    if len(valid_pvals) > 0:
        from statsmodels.stats.multitest import multipletests
        reject, fdr, _, _ = multipletests(valid_pvals, method='fdr_bh')
        variant_stats.loc[valid_pvals.index, 'fdr'] = fdr
    else:
        variant_stats['fdr'] = np.nan
    
    print(f"  Variants tested: {variant_stats['p_value'].notna().sum():,}")
    print(f"  Significant (FDR<0.05): {(variant_stats['fdr'] < 0.05).sum():,}")
    
    return variant_stats

def classify_variants(variant_stats, fc_threshold=0.5, fdr_threshold=0.05):
    """Classify variants as enhancers, silencers, etc."""
    print("\nClassifying variants...")
    
    def classify(row):
        if pd.isna(row['fdr']):
            return 'Not tested'
        elif row['fdr'] >= fdr_threshold:
            return 'Not significant'
        elif row['log2FC'] > fc_threshold:
            return 'Enhancer'
        elif row['log2FC'] < -fc_threshold:
            return 'Silencer'
        else:
            return 'Neutral'
    
    variant_stats['classification'] = variant_stats.apply(classify, axis=1)
    
    # Summary
    class_summary = variant_stats['classification'].value_counts()
    print("\nClassification summary:")
    for class_name, count in class_summary.items():
        print(f"  {class_name}: {count:,}")
    
    return variant_stats

def main():
    """Main function"""
    print("="*60)
    print("MPRA Analysis - Statistical Analysis")
    print("="*60)
    
    # Load data
    data, variant_col = load_data()
    
    # Calculate activities
    data = calculate_activities(data, variant_col)
    
    # Aggregate by variant
    variant_stats = aggregate_by_variant(data, variant_col)
    
    # Statistical testing
    variant_stats = statistical_testing(variant_stats)
    
    # Classify variants
    variant_stats = classify_variants(variant_stats)
    
    # Save results
    print("\nSaving results...")
    variant_stats.to_csv(OUTPUT_DIR / "variant_activities.csv", index=False)
    variant_stats.to_pickle(OUTPUT_DIR / "variant_activities.pkl")
    
    # Save top hits
    enhancers = variant_stats[variant_stats['classification'] == 'Enhancer'].sort_values('log2FC', ascending=False)
    silencers = variant_stats[variant_stats['classification'] == 'Silencer'].sort_values('log2FC')
    
    enhancers.head(50).to_csv(OUTPUT_DIR / "top_enhancers.csv", index=False)
    silencers.head(50).to_csv(OUTPUT_DIR / "top_silencers.csv", index=False)
    
    print(f"\n✓ Statistical analysis complete!")
    print(f"  Results saved to: {OUTPUT_DIR}")
    print(f"\nTop 5 Enhancers:")
    print(enhancers.head(5)[[variant_col, 'log2FC', 'fdr', 'n_barcodes']].to_string(index=False))
    print(f"\nTop 5 Silencers:")
    print(silencers.head(5)[[variant_col, 'log2FC', 'fdr', 'n_barcodes']].to_string(index=False))
    
    print(f"\nNext step: Run 02_custom_plots.py for plots")

if __name__ == "__main__":
    main()
