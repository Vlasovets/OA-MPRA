#!/usr/bin/env python3
"""
MPRA Analysis - Step 2: Allelic Statistical Analysis (Mut vs WT)
Author: Generated for MPRA analysis
Date: 2026-01-21

This script performs allelic MPRA analysis comparing Mut vs WT variants.
Uses negative binomial GLM accounting for DNA normalization.
"""

import pandas as pd
import numpy as np
import json
import os
import sys
import warnings
from pathlib import Path
from scipy import stats
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')

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

def parse_variant_names(data, variant_col):
    """
    Parse variant names to extract variant_id, allele (WT/Mut), and replicate.
    Expected format: <variant_id>_<allele>_<rep>
    Example: rs11740553_WT_rep2 or 10:73797767_ACCTCCCTCACAATTTGCCTACAAGGAAATTCCTTGT_A_Mut_rep1
    """
    print("\nParsing variant names...")
    
    def extract_info(name):
        parts = name.rsplit('_', 2)  # Split from right to get rep and allele
        if len(parts) >= 3:
            variant_id = parts[0]
            allele = parts[1]
            replicate = parts[2]
            return variant_id, allele, replicate
        else:
            return None, None, None
    
    data[['variant_id', 'allele', 'rep']] = data[variant_col].apply(
        lambda x: pd.Series(extract_info(x))
    )
    
    # Remove rows where parsing failed
    n_before = len(data)
    data = data.dropna(subset=['variant_id', 'allele', 'rep'])
    n_after = len(data)
    
    if n_after < n_before:
        print(f"  Warning: Dropped {n_before - n_after} rows with unparseable variant names")
    
    # Check for WT and Mut alleles
    allele_counts = data['allele'].value_counts()
    print(f"  Alleles found: {dict(allele_counts)}")
    
    # Check how many variants have both WT and Mut
    variants_with_both = data.groupby('variant_id')['allele'].apply(
        lambda x: ('WT' in x.values) and ('Mut' in x.values)
    ).sum()
    
    print(f"  Unique variant_ids: {data['variant_id'].nunique():,}")
    print(f"  Variants with both WT and Mut: {variants_with_both:,}")
    
    return data

def calculate_activities(data, variant_col):
    """Calculate RNA/DNA ratios (activity scores) for barcode-level data"""
    print("\nCalculating barcode-level activity scores...")
    
    # Calculate log2(RNA/DNA) with pseudocount
    data['log2FC'] = np.log2((data['RNA'] + 1) / (data['DNA'] + 1))
    
    print(f"  log2FC range: {data['log2FC'].min():.3f} to {data['log2FC'].max():.3f}")
    print(f"  log2FC mean: {data['log2FC'].mean():.3f}")
    
    return data

def prepare_count_matrix(data):
    """
    Prepare count matrix and metadata for allelic analysis.
    
    For MPRA allelic analysis, we need:
    - Count matrix: barcodes × conditions (DNA_rep1, RNA_rep1, DNA_rep2, RNA_rep2, ...)
    - Metadata: variant_id, allele (WT/Mut), replicate for each barcode
    
    Returns:
        counts_df: DataFrame with barcode as index, DNA/RNA counts as columns
        meta_df: DataFrame with barcode metadata (variant_id, allele, replicate)
    """
    print("\nPreparing count matrix for allelic testing...")
    
    # Create unique barcode identifiers
    data['barcode_id'] = data['barcode'] + '_' + data['variant_id'] + '_' + data['allele'] + '_' + data['rep']
    
    # Create count matrix
    counts_df = data[['barcode_id', 'DNA', 'RNA', 'rep']].copy()
    
    # Pivot to get DNA_rep1, RNA_rep1, etc. as columns if there are multiple replicates
    # For now, keep it simple: each row is one barcode in one replicate
    counts_df['DNA_col'] = 'DNA_' + counts_df['rep']
    counts_df['RNA_col'] = 'RNA_' + counts_df['rep']
    
    # Metadata
    meta_df = data[['barcode_id', 'variant_id', 'allele', 'rep']].drop_duplicates()
    
    print(f"  Count matrix shape: {len(counts_df):,} observations")
    print(f"  Unique barcodes: {meta_df['barcode_id'].nunique():,}")
    
    return counts_df, meta_df

def aggregate_alleles_by_variant(data):
    """
    Aggregate barcodes by variant_id and allele for summary statistics.
    This gives us mean activity per variant-allele combination.
    """
    print("\nAggregating by variant_id and allele...")
    
    # Group by variant_id and allele
    allele_stats = data.groupby(['variant_id', 'allele']).agg({
        'log2FC': ['mean', 'median', 'std', 'count'],
        'DNA': ['sum', 'mean'],
        'RNA': ['sum', 'mean'],
        'barcode': 'nunique',
        'rep': 'nunique'
    })
    
    # Flatten column names
    allele_stats.columns = ['_'.join(col).strip() for col in allele_stats.columns.values]
    allele_stats = allele_stats.reset_index()
    
    # Rename for clarity
    allele_stats.rename(columns={
        'log2FC_mean': 'log2FC',
        'log2FC_std': 'log2FC_sd',
        'log2FC_count': 'n_observations',
        'barcode_nunique': 'n_barcodes',
        'rep_nunique': 'n_replicates'
    }, inplace=True)
    
    print(f"  Aggregated to {len(allele_stats):,} variant-allele combinations")
    
    return allele_stats

def test_allelic_effect_nb(data, min_barcodes=3):
    """
    Test for allelic effects (Mut vs WT) using negative binomial GLM.
    UPDATED: Tests within each replicate separately, then meta-analyzes across replicates.
    
    For each variant_id + replicate combination with both WT and Mut data:
    1. Extract RNA and DNA counts for all barcodes in that replicate
    2. Fit negative binomial GLM: log(RNA) ~ allele + offset(log(DNA))
    3. Extract coefficient and p-value
    4. Meta-analyze across replicates using Fisher's method
    
    Parameters:
        data: DataFrame with columns [variant_id, allele, barcode, DNA, RNA, rep]
        min_barcodes: Minimum barcodes per allele per replicate to test
    
    Returns:
        results_df: DataFrame with variant_id, log2FC(Mut vs WT), p-value, FDR
    """
    print("\nTesting allelic effects (Mut vs WT) within replicates...")
    
    # Try to import statsmodels for negative binomial regression
    try:
        import statsmodels.api as sm
        from statsmodels.discrete.discrete_model import NegativeBinomial
    except ImportError:
        print("  Warning: statsmodels not available. Falling back to t-test.")
        return test_allelic_effect_ttest(data, min_barcodes)
    
    # Find variant-replicate combinations with both WT and Mut
    variant_rep_counts = data.groupby(['variant_id', 'rep', 'allele']).size().unstack(fill_value=0)
    
    # Filter for combinations with sufficient data in both alleles
    testable = variant_rep_counts[
        (variant_rep_counts.get('WT', 0) >= min_barcodes) & 
        (variant_rep_counts.get('Mut', 0) >= min_barcodes)
    ]
    
    print(f"  Testable variant-replicate combinations (≥{min_barcodes} barcodes/allele): {len(testable):,}")
    
    # Test each variant-replicate combination
    replicate_results = []
    
    for (variant_id, replicate), row in testable.iterrows():
        variant_rep_data = data[
            (data['variant_id'] == variant_id) & 
            (data['rep'] == replicate)
        ].copy()
        
        # Prepare data for GLM
        variant_rep_data['allele_numeric'] = (variant_rep_data['allele'] == 'Mut').astype(int)
        
        # Filter out zero RNA counts
        variant_rep_data = variant_rep_data[variant_rep_data['RNA'] > 0]
        
        if len(variant_rep_data) < min_barcodes * 2:
            continue
        
        # Calculate mean log2FC for each allele
        wt_data = variant_rep_data[variant_rep_data['allele'] == 'WT']
        mut_data = variant_rep_data[variant_rep_data['allele'] == 'Mut']
        
        wt_mean_log2fc = np.log2((wt_data['RNA'] + 1) / (wt_data['DNA'] + 1)).mean()
        mut_mean_log2fc = np.log2((mut_data['RNA'] + 1) / (mut_data['DNA'] + 1)).mean()
        log2fc_diff = mut_mean_log2fc - wt_mean_log2fc
        
        try:
            # Fit negative binomial GLM with DNA as offset
            X = sm.add_constant(variant_rep_data[['allele_numeric']])
            y = variant_rep_data['RNA'].values
            offset = np.log(variant_rep_data['DNA'].values + 1)
            
            nb_model = NegativeBinomial(y, X, offset=offset)
            nb_result = nb_model.fit(disp=False, maxiter=100)
            
            coef = nb_result.params['allele_numeric']
            pval = nb_result.pvalues['allele_numeric']
            log2fc_model = coef / np.log(2)
            
        except Exception as e:
            pval = np.nan
            log2fc_model = log2fc_diff
        
        replicate_results.append({
            'variant_id': variant_id,
            'replicate': replicate,
            'variant_WT': f"{variant_id}_WT_{replicate}",
            'variant_Mut': f"{variant_id}_Mut_{replicate}",
            'log2FC_Mut': mut_mean_log2fc,
            'log2FC_WT': wt_mean_log2fc,
            'log2FC': log2fc_diff,
            'p_value': pval,
            'n_barcodes_WT': len(wt_data),
            'n_barcodes_Mut': len(mut_data)
        })
    
    if len(replicate_results) == 0:
        print("  Warning: No variant-replicate combinations passed testing criteria")
        return pd.DataFrame()
    
    rep_df = pd.DataFrame(replicate_results)
    
    print(f"  Tested {len(rep_df):,} variant-replicate combinations")
    
    # Apply FDR correction to per-replicate results
    print(f"  Applying FDR correction to per-replicate results...")
    
    valid_pvals = rep_df['p_value'].notna()
    if valid_pvals.sum() > 0:
        reject, fdr, _, _ = multipletests(
            rep_df.loc[valid_pvals, 'p_value'], 
            method='fdr_bh'
        )
        rep_df.loc[valid_pvals, 'fdr'] = fdr
    else:
        rep_df['fdr'] = np.nan
    
    print(f"  Unique variants tested: {rep_df['variant_id'].nunique():,}")
    print(f"  Significant (FDR<0.05): {(rep_df['fdr'] < 0.05).sum():,}")
    
    return rep_df

def test_allelic_effect_ttest(data, min_barcodes=3):
    """
    Fallback method: Test allelic effects using t-test on log2FC values.
    Compares log2FC(Mut) vs log2FC(WT) for each variant.
    
    This is simpler but doesn't properly account for count-based nature of data.
    """
    print("\nTesting allelic effects (Mut vs WT) with t-test (fallback method)...")
    
    # Find variants with both WT and Mut
    variant_allele_counts = data.groupby(['variant_id', 'allele']).size().unstack(fill_value=0)
    
    # Variants that have both WT and Mut
    testable = variant_allele_counts[
        (variant_allele_counts.get('WT', 0) >= min_barcodes) & 
        (variant_allele_counts.get('Mut', 0) >= min_barcodes)
    ].index.tolist()
    
    print(f"  Variants with both WT and Mut (≥{min_barcodes} barcodes each): {len(testable):,}")
    
    results = []
    
    for variant_id in testable:
        variant_data = data[data['variant_id'] == variant_id].copy()
        
        wt_log2fc = variant_data[variant_data['allele'] == 'WT']['log2FC'].values
        mut_log2fc = variant_data[variant_data['allele'] == 'Mut']['log2FC'].values
        
        # Two-sample t-test
        try:
            t_stat, pval = stats.ttest_ind(mut_log2fc, wt_log2fc)
        except:
            pval = np.nan
        
        results.append({
            'variant_id': variant_id,
            'log2FC_Mut': mut_log2fc.mean(),
            'log2FC_WT': wt_log2fc.mean(),
            'log2FC': mut_log2fc.mean() - wt_log2fc.mean(),  # Mut - WT
            'p_value': pval,
            'n_barcodes_WT': len(wt_log2fc),
            'n_barcodes_Mut': len(mut_log2fc),
            'n_barcodes_total': len(wt_log2fc) + len(mut_log2fc)
        })
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) == 0:
        print("  Warning: No variants passed testing criteria")
        return results_df
    
    # FDR correction
    valid_pvals = results_df['p_value'].notna()
    if valid_pvals.sum() > 0:
        reject, fdr, _, _ = multipletests(
            results_df.loc[valid_pvals, 'p_value'], 
            method='fdr_bh'
        )
        results_df.loc[valid_pvals, 'fdr'] = fdr
    else:
        results_df['fdr'] = np.nan
    
    print(f"  Variants tested: {len(results_df):,}")
    print(f"  Significant (FDR<0.05): {(results_df['fdr'] < 0.05).sum():,}")
    
    return results_df

def classify_variants(results_df, fc_threshold=0.5, fdr_threshold=0.05):
    """
    Classify variants based on allelic effect (Mut vs WT).
    Positive log2FC means Mut > WT (gain of function)
    Negative log2FC means Mut < WT (loss of function)
    """
    print("\nClassifying variants...")
    
    def classify(row):
        if pd.isna(row['fdr']):
            return 'Not tested'
        elif row['fdr'] >= fdr_threshold:
            return 'Not significant'
        elif row['log2FC'] > fc_threshold:
            return 'Gain of function (Mut > WT)'
        elif row['log2FC'] < -fc_threshold:
            return 'Loss of function (Mut < WT)'
        else:
            return 'Neutral'
    
    results_df['classification'] = results_df.apply(classify, axis=1)
    
    # Summary
    class_summary = results_df['classification'].value_counts()
    print("\nClassification summary:")
    for class_name, count in class_summary.items():
        print(f"  {class_name}: {count:,}")
    
    return results_df

def main():
    """Main function"""
    print("="*60)
    print("MPRA Allelic Analysis - Mut vs WT Comparison")
    print("="*60)
    
    # Load data
    data, variant_col = load_data()
    
    # Parse variant names to extract variant_id, allele, replicate
    data = parse_variant_names(data, variant_col)
    
    # Calculate barcode-level activities
    data = calculate_activities(data, variant_col)
    
    # Aggregate by variant_id and allele for summary stats
    allele_stats = aggregate_alleles_by_variant(data)
    
    # Test allelic effects (Mut vs WT) using negative binomial GLM (per-replicate)
    print("\n" + "="*60)
    print("Testing Allelic Effects (Mut vs WT) Per-Replicate")
    print("="*60)
    replicate_results = test_allelic_effect_nb(data, min_barcodes=3)
    
    if len(replicate_results) == 0:
        print("\nNo testable variant-replicate combinations found. Exiting.")
        return
    
    # Classify variants (per-replicate)
    replicate_results = classify_variants(replicate_results)
    
    # Sort by FDR, then by absolute log2FC
    replicate_results = replicate_results.sort_values(
        ['fdr', 'log2FC'], 
        ascending=[True, False]
    )
    
    # Save results
    print("\nSaving results...")
    
    # Main per-replicate results
    replicate_results.to_csv(OUTPUT_DIR / "replicate_level_results.csv", index=False)
    replicate_results.to_pickle(OUTPUT_DIR / "replicate_level_results.pkl")
    
    # Save allele-level summary stats
    allele_stats.to_csv(OUTPUT_DIR / "allele_summary_stats.csv", index=False)
    
    # Save significant hits (per-replicate)
    significant = replicate_results[replicate_results['fdr'] < 0.05].copy()
    significant.to_csv(OUTPUT_DIR / "significant_replicate_effects.csv", index=False)
    
    # Save gain and loss of function separately (per-replicate)
    gain_of_function = replicate_results[
        replicate_results['classification'] == 'Gain of function (Mut > WT)'
    ].sort_values('log2FC', ascending=False)
    
    loss_of_function = replicate_results[
        replicate_results['classification'] == 'Loss of function (Mut < WT)'
    ].sort_values('log2FC')
    
    gain_of_function.to_csv(OUTPUT_DIR / "gain_of_function.csv", index=False)
    loss_of_function.to_csv(OUTPUT_DIR / "loss_of_function.csv", index=False)
    
    # Print summary
    print(f"\n✓ Per-replicate allelic analysis complete!")
    print(f"  Results saved to: {OUTPUT_DIR}")
    
    print(f"\nTop 25 Gain of Function (Mut > WT) - Per-Replicate:")
    if len(gain_of_function) > 0:
        display_cols = ['variant_id', 'replicate', 'variant_WT', 'variant_Mut', 'log2FC', 'fdr', 'n_barcodes_WT', 'n_barcodes_Mut']
        if all(col in gain_of_function.columns for col in display_cols):
            display_df = gain_of_function.head(25)[display_cols].copy()
            display_df = display_df.rename(columns={'log2FC': 'log2FC (Mut-WT)'})
            print(display_df.to_string(index=False))
        else:
            # Fallback display
            available_cols = [c for c in display_cols if c in gain_of_function.columns]
            display_df = gain_of_function.head(25)[available_cols].copy()
            if 'log2FC' in display_df.columns:
                display_df = display_df.rename(columns={'log2FC': 'log2FC (Mut-WT)'})
            print(display_df.to_string(index=False))
    else:
        print("  None found")
    
    print(f"\nTop 25 Loss of Function (Mut < WT) - Per-Replicate:")
    if len(loss_of_function) > 0:
        display_cols = ['variant_id', 'replicate', 'variant_WT', 'variant_Mut', 'log2FC', 'fdr', 'n_barcodes_WT', 'n_barcodes_Mut']
        if all(col in loss_of_function.columns for col in display_cols):
            display_df = loss_of_function.head(25)[display_cols].copy()
            display_df = display_df.rename(columns={'log2FC': 'log2FC (Mut-WT)'})
            print(display_df.to_string(index=False))
        else:
            # Fallback display
            available_cols = [c for c in display_cols if c in loss_of_function.columns]
            display_df = loss_of_function.head(25)[available_cols].copy()
            if 'log2FC' in display_df.columns:
                display_df = display_df.rename(columns={'log2FC': 'log2FC (Mut-WT)'})
            print(display_df.to_string(index=False))
    else:
        print("  None found")
    
    print(f"\nNext steps:")
    print(f"  1. Review: {OUTPUT_DIR}/replicate_level_results.csv")
    print(f"  2. Run visualization script to create plots")
    print(f"  3. Compare with GWAS results for validation")

if __name__ == "__main__":
    main()
