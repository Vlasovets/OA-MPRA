#!/usr/bin/env python3
"""
MPRA Analysis - Step 3: Visualization
Author: Generated for MPRA analysis
Date: 2026-01-20
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
import os
from pathlib import Path

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

# Check for analysis suffix from environment
ANALYSIS_SUFFIX = os.environ.get('ANALYSIS_SUFFIX', '')

# Paths
if ANALYSIS_SUFFIX:
    RESULTS_DIR = Path(f"../results/{ANALYSIS_SUFFIX}/mpra_analysis")
    DATA_DIR = Path(f"../results/{ANALYSIS_SUFFIX}/tables")
    PLOT_DIR_QC = Path(f"../plots/{ANALYSIS_SUFFIX}/qc")
    PLOT_DIR_ACT = Path(f"../plots/{ANALYSIS_SUFFIX}/activities")
else:
    RESULTS_DIR = Path("../results/mpra_analysis")
    DATA_DIR = Path("../results/tables")
    PLOT_DIR_QC = Path("../plots/qc")
    PLOT_DIR_ACT = Path("../plots/activities")

# Create directories
PLOT_DIR_QC.mkdir(parents=True, exist_ok=True)
PLOT_DIR_ACT.mkdir(parents=True, exist_ok=True)

def load_data():
    """Load processed data"""
    print("Loading data...")
    
    # Load variant activities
    results = pd.read_pickle(RESULTS_DIR / "variant_activities.pkl")
    
    # Load barcode-level data
    try:
        barcodes = pd.read_pickle(DATA_DIR / "data_filtered.pkl")
    except:
        barcodes = None
        print("  Warning: Barcode data not available")
    
    # Load metadata
    with open(DATA_DIR / "metadata.json", 'r') as f:
        metadata = json.load(f)
    
    print(f"  Loaded {len(results):,} variants")
    if barcodes is not None:
        print(f"  Loaded {len(barcodes):,} barcode observations")
    
    return results, barcodes, metadata

def plot_activity_distribution(df, output_dir):
    """Plot distribution of activity scores"""
    print("\nCreating activity distribution plots...")
    
    if 'log2FC' not in df.columns:
        print("  Skipping: log2FC not found")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Histogram
    axes[0, 0].hist(df['log2FC'].dropna(), bins=50, edgecolor='black', alpha=0.7)
    axes[0, 0].axvline(0, color='red', linestyle='--', linewidth=1)
    axes[0, 0].set_xlabel('log2(RNA/DNA)')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Activity Distribution')
    
    # Box plot by classification
    if 'classification' in df.columns:
        df_plot = df[df['classification'].isin(['Enhancer', 'Silencer', 'Neutral', 'Not significant'])]
        sns.boxplot(data=df_plot, x='classification', y='log2FC', ax=axes[0, 1])
        axes[0, 1].set_xticklabels(axes[0, 1].get_xticklabels(), rotation=45, ha='right')
        axes[0, 1].set_title('Activity by Classification')
    
    # Volcano plot (if p-values available)
    if 'fdr' in df.columns:
        df_vol = df.dropna(subset=['log2FC', 'fdr']).copy()
        df_vol['neg_log10_fdr'] = -np.log10(df_vol['fdr'] + 1e-300)
        
        colors = ['gray'] * len(df_vol)
        if 'classification' in df.columns:
            colors = df_vol['classification'].map({
                'Enhancer': 'red',
                'Silencer': 'blue',
                'Neutral': 'orange',
                'Not significant': 'gray'
            }).fillna('gray')
        
        axes[1, 0].scatter(df_vol['log2FC'], df_vol['neg_log10_fdr'], 
                          c=colors, alpha=0.5, s=20)
        axes[1, 0].axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=1)
        axes[1, 0].axvline(0, color='black', linestyle='--', linewidth=1)
        
        # Add labels for enhancers (red dots)
        if 'classification' in df.columns:
            enhancers = df_vol[df_vol['classification'] == 'Enhancer']
            for idx, row in enhancers.iterrows():
                variant_name = row.iloc[0]  # First column is variant name
                # Truncate long names
                display_name = variant_name[:20] + '...' if len(variant_name) > 20 else variant_name
                axes[1, 0].annotate(display_name, 
                                   xy=(row['log2FC'], row['neg_log10_fdr']),
                                   xytext=(5, 5), textcoords='offset points',
                                   fontsize=6, alpha=0.8,
                                   bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3))
        
        axes[1, 0].set_xlabel('log2(Fold Change)')
        axes[1, 0].set_ylabel('-log10(FDR)')
        axes[1, 0].set_title('Volcano Plot')
    
    # CDF plot
    sorted_fc = np.sort(df['log2FC'].dropna())
    cdf = np.arange(1, len(sorted_fc) + 1) / len(sorted_fc)
    axes[1, 1].plot(sorted_fc, cdf, linewidth=2)
    axes[1, 1].axvline(0, color='red', linestyle='--', linewidth=1)
    axes[1, 1].set_xlabel('log2(RNA/DNA)')
    axes[1, 1].set_ylabel('Cumulative Probability')
    axes[1, 1].set_title('Cumulative Distribution')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / "activity_distributions.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_dir / 'activity_distributions.png'}")

def plot_qc_metrics(df, barcodes, output_dir):
    """Plot QC metrics"""
    print("\nCreating QC plots...")
    
    if barcodes is None:
        print("  Skipping: barcode data not available")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # DNA vs RNA counts
    axes[0, 0].scatter(np.log10(barcodes['DNA'] + 1), 
                       np.log10(barcodes['RNA'] + 1),
                       alpha=0.3, s=10)
    axes[0, 0].set_xlabel('log10(DNA counts + 1)')
    axes[0, 0].set_ylabel('log10(RNA counts + 1)')
    axes[0, 0].set_title('DNA vs RNA Counts')
    
    # Barcode count distributions
    axes[0, 1].hist(np.log10(barcodes['DNA'] + 1), bins=50, alpha=0.5, label='DNA')
    axes[0, 1].hist(np.log10(barcodes['RNA'] + 1), bins=50, alpha=0.5, label='RNA')
    axes[0, 1].set_xlabel('log10(counts + 1)')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].set_title('Count Distributions')
    axes[0, 1].legend()
    
    # Barcodes per variant
    if 'n_barcodes' in df.columns:
        axes[1, 0].hist(df['n_barcodes'].dropna(), bins=30, edgecolor='black')
        axes[1, 0].set_xlabel('Number of Barcodes per Variant')
        axes[1, 0].set_ylabel('Frequency')
        axes[1, 0].set_title('Barcode Coverage per Variant')
    
    # Replicate comparison
    if 'replicate' in barcodes.columns:
        rep_counts = barcodes.groupby('replicate').size()
        axes[1, 1].bar(range(len(rep_counts)), rep_counts.values)
        axes[1, 1].set_xticks(range(len(rep_counts)))
        axes[1, 1].set_xticklabels(rep_counts.index, rotation=45)
        axes[1, 1].set_ylabel('Number of Observations')
        axes[1, 1].set_title('Observations per Replicate')
    
    plt.tight_layout()
    plt.savefig(output_dir / "qc_metrics.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_dir / 'qc_metrics.png'}")

def plot_top_variants(df, output_dir, n=20):
    """Plot top enhancers and silencers"""
    print(f"\nCreating top {n} variants plot...")
    
    if 'log2FC' not in df.columns or 'classification' not in df.columns:
        print("  Skipping: required columns not found")
        return
    
    # Get top enhancers and silencers (sorted by log2FC for effect size)
    enhancers = df[df['classification'] == 'Enhancer'].nlargest(n, 'log2FC')
    silencers = df[df['classification'] == 'Silencer'].nsmallest(n, 'log2FC')
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 8))
    
    # Top enhancers
    if len(enhancers) > 0:
        y_pos = np.arange(len(enhancers))
        axes[0].barh(y_pos, enhancers['log2FC'].values, color='red', alpha=0.7)
        axes[0].set_yticks(y_pos)
        axes[0].set_yticklabels(enhancers.iloc[:, 0].values, fontsize=8)
        axes[0].set_xlabel('log2(Fold Change)')
        axes[0].set_title(f'Top {len(enhancers)} Enhancers')
        axes[0].invert_yaxis()
    
    # Top silencers
    if len(silencers) > 0:
        y_pos = np.arange(len(silencers))
        axes[1].barh(y_pos, silencers['log2FC'].values, color='blue', alpha=0.7)
        axes[1].set_yticks(y_pos)
        axes[1].set_yticklabels(silencers.iloc[:, 0].values, fontsize=8)
        axes[1].set_xlabel('log2(Fold Change)')
        axes[1].set_title(f'Top {len(silencers)} Silencers')
        axes[1].invert_yaxis()
    
    plt.tight_layout()
    plt.savefig(output_dir / "top_variants.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_dir / 'top_variants.png'}")

def plot_dna_rna_heatmap(barcodes, output_dir):
    """Plot heatmap showing DNA vs RNA count patterns"""
    print("\nCreating DNA vs RNA count heatmap...")
    
    if barcodes is None:
        print("  Skipping: barcode data not available")
        return
    
    # Aggregate by replicate
    replicate_counts = barcodes.groupby('replicate')[['DNA', 'RNA']].agg(['sum', 'mean', 'median']).reset_index()
    
    # Create bins for DNA and RNA counts
    dna_bins = [0, 10, 20, 50, 100, 500, np.inf]
    rna_bins = [0, 3, 10, 50, 100, 500, np.inf]
    dna_labels = ['0-10', '10-20', '20-50', '50-100', '100-500', '>500']
    rna_labels = ['0-3', '3-10', '10-50', '50-100', '100-500', '>500']
    
    barcodes_copy = barcodes.copy()
    barcodes_copy['DNA_bin'] = pd.cut(barcodes_copy['DNA'], bins=dna_bins, labels=dna_labels, include_lowest=True)
    barcodes_copy['RNA_bin'] = pd.cut(barcodes_copy['RNA'], bins=rna_bins, labels=rna_labels, include_lowest=True)
    
    # Create contingency table
    heatmap_data = pd.crosstab(barcodes_copy['RNA_bin'], barcodes_copy['DNA_bin'])
    
    # Identify problematic cases (low DNA but high RNA)
    low_dna_high_rna = barcodes_copy[(barcodes_copy['DNA'] < 10) & (barcodes_copy['RNA'] >= 50)]
    n_problematic = len(low_dna_high_rna)
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Heatmap 1: Count frequency
    sns.heatmap(heatmap_data, annot=True, fmt='d', cmap='YlOrRd', ax=axes[0], 
                cbar_kws={'label': 'Number of Observations'})
    axes[0].set_xlabel('DNA Count Bins')
    axes[0].set_ylabel('RNA Count Bins')
    axes[0].set_title('DNA vs RNA Count Distribution\n(All Observations)')
    
    # Heatmap 2: Log-scaled percentage
    heatmap_pct = (heatmap_data / heatmap_data.sum().sum() * 100)
    sns.heatmap(heatmap_pct, annot=True, fmt='.2f', cmap='YlOrRd', ax=axes[1],
                cbar_kws={'label': 'Percentage (%)'})
    axes[1].set_xlabel('DNA Count Bins')
    axes[1].set_ylabel('RNA Count Bins')
    axes[1].set_title(f'DNA vs RNA Count Distribution (%)\nProblematic cases (DNA<10, RNA≥50): {n_problematic:,}')
    
    # Highlight problematic region
    # RNA≥50 is bins 3-5 (indices 3,4,5), DNA<10 is bin 0
    for rna_idx in [3, 4, 5]:
        if rna_idx < len(heatmap_pct):
            rect = plt.Rectangle((0, rna_idx), 1, 1, fill=False, 
                                edgecolor='red', lw=3, linestyle='--')
            axes[1].add_patch(rect)
    
    plt.tight_layout()
    plt.savefig(output_dir / "dna_rna_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_dir / 'dna_rna_heatmap.png'}")
    print(f"  Problematic observations (DNA<10, RNA≥50): {n_problematic:,} ({n_problematic/len(barcodes_copy)*100:.2f}%)")
    
    # Save problematic cases
    if n_problematic > 0:
        problematic_file = output_dir / "problematic_low_dna_high_rna.csv"
        low_dna_high_rna.to_csv(problematic_file, index=False)
        print(f"  Saved problematic cases to: {problematic_file}")

def main():
    """Main function"""
    print("="*60)
    print("MPRA Analysis - Visualization")
    print("="*60)
    
    # Load data
    results, barcodes, metadata = load_data()
    
    # Generate plots
    plot_activity_distribution(results, PLOT_DIR_ACT)
    plot_qc_metrics(results, barcodes, PLOT_DIR_QC)
    plot_top_variants(results, PLOT_DIR_ACT, n=20)  # Sorted by log2FC (effect size)
    plot_dna_rna_heatmap(barcodes, PLOT_DIR_QC)
    
    print("\n✓ Visualization complete!")
    print(f"  QC plots: {PLOT_DIR_QC}")
    print(f"  Activity plots: {PLOT_DIR_ACT}")

if __name__ == "__main__":
    main()
