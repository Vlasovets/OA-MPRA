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
    
    # Try to load per-replicate results first
    try:
        results = pd.read_pickle(RESULTS_DIR / "replicate_level_results.pkl")
        analysis_type = "per_replicate"
        print(f"  Loaded per-replicate results")
    except:
        try:
            results = pd.read_pickle(RESULTS_DIR / "variant_activities.pkl")
            analysis_type = "single"
            print(f"  Loaded variant activities")
        except Exception as e:
            print(f"  Error loading results: {e}")
            results = None
            analysis_type = None
    
    # Load barcode-level data
    try:
        barcodes = pd.read_pickle(DATA_DIR / "data_filtered.pkl")
    except:
        barcodes = None
        print("  Warning: Barcode data not available")
    
    # Load metadata
    try:
        with open(DATA_DIR / "metadata.json", 'r') as f:
            metadata = json.load(f)
    except:
        metadata = {}
        print("  Warning: Metadata not available")
    
    if results is not None:
        print(f"  Loaded {len(results):,} variants/comparisons")
    if barcodes is not None:
        print(f"  Loaded {len(barcodes):,} barcode observations")
    
    return results, barcodes, metadata, analysis_type

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
    axes[0, 0].set_xlabel('log2FC (Mut-WT)')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Allelic Effect Distribution (Per-Replicate)')
    
    # Box plot by classification
    if 'classification' in df.columns:
        df_plot = df[df['classification'].isin(['Gain of function (Mut > WT)', 'Loss of function (Mut < WT)', 'Neutral', 'Not significant'])]
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
                'Gain of function (Mut > WT)': 'red',
                'Loss of function (Mut < WT)': 'blue',
                'Neutral': 'orange',
                'Not significant': 'gray'
            }).fillna('gray')
        
        axes[1, 0].scatter(df_vol['log2FC'], df_vol['neg_log10_fdr'], 
                          c=colors, alpha=0.5, s=20)
        axes[1, 0].axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=1)
        axes[1, 0].axvline(0, color='black', linestyle='--', linewidth=1)
        
        # Add labels for top gain-of-function (red dots) and loss-of-function (blue dots)
        if 'classification' in df.columns:
            # Top gain-of-function - only label top 5 to avoid overlap
            gain_func = df_vol[df_vol['classification'] == 'Gain of function (Mut > WT)'].nlargest(5, 'neg_log10_fdr')
            for i, (idx, row) in enumerate(gain_func.iterrows()):
                variant_name = row['variant_id'] if 'variant_id' in row else str(row.iloc[0])
                # Truncate long names
                display_name = variant_name[:12] + '...' if len(variant_name) > 12 else variant_name
                # Stagger y offset to reduce overlap
                y_offset = 5 + (i % 3) * 8
                axes[1, 0].annotate(display_name, 
                                   xy=(row['log2FC'], row['neg_log10_fdr']),
                                   xytext=(5, y_offset), textcoords='offset points',
                                   fontsize=5, alpha=0.9,
                                   bbox=dict(boxstyle='round,pad=0.2', facecolor='yellow', alpha=0.3),
                                   arrowprops=dict(arrowstyle='-', lw=0.5, alpha=0.5))
            
            # Top loss-of-function - only label top 5 to avoid overlap
            loss_func = df_vol[df_vol['classification'] == 'Loss of function (Mut < WT)'].nlargest(5, 'neg_log10_fdr')
            for i, (idx, row) in enumerate(loss_func.iterrows()):
                variant_name = row['variant_id'] if 'variant_id' in row else str(row.iloc[0])
                # Truncate long names
                display_name = variant_name[:12] + '...' if len(variant_name) > 12 else variant_name
                # Stagger y offset to reduce overlap
                y_offset = 5 + (i % 3) * 8
                axes[1, 0].annotate(display_name, 
                                   xy=(row['log2FC'], row['neg_log10_fdr']),
                                   xytext=(-5, y_offset), textcoords='offset points',
                                   fontsize=5, alpha=0.9,
                                   bbox=dict(boxstyle='round,pad=0.2', facecolor='lightblue', alpha=0.3),
                                   ha='right',
                                   arrowprops=dict(arrowstyle='-', lw=0.5, alpha=0.5))
        
        axes[1, 0].set_xlabel('log2FC (Mut-WT)')
        axes[1, 0].set_ylabel('-log10(FDR)')
        axes[1, 0].set_title('Volcano Plot (Per-Replicate Comparisons)')
    
    # CDF plot
    sorted_fc = np.sort(df['log2FC'].dropna())
    cdf = np.arange(1, len(sorted_fc) + 1) / len(sorted_fc)
    axes[1, 1].plot(sorted_fc, cdf, linewidth=2)
    axes[1, 1].axvline(0, color='red', linestyle='--', linewidth=1)
    axes[1, 1].set_xlabel('log2FC (Mut-WT)')
    axes[1, 1].set_ylabel('Cumulative Probability')
    axes[1, 1].set_title('Cumulative Distribution of Allelic Effects')
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
    
    # Barcodes per variant-replicate
    if 'n_barcodes_WT' in df.columns and 'n_barcodes_Mut' in df.columns:
        all_barcodes = pd.concat([df['n_barcodes_WT'].dropna(), df['n_barcodes_Mut'].dropna()])
        axes[1, 0].hist(all_barcodes, bins=30, edgecolor='black', alpha=0.7)
        axes[1, 0].set_xlabel('Number of Barcodes per Variant-Replicate-Allele')
        axes[1, 0].set_ylabel('Frequency')
        axes[1, 0].set_title('Barcode Coverage Distribution')
    
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
    """Plot top gain and loss of function variants"""
    print(f"\nCreating top {n} variants plot...")
    
    if 'log2FC' not in df.columns or 'classification' not in df.columns:
        print("  Skipping: required columns not found")
        return
    
    # Get top gain and loss of function (sorted by log2FC for effect size)
    gain_func = df[df['classification'] == 'Gain of function (Mut > WT)'].nlargest(n, 'log2FC')
    loss_func = df[df['classification'] == 'Loss of function (Mut < WT)'].nsmallest(n, 'log2FC')
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 8))
    
    # Top gain of function
    if len(gain_func) > 0:
        y_pos = np.arange(len(gain_func))
        # Create labels with variant_WT and variant_Mut if available
        if 'variant_WT' in gain_func.columns and 'variant_Mut' in gain_func.columns:
            labels = [f"{row['variant_id']} ({row['replicate']})" for idx, row in gain_func.iterrows()]
        else:
            labels = gain_func['variant_id'].values if 'variant_id' in gain_func.columns else gain_func.iloc[:, 0].values
        
        axes[0].barh(y_pos, gain_func['log2FC'].values, color='red', alpha=0.7)
        axes[0].set_yticks(y_pos)
        axes[0].set_yticklabels(labels, fontsize=7)
        axes[0].set_xlabel('log2FC (Mut-WT)')
        axes[0].set_title(f'Top {len(gain_func)} Gain of Function')
        axes[0].invert_yaxis()
    
    # Top loss of function
    if len(loss_func) > 0:
        y_pos = np.arange(len(loss_func))
        # Create labels with variant_WT and variant_Mut if available
        if 'variant_WT' in loss_func.columns and 'variant_Mut' in loss_func.columns:
            labels = [f"{row['variant_id']} ({row['replicate']})" for idx, row in loss_func.iterrows()]
        else:
            labels = loss_func['variant_id'].values if 'variant_id' in loss_func.columns else loss_func.iloc[:, 0].values
        
        axes[1].barh(y_pos, loss_func['log2FC'].values, color='blue', alpha=0.7)
        axes[1].set_yticks(y_pos)
        axes[1].set_yticklabels(labels, fontsize=7)
        axes[1].set_xlabel('log2FC (Mut-WT)')
        axes[1].set_title(f'Top {len(loss_func)} Loss of Function')
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

def plot_allelic_comparison(allelic_results, barcodes, metadata, output_dir):
    """
    Create boxplot comparing WT vs Mut activity across all variants.
    Highlights statistically significant variants as red dots.
    
    Parameters:
        allelic_results: DataFrame with allelic comparison results (variant_id, log2FC, fdr, etc.)
        barcodes: DataFrame with barcode-level data (including variant_id, allele, log2FC)
        metadata: Dictionary with variant column name
        output_dir: Path to save plots
    """
    print("\nCreating allelic comparison boxplot...")
    
    if barcodes is None:
        print("  Skipping: barcode data not available")
        return
    
    # Get variant column name from metadata
    variant_col = metadata.get('variant_column', None)
    if variant_col is None or variant_col not in barcodes.columns:
        print(f"  Error: Could not identify variant column")
        print(f"  Metadata keys: {metadata.keys()}")
        print(f"  Available columns: {list(barcodes.columns)}")
        return
    
    print(f"  Parsing variant names from column: {variant_col}")
    
    # Parse variant_id, allele, and replicate from variant names
    # Expected format: <variant_id>_<allele>_<rep>
    def extract_info(name):
        parts = str(name).rsplit('_', 2)  # Split from right
        if len(parts) >= 3:
            variant_id = parts[0]
            allele = parts[1]
            replicate = parts[2]
            return variant_id, allele, replicate
        return None, None, None
    
    parsed = barcodes[variant_col].apply(lambda x: pd.Series(extract_info(x)))
    barcodes = barcodes.copy()
    barcodes['variant_id'] = parsed[0]
    barcodes['allele'] = parsed[1]
    barcodes['rep_parsed'] = parsed[2]
    
    # Remove rows where parsing failed
    barcodes = barcodes.dropna(subset=['variant_id', 'allele'])
    
    print(f"  Parsed {len(barcodes)} observations")
    print(f"  Alleles found: {barcodes['allele'].value_counts().to_dict()}")
    
    # Filter for WT and Mut only
    allelic_data = barcodes[barcodes['allele'].isin(['WT', 'Mut'])].copy()
    
    if len(allelic_data) == 0:
        print("  No WT/Mut data found after filtering")
        return
    
    print(f"  Filtered to {len(allelic_data)} WT/Mut observations")
    
    # Calculate log2FC if not present
    if 'log2FC' not in allelic_data.columns:
        allelic_data['log2FC'] = np.log2((allelic_data['RNA'] + 1) / (allelic_data['DNA'] + 1))
    
    # Find variants with both WT and Mut
    variant_counts = allelic_data.groupby(['variant_id', 'allele']).size().unstack(fill_value=0)
    variants_both = variant_counts[(variant_counts.get('WT', 0) > 0) & 
                                   (variant_counts.get('Mut', 0) > 0)].index.tolist()
    
    allelic_data_both = allelic_data[allelic_data['variant_id'].isin(variants_both)]
    
    print(f"  Found {len(variants_both)} variants with both WT and Mut data")
    print(f"  Total observations: {len(allelic_data_both)} ({allelic_data_both['allele'].value_counts().to_dict()})")
    
    # Get replicate information
    if 'replicate' not in allelic_data_both.columns:
        # Try to get from rep_parsed
        allelic_data_both['replicate'] = allelic_data_both['rep_parsed']
    
    replicates = sorted(allelic_data_both['replicate'].dropna().unique())
    print(f"  Replicates found: {replicates}")
    
    # Create figure with two panels
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Panel 1: Boxplots per replicate (4 pairs of WT/Mut)
    ax1 = axes[0]
    
    # Prepare data for boxplot - one pair (WT, Mut) per replicate
    box_data = []
    box_positions = []
    box_labels = []
    box_colors = []
    
    for i, rep in enumerate(replicates):
        rep_data = allelic_data_both[allelic_data_both['replicate'] == rep]
        
        wt_vals = rep_data[rep_data['allele'] == 'WT']['log2FC'].dropna()
        mut_vals = rep_data[rep_data['allele'] == 'Mut']['log2FC'].dropna()
        
        # Add to box data
        box_data.extend([wt_vals, mut_vals])
        box_positions.extend([i*3, i*3 + 1])  # Space out replicates
        box_labels.extend([f'{rep}\nWT', f'{rep}\nMut'])
        box_colors.extend(['lightblue', 'lightcoral'])
    
    # Create boxplots with outliers visible
    bp = ax1.boxplot(box_data,
                     positions=box_positions,
                     widths=0.6,
                     patch_artist=True,
                     showfliers=True,  # Show outliers
                     flierprops=dict(marker='o', markerfacecolor='gray', 
                                    markersize=4, linestyle='none', 
                                    markeredgecolor='darkgray', alpha=0.5))
    
    # Color boxes
    for patch, color in zip(bp['boxes'], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # Overlay red dots ONLY for significant variant-replicate combinations
    if allelic_results is not None and 'fdr' in allelic_results.columns:
        # Get significant variant-replicate combinations
        sig_results = allelic_results[allelic_results['fdr'] < 0.05].copy()
        
        if len(sig_results) > 0:
            print(f"  Overlaying {len(sig_results)} significant variant-replicate combinations")
            
            for i, rep in enumerate(replicates):
                # Get significant variants for this replicate
                sig_rep = sig_results[sig_results['replicate'] == rep]
                
                if len(sig_rep) > 0:
                    sig_variant_ids = sig_rep['variant_id'].tolist()
                    
                    # Get barcode-level data for significant variants in this replicate
                    sig_data = allelic_data_both[
                        (allelic_data_both['replicate'] == rep) & 
                        (allelic_data_both['variant_id'].isin(sig_variant_ids))
                    ]
                    
                    # WT significant
                    wt_sig = sig_data[sig_data['allele'] == 'WT']['log2FC'].dropna()
                    if len(wt_sig) > 0:
                        wt_x = np.random.normal(i*3, 0.08, size=len(wt_sig))
                        ax1.scatter(wt_x, wt_sig, color='red', alpha=0.6, s=15, zorder=3)
                    
                    # Mut significant
                    mut_sig = sig_data[sig_data['allele'] == 'Mut']['log2FC'].dropna()
                    if len(mut_sig) > 0:
                        mut_x = np.random.normal(i*3 + 1, 0.08, size=len(mut_sig))
                        ax1.scatter(mut_x, mut_sig, color='red', alpha=0.6, s=15, zorder=3)
    
    ax1.set_xticks(box_positions)
    ax1.set_xticklabels(box_labels, fontsize=9)
    ax1.set_ylabel('log2(RNA/DNA)', fontsize=11)
    ax1.set_xlabel('Replicate and Allele', fontsize=11)
    ax1.set_title('Activity Comparison: WT vs Mut by Replicate\n(Red dots = significant variants, FDR<0.05)', fontsize=12)
    ax1.axhline(0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Panel 2: Paired difference plot (log2FC_Mut - log2FC_WT for each variant)
    ax2 = axes[1]
    
    if allelic_results is not None and 'log2FC' in allelic_results.columns:
        # Sort by log2FC difference
        allelic_sorted = allelic_results.sort_values('log2FC')
        
        # Separate significant and non-significant
        sig_mask = allelic_sorted['fdr'] < 0.05 if 'fdr' in allelic_sorted.columns else [False] * len(allelic_sorted)
        
        colors = ['red' if sig else 'gray' for sig in sig_mask]
        alphas = [0.8 if sig else 0.3 for sig in sig_mask]
        
        # Plot
        y_pos = np.arange(len(allelic_sorted))
        
        for i, (idx, row) in enumerate(allelic_sorted.iterrows()):
            ax2.scatter(row['log2FC'], i, color=colors[i], alpha=alphas[i], s=30)
        
        ax2.axvline(0, color='black', linestyle='--', linewidth=1)
        ax2.set_xlabel('log2FC (Mut - WT)')
        ax2.set_ylabel('Variant Index (sorted by effect)')
        ax2.set_title(f'Allelic Effect Size Distribution\n(Red = FDR<0.05, n={sum(sig_mask)})')
        ax2.grid(True, alpha=0.3, axis='x')
        
        # Add labels for top gain-of-function (red dots on right) - only top 5
        gain_sig = allelic_sorted[(allelic_sorted['log2FC'] > 0) & sig_mask].nlargest(5, 'log2FC')
        for i, (idx, row) in enumerate(gain_sig.iterrows()):
            y_idx = allelic_sorted.index.get_loc(idx)
            variant_name = row['variant_id'] if 'variant_id' in row else str(row.iloc[0])
            replicate = row['replicate'] if 'replicate' in row else ''
            display_name = f"{variant_name[:10]}({replicate})" if replicate else (variant_name[:10] + '...' if len(variant_name) > 10 else variant_name)
            # Alternate y position to reduce overlap
            y_offset = (i % 2) * 10 - 5
            ax2.annotate(display_name,
                        xy=(row['log2FC'], y_idx),
                        xytext=(8, y_offset), textcoords='offset points',
                        fontsize=5, alpha=0.9,
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='yellow', alpha=0.3),
                        arrowprops=dict(arrowstyle='-', lw=0.5, alpha=0.5))
        
        # Add labels for top loss-of-function (red dots on left) - only top 5
        loss_sig = allelic_sorted[(allelic_sorted['log2FC'] < 0) & sig_mask].nsmallest(5, 'log2FC')
        for i, (idx, row) in enumerate(loss_sig.iterrows()):
            y_idx = allelic_sorted.index.get_loc(idx)
            variant_name = row['variant_id'] if 'variant_id' in row else str(row.iloc[0])
            replicate = row['replicate'] if 'replicate' in row else ''
            display_name = f"{variant_name[:10]}({replicate})" if replicate else (variant_name[:10] + '...' if len(variant_name) > 10 else variant_name)
            # Alternate y position to reduce overlap
            y_offset = (i % 2) * 10 - 5
            ax2.annotate(display_name,
                        xy=(row['log2FC'], y_idx),
                        xytext=(-8, y_offset), textcoords='offset points',
                        fontsize=5, alpha=0.9,
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='lightblue', alpha=0.3),
                        ha='right',
                        arrowprops=dict(arrowstyle='-', lw=0.5, alpha=0.5))
        
        # Add text annotations for statistics
        n_gain = sum((allelic_sorted['log2FC'] > 0.5) & sig_mask)
        n_loss = sum((allelic_sorted['log2FC'] < -0.5) & sig_mask)
        
        textstr = f'Gain of function: {n_gain}\nLoss of function: {n_loss}'
        ax2.text(0.02, 0.98, textstr, transform=ax2.transAxes, 
                fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(output_dir / "allelic_comparison_boxplot.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_dir / 'allelic_comparison_boxplot.png'}")
    
    # Statistical summary
    wt_mean = allelic_data_both[allelic_data_both['allele'] == 'WT']['log2FC'].mean()
    mut_mean = allelic_data_both[allelic_data_both['allele'] == 'Mut']['log2FC'].mean()
    wt_median = allelic_data_both[allelic_data_both['allele'] == 'WT']['log2FC'].median()
    mut_median = allelic_data_both[allelic_data_both['allele'] == 'Mut']['log2FC'].median()
    
    print(f"  Summary statistics:")
    print(f"    WT:  mean={wt_mean:.3f}, median={wt_median:.3f}, n={len(wt_vals)}")
    print(f"    Mut: mean={mut_mean:.3f}, median={mut_median:.3f}, n={len(mut_vals)}")
    print(f"    Difference (Mut-WT): {mut_mean - wt_mean:.3f}")


def main():
    """Main function"""
    print("="*60)
    print("MPRA Analysis - Visualization")
    print("="*60)
    
    # Load data
    results, barcodes, metadata, analysis_type = load_data()
    
    # Generate plots
    plot_activity_distribution(results, PLOT_DIR_ACT)
    plot_qc_metrics(results, barcodes, PLOT_DIR_QC)
    plot_top_variants(results, PLOT_DIR_ACT, n=20)  # Sorted by log2FC (effect size)
    plot_dna_rna_heatmap(barcodes, PLOT_DIR_QC)
    
    # Add allelic comparison plot if we have per-replicate results
    if analysis_type == "per_replicate":
        plot_allelic_comparison(results, barcodes, metadata, PLOT_DIR_ACT)
    
    print("\n✓ Visualization complete!")
    print(f"  QC plots: {PLOT_DIR_QC}")
    print(f"  Activity plots: {PLOT_DIR_ACT}")

if __name__ == "__main__":
    main()
