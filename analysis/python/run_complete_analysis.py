#!/usr/bin/env python3
"""
MPRA Analysis - Complete Pipeline with Before/After Filtering Comparison
Author: Generated for MPRA analysis
Date: 2026-01-20

This script runs the complete analysis twice:
1. Without artifact filtering (baseline)
2. With artifact filtering (cleaned)

Results are saved in separate directories for comparison.
"""

import subprocess
import sys
from pathlib import Path
import shutil

def create_directories(base_dir, suffix):
    """Create organized directory structure"""
    dirs = {
        'results': base_dir / 'results' / suffix,
        'results_mpra': base_dir / 'results' / suffix / 'mpra_analysis',
        'results_tables': base_dir / 'results' / suffix / 'tables',
        'plots_qc': base_dir / 'plots' / suffix / 'qc',
        'plots_activities': base_dir / 'plots' / suffix / 'activities',
    }
    
    for dir_path in dirs.values():
        dir_path.mkdir(parents=True, exist_ok=True)
    
    return dirs

def run_analysis(filter_artifacts, suffix):
    """Run complete analysis pipeline"""
    print("\n" + "="*70)
    print(f"Running Analysis: {suffix.upper()}")
    print("="*70)
    
    base_dir = Path("..").resolve()
    
    # Create directory structure
    print(f"\nCreating directory structure for {suffix}...")
    dirs = create_directories(base_dir, suffix)
    
    # Temporarily modify the scripts to use correct output paths
    # We'll pass the suffix as an environment variable
    import os
    os.environ['ANALYSIS_SUFFIX'] = suffix
    os.environ['FILTER_ARTIFACTS'] = str(filter_artifacts)
    
    # Run each script
    scripts = [
        ('01_load_data.py', 'Data Loading & Filtering'),
        ('02_mpra_analysis.py', 'Statistical Analysis'),
        ('03_visualization.py', 'Visualization'),
    ]
    
    for script, description in scripts:
        print(f"\n{'─'*70}")
        print(f"Step: {description}")
        print(f"{'─'*70}")
        
        result = subprocess.run(
            [sys.executable, script],
            capture_output=False,
            text=True
        )
        
        if result.returncode != 0:
            print(f"\n❌ Error in {script}")
            return False
    
    print(f"\n✓ {suffix.upper()} analysis complete!")
    return True

def create_comparison_summary():
    """Create a summary comparing both analyses"""
    print("\n" + "="*70)
    print("Creating Comparison Summary")
    print("="*70)
    
    import pandas as pd
    import json
    
    summary = {}
    
    for suffix in ['before_filtering', 'after_filtering']:
        try:
            # Load metadata
            meta_file = Path(f"../results/{suffix}/tables/metadata.json")
            if meta_file.exists():
                with open(meta_file, 'r') as f:
                    meta = json.load(f)
                summary[suffix] = meta
            
            # Load variant results
            results_file = Path(f"../results/{suffix}/mpra_analysis/variant_activities.csv")
            if results_file.exists():
                results = pd.read_csv(results_file)
                summary[suffix]['n_significant'] = int((results['fdr'] < 0.05).sum())
                summary[suffix]['n_enhancers'] = int(((results['fdr'] < 0.05) & (results['log2FC'] > 0.5)).sum())
                summary[suffix]['n_silencers'] = int(((results['fdr'] < 0.05) & (results['log2FC'] < -0.5)).sum())
        except Exception as e:
            print(f"  Warning: Could not load {suffix} results: {e}")
    
    # Save comparison
    comparison_file = Path("../results/filtering_comparison_summary.json")
    with open(comparison_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n✓ Comparison summary saved to: {comparison_file}")
    
    # Print comparison
    print("\n" + "="*70)
    print("FILTERING COMPARISON")
    print("="*70)
    
    for suffix in ['before_filtering', 'after_filtering']:
        if suffix in summary:
            print(f"\n{suffix.upper().replace('_', ' ')}:")
            print(f"  Total observations: {summary[suffix].get('total_observations', 'N/A'):,}")
            print(f"  Unique barcodes: {summary[suffix].get('unique_barcodes', 'N/A'):,}")
            print(f"  Unique variants: {summary[suffix].get('unique_variants', 'N/A'):,}")
            print(f"  Significant variants: {summary[suffix].get('n_significant', 'N/A')}")
            print(f"  Enhancers: {summary[suffix].get('n_enhancers', 'N/A')}")
            print(f"  Silencers: {summary[suffix].get('n_silencers', 'N/A')}")

def main():
    """Main function"""
    print("="*70)
    print("MPRA COMPLETE ANALYSIS WITH FILTERING COMPARISON")
    print("="*70)
    print("\nThis script will run the complete analysis pipeline twice:")
    print("  1. WITHOUT artifact filtering (baseline)")
    print("  2. WITH artifact filtering (cleaned)")
    print("\nResults will be saved in separate directories for comparison.")
    
    input("\nPress Enter to continue or Ctrl+C to cancel...")
    
    # Run analysis without filtering
    success1 = run_analysis(filter_artifacts=False, suffix='before_filtering')
    
    if not success1:
        print("\n❌ Analysis without filtering failed. Stopping.")
        return
    
    # Run analysis with filtering
    success2 = run_analysis(filter_artifacts=True, suffix='after_filtering')
    
    if not success2:
        print("\n❌ Analysis with filtering failed.")
        return
    
    # Create comparison summary
    create_comparison_summary()
    
    print("\n" + "="*70)
    print("✓ ALL ANALYSES COMPLETE!")
    print("="*70)
    print("\nResults organized as:")
    print("  📊 results/before_filtering/  - Analysis without artifact filtering")
    print("  📊 results/after_filtering/   - Analysis with artifact filtering")
    print("  📈 plots/before_filtering/    - Plots without filtering")
    print("  📈 plots/after_filtering/     - Plots with filtering")
    print("  📄 results/filtering_comparison_summary.json - Comparison summary")

if __name__ == "__main__":
    main()
