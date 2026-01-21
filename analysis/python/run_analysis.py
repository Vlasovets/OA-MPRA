#!/usr/bin/env python3
"""
MPRA Analysis - Complete Pipeline Runner
Author: Generated for MPRA analysis
Date: 2026-01-21

Runs the complete per-replicate allelic analysis pipeline.
"""

import subprocess
import sys
from pathlib import Path

def run_script(script_path, description):
    """Run a single analysis script"""
    print("\n" + "="*70)
    print(f"Step: {description}")
    print("="*70 + "\n")
    
    result = subprocess.run(
        [sys.executable, script_path],
        capture_output=False,
        text=True
    )
    
    if result.returncode != 0:
        print(f"\n❌ Error in {script_path}")
        return False
    
    print(f"\n✓ {description} complete!")
    return True

def main():
    """Run complete analysis pipeline"""
    print("="*70)
    print("MPRA PER-REPLICATE ALLELIC ANALYSIS PIPELINE")
    print("="*70)
    print("\nThis will run:")
    print("  1. Data Loading & Filtering (01_load_data.py)")
    print("  2. Per-Replicate Statistical Analysis (02_mpra_analysis.py)")
    print("  3. Visualization (03_visualization.py)")
    
    # Get current directory
    script_dir = Path(__file__).parent
    
    # Define pipeline steps
    steps = [
        (script_dir / "01_load_data.py", "Data Loading & Filtering"),
        (script_dir / "02_mpra_analysis.py", "Per-Replicate Statistical Analysis"),
        (script_dir / "03_visualization.py", "Visualization"),
    ]
    
    # Run each step
    for script_path, description in steps:
        if not script_path.exists():
            print(f"\n❌ Script not found: {script_path}")
            return
        
        success = run_script(script_path, description)
        if not success:
            print("\n❌ Pipeline failed. Stopping.")
            return
    
    # Summary
    print("\n" + "="*70)
    print("✓ COMPLETE ANALYSIS PIPELINE FINISHED!")
    print("="*70)
    print("\nResults saved to:")
    print("  📊 ../results/mpra_analysis/")
    print("     - replicate_level_results.csv")
    print("     - significant_replicate_effects.csv")
    print("     - gain_of_function.csv")
    print("     - loss_of_function.csv")
    print("  📈 ../plots/qc/")
    print("  📈 ../plots/activities/")
    print("\nNext steps:")
    print("  1. Review significant variants")
    print("  2. Validate against GWAS results")
    print("  3. Update analysis report")

if __name__ == "__main__":
    main()
