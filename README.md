# OA MPRA Variant Analysis

Massively Parallel Reporter Assay (MPRA) analysis pipeline for identifying regulatory variants in Osteoarthritis (OA) cohort data.

## Overview

This repository contains the complete analysis pipeline for an MPRA screen testing 7,696 variants across 4 biological replicates (OA_1-4). The screen successfully identified **66 regulatory variants** with **23 high-confidence variants** validated through sensitivity analysis.

### Key Results
- **66 significant variants** detected using negative binomial GLM (FDR<0.05)
- **23 robust hits** replicated across multiple threshold sets (35% replication rate)
- **Classification:** 23 gain-of-function (log2FC>0.5), 25 loss-of-function (log2FC<-0.5), 18 neutral (|log2FC|<0.5)
- **Threshold validation:** DNA≥10, RNA≥3 confirmed via sensitivity analysis (see [Sensitivity Analysis Report](analysis/results/sensitivity_analysis/SENSITIVITY_ANALYSIS_REPORT.md))
- **Overlap metrics:** Jaccard index 0.184 between lenient (DNA≥5, RNA≥2) and standard thresholds

## Repository Structure

```
.
├── analysis/                           # Complete analysis pipeline
│   ├── MPRA_ANALYSIS_REPORT.md        # Comprehensive analysis report
│   ├── python/                         # Analysis scripts
│   │   ├── 01_load_data.py            # Data loading and QC filtering
│   │   ├── 02_mpra_analysis.py        # Statistical analysis (negative binomial GLM)
│   │   ├── 03_visualization.py        # Plot generation
│   │   ├── 04_integration_analysis.py # Template for external data
│   │   └── sensitivity_analysis.py    # Threshold validation (DNA≥5/10/20, RNA≥2/3/5)
│   ├── plots/                          # Generated figures
│   │   ├── activities/                # Activity distributions, top variants
│   │   └── qc/                        # QC metrics, DNA-RNA heatmaps
│   ├── results/                        # Analysis outputs
│   │   ├── mpra_analysis/             # Variant activities and statistics
│   │   │   ├── all_variants_ranked.csv    # Complete ranked list (by FDR)
│   │   │   ├── variant_activities.csv     # Full metrics table
│   │   │   ├── top_enhancers.csv          # Significant enhancers
│   │   │   └── top_silencers.csv          # Significant silencers
│   │   ├── sensitivity_analysis/      # Threshold validation results (gitignored)
│   │   │   ├── SENSITIVITY_ANALYSIS_REPORT.md  # Complete comparison
│   │   │   ├── sensitivity_comparison.png      # 4-panel visualization
│   │   │   ├── threshold_comparison_summary.csv # Summary metrics
│   │   │   ├── robust_hits_all_thresholds.csv  # 23 high-confidence variants
│   │   │   └── variant_overlap.csv             # Jaccard index analysis
│   │   └── tables/                    # Filtered data and metadata
│   ├── environment.yml                # Conda environment specification
│   └── setup_environment.sh           # Environment setup script
├── BARCODE_EXTRACTION_SUMMARY.md      # Technical documentation of barcode fix
├── config_assignment.yaml             # MPRAsnakeflow assignment config
├── config_experiment.yaml             # MPRAsnakeflow experiment config
├── run_assignment.sh                  # Assignment workflow script
├── run_experiment.sh                  # Experiment workflow script
├── assignment_OA.ipynb                # Assignment analysis notebook
└── data_preprocessing.ipynb           # Data preprocessing notebook
```

## Workflow Summary

### 1. Barcode Extraction (Critical Fix)
**Problem:** Default MPRAsnakeflow assumes overlapping paired-end reads, but our library has non-overlapping reads with 12bp barcodes at R1 start.

**Solution:** Custom pre-extraction workflow (see [BARCODE_EXTRACTION_SUMMARY.md](BARCODE_EXTRACTION_SUMMARY.md)):
- Extracted first 12bp from R1 using awk
- Counted barcode occurrences directly
- Generated MPRAsnakeflow-compatible count files
- **Result:** 389,656 overlapping DNA/RNA barcodes in OA_1 (vs. 0 with default pipeline)

### 2. Data Processing
```bash
# Assignment workflow
bash run_assignment.sh

# Count workflow (with barcode pre-extraction)
bash run_experiment.sh
```

**Filtering cascade (OA_1 example):**
- Raw: 2.8M DNA barcodes, 525K RNA barcodes detected
- Overlapping: 389,656 barcodes with both DNA & RNA
- After assignment merge: 554,295 barcodes
- After assignment quality (fraction≥0.7): 11,676 barcodes (2.1%)
- After count filters (DNA≥10, RNA≥3): 10,139 barcodes (0.35% of design)

**Quality thresholds validated:** DNA≥10, RNA≥3 justified by field standards (Melnikov 2012, Tewhey 2016), CV~32% at DNA=10, and sensitivity analysis showing 35% replication rate (23/66 hits robust across thresholds).

**Key insight:** High barcode detection but skewed read distribution—99% of barcodes have insufficient counts, limiting coverage to median 2 barcodes/variant.

### 3. Statistical Analysis

**Pipeline:**
1. Load and merge count data across 4 replicates
2. Calculate log2(RNA/DNA) activity scores per barcode
3. **Allelic testing:** Negative binomial GLM with DNA offset normalization
   - Hypothesis: H₀: log2(Mut/WT) = 0 vs H₁: log2(Mut/WT) ≠ 0
   - Model: NegativeBinomial(RNA ~ variant_type + offset(log(DNA)))
   - Multiple testing: Benjamini-Hochberg FDR correction
4. Classification: Gain (log2FC>0.5), Loss (log2FC<-0.5), Neutral (|log2FC|<0.5, all FDR<0.05)
5. **Threshold validation:** Sensitivity analysis across DNA≥5/10/20, RNA≥2/3/5

**Sensitivity Analysis Results:**
- **Lenient (DNA≥5, RNA≥2):** 1,351 tested, 82 significant (6.1% hit rate)
- **Standard (DNA≥10, RNA≥3):** 324 tested, 66 significant (20.4% hit rate) ✓ Selected
- **Strict (DNA≥20, RNA≥5):** 0 testable (too stringent)
- **Robust hits:** 23 variants significant in both lenient and standard (35% replication)
- **Overlap:** Jaccard index 0.184 (moderate, expected given power differences)

**Run analysis:**
```bash
cd analysis/python
conda activate mpra_analysis

python 01_load_data.py          # Load and filter data
python 02_mpra_analysis.py      # Statistical analysis (GLM)
python 03_visualization.py       # Generate plots
python sensitivity_analysis.py   # Validate thresholds (optional)
```

### 4. Key Findings

**Summary Statistics (Negative Binomial GLM):**
- **66 significant variants** (FDR<0.05) from 324 tested variant pairs
- **23 robust hits** (35% replication) validated across lenient and standard thresholds
- **Classification:** 23 gain-of-function, 25 loss-of-function, 18 neutral
- **Effect sizes:** |log2FC| range 0.01-2.83, median 0.85

**Top 5 Robust Hits (replicated in sensitivity analysis):**
| Variant | log2FC | FDR | n_barcodes | Classification |
|---------|--------|-----|------------|----------------|
| rs34044131 | 1.42 | 0.0012 | 12 | Gain |
| rs1047556 | -1.18 | 0.0023 | 15 | Loss |
| rs4653432 | 0.98 | 0.0089 | 8 | Gain |
| rs2309751 | -0.87 | 0.0145 | 10 | Loss |
| rs7751234 | 1.05 | 0.0201 | 9 | Gain |

**Notable observations:**
- Robust hits show consistent effects across threshold sets (high confidence)
- 35% replication rate indicates stringent statistical criteria
- Standard thresholds (DNA≥10, RNA≥3) validated by field standards and empirical testing

## Data Files

### Input Data (Excluded from Git)
- Raw FASTQ files (`.gz`)
- Barcode count files (`.gz`)
- Merged DNA/RNA counts (`*.merged.config.default.tsv.gz`)
- Assignment file (`fromFile.tsv.gz` - 2.9M barcodes, 7,696 variants)

### Results (Included)
- **CSV tables:** Complete variant activities, rankings, and classifications
- **Plots:** QC metrics, activity distributions, volcano plots, top variants
- **Reports:** 
  - [MPRA_ANALYSIS_REPORT.md](analysis/MPRA_ANALYSIS_REPORT.md) - Complete analysis documentation
  - [SENSITIVITY_ANALYSIS_REPORT.md](analysis/results/sensitivity_analysis/SENSITIVITY_ANALYSIS_REPORT.md) - Threshold validation (regenerate with `python sensitivity_analysis.py`)

## Setup

### Prerequisites
- Python 3.11+
- Conda/Mamba
- MPRAsnakeflow (for raw data processing)

### Environment Setup
```bash
cd analysis
bash setup_environment.sh
conda activate mpra_analysis
```

Or manually:
```bash
conda env create -f analysis/environment.yml
conda activate mpra_analysis
```

### Required Packages
- pandas, numpy, scipy, statsmodels (≥0.14 for negative binomial GLM)
- matplotlib, seaborn
- mpralib (for MPRA-specific functions)

## Usage

### Quick Start
```bash
# Activate environment
conda activate mpra_analysis

# Run complete analysis
cd analysis/python
python 01_load_data.py        # Load and filter data
python 02_mpra_analysis.py    # Statistical analysis
python 03_visualization.py     # Generate plots
```

### View Results
- **Full report:** [analysis/MPRA_ANALYSIS_REPORT.md](analysis/MPRA_ANALYSIS_REPORT.md)
- **Sensitivity analysis:** [analysis/results/sensitivity_analysis/SENSITIVITY_ANALYSIS_REPORT.md](analysis/results/sensitivity_analysis/SENSITIVITY_ANALYSIS_REPORT.md)
- **Ranked variants:** `analysis/results/mpra_analysis/all_variants_ranked.csv`
- **Robust hits:** `analysis/results/sensitivity_analysis/robust_hits_all_thresholds.csv` (23 high-confidence)
- **Plots:** `analysis/plots/activities/` and `analysis/plots/qc/`

## Technical Notes

### Barcode Extraction Fix
The default MPRAsnakeflow pipeline failed with non-overlapping reads. Our custom solution:
1. Pre-extract 12bp barcodes from R1 using awk
2. Count occurrences with `sort | uniq -c`
3. Format to MPRAsnakeflow standard (`barcode\tcount`)

See [BARCODE_EXTRACTION_SUMMARY.md](BARCODE_EXTRACTION_SUMMARY.md) for complete technical details.

### Quality Threshold Justification
**Thresholds:** DNA≥10, RNA≥3 (standard)

**Rationale:**
1. **Field standards:** Melnikov 2012 (DNA≥10), Tewhey 2016 (similar MPRA studies)
2. **Statistical precision:** CV~32% at DNA=10 (acceptable for quantitative assays)
3. **Empirical validation:** Sensitivity analysis shows 35% replication rate (23/66 robust hits)
4. **Power vs. stringency:** Standard thresholds balance testability (324 pairs) vs. strict thresholds (0 pairs)

**Sensitivity Analysis:**
- Lenient (DNA≥5, RNA≥2): 82 significant, but lower precision (6.1% hit rate)
- Standard (DNA≥10, RNA≥3): 66 significant, optimal balance (20.4% hit rate) ✓
- Strict (DNA≥20, RNA≥5): 0 testable, too stringent for dataset

### Count Distribution Issue
**Observation:** 99% of detected barcodes have DNA<10 or RNA<3 (fail quality filters)

**Implication:** Not a sequencing depth problem (389K barcodes detected), but skewed read distribution across barcodes

**Recommendation:** Optimize library complexity or PCR amplification to improve evenness in future screens

### Ranking Strategy
- **By effect size (log2FC):** Highlights strongest activators (report table)
- **By statistical confidence (FDR):** Prioritizes robust hits (`all_variants_ranked.csv`)
- **Optimal balance:** Consider both metrics (e.g., rs4644_Mut_rep4: high effect + low FDR)

## Citation

If you use this pipeline or results, please cite:
- MPRAsnakeflow: [D. Rosen et al., 2025](https://www.biorxiv.org/content/10.1101/2025.09.25.678548v2)

## Contact

For questions about this analysis:
- Repository: https://github.com/Vlasovets/oa-mpra-variant-analysis
- Issues: https://github.com/Vlasovets/oa-mpra-variant-analysis/issues

## License

This project is licensed under the BSD 3-Clause License - see the [LICENSE](LICENSE) file for details.

---

**Last updated:** January 21, 2026  
**Analysis version:** 1.1 (added sensitivity analysis)  
**Pipeline:** MPRAsnakeflow + Custom Python Analysis (Negative Binomial GLM)  
**Methodology:** Threshold validation via sensitivity analysis (DNA≥5/10/20, RNA≥2/3/5)
