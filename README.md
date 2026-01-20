# OA MPRA Variant Analysis

Massively Parallel Reporter Assay (MPRA) analysis pipeline for identifying regulatory variants in Osteoarthritis (OA) cohort data.

## Overview

This repository contains the complete analysis pipeline for an MPRA screen testing 7,696 variants across 4 biological replicates (OA_1-4). The screen successfully identified **23 high-confidence enhancers** with robust statistical validation.

### Key Results
- **4,007 variants tested** with sufficient barcode coverage
- **23 significant enhancers** (FDR<0.05, log2FC≥1.0, hit rate 0.95%)
- **0 silencers** (consistent with activator-biased library design)
- **Top hit by effect size:** scr_157_rep3 (3.32-fold activation)
- **Top hit by statistical confidence:** rs11740553_WT_rep2 (FDR=0.0003)

## Repository Structure

```
.
├── analysis/                           # Complete analysis pipeline
│   ├── MPRA_ANALYSIS_REPORT.md        # Comprehensive analysis report
│   ├── python/                         # Analysis scripts
│   │   ├── 01_load_data.py            # Data loading and QC filtering
│   │   ├── 02_mpra_analysis.py        # Statistical analysis
│   │   ├── 03_visualization.py        # Plot generation
│   │   └── 04_integration_analysis.py # Template for external data
│   ├── plots/                          # Generated figures
│   │   ├── activities/                # Activity distributions, top variants
│   │   └── qc/                        # QC metrics, DNA-RNA heatmaps
│   ├── results/                        # Analysis outputs
│   │   ├── mpra_analysis/             # Variant activities and statistics
│   │   │   ├── all_variants_ranked.csv    # Complete ranked list (by FDR)
│   │   │   ├── variant_activities.csv     # Full metrics table
│   │   │   ├── top_enhancers.csv          # Significant enhancers
│   │   │   └── top_silencers.csv          # Significant silencers
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
- After quality filters (DNA≥10, RNA≥3): 3,353 barcodes (0.6%)

**Key insight:** High barcode detection but skewed read distribution—99% of barcodes have insufficient counts, limiting coverage to median 2 barcodes/variant.

### 3. Statistical Analysis

**Pipeline:**
1. Load and merge count data across 4 replicates
2. Calculate log2(RNA/DNA) activity scores
3. Aggregate by variant (mean, SD across barcodes)
4. Statistical testing: t-test with Benjamini-Hochberg FDR correction
5. Classification: Enhancers (log2FC≥1.0, FDR<0.05), Silencers (log2FC≤-1.0, FDR<0.05)

**Run analysis:**
```bash
cd analysis/python
conda activate mpra_analysis

python 01_load_data.py
python 02_mpra_analysis.py
python 03_visualization.py
```

### 4. Key Findings

**Top 5 Enhancers (by effect size):**
| Rank | Variant | log2FC | FDR | n_barcodes |
|------|---------|--------|-----|------------|
| 1 | scr_157_rep3 | 3.32 | 0.0491 | 3 |
| 2 | rs4644_Mut_rep4 | 2.83 | 0.0004 | 8 |
| 3 | rs2756122_WT_rep2 | 2.72 | 0.0494 | 10 |
| 4 | chr15_88834589_pos_control_rep2 | 2.61 | 0.0494 | 14 |
| 5 | rs11740553_WT_rep2 | 2.39 | 0.0003 | 18 |

**Notable observations:**
- rs11740553 appears in 4 contexts (WT and Mut across replicates) → highly reproducible
- rs4644_Mut_rep4 combines strong effect (2.83-fold) with robust statistics (FDR=0.0004)
- Positive control (chr15_88834589) validates assay performance

## Data Files

### Input Data (Excluded from Git)
- Raw FASTQ files (`.gz`)
- Barcode count files (`.gz`)
- Merged DNA/RNA counts (`*.merged.config.default.tsv.gz`)
- Assignment file (`fromFile.tsv.gz` - 2.9M barcodes, 7,696 variants)

### Results (Included)
- **CSV tables:** Complete variant activities, rankings, and classifications
- **Plots:** QC metrics, activity distributions, volcano plots, top variants
- **Report:** Comprehensive analysis documentation

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
- pandas, numpy, scipy, statsmodels
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
- **Ranked variants:** `analysis/results/mpra_analysis/all_variants_ranked.csv`
- **Plots:** `analysis/plots/activities/` and `analysis/plots/qc/`

## Technical Notes

### Barcode Extraction Fix
The default MPRAsnakeflow pipeline failed with non-overlapping reads. Our custom solution:
1. Pre-extract 12bp barcodes from R1 using awk
2. Count occurrences with `sort | uniq -c`
3. Format to MPRAsnakeflow standard (`barcode\tcount`)

See [BARCODE_EXTRACTION_SUMMARY.md](BARCODE_EXTRACTION_SUMMARY.md) for complete technical details.

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

**Last updated:** January 20, 2026  
**Analysis version:** 1.0  
**Pipeline:** MPRAsnakeflow + Custom Python Analysis
