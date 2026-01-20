# MPRA Count Experiment - Barcode Extraction Fix

## Problem Identified

The MPRAsnakeflow count workflow (`counts_noUMI.smk`) was designed for **overlapping paired-end reads** but your MPRA library has:
- Barcodes at the **start of R1** (first 12bp)
- R1 and R2 **don't overlap** (construct ~151bp, longer than read length)
- The workflow's `--mergeoverlap` parameter only merged 249 reads out of 9.4M

## Solution Implemented

### 1. Pre-extracted Barcodes from R1
```bash
cd /lustre/.../real_data/counts
for f in *_R1_001.fastq.gz; do 
  zcat "$f" | awk 'NR%4==1 {print} NR%4==2 {print substr($0,1,12)} NR%4==3 {print} NR%4==0 {print substr($0,1,12)}' | \
  gzip > "${f%.fastq.gz}_first12bp.fastq.gz"
done
```

### 2. Counted Barcode Occurrences
```bash
for f in *_first12bp.fastq.gz; do 
  sample=$(echo $f | sed 's/_first12bp.fastq.gz//')
  zcat "$f" | awk 'NR%4==2' | sort | uniq -c | awk '{print $2"\t"$1}' | \
  gzip > "${sample}_barcode_counts.tsv.gz"
done
```

### 3. Created Merged DNA/RNA Count Files
```bash
cd results/experiments/exampleCount/counts_manual
for rep in 1 2 3 4; do
  join -t $'\t' -a 1 -a 2 -e 0 -o 0,1.2,2.2 \
    <(zcat OA_${rep}_DNA_final_counts.tsv.gz | sort -k1,1) \
    <(zcat OA_${rep}_RNA_final_counts.tsv.gz | sort -k1,1) | \
    gzip > OA_${rep}.merged.config.default.tsv.gz
done
```

## Results

### Before Fix (Using workflow's merge logic):
- OA_1: **0 overlapping barcodes** between DNA and RNA
- OA_2: **0 overlapping barcodes**
- OA_3: **0 overlapping barcodes**
- OA_4: **0 overlapping barcodes**
- Merged files: 20 bytes (empty)

### After Fix (Pre-extracted 12bp barcodes):
- **OA_1 DNA:** 2,799,363 unique barcodes (531,625 in assignment)
- **OA_1 RNA:** 524,996 unique barcodes (99,912 in assignment)
- **OA_1 Overlap:** **389,656 barcodes** ✅
- Merged files: 7.5-11 MB with millions of entries

### File Locations

**Barcode count files:**
- DNA: `real_data/counts/24L01206[4-7]_S[1-4]_L001_R1_001_barcode_counts.tsv.gz`
- RNA: `real_data/counts/24L00751[2-4]_S[1-3]_L001_R1_001_barcode_counts.tsv.gz`
- RNA: `real_data/counts/24L007692_S5_L001_R1_001_barcode_counts.tsv.gz`

**Final count files (MPRAsnakeflow format):**
- `results/experiments/exampleCount/counts/OA_[1-4]_DNA_final_counts.tsv.gz`
- `results/experiments/exampleCount/counts/OA_[1-4]_RNA_final_counts.tsv.gz`

**Merged files:**
- `results/experiments/exampleCount/counts/OA_[1-4].merged.config.default.tsv.gz`

## Next Steps

The workflow can now continue with:
1. R correlation scripts (should now work with non-empty merged files)
2. Assigned counts (merging with assignment_barcodes)
3. Statistical analysis
4. QC reports

## Technical Notes

- The `feat/allow_single_end_reads` branch has single-end support for **assignment only**, not counts
- For future experiments with non-overlapping reads, consider pre-extracting barcodes before running MPRAsnakeflow
- Runtime was increased from 360 to 1440 minutes for full-length FASTQ processing
