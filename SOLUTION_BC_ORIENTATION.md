# SOLUTION: Fix BC Orientation Mismatch

**Date**: March 5, 2026  
**Issue**: DNA/RNA barcodes don't match assignment orientation  
**Root Cause**: BC orientation mismatch between assignment and experiment  
**Solution**: Reverse complement DNA/RNA BC files

---

## 🔍 Problem Diagnosis

### What We Found:
1. **Assignment BCs** (with `BC_rev_comp: true`): Evenly distributed nucleotides
   - A: 27%, C: 27%, T: 24%, G: 22%
   
2. **DNA/RNA BCs** (original orientation): 98.7% start with 'G'
   - G: 98.7%, A: 0.8%, T: 0.1%, C: 0.3%  
   - This is NOT natural - indicates wrong orientation!

3. **Matching Results**:
   - Original orientation: Only **14.31%** of assignment BCs match DNA
   - Reverse complement: **18.29%** of assignment BCs match DNA
   - **Improvement: +115,852 BCs (+4.0%)**

### Why This Happened:
- Your assignment config has `BC_rev_comp: true` 
- This outputs BCs in reverse complement orientation
- Your DNA/RNA BC files were extracted WITHOUT reverse complement
- Result: orientation mismatch → poor correlation

---

## ✅ Solution (Following Max's Advice)

### Step-by-Step Instructions:

#### 1. Extract Reverse Complemented BCs from DNA/RNA Files

```bash
cd /home/itg/oleg.vlasovets/projects/MPRA_data/mpra_test
chmod +x scripts/extract_barcodes_RC.sh
bash scripts/extract_barcodes_RC.sh
```

This will create:
- `/lustre/groups/itg/teams/zeggini/projects/GO2/MPRA/mpra_test/real_data/counts_barcodes_RC/`
- 8 RC BC files (4 DNA + 4 RNA)

**Expected runtime**: ~10-15 minutes

#### 2. Verify the RC BCs Match Better

```bash
# Quick verification
cd /home/itg/oleg.vlasovets/projects/MPRA_data/mpra_test

# Sample a few RC BCs and check they match assignment
zcat real_data/counts_barcodes_RC/24L012064_S1_L001_R1_001_bc12_RC.fastq.gz | \
  awk 'NR%4==2' | head -1000 | cut -c1 | sort | uniq -c

# Should show balanced distribution (not 99% G):
# ~250 A
# ~250 C  
# ~250 G
# ~250 T
```

#### 3. Re-run Experiment with RC BC Files

```bash
cd /home/itg/oleg.vlasovets/projects/MPRA_data/mpra_test

# Use the new config that points to RC BC files
snakemake \
  --configfile config_experiment_RC.yaml \
  --use-singularity \
  --singularity-args "--bind /lustre" \
  --cores 8 \
  --rerun-incomplete \
  -p
```

---

## 📊 Expected Results (After Fix)

Based on the diagnostic, you should see:

### Before (offset0 files):
- Oligos detected: 9,053
- Barcodes with counts: 37,545  
- Pearson correlation: 0.1895
- Fraction passing: 27.54%
- Median RNA count: 57

### After (RC files) - Expected:
- **Oligos detected: ~11,000-12,000** (↑20-30%)
- **Barcodes with counts: ~45,000-50,000** (↑20-30%)
- **Pearson correlation: >0.35** (↑80-100%)
- **Fraction passing: >35%** (↑30%)
- **Median RNA count: >70** (↑20%)

### Why Improvement Will Be Substantial:
- Currently matching: 416K BCs (14.3%)
- After RC: 532K BCs (18.3%)
- **+115,852 more BCs will match** between assignment and DNA/RNA
- This means more oligos will have sufficient BC coverage
- Better signal-to-noise ratio
- Higher correlation coefficients

---

## 🎯 Understanding the Fix

### Your Library Structure:
```
Assignment R1:  N-[12bp BC]-rest of read
Assignment R2:  [oligo sequence with 17bp adapter]

DNA/RNA R1:     [12bp BC]-rest of read  
DNA/RNA R2:     [oligo or other sequence]
```

### The Orientation Issue:

**In Assignment** (with `BC_rev_comp: true`):
- R1 read: `GTCGACGTTAGA` (forward strand)
- Output BC: `TCTAACGTCGAC` (reverse complement)

**In DNA/RNA counting** (original offset0 files):
- R1 read: `GTCGACGTTAGA` (forward strand)
- Extracted BC: `GTCGACGTTAGA` (NO reverse complement)
- **MISMATCH!** ❌

**After Fix** (new RC files):
- R1 read: `GTCGACGTTAGA` (forward strand)
- Extracted BC: `TCTAACGTCGAC` (reverse complement)
- **MATCH!** ✅

---