# Barcode Overlap Analysis (Raw Counts, Before Assignment)

**Source:** `results/experiments/examplenoMinCount5/statistic/bc_overlap.counts.noMinCount.tsv`

Max requested: "It would be great to see the overlap of barcodes for DNA and RNA separately *before* assignment."

**Note:** Raw barcode overlap statistics are identical across all filtering configurations (noMinCount1, noMinCount3, noMinCount5) because they are calculated from the raw barcode files *before* any filtering is applied. The filtering parameters (bc_threshold, min_dna_counts, min_rna_counts) only affect downstream assignment and oligo-level statistics.

---

## DNA Replicates - Barcode Overlap

| Replicate Pair | BCs Rep A | BCs Rep B | BCs Overlap | Overlap % (A) | Overlap % (B) |
|----------------|-----------|-----------|-------------|---------------|---------------|
| **Rep 1 vs 2** | 2,799,319 | 2,920,524 | 2,314,674   | **82.7%**     | **79.3%**     |
| **Rep 1 vs 3** | 2,799,319 | 2,806,296 | 2,242,575   | **80.1%**     | **79.9%**     |
| **Rep 1 vs 4** | 2,799,319 | 2,391,239 | 1,955,630   | **69.9%**     | **81.8%**     |
| **Rep 2 vs 3** | 2,920,524 | 2,806,296 | 2,318,418   | **79.4%**     | **82.6%**     |
| **Rep 2 vs 4** | 2,920,524 | 2,391,239 | 2,016,998   | **69.1%**     | **84.3%**     |
| **Rep 3 vs 4** | 2,806,296 | 2,391,239 | 1,957,965   | **69.8%**     | **81.9%**     |

**Summary:** DNA replicates show **good consistency** with 70-84% barcode overlap across all pairs. Replicate 4 has fewer unique barcodes (~2.4M vs ~2.8M) but still maintains reasonable overlap.

---

## RNA Replicates - Barcode Overlap

| Replicate Pair | BCs Rep A | BCs Rep B | BCs Overlap | Overlap % (A) | Overlap % (B) |
|----------------|-----------|-----------|-------------|---------------|---------------|
| **Rep 1 vs 2** | 524,996   | 443,750   | 142,320     | **27.1%**     | **32.1%**     |
| **Rep 1 vs 3** | 524,996   | 1,168,849 | 270,135     | **51.5%**     | **23.1%**     |
| **Rep 1 vs 4** | 524,996   | 410,043   | 269,035     | **51.2%**     | **65.6%**     |
| **Rep 2 vs 3** | 443,750   | 1,168,849 | 313,424     | **70.6%**     | **26.8%**     |
| **Rep 2 vs 4** | 443,750   | 410,043   | 118,419     | **26.7%**     | **28.9%**     |
| **Rep 3 vs 4** | 1,168,849 | 410,043   | 220,215     | **18.8%**     | **53.7%**     |

**Summary:** RNA replicates show **VERY POOR consistency** with only 19-71% overlap. The asymmetry in overlap percentages is striking:
- Rep 3 has 2-3x more unique barcodes (1.17M) than other replicates (410K-525K)
- No clear pairing pattern (unlike Max observed: Rep 1+4 vs Rep 2+3)
- Suggests different library fractions or treatment conditions

---

## Key Observations

### 1. **DNA vs RNA Barcode Counts**
- **DNA:** ~2.4-2.9 million unique barcodes per replicate
- **RNA:** ~410K-1.2 million unique barcodes per replicate
- **Ratio:** DNA has **2-7x more unique barcodes** than RNA

This confirms Max's point #2: "Many barcodes are not present in the RNA, and your small number of RNA barcodes have very high counts."

### 2. **Replicate Inconsistency**
- **DNA:** Relatively consistent across replicates (70-84% overlap)
- **RNA:** Highly inconsistent (19-71% overlap)
- **Replicate 3 (RNA):** Has far more unique barcodes than other RNA replicates

This supports Max's point #1: "Replicates 1 and 4 match and also slightly replicates 2 and 3, but not consistently across all samples. In particular, the RNA seems to be off."

### 3. **Barcode Dropout**
The poor RNA overlap suggests:
- **Bottleneck/dropout:** Barcodes were lost during library prep or RNA extraction
- **PCR amplification bias:** Only a subset of the library was efficiently amplified
- **Different treatment:** RNA samples may have been processed differently

### 4. **Estimated Total Unique Barcodes**
Using Lincoln-Peterson estimator:
- **DNA:** ~3.48 million total unique barcodes (mean estimate)
- **RNA:** ~1.68 million total unique barcodes (mean estimate)
- **Expected (design):** 17,997 oligos × ~230 BCs/oligo = ~4.1M theoretical maximum

---

## Raw Read Depth Statistics

| Sample | Type | Total Reads | Unique BCs | Reads/BC | BCs Shared (DNA+RNA) |
|--------|------|-------------|------------|----------|----------------------|
| OA.1   | DNA  | 9,622,075   | 2,799,363  | 3.4      | 389,656 (13.9%)      |
| OA.2   | DNA  | 10,074,852  | 2,920,566  | 3.4      | 332,990 (11.4%)      |
| OA.3   | DNA  | 8,849,208   | 2,806,329  | 3.2      | 867,922 (30.9%)      |
| OA.4   | DNA  | 6,073,076   | 2,391,262  | 2.5      | 278,225 (11.6%)      |
| OA.1   | RNA  | 16,410,470  | 524,996    | 31.3     | 389,656 (74.2%)      |
| OA.2   | RNA  | 20,640,344  | 443,750    | 46.5     | 332,990 (75.0%)      |
| OA.3   | RNA  | 35,791,170  | 1,168,849  | 30.6     | 867,922 (74.3%)      |
| OA.4   | RNA  | 9,438,216   | 410,043    | 23.0     | 278,225 (67.9%)      |

**Key findings:**
- **DNA:** 2.5-3.4 reads per barcode (low depth, typical for DNA)
- **RNA:** 23-46.5 reads per barcode (much higher depth, as expected)
- **Sequencing ratio (RNA:DNA):** ~3:1 total reads (appropriate)
- **Barcode ratio (DNA:RNA):** ~3-7:1 (INVERTED - DNA has more BCs!)

---

## Barcode-Level Correlations (Raw Counts)

**DNA Spearman correlations:**
- Rep 1 vs 2: 0.410
- Rep 1 vs 3: 0.408
- Rep 1 vs 4: 0.352
- Rep 2 vs 3: 0.434
- Rep 2 vs 4: 0.352
- Rep 3 vs 4: 0.355

**RNA Spearman correlations:**
- Rep 1 vs 2: 0.259
- Rep 1 vs 3: 0.237
- Rep 1 vs 4: **0.696** ← BEST RNA pair
- Rep 2 vs 3: **0.421** ← Second best RNA pair
- Rep 2 vs 4: 0.253
- Rep 3 vs 4: 0.231

**Observations:**
- DNA correlations are moderate (0.35-0.43) - acceptable for raw barcode counts
- RNA correlations are mostly VERY POOR (0.23-0.26), except:
  - **Rep 1 + Rep 4:** Strong correlation (0.696) - these match well
  - **Rep 2 + Rep 3:** Moderate correlation (0.421) - these match somewhat
- This confirms Max's observation: "Replicates 1 and 4 match and also slightly replicates 2 and 3"

---

## Conclusion

This raw barcode overlap data (consistent across all filtering configurations) confirms Max's diagnosis of experimental failure:
- ✅ DNA replicates are consistent (good technical quality)
- ❌ RNA replicates are highly inconsistent (poor technical quality)
- ❌ Major barcode dropout in RNA compared to DNA
- ❌ RNA barcode diversity varies wildly between replicates (410K to 1.2M)
- ✅ **Replicate pairing exists** (Rep 1+4 and Rep 2+3 cluster), but overall RNA quality is poor

**The experiment appears to have a serious RNA library issue** - either bottleneck, PCR bias, or differential treatment of replicates.

**Important:** These poor overlap statistics appear in the raw data *before* any filtering, so they cannot be improved by adjusting filtering parameters (bc_threshold, min_counts). The problem is fundamental to the experimental library preparation.

---

## File Locations

Raw barcode overlap statistics (identical across all experiments):
```
/lustre/groups/itg/olegv/mpra_analysis/results/experiments/examplenoMinCount*/statistic/
├── bc_overlap.counts.noMinCount.tsv              # Raw BC overlap (this analysis)
├── bc_overlap.assigned_counts.noMinCount.fromFile.tsv  # Assigned BC overlap
├── counts.raw.tsv                                # Read depth statistics
├── statistic_bc_correlation_merged.noMinCount.tsv  # BC-level correlations
└── statistic_oligo_correlation_merged.fromFile.noMinCount.tsv  # Final oligo correlations
```

---

## Final Filtering Results Comparison

**Impact of bc_threshold on final oligo-level results** (all with min_dna_counts=0, min_rna_counts=0):

| Config       | bc_threshold | Correlation | Passing | Median BCs | Median RNA |
|--------------|--------------|-------------|---------|------------|------------|
| noMinCount1  | 1            | **0.3316**  | **97.9%** | 32       | 98         |
| noMinCount3  | 3            | 0.3214      | 93.1%   | 32         | 98         |
| noMinCount5  | 5            | 0.3071      | 88.5%   | 32         | 98         |

**Best result:** noMinCount1 (bc_threshold=1) with 0.3316 correlation

While bc_threshold=1 gives the best correlation, all three show similar patterns in the **raw barcode overlap** (analyzed above), confirming that the poor RNA replicate quality is a fundamental experimental issue, not a filtering artifact.
