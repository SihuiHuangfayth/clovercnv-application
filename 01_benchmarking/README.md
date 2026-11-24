# Benchmarking (BAMSurgeon-based simulation)

Benchmarking of CLOverCNV against CNVnator and ExomeDepth using synthetic CNV/ecDNA spiked into a *P. falciparum* Dd2 baseline BAM. This section contributes to **Figure 2 (Panel F)** of the CLOver manuscript.

---

## Quick facts

- **Baseline**: `/home/gaz6bs/data/row_Dd2/D2B.clean.bam` (Dd2 strain, pre-characterized with 5 endogenous CNVs)
- **Reference**: `/home/gaz6bs/data/row_Dd2/Dd2.fasta`
- **Spike-in tool**: BAMSurgeon (with `--maxdepth 50000` to accommodate ecDNA-level fold changes)
- **Simulated samples requested**: 100 (SLURM array 51–100 in V2 + earlier V1 runs)
- **Simulated samples surviving QC**: **18** (numbered 13, 14, 17, 19, 29, 37, 39, 48, 50, 51, 58, 62, 64, 66, 71, 78, 86, 87)
- **Truth-set variants surviving spike + verification**: **69 synthetic + 5 endogenous = 74 total** (reverse-inferred from V3 evaluation: TP+FN=74 for CLOverCNV)
- **Compared tools**: CLOverCNV vs CNVnator vs ExomeDepth
- **Truth-set composition**: 70% standard DUP/DEL (CN=2 or CN=0) + 30% ecDNA-like (CN=5-20)

---

## Data flow

```
[Stage 1] generate_safe_sv_list.py
  ├─ Input: /home/gaz6bs/data/row_Dd2/Dd2.fasta.fai (chromosome sizes)
  └─ Output: safe_list_<ID>.txt (per-sample SV request list, 5 SV/chrom max)
     - 3-50 kb size range
     - 50 kb buffer from chromosome ends
     - 30% of DUPs are ecDNA-like (CN 5-20); 70% are standard (CN=1)
     - "Safe mode": avoids AT-extreme, low-complexity, and repetitive regions
       to prevent BAMSurgeon local reassembly failures ("silent drop")

     ↓

[Stage 2] BAMSurgeon addsv.py (via run_safe_mimic_V2.slurm)
  ├─ Input: safe_list_<ID>.txt + D2B.clean.bam + Dd2.fasta
  └─ Output: simulated_ecDNA_<ID>.bam (in final_bams_V2/)
     - Local reassembly-based synthetic reads
     - Learns error profile + fragment size from baseline
     - --maxdepth 50000 (modified from default to accommodate CN=20 events)

     ↓

[Stage 3] check_mimic_signal.py (silent-drop QC gate)
  ├─ Compares avg depth (original BAM vs simulated BAM) at each ecDNA locus
  ├─ Threshold: sim_depth / orig_depth ≥ target_CN × 0.3 AND depth_delta > 10x
  └─ Fails SLURM job (via set -e) if any ecDNA event lacks sufficient signal
     Result: only samples with all ecDNA events successfully spiked reach evaluation

     ↓

[Stage 4] Truth-set curation (manual/scripted, produces table.tsv + Mask_Table.csv)
  ├─ table.tsv: verified true variants (post-BAMSurgeon; only survivors)
  ├─ background.tsv: 5 endogenous Dd2 CNVs (natively present in baseline)
  └─ Mask_Table.csv: per-sample repetitive-region masking (suppresses FP)

     ↓

[Stage 5] Tool runs (external — outputs collected as *.tsv)
  ├─ CLOverCNV.tsv (58 rows of calls across the 18 simulated samples)
  ├─ CNVnator.tsv (many rows — high FP baseline)
  └─ exomedepth.tsv

     ↓

[Stage 6] benchmark_eval_v3.py
  ├─ Inputs: CLOverCNV.tsv + CNVnator.tsv + exomedepth.tsv + table.tsv + Mask_Table.csv
  ├─ Standardizes: DUP/DEL vocabulary, sample ID format
  ├─ Applies mask filter to predictions
  ├─ Reciprocal overlap threshold sweep: 0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.01, 0.001
  └─ Output: Multi_Threshold_Report_Full.csv (per-tool per-threshold TP/FP/FN/Prec/Rec/F1)

     ↓

[Stage 7] analyze_clover_fn_v3.py (post-hoc FN diagnostics)
  ├─ Compares FN sets between threshold=0.5 and threshold=0.1
  └─ Output: CLOverCNV_Persistent_FN_0.1.csv (FN persistent even at loose threshold, by type × length)
```

---

## Key parameter choices

### Why `--maxdepth 50000`?
BAMSurgeon's default (`--maxdepth 8000`) caps at ~8000× depth, insufficient for ecDNA events at CN=20 in a 60× baseline (which would produce depth ≥ 1200× × 20 = 24000×). Setting `--maxdepth 50000` gives headroom while remaining well below the practical limit.

**Note discrepancy with manuscript methods**: The methods draft states `--maxdepth 20000`. The V2 SLURM script uses `50000`. Reconcile before submission — either update methods to match V2 or verify which was actually used for the paper V3 evaluation.

### Why `EC_DNA_RATIO = 0.3`?
Reflects the expected stochastic distribution of ecDNA events in parasite populations under drug selection. 30% ecDNA-like (high-copy) + 70% standard CNV in the truth set is intended to represent a realistic mix, not pure ecDNA enrichment.

### Why `SIGNAL_RELAX_RATIO = 0.3` in check_mimic_signal.py?
BAMSurgeon does not always achieve exact target CN due to local reassembly variance. Accepting events where measured depth ≥ 30% of target CN allows for typical BAMSurgeon "undershoot" while still catching cases where the spike failed entirely (silent drop).

### Why `VARIANTS_PER_CHROM = 5`?
Lower density reduces the probability of accidental overlap between simulated SVs (which would corrupt the truth set), and allows BAMSurgeon more per-variant compute headroom, reducing silent-drop rate.

### Why reciprocal overlap threshold 0.1 as primary criterion?
See Multi_Threshold_Report_Full.csv: CLOverCNV F1 climbs from 0.338 (threshold 0.5) → 0.769 (threshold 0.1) → 0.785 (threshold 0.01, plateau). The plateau at ≤ 0.1 indicates the tool's segmentation-boundary logic gives calls that overlap the true region but don't achieve tight reciprocal overlap at strict thresholds. 0.1 captures the tool's true sensitivity while remaining defensible.

---

## Files in this directory

### Scripts (`scripts/`)
| File | Stage | Purpose |
|---|---|---|
| `generate_safe_sv_list.py` | 1 | Generate per-sample SV request list with safe-mode constraints |
| `run_safe_mimic_V2.slurm` | 2 | SLURM array job: BAMSurgeon spike + auto-QC |
| `check_mimic_signal.py` | 3 | Silent-drop verification (depth fold-change check) |
| `rescue_sort_index.slurm` | 2 (patch) | Rescue sort/index for BAMs that failed initial finalize step |
| `benchmark_eval_v3.py` | 6 | Multi-threshold tool evaluation → F1/Precision/Recall |
| `analyze_clover_fn_v3.py` | 7 | Persistent false-negative diagnostics for CLOverCNV |

### Truth data (`truth/`)
| File | Contents |
|---|---|
| `table.tsv` | Verified true positive variants (post-spike; only survivors of BAMSurgeon + check_mimic) |
| `background.tsv` | 5 endogenous Dd2 CNVs (natively present in baseline BAM) |
| `Mask_Table.csv` | Per-sample repetitive-region mask (used to filter tool predictions in evaluation) |

### Tool outputs (`tool_outputs/`)
| File | Contents |
|---|---|
| `CLOverCNV.tsv` | 58 raw CLOverCNV calls across 18 samples |
| `CNVnator.tsv` | Raw CNVnator calls (~9700 predictions — many FPs) |
| `exomedepth.tsv` | Raw ExomeDepth calls |

### Evaluation results (`results/`)
| File | Contents |
|---|---|
| `Multi_Threshold_Report_Full.csv` | Full multi-threshold report (9 thresholds × 3 tools = 27 rows) |
| `Multi_Threshold_Report.csv` | Reduced report (5 thresholds) |
| `CLOverCNV_Persistent_FN_0.1.csv` | FNs still missing at threshold=0.1, stratified by type × length |

---

## Headline results (V3, from Multi_Threshold_Report_Full.csv)

At reciprocal overlap threshold = 0.1 (primary criterion):

| Tool | TP | FP | FN | Precision | Recall | F1 |
|---|---|---|---|---|---|---|
| **CLOverCNV** | **50** | **6** | 24 | **0.893** | 0.676 | **0.769** |
| CNVnator | 41 | 9673 | 33 | 0.004 | 0.554 | 0.008 |
| ExomeDepth | 53 | 14 | 21 | 0.791 | 0.716 | 0.752 |

**Interpretation notes for manuscript**:
- CLOverCNV has the **highest precision (0.89)** — 2× fewer FPs than ExomeDepth
- ExomeDepth is a close F1 competitor (0.75 vs 0.77) — the differentiator is precision, not F1
- CNVnator produces ~9700 FPs — two orders of magnitude worse. Not a viable tool for *P. falciparum* CNV calling with default parameters.
- **Not yet stratified**: performance on ecDNA-like (CN>1) subset vs standard CNVs. Stratification likely reveals CLOverCNV's unique advantage. Consider adding this analysis before submission.

Plateau behavior: CLOverCNV F1 plateaus at threshold ≤ 0.05 (F1=0.785); ExomeDepth continues improving marginally to 0.780 at threshold=0.001.

---

## Silent drop: methodological transparency

BAMSurgeon fails silently in low-complexity regions of the AT-rich *P. falciparum* genome, dropping variants without raising errors. The pipeline addresses this by:

1. **Safe-mode generation** (Stage 1) — avoids regions likely to fail
2. **Post-simulation depth verification** (Stage 3) — via `check_mimic_signal.py`, which fails the SLURM job if any ecDNA event lacks expected fold-change signal
3. **Empirical survival rate** — of 100 requested samples, 18 survived full QC (survivals numbered 13-87). Per-sample, requesting ~14 SVs typically yielded ~4 spiked survivors. Total truth-set count: 74 (69 synthetic + 5 endogenous).

**For manuscript methods**: Only 18 samples because BAMSurgeon's silent-drop rate on this genome required post-hoc verification, and we retained only variants whose spiking was empirically confirmed.

---

## Open questions

- [ ] **Stratify results by variant class**: F1/Precision/Recall separately for standard CNV (CN=2, CN=0) vs ecDNA-like (CN≥5). Truth set has both; current eval collapses them.
- [ ] **Add per-tool ecDNA-specific classification metric**: CLOverCNV outputs a `direction` and support-tier column that ExomeDepth/CNVnator lack — this is CLOverCNV's differentiating capability but is not exercised in the current F1 comparison.

---

## Data paths (Rivanna)

```
/home/gaz6bs/archive/mimic_data/
├── ground_truth_V2/         # safe_list_<ID>.txt files (per-sample SV requests)
├── final_bams_V2/           # BAMSurgeon output BAMs (simulated_ecDNA_<ID>.bam)
├── benchmarking/version_3/  # V3 evaluation (source of files in results/)
├── benchmarking/CLOverCNV/  # CLOverCNV raw output per sample (segments.calls.tsv.gz)
└── benchmarking/cnvnator/   # CNVnator raw output
```

Large data files (BAMs, ground-truth intermediate files) are NOT copied to this repo per .gitignore. Rivanna remains the source of truth for raw data.
