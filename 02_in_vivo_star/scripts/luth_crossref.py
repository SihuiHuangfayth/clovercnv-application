#!/usr/bin/env python
"""
luth_crossref.py — Cross-reference CLOverCNV calls with Luth 2024 Data S4 CNV list.

Purpose
-------
Jenny's meeting (2026-08-13) asked: does CLOverCNV output correspond to what Luth 2024 reported?
"50/50 correspondence granularity is sufficient."

This script uses a RATIO-BASED filter (not keep_suggested=TRUE alone), because our validation
showed CLOverCNV's default confidence thresholds are conservative and filter out real biological
signals — e.g., pfabci3 in MMV006901: 8/8 samples show elevated ratio, but only 1/8 passes
keep_suggested=TRUE.

Data source
-----------
Luth 2024 Supplementary Data 4 (NIHMS2048195), filtered to 42 CNVs across our 4 cohorts:
  - Suloctidil (13 CNVs, 10 clones)
  - MMV006901 (14 CNVs, 8 clones)
  - DSM267 (1 CNV, 1 clone: F2-2-E8)
  - MMV019662 (14 CNVs, 13 clones; 1 involves excluded sample SRR5161093)

Loaded from: 02_in_vivo_star/luth_data_s4_filtered.tsv

Approach
--------
For each Luth CNV entry (which is per-clone), check whether the specific matched CLOverCNV
sample has an overlapping segment.
  - Overlap = same chromosome, coordinates intersect.
  - Ratio-based detection: gain if ratio > GAIN_RATIO_THRESHOLD; loss if ratio < LOSS_RATIO_THRESHOLD.
  - Also report keep_suggested status for informational purposes.

Usage
-----
    conda activate CLOverCNV
    python luth_crossref.py

Prerequisites:
    - CLOverCNV outputs at /scratch/gaz6bs/clovercnv/SRR*/SRR*/*.segments.calls.final2.tsv.gz
    - Filtered Luth Data S4 at LUTH_DATA_S4_TSV path below

Output:
    - luth_correspondence_per_cnv.tsv     : one row per Luth CNV × matched sample
    - luth_correspondence_by_compound.tsv : summary per compound cohort
    - console: human-readable summary for email to Jenny
"""

import pandas as pd
import glob
import os
from collections import defaultdict


# ============================================================
# CONFIGURATION
# ============================================================

CLOVER_DIR = "/scratch/gaz6bs/clovercnv"
LUTH_DATA_S4_TSV = "/scratch/gaz6bs/clovercnv-application/02_in_vivo_star/luth_data_s4_filtered.tsv"
OUTPUT_DIR = "/scratch/gaz6bs/clovercnv/analysis"

# Ratio thresholds — relaxed vs default because default keep_suggested is too conservative
# (pfabci3 shows 1.6-2.9x ratio in 8/8 MMV006901 samples but only 1/8 kept).
GAIN_RATIO_THRESHOLD = 1.5
LOSS_RATIO_THRESHOLD = 0.7


# ============================================================
# HELPERS
# ============================================================

def load_luth_s4(path):
    """Load pre-filtered Luth Data S4 (42 CNVs across our 4 cohorts)."""
    df = pd.read_csv(path, sep='\t')
    # Skip entries involving excluded samples
    excluded_mask = df['srr'].astype(str).str.startswith('EXCLUDED')
    if excluded_mask.any():
        n_excluded = excluded_mask.sum()
        excluded_clones = df.loc[excluded_mask, 'clone_name'].unique()
        print(f"Note: Skipping {n_excluded} Luth CNV(s) involving excluded samples: {list(excluded_clones)}")
        df = df.loc[~excluded_mask].copy()
    return df


def load_sample_calls(srr):
    """Load CLOverCNV final2 output for one sample. Returns None if not found."""
    pattern = f"{CLOVER_DIR}/{srr}/{srr}/*.segments.calls.final2.tsv.gz"
    matches = glob.glob(pattern)
    if not matches:
        return None
    return pd.read_csv(matches[0], sep='\t', compression='gzip')


def check_cnv_in_sample(luth_cnv, sample_calls):
    """
    For one Luth CNV row and the matched sample's CLOverCNV output,
    check for overlap and ratio signal at that locus.
    """
    chrom = luth_cnv['chrom']
    start = int(luth_cnv['cnv_start'])
    end = int(luth_cnv['cnv_end'])
    direction = luth_cnv['direction']

    # Same-chromosome segments
    same_chr = sample_calls[sample_calls['chrom'] == chrom]
    overlapping = same_chr[(same_chr['end'] > start) & (same_chr['start'] < end)]

    if len(overlapping) == 0:
        return {
            'overlapping_segments': 0,
            'ratio_signal_detected': False,
            'keep_suggested_any': False,
            'max_ratio_in_region': None,
            'best_segment_ratio': None,
            'best_segment_call': None,
            'best_segment_kept': False,
        }

    # Ratio-based detection based on Luth's stated direction
    if direction == 'gain':
        ratio_hits = overlapping[overlapping['ratio'] > GAIN_RATIO_THRESHOLD]
        best_idx = overlapping['ratio'].idxmax()
    elif direction == 'loss':
        ratio_hits = overlapping[overlapping['ratio'] < LOSS_RATIO_THRESHOLD]
        best_idx = overlapping['ratio'].idxmin()
    else:  # 'both' or unknown
        ratio_hits = overlapping[
            (overlapping['ratio'] > GAIN_RATIO_THRESHOLD) |
            (overlapping['ratio'] < LOSS_RATIO_THRESHOLD)
        ]
        best_idx = overlapping['ratio'].idxmax()

    best = overlapping.loc[best_idx]
    return {
        'overlapping_segments': len(overlapping),
        'ratio_signal_detected': len(ratio_hits) > 0,
        'keep_suggested_any': bool((overlapping['keep_suggested'] == True).any()),
        'max_ratio_in_region': float(overlapping['ratio'].max()),
        'best_segment_ratio': float(best['ratio']),
        'best_segment_call': str(best.get('call', 'N/A')),
        'best_segment_kept': bool(best.get('keep_suggested', False)),
    }


# ============================================================
# MAIN
# ============================================================

def main():
    luth = load_luth_s4(LUTH_DATA_S4_TSV)
    print(f"Loaded {len(luth)} Luth CNVs to check (post-exclusion).")

    results = []
    for _, cnv in luth.iterrows():
        srr = cnv['srr']
        calls = load_sample_calls(srr)
        if calls is None:
            print(f"⚠ No CLOverCNV output found for {srr} ({cnv['clone_name']})")
            continue
        r = check_cnv_in_sample(cnv, calls)
        results.append({
            'srr': srr,
            'clone_name': cnv['clone_name'],
            'compound': cnv['compound'],
            'chrom': cnv['chrom'],
            'cnv_start': int(cnv['cnv_start']),
            'cnv_end': int(cnv['cnv_end']),
            'cnv_size_genes': cnv['cnv_size_genes'],
            'luth_direction': cnv['direction'],
            'luth_avg_ratio': cnv['avg_copy_ratio'],
            'first_gene': str(cnv['gene_list']).split(',')[0] if cnv['gene_list'] else '',
            **r,
        })

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    per_cnv_df = pd.DataFrame(results)
    per_cnv_df.to_csv(f"{OUTPUT_DIR}/luth_correspondence_per_cnv.tsv", sep='\t', index=False)

    # Per-compound summary
    summary_rows = []
    for compound in sorted(per_cnv_df['compound'].unique()):
        sub = per_cnv_df[per_cnv_df['compound'] == compound]
        n = len(sub)
        summary_rows.append({
            'compound': compound,
            'n_luth_cnvs_checked': n,
            'n_ratio_detected': sub['ratio_signal_detected'].sum(),
            'n_keep_suggested_any': sub['keep_suggested_any'].sum(),
            'ratio_recall_rate': round(sub['ratio_signal_detected'].mean(), 3),
            'keep_suggested_recall_rate': round(sub['keep_suggested_any'].mean(), 3),
            'median_max_ratio_in_luth_region': round(sub['max_ratio_in_region'].median(), 3)
                if sub['max_ratio_in_region'].notna().any() else None,
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(f"{OUTPUT_DIR}/luth_correspondence_by_compound.tsv", sep='\t', index=False)

    # Console summary — Jenny-ready
    n_total = len(per_cnv_df)
    n_ratio = per_cnv_df['ratio_signal_detected'].sum()
    n_kept = per_cnv_df['keep_suggested_any'].sum()

    print("\n" + "=" * 70)
    print("Luth 2024 Data S4 correspondence — SUMMARY")
    print("=" * 70)
    print(f"\nTotal Luth CNVs checked: {n_total}")
    print(f"  Ratio-level correspondence (any signal at region): {n_ratio}/{n_total} ({100*n_ratio/n_total:.1f}%)")
    print(f"  Pipeline-level (keep_suggested=TRUE): {n_kept}/{n_total} ({100*n_kept/n_total:.1f}%)")

    print("\n--- Per-compound breakdown ---")
    for _, r in summary_df.iterrows():
        print(f"  {r['compound']:12} n={r['n_luth_cnvs_checked']:3}  "
              f"ratio_recall={r['n_ratio_detected']}/{r['n_luth_cnvs_checked']} "
              f"({r['ratio_recall_rate']*100:.1f}%);  "
              f"keep_suggested={r['n_keep_suggested_any']}/{r['n_luth_cnvs_checked']}")

    # Flag CNVs missed even at ratio level
    missed_ratio = per_cnv_df[~per_cnv_df['ratio_signal_detected']]
    if len(missed_ratio) > 0:
        print(f"\n--- CNVs NOT recovered at ratio level ({len(missed_ratio)}) ---")
        for _, r in missed_ratio.iterrows():
            max_r = r['max_ratio_in_region']
            max_str = f"max_ratio={max_r:.2f}" if pd.notna(max_r) else "no overlap"
            print(f"  {r['srr']} {r['compound']} {r['chrom']} {r['first_gene']}: {max_str}")

    print(f"\nOutputs:")
    print(f"  - {OUTPUT_DIR}/luth_correspondence_per_cnv.tsv  (per-CNV × per-sample detail)")
    print(f"  - {OUTPUT_DIR}/luth_correspondence_by_compound.tsv  (per-compound summary)")


if __name__ == "__main__":
    main()
