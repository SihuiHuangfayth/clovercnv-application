#!/usr/bin/env python
"""
luth_crossref.py — Cross-reference CLOverCNV calls with published CNV list from Luth 2024 Data S4.

Purpose
-------
Jenny's meeting (2026-08-13) asked: does CLOverCNV output correspond to what Luth 2024 reported?
"50/50 correspondence granularity is sufficient."

This script uses a ratio-based filter (NOT keep_suggested=TRUE alone) because our validation
showed CLOverCNV's default confidence thresholds are conservative and filter out real biological
signals (e.g., pfabci3 in MMV006901 — 8/8 samples show elevated ratio but only 1/8 passes keep_suggested).

Approach
--------
1. Load Luth Data S4 CNV list (42 CNVs across 4 cohorts).
2. For each Luth CNV, check each cohort sample's CLOverCNV output for overlapping segments.
3. Use RATIO-BASED filter: ratio > 1.5 (gain) or < 0.7 (loss) — regardless of keep_suggested.
4. Also report keep_suggested status for informational purposes.
5. Cross-reference at GENE-LEVEL (not exact coordinates) — Luth's CNVs are typically gene-annotated.
6. Output: correspondence table + summary statistics.

Usage
-----
Prerequisites:
    - Populate LUTH_DATA_S4 with actual coordinates from Luth 2024 Supplementary Data S4.
    - CLOverCNV outputs at /scratch/gaz6bs/clovercnv/SRR*/SRR*/*.segments.calls.final2.tsv.gz.
    - Sample manifest at ~/repo/02_in_vivo_star/sample_manifest_cohort.tsv (SRR → compound mapping).

Run:
    python luth_crossref.py

Output:
    - luth_correspondence_per_sample.tsv    : per-sample × per-CNV recall matrix
    - luth_correspondence_summary.tsv       : per-CNV recall rate across cohort
    - console: human-readable summary for email to Jenny
"""

import pandas as pd
import glob
import os
from collections import defaultdict


# ============================================================
# CONFIGURATION — populate before running
# ============================================================

CLOVER_DIR = "/scratch/gaz6bs/clovercnv"
MANIFEST_PATH = "/scratch/gaz6bs/clovercnv-application/02_in_vivo_star/sample_manifest_cohort.tsv"
OUTPUT_DIR = "/scratch/gaz6bs/clovercnv/analysis"

# Ratio thresholds (relaxed, per interim finding: pfabci3 is 1.6–2.9× but 7/8 fail default thresholds)
GAIN_RATIO_THRESHOLD = 1.5
LOSS_RATIO_THRESHOLD = 0.7

# Minimum reciprocal overlap to consider "same locus" (loose for gene-level correspondence)
MIN_OVERLAP_BP = 1  # Any overlap counts as "same locus" for coarse correspondence


# ============================================================
# LUTH DATA S4 — TEMPLATE (populate with actual values from paper)
# ============================================================
# Each row: (chrom, start, end, gene_or_locus, cohort, direction, notes)
# Extract from Luth 2024 Supplementary Data S4 (42 CNVs across 4 cohorts studied here).
#
# Example row structure (REPLACE WITH REAL DATA):
LUTH_DATA_S4 = [
    # -------- MMV006901 (known: pfabci3) --------
    # ("Pf3D7_03_v3", 810000, 887000, "pfabci3", "MMV006901", "gain", "Luth Data S4 row X"),
    #
    # -------- Suloctidil --------
    # ("Pf3D7_XX_v3", start, end, "geneX", "Suloctidil", "gain/loss", "..."),
    #
    # -------- DSM267 (known: dhodh, dhfr-ts) --------
    # ("Pf3D7_06_v3", start, end, "dhodh", "DSM267", "gain", "..."),
    # ("Pf3D7_04_v3", start, end, "dhfr-ts", "DSM267", "gain", "Mandt2019 Table 1 also"),
    #
    # -------- MMV019662 (known: pfncr1 as SNV target; CNV is discovery goal) --------
    # ("Pf3D7_01_v3", start, end, "pfncr1", "MMV019662", "gain/loss", "..."),
]


# ============================================================
# HELPERS
# ============================================================

def load_manifest(path):
    """Load sample manifest, return dict: SRR -> compound."""
    df = pd.read_csv(path, sep='\t')
    return dict(zip(df['srr'], df['cohort']))


def load_sample_calls(srr):
    """Load CLOverCNV final2 output for one sample."""
    pattern = f"{CLOVER_DIR}/{srr}/{srr}/*.segments.calls.final2.tsv.gz"
    matches = glob.glob(pattern)
    if not matches:
        return None
    return pd.read_csv(matches[0], sep='\t', compression='gzip')


def check_luth_cnv_in_sample(luth_row, sample_calls):
    """
    For one Luth CNV and one sample's CLOverCNV calls, check for overlap.
    Returns: dict with match info.
    """
    chrom, start, end, gene, cohort, direction, notes = luth_row
    # Filter sample calls to same chromosome
    same_chr = sample_calls[sample_calls['chrom'] == chrom].copy()
    # Find overlapping segments
    overlapping = same_chr[
        (same_chr['end'] > start) & (same_chr['start'] < end)
    ]
    if len(overlapping) == 0:
        return {
            'overlapping_segments': 0,
            'ratio_signal_detected': False,
            'keep_suggested_any': False,
            'max_ratio_in_region': None,
            'best_segment_call': None,
        }
    # Apply ratio-based filter
    if direction == 'gain':
        ratio_hits = overlapping[overlapping['ratio'] > GAIN_RATIO_THRESHOLD]
    elif direction == 'loss':
        ratio_hits = overlapping[overlapping['ratio'] < LOSS_RATIO_THRESHOLD]
    else:  # Unknown direction, check both
        ratio_hits = overlapping[
            (overlapping['ratio'] > GAIN_RATIO_THRESHOLD) |
            (overlapping['ratio'] < LOSS_RATIO_THRESHOLD)
        ]
    kept_any = (overlapping['keep_suggested'] == True).any()
    best_call = overlapping.sort_values('ratio', ascending=(direction == 'loss')).iloc[0]
    return {
        'overlapping_segments': len(overlapping),
        'ratio_signal_detected': len(ratio_hits) > 0,
        'keep_suggested_any': bool(kept_any),
        'max_ratio_in_region': float(best_call['ratio']),
        'best_segment_call': str(best_call.get('call', 'N/A')),
    }


# ============================================================
# MAIN
# ============================================================

def main():
    if not LUTH_DATA_S4:
        print("=" * 70)
        print("⚠ LUTH_DATA_S4 is empty. Populate it before running.")
        print("=" * 70)
        print("\nTo populate:")
        print("  1. Open Luth 2024 Supplementary Data S4 (Excel/PDF).")
        print("  2. For each CNV row in cohorts Suloctidil/MMV006901/DSM267/MMV019662,")
        print("     add a tuple to LUTH_DATA_S4 above:")
        print("     (chrom, start, end, gene_or_locus, cohort, direction, notes)")
        print("  3. Direction should be 'gain' or 'loss' (or 'both' if unclear).")
        print("  4. Save and re-run this script.")
        return

    print(f"Loaded {len(LUTH_DATA_S4)} Luth CNVs to check.")
    manifest = load_manifest(MANIFEST_PATH)
    print(f"Loaded manifest: {len(manifest)} cohort samples.")

    per_sample_results = []
    per_cnv_summary = defaultdict(lambda: {
        'gene': '', 'cohort': '', 'chrom': '', 'start': 0, 'end': 0,
        'n_samples_checked': 0,
        'n_samples_with_overlap': 0,
        'n_samples_ratio_detected': 0,
        'n_samples_kept_suggested': 0,
        'sample_ratios': [],
    })

    for luth_row in LUTH_DATA_S4:
        chrom, start, end, gene, cohort, direction, notes = luth_row
        cnv_key = f"{gene}_{chrom}:{start}-{end}"
        per_cnv_summary[cnv_key].update({
            'gene': gene, 'cohort': cohort, 'chrom': chrom,
            'start': start, 'end': end, 'direction': direction,
        })

        # Find cohort samples for this Luth CNV
        cohort_samples = [srr for srr, comp in manifest.items() if comp == cohort]
        for srr in cohort_samples:
            calls = load_sample_calls(srr)
            if calls is None:
                continue
            result = check_luth_cnv_in_sample(luth_row, calls)
            per_sample_results.append({
                'luth_cnv': cnv_key,
                'gene': gene,
                'cohort': cohort,
                'srr': srr,
                'expected_direction': direction,
                **result,
            })
            summary = per_cnv_summary[cnv_key]
            summary['n_samples_checked'] += 1
            if result['overlapping_segments'] > 0:
                summary['n_samples_with_overlap'] += 1
            if result['ratio_signal_detected']:
                summary['n_samples_ratio_detected'] += 1
            if result['keep_suggested_any']:
                summary['n_samples_kept_suggested'] += 1
            if result['max_ratio_in_region'] is not None:
                summary['sample_ratios'].append(result['max_ratio_in_region'])

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    per_sample_df = pd.DataFrame(per_sample_results)
    per_sample_df.to_csv(f"{OUTPUT_DIR}/luth_correspondence_per_sample.tsv", sep='\t', index=False)

    summary_rows = []
    for cnv_key, s in per_cnv_summary.items():
        n = s['n_samples_checked']
        summary_rows.append({
            'luth_cnv': cnv_key, 'gene': s['gene'], 'cohort': s['cohort'],
            'direction': s['direction'], 'n_samples': n,
            'ratio_recall_rate': (s['n_samples_ratio_detected'] / n) if n else 0,
            'keep_suggested_recall_rate': (s['n_samples_kept_suggested'] / n) if n else 0,
            'ratio_recall_count': f"{s['n_samples_ratio_detected']}/{n}",
            'keep_suggested_recall_count': f"{s['n_samples_kept_suggested']}/{n}",
            'median_max_ratio': (sorted(s['sample_ratios'])[len(s['sample_ratios'])//2]
                                 if s['sample_ratios'] else None),
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(f"{OUTPUT_DIR}/luth_correspondence_summary.tsv", sep='\t', index=False)

    # Console summary for email to Jenny
    print("\n" + "=" * 70)
    print("Luth 2024 Data S4 correspondence — SUMMARY")
    print("=" * 70)
    print(f"\nTotal Luth CNVs checked: {len(per_cnv_summary)}")
    ratio_recovered = sum(1 for s in per_cnv_summary.values()
                          if s['n_samples_ratio_detected'] > 0)
    kept_recovered = sum(1 for s in per_cnv_summary.values()
                         if s['n_samples_kept_suggested'] > 0)
    print(f"Recovered at ratio level (any sample): {ratio_recovered}/{len(per_cnv_summary)}")
    print(f"Recovered at keep_suggested level (any sample): {kept_recovered}/{len(per_cnv_summary)}")

    print("\n--- By cohort ---")
    for cohort in sorted(set(s['cohort'] for s in per_cnv_summary.values())):
        cohort_cnvs = [s for s in per_cnv_summary.values() if s['cohort'] == cohort]
        n = len(cohort_cnvs)
        ratio_recov = sum(1 for s in cohort_cnvs if s['n_samples_ratio_detected'] > 0)
        kept_recov = sum(1 for s in cohort_cnvs if s['n_samples_kept_suggested'] > 0)
        print(f"  {cohort}: {ratio_recov}/{n} at ratio level; {kept_recov}/{n} at keep_suggested level")

    print(f"\nOutputs saved to {OUTPUT_DIR}/")
    print("  - luth_correspondence_per_sample.tsv  (per-sample × per-CNV recall matrix)")
    print("  - luth_correspondence_summary.tsv     (per-CNV recall stats)")


if __name__ == "__main__":
    main()
