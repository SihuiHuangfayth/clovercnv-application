#!/usr/bin/env python
"""
Aggregate CLOverCNV keep_suggested=TRUE calls across all 64 samples.

Output:
  - all_calls.tsv: every keep_suggested=TRUE call from every sample
  - per_sample_summary.tsv: N calls per sample by direction/tier
"""

import pandas as pd
import glob
import os
import gzip

INPUT_DIR = "/scratch/gaz6bs/clovercnv/"
OUTPUT_DIR = "/scratch/gaz6bs/clovercnv/analysis/"

# Find all final2 files (handles the nested SRR/SRR/ structure)
final2_files = sorted(glob.glob(f"{INPUT_DIR}/SRR*/SRR*/*.segments.calls.final2.tsv.gz"))
print(f"Found {len(final2_files)} final2 files")

all_calls = []
per_sample_stats = []

for fp in final2_files:
    srr = os.path.basename(fp).split('.')[0]
    df = pd.read_csv(fp, sep='\t', compression='gzip')
    
    kept = df[df['keep_suggested'] == True].copy()
    
    # Aggregate stats regardless of whether kept is empty
    stats = {
        'srr': srr,
        'total_segments': len(df),
        'kept_calls': len(kept),
        'gain_calls': (kept['direction'] == 'gain').sum() if len(kept) else 0,
        'loss_calls': (kept['direction'] == 'loss').sum() if len(kept) else 0,
        'strong_effect': kept['call'].str.contains('strongEffect', na=False).sum() if len(kept) else 0,
        'weak_effect': kept['call'].str.contains('weakEffect', na=False).sum() if len(kept) else 0,
        'max_ratio': kept['ratio'].max() if len(kept) else 0,
        'high_cn_calls_ge5': (kept['ratio'] >= 5).sum() if len(kept) else 0,   # ecDNA-candidate calls
    }
    per_sample_stats.append(stats)
    
    if len(kept):
        all_calls.append(kept)

# Write outputs
all_calls_df = pd.concat(all_calls, ignore_index=True) if all_calls else pd.DataFrame()
summary_df = pd.DataFrame(per_sample_stats)

os.makedirs(OUTPUT_DIR, exist_ok=True)
all_calls_df.to_csv(f"{OUTPUT_DIR}/all_kept_calls.tsv", sep='\t', index=False)
summary_df.to_csv(f"{OUTPUT_DIR}/per_sample_summary.tsv", sep='\t', index=False)

# Console summary
print(f"\nTotal kept calls across all samples: {len(all_calls_df)}")
print(f"Samples with any kept call: {(summary_df['kept_calls'] > 0).sum()} / {len(summary_df)}")
print(f"Samples with high-CN calls (ratio >= 5): {(summary_df['high_cn_calls_ge5'] > 0).sum()}")
print(f"\nTop 10 samples by kept_calls:")
print(summary_df.nlargest(10, 'kept_calls')[['srr', 'kept_calls', 'gain_calls', 'loss_calls', 'high_cn_calls_ge5', 'max_ratio']].to_string(index=False))

print(f"\n\nOutputs written to {OUTPUT_DIR}")
