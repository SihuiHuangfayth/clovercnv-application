import pandas as pd
import glob
import os
import re

def extract_sample_id(x):
    match = re.search(r'(\d+)', str(x))
    if match:
        return f"simulated_ecDNA_{match.group(1)}"
    return str(x)

def standardize_type(t):
    t = str(t).lower()
    if 'dup' in t or 'gain' in t: return 'DUP'
    if 'del' in t or 'loss' in t: return 'DEL'
    return 'UNKNOWN'

def check_reciprocal_overlap(start1, end1, start2, end2, threshold):
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    overlap_len = overlap_end - overlap_start
    if overlap_len <= 0: return False

    len1, len2 = end1 - start1, end2 - start2
    if len1 <= 0 or len2 <= 0: return False

    ratio1 = overlap_len / len1
    ratio2 = overlap_len / len2
    return (ratio1 >= threshold) and (ratio2 >= threshold)

def filter_masked_predictions(pred_df, mask_df):
    filtered_preds = []
    for _, pred in pred_df.iterrows():
        is_masked = False
        relevant_masks = mask_df[(mask_df['Sample'] == pred['Sample']) & 
                                 (mask_df['Chrom'] == pred['Chrom'])]
        for _, mask in relevant_masks.iterrows():
            if min(pred['End'], mask['End']) > max(pred['Start'], mask['Start']):
                is_masked = True
                break
        if not is_masked:
            filtered_preds.append(pred)
    return pd.DataFrame(filtered_preds)

def evaluate_tool(pred_df, truth_df, threshold):
    tp = fp = 0
    truth_matched = set()
    if len(pred_df) > 0:
        for _, pred in pred_df.iterrows():
            matched = False
            relevant_truths = truth_df[(truth_df['Sample'] == pred['Sample']) &
                                       (truth_df['Chrom'] == pred['Chrom']) &
                                       (truth_df['Type'] == pred['Type'])]
            for idx, truth in relevant_truths.iterrows():
                if check_reciprocal_overlap(pred['Start'], pred['End'], truth['Start'], truth['End'], threshold):
                    matched = True
                    truth_matched.add(idx)
                    break
            if matched: tp += 1
            else: fp += 1
    fn = len(truth_df) - len(truth_matched)
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    return {'TP': tp, 'FP': fp, 'FN': fn, 'Precision': precision, 'Recall': recall, 'F1_Score': f1}

def main():
    print("1. Loading Mask Table...")
    mask_df = pd.read_csv("Mask_Table.csv")

    print("2. Preparing Data...")
    t2_raw = pd.read_csv("table.tsv", sep="\t")
    truth_df = pd.DataFrame({
        'Sample': t2_raw['file'].apply(extract_sample_id), 'Chrom': t2_raw['chrm'],
        'Start': t2_raw['start'], 'End': t2_raw['stop'], 'Type': t2_raw['type'].apply(standardize_type)
    })

    # Load Tools
    clover_raw = pd.read_csv("CLOverCNV.tsv", sep="\t")
    clover_df = pd.DataFrame({'Sample': clover_raw['sample_name'].apply(extract_sample_id), 'Chrom': clover_raw['chrom'], 'Start': clover_raw['start'], 'End': clover_raw['end'], 'Type': clover_raw['direction'].apply(standardize_type)})
    
    cnv_raw = pd.read_csv("CNVnator.tsv", sep="\t")
    cnv_df = pd.DataFrame({'Sample': cnv_raw['file'].apply(extract_sample_id), 'Chrom': cnv_raw['chrom'], 'Start': cnv_raw['start'], 'End': cnv_raw['end'], 'Type': cnv_raw['type'].apply(standardize_type)})
    
    exome_raw = pd.read_csv("exomedepth.tsv", sep="\t")
    exome_df = pd.DataFrame({'Sample': exome_raw['sample_name'].apply(extract_sample_id), 'Chrom': exome_raw['chromosome'], 'Start': exome_raw['start'], 'End': exome_raw['end'], 'Type': exome_raw['type'].apply(standardize_type)})

    tools = {'CLOverCNV': clover_df, 'CNVnator': cnv_df, 'ExomeDepth': exome_df}
    thresholds = [0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.01, 0.001]
    all_results = []

    print("\n3. Running Multi-Threshold Evaluation...")
    for t in thresholds:
        for name, df in tools.items():
            filtered = filter_masked_predictions(df, mask_df)
            metrics = evaluate_tool(filtered, truth_df, t)
            metrics['Tool'] = name
            metrics['Threshold'] = t
            all_results.append(metrics)

    report_df = pd.DataFrame(all_results)
    display_table = report_df.pivot_table(
        index='Tool',
        columns='Threshold',
        values=['Recall', 'F1_Score']
    )

    print("\n" + "="*80)
    print("Threshold test (Recall & F1_Score)")
    print(display_table.round(4))
    print("="*80)

    report_df.to_csv("Multi_Threshold_Report_Full.csv", index=False)
    print("\nFull results saved to: Multi_Threshold_Report_Full.csv")

if __name__ == "__main__":
    main()
