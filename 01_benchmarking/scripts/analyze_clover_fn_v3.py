import pandas as pd
import re

def extract_sample_id(x):
    match = re.search(r'(\d+)', str(x))
    if match: return f"simulated_ecDNA_{match.group(1)}"
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
    l1, l2 = end1 - start1, end2 - start2
    return (overlap_len / l1 >= threshold) and (overlap_len / l2 >= threshold)

def get_fn_list(truth_df, pred_df, threshold):
    fn_list = []
    for idx, truth in truth_df.iterrows():
        matched = False
        relevant = pred_df[(pred_df['Sample'] == truth['Sample']) & 
                           (pred_df['Chrom'] == truth['Chrom']) & 
                           (pred_df['Type'] == truth['Type'])]
        for _, pred in relevant.iterrows():
            if check_reciprocal_overlap(truth['Start'], truth['End'], pred['Start'], pred['End'], threshold):
                matched = True
                break
        if not matched:
            fn_list.append(truth.to_dict())
    return fn_list

def main():
    mask_df = pd.read_csv("Mask_Table.csv")
    t2_raw = pd.read_csv("table.tsv", sep="\t")
    truth_df = pd.DataFrame({
        'ID': range(len(t2_raw)), 'Sample': t2_raw['file'].apply(extract_sample_id), 
        'Chrom': t2_raw['chrm'], 'Start': t2_raw['start'], 'End': t2_raw['stop'], 
        'Type': t2_raw['type'].apply(standardize_type), 'Signal': t2_raw['signal']
    })
    truth_df['Length'] = truth_df['End'] - truth_df['Start']

    clover_raw = pd.read_csv("CLOverCNV.tsv", sep="\t")
    clover_df = pd.DataFrame({'Sample': clover_raw['sample_name'].apply(extract_sample_id), 'Chrom': clover_raw['chrom'], 'Start': clover_raw['start'], 'End': clover_raw['end'], 'Type': clover_raw['direction'].apply(standardize_type)})

    print("Comparing FN between 0.5 and 0.1 thresholds...")
    fn_strict = get_fn_list(truth_df, clover_df, 0.5)
    fn_relaxed = get_fn_list(truth_df, clover_df, 0.1)

    strict_ids = set(f['ID'] for f in fn_strict)
    relaxed_ids = set(f['ID'] for f in fn_relaxed)
    recovered_ids = strict_ids - relaxed_ids

    print(f"\n[Summary]")
    print(f"FNs at 0.5 threshold: {len(fn_strict)}")
    print(f"FNs at 0.1 threshold: {len(fn_relaxed)}")
    print(f"Variants 'Recovered' by lowering threshold: {len(recovered_ids)}")

    # Analyze Persistent FNs (the ones still missing at 0.1)
    if fn_relaxed:
        fn_df = pd.DataFrame(fn_relaxed)
        print("\n[Persistent FN Analysis (at 0.1 threshold)]")
        print(">>> By Type:")
        print(fn_df['Type'].value_counts())
        
        bins = [0, 1000, 2000, 5000, float('inf')]
        labels = ['<1kb', '1-2kb', '2-5kb', '>5kb']
        fn_df['LenCat'] = pd.cut(fn_df['Length'], bins=bins, labels=labels)
        print("\n>>> By Length & Type:")
        print(fn_df.groupby(['Type', 'LenCat']).size().unstack(fill_value=0))

        fn_df.to_csv("CLOverCNV_Persistent_FN_0.1.csv", index=False)
        print("\nDetails saved to CLOverCNV_Persistent_FN_0.1.csv")

if __name__ == "__main__":
    main()
