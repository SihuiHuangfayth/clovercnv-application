import sys
import subprocess
import os

# 1. 接收 SLURM 传进来的 3 个参数
if len(sys.argv) != 4:
    print("Usage: python check_mimic_signal.py <original_bam> <simulated_bam> <ground_truth_txt>")
    sys.exit(1)

ORIGINAL_BAM = sys.argv[1]
SIMULATED_BAM = sys.argv[2]
GROUND_TRUTH = sys.argv[3]

# 设定宽容的信号阈值（实际 Fold Change 达到目标 CN 的 30% 即认为有足够信号）
SIGNAL_RELAX_RATIO = 0.3 

def get_avg_depth(bam_file, chrom, start, end):
    """调用 samtools 获取区域平均深度"""
    region = f"{chrom}:{start}-{end}"
    cmd = f"samtools depth -r {region} {bam_file} | awk '{{sum+=$3}} END {{if (NR>0) print sum/NR; else print 0}}'"
    try:
        output = subprocess.check_output(cmd, shell=True, text=True).strip()
        return float(output) if output else 0.0
    except Exception as e:
        print(f"Error calculating depth for {region}: {e}")
        return 0.0

def main():
    print(f"\n============================================")
    print(f"自动化质检 (QC Check): {os.path.basename(SIMULATED_BAM)}")
    print(f"============================================")
    
    failed_events = 0
    total_highcopy_events = 0

    print(f"{'Region':<25} | {'Target CN':<10} | {'Orig Depth':<12} | {'Sim Depth':<12} | {'Fold Change':<12} | {'Status'}")
    print("-" * 105)

    with open(GROUND_TRUTH, 'r') as f:
        for line in f:
            parts = line.strip().split()
            # 识别真值文件中的 DUP 事件
            if len(parts) >= 5 and parts[3] == "DUP":
                target_cn = int(parts[4])
                # 只质检 ecDNA 高拷贝事件 (CN > 1)
                if target_cn > 1:
                    total_highcopy_events += 1
                    chrom, start, end = parts[0], parts[1], parts[2]
                    
                    orig_depth = get_avg_depth(ORIGINAL_BAM, chrom, start, end)
                    sim_depth = get_avg_depth(SIMULATED_BAM, chrom, start, end)
                    
                    # 计算 Fold Change
                    fold_change = sim_depth / orig_depth if orig_depth > 0 else sim_depth
                    
                    # 判断逻辑：倍数达标，且绝对深度有明显增加 (>10x)
                    is_passed = fold_change >= (target_cn * SIGNAL_RELAX_RATIO) and (sim_depth - orig_depth > 10)
                    status = "PASS" if is_passed else "FAIL"
                    
                    if not is_passed:
                        failed_events += 1
                        
                    print(f"{chrom}:{start}-{end:<10} | {target_cn:<10} | {orig_depth:<12.2f} | {sim_depth:<12.2f} | {fold_change:<12.2f} | {status}")

    print("-" * 105)
    
    # 边缘情况：如果这个样本本来就没有随机到高拷贝变异
    if total_highcopy_events == 0:
        print("QC Note: No HighCopy events (CN > 1) in this sample. Skipping depth check.")
        sys.exit(0)

    # 最终裁决
    if failed_events > 0:
        print(f"QC FAILED: {failed_events}/{total_highcopy_events} events did not produce sufficient signal.")
        sys.exit(1) # 触发 SLURM 的 set -e，终止任务
    else:
        print(f"QC PASSED: All {total_highcopy_events} ecDNA events successfully produced signals in the BAM!")
        sys.exit(0)

if __name__ == "__main__":
    main()
