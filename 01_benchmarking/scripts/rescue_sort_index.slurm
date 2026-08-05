#!/bin/bash
#SBATCH --job-name=rescue_bams
#SBATCH --output=rescue_bams.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=05:00:00
#SBATCH --partition=standard

set -e
module load samtools

BAM_DIR="/home/gaz6bs/mimic_data/final_bams_V2"
cd $BAM_DIR

echo "开始批量排序和建索引..."

for BAM in simulated_ecDNA_*.bam; do
    # 跳过已经包含 .sorted.bam 的文件
    if [[ "$BAM" == *".sorted.bam"* ]]; then
        continue
    fi

    SORTED_BAM="${BAM%.bam}.sorted.bam"
    
    # 检查是否已经存在对应的 sorted 文件
    if [ ! -f "$SORTED_BAM" ]; then
        echo "处理文件: $BAM -> $SORTED_BAM"
        samtools sort -@ 4 -o "$SORTED_BAM" "$BAM"
        samtools index "$SORTED_BAM"
        
        # 排序完成后，删除原始未排序的 BAM 以节省集群存储空间
        rm "$BAM"
    fi
done

echo "所有 BAM 文件排序和索引完成！"
