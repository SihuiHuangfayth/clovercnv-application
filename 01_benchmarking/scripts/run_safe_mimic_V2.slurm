#!/bin/bash
#SBATCH --job-name=ecDNA_safe
#SBATCH --output=ecDNA_safe_%A_%a.log
#SBATCH --error=ecDNA_safe_%A_%a.err
#SBATCH --array=51-100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=30:00:00
#SBATCH --partition=standard

set -e
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate bamsurgeon

ID=${SLURM_ARRAY_TASK_ID}

LIST_DIR="../ground_truth_V2"
OUT_DIR="../final_bams_V2"
TMP_DIR="../tmp_mimic_${ID}"

mkdir -p $LIST_DIR $OUT_DIR $TMP_DIR

LIST_FILE="${LIST_DIR}/safe_list_${ID}.txt"
ORIG_BAM="/home/gaz6bs/data/row_Dd2/D2B.clean.bam" 
REF_FASTA="/home/gaz6bs/data/row_Dd2/Dd2.fasta"

echo "============================================"
echo "Sample: ${ID} Begin"

python generate_safe_sv_list.py $LIST_FILE

rm -rf $TMP_DIR
mkdir -p $TMP_DIR

echo "Simulating..."
addsv.py \
  -p 4 \
  -v $LIST_FILE \
  -f $ORIG_BAM \
  -r $REF_FASTA \
  -o $OUT_DIR/simulated_ecDNA_${ID}.bam \
  --tmpdir $TMP_DIR \
  --maxdepth 50000

echo "Build Index..."
samtools index $OUT_DIR/simulated_ecDNA_${ID}.bam

echo "AUTO - QC Check..."

python check_mimic_signal.py $ORIG_BAM $OUT_DIR/simulated_ecDNA_${ID}.bam $LIST_FILE

echo "Sample ${ID} Passed，All done！"
