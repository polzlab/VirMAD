#!/usr/bin/env bash
#SBATCH --job-name=clust_vib
#SBATCH --cpus-per-task=20
#SBATCH -N 1
#SBATCH --mem=50G
#SBATCH --time=2-00:00:00
#SBATCH --array=1
#SBATCH --constraint=avx2

organism="vibrio" # {vibrio, phages}
nt=20 # number of cpus
i=${SLURM_ARRAY_TASK_ID}
# IDS=("0" "30" "50" "70" "80" "90" "95" "99" "99.9")
IDS=("99" "99.9")
ID=${IDS[$i]}
$HOME/mm/bin/bash
pcd="../../output/protclu"
prots_path="$pcd/$organism/prots.faa"
clu_path="$TMPDIR/no_partial_$ID.clu"
final_path="$pcd/$organism/prots_$ID.clu"

diamond deepclust \
  --approx-id $ID \
  --member-cover 60 \
  --masking none \
  --motif-masking 0 \
  --evalue 10 \
  --db "$prots_path" \
  --out "${clu_path}_pre" \
  --threads $nt \
  --tmpdir "$TMPDIR"

diamond recluster \
  --approx-id $ID \
  --member-cover 60 \
  --masking none \
  --motif-masking 0 \
  --evalue 10 \
  --db "$prots_path" \
  --clusters "${clu_path}_pre" \
  --out "${final_path}" \
  --threads $nt \
  --tmpdir "$TMPDIR"
rm "${clu_path}_pre"
