#!/usr/bin/env bash
#SBATCH --job-name=defensefinder
#SBATCH --cpus-per-task=5
#SBATCH -N 1
#SBATCH --mem=5G
#SBATCH --time=1:00:00
#SBATCH --array=2-800
#SBATCH --output=/dev/null

. "$CONDA_PREFIX/etc/profile.d/conda.sh"
conda activate defensefinder

models_dir="$HOME/scratch/dbs/defensefinder"
gendir="../../data/vibrio/annotations/faa"
ind=${SLURM_ARRAY_TASK_ID}
file=$(ls ${gendir}/*.faa | sed -n ${ind}p)
bn1=$(basename $file)
bn=${bn1:0:-4}
outdir="../../output/defensefinder/$bn"

# run
defense-finder run --models-dir $models_dir -w 5 -o $outdir $file
