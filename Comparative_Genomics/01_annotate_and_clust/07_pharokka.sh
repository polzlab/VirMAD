#!/usr/bin/env bash
#SBATCH --job-name=pharokka
#SBATCH --cpus-per-task=20
#SBATCH -N 1
#SBATCH --mem=5G
#SBATCH --time=10:00:00

. "$CONDA_PREFIX/etc/profile.d/conda.sh"
conda activate pharokka

dbdir="$HOME/scratch/dbs/pharokka"
protpath="../../output/protclu/phages/prots.faa"
outdir="../../output/pharokka"

pharokka_proteins.py -i $protpath -o $outdir -d $dbdir -t 20
