#!/usr/bin/env bash
#SBATCH --job-name=bakta
#SBATCH -c 5
#SBATCH -N 1
#SBATCH --mem=30G
#SBATCH --time=5:00:00
#SBATCH --array=1-900

# get slurm array ind
ind=$SLURM_ARRAY_TASK_ID

# define input output variables
indir="/lisc/scratch/dome/pollak/vibrio_matrix/data/vibrio/contigs_fixed_hdrs"
nfiles=$(ls $indir/*.fa | wc -l)
outdir="/lisc/scratch/dome/pollak/vibrio_matrix/data/vibrio/annotations"
if [ "$ind" -le "$nfiles" ]; then

  # check that output file does not exist
  file=$(ls $indir/*.fa | sed -n ${ind}p)
  bn1=$(basename $file)
  bn=${bn1:0:-3}

  if [ ! -f $outdir/faa/$bn.faa ]; then

    # make tmp dir
    otd="$TMPDIR/$ind"
    [[ ! -d $otd ]] && mkdir $otd

    # activate conda environment
    . "$CONDA_PREFIX/etc/profile.d/conda.sh"
    conda activate bakta


    # copy database to $TMPDIR
    db1="/lisc/user/pollak/scratch/dbs/bakta/bakta.db.zst"
    dbdir="$TMPDIR/db"
    [[ ! -d $dbdir ]] && mkdir $dbdir
    echo "Copying database to $TMPDIR"
    cp $db1 $TMPDIR
    echo "extracting archive"
    cd $TMPDIR
    tar -I"unzstd" -xf bakta.db.zst
    rm bakta.db.zst
    cd -


    # run bakta
    bakta --force --skip-crispr --db $dbdir --output $otd --threads 5 $file
    mv $otd/*hypothetical* $outdir/hypotheticals
    mv $otd/*.embl $outdir/embl
    mv $otd/*.faa $outdir/faa
    mv $otd/*.ffn $outdir/ffn
    mv $otd/*.gbff $outdir/gbff
    mv $otd/*.gff3 $outdir/gff
    mv $otd/*.json $outdir/json
    mv $otd/*.log $outdir/logs
    mv $otd/*.tsv $outdir/tsv
    rm -rf $otd
    rm -rf $TMPDIR

  fi
fi
