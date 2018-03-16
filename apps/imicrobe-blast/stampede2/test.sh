#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -p normal
#SBATCH -J imcrblst
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user jklynch@email.arizona.edu

OUT_DIR="$SCRATCH/imicrobe-blast/test"

[[ -d "$OUT_DIR" ]] && rm -rf "$OUT_DIR"

DATA_DIR=/work/05066/imicrobe/iplantc.org/data/mmetsp/fasta

# these data files are large
# -rw------- 1 kyclark G-818024   54525952 Oct 10 12:36 1675-MMETSP0011_2-Rhodosorus-marinus-769.1.fa
# -rw------- 1 kyclark G-818024    1638400 Oct  9 19:28 2283-MMETSP1452-Synchroma-pusillum-CCMP3072.1.fa

sh run.sh -q "$DATA_DIR/2283-MMETSP1452-Synchroma-pusillum-CCMP3072.1.fa" -q "$DATA_DIR/1675-MMETSP0011_2-Rhodosorus-marinus-769.1.fa" -o $OUT_DIR
