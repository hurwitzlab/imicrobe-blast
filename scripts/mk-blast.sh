#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -J mkblast
#SBATCH --mail-user=kyclark@email.arizona.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -u

function lc() {
  wc -l $1 | awk '{print $1}'
}

module load blast

BASE_DIR="$WORK/iplantc.org/data/imicrobe"
IN_DIR="$BASE_DIR/projects"
OUT_DIR="$BASE_DIR/blast2"

[[ ! -d "$OUT_DIR" ]] && mkdir -p "$OUT_DIR"

PROJECT_DIRS=$(mktemp)
find "$IN_DIR" -mindepth 1 -maxdepth 1 -type d > "$PROJECT_DIRS"
NUM_PROJECT_DIRS=$(lc "$PROJECT_DIRS")

if [[ $NUM_PROJECT_DIRS -lt 1 ]]; then
    echo "Found 0 project dirs in IN_DIR \"$IN_DIR\""
    exit 1
fi

PARAM="$$.param"
cat /dev/null > "$PARAM"

while read -r PROJECT_DIR; do
    PROJECT_ID=$(basename "$PROJECT_DIR")
    FILES=$(mktemp)
    find "$PROJECT_DIR" -type f -size +0c \( -name \*.fa -o -name \*.fna -o -name \*.fasta \) > "$FILES"
    NUM_FILES=$(lc "$FILES")

    if [[ $NUM_FILES -lt 1 ]]; then
        echo "Found no files in PROJECT_DIR \"$PROJECT_DIR\""
        continue
    fi

    echo "PROJECT_ID \"$PROJECT_ID\" NUM_FILES \"$NUM_FILES\""

    i=0
    while read -r FILE; do
        BASENAME=$(basename $FILE)
        EXT="${BASENAME##*.}"
        SAMPLE_ID=$(basename $(dirname "$FILE"))
        SAMPLE_DIR="$OUT_DIR/$PROJECT_ID/$SAMPLE_ID"
  
        [[ ! -d "$SAMPLE_DIR" ]] && mkdir -p "$SAMPLE_DIR"
  
        let i++
        printf "%5d: %s => %s\n" $i $BASENAME $SAMPLE_DIR
  
        TITLE="iMicrobe Sample $SAMPLE_ID"
        ARGS="-dbtype nucl -title \"$TITLE\" -out $SAMPLE_DIR/$BASENAME"
  
        if [[ $EXT == 'gz' ]]; then
            echo "gunzip -c $FILE | makeblastdb -in - $ARGS" >> $PARAM
        else
            echo "makeblastdb -in $FILE $ARGS" >> $PARAM
        fi
    done < "$FILES"
    rm "$FILES"
done < "$PROJECT_DIRS"

NJOBS=$(lc "$PARAM")

echo "Starting $NJOBS parallel jobs... $(date)"
export LAUNCHER_DIR="$HOME/src/launcher"
export LAUNCHER_WORKDIR="$PWD"
export LAUNCHER_PPN=8
export LAUNCHER_JOB_FILE="$PARAM"
export LAUNCHER_RMI=SLURM
export LAUNCHER_NJOBS=$NJOBS
time $LAUNCHER_DIR/paramrun
echo "Finished parallel $(date)"

echo "Done."
