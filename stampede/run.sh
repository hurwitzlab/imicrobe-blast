#!/bin/bash

set -u

function lc() {
  wc -l "$1" | cut -d ' ' -f 1
}

function HELP() {
  printf "Usage:\n  %s -q QUERY_DIR -o OUT_DIR\n\n" $(basename $0)

  echo "Required arguments:"
  echo " -q QUERY_DIR"
  echo " -o OUT_DIR"
  echo ""
  exit 0
}

if [[ $# -eq 0 ]]; then
  HELP
fi

BIN=$( cd "$( dirname "$0" )" && pwd )
chmod +x *.pl

echo "----------"
echo "BIN \"$BIN\""
echo "Contents of $BIN"
ls -lh $BIN
echo "----------"

QUERY_DIR=""
OUT_DIR=$BIN

while getopts :o:q:h OPT; do
  case $OPT in
    h)
      HELP
      ;;
    o)
      OUT_DIR="$OPTARG"
      ;;
    q)
      QUERY_DIR="$OPTARG"
      ;;
    :)
      echo "Error: Option -$OPTARG requires an argument."
      exit 1
      ;;
    \?)
      echo "Error: Invalid option: -${OPTARG:-""}"
      exit 1
  esac
done

if [[ ! -d $QUERY_DIR ]]; then
  echo QUERY_DIR \"$QUERY_DIR\" does not exist.
  exit 1
fi

if [[ ! -d $OUT_DIR ]]; then
  mkdir -p "$OUT_DIR"
fi

BLAST_OUT_DIR="$OUT_DIR/blast-out"
if [[ ! -d $BLAST_OUT_DIR ]]; then
  mkdir -p "$BLAST_OUT_DIR"
fi

#
# Convert BAM files to FASTA if necessary
#
BAM_FILES=$(mktemp)
find "$QUERY_DIR" -name \*.bam > "$BAM_FILES"
NUM_BAM=$(lc "$BAM_FILES")

if [[ $NUM_BAM -gt 0 ]]; then
  while read BAM_FILE; do
    BASENAME=$(basename $BAM_FILE '.bam')
    FASTA="$QUERY_DIR/${BASENAME}.fa"

    if [[ ! -s $FASTA ]]; then
      echo "Converting BAM_FILE \"$BASENAME\""
      samtools fasta -0 "$FASTA" "$BAM_FILE"
    fi
  done < $BAM_FILES
fi
rm "$BAM_FILES"

INPUT_FILES=$(mktemp)
find "$QUERY_DIR" -type f -not -name \*.bam > "$INPUT_FILES"
NUM_INPUT=$(lc "$INPUT_FILES")

if [[ $INPUT_FILES -lt 1 ]]; then
  echo No files found in QUERY_DIR \"$QUERY_DIR\"
  exit 1
fi

BLAST_DIR="$SCRATCH/imicrobe/blast"

if [[ ! -d $BLAST_DIR ]]; then
  echo BLAST_DIR \"$BLAST_DIR\" does not exist.
  exit 1
fi

BLAST_DBS=$(mktemp)
find "$BLAST_DIR" -type f -name | sed "s/\.fa..*/.fa/" | uniq > "$BLAST_DBS"
NUM_BLAST=$(lc "$BLAST_DBS")

if [[ $NUM_BLAST -lt 1 ]]; then
  echo No files found in BLAST_DIR \"$BLAST_DIR\"
  exit 1
fi

BLAST_ARGS="-perc_identity .98 -outfmt 6"
PARAM="$$.param"

i=0
while read INPUT_FILE; do
  BASENAME=$(basename $INPUT_FILE)

  let i++
  printf "%3d: %s\n" $i $BASENAME

  while read BLAST_DB; do
    echo "blastn $BLAST_ARGS -db $BLAST_DB -query $INPUT_FILE -out $BLAST_OUT_DIR/$BASENAME" >> $PARAM
  done < $BLAST_DBS
done < $INPUT_FILES

echo "Starting launcher $(date)"
export LAUNCHER_DIR="$HOME/src/launcher"
export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=$PARAM
export LAUNCHER_PPN=2

$LAUNCHER_DIR/paramrun
echo "Ended launcher $(date)"

rm $INPUT_FILES
rm $BLAST_DBS
rm $PARAM
