#!/bin/bash

# Author: Ken Youens-Clark <kyclark@email.arizona.edu>

set -u

BIN=$( cd "$( dirname "$0" )" && pwd )
QUERY=""
PCT_ID=".98"
OUT_DIR="$BIN"
NUM_THREADS=12
MAX_SEQS_PER_FILE=25

module load blast

function lc() {
  wc -l "$1" | awk '{print $1}'
}

function HELP() {
  printf "Usage:\n  %s -q QUERY -o OUT_DIR\n\n" $(basename $0)

  echo "Required arguments:"
  echo " -q QUERY"
  echo
  echo "Options:"
  echo
  echo " -p PCT_ID ($PCT_ID)"
  echo " -o OUT_DIR ($OUT_DIR)"
  echo " -n NUM_THREADS ($NUM_THREADS)"
  echo 
  exit 0
}

if [[ $# -eq 0 ]]; then
  HELP
fi

while getopts :o:n:p:q:h OPT; do
  case $OPT in
    h)
      HELP
      ;;
    n)
      NUM_THREADS="$OPTARG"
      ;;
    o)
      OUT_DIR="$OPTARG"
      ;;
    p)
      PCT_ID="$OPTARG"
      ;;
    q)
      QUERY="$OPTARG"
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

# 
# TACC docs recommend tar'ing a "bin" dir of scripts in order 
# to maintain file permissions such as the executable bit; 
# otherwise, you would need to "chmod +x" the files or execute
# like "python script.py ..."
#
SCRIPTS="scripts.tgz"
if [[ -e $SCRIPTS ]]; then
  echo "Untarring $SCRIPTS to bin"
  if [[ ! -d bin ]]; then
    mkdir bin
  fi
  tar -C bin -xvf $SCRIPTS
fi

if [[ -e "$BIN/bin" ]]; then
  PATH="$BIN/bin:$PATH"
fi

if [[ $NUM_THREADS -lt 1 ]]; then
  echo "NUM_THREADS \"$NUM_THREADS\" cannot be less than 1"
  exit 1
fi

if [[ -d "$OUT_DIR" ]]; then
  mkdir -p "$OUT_DIR"
fi

BLAST_OUT_DIR="$OUT_DIR/blast-out"
if [[ ! -d "$BLAST_OUT_DIR" ]]; then
  mkdir -p "$BLAST_OUT_DIR"
fi

INPUT_FILES=$(mktemp)
if [[ -d $QUERY ]]; then
  find "$QUERY" -type f > "$INPUT_FILES"
else
  echo "$QUERY" > $INPUT_FILES
fi
NUM_INPUT=$(lc "$INPUT_FILES")

if [[ $NUM_INPUT -lt 1 ]]; then
  echo "No input files found"
  exit 1
fi

SPLIT_DIR="$OUT_DIR/split"
if [[ ! -d $SPLIT_DIR ]]; then
  mkdir -p $SPLIT_DIR
fi

while read FILE; do
  fasplit.py -f $FILE -o $SPLIT_DIR/$(basename $FILE) -n $MAX_SEQS_PER_FILE
done < $INPUT_FILES

SPLIT_FILES=$(mktemp)
find $SPLIT_DIR -type f -size +0c > $SPLIT_FILES
NUM_SPLIT=$(lc $SPLIT_FILES)

echo "After splitting, there are NUM_SPLIT \"$NUM_SPLIT\""

BLAST_DIR="$WORK/imicrobe/blast"

if [[ ! -d "$BLAST_DIR" ]]; then
  echo "BLAST_DIR \"$BLAST_DIR\" does not exist."
  exit 1
fi

BLAST_ARGS="-outfmt 6 -num_threads $NUM_THREADS"
BLAST_PARAM="$$.blast.param"
BLAST_DBS=$(mktemp)

find $BLAST_DIR -type f -size +0c | sed "s/\..*//" | sort | uniq > $BLAST_DBS

cat /dev/null > $BLAST_PARAM # make sure it's empty

i=0
while read INPUT_FILE; do
  SAMPLE_NAME=$(basename $(dirname $INPUT_FILE))
  BASENAME=$(basename "$INPUT_FILE")
  SAMPLE_DIR="${BLAST_OUT_DIR}/${SAMPLE_NAME}"

  echo "SAMPLE_DIR \"$SAMPLE_DIR\""

  if [[ ! -e $SAMPLE_DIR ]]; then
    mkdir -p "$SAMPLE_DIR"
  fi

  let i++
  printf "%3d: %s\n" "$i" "$BASENAME"
  EXT="${BASENAME##*.}"
  TYPE="unknown"
  if [[ $EXT == 'fa'    ]] || \
     [[ $EXT == 'fna'   ]] || \
     [[ $EXT == 'fas'   ]] || \
     [[ $EXT == 'fasta' ]] || \
     [[ $EXT == 'ffn'   ]];
  then
    TYPE="dna"
  elif [[ $EXT == 'faa' ]]; then
    TYPE="prot"
  elif [[ $EXT == 'fra' ]]; then
    TYPE="rna"
  fi

  BLAST_TO_DNA=""
  if [[ $TYPE == 'dna' ]]; then 
    BLAST_TO_DNA='blastn'
  elif [[ $TYPE == 'prot' ]]; then
    BLAST_TO_DNA='tblastn'
  else
    echo "Cannot BLAST $BASENAME to DNA (not DNA or prot)"
  fi

  if [[ ${#BLAST_TO_DNA} -gt 0 ]]; then
    while read DB; do
      BASE_DB=$(basename $DB)
      DB_DIR="${SAMPLE_DIR}/${BASE_DB}"

      if [[ ! -e $DB_DIR ]]; then
        mkdir -p "$DB_DIR"
      fi

      echo "$BLAST_TO_DNA $BLAST_ARGS -perc_identity $PCT_ID -db \"$BLAST_DIR/$BASE_DB\" -query \"$INPUT_FILE\" -out \"${DB_DIR}/${BASENAME}\"" >> $BLAST_PARAM
    done < $BLAST_DBS
  fi

#  BLAST_TO_PROT=""
#  if [[ $TYPE == 'dna' ]]; then 
#    BLAST_TO_PROT='blastx'
#  elif [[ $TYPE == 'prot' ]]; then
#    BLAST_TO_PROT='blastp'
#  else
#    echo "Cannot BLAST $BASENAME to PROT (not DNA or prot)"
#  fi
#
#  if [[ ${#BLAST_TO_PROT} -gt 0 ]]; then
#    echo "$BLAST_TO_PROT $BLAST_ARGS -db $BLAST_DIR/proteins -query $INPUT_FILE -out $BLAST_OUT_DIR/$BASENAME-proteins.tab" >> $BLAST_PARAM
#  fi
done < "$SPLIT_FILES"
rm "$SPLIT_FILES"

NUM_JOBS=$(lc $PARAM)

echo "Starting launcher NUM_JOBS \"$NUM_JOBS\" for BLAST"
export LAUNCHER_DIR="$HOME/src/launcher"
export LAUNCHER_NJOBS=$NUM_JOBS
export LAUNCHER_NHOSTS=4
export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_WORKDIR=$BIN
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=$BLAST_PARAM
export LAUNCHER_PPN=4
export LAUNCHER_SCHED=interleaved
$LAUNCHER_DIR/paramrun
echo "Ended launcher for BLAST"
rm $BLAST_PARAM

#
# Remove the empty files, directories
#
find $BLAST_OUT_DIR -type f -size 0 -exec rm {} \;
find $BLAST_OUT_DIR -type d -empty -exec rmdir {} \;

#
# Rollup the 1/2/... files into one
# Add annotations
#
SAMPLE_DIRS=$(mktemp)
find $BLAST_OUT_DIR -maxdepth 1 -mindepth 1 -type d > $SAMPLE_DIRS
ANNOT_PARAM="$$.annot.param"

while read SAMPLE_DIR; do
  echo "SAMPLE_DIR $SAMPLE_DIR"
  DB_DIRS=$(mktemp)
  find $SAMPLE_DIR -maxdepth 1 -mindepth 1 -type d > $DB_DIRS

  while read DB_DIR; do
    echo "DB_DIR $DB_DIR"
    DB_DIR_BASE=$(basename $DB_DIR)
    SUMMARY="${SAMPLE_DIR}/${DB_DIR_BASE}.tab"
    echo "SUMMARY $SUMMARY"
    cat $DB_DIR/* > $SUMMARY
    rm -rf $DB_DIR
  done < $DB_DIRS
  rm $DB_DIRS

  echo "annotate.py -b \"$SAMPLE_DIR\" -a \"${WORK}/imicrobe/annotations/annots.db\" -o \"${OUT_DIR}/annotations/$(basename $SAMPLE_DIR)\"" > $ANNOT_PARAM
done < $SAMPLE_DIRS

rm $SAMPLE_DIRS

NUM_JOBS=$(lc $ANNOT_PARAM)
echo "Starting for annotation"
export LAUNCHER_NHOSTS=1
export LAUNCHER_NJOBS=$(lc $ANNOT_PARAM)
export LAUNCHER_JOB_FILE=$ANNOT_PARAM
export LAUNCHER_PPN=4
$LAUNCHER_DIR/paramrun
echo "Ended launcher for annotation"

rm -rf "$SPLIT_DIR"
rm "$ANNOT_PARAM"
