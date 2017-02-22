#!/bin/bash

# Author: Ken Youens-Clark <kyclark@email.arizona.edu>

set -u

BIN=$( cd "$( dirname "$0" )" && pwd )
QUERY=""
PCT_ID=".98"
OUT_DIR="$BIN"
NUM_THREADS=12

module load blast

function lc() {
  wc -l "$1" | cut -d ' ' -f 1
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

# Here is a place for fasplit.py to ensure not too 
# many sequences in each query.

BLAST_DIR="$WORK/ohana/blast"

if [[ ! -d "$BLAST_DIR" ]]; then
  echo "BLAST_DIR \"$BLAST_DIR\" does not exist."
  exit 1
fi

BLAST_DIR="$WORK/ohana/blast"
BLAST_ARGS="-outfmt 6 -num_threads $NUM_THREADS"
BLAST_PARAM="$$.blast.param"

cat /dev/null > $BLAST_PARAM # make sure it's empty

i=0
while read INPUT_FILE; do
  BASENAME=$(basename "$INPUT_FILE")

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
    echo "$BLAST_TO_DNA $BLAST_ARGS -perc_identity $PCT_ID -db $BLAST_DIR/contigs -query $INPUT_FILE -out $BLAST_OUT_DIR/$BASENAME-contigs.tab" >> $BLAST_PARAM
    echo "$BLAST_TO_DNA $BLAST_ARGS -perc_identity $PCT_ID -db $BLAST_DIR/genes -query $INPUT_FILE -out $BLAST_OUT_DIR/$BASENAME-genes.tab" >> $BLAST_PARAM
  fi

  BLAST_TO_PROT=""
  if [[ $TYPE == 'dna' ]]; then 
    BLAST_TO_PROT='blastx'
  elif [[ $TYPE == 'prot' ]]; then
    BLAST_TO_PROT='blastp'
  else
    echo "Cannot BLAST $BASENAME to PROT (not DNA or prot)"
  fi

  if [[ ${#BLAST_TO_PROT} -gt 0 ]]; then
    echo "$BLAST_TO_PROT $BLAST_ARGS -db $BLAST_DIR/proteins -query $INPUT_FILE -out $BLAST_OUT_DIR/$BASENAME-proteins.tab" >> $BLAST_PARAM
  fi
done < "$INPUT_FILES"
rm "$INPUT_FILES"

echo "Starting launcher for BLAST"
export LAUNCHER_DIR="$HOME/src/launcher"
export LAUNCHER_NJOBS=$(lc $BLAST_PARAM)
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
# Now we need to add Eggnog (and eventually Pfam, KEGG, etc.)
# annotations to the "*-genes.tab" files.
# 
ANNOT_PARAM="$$.annot.param"
cat /dev/null > $ANNOT_PARAM

GENE_HITS=$(mktemp)
find $BLAST_OUT_DIR -size +0c -name \*-genes.tab > $GENE_HITS
while read FILE; do
  BASENAME=$(basename $FILE '.tab')
  echo "Annotating $FILE"
  echo "annotate.py -b \"$FILE\" -a \"${WORK}/ohana/sqlite\" -o \"${OUT_DIR}/annotations\"" >> $ANNOT_PARAM
done < $GENE_HITS

# Probably should run the above annotation with launcher, but I was 
# having problems with this.
echo "Starting launcher for annotation"
export LAUNCHER_NHOSTS=1
export LAUNCHER_NJOBS=$(lc $ANNOT_PARAM)
export LAUNCHER_JOB_FILE=$ANNOT_PARAM
$LAUNCHER_DIR/paramrun
echo "Ended launcher for annotation"
rm "$ANNOT_PARAM"
