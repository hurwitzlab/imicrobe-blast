#!/bin/bash

# Author: Ken Youens-Clark <kyclark@email.arizona.edu>

set -u

QUERY=""
PCT_ID=".98"
OUT_DIR="$PWD/blast-out"
NUM_THREADS=12
MAX_SEQS_PER_FILE=1000
IMG="imicrobe-blast-0.0.3.img"
BLAST_DB_DIR="/work/05066/imicrobe/iplantc.org/data/blast/dbs"
ANNOT_DB="/work/05066/imicrobe/iplantc.org/data/imicrobe-annotdb/annots.db"
BLAST_DB_LIST=""

export LAUNCHER_DIR="$HOME/src/launcher"
export LAUNCHER_PLUGIN_DIR="$LAUNCHER_DIR/plugins"
export LAUNCHER_WORKDIR="$PWD"
export LAUNCHER_RMI=SLURM
export LAUNCHER_SCHED=interleaved

PATH="$LAUNCHER_DIR:$PATH"

module load blast

function lc() {
  wc -l "$1" | awk '{print $1}'
}

function HELP() {
  printf "Usage:\n  %s -q QUERY -o OUT_DIR\n\n" "$(basename "$0")"

  echo "Required arguments:"
  echo " -q QUERY"
  echo
  echo "Options:"
  echo
  echo " -p PCT_ID ($PCT_ID)"
  echo " -o OUT_DIR ($OUT_DIR)"
  echo " -n NUM_THREADS ($NUM_THREADS)"
  echo " -b BLAST_DB_LIST"
  echo 
  exit 0
}

if [[ $# -eq 0 ]]; then
  HELP
fi

while getopts :b:o:n:p:q:h OPT; do
  case $OPT in
    h)
      HELP
      ;;
    b)
      BLAST_DB_LIST="$OPTARG"
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
      QUERY="$QUERY $OPTARG"
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

if [[ ! -d "$BLAST_DB_DIR" ]]; then
    echo "BLAST_DB_DIR \"$BLAST_DB_DIR\" does not exist."
    exit 1
fi

if [[ $NUM_THREADS -lt 1 ]]; then
    echo "NUM_THREADS \"$NUM_THREADS\" cannot be less than 1"
    exit 1
fi

[[ -d "$OUT_DIR" ]] && mkdir -p "$OUT_DIR"

BLAST_OUT_DIR="$OUT_DIR/blast-out"
[[ ! -d "$BLAST_OUT_DIR" ]] && mkdir -p "$BLAST_OUT_DIR"

INPUT_FILES=$(mktemp)
for QRY in $QUERY; do
    if [[ -d "$QRY" ]]; then
        find "$QRY" -type f >> "$INPUT_FILES"
    elif [[ -f "$QRY" ]]; then
        echo "$QRY" >> "$INPUT_FILES"
    else
        echo "$QRY neither file nor directory"
    fi
done

NUM_INPUT=$(lc "$INPUT_FILES")

if [[ $NUM_INPUT -lt 1 ]]; then
    echo "No input files found"
    exit 1
else
    echo "Found $NUM_INPUT input files"
fi

#
# Split the input files
#
SPLIT_DIR="$OUT_DIR/split"
echo "SPLIT_DIR \"$SPLIT_DIR\""
[[ ! -d "$SPLIT_DIR" ]] && mkdir -p "$SPLIT_DIR"

SPLIT_PARAM="$$.split.param"
while read -r FILE; do
    BASENAME=$(basename "$FILE")
    echo "singularity exec $IMG fasplit.py -f $FILE -o $SPLIT_DIR/$BASENAME -n $MAX_SEQS_PER_FILE" >> "$SPLIT_PARAM"
done < "$INPUT_FILES"

NUM_SPLIT=$(lc "$SPLIT_PARAM")
echo "Starting launcher NUM_SPLIT \"$NUM_SPLIT\" for split"
LAUNCHER_JOB_FILE="$SPLIT_PARAM"
export LAUNCHER_JOB_FILE
paramrun
echo "Ended launcher for split"

SPLIT_FILES=$(mktemp)
find "$SPLIT_DIR" -type f -size +0c > "$SPLIT_FILES"
NUM_SPLIT=$(lc "$SPLIT_FILES")

echo "After splitting, there are NUM_SPLIT \"$NUM_SPLIT\""

if [[ "$NUM_SPLIT" -lt 1 ]]; then
    echo "Something went wrong with splitting."
    exit 1
fi

# 
# Run BLAST
#
BLAST_ARGS="-outfmt 6 -num_threads $NUM_THREADS"
BLAST_PARAM="$$.blast.param"
cat /dev/null > "$BLAST_PARAM" 
#BLAST_DBS=$(mktemp)
BLAST_DBS='blast-dbs'

if [[ -z "$BLAST_DB_LIST" ]] && [[ -f "$BLAST_DB_LIST" ]]; then
    TMP=$(mktemp)
    while read -r DB_NAME; do
        echo "Looking for DB_NAME \"$DB_NAME\""
        find "$BLAST_DB_DIR" -name $DB_NAME\\* -type f -size +0c >> "$TMP"
    done < "$BLAST_DB_LIST"
    sed "s/\..*//" "$TMP" | sort | uniq > "$BLAST_DBS"
    rm "$TMP"
else
    find "$BLAST_DB_DIR" -type f -size +0c | perl -pe "s/\.\w+$//d" | sort | uniq > "$BLAST_DBS"
fi

NUM_BLAST=$(lc "$BLAST_DBS")
echo "Found NUM_BLAST \"$NUM_BLAST\" dbs"

if [[ $NUM_BLAST -lt 1 ]]; then
    echo "Cannot continue without BLAST dbs"
    exit 1
fi

i=0
while read -r SPLIT_FILE; do
    SPLIT_NAME=$(basename "$SPLIT_FILE")
    QUERY_NAME=$(basename $(dirname "$SPLIT_FILE"))
    QUERY_OUT_DIR="$BLAST_OUT_DIR/$QUERY_NAME"
  
    [[ ! -e "$QUERY_OUT_DIR" ]] && mkdir -p "$QUERY_OUT_DIR"
    echo "QUERY_OUT_DIR \"$QUERY_OUT_DIR\""
  
    let i++
    printf "%5d: QUERY %s\n" "$i" "$SPLIT_NAME"
    EXT="${QUERY_NAME##*.}"
    TYPE="unknown"
    if [[ $EXT == "fa"    ]] || \
       [[ $EXT == "fna"   ]] || \
       [[ $EXT == "fas"   ]] || \
       [[ $EXT == "fasta" ]] || \
       [[ $EXT == "ffn"   ]];
    then
        TYPE="dna"
    elif [[ $EXT == "faa" ]]; then
        TYPE="prot"
    elif [[ $EXT == "fra" ]]; then
        TYPE="rna"
    fi
  
    BLAST_TO_DNA=""
    if [[ $TYPE == "dna" ]]; then 
        BLAST_TO_DNA="blastn"
    elif [[ $TYPE == "prot" ]]; then
        BLAST_TO_DNA="tblastn"
    else
        echo "Cannot BLAST \"$QUERY_NAME\" to DNA (not DNA or prot)"
    fi
  
    if [[ ${#BLAST_TO_DNA} -gt 0 ]]; then
        i=0
        while read -r BLAST_DB; do
            SAMPLE_ID=$(basename $(dirname "$BLAST_DB")) 
            SAMPLE_NAME=$(basename "$BLAST_DB") 
            HITS_DIR="$QUERY_OUT_DIR/$SAMPLE_ID"
      
            [[ ! -e "$HITS_DIR" ]] && mkdir -p "$HITS_DIR"
      
            echo "singularity exec $IMG $BLAST_TO_DNA $BLAST_ARGS -perc_identity $PCT_ID -db \"$BLAST_DB\" -query \"$SPLIT_FILE\" -out \"$HITS_DIR/$SAMPLE_NAME-$SPLIT_NAME\"" >> "$BLAST_PARAM"
        done < "$BLAST_DBS"
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
#    echo "$BLAST_TO_PROT $BLAST_ARGS -db $BLAST_DIR/proteins -query $QUERY_FILE -out $BLAST_OUT_DIR/$BASENAME-proteins.tab" >> $BLAST_PARAM
#  fi
done < "$SPLIT_FILES"
rm "$SPLIT_FILES"

NUM_JOBS=$(lc "$BLAST_PARAM")

echo "Starting launcher NUM_JOBS \"$NUM_JOBS\" for BLAST"
LAUNCHER_JOB_FILE="$BLAST_PARAM"
export LAUNCHER_JOB_FILE
paramrun
echo "Ended launcher for BLAST"
rm "$BLAST_PARAM"

#
# Remove the empty files, directories
#
echo "Removing empty files/dirs from BLAST_OUT_DIR \"$BLAST_OUT_DIR\""
find "$BLAST_OUT_DIR" -type f -size 0 -exec rm {} \; 2>/dev/null
find "$BLAST_OUT_DIR" -type d -empty -exec rmdir {} \; 2>/dev/null

#
# Rollup the 1/2/... files into one
# Add annotations
#
QUERY_DIRS=$(mktemp)
find "$BLAST_OUT_DIR" -maxdepth 1 -mindepth 1 -type d > "$QUERY_DIRS"

ANNOT_PARAM="$$.annot.param"
cat /dev/null > "$ANNOT_PARAM"

i=0
while read -r QUERY_DIR; do
    let i++
    QUERY_FILE=$(basename "$QUERY_DIR")
    printf "%5d: QUERY %s\n" $i "$QUERY_FILE"
    HITS=$(mktemp)
    find "$QUERY_DIR" -maxdepth 1 -mindepth 1 -type d > "$HITS"
    TOTAL_HITS=0
  
    while read -r HIT_DIR; do
        HIT_ID=$(basename "$HIT_DIR")
        ALL_HITS="$HIT_DIR.tab"
        cat $HIT_DIR/* > "$ALL_HITS"
        rm -rf "$HIT_DIR"
        NUM_HITS=$(lc "$ALL_HITS")
        echo "        => NUM_HITS \"$NUM_HITS\" to SAMPLE_ID \"$HIT_ID\""
        TOTAL_HITS=$((TOTAL_HITS + NUM_HITS))
    done < "$HITS"
    rm "$HITS"
    echo ">>>>> TOTAL HITS \"$TOTAL_HITS\" <<<<<"

    echo "singularity exec $IMG annotate.py -b $QUERY_DIR -a $ANNOT_DB -o $QUERY_DIR/annots.tab" >> "$ANNOT_PARAM"
done < "$QUERY_DIRS"
rm "$QUERY_DIRS"

#
# Annotate the output
#
NUM_JOBS=$(lc "$ANNOT_PARAM")
echo "Starting launcher NUM_JOBS \"$NUM_JOBS\" for annotation"
LAUNCHER_JOB_FILE="$ANNOT_PARAM"
export LAUNCHER_JOB_FILE
paramrun
echo "Ended launcher for annotation"

# 
# Clean up on aisle 3
# 
rm "$ANNOT_PARAM"
rm -rf "$SPLIT_DIR"

echo "Done."
echo "Comments to Ken Youens-Clark kyclark@email.arizona.edu"
