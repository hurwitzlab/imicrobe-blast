#!/bin/bash

# Author: Ken Youens-Clark <kyclark@email.arizona.edu>
# Author: Joshua Lynch <jklynch@email.arizona.edu>

module load tacc-singularity
module load launcher
module load blast

set -u
set -e

QUERY=""
PCT_ID=".98"
OUT_DIR="$PWD/blast-out"
NUM_THREADS=48
IMG="imicrobe-blast-0.0.5.img"
##BLAST_DB_DIR="/work/05066/imicrobe/iplantc.org/data/blast/one-db"
BLAST_DB="/work/05066/imicrobe/iplantc.org/data/blast/imicrobe"
ANNOT_DB="/work/05066/imicrobe/iplantc.org/data/imicrobe-annotdb/annots.db"
##BLAST_DB_LIST=""

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
    exit ${1:-0}
}

[[ $# -eq 0 ]] && HELP 1

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

#if [[ ! -d "$BLAST_DB_DIR" ]]; then
#    echo "BLAST_DB_DIR \"$BLAST_DB_DIR\" does not exist."
#    exit 1
#fi

if [[ ! -f "$BLAST_DB.nal" ]]; then
    echo "BLAST_DB \"$BLAST_DB\" does not exist."
    exit 1
fi

if [[ $NUM_THREADS -lt 1 ]]; then
    echo "NUM_THREADS \"$NUM_THREADS\" cannot be less than 1"
    exit 1
fi

[[ -d "$OUT_DIR" ]] && mkdir -p "$OUT_DIR"

BLAST_OUT_DIR="$OUT_DIR/blast"
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
    echo "Found $NUM_INPUT input file(s)"
fi

#
# Run BLAST
#
BLAST_ARGS="-outfmt 6 -num_threads $NUM_THREADS"

FILE_NUM=0
##while read -r SPLIT_FILE; do
while read -r SPLIT_FILE; do
    SPLIT_NAME=$(basename "$SPLIT_FILE")
    #QUERY_NAME=$(basename $(dirname "$SPLIT_FILE"))
    QUERY_OUT_DIR="$BLAST_OUT_DIR"
  
    [[ ! -d "$QUERY_OUT_DIR" ]] && mkdir -p "$QUERY_OUT_DIR"
  
    FILE_NUM=$((FILE_NUM + 1))
    printf "%5d: QUERY %s\n" "$FILE_NUM" "$SPLIT_NAME"
    ##EXT="${QUERY_NAME##*.}"
    EXT="${SPLIT_NAME##*.}"
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
        echo "Cannot BLAST \"$SPLIT_FILE\" to DNA (not DNA or prot)"
    fi
  
    if [[ ${#BLAST_TO_DNA} -gt 0 ]]; then
        HITS_DIR="$QUERY_OUT_DIR"
        [[ ! -d "$HITS_DIR" ]] && mkdir -p "$HITS_DIR"
      
        ##echo "singularity exec $IMG $BLAST_TO_DNA $BLAST_ARGS -perc_identity $PCT_ID -db \"$BLAST_DB\" -query \"$SPLIT_FILE\" -out \"$HITS_DIR/$SAMPLE_NAME-$SPLIT_NAME\"" >> "$BLAST_PARAM"
        echo "$BLAST_TO_DNA $BLAST_ARGS -perc_identity $PCT_ID -db $BLAST_DB -query $SPLIT_FILE -out $HITS_DIR/$SPLIT_NAME"
        $BLAST_TO_DNA $BLAST_ARGS -perc_identity $PCT_ID -db $BLAST_DB -query $SPLIT_FILE -out $HITS_DIR/$SPLIT_NAME
    fi
done < "$INPUT_FILES"
##done < "$SPLIT_FILES"

echo "Done."
echo "Comments to Ken Youens-Clark kyclark@email.arizona.edu"
