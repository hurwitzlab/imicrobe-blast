#!/bin/bash

# Author: Ken Youens-Clark <kyclark@email.arizona.edu>
# Author: Joshua Lynch <jklynch@email.arizona.edu>

module load tacc-singularity

# the next two lines load 'testing' version
# of the TACC launcher
module use /scratch/01255/siliu/modulefiles
module load launcher/3.2

module load blast

module list

set -u
set -e

QUERY=""
PCT_ID=".98"
OUT_DIR="$PWD/ohana-blast-out"
IMG="imicrobe-blast-0.0.5.img"
BLAST_DB_DIR="/work/05066/imicrobe/iplantc.org/data/ohana/blast/db24/genes"
# this is a BLAST environment variable
export BLASTDB=$BLAST_DB_DIR

PARAMRUN="$TACC_LAUNCHER_DIR/paramrun"
echo "PARAMRUN: ${PARAMRUN}"

export LAUNCHER_WORKDIR="$PWD"
export LAUNCHER_SCHED=dynamic
export LAUNCHER_PPN=24

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
    #echo " -n NUM_THREADS ($NUM_THREADS)"
    echo " -b BLAST_DB_LIST"
    echo
    exit ${1:-0}
}

[[ $# -eq 0 ]] && HELP 1

while getopts :b:n:o:p:q:h OPT; do
    case $OPT in
        h)
            HELP
            ;;
        b)
            BLAST_DB_LIST="$OPTARG"
            ;;
        #n)
        #    NUM_THREADS="$OPTARG"
        #    ;;
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

#if [[ ! -f "$BLAST_DB_DIR/$BLAST_DB.nal" ]]; then
#    echo "BLAST_DB \"$BLAST_DB\" does not exist."
#    exit 1
#fi

#if [[ $NUM_THREADS -lt 1 ]]; then
#    echo "NUM_THREADS \"$NUM_THREADS\" cannot be less than 1"
#    exit 1
#fi

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
# Split the input files
#
#SPLIT_DIR="$OUT_DIR/split"
#echo "SPLIT_DIR \"$SPLIT_DIR\""
#[[ ! -d "$SPLIT_DIR" ]] && mkdir -p "$SPLIT_DIR"
#
#SPLIT_PARAM="$$.split.param"
#while read -r FILE; do
#    BASENAME=$(basename "$FILE")
#    echo "singularity exec $IMG fasplit -f $FILE -o $SPLIT_DIR -n 24" >> "$SPLIT_PARAM"
#done < "$INPUT_FILES"
#
#NUM_SPLIT=$(lc "$SPLIT_PARAM")
#echo "Starting launcher NUM_SPLIT \"$NUM_SPLIT\" for split"
#export LAUNCHER_JOB_FILE="$SPLIT_PARAM"
#$PARAMRUN
#echo "Ended launcher for split"
#
#SPLIT_FILES=$(mktemp)
#find "$SPLIT_DIR" -type f -size +0c > "$SPLIT_FILES"
#NUM_SPLIT=$(lc "$SPLIT_FILES")
#
#echo "After splitting, there are NUM_SPLIT \"$NUM_SPLIT\""
#
#if [[ "$NUM_SPLIT" -lt 1 ]]; then
#    echo "Something went wrong with splitting."
#    exit 1
#fi

#
# Run BLAST
#
BLAST_PARAM="$$.blast.param"
BLAST_ARGS="-outfmt 6 -num_threads 2"

FILE_NUM=0
while read -r INPUT_FILE; do
    INPUT_NAME=$(basename "$INPUT_FILE")
    QUERY_OUT_DIR="$BLAST_OUT_DIR"

    [[ ! -d "$QUERY_OUT_DIR" ]] && mkdir -p "$QUERY_OUT_DIR"

    FILE_NUM=$((FILE_NUM + 1))
    printf "%5d: QUERY %s\n" "$FILE_NUM" "$INPUT_NAME"
    ##EXT="${QUERY_NAME##*.}"
    EXT="${INPUT_NAME##*.}"
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
        echo "Cannot BLAST \"$INPUT_FILE\" to DNA (not DNA or prot)"
    fi

    if [[ ${#BLAST_TO_DNA} -gt 0 ]]; then
        HITS_DIR="$QUERY_OUT_DIR"
        [[ ! -d "$HITS_DIR" ]] && mkdir -p "$HITS_DIR"

        while read -r BLAST_DB; do
            echo "$BLAST_TO_DNA $BLAST_ARGS -perc_identity $PCT_ID -db \"$BLAST_DB\" -query $INPUT_FILE -out $HITS_DIR/$INPUT_NAME-$BLAST_DB" >> "$BLAST_PARAM"
        done < blast_dbs_24.txt

    fi
done < "$INPUT_FILES"

NUM_SPLIT=$(lc "$BLAST_PARAM")
echo "Starting $LAUNCHER_PPN launcher tasks"
echo "Starting $NUM_SPLIT launcher jobs"
export LAUNCHER_JOB_FILE="$BLAST_PARAM"
$PARAMRUN
echo "Ended launcher for BLAST"


echo "Done."
echo "Comments to Ken Youens-Clark kyclark@email.arizona.edu"
