#!/bin/bash

while read -r BLAST_DB; do
    echo "BLAST DB: $BLAST_DB"
done < blast_dbs_24.txt