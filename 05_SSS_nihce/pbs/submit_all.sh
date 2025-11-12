#!/bin/bash

GENE_LIST_DIR="/dellstorage09/quj_lab/yanghang/spatial/ref"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

N_GENE_LISTS=$(ls "$GENE_LIST_DIR"/*.txt 2>/dev/null | wc -l)

if [[ $N_GENE_LISTS -eq 0 ]]; then
    echo "Error: No .txt files found in $GENE_LIST_DIR"
    exit 1
fi

echo "Found $N_GENE_LISTS gene list files"
echo "Submitting array job: 1-${N_GENE_LISTS}"

qsub -t 1-${N_GENE_LISTS} "$SCRIPT_DIR/submit_niche_analysis.pbs"

echo "Jobs submitted successfully"