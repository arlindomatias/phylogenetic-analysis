#!/usr/bin/env bash
# Align one marker gene and trim poorly aligned columns.
#
# Usage: align.sh <marker.faa> <marker.aln.faa> <marker.trim.faa> <gap_threshold> [threads]
set -euo pipefail

in=$1
aln=$2
trimmed=$3
gap_threshold=$4
threads=${5:-1}

mafft --auto --anysymbol --thread "$threads" "$in" > "$aln"

# -gt keeps columns present in at least this fraction of taxa; -cons guarantees that
# trimAl never returns fewer than 20% of the columns on a divergent alignment.
trimal -in "$aln" -out "$trimmed" -gt "$gap_threshold" -cons 20

cols_in=$(awk '/^>/{if (n++) exit; next} {l += length($0)} END {print l}' "$aln")
cols_out=$(awk '/^>/{if (n++) exit; next} {l += length($0)} END {print l}' "$trimmed")
echo "$(basename "$in"): $(grep -c '^>' "$in") taxa, ${cols_in} -> ${cols_out} columns"
