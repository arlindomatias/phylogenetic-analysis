#!/usr/bin/env bash
# Search the marker queries against every protein of the family with DIAMOND blastp.
#
# Usage: homology.sh <queries.faa> <protein.faa> <hits.tsv> <evalue> <sensitivity>
# The subject database is the whole family proteome; hits are filtered later, per
# marker and per genome, by scripts/select_orthologs.py.
set -euo pipefail

queries=$1
subjects=$2
out=$3
evalue=$4
sensitivity=$5

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

diamond makedb --in "$subjects" --db "$tmp/db" --quiet

diamond blastp \
    --query "$queries" \
    --db "$tmp/db" \
    --out "$out" \
    --outfmt 6 qseqid sseqid pident length qlen slen qstart qend evalue bitscore \
    --evalue "$evalue" \
    --"$sensitivity" \
    --max-target-seqs 5000 \
    --quiet

echo "hits: $(wc -l < "$out")"
