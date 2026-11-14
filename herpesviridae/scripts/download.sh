#!/usr/bin/env bash
# Download the RefSeq genomes and proteomes of a viral taxon from NCBI Virus.
#
# Usage: download.sh <taxon> <include> <refseq_only:true|false> <out_dir>
# Writes <out_dir>/{genomic.fna,protein.faa,data_report.jsonl}.
set -euo pipefail

taxon=$1
include=$2
refseq_only=$3
out_dir=$4

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

args=(virus genome taxon "$taxon" --include "$include" --filename "$tmp/dataset.zip")
if [ "$refseq_only" = "true" ]; then
    args+=(--refseq)
fi

echo "== datasets download ${args[*]}"
datasets download "${args[@]}"
unzip -oq "$tmp/dataset.zip" -d "$tmp/unpacked"

mkdir -p "$out_dir"
src="$tmp/unpacked/ncbi_dataset/data"
for f in genomic.fna protein.faa data_report.jsonl; do
    if [ -f "$src/$f" ]; then
        cp "$src/$f" "$out_dir/$f"
    else
        echo "ERROR: $f missing from the downloaded package" >&2
        exit 1
    fi
done

echo "== downloaded"
echo "genomes:  $(grep -c '^>' "$out_dir/genomic.fna")"
echo "proteins: $(grep -c '^>' "$out_dir/protein.faa")"
echo "records:  $(wc -l < "$out_dir/data_report.jsonl")"
