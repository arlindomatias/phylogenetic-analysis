"""Concatenate the trimmed marker alignments into a partitioned supermatrix.

Usage: concat.py <config.yaml> <supermatrix.fasta> <partitions.nex> <stats.csv> <trimmed.faa>...

A virus enters the supermatrix when it has at least `supermatrix.min_markers` markers;
markers it lacks are filled with gaps, which IQ-TREE treats as missing data.
"""
import csv
import os
import sys

import yaml

import vphylo


def main(config_path, matrix_path, partition_path, stats_path, alignment_paths):
    config = yaml.safe_load(open(config_path))
    min_markers = int(config["supermatrix"]["min_markers"])

    blocks = {}   # marker -> {label: aligned sequence}
    widths = {}
    for path in alignment_paths:
        marker = os.path.basename(path).split(".")[0]
        block = {name.split()[0]: seq for name, seq in vphylo.read_fasta(path)}
        if not block:
            raise SystemExit(f"ERROR: {path} is empty")
        lengths = {len(s) for s in block.values()}
        if len(lengths) != 1:
            raise SystemExit(f"ERROR: {path} is not an alignment (column counts {sorted(lengths)})")
        blocks[marker] = block
        widths[marker] = lengths.pop()

    markers = sorted(blocks)
    all_labels = sorted({lab for b in blocks.values() for lab in b})
    present = {lab: sum(1 for m in markers if lab in blocks[m]) for lab in all_labels}
    labels = [lab for lab in all_labels if present[lab] >= min_markers]
    excluded = [lab for lab in all_labels if present[lab] < min_markers]

    records, rows = [], []
    for lab in labels:
        seq = "".join(blocks[m].get(lab, "-" * widths[m]) for m in markers)
        records.append((lab, seq))
        gaps = seq.count("-") + seq.count("X")
        rows.append({"label": lab, "n_markers": present[lab],
                     "missing_markers": ";".join(m for m in markers if lab not in blocks[m]),
                     "columns": len(seq), "pct_missing": round(100.0 * gaps / len(seq), 1)})
    vphylo.write_fasta(matrix_path, records)

    with open(partition_path, "w") as out:
        out.write("#nexus\nbegin sets;\n")
        start = 1
        for m in markers:
            end = start + widths[m] - 1
            out.write(f"    charset {m} = {start}-{end};\n")
            start = end + 1
        out.write("    charpartition markers = " +
                  ", ".join(f"{config['tree']['model']}:{m}" for m in markers) + ";\nend;\n")

    with open(stats_path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (-r["n_markers"], r["label"])))

    total = sum(widths[m] for m in markers)
    print(f"supermatrix: {len(records)} taxa x {total} columns, {len(markers)} partitions")
    for m in markers:
        print(f"  {m:22s} {widths[m]:5d} columns, {len(blocks[m]):4d} taxa")
    if excluded:
        print(f"excluded (< {min_markers} markers): {len(excluded)}")
        for lab in excluded:
            print(f"  {lab} ({present[lab]} markers)")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5:])
