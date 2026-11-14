"""Build the taxon table from the NCBI data report.

Usage: metadata.py <data_report.jsonl> <taxa_metadata.csv>

One row per genome: accession, organism, ICTV ranks (family/subfamily/genus/species),
host, genome size and the tree label used by every downstream step.
"""
import csv
import sys

import vphylo

FIELDS = ["label", "organism", "accession", "subfamily", "genus", "species", "tax_id",
          "host", "genome_length", "gene_count", "protein_count", "isolate",
          "release_date", "family", "lineage"]


def main(report_path, out_path):
    rows = vphylo.read_report(report_path)
    for r in rows:
        r["label"] = vphylo.label(r["organism"])
    rows.sort(key=lambda r: (r["subfamily"], r["genus"], r["organism"]))

    with open(out_path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)

    by_sub = {}
    for r in rows:
        by_sub[r["subfamily"]] = by_sub.get(r["subfamily"], 0) + 1
    print(f"genomes: {len(rows)}  organisms: {len({r['organism'] for r in rows})}")
    for sub, n in sorted(by_sub.items()):
        print(f"  {sub}: {n}")
    labels = [r["label"] for r in rows]
    dups = sorted({r["label"] for r in rows if labels.count(r["label"]) > 1})
    if dups:
        print(f"viruses with more than one RefSeq genome ({len(dups)}), merged into one taxon "
              f"downstream: {', '.join(dups)}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
