"""Measure how far the RefSeq annotation text alone would have got us.

Usage: annotation_check.py <config.yaml> <protein.faa> <annotation_check.csv>

For every marker, counts the viruses that carry a protein whose product name is literally the
name of the HSV-1 query protein. Comparing that count with the homology-based recovery in
`marker_coverage.csv` is the evidence for building the orthologue sets by homology instead of
by annotation text.
"""
import csv
import sys

import yaml

import vphylo


def main(config_path, faa_path, out_path):
    config = yaml.safe_load(open(config_path))
    by_accession = {acc: marker for marker, acc in config["markers"].items()}

    query_name = {}
    taxa = set()
    names_per_taxon = {}
    for header, _ in vphylo.read_fasta(faa_path):
        info = vphylo.parse_header(header)
        if not info["organism"]:
            continue
        taxon = vphylo.label(info["organism"])
        taxa.add(taxon)
        names_per_taxon.setdefault(taxon, set()).add(info["name"].strip().lower())
        marker = by_accession.get(vphylo.protein_accession(info["id"]))
        if marker:
            query_name[marker] = info["name"].strip()

    rows = []
    for marker in sorted(config["markers"]):
        name = query_name.get(marker, "")
        n = sum(1 for t in taxa if name.lower() in names_per_taxon[t]) if name else 0
        rows.append({"marker": marker, "query_product_name": name,
                     "n_taxa_with_that_exact_name": n,
                     "pct_taxa": round(100.0 * n / len(taxa), 1) if taxa else 0.0})

    with open(out_path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print(f"viruses in the proteome file: {len(taxa)}")
    for r in rows:
        print(f"  {r['marker']:22s} {r['n_taxa_with_that_exact_name']:4d} viruses annotate a protein "
              f"named '{r['query_product_name']}' ({r['pct_taxa']}%)")


if __name__ == "__main__":
    main(*sys.argv[1:4])
