"""Pull the HSV-1 marker proteins out of the downloaded proteome, to use as BLAST queries.

Usage: extract_queries.py <config.yaml> <protein.faa> <queries.faa>

Each query is renamed to the marker name from the config, so DIAMOND reports the marker
directly in the `qseqid` column.
"""
import sys

import yaml

import vphylo


def main(config_path, faa_path, out_path):
    config = yaml.safe_load(open(config_path))
    wanted = {acc: name for name, acc in config["markers"].items()}

    found = {}
    for header, seq in vphylo.read_fasta(faa_path):
        info = vphylo.parse_header(header)
        acc = vphylo.protein_accession(info["id"])
        if acc in wanted and acc not in found:
            found[acc] = (wanted[acc], seq, info)

    missing = sorted(set(wanted) - set(found))
    if missing:
        raise SystemExit(f"ERROR: query proteins absent from {faa_path}: {', '.join(missing)}")

    records = [(name, seq) for name, seq, _ in
               sorted(found.values(), key=lambda t: t[0])]
    vphylo.write_fasta(out_path, records)

    print(f"{len(records)} marker queries from {faa_path}")
    for acc, (name, seq, info) in sorted(found.items(), key=lambda t: t[1][0]):
        print(f"  {name:22s} {acc:16s} {len(seq):5d} aa  "
              f"{info['name'] or '(unnamed CDS)'} [{info['organism']}]")


if __name__ == "__main__":
    main(*sys.argv[1:4])
