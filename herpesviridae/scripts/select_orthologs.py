"""Turn DIAMOND hits into one orthologue set per marker gene.

Usage: select_orthologs.py <config.yaml> <protein.faa> <hits.tsv> <out_dir> <coverage.csv> <presence.csv>

For every (marker, virus) pair the best-scoring hit that passes the identity and query
coverage floors is taken as the orthologue. If one protein is the best hit of two markers
in the same virus (paralogy or a fused annotation), it is kept only for the marker with the
higher bit score, so no sequence is used twice in the supermatrix.
"""
import csv
import os
import sys

import yaml

import vphylo

COLS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen",
        "qstart", "qend", "evalue", "bitscore"]


def main(config_path, faa_path, hits_path, out_dir, coverage_path, presence_path):
    config = yaml.safe_load(open(config_path))
    markers = sorted(config["markers"])
    hp = config["homology"]
    min_id, min_qcov = float(hp["min_identity"]), float(hp["min_query_coverage"])

    # viruses are keyed by tree label, which merges NCBI organism names that differ only in
    # capitalisation (two RefSeq genomes of the same virus)
    seqs, taxon_of, organism_of = {}, {}, {}
    for header, seq in vphylo.read_fasta(faa_path):
        info = vphylo.parse_header(header)
        if not info["organism"]:
            continue
        taxon = vphylo.label(info["organism"])
        seqs[info["id"]] = seq
        taxon_of[info["id"]] = taxon
        organism_of.setdefault(taxon, info["organism"])

    # best hit per (marker, organism)
    best = {}   # (marker, taxon label) -> (bitscore, subject id, identity, query coverage)
    kept = dropped_id = dropped_cov = 0
    with open(hits_path) as fh:
        for row in csv.DictReader(fh, fieldnames=COLS, delimiter="\t"):
            marker, sid = row["qseqid"], row["sseqid"]
            taxon = taxon_of.get(sid)
            if marker not in config["markers"] or not taxon:
                continue
            if float(row["pident"]) < min_id:
                dropped_id += 1
                continue
            qcov = 100.0 * (int(row["qend"]) - int(row["qstart"]) + 1) / int(row["qlen"])
            if qcov < min_qcov:
                dropped_cov += 1
                continue
            kept += 1
            score = float(row["bitscore"])
            key = (marker, taxon)
            if key not in best or score > best[key][0]:
                best[key] = (score, sid, float(row["pident"]), qcov)

    # one protein cannot serve two markers in the same virus
    by_subject = {}
    for (marker, taxon), (score, sid, pid, qcov) in best.items():
        by_subject.setdefault((taxon, sid), []).append((score, marker))
    conflicts = 0
    for (taxon, sid), entries in by_subject.items():
        if len(entries) > 1:
            entries.sort(reverse=True)
            for _, marker in entries[1:]:
                del best[(marker, taxon)]
                conflicts += 1

    taxa = sorted({taxon for _, taxon in best})
    os.makedirs(out_dir, exist_ok=True)
    stats = []
    for marker in markers:
        records = []
        for taxon in taxa:
            hit = best.get((marker, taxon))
            if hit:
                records.append((taxon, seqs[hit[1]]))
        vphylo.write_fasta(os.path.join(out_dir, f"{marker}.faa"), records)
        lengths = [len(s) for _, s in records]
        ids = [best[(marker, t)][2] for t in taxa if (marker, t) in best]
        stats.append({
            "marker": marker,
            "query": config["markers"][marker],
            "n_taxa": len(records),
            "pct_taxa": round(100.0 * len(records) / len(taxa), 1) if taxa else 0.0,
            "min_length": min(lengths) if lengths else 0,
            "median_length": sorted(lengths)[len(lengths) // 2] if lengths else 0,
            "max_length": max(lengths) if lengths else 0,
            "median_identity_to_query": round(sorted(ids)[len(ids) // 2], 1) if ids else 0.0,
            "min_identity_to_query": round(min(ids), 1) if ids else 0.0,
        })

    with open(coverage_path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=list(stats[0]))
        w.writeheader()
        w.writerows(stats)

    with open(presence_path, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["label", "organism", "n_markers"] + markers)
        for taxon in taxa:
            flags = [1 if (m, taxon) in best else 0 for m in markers]
            w.writerow([taxon, organism_of.get(taxon, taxon), sum(flags)] + flags)

    print(f"hits kept {kept}; dropped by identity {dropped_id}, by coverage {dropped_cov}; "
          f"paralogy conflicts resolved {conflicts}")
    print(f"viruses with at least one marker: {len(taxa)}")
    for s in stats:
        print(f"  {s['marker']:22s} {s['n_taxa']:4d} taxa ({s['pct_taxa']:5.1f}%)  "
              f"median identity to query {s['median_identity_to_query']:5.1f}%")


if __name__ == "__main__":
    main(*sys.argv[1:7])
