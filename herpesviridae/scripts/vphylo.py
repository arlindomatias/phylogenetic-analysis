"""Shared helpers for the Orthoherpesviridae phylogeny workflow (standard library only).

Covers the three things every step needs: reading the NCBI FASTA/JSONL files, turning an
NCBI lineage into ICTV ranks (subfamily / genus / species), and turning an organism name
into a label that survives Newick and PHYLIP.
"""
import json
import re

# `>YP_009137105.1:1-1235 DNA polymerase catalytic subunit [organism=X] [isolate=Y]`
_HEADER = re.compile(
    r"^(?P<id>\S+)(?:\s+(?P<name>.*?))?\s*(?:\[organism=(?P<organism>[^\]]*)\])?"
    r"\s*(?:\[isolate=(?P<isolate>[^\]]*)\])?\s*$"
)


def read_fasta(path):
    """Yield (header, sequence) pairs; header is the '>' line without the '>'."""
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            elif line:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def write_fasta(path, records, width=60):
    """Write (name, sequence) pairs, wrapped at `width` characters."""
    with open(path, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")


def parse_header(header):
    """Split an NCBI protein FASTA header into id, product name, organism and isolate."""
    m = _HEADER.match(header)
    if not m:
        return {"id": header.split()[0], "name": "", "organism": "", "isolate": ""}
    d = m.groupdict()
    name = (d["name"] or "").strip()
    if name.startswith("CDS:"):            # unnamed CDS entries: `CDS:` with nothing after it
        name = name[4:].strip()
    return {
        "id": d["id"],
        "name": name,
        "organism": (d["organism"] or "").strip(),
        "isolate": (d["isolate"] or "").strip(),
    }


def protein_accession(seq_id):
    """`YP_009137105.1:1-1235` -> `YP_009137105.1`."""
    return seq_id.split(":", 1)[0]


def ictv_ranks(lineage_names):
    """Pull family / subfamily / genus / species out of an NCBI lineage name list.

    ICTV names are recognised by their suffixes (-viridae, -virinae) and by position:
    the genus is the rank right below the subfamily and the species right below the genus.
    Unassigned ranks come back as 'unclassified'.
    """
    ranks = {"family": "unclassified", "subfamily": "unclassified",
             "genus": "unclassified", "species": "unclassified"}
    # NCBI inserts holding nodes such as 'unclassified Gammaherpesvirinae'; they are not ranks,
    # and dropping them leaves the genus of those viruses correctly unassigned.
    lineage_names = [n for n in lineage_names
                     if not n.lower().startswith(("unclassified", "unassigned"))]
    for i, name in enumerate(lineage_names):
        low = name.lower()
        if low.endswith("viridae"):
            ranks["family"] = name
        elif low.endswith("virinae"):
            ranks["subfamily"] = name
            if i + 1 < len(lineage_names) and " " not in lineage_names[i + 1]:
                ranks["genus"] = lineage_names[i + 1]          # genus: one word, below subfamily
                if i + 2 < len(lineage_names):
                    ranks["species"] = lineage_names[i + 2]    # species: binomial, below genus
    return ranks


def read_report(path):
    """Read data_report.jsonl into a list of flat dicts, one per genome."""
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            virus = r.get("virus", {})
            lineage = [x["name"] for x in virus.get("lineage", [])]
            ranks = ictv_ranks(lineage)
            rows.append({
                "accession": r.get("accession", ""),
                "organism": virus.get("organismName", ""),
                "tax_id": virus.get("taxId", ""),
                "family": ranks["family"],
                "subfamily": ranks["subfamily"],
                "genus": ranks["genus"],
                "species": ranks["species"],
                "host": (r.get("host") or {}).get("organismName", ""),
                "genome_length": r.get("length", ""),
                "gene_count": r.get("geneCount", ""),
                "protein_count": r.get("proteinCount", ""),
                "isolate": (r.get("isolate") or {}).get("name", ""),
                "release_date": (r.get("releaseDate") or "")[:10],
                "lineage": "; ".join(lineage),
            })
    return rows


def label(organism):
    """Organism name -> tree label safe for Newick/PHYLIP (no spaces, commas, colons).

    The first letter is upper-cased so that NCBI entries differing only in capitalisation
    ('Anatid alphaherpesvirus 1' and 'anatid alphaherpesvirus 1' are two RefSeq genomes of
    the same virus, tax id 104388) collapse onto one taxon. Newick readers compare taxon
    labels case-insensitively anyway, so leaving them apart would break the tree parsing.
    """
    out = re.sub(r"[^A-Za-z0-9._-]+", "_", organism.strip())
    out = re.sub(r"_+", "_", out).strip("_")
    return out[:1].upper() + out[1:]
