"""Confront the ML tree with the ICTV taxonomy of the family.

Usage: taxonomy_check.py <config.yaml> <tree.treefile> <taxa_metadata.csv> \
           <monophyly_subfamily.csv> <monophyly_genus.csv> <human_herpesviruses.csv> \
           <rooted.nwk> <taxonomy_summary.csv>

For every subfamily and every genus with at least two sampled viruses, the script asks
whether its members form a clade. The test is done on bipartitions (splits), so the answer
does not depend on where the tree is rooted: a group is monophyletic when the set of its
members, or its complement, is exactly the leaf set on one side of some branch.

Branch support is read from the IQ-TREE node labels (SH-aLRT / ultrafast bootstrap).
Groups that fail the test get the smallest clade that contains them, with the intruding
viruses listed, which is what makes a taxonomic conflict interpretable.
"""
import csv
import sys

import dendropy
import yaml

RANKS = ["subfamily", "genus"]


def read_metadata(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def split_support(label):
    """IQ-TREE internal label 'SH-aLRT/UFboot' -> (alrt, ufboot) as floats or ''."""
    if not label:
        return "", ""
    parts = str(label).split("/")
    out = []
    for p in parts[:2]:
        try:
            out.append(float(p))
        except ValueError:
            out.append("")
    while len(out) < 2:
        out.append("")
    return out[0], out[1]


def leaf_labels(node):
    return {lf.taxon.label for lf in node.leaf_iter()}


def test_groups(tree, all_taxa, groups):
    """For each group, say whether it is a clade on the unrooted tree and with what support."""
    edges = []
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            continue
        leaves = leaf_labels(node)
        if 1 < len(leaves) < len(all_taxa):
            edges.append((leaves, node))

    rows = []
    for name, members in sorted(groups.items()):
        members = set(members)
        row = {"taxon": name, "n_taxa": len(members)}
        if len(members) < 2:
            row.update({"monophyletic": "single taxon", "support_alrt": "", "support_ufboot": "",
                        "smallest_clade_size": len(members), "n_intruders": 0, "intruders": ""})
            rows.append(row)
            continue

        match = None
        for leaves, node in edges:
            if leaves == members or (all_taxa - leaves) == members:
                match = node
                break
        if match is not None:
            alrt, ufboot = split_support(match.label)
            row.update({"monophyletic": "yes", "support_alrt": alrt, "support_ufboot": ufboot,
                        "smallest_clade_size": len(members), "n_intruders": 0, "intruders": ""})
        else:
            mrca = tree.mrca(taxon_labels=sorted(members))
            clade = leaf_labels(mrca) if mrca is not None else all_taxa
            intruders = sorted(clade - members)
            alrt, ufboot = split_support(mrca.label if mrca is not None else "")
            row.update({"monophyletic": "no", "support_alrt": alrt, "support_ufboot": ufboot,
                        "smallest_clade_size": len(clade), "n_intruders": len(intruders),
                        "intruders": "; ".join(intruders)})
        rows.append(row)
    return rows


def sister_taxa(tree, label):
    """Closest relative(s) of a tip: the other tips of its parent node."""
    node = tree.find_node_with_taxon_label(label)
    if node is None or node.parent_node is None:
        return ""
    sisters = leaf_labels(node.parent_node) - {label}
    if not sisters:                            # polytomy-free trees always have a sister, but be safe
        gp = node.parent_node.parent_node
        sisters = (leaf_labels(gp) - {label}) if gp is not None else set()
    return "; ".join(sorted(sisters))


def main(config_path, tree_path, metadata_path, sub_path, genus_path, humans_path,
         rooted_path, summary_path):
    config = yaml.safe_load(open(config_path))
    highlight = config.get("highlight", {})
    meta = read_metadata(metadata_path)

    tree = dendropy.Tree.get(path=tree_path, schema="newick", preserve_underscores=True)
    tree.encode_bipartitions()
    all_taxa = {lf.taxon.label for lf in tree.leaf_node_iter()}

    # one row per taxon; when a virus has several RefSeq genomes, prefer the record that
    # carries host information
    by_label = {}
    for r in meta:
        kept = by_label.get(r["label"])
        if kept is None or (not kept["host"] and r["host"]):
            by_label[r["label"]] = r
    missing = sorted(all_taxa - set(by_label))
    if missing:
        raise SystemExit(f"ERROR: tips absent from the metadata table: {', '.join(missing)}")

    outputs = {}
    counts = {}
    for rank, out_path in zip(RANKS, [sub_path, genus_path]):
        groups = {}
        for label in all_taxa:
            name = by_label[label][rank]
            if name and name != "unclassified":
                groups.setdefault(name, set()).add(label)
        rows = test_groups(tree, all_taxa, groups)
        outputs[rank] = rows
        with open(out_path, "w", newline="") as out:
            w = csv.DictWriter(out, fieldnames=["rank"] + list(rows[0]))
            w.writeheader()
            for r in rows:
                w.writerow({"rank": rank, **r})
        testable = [r for r in rows if r["monophyletic"] in ("yes", "no")]
        counts[rank] = {
            "groups": len(rows),
            "testable": len(testable),
            "monophyletic": sum(1 for r in testable if r["monophyletic"] == "yes"),
            "singletons": sum(1 for r in rows if r["monophyletic"] == "single taxon"),
        }
        print(f"{rank}: {counts[rank]['monophyletic']}/{counts[rank]['testable']} groups with "
              f"2+ sampled viruses are monophyletic "
              f"({counts[rank]['singletons']} single-virus groups not testable)")
        for r in rows:
            if r["monophyletic"] == "no":
                shown = r["intruders"].split("; ")
                more = f" (+{len(shown) - 8} more)" if len(shown) > 8 else ""
                print(f"  NOT monophyletic: {r['taxon']} ({r['n_taxa']} taxa) — smallest clade has "
                      f"{r['smallest_clade_size']}, intruders: "
                      f"{'; '.join(shown[:8])}{more}")

    # rooted copy for plotting: midpoint rooting, a neutral choice for a family whose
    # root position is not settled by an outgroup here
    rooted = dendropy.Tree.get(path=tree_path, schema="newick", preserve_underscores=True)
    rooted.reroot_at_midpoint(update_bipartitions=True, suppress_unifurcations=True)
    rooted.ladderize(ascending=False)
    rooted.write(path=rooted_path, schema="newick", suppress_rooting=True,
                 unquoted_underscores=True, suppress_internal_node_labels=False)

    # human herpesviruses: the viruses of this family isolated from Homo sapiens
    human_rows = []
    for label in sorted(all_taxa):
        r = by_label[label]
        if r["host"] == "Homo sapiens" or r["organism"].startswith("Human herpesvirus") \
                or r["organism"].startswith("Human "):
            human_rows.append({
                "label": label,
                "organism": r["organism"],
                "accession": r["accession"],
                "subfamily": r["subfamily"],
                "genus": r["genus"],
                "species": r["species"],
                "genome_length": r["genome_length"],
                "cognitive_decline_note": highlight.get(r["organism"], ""),
                "closest_relative_in_tree": sister_taxa(rooted, label),
            })
    with open(humans_path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=list(human_rows[0]))
        w.writeheader()
        w.writerows(human_rows)
    print(f"human-derived viruses in the tree: {len(human_rows)} "
          f"({sum(1 for r in human_rows if r['cognitive_decline_note'])} highlighted)")

    summary = [
        {"metric": "taxa in tree", "value": len(all_taxa)},
        {"metric": "marker genes", "value": len(config["markers"])},
        {"metric": "subfamilies sampled", "value": counts["subfamily"]["groups"]},
        {"metric": "subfamilies monophyletic", "value":
            f"{counts['subfamily']['monophyletic']}/{counts['subfamily']['testable']}"},
        {"metric": "genera sampled", "value": counts["genus"]["groups"]},
        {"metric": "genera testable (2+ viruses)", "value": counts["genus"]["testable"]},
        {"metric": "genera monophyletic", "value":
            f"{counts['genus']['monophyletic']}/{counts['genus']['testable']}"},
        {"metric": "human-derived viruses", "value": len(human_rows)},
    ]
    with open(summary_path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=["metric", "value"])
        w.writeheader()
        w.writerows(summary)


if __name__ == "__main__":
    main(*sys.argv[1:9])
