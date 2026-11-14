"""Figures: annotated ML tree, marker coverage, genome size by subfamily.

Usage: plot_tree.py <config.yaml> <rooted.nwk> <taxa_metadata.csv> <marker_coverage.csv> <plots_dir>

The tree is drawn directly with matplotlib (rectangular phylogram): tip labels and branches
coloured by ICTV subfamily, filled dots on branches with ultrafast bootstrap >= 95, and the
human herpesviruses linked to cognitive decline in ageing called out on the right margin.
"""
import csv
import os
import sys

import dendropy
import matplotlib
import yaml

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

UFBOOT_CUTOFF = 95


def save(fig, plots_dir, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(plots_dir, f"{name}.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {name}.png / {name}.pdf")


def read_csv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def layout(tree):
    """Assign (x = root distance, y = tip order) to every node of a rooted tree."""
    x, y = {}, {}
    tip_counter = [0]
    for node in tree.preorder_node_iter():
        edge = node.edge_length or 0.0
        x[node] = (x[node.parent_node] + edge) if node.parent_node is not None else 0.0
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            y[node] = tip_counter[0]
            tip_counter[0] += 1
        else:
            ys = [y[c] for c in node.child_node_iter()]
            y[node] = sum(ys) / len(ys)
    return x, y


# <figure name="tree_subfamily">
def plot_tree(tree, meta_by_label, config, plots_dir):
    colors = config["colors"]
    highlight = config.get("highlight", {})
    x, y = layout(tree)
    n_tips = sum(1 for _ in tree.leaf_node_iter())

    fig, ax = plt.subplots(figsize=(13, max(8, 0.19 * n_tips)))

    def subfamily_of(node):
        subs = {meta_by_label[lf.taxon.label]["subfamily"] for lf in node.leaf_iter()}
        return subs.pop() if len(subs) == 1 else None

    for node in tree.preorder_node_iter():
        sub = subfamily_of(node)
        color = colors.get(sub, "#555555") if sub else "#999999"
        if node.parent_node is not None:
            ax.plot([x[node.parent_node], x[node]], [y[node], y[node]],
                    color=color, lw=1.0, solid_capstyle="round", zorder=2)
        children = list(node.child_node_iter())
        if children:
            ys = [y[c] for c in children]
            ax.plot([x[node], x[node]], [min(ys), max(ys)], color=color, lw=1.0, zorder=2)
            try:                                  # IQ-TREE label: 'SH-aLRT/UFboot'
                ufboot = float(str(node.label).split("/")[1])
            except (AttributeError, IndexError, ValueError):
                ufboot = None
            if ufboot is not None and ufboot >= UFBOOT_CUTOFF:
                ax.plot([x[node]], [y[node]], marker="o", ms=2.6, color="#333333", zorder=3)

    xmax = max(x.values())
    for leaf in tree.leaf_node_iter():
        label = leaf.taxon.label
        row = meta_by_label[label]
        organism = row["organism"][:1].upper() + row["organism"][1:]  # NCBI mixes capitalisation
        note = highlight.get(organism)
        color = colors.get(row["subfamily"], colors.get("unclassified", "#777777"))
        ax.plot([x[leaf], xmax * 1.015], [y[leaf], y[leaf]], color="#dddddd", lw=0.4, zorder=1)
        ax.text(xmax * 1.02, y[leaf], organism, fontsize=5.6, va="center",
                color=colors["highlight"] if note else color,
                fontweight="bold" if note else "normal")
        if note:
            ax.text(xmax * 1.52, y[leaf], "◀ " + note, fontsize=5.6, va="center",
                    color=colors["highlight"], fontweight="bold")

    bar = round(xmax / 10, 2) or xmax / 10
    ax.plot([0, bar], [-2.2, -2.2], color="#333333", lw=1.4)
    ax.text(bar / 2, -3.4, f"{bar} subst./site", ha="center", fontsize=7, color="#333333")

    handles = [Line2D([], [], color=colors[s], lw=2, label=s)
               for s in ("Alphaherpesvirinae", "Betaherpesvirinae", "Gammaherpesvirinae")]
    handles.append(Line2D([], [], color="#333333", marker="o", ms=4, lw=0,
                          label=f"ultrafast bootstrap ≥ {UFBOOT_CUTOFF}%"))
    handles.append(Line2D([], [], color=colors["highlight"], lw=2,
                          label="linked to cognitive decline in ageing"))
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=7.5)

    ax.set_xlim(-xmax * 0.02, xmax * 2.25)
    ax.set_ylim(-5, n_tips + 1)
    ax.axis("off")
    ax.set_title("Orthoherpesviridae: maximum likelihood tree of 8 core genes\n"
                 "midpoint-rooted; tips coloured by ICTV subfamily", fontsize=11, loc="left")
    save(fig, plots_dir, "tree_subfamily")
# </figure>


# <figure name="marker_coverage">
def plot_coverage(coverage, plots_dir):
    rows = sorted(coverage, key=lambda r: int(r["n_taxa"]))
    names = [r["marker"] for r in rows]
    taxa = [int(r["n_taxa"]) for r in rows]
    ident = [float(r["median_identity_to_query"]) for r in rows]

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
    axes[0].barh(names, taxa, color="#3399FF")
    for i, v in enumerate(taxa):
        axes[0].text(v + 1, i, str(v), va="center", fontsize=8)
    axes[0].set_xlabel("viruses with an orthologue")
    axes[0].set_title("Marker coverage", fontsize=10, loc="left")

    axes[1].barh(names, ident, color="#009933")
    for i, v in enumerate(ident):
        axes[1].text(v + 0.4, i, f"{v:.0f}%", va="center", fontsize=8)
    axes[1].set_xlabel("median identity to the HSV-1 query (%)")
    axes[1].set_title("Divergence from the query", fontsize=10, loc="left")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=8)
    fig.tight_layout()
    save(fig, plots_dir, "marker_coverage")
# </figure>


# <figure name="genome_size_by_subfamily">
def plot_genome_size(meta_by_label, config, plots_dir):
    colors = config["colors"]
    groups = {}
    for row in meta_by_label.values():
        try:
            kb = int(row["genome_length"]) / 1000.0
        except (TypeError, ValueError):
            continue
        groups.setdefault(row["subfamily"], []).append((kb, row["organism"]))
    order = [s for s in ("Alphaherpesvirinae", "Betaherpesvirinae", "Gammaherpesvirinae")
             if s in groups] + [s for s in sorted(groups) if s not in
                                ("Alphaherpesvirinae", "Betaherpesvirinae", "Gammaherpesvirinae")]
    highlight = config.get("highlight", {})

    fig, ax = plt.subplots(figsize=(8, 4))
    rng = __import__("random").Random(config["seed"])
    for i, sub in enumerate(order):
        values = [v for v, _ in groups[sub]]
        ax.boxplot(values, positions=[i], widths=0.55, showfliers=False,
                   medianprops={"color": "#333333"})
        for kb, organism in groups[sub]:
            is_hl = organism in highlight
            ax.plot(i + rng.uniform(-0.17, 0.17), kb, marker="o",
                    ms=6 if is_hl else 3.4,
                    color=colors["highlight"] if is_hl else colors.get(sub, "#777777"),
                    alpha=1.0 if is_hl else 0.55, zorder=3 if is_hl else 2)
            if is_hl:
                ax.annotate(organism.replace("Human ", "H. "), (i, kb), fontsize=6.5,
                            xytext=(10, 0), textcoords="offset points", va="center",
                            color=colors["highlight"], fontweight="bold")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([s.replace("herpesvirinae", "-") for s in order], fontsize=9)
    ax.set_ylabel("genome size (kb)")
    ax.set_title("Genome size by subfamily (RefSeq genomes in the tree)", fontsize=10, loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, plots_dir, "genome_size_by_subfamily")
# </figure>


def main(config_path, tree_path, metadata_path, coverage_path, plots_dir):
    config = yaml.safe_load(open(config_path))
    os.makedirs(plots_dir, exist_ok=True)
    tree = dendropy.Tree.get(path=tree_path, schema="newick", preserve_underscores=True)
    meta = read_csv(metadata_path)
    tips = {lf.taxon.label for lf in tree.leaf_node_iter()}
    meta_by_label = {}
    for r in meta:
        if r["label"] in tips:
            meta_by_label.setdefault(r["label"], r)
    plot_tree(tree, meta_by_label, config, plots_dir)
    plot_coverage(read_csv(coverage_path), plots_dir)
    plot_genome_size(meta_by_label, config, plots_dir)


if __name__ == "__main__":
    main(*sys.argv[1:6])
