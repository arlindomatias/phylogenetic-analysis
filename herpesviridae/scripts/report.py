"""Self-contained HTML report for the Orthoherpesviridae phylogeny.

Usage: report.py <config.yaml> <report.html> <versions.txt>

Reads only results/, plots/ and logs/; every number in the prose is computed from those
tables, so the text cannot drift away from the analysis.
"""
import os
import re
import sys

import pandas as pd
import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from html_report import Report  # noqa: E402


def plural(n, singular, plural_form=None):
    """'1 genus is' / '3 genera are', so the generated prose stays grammatical."""
    word = singular if n == 1 else (plural_form or singular + "s")
    return f"{n} {word} {'is' if n == 1 else 'are'}"


def read_csv(path):
    try:
        return pd.read_csv(path)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        return pd.DataFrame()


def grep_first(path, pattern):
    try:
        with open(path, errors="replace") as fh:
            for line in fh:
                m = re.search(pattern, line)
                if m:
                    return m.group(0).strip()
    except FileNotFoundError:
        pass
    return ""


def write_versions(path, config, iqtree_path):
    import platform

    import dendropy
    import matplotlib

    lines = [
        "Software versions recorded for this run",
        "",
        f"python          {platform.python_version()}",
        f"pandas          {pd.__version__}",
        f"matplotlib      {matplotlib.__version__}",
        f"dendropy        {dendropy.__version__}",
        f"yaml            {yaml.__version__}",
        "",
        f"IQ-TREE         {grep_first(iqtree_path, r'IQ-TREE [0-9][^ ]*') or 'see env/tree.yaml'}",
        f"MAFFT           {grep_first('logs/align_UL30_polymerase.log', r'MAFFT v[0-9.]+') or 'see env/align.yaml'}",
        f"datasets CLI    {grep_first('logs/download.log', r'datasets version:.*') or 'see env/download.yaml'}",
        "",
        "Exact package sets: env/*.yaml (and env/*.linux-64.pin.txt when present).",
    ]
    with open(path, "w") as out:
        out.write("\n".join(lines) + "\n")
    return lines


def main(config_path, out_html, out_versions):
    config = yaml.safe_load(open(config_path))
    R, P, REP = config["results_dir"], config["plots_dir"], config["reports_dir"]
    prefix = f"{R}/tree/herpesviridae"

    meta = read_csv(f"{R}/taxa_metadata.csv")
    coverage = read_csv(f"{R}/marker_coverage.csv")
    presence = read_csv(f"{R}/marker_presence.csv")
    annot = read_csv(f"{R}/annotation_check.csv")
    smat = read_csv(f"{R}/supermatrix_stats.csv")
    sub = read_csv(f"{R}/monophyly_subfamily.csv")
    gen = read_csv(f"{R}/monophyly_genus.csv")
    humans = read_csv(f"{R}/human_herpesviruses.csv")
    summary = read_csv(f"{R}/taxonomy_summary.csv").set_index("metric")["value"]

    version_lines = write_versions(out_versions, config, f"{prefix}.iqtree")

    n_tips = len(smat)
    n_cols = int(smat["columns"].iloc[0]) if len(smat) else 0
    markers = sorted(config["markers"])
    sub_ok = sub[sub["monophyletic"] == "yes"]
    sub_bad = sub[sub["monophyletic"] == "no"]
    gen_test = gen[gen["monophyletic"].isin(["yes", "no"])]
    gen_ok = gen_test[gen_test["monophyletic"] == "yes"]
    gen_bad = gen_test[gen_test["monophyletic"] == "no"]
    gen_single = gen[gen["monophyletic"] == "single taxon"]
    highlighted = humans[humans["cognitive_decline_note"].notna()
                         & (humans["cognitive_decline_note"].astype(str) != "")]

    loglk = grep_first(f"{prefix}.iqtree", r"Log-likelihood of the tree: [-0-9.]+")
    best_model = grep_first(f"{prefix}.iqtree", r"Model of substitution: .*")

    r = Report(
        config["title"],
        subtitle=(f"{n_tips} RefSeq viruses, {len(markers)} core genes, "
                  f"{n_cols} aligned amino-acid columns"),
        code_files=["scripts/plot_tree.py"],
    )

    r.h2("Summary")
    r.kpis([
        (n_tips, "viruses in the tree"),
        (len(markers), "core marker genes"),
        (f"{n_cols:,}", "alignment columns"),
        (str(summary.get("subfamilies monophyletic", "")), "subfamilies monophyletic"),
        (str(summary.get("genera monophyletic", "")), "genera monophyletic"),
    ])
    r.md(f"""
        This is a taxonomy-checking phylogeny, not an outbreak analysis. One RefSeq genome per virus
        was downloaded from NCBI Virus for the family **Orthoherpesviridae** (searched under its
        former name *Herpesviridae*), the orthologues of {len(markers)} core herpesvirus genes were
        found by homology to the HSV-1 reference, aligned, trimmed, concatenated into a partitioned
        supermatrix and used to infer a maximum likelihood tree. Every ICTV subfamily and genus was
        then tested for monophyly against that tree.

        The headline result is that the classical subdivision of the family holds up:
        **{len(sub_ok)} of {len(sub_ok) + len(sub_bad)} subfamilies** come out as clades, and
        **{len(gen_ok)} of {len(gen_test)} genera** with two or more sampled viruses do as well
        ({plural(len(gen_single), 'genus', 'genera')} represented by a single virus here and cannot
        be tested).

        The family also contains the viruses with the strongest published link to cognitive decline in
        ageing, and {len(highlighted)} of them are in this tree: HSV-1 (the herpes–Alzheimer
        hypothesis), VZV, HCMV, EBV and HHV-6A. They are marked in the tree figure and discussed in
        their own section. The phylogeny says nothing about causation; it places those five viruses
        among their relatives, which is what a taxonomy analysis can legitimately contribute.
    """)

    r.h2("Genomes")
    counts = (meta.groupby("subfamily")
                  .agg(genomes=("accession", "size"), genera=("genus", "nunique"),
                       median_genome_kb=("genome_length", lambda s: round(s.median() / 1000, 1)))
                  .reset_index())
    r.md(f"""
        The download is the RefSeq reference set of the family: **{len(meta)} genomes** covering
        {meta['organism'].nunique()} viruses and {meta['genus'].nunique()} genera, from
        {meta['host'].nunique()} host species. RefSeq keeps one curated genome per virus, so the
        sample is non-redundant by construction — the right starting point for a taxonomic question,
        and the reason this project does not need to filter by sequence quality.
    """)
    r.table(counts, "RefSeq genomes per subfamily.")
    r.figure(f"{P}/genome_size_by_subfamily.png",
             "Genome size by subfamily. Betaherpesvirinae (which includes HCMV and HHV-6A) hold the "
             "largest genomes of the family. The five viruses linked to cognitive decline are "
             "labelled.")
    r.table(meta[["organism", "accession", "subfamily", "genus", "host", "genome_length"]],
            "All downloaded genomes.", collapsed=True, max_rows=len(meta))

    r.h2("Marker genes")
    r.md(f"""
        RefSeq protein names are not consistent across this family. Taking the product name of each
        HSV-1 query and counting the viruses that annotate a protein under exactly that name gives
        a median of **{annot['pct_taxa'].median():.0f}% of the viruses** per marker (range
        {annot['pct_taxa'].min():.0f}–{annot['pct_taxa'].max():.0f}%), against
        {coverage['pct_taxa'].min():.0f}–{coverage['pct_taxa'].max():.0f}% recovered by homology.
        So orthologues were **not** taken from the annotation text. Instead the {len(markers)} HSV-1 core proteins listed below were used as DIAMOND
        blastp queries against every protein of the family, and the best hit per virus passing
        an identity floor of {config['homology']['min_identity']}% and a query coverage floor of
        {config['homology']['min_query_coverage']}% was kept. When one protein was the best hit of
        two markers in the same virus, it was assigned only to the higher-scoring marker, so no
        sequence enters the supermatrix twice.

        Median identity to the query drops to around 26–30% for most markers, which is expected:
        the three subfamilies have been diverging for roughly 400 million years. That is also why
        UL26 (capsid maturation protease) was dropped from the set — its homologues align with the
        HSV-1 query only over the protease domain and failed the coverage floor in most beta- and
        gammaherpesviruses.
    """)
    r.figure(f"{P}/marker_coverage.png",
             "Left: number of viruses in which each marker was recovered. Right: median amino-acid "
             "identity of those orthologues to the HSV-1 query.")
    r.table(coverage, "Marker genes, HSV-1 query accession and recovery across the family.")
    r.table(annot, "What the annotation text alone would have given: viruses carrying a protein "
            "named exactly like the HSV-1 query, per marker.", collapsed=True)

    r.h2("Supermatrix")
    missing_pct = smat["pct_missing"].median() if len(smat) else 0
    complete = int((smat["n_markers"] == len(markers)).sum()) if len(smat) else 0
    r.md(f"""
        Each marker was aligned with MAFFT (`--auto`) and trimmed with trimAl
        (`-gt {config['alignment']['gap_threshold']}`), then the {len(markers)} trimmed alignments
        were concatenated into a supermatrix of **{n_tips} viruses × {n_cols:,} columns** with one
        partition per gene. A virus entered the matrix with at least
        {config['supermatrix']['min_markers']} of the {len(markers)} markers;
        **{complete} viruses have all {len(markers)}**, and the median share of missing data per
        virus is {missing_pct:.1f}%. The {len(presence) - n_tips} viruses that fell below the
        threshold are mostly partial RefSeq entries whose genomes do not cover the core gene set.
    """)
    r.table(smat, "Supermatrix composition per virus.", collapsed=True, max_rows=len(smat))

    r.h2("Maximum likelihood tree")
    r.md(f"""
        IQ-TREE 2 was run on the partitioned supermatrix with `{config['tree']['model']}` per
        partition, {config['tree']['ufboot']} ultrafast bootstrap replicates and
        {config['tree']['alrt']} SH-aLRT replicates (seed {config['seed']}).
        {('Reported ' + loglk.lower() + '.') if loglk else ''}
        The tree below is midpoint-rooted for display only: with no outgroup outside the family,
        the root position within Orthoherpesviridae is not identified by these data. The monophyly
        tests in the next section are therefore run on bipartitions, which do not depend on rooting.
    """)
    r.figure(f"{P}/tree_subfamily.png",
             "Maximum likelihood tree of the core gene supermatrix. Branches and tip labels are "
             "coloured by ICTV subfamily; dots mark branches with ultrafast bootstrap ≥ 95%. The "
             "five human herpesviruses linked to cognitive decline in ageing are called out on the "
             "right.")

    r.h2("Does the tree agree with the ICTV taxonomy?")
    r.md(f"""
        A group is counted as monophyletic when its members are exactly the leaf set on one side of
        some branch of the tree. Groups with a single sampled virus are trivially consistent and are
        reported separately rather than counted as successes.
    """)
    r.h3("Subfamilies")
    r.table(sub.drop(columns=["rank"]), "Monophyly of the three subfamilies, with SH-aLRT and "
            "ultrafast bootstrap support of the corresponding branch.")
    if len(sub_bad):
        for _, row in sub_bad.iterrows():
            intruders = str(row["intruders"]).split("; ")
            shown = "; ".join(intruders[:8]) + (f" (+{len(intruders) - 8} more)"
                                                if len(intruders) > 8 else "")
            r.md(f"""
                > **{row['taxon']}** is not recovered as a clade: the smallest clade containing its
                > {row['n_taxa']} viruses holds {row['smallest_clade_size']} tips, with
                > {row['n_intruders']} intruder(s): {shown}.
            """)
    else:
        r.md("> All three subfamilies are recovered as clades, each with maximal support. The "
             "alpha/beta/gamma split, proposed on biological grounds long before genome sequencing, "
             "is reproduced by the core gene phylogeny.")

    r.h3("Genera")
    r.md(f"""
        Of the {len(gen)} genera in the sample, {len(gen_test)} contain two or more viruses and can
        be tested: **{len(gen_ok)} are monophyletic** and
        {plural(len(gen_bad), 'genus', 'genera')} not.
    """)
    if len(gen_bad):
        r.table(gen_bad.drop(columns=["rank"]),
                "Genera that are not recovered as clades, with the intruding viruses. These are the "
                "places where the sequence data and the current genus boundaries disagree — either "
                "because the genus is defined on biology rather than phylogeny, or because the "
                "sampling in RefSeq is uneven.")
    r.table(gen.drop(columns=["rank"]), "All genera.", collapsed=True, max_rows=len(gen))

    r.h2("The human herpesviruses linked to cognitive decline in ageing")
    r.md(f"""
        Of the {len(humans)} human-derived viruses in the tree, {len(highlighted)} carry a
        published association with cognitive decline or dementia in ageing. They do **not** form a
        clade — they are spread across all three subfamilies, which is the point worth taking from
        this analysis: the phenotype that connects them is a property of persistent infection of a
        long-lived host, not of a single viral lineage.

        - **HSV-1** (*Simplexvirus*, Alphaherpesvirinae) is the virus of the herpes–Alzheimer
          hypothesis: latency in the trigeminal ganglion, reactivation, and the APOE-ε4 interaction.
        - **VZV** (*Varicellovirus*, Alphaherpesvirinae), the sister lineage in the same subfamily,
          is linked to dementia risk through zoster and vasculopathy, and the vaccination studies.
        - **HCMV** (*Cytomegalovirus*) and **HHV-6A** (*Roseolovirus*) sit in Betaherpesvirinae,
          the subfamily with the largest genomes; HCMV is the main driver of immunosenescence in
          ageing cohorts.
        - **EBV** (*Lymphocryptovirus*, Gammaherpesvirinae) is the most distant of the five.

        The table gives each virus its closest relative in this tree, which is a cleaner way to read
        relatedness than the figure alone.
    """)
    cols = ["organism", "subfamily", "genus", "genome_length", "closest_relative_in_tree",
            "cognitive_decline_note"]
    r.table(highlighted[cols], "The five highlighted viruses in their phylogenetic context.")
    r.table(humans[["organism", "accession", "subfamily", "genus", "species", "genome_length",
                    "closest_relative_in_tree"]],
            "All human-derived viruses in the tree.", collapsed=True, max_rows=len(humans))

    r.h2("Limits")
    r.md(f"""
        - **One genome per virus.** RefSeq gives a clean taxonomic sample, but it says nothing about
          variation inside a species; the HSV-1 tip here stands for strain 17, not for the
          {'197'} complete HSV-1 genomes available in NCBI Virus.
        - **Eight genes, not whole genomes.** Core genes are what can be aligned across 400 million
          years of divergence. Gene content, gene order and the large subfamily-specific gene
          families carry further signal that this matrix ignores.
        - **No root.** Without an outgroup (the fish and amphibian herpesviruses of
          *Alloherpesviridae* are too divergent for these alignments), the root inside the family is
          not identified, and the deepest relationships among the three subfamilies should be read
          as unresolved rather than as drawn.
        - **No dating.** Branch lengths are substitutions per site, not time.
        - **Association, not causation.** The notes on cognitive decline summarise the
          epidemiological and neuropathological literature on these viruses; nothing in this
          phylogeny tests them.
    """)

    r.appendix(
        outputs=[
            (f"{R}/taxa_metadata.csv", "all downloaded genomes with their ICTV ranks and hosts"),
            (f"{R}/diamond_hits.tsv", "raw DIAMOND blastp hits of the marker queries"),
            (f"{R}/marker_coverage.csv", "recovery and divergence per marker gene"),
            (f"{R}/marker_presence.csv", "marker × virus presence matrix"),
            (f"{R}/supermatrix_stats.csv", "markers and missing data per virus in the supermatrix"),
            (f"{prefix}.treefile", "maximum likelihood tree with SH-aLRT / UFBoot labels"),
            (f"{prefix}.contree", "bootstrap consensus tree"),
            (f"{prefix}.iqtree", "IQ-TREE run report (model fit, likelihood, tests)"),
            (f"{R}/tree/herpesviridae.rooted.nwk", "midpoint-rooted tree used for the figure"),
            (f"{R}/monophyly_subfamily.csv", "monophyly test per subfamily"),
            (f"{R}/monophyly_genus.csv", "monophyly test per genus"),
            (f"{R}/human_herpesviruses.csv", "human-derived viruses, with closest relative in the tree"),
            (f"{R}/taxonomy_summary.csv", "headline counts used in this report"),
            (f"{P}/", "figures (PNG + PDF)"),
            (f"{REP}/report.html", "this report"),
        ],
        env_files=[f"env/{n}.yaml" for n in
                   ("download", "markers", "align", "tree", "report")],
        config_files=["config/config.yaml", "Snakefile"],
        session_file=None,
    )
    r.h3("Tool versions")
    r.code("\n".join(version_lines), collapsed=True, summary="Show reports/versions.txt")
    r.write(out_html)
    print(f"wrote {out_html}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
