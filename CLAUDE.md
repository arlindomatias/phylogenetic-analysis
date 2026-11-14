# phylogenetic-analysis

Repository `arlindomatias/phylogenetic-analysis` (public; local clone `~/github/phylogenetic-analysis`).
Phylogenetic analyses of public sequence data, **one folder per dataset**. It started as the folder
`viral_phylogeny` of the private repository `~/github/private`, with history kept, and follows the
conventions of `~/github/proteomics` and `~/github/bulk-rna-seq`.

## Layout

```
phylogenetic-analysis/
├── CLAUDE.md                 # this file: repository-wide conventions
├── README.md                 # index of the analyses
├── .gitignore                # applies to every analysis folder (data/, .snakemake/, logs/, …)
└── <dataset>/                # one self-contained analysis
    ├── CLAUDE.md             #   decisions specific to that dataset
    ├── README.md
    ├── Snakefile             #   full workflow, target: reports/report.html
    ├── config/config.yaml    #   every parameter
    ├── env/                  #   one conda environment per rule (--use-conda)
    ├── scripts/              #   all code
    ├── data/                 #   downloads and intermediates (not versioned)
    └── results/ plots/ reports/
```

Current analyses: `herpesviridae/` (taxonomy of *Orthoherpesviridae* from RefSeq genomes).

## Conventions

- **Run Snakemake from inside the analysis folder**, never from the repository root; all paths in a
  `config.yaml` are relative to its own folder.
- **Languages:** bash for download and tool invocation, Python for analysis and figures. Avoid R
  unless a step genuinely needs a Bioconductor package.
- **Scripts take explicit CLI arguments** (never the `snakemake` object), so each one can be run by
  hand. Every rule redirects stdout/stderr to `logs/<rule>.log`.
- **`scripts/html_report.py` is a shared file**, byte-identical to the copies in `~/github/proteomics`
  and `~/github/bulk-rna-seq`. Do not edit it for a project-specific need.
- **Reports compute their own numbers** from the tables in `results/`; nothing is hard-coded in the
  prose, so the text cannot drift away from the analysis.
- **Writing language is English** everywhere in the repository, including commit messages.

## Adding an analysis

Copy the structure of an existing folder, keep the rule naming (`download → … → report`), add a row
to the table in `README.md`, and record in the folder's own `CLAUDE.md` the decisions that the code
alone does not explain — why a marker, threshold or taxon was excluded, and what the known
limitations are.
