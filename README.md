# phylogenetic-analysis

Phylogenetic analyses of public sequence data, one folder per dataset. Each folder is a
self-contained Snakemake workflow — download, alignment, tree inference, interpretation and a
single HTML report — that can be run from scratch with nothing but a network connection.

## Analyses

| Folder | Question | Data | Result |
|---|---|---|---|
| [`herpesviridae/`](herpesviridae/) | Does a core gene phylogeny reproduce the ICTV taxonomy of the herpesvirus family, and where do the human herpesviruses linked to cognitive decline in ageing sit in it? | 138 RefSeq genomes of *Orthoherpesviridae* (NCBI Virus) | 3/3 subfamilies and 14/15 testable genera recovered as clades; the five viruses of interest are spread across all three subfamilies |

## Running an analysis

```bash
cd herpesviridae
snakemake --cores 6 --use-conda      # everything, ending in reports/report.html
snakemake --cores 6 -n               # dry run
```

Every analysis builds its own conda environments from `env/*.yaml`, so the only prerequisites are
`snakemake` and `conda`.

## Conventions

- One folder per dataset, each with its own `Snakefile`, `config/config.yaml`, `env/`, `scripts/`,
  `results/`, `plots/` and `reports/`, plus a `CLAUDE.md` recording the analysis decisions.
- Bash for downloads and tool invocation, Python for analysis and figures.
- `results/`, `plots/` and `reports/` are versioned; `data/` is not — the workflows regenerate it.
- Every number in a report is computed from the tables in `results/`, never written by hand.
