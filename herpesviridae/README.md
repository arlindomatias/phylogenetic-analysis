# herpesviridae

Phylogenetic taxonomy of the herpesvirus family **Orthoherpesviridae** (formerly *Herpesviridae*) from
NCBI RefSeq genomes, with a focus on the five human herpesviruses that carry a published association
with cognitive decline and dementia in ageing: **HSV-1**, **VZV**, **HCMV**, **EBV** and **HHV-6A**.

The question is taxonomic, not clinical: does a maximum likelihood tree built from conserved core
genes reproduce the ICTV subfamilies and genera of the family, and where do those five viruses sit in
it? The answer, for this dataset, is that the three subfamilies are recovered as clades with maximal
support, most genera are monophyletic, and the five viruses of interest are spread across all three
subfamilies rather than forming a group of their own.

## Pipeline

```
download  →  metadata                                           NCBI Virus RefSeq genomes + proteomes
          →  queries  →  homology  →  orthologs                 8 HSV-1 core proteins as DIAMOND queries
                                   →  align (per marker)        MAFFT --auto, then trimAl
                                   →  concat                    partitioned supermatrix
                                   →  tree                      IQ-TREE 2, UFBoot + SH-aLRT
                                   →  taxonomy  →  plots  →  report
```

Orthologues are found by **homology**, not by annotation text: RefSeq protein names are inconsistent
across the family (only 59 of 126 viruses have a protein literally named *DNA polymerase catalytic
subunit*), so the HSV-1 reference proteins are used as DIAMOND blastp queries and the best hit per
virus is taken, subject to identity and query coverage floors.

Monophyly is tested on **bipartitions**, so the result does not depend on where the tree is rooted —
which matters here, because there is no usable outgroup inside the alignable range.

## Running

```bash
snakemake --cores 6 --use-conda      # everything, from download to reports/report.html
snakemake --cores 6 -n               # dry run
```

The only input is the network: `scripts/download.sh` fetches the data from NCBI Virus with the
`datasets` CLI. All parameters are in `config/config.yaml`.

## Layout

```
herpesviridae/
├── Snakefile                 # download → … → report (target of rule all)
├── config/config.yaml        # taxon, marker genes, thresholds, model, highlighted viruses, colors
├── env/                      # one conda environment per rule (--use-conda)
├── scripts/                  # all code (bash for I/O, python for analysis and figures)
├── data/                     # download, marker FASTAs, alignments, supermatrix (not versioned)
├── results/                  # tables, DIAMOND hits, tree files, monophyly tests
├── plots/                    # figures (PNG + PDF)
└── reports/report.html       # self-contained report
```

## Data

NCBI Virus, RefSeq reference genomes of the family (138 genomes, 126 viruses), downloaded with
`datasets download virus genome taxon Herpesviridae --refseq`. Nothing is redistributed in this
repository; the workflow re-downloads it.
