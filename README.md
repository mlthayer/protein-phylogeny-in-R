Workflow A: Protein phylogeny in R (generic template)

This repository is a reusable template to build a protein phylogenetic tree in R
from a pre-aligned amino-acid FASTA.

Input
Place your aligned FASTA here:

- `Data/aligned_proteins.fa`

All sequences must be the same length (aligned).

Run
Open the repository in RStudio and run:

```r
source("scripts/01_build_tree_from_aligned_fasta.R")

Workflow B: Start from an unaligned FASTA (alignment + tree)

Use this workflow if you have **unaligned protein sequences** and want the repository to:
1) align them, then
2) build a phylogenetic tree.

Input
Place your unaligned protein FASTA here:

- `Data/proteins.fa`

FASTA headers can be anything, but they should be unique and descriptive.

Run
From the repository root in RStudio:

```r
source("scripts/00_align_then_build_tree.R")
