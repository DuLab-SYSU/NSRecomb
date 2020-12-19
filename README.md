# Recombination Analysis for SARS-CoV-2

This project was designed to explore the recombintaion events during the evolution of *Coronavirus*.

The final purpose is to generate the ARG process.

## Framework

- Main.py
- Decorator.py
- Prepare.py
- BlastTools.py
- Verify.py
- Iteration.py
- Merge.py
- PValue.py
- Plot.py

### BlastTools

This module contains 3 functions, `_get_backbone`, `get_backbone` and `get_hits`, respectively.
`get_backbone` was used to get the corresponding backbone fragment from the backbone sequence.
`get_hits` was used to get these hits with the highest `score`

### Verify


### Iteration


### Merge


### PValue


### Plot


## Design

输入 Query 序列

## Dependency

- blast+ 2.11.0
  - query, fasta file
  - database, local fasta file
  - negative sequence id lists
  - backbone file, fasta
- KaKs_Calculator 2.0
- Mafft 7.475
- Biopython
- numpy
- matplotlib
