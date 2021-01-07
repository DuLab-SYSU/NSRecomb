# Recombination Analysis for SARS-CoV-2

This project was designed to explore the recombintaion events during the evolution of *Coronavirus*.

## Dependency

- blast+ 2.11.0
  - query, fasta file, in data directory
  - database, local fasta file, in data directory --> db directory
  - negative sequence id lists, in data directory
  - backbone file, fasta, in data directory
- KaKs_Calculator 2.0, intergrated in program
- Mafft 7.475, need to install and specify in the PATH
- Biopython
- numpy
- pandas
- matplotlib

## Modules

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

This module was designed to identify the potential recombintaion regions based on the similarity analysis using the representative strain as the query against all other candidate strains.  

Taking the closest lineage from the genomic tree as the backbone (user given or auto set with `get_gene_backbone` function).  

For a specified window covering certain numbers codons, the most similar lineage was identified based on blastn algorithm (`get_hits` function). The coressponding backbone fragments was identified with `get_window_backbone` function.  

### Verify

This module was designed to test the robustness of the recombination and to exclude the possibility of noise and convergent evolution based on synomous mutation distance and bootstrap method.  

Bootstapping was performed by constructing random sequences (query-backbone-similar triple (`multiseqalign` function)) with the same number of codons through randomly selecting nucleotide triples from the real sequences, and repeated 1000 times. (`bootstrap2` function)  

### Iteration

This module was the main module to identify recombination regions. Intergrating search phage and verify phage in a `window` function and able to utilize multi cpus.

### Merge

This module was designed to analysis the log file from the recombination analysis, and merge series of confirmed overlapping recombiation windows to longer regions.  

ATTEINTION: This is only experimental code and can only be used to reproduce the reslts of this study. The general verison will be released in the future verison.  

### PValue

This module will be designed to given the probaility of the identified recombination region. And will be released in the future verison.  

### Main

Main program body.  

### Interface

GUI interface.  

## Install & Usage

``` shell
git clone https://github.com/DuLab-SYSU/NSRcomb.git
cd NSRcomb/bin
# CLI
python Main.py -db all-cov -q ../data/SARS-CoV-2.ORF1ab.fasta -n ../data/SARS-CoV-2.acc -b ../data/SARS-CoV-2.backbone -o ../result/SARS-CoV-2.ORF1ab.result
#GUI
python Interface.py
```
