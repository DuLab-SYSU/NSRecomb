[toc]

# Protocal for simulation

1. `ms` was used to generate evolutionary histories with and without recombination.
2. `seqgen`
3. apply our pipeline


# Convergent evolution

For each genealogy, 100 replicates was generated with sample size as 10. Under no recomination scenario, randomly select 2 sequences from a replicate, then calculates the number of different nts, making differing amino acid changing to another same amino acid within one substitiution. Calculate the decreased number of different nts and make the same number synonymou muataions in other nt to matain the level of divergence.

1. calculate the number of different nts.
2. chage different aa into another same aa with one substitution.
3. compensate the decreased number of nt substitutions in other codon with synonymous substitution.
