# protocal for simulation

1. Simulate a tree (genealogy) under the coalescent (without recombination) using `ms`
2. Evolve amino acid sequences from a common ancestor along the tree using `seq-gen`. If insertions and/or deletions are required, we use `INDELible` instead.
3. Denerate recombinant sequences from two or more randomly chosen sequences in the dateset, with breakpoint chosen uniformly at random along the genome. The parent sequences are removerd from the dataset.
