import random
from io import StringIO
from Bio import SeqIO
# from Bio import Phylo
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import simulation_paml


def full_mutate(codon, base_table=None, codon_table=None):
    if codon_table is None:
        codon_table = CodonTable.unambiguous_dna_by_name["Standard"]
    if base_table is None:
        base_table = ['A', 'T', 'C', 'G']

    if not isinstance(codon, str):
        raise TypeError("codon must be a string")

    try:
        amino_acid = codon_table.forward_table[codon]
    except KeyError:
        raise RuntimeError(
            "Unknown codon detected ({codon})."
        )
    # print("Original codon: ", codon)
    # print("Original amino acid: ", amino_acid)
    target = set()
    for idx, base_i in enumerate(codon):
        for base_j in base_table:
            if base_i != base_j:
                # print("Mutate: %s -> %s" % (base_i, base_j))
                tmp_ = list(codon)
                tmp_[idx] = base_j
                new_codon = "".join(tmp_)
                if new_codon not in codon_table.stop_codons:
                    new_amino_acid = codon_table.forward_table[new_codon]
                    if new_amino_acid != amino_acid:
                        target.add((new_codon, new_amino_acid))
    return target


def _cal_full_mutate(seq):
    if not isinstance(seq, Seq):
        raise TypeError("seq must be Seq or SeqRecord")
    if len(seq) % 3 != 0:
        raise ValueError("seq must be codon sequence")

    rf_table = []
    for i in range(0, len(seq), 3):
        # rf = i // 3
        targets = full_mutate(str(seq[i:i+3]))
        rf_table.append((i, targets))
    return rf_table


def find_mutatble_site(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("The length of seq1 and seq2 must be same.")

    rf_table1 = _cal_full_mutate(seq1.seq)
    rf_table2 = _cal_full_mutate(seq2.seq)

    # intersections = []
    targets = []
    for i, j in zip(rf_table1, rf_table2):
        intersection = i[1].intersection(j[1])
        if intersection:
            # intersections.append((i[0], intersection))
            targets.append((i[0], random.choice(list(intersection))))
    return targets


def target_mutate(seq, mutataion):
    mutable_seq = seq.seq.tomutable()
    for idx, (codon, _) in mutataion:
        mutable_seq[idx: idx+3] = codon
    new_seq = mutable_seq.toseq()
    return SeqRecord(new_seq, id=seq.id)


def covergent_evolution(seq1, seq2):
    candicant_targets = find_mutatble_site(seq1, seq2)
    # print("There %s candicant sites." % len(candicant_targets))
    # print(candicant_targets)

    target_num = len(candicant_targets) // 5
    targets = random.sample(candicant_targets, k=target_num)

    new_seq1 = target_mutate(seq1, targets)
    new_seq2 = target_mutate(seq2, targets)
    return new_seq1, new_seq2


if __name__ == '__main__':
    replicates_1 = simulation_paml.parser_paml(
        "/home/zeng/python_work/nCov/simulation/mc3_1.paml")
    replicates_2 = simulation_paml.parser_paml(
        "/home/zeng/python_work/nCov/simulation/mc3_2.paml")
    replicates_3 = simulation_paml.parser_paml(
        "/home/zeng/python_work/nCov/simulation/mc3_3.paml")

    print(len(replicates_2), len(replicates_2[0]), len(replicates_2[0][0]))

    replicates = []
    for x, y, z in zip(replicates_1, replicates_2, replicates_3):
        y1, y2, y3, y4 = y
        new_y1, new_y3 = covergent_evolution(y1, y3)
        full_seq1 = x[0] + new_y1 + z[0]
        full_seq2 = x[1] + y2 + z[1]
        full_seq3 = x[2] + new_y3 + z[2]
        full_seq4 = x[3] + y4 + z[3]
        replicates.append((full_seq1, full_seq2, full_seq3, full_seq4))

    # q, b, s, _ = replicates[0]
    # n_step = (len(q) - 501) // 3
    # for i in range(n_step):
    #     print(simulation_paml.one_window(i*3, i*3+501, q, b, s))
    # print(simulation_paml.one_window(3429, 3429+501, q, b, s))
    for idx, replicate_i in enumerate(replicates):
        query, backbone, subject, _ = replicate_i
        n_step = (len(query) - 501) // 3

        simulation_paml.main(n_step, query, backbone, subject, log='../simulation/s3/test2_log_%s' % idx)
