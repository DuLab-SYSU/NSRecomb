import random
from pathlib import Path
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.codonalign.codonseq import CodonSeq


CODON_TABLE = CodonTable.unambiguous_dna_by_name["Standard"]


def diff_nt(seq1, seq2):
    idxs = []
    for idx, (i, j) in enumerate(zip(seq1, seq2)):
        if i != j:
            idxs.append(idx)
    return idxs


def full_mutate(codon):
    base_table = ['T', 'C', 'G', 'A']
    amino_acid = CODON_TABLE.forward_table[codon]

    c1, c2, c3 = list(codon)
    new_codon1, new_codon2, new_codon3 = [], [], []
    for base in base_table:
        if c1 != base:
            new_codon1_i = ''.join([base, c2, c3])
            if new_codon1_i not in CODON_TABLE.stop_codons:
                new_codon1.append(new_codon1_i)
        if c2 != base:
            new_codon2_i = ''.join([c1, base, c3])
            if new_codon2_i not in CODON_TABLE.stop_codons:
                new_codon2.append(new_codon2_i)
        if c3 != base:
            new_codon3_i = ''.join([c1, c2, base])
            if new_codon3_i not in CODON_TABLE.stop_codons:
                new_codon3.append(new_codon3_i)

    return set(new_codon1 + new_codon2 + new_codon3)


def common_codon(codon1, codon2):
    s1 = full_mutate(codon1)
    s2 = full_mutate(codon2)
    return list(s1.intersection(s2))


def diff_codon(seq1, seq2):
    codon_seq1 = CodonSeq.from_seq(seq1)
    codon_seq2 = CodonSeq.from_seq(seq2)
    
    codon_len = codon_seq1.get_codon_num()
    idxs_synon, idxs_nonsy, idxs_same = [], [], []
    for idx in range(codon_len):
        codon_1 = codon_seq1.get_codon(idx)
        codon_2 = codon_seq2.get_codon(idx)
        aa_1 = CODON_TABLE.forward_table[codon_1]
        aa_2 = CODON_TABLE.forward_table[codon_2]
        if codon_1 != codon_2 and aa_1 != aa_2:
            idxs_nonsy.append(idx)
        elif codon_1 != codon_2 and aa_1 == aa_2:
            idxs_synon.append(idx)
        else:
            idxs_same.append(idx)
    return idxs_synon, idxs_nonsy, idxs_same


def convergent_mutate(seq1, seq2):
    idxs_synon, idxs_nonsy, idxs_same = diff_codon(seq1, seq2)
    codon_seq1 = CodonSeq.from_seq(seq1)
    codon_seq2 = CodonSeq.from_seq(seq2)

    mutableseq1 = seq1.tomutable()
    mutableseq2 = seq2.tomutable()
    for idx in idxs_nonsy:
        c1 = codon_seq1.get_codon(idx)
        c2 = codon_seq2.get_codon(idx)
        candi = common_codon(c1, c2)
        if candi:
            target_codon = random.choice(candi)
            mutableseq1[idx*3: idx*3+3] = target_codon
            mutableseq2[idx*3: idx*3+3] = target_codon
        else:
            continue

    nt_diff1 = diff_nt(seq1, seq2)
    nt_diff2 = diff_nt(mutableseq1.toseq(), mutableseq2.toseq())

    delta_nt = len(nt_diff1) - len(nt_diff2)
    # print(delta_nt)
    conpensate_idx = random.sample(idxs_same, delta_nt)
    for idx in conpensate_idx:
        codon_i = str(mutableseq1[idx*3: idx*3+3])
        aa_i = CODON_TABLE.forward_table[codon_i]
        syno_codon = [x for x in full_mutate(codon_i) if CODON_TABLE.forward_table[x] == aa_i]
        if syno_codon:
            mutableseq1[idx*3: idx*3+3] = random.choice(syno_codon)
        else:
            continue

    new_seq1 = mutableseq1.toseq()
    new_seq2 = mutableseq2.toseq()
    return new_seq1, new_seq2


def convert_replicates(in_path, out_path=None):
    replicate = list(SeqIO.parse(in_path, 'fasta'))
    idx1, idx2 = sorted(random.sample(range(len(replicate)), 2))

    seq1 = replicate.pop(idx2)
    seq2 = replicate.pop(idx1)
    
    seq1_1, seq1_2, seq1_3 = seq1[0: 900], seq1[900: 1800], seq1[1800: 3000]
    seq2_1, seq2_2, seq2_3 = seq2[0: 900], seq2[900: 1800], seq2[1800: 3000]

    new_seq1_2, new_seq2_2 = convergent_mutate(seq1_2.seq, seq2_2.seq)    
    new_seq1 = seq1_1.seq + new_seq1_2 + seq1_3.seq
    new_seq2 = seq2_1.seq + new_seq2_2 + seq2_3.seq
    new_seq1 = SeqRecord(new_seq1, id=seq1.id)
    new_seq2 = SeqRecord(new_seq2, id=seq2.id)

    replicate.extend([new_seq1, new_seq2])

    SeqIO.write(replicate, out_path, "fasta")
    

def generate_convergent_dataset():
    original_path = Path.cwd()
    norecomb_dir = original_path.joinpath('norecombinate')
    for target_dir_i in [x for x in norecomb_dir.iterdir() if x.is_dir()]:
        in_dir = target_dir_i.joinpath('Results')
        replicates = list(in_dir.glob('seque*'))
        for replicate in replicates:
            out_dir = original_path.joinpath('convergent2', target_dir_i.name)
            if not out_dir.exists():
                out_dir.mkdir(parents=True)
            out_path = str(out_dir.joinpath(replicate.name))
            convert_replicates(str(replicate), out_path=out_path)

generate_convergent_dataset()
