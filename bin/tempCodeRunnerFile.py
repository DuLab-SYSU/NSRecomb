

def cal_ds(seq1, seq2):
    new_seq1 = CodonSeq.from_seq(seq1)
    new_seq2 = CodonSeq.from_seq(seq2)
    ds = cal_dn_ds(new_seq1, new_seq2, method='LWL85')[1]
    return ds