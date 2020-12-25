import json
from io import StringIO
from Bio import SeqIO
from Bio.codonalign.codonseq import CodonSeq
from Bio.codonalign.codonseq import cal_dn_ds
import Verify
from Decorator import execute_time
from multiprocessing import Pool


def parser_paml(_in, num_seq=4):
    replicates = []
    with open(_in) as f:
        count = 0
        sample = []
        for line in f:
            if line != '\n':
                sample.append(line)
                count += 1
            if count == num_seq+1:
                sample_ = SeqIO.parse(StringIO("".join(sample)), 'phylip')
                replicates.append(list(sample_))
                count = 0
                sample = []
    return replicates


def cal_ds(seq1, seq2):
    new_seq1 = CodonSeq.from_seq(seq1)
    new_seq2 = CodonSeq.from_seq(seq2)
    ds = cal_dn_ds(new_seq1, new_seq2, method='LWL85')[1]
    return ds


def one_window(start, end, query, backbone, subject):
    query_i = query[start: end]
    backbone_i = backbone[start: end]
    subject_i = subject[start: end]
    try:
        ds1 = cal_ds(query_i.seq, backbone_i.seq)
        ds2 = cal_ds(query_i.seq, subject_i.seq)
    except Exception:
        print("ERROR,\n%s\n%s\n%s\n" % (query_i.seq, backbone_i.seq, subject_i.seq))
    bs_value = Verify.bootstrap(query, backbone, subject)
    return start, [query_i.id, backbone_i.id, subject_i.id, bs_value, ds1, ds2]


@execute_time
def main(n_step, query, backbone, subject, log='test_log'):
    p = Pool(10)
    res_l = []
    for i in range(0, n_step):
        res = p.apply_async(one_window, args=(i*3, i*3+501, query, backbone, subject, ))
        res_l.append(res)
    p.close()
    p.join()
    results = {item.get()[0]: item.get()[1] for item in res_l}
    with open(log, 'w') as w:
        w.write(json.dumps(results, indent=4, sort_keys=True))


if __name__ == '__main__':
    replicates_1 = parser_paml("/home/zeng/python_work/nCov/simulation/mc2_1.paml")
    replicates_2 = parser_paml("/home/zeng/python_work/nCov/simulation/mc2_2.paml")
    replicates_3 = parser_paml("/home/zeng/python_work/nCov/simulation/mc2_3.paml")
    replicates = []
    for x, y, z in zip(replicates_1, replicates_2, replicates_3):
        sample = []
        for x_i, y_i, z_i in zip(x, y, z):
            sample.append(x_i+y_i+z_i)
        else:
            replicates.append(sample)

    print(len(replicates), len(replicates[0]), len(replicates[0][0]))

    # print(one_window(3, 504, query, backbone, subject))
    # main((9000-501)//3, *replicates[0][:3], log='test_log')
    for idx, replicate_i in enumerate(replicates):
        query, backbone, subject, _ = replicate_i
        n_step = (len(query) - 501) // 3

        main(n_step, query, backbone, subject, log='../simulation/s2/test1_log_%s' % idx)
