import tempfile
import subprocess
import numpy as np
from io import StringIO
from Bio import SeqIO


def multiseqalign(query_i, backbone_i, sbjct_i):
    ref_seq = tempfile.NamedTemporaryFile(mode='w+t')
    other_seq = tempfile.NamedTemporaryFile(mode='w+t')
    SeqIO.write(query_i, ref_seq, "fasta")
    SeqIO.write([backbone_i, sbjct_i], other_seq, "fasta")
    ref_seq.flush()
    ref_seq.seek(0)
    other_seq.flush()
    other_seq.seek(0)

    mafft_client = subprocess.Popen(
        'mafft --auto --keeplength --addfragments %s %s'
        % (other_seq.name, ref_seq.name),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    stdout, _ = mafft_client.communicate()
    msa = SeqIO.parse(StringIO(stdout.decode('utf-8')), "fasta")
    ref_seq.close()
    other_seq.close()
    align_query_i, align_backbone_i, align_sbjct_i = list(msa)
    align_query_i.annotations = query_i.annotations
    align_backbone_i.annotations = backbone_i.annotations
    align_sbjct_i.annotations = sbjct_i.annotations
    return align_query_i, align_backbone_i, align_sbjct_i


def _deal_nan(x):
    if x == 'NA':
        return 0
    elif x == '-nan':
        return 1
    else:
        return float(x)


def _parse_kaks_result(file_contents):
    row1, row2 = file_contents.split('\n')[1: 3]
    _, ks1 = row1.split('\t')[2: 4]
    _, ks2 = row2.split('\t')[2: 4]
    trans_ = map(_deal_nan, [ks1, ks2])
    return list(trans_)


def kaks(query_i, backbone_i, sbjct_i, method='LWL'):
    axt = tempfile.NamedTemporaryFile('w+t')
    axt.write('%s-%s\n' % (query_i.id, backbone_i.id))
    axt.write('%s\n%s\n\n' % (str(query_i.seq), str(backbone_i.seq)))
    axt.write('%s-%s\n' % (query_i.id, sbjct_i.id))
    axt.write('%s\n%s\n' % (str(query_i.seq), str(sbjct_i.seq)))
    axt.flush()
    axt.seek(0)
    # print(axt.read())
    # axt.seek(0)

    kk = tempfile.NamedTemporaryFile('w+t')

    kaks_client = subprocess.Popen(
        'KaKs_Calculator -i %s -o %s -m %s' % (axt.name, kk.name, method),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    _ = kaks_client.communicate()
    # print(stdout.decode('utf-8'), stderr.decode('utf-8'))
    kk.flush()
    kaks_result = kk.read()
    kk.close()
    print(kaks_result)
    return _parse_kaks_result(kaks_result)


def _get_identity(seq1, seq2):
    count_ident = sum([1 for i, j in zip(seq1, seq2) if i == j])
    return count_ident / len(seq1)


def bootstrap(query_i, backbone_i, sbjct_i, sample_size=1000):
    length_ = len(query_i)
    idx = np.random.randint(0, length_, size=(sample_size, length_))
    bs_q = np.array(list(str(query_i.seq)))[idx]
    bs_b = np.array(list(str(backbone_i.seq)))[idx]
    bs_s = np.array(list(str(sbjct_i.seq)))[idx]
    count = 0
    for bs_qi, bs_bi, bs_si in zip(bs_q, bs_b, bs_s):
        qb_ident = _get_identity(bs_qi, bs_bi)
        qs_ident = _get_identity(bs_qi, bs_si)
        if qs_ident > qb_ident:
            count += 1
    return count / sample_size


if __name__ == "__main__":
    import BlastTools

    query = SeqIO.read('../data/RmYN02_S.fasta', 'fasta')
    backbone = BlastTools._get_backbone(query, db='../db/all-cov',
                                        negative_seqidlist='../data/RmYN02.acc')

    query_i = query[1002: 1503]
    backbone_i = BlastTools.get_backbone(query_i, backbone)
    hits = BlastTools.get_hits(query_i,
                               db='../db/all-cov',
                               negative_seqidlist='../data/RmYN02.acc',
                               outfmt=5,
                               task='blastn',
                               qcov_hsp_perc=70)
    print("*************** Preparation *********************")
    print("Query:\n%s\n" % query)
    print('Backbone: ')
    print(backbone, '\n')
    print("*************** Iteration *********************")
    print("query_i:\n%s\n" % query_i)
    print("Backbone_i: ")
    print(backbone_i, '\n')
    print('# of best hits: ', len(hits))
    for hit in hits:
        print(hit)

    align_seqs = multiseqalign(query_i, backbone_i, hits[0])
    print('\n')
    print("****************** MSA *************************")
    for msa_i in align_seqs:
        print(msa_i, len(msa_i))

    print('\n')
    print("************ KaKs_Calculator ****************")
    ds1, ds2 = kaks(*align_seqs)
    print('dS of query_i and backbone_i: ', ds1)
    print('dS of query_i and subject_i: ', ds2)
    # kaks_bio(*align_seqs)

    print('\n')
    print("************ BootStrap ****************")
    bs_value = bootstrap(*align_seqs)
    print('BootStrap Value: ', bs_value)
