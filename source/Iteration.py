import json
import os

from source import BlastTools
from source import Verify
from source.Decorator import execute_time

from Bio import SeqIO
from multiprocessing import Pool


def window_codon_bootstrap(start, end, query, backbone, db, negative_seqidlist, qcov_hsp_perc=70, bootstrap_threshold=0.8):
    query_i = query[start: end]
    backbone_i = BlastTools.get_window_backbone(query_i, backbone, qcov_hsp_perc)

    result = []
    hits = BlastTools.get_hits(query_i, db, negative_seqidlist, qcov_hsp_perc)
    for hit in hits:
        triple_align = Verify.multiseqalign(query_i, backbone_i, hit)
        bs_ident, bs_ds = Verify.bootstrap2(*triple_align)
        
        result.append([query_i.id, backbone_i.id, hit.id, bs_ident, bs_ds])
        if (bs_ident > bootstrap_threshold) and (bs_ds > bootstrap_threshold):
            print("Fragment: %s-%s, (%s-%s-%s), BootStrap Value:  %s, %s" % (start, end, query_i.id, backbone_i.id, hit.id, bs_ident, bs_ds))

    return (start, result)


def window_site_bootstrap(start, end, query, backbone, db, negative_seqidlist, qcov_hsp_perc=70, bootstrap_threshold=0.8):
    query_i = query[start: end]
    backbone_i = BlastTools.get_window_backbone(query_i, backbone, qcov_hsp_perc)

    result = []
    hits = BlastTools.get_hits(query_i, db, negative_seqidlist, qcov_hsp_perc)
    for hit in hits:
        triple_align = Verify.multiseqalign(query_i, backbone_i, hit)
        bs_ident = Verify.bootstrap(*triple_align)
        ds1, ds2 = Verify._kaks(*triple_align)
        result.append([query_i.id, backbone_i.id, hit.id, bs_ident, ds1, ds2])
        if (bs_ident > bootstrap_threshold) and (ds1 > ds2):
            print("Fragment: %s-%s, (%s-%s-%s), BootStrap Value:  %s, dS(Q-B): %s, dS(Q-S): %s" % (start, end, query_i.id, backbone_i.id, hit.id, bs_ident, ds1, ds2))

    return (start, result)


@execute_time
def parallel(func, cpu_num, tasks, window_size, query, backbone, db, negative_seqidlist):
    p = Pool(cpu_num)
    res_l = []
    for start, end in tasks:
        res = p.apply_async(func, args=(start, end, query, backbone, db, negative_seqidlist, ))
        res_l.append(res)
    p.close()
    p.join()
    results = {int(res.get()[0]): res.get()[1] for res in res_l}
    return results


def _get_bounds(sequence_length, window_size, step):
    if sequence_length < window_size:
        n_step =  1
    else:
        n_step = (sequence_length - window_size) // step + 1

    tasks = []
    for i in range(0, n_step):
        start = i * step
        end = i * step + window_size
        tasks.append([start, end])
    return tasks


def main(out, query, backbone, db, negative_seqidlist, cpu_num, window_size, step, codon_bootstrap=False):
    sequence_length = len(query)
    tasks = _get_bounds(sequence_length, window_size, step)
    print("Sequence Length: %s, # of Tasks: %s" % (sequence_length, len(tasks)))
    print('# of CPUs:', cpu_num)
    
    if codon_bootstrap:
        func = window_codon_bootstrap
    else:
        func = window_site_bootstrap
    results = parallel(func, cpu_num, tasks, window_size, query, backbone, db, negative_seqidlist)
    with open(out, 'w') as w:
        w.write(json.dumps(results, indent=4, sort_keys=True))


if __name__ == '__main__':
    import BlastTools
    import Verify
    from Decorator import execute_time

    BlastTools.BLAST_PATH = '../package/ncbi-blast-2.11.0+/bin'
    Verify.MAFFT_PATH = '../package/mafft-linux64'
    Verify.KAKS_PATH = '../package/KaKs_Calculator2.0/bin'
    
    QUERY = SeqIO.read('../example/SARS-CoV-2.S.fasta', 'fasta')
    NEGATIVE_SEQIDLIST  = '../example/SARS-CoV-2.acc'
    DATABASE = '../db/all_cov_unify'

    # BACKBONE = SeqIO.read('../data/SARS-CoV-2.backbone', 'fasta')
    BACKBONE = BlastTools.get_gene_backbone(QUERY, DATABASE, NEGATIVE_SEQIDLIST)
    # print(BACKBONE)

    
    WINDOW_SIZA = 501
    STEP = 6
    NUM_CPUs = 10
    OUTPUT = '../result/' + "log"
    CODON = False

    if not os.path.exists("../result"):
        os.makedirs("../result")

    start = 2160
    res = window_site_bootstrap(start, start+WINDOW_SIZA, QUERY, BACKBONE, DATABASE, NEGATIVE_SEQIDLIST)
    print(res)


    # main(OUTPUT, QUERY, BACKBONE, DATABASE, NEGATIVE_SEQIDLIST, NUM_CPUs, WINDOW_SIZA, STEP, codon_bootstrap=CODON)

    