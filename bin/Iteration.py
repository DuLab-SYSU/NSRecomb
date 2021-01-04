import json
import os
import argparse

import BlastTools
import Verify
from Decorator import execute_time

from Bio import SeqIO
from multiprocessing import Pool


def window(start, end, query, backbone, db, negative_seqidlist, qcov_hsp_perc=70, bootstrap_threshold=0.8):
    query_i = query[start: end]
    backbone_i = BlastTools.get_window_backbone(query_i, backbone, qcov_hsp_perc)

    result = []
    try:
        hits = BlastTools.get_hits(query_i, db, negative_seqidlist, qcov_hsp_perc)
        for hit in hits:
            triple_align = Verify.multiseqalign(query_i, backbone_i, hit)
            bs_ident, bs_ds = Verify.bootstrap2(*triple_align)
            
            result.append([query_i.id, backbone_i.id, hit.id, bs_ident, bs_ds])
            if (bs_ident > bootstrap_threshold) and (bs_ds > bootstrap_threshold):
                print("Fragment: %s-%s, (%s-%s-%s), BootStrap Value:  %s, %s"
                      % (start, end, query_i.id, backbone_i.id, hit.id, bs_ident, bs_ds))

        return (start, result)
    except Exception as e:
        print(str(e))
        print("Start site ", start)
        return (start, ['NA']*5)


@execute_time
def parallel(func, cpu_num, n_step, window_size,
             query, backbone, db, negative_seqidlist,
             qcov_hsp_perc=70, bootstrap_threshold=0.8):
    print('# Number of CPUs:', cpu_num)
    p = Pool(cpu_num)
    res_l = []
    for i in range(0, n_step):
        start = i * 3
        end = i * 3 + window_size
        res = p.apply_async(func, args=(start, end, query, backbone, db, negative_seqidlist, qcov_hsp_perc, bootstrap_threshold, ))
        res_l.append(res)
    p.close()
    p.join()
    results = {res.get()[0]: res.get()[1] for res in res_l}
    return results


def _get_bounds(sequence_length, window_size, step):
    if sequence_length < window_size:
        return 1
    else:
        n_step = (sequence_length - window_size) // step
        return n_step + 1


def main(out='log', _backbone=None):
    query = SeqIO.read(QUERY, 'fasta')
    sequence_length = len(query)
    n_step = _get_bounds(sequence_length, WINDOW_SIZA, STEP)
    print("# Sequence Length: %s, Task counts: %s" % (sequence_length, n_step))
    
    if _backbone:
        backbone = SeqIO.read(_backbone, 'fasta')
    else:
        try:
            backbone = BlastTools.get_gene_backbone(query, db=DATABASE, negative_seqidlist=NEGATIVE_SEQIDLIST)
        except Exception:
            raise("""
                  No backbone sequences in the database, which may due to the
                  long divergent time between your query sequence and the target
                  database.\n
                  Please change your database or specify the backbone
                  by your own with switch -b.
                  """)

    results = parallel(window, NUM_CPUs, n_step, WINDOW_SIZA,
                       query, backbone, db=DATABASE, negative_seqidlist=NEGATIVE_SEQIDLIST, qcov_hsp_perc=70, bootstrap_threshold=0.8)
    with open(out, 'w') as w:
        w.write(json.dumps(results, indent=4, sort_keys=True))


if __name__ == '__main__':    
    QUERY = SeqIO.read('../data/SARS-CoV-2.ORF1ab.fasta', 'fasta')
    BACKBONE = SeqIO.read('../data/SARS-CoV-2.backbone', 'fasta')
    DATABASE = '../db/all-cov'
    NEGATIVE_SEQIDLIST  = '../data/SARS-CoV-2.acc'
    WINDOW_SIZA = 501
    STEP = 3
    NUM_CPUs = 10
    
    
    start = 3
    res = window(start, start+WINDOW_SIZA, QUERY, BACKBONE, DATABASE, NEGATIVE_SEQIDLIST)
    print(res)
