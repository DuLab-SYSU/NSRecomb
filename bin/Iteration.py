import BlastTools
import Verify
from Bio import SeqIO
from Decorator import execute_time
from multiprocessing import Pool
import json
import os
import argparse


def window_i(start, end, query, backbone,
             db, negative_seqidlist,
             bootstrap_threshold=0.8):
    query_i = query[start: end]
    backbone_i = BlastTools.get_backbone(query_i, backbone)
    result = []
    try:
        hits = BlastTools.get_hits(query_i,
                                   db=db,
                                   negative_seqidlist=negative_seqidlist)
        for hit in hits:
            align_seqs = Verify.multiseqalign(query_i, backbone_i, hit)
            ds_qb, ds_qs = Verify.kaks(*align_seqs)
            bs_value = Verify.bootstrap(*align_seqs)
            result.append([query_i.id, backbone_i.id, hit.id, bs_value, ds_qb, ds_qs])
            if (bs_value > bootstrap_threshold) and (ds_qb > ds_qs):
                print("Fragment: %s-%s, (%s-%s-%s),\
                       BootStrap Value:  %s,\
                       dS(Q-B): %s,\
                       dS(Q-S): %s"
                      % (start, end,
                         query_i.id, backbone_i.id, hit.id,
                         bs_value, ds_qb, ds_qs))
        return (start, result)
    except Exception as e:
        print(str(e))
        print("Start site ", start)
        return (start, ['NA'])


def get_bounds(sequence_length, window_size, step):
    if sequence_length < window_size:
        return 1
    else:
        n_step = (sequence_length - window_size) // step
        return n_step


@execute_time
def parallel(cpu_num, func, n_step, query, backbone,
             window_size=501,
             db='db/test',
             negative_seqidlist='q.acc'):
    print('# Number of CPUs:', cpu_num)
    p = Pool(cpu_num)
    res_l = []
    for i in range(0, n_step):
        res = p.apply_async(func, args=(i*3, i*3+window_size, query, backbone,
                                        db, negative_seqidlist, ))
        res_l.append(res)
    p.close()
    p.join()
    results = {item.get()[0]: item.get()[1] for item in res_l}
    return results


def main(out='log', _backbone=None):
    query = SeqIO.read(QUERY, 'fasta')
    if not _backbone:
        try:
            backbone = BlastTools._get_backbone(query,
                                                db=DATABASE,
                                                negative_seqidlist=NEGATIVE_SEQIDLIST)
        except Exception:
            raise("No backbone sequences in the database, which may due to the\
                  long divergent time between your query sequence and the target\
                  database.\nPlease change your database or specify the backbone\
                  by your own with switch -b.")
    else:
        backbone = SeqIO.read(_backbone, 'fasta')

    sequence_length = len(query)
    n_step = get_bounds(sequence_length, WINDOW_SIZA, STEP)
    print("# Sequence Length: %s, Task counts: %s" % (sequence_length, n_step))

    results = parallel(NUM_CPUs, window_i, n_step+1, query, backbone,
                       window_size=WINDOW_SIZA,
                       db=DATABASE,
                       negative_seqidlist=NEGATIVE_SEQIDLIST)
    with open(out, 'w') as w:
        w.write(json.dumps(results, indent=4, sort_keys=True))


def unit_test(start, _backbone=None):
    query = SeqIO.read(QUERY, 'fasta')
    if not _backbone:
        backbone = BlastTools._get_backbone(query,
                                            db=DATABASE,
                                            negative_seqidlist=NEGATIVE_SEQIDLIST)
    else:
        backbone = SeqIO.read(_backbone, 'fasta')

    sequence_length = len(query)
    n_step = get_bounds(sequence_length, WINDOW_SIZA, STEP)
    print("# Sequence Length: %s, Task counts: %s" % (sequence_length, n_step))

    test_start_point = start
    window_i(test_start_point, test_start_point+WINDOW_SIZA,
             query, backbone,
             DATABASE, NEGATIVE_SEQIDLIST,
             bootstrap_threshold=0.8)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Software for identifying recombination.")
    parser.add_argument('--query', '-q', help='Query file path', required=True)
    parser.add_argument('--backbone', '-b', help='Backbone file path', required=False)
    parser.add_argument('--negative_seqidlist', '-n', help='Negative sequence id file, binary format', required=True)
    parser.add_argument('--database', '-db', help='Database generate by makeblastdb', required=True)
    parser.add_argument('--window_size', '-ws', help='Window size', required=False, default=501, type=int)
    parser.add_argument('--num_cpus', '-nc', help='Number of CPUs', required=False, default=10, type=int)
    parser.add_argument('--output', '-o', help='Output file path', required=False, default='log')

    args = parser.parse_args()

    NUM_CPUs = args.num_cpus

    QUERY = args.query
    BACKBONE = args.backbone
    NEGATIVE_SEQIDLIST = args.negative_seqidlist
    DATABASE = args.database
    WINDOW_SIZA = args.window_size
    STEP = 3
    OUTPUT = args.output

    if not os.path.exists("../result"):
        os.makedirs("../result")

    if not BACKBONE:
        main(OUTPUT)
    else:
        main(OUTPUT, _backbone=BACKBONE)

    # unit_test(1164)
