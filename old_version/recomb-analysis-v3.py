# -*- coding: utf-8 -*-
import os, sys
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from collections import Counter, defaultdict
import numpy as np
from multiprocessing import Pool, Process


def split_sequence(_in, _out, _length):
    fr = open(_in, 'r')
    fw = open(_out, 'w')
    seq = ''
    for line in fr:
        if not line.startswith(">"):
            seq += line.strip()
    for i in range(0, len(seq)-_length):
        fragment = seq[i:i+_length]
        fw.write(">fragment_%s_%s\n%s\n" % (i+1, i+_length, fragment))


def blastn(_query, _mask, _out):
    command = "blastn -task blastn -max_hsps 1 -num_alignments 1 -num_threads 4 \
            -db all-cov/all-cov -query %s \
            -negative_seqidlist %s \
            -outfmt '6 qseqid sacc sstart send pident length mismatch qcovhsp stitle'" % (_query, _mask)
    result = os.popen(command)  
    res = result.read()
    with open(_out, "w") as f:
        f.write(res)


def read_log(_in):
    with open(_in) as f:
        res = f.read()
    qseqids, saccs, sstarts, sends, pidents, _, _, _, stitles = zip(*[line.split('\t') for line in res.split("\n") if line != ''])
    qstarts = list(map(lambda x: int(x.split("_")[1]), qseqids))
    qends = list(map(lambda x: int(x.split("_")[2]), qseqids))
    sstarts = list(map(lambda x :int(x), qstarts))
    sends = list(map(lambda x :int(x), qends))
    pidents = list(map(lambda x :float(x), pidents))
    return qstarts, qends, saccs, sstarts, sends, pidents, stitles


def read_fragments(qstarts, pidents, stitles):
    x_list = defaultdict(list)
    y_list = defaultdict(list)
    for qstart, pident, stitle in zip(qstarts, pidents, stitles):
        x_list[stitle].append(qstart)
        y_list[stitle].append(pident)

    fragments = defaultdict(list)
    start = qstarts[0]
    title = stitles[0]
    for i in range(1, len(stitles)):
        if stitles[i] != stitles[i-1]:
            end = qstarts[i-1]
            fragments[title].append((start, end))
            start = qstarts[i]
            title = stitles[i]
    fragments[title].append((start, qstarts[-1]))
    return fragments, x_list, y_list


def plot_fig(qstarts, stitles, pidents, colormap, out_name):
    fragments, x_list, y_list = read_fragments(qstarts, pidents, stitles)

    fig, ax = plt.subplots(figsize=(16, 4))
    fig.subplots_adjust(0.025, 0.1, 0.85, 0.9)
    for k, v in fragments.items():
        print("%s: %s" % (k, v))
        hl = defaultdict(list)
        for item in v:
            hl['xmin'].append(item[0])
            hl['xmax'].append(item[1])
            hl['y'].append(105)
        plt.hlines(y=hl['y'], xmin=hl['xmin'], xmax=hl['xmax'], \
                   label=k, color=colormap[k], linewidth=8)
        plt.scatter(x=x_list[k], y=y_list[k] , color=colormap[k], s=2)
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)
    fig.savefig("%s" % out_name)


def main(_query, _mask, window_size):
    infile = _query
    split = "%s.%s.split" % (_query, window_size)
    log =  "%s.%s.log" % (_query, window_size)
    jpg =  "%s.%s.jpg" % (_query, window_size)
    split_sequence(infile, split, window_size)
    blastn(split, _mask, log)
    qstarts, qends, saccs, sstarts, sends, pidents, stitles = read_log(log)
    qstarts = [(qstart+qend)//2 for qstart, qend in zip(qstarts, qends)]
    colormap = {}
    colors = list(mcolors.TABLEAU_COLORS.values())
    top_9 = [i[0] for i in sorted(list(Counter(stitles).most_common(9)), key=lambda x: x[1], reverse=True)]

    stitles2 = list(map(lambda x: 'other' if x not in top_9 else x, stitles))
    colormap = {top_9[i]: colors[i] for i in range(len(top_9))}
    colormap['other'] = colors[9]

    plot_fig(qstarts, stitles2, pidents, colormap, jpg)
    print("\n".join(set(saccs)))


if __name__ == '__main__':
    params = sys.argv
    workspace = '/home/zeng/Desktop/'
    os.chdir(workspace)
#    main("2019-ncov.fasta", "2019-ncov.idlist", 1000)
    main(params[1], params[2], int(params[3]))

