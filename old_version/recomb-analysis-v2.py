# -*- coding: utf-8 -*-

import os, sys
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from collections import Counter, defaultdict
import numpy as np
from multiprocessing import Pool, Process


class Blastor():
    def __init__(self, _query,  _db, _mask):
        self.query = _query
        self.db = _db
        self.mask = _mask
        self.query_seq = self._read_query(_query)
        self.query_len = len(self.query_seq)

    def _read_query(self, _in):
        with open(_in)as f:
            seq = ''
            for line in f:
                if not line.startswith(">"):
                    seq += line.strip()
        return seq

    def _blastn(self, query_loc):
        command = "blastn -query %s -query_loc %s -db %s \
                   -num_alignments 10 -task blastn -max_hsps 1\
                   -outfmt '6 qseqid qstart qend sseqid staxid sacc stitle scomname sstart send pident' \
                   -negative_seqidlist %s"  % (self.query, query_loc, self.db, self.mask)
        result = os.popen(command)  
        res = result.read()
        return res

    def one_blast(self, query_loc='all'):
        if query_loc == 'all':
            result = self._blastn('1-%s' % self.query_len)
        else:
            result = self._blastn(query_loc)
        print('The best alignment for %s is %s' % (query_loc, result.split("\n")[0].split("\t")))
        return result.split("\n")[0].split("\t")

    def window_blast(self, log_file, window_size=1000):
        fw = open("%s" % log_file, "w")
        results = []
        for i in range(1, self.query_len - window_size):
            query_loc = "%s-%s" % (i, i+window_size)
            result = self.one_blast(query_loc)
            if result == ['']:
                results.append([None]*11)
                fw.write("\t".join(['None']*11) + "\n")
            else:
                results.append(result)
                fw.write("\t".join(result) + "\n")
        fw.close()
        return results


def generate_fragments(qstarts, stitles, pidents):
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
            print("%s : %s-%s, length: %s" % (title, start, end, end-start))
            start = qstarts[i]
            title = stitles[i]
    print("%s : %s-%s, length: %s" % (title, start, qstarts[-1], qstarts[-1]-start))
    fragments[title].append((start, qstarts[-1]))
    return fragments, x_list, y_list


def plot_fig(qstarts, stitles, pidents, colormap, out_name):
    fragments, x_list, y_list = generate_fragments(qstarts, stitles, pidents)

    fig, ax = plt.subplots(figsize=(16, 4))
    fig.subplots_adjust(0.025, 0.1, 0.85, 0.9)
    #ax.yaxis.set_visible(False)
    for k, v in fragments.items():
        hl = defaultdict(list)
        for item in v:
            hl['xmin'].append(item[0])
            hl['xmax'].append(item[1])
            hl['y'].append(105)
        plt.hlines(y=hl['y'], xmin=hl['xmin'], xmax=hl['xmax'], \
                   label=k, color=colormap[k], linewidth=8)
        plt.scatter(x=x_list[k], y=y_list[k] , color=colormap[k], s=2)
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)
    fig.savefig("/home/zeng/python_work/bioinfo/nCov/result/%s.jpg" % out_name)


def main(_in_file, _mask, window_size, outplot=False):
    out_name = "%s-%s" % (_in_file.split('.')[0], window_size)
    test = Blastor('/home/zeng/Desktop/%s' % _in_file, 'all-cov/all-cov', '/home/zeng/Desktop/%s' % _mask)
    b = test.window_blast("result/%s.log" % out_name, window_size=int(window_size))
    if outplot:
        _, qstarts, _, _, _, _, stitles, snames, _, _, pidents = zip(*b)
        qstarts = [int(x) for x in qstarts if x != 'None']
        pidents = [float(x) for x in pidents if x != 'None']
        stitles = [x for x in stitles if x != 'None']
        snames = [x for x in snames if x != 'None']

        colormap = {}
        colors = list(mcolors.TABLEAU_COLORS.values())
        top_9 = [i[0] for i in sorted(list(Counter(stitles).most_common(9)), key=lambda x: x[1], reverse=True)]

        stitles2 = list(map(lambda x: 'other' if x not in top_9 else x, stitles))
        colormap = {top_9[i]: colors[i] for i in range(len(top_9))}
        colormap['other'] = colors[9]

        plot_fig(qstarts, stitles2, pidents, colormap, out_name)


if __name__ == "__main__":
    _, _in_file, _mask = sys.argv
    # _in_file = 'RaTG13.fasta'
    # _mask = '2019-ncov-ratg.idlist'
    # window_size = 500

    # qaccs, qstarts, qends, sseqids, staxids, saccs, stitles, snames, sstarts, sends, pidents

    process_list = []
    tasks = [1000]
    for task in tasks:
        p = Process(target=main, args=(_in_file, _mask, task, False, ))
        p.start()
        process_list.append(p)
    
    for p in process_list:
        p.join()

#    main(_in_file, _mask, task, 500)