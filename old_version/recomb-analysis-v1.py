# -*- coding: utf-8 -*-

import sys, os
import random
import re
import argparse
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

class RandomCut():
    def __init__(self, in_file='', length=0, mode='all', seg_length='9-100'):
        self.full_sequence = self.read_sequence(in_file) # 读取序列
        self.full_length = len(self.full_sequence) # 序列长度
        self.mode = mode
        self.range = self.mode_select(mode) # 读取切片
        self.sequence = self.full_sequence[self.range[0]-1: self.range[1]] # 切取序列
        self.length = len(self.sequence) # 子序列长度
        self.seg_length = self.extract_seg(seg_length) # 获取随机切片长度范围
    
    def __len__(self):
        return self.length
    
    def extract_seg(self, seg_length):
        return tuple(map(int, seg_length.split('-')))
    
    def mode_select(self, mode):
        '''mode in ['s', 'e', 'm', 'n'] or range xx-xx'''
        if mode == 's':
            return (21563, 25385)
        elif mode == 'e':
            return (26245, 26472)
        elif mode == 'm':
            return (26523, 27191)
        elif mode == 'n':
            return (28274, 29533)
        elif mode == 'orf1ab':
            return (266, 21555)
        elif mode == 'all':
            return (1, self.full_length)
        else:
            return tuple(map(int, mode.split('-')))

    def read_sequence(self, _in):
        seq = ''
        with open(_in) as f:
            for line in f:
                if not line.startswith(">"):
                    seq += line.replace("\n", "")     
        return seq

    def random_cut(self):
        len_l, len_u = self.seg_length
        len_of_seg = random.randrange(len_l, len_u+1)
        start = random.randrange(0, self.length - len_of_seg)
        end = start + len_of_seg
        return start, end, self.range[0]

    def random_cut_write(self):
        start, end, _ = self.random_cut()
        fragment, _range = self.sequence[start: end], (start, end-1)   
        with open('fragment-%s.fasta' % self.mode, 'w') as f:
            query = ">fragment[%s..%s]\n%s\n" % (_range[0], _range[1], fragment)
            f.write(query)


class Blastor():
    def __init__(self, _query, _blastdb, _mask='2697049'):
        self.query = _query
        self.blastdb = _blastdb
        self.mask = _mask

    def blastn(self, query_loc, quiet=False):
        command = "blastn -query %s -query_loc %s -db %s -task blastn -negative_taxids %s -outfmt '6 qseqid saccver stitle pident sstart send length qstart qend' -num_alignments 10" % (self.query, query_loc, self.blastdb, self.mask)
        result = os.popen(command)  
        res = result.read()
        if not quiet:
            if res:
                print(res)
            else:
                print("No hits.")
        return res
    
    def extract_result(self, query_loc, quiet=False):
        res_blast = self.blastn(query_loc, quiet=quiet)
        if res_blast:
            return res_blast.split('\n')[0]


class Leaderboard():
    def __init__(self):
        pass


def read_log(_in):
    sacs = []
    snames = []
    pidents = []
    qstarts = []
    qends = []
    with open(_in) as f:
        for line in f:
            row = line.split('\t')
            sacs.append(row[1])
            snames.append(row[2])
            pidents.append(row[3])
            qstarts.append(row[7])
            qends.append(row[8])
    return sacs, snames, pidents, qstarts, qends 


def judge_exist(_in):
    status = os.path.exists(_in)
    if status:
        os.system('rm -rf %s' % _in) 


def randomcolor(num_of_color):
    color_B = sorted(random.sample(range(0, 256), num_of_color))
    color_hex = list(map(lambda x: "#%s00FA" % str(hex(x))[-2:].replace('x', '0').upper(), color_B))
    return color_hex



def path_proc():
    current_path = os.getcwd()
    judge_exist(os.path.join(current_path, 'segment-%s.fasta' % mode))
    # os.chdir(current_path + "/tmp/")
    return current_path


def main():
    current_path = path_proc()

    nCOV = RandomCut(in_file=in_file, mode=mode, seg_length=seg_length)
    blastor = Blastor('%s' % in_file, 'coronavrius/coronavrius', mask)

    def get_sseq_id():
        start, end, offset = nCOV.random_cut()
        ac_name_pident_s_e  = blastor.extract_result('%s-%s' % (start+offset, end+offset), quiet=True)
        if ac_name_pident_s_e:
            return ac_name_pident_s_e.split("\t")
        else:
            return get_sseq_id()

    results = []
    judge_exist(os.path.join(current_path, 'log-%s' % task))
    with open(os.path.join(current_path, 'log-%s' % task), 'a') as f:
        for _ in range(iter_nums):
            a = get_sseq_id()
            results.append(a)
            print("\t".join(a))
            f.write("\t".join(a) + "\n")
    

    _, _, stitles, pidents, sstarts, sends, lengths, qstarts, qends = list(zip(*results))

    pidents = list(map(float, pidents))
    sstarts = list(map(int, sstarts))
    sends = list(map(int, sends))
    lengths = list(map(int, lengths))
    qstarts = list(map(int, qstarts))
    qends = list(map(int, qends))

    # sacs, snames, pidents, pstarts, pends = read_log(os.path.join(current_path, 'log-%s' % task))

    cmap = {}
    top_4 = [i[0] for i in Counter(stitles).most_common(4)]
    top_4_color = ['#ff0000', '#ff7400', '#009999', '#00cc00']
    for i in range(4):
        cmap[top_4[i]] = top_4_color[i]
    cmap['other'] = '#67e667'

    print(list(cmap.items()))

    stitles2 = list(map(lambda x: x if x in top_4 else 'other', stitles))

    # color_sp = randomcolor(len(set(snames)))
    # cmap = {n: color_sp[i] for i, n in enumerate(set(snames))}
    # cmap['RaTG13'] = '#FF3030'

    plt.figure(figsize=(16, 9))
    for sname, pident, qstart, qend in zip(stitles2, pidents, qstarts, qends):
        plt.hlines(pident, qstart, qend, label=sname, color=cmap[sname])
    plt.ylabel('% Identity')
    plt.xlim(xmin=min(qstarts), xmax=max(qends))
    plt.title("%s, Random Times: %s, Fragment length: %s" % (mode.upper(), iter_nums, seg_length))

    annoation_line_length = nCOV.length // 100 *2
    annoation_line_start = nCOV.length + 100
    annoation_line_end = annoation_line_length + annoation_line_start
    for i, s in enumerate(list(cmap.items())):
        annoation_line_y = 100 - i
        plt.hlines(annoation_line_y, annoation_line_start, annoation_line_end, colors=s[1])
        plt.text(annoation_line_end+50 ,annoation_line_y-0.25, s[0], color=s[1])


    plt.savefig(os.path.join(current_path, "result-%s.jpg" % task))

    #plot_colortable(cmap, "Color Table of %s" % task, emptycols=0)

if __name__ == '__main__':
    params = sys.argv[1:]
    task = params[0]
    in_file = params[1]
    mode = params[2]
    iter_nums = int(params[3])
    seg_length = params[4]
    mask = params[5]


    main()

    os.system('rm -rf fragment-*')
