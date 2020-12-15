# -*- coding: utf-8 -*-
# Requirement:
#   blsatdb: all-cov/all-cov
#   blastdb_version: 5 
#   query: xxxx.fasta
#   mask information: idlist.xxxx (one id one line)

import os, sys
import json
import time
import numpy as np
from collections import Counter, defaultdict
from multiprocessing import Pool, Process
from multiprocessing.dummy import Pool as ThreadPool


is_quiet = True
#window_length = 251
step = 3
background_length = 750
pidebt_threshold = 80
ks_threshold = 0.17

# Query sequences class
# attr  : acc, title, sequence, length
# method: get_partial
class Query():
    def __init__(self, _in=None, acc=None, title=None, sequence=None, file=None):
        if _in:
            self.acc, self.title, self.sequence = self._read_seq(_in)
            self.file = _in
        else:
            self.acc = acc
            self.title = title
            self.sequence = sequence
            self.file = file

    def _read_seq(self, _in):
        with open(_in) as f:
            for line in f:
                if line.startswith('>'):
                    des = line[1:].rstrip()
                    acc, *title = des.split()
                    seq = ''
                else:
                    seq += line.rstrip()
        return acc, ' '.join(title), seq
    
    def get_partial(self, start, end):
        start2 = 1 if start < 1 else start
        end2 = len(self) if end > len(self) else end
        seq = self.sequence[start2-1: end2]
        file = '%s_%s_%s.fasta' % (self.acc, start, end)
        with open(file, 'w') as f:
            f.write(">%s:%s..%s %s\n%s\n" % (self.acc, start, end, self.title, seq))
        return Query(acc='%s:%s..%s' % (self.acc, start, end), title=self.title, sequence=seq, file=file)
    
    def __len__(self):
        return len(self.sequence)
                    

# result class
# attr: acc, title, start, end, pident, sequence
class Subject():
    """the class of one hsp in blastn result"""
    def __init__(self, acc=None, title=None, start=None, end=None, pident=None, score=None, sequence=None):
        self.acc = acc
        self.title = title
        self.start = start
        self.end = end
        self. pident = pident
        self.score = score
        self.sequence = sequence
    
    def to_file(self):
        self.file = '_'.join([self.acc, self.start, self.end]) + ".fasta"
        with open(self.file, 'w') as f:
            f.write(">%s:%s..%s %s\n%s\n" % (self.acc, self.start, self.end, self.title, self.sequence))

    def to_background_file(self, length=background_length):
        start = self.start - length
        end = self.end + length
        if start < 1:
            start = 1
        _range = '%s-%s' % (start, end)
        file = '_'.join([self.acc, str(self.start-length), str(self.end+length)]) + '.fasta'
        self.background_file = file
        commond = "blastdbcmd -db all-cov/all-cov -entry %s -range %s -out %s" % (self.acc, _range, file)
        os.system(commond)


class Result():
    """list of Subject class"""
    def __init__(self, _in):
        self.result = self._parse(_in)

    def _parse(self, _in):
        subjects = []
        for i in _in:
            _, _, _, sacc, sstart, send, pident, score, sseq, stitle = i
            subject = Subject(sacc, stitle, int(sstart), int(send), float(pident), int(score), sseq)
            subjects.append(subject)
        return subjects

    def top1(self, quiet=is_quiet):
        result = sorted(self.result, key=lambda x: x.score, reverse=True)
        top_score = result[0].score
        top1 = [i for i in result if i.score == top_score]
        if not quiet:
            print("Top HSPs:")
            for i in top1:
                print('\t', i.acc, i.pident, i.title)
        return top1
    
    def extract(self, acc):
        return [hsp for hsp in self.result if hsp.acc == acc]


# blastn
def blastn(query, mask):
    """
    accept a query fasta file, a mask file
    return Result class => list of Subject class
    """
    commond = "blastn -query %s -db all-cov/all-cov -max_hsps 1 \
                      -num_alignments 100 -outfmt '6 qacc qstart qend sacc sstart send pident score sseq stitle' \
                      -negative_seqidlist %s -task blastn \
                      2> /dev/null" % (query.file, mask)

    res = [x.rstrip().split('\t') for x in os.popen(commond).readlines()]
    result = Result(res)
    return result


# blastdbcmd
def background_identity(query, subject, length=background_length, quiet=is_quiet):
    # subject.to_background_file()
    start = subject.start-length if subject.start-length > 0 else 1
    end = subject.end+length
    commond ="blastdbcmd -db all-cov/all-cov -entry %s -range %s-%s | blastn -subject %s -outfmt '6 pident' -max_hsps 1" % (subject.acc, start, end, query.file)
    try:
        pident = float(os.popen(commond).read().strip())
    except:
        pident = 0
    if not quiet:
        print('Background identity: ', pident)
    return pident


# KaKs
# two sequece alignment
def align(query, subject):
    with open('%s.to_align' % (query.file), 'w') as f:
        f.write(">%s\n%s\n" % (query.acc, query.sequence))
        f.write(">%s:%s..%s\n%s\n" % (subject.acc, subject.start, subject.end, subject.sequence))
    os.system('mafft --quiet %s.to_align > %s.align' % (query.file, query.file))
    os.remove('%s.to_align' % (query.file))

# fasta => AXT
def fasta2AXT(_in):
    idlist = []
    data = {}
    with open(_in) as f:
        for line in f:
            if line.startswith(">"):
                des = line[1:].rstrip()
                acc = des.split()[0]
                idlist.append(acc)
                data[acc] = ''
            else:
                data[acc] += line.rstrip()
    os.remove(_in)
    with open(_in.replace('.align', '.axt'), 'w') as f:
        f.write('-'.join(idlist) + '\n')
        for i in idlist:
            f.write(data[i] + '\n')

def kaks(query, subject, quiet=is_quiet):
    align(query, subject)
    fasta2AXT('%s.align' % query.file)

    # KaKs_Calculator
    commond = "KaKs_Calculator -i %s.axt -o %s.axt.kaks -m LWL >> /dev/null 2>&1" % (query.file, query.file)
    os.system(commond)
    try:
        with open('%s.axt.kaks' % query.file) as f:
            f.readline()
            ks = f.readline().rstrip().split()[3]
            if ks == 'NA':
                ks = 0
    except:
        ks = 'ERROR'
    finally:
        os.remove('%s.axt' % query.file)
        os.remove('%s.axt.kaks' % query.file)
    if not quiet:
        print("Ks: ", ks)
    return float(ks) if not ks=='ERROR' else ks


def one_blast(start):
    end = start + window_length
    query_partial = query.get_partial(start, end)
    query_background = query.get_partial(start-background_length, end+background_length)
    try:
        query_result = blastn(query_partial, mask)
        top_query = query_result.top1()
        top_background = blastn(query_background, mask).top1()
    except:
        os.remove(query_partial.file)
        os.remove(query_background.file)
        return (start, ['NA'])

    # if top 1 only 1 and is background target
    # if top 1 not only 1 but with background target  
    set_query = set(x.title for x in top_query)
    set_background = set(x.title for x in top_background)
    insert = set_query.intersection(set_background)
    if insert:
        best_source = list(insert)
    # if top 1 only 1 but not background target
    # if top 1 not only 1 and without background target
    else:
        query_background_partial = query_result.extract(top_background[0].acc)
        if query_background_partial:
            ks2 = kaks(query_partial, query_background_partial[0])
            if ks2 == 'ERROR':
                ks2 = ks_threshold
        else:
            ks2 = ks_threshold

        candidate_source = []
        for hsp in top_query:
            b_ident = background_identity(query, hsp)
            ks1 = kaks(query_partial, hsp)
            if ks1 != 'ERROR':
                if b_ident > pidebt_threshold and ks1 < ks2:
                    candidate_source.append(hsp)
        if not candidate_source:
            candidate_source.extend(top_background)
        best_source = list(set([x.title for x in candidate_source]))
    
    os.remove(query_partial.file)
    os.remove(query_background.file)

    print('Background source for %s - %s:' % (start, end))
    for x in top_background:
        print("\t%s\t%s\t%s" % (x.acc, x.pident, x.title))
    print("Candidate source for %s - %s:" % (start, end))
    for x in top_query:
        print("\t%s\t%s\t%s" % (x.acc, x.pident, x.title))
    print('Best source: %s' % ', '.join(best_source))
    # print("\n")

    return (start, best_source)


def main(stop):
    # 创建一个工作者线程池
    pool = Pool(10) 
    # 在各个线程中打开url，并返回结果
    results = pool.map(one_blast, range(1, stop, step))
    #close the pool and wait for the work to finish
    pool.close() 
    pool.join()
    results = {x: y for x, y in results}
    return results


if __name__ == '__main__':
    os.chdir('/home/zeng/Desktop/recombination-analysis-0312')
    mask = 'idlist.2019-ncov-ratg-batlike-pangolin'
    query = Query('bat_SL_ZC45_cds.fasta')
    
    tasks = [50, 101, 152, 200, 251, 500]
    for window_length in tasks:
        stop = len(query) - window_length

        # one_blast(6603)

        s = time.time()
        results = main(stop)
        e = time.time()
        print("Total cost %s seconds.\n" % (e-s))

        with open('bat_SL_ZC45_80_%s_cds.json' % window_length, 'w') as f:
            f.write(json.dumps(results, indent='\t', sort_keys=True, separators=(',', ': ')))
