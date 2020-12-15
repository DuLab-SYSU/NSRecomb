# requirements: blastn, mafft, KaKs_calculator
# input: query, blast_database, bnegative_seqidlist, backbone

import os, sys 
import numpy as np
import re, copy, time, json, random
from io import StringIO
from collections import Counter, defaultdict
from multiprocessing import Pool, Process

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Blast.Applications import NcbiblastnCommandline


BLASTDB = "all-cov/all-cov"

class Sequence:
    '''
    Sequence 类, 序列对象
    存储序列 acc title sequence 相对位置 等信息
    可以构建时输入 或初始化空对象,随后从fasta文件载入
    具有extract_ fasta_ 方法,分别可以提取部分序列 及格式化输出
    '''
    def __init__(self, acc=None, title=None, sequence=None, start=None, end=None):
        self.acc = acc
        self.title = title
        self.sequence = sequence
        self.start = start
        self.end = end
    
    def __len__(self):
        return len(self.sequence)

    def __setattr__(self, name, value):
        self.__dict__[name] = value

    def read_seqs(self, file):
        count = 0
        sequence = ''
        with open(file) as f:
            for line in f:
                if line.startswith(">"):
                    count += 1
                    if count == 2:
                        break
                    self.acc, *title = line[1:-1].split()
                    self.title = ' '.join(title)
                else:
                    sequence += line.strip()
            self.sequence = sequence
            self.file = file
            self.start = 1
            self.end = len(self.sequence)
    
    def extract_(self, start, end):
        sequence = self.sequence[start-1: end]
        return Sequence(self.acc, self.title, sequence, start, end)

    def fasta_(self):
        return ">%s %s_%s_%s\n%s\n" % (self.acc, self.start, self.end, self.title, self.sequence)


class Subject(Sequence):
    '''
    Subject 类
    Sequence 子类, 存储blastn结果,如 pident score
    '''
    def __init__(self, acc, title, sequence, start, end, pident=0, score=0):
        super().__init__(acc, title, sequence, start, end)
        self.pident = pident
        self.score = score
    
    def get_bident(self, query, length=1000, db=BLASTDB):
    #TODO 优化获取背景相似度
        start = 1 if self.start - length < 1 else self.start - length
        end = len(self) if self.end + length > len(self) else self.end + length
        commond ="blastdbcmd -db %s -entry %s -range %s-%s | \
                  blastn -subject %s -outfmt '6 pident' \
                  -max_hsps 1" % (db, self.acc, start, end, query.file)
        try:
            bident = float(os.popen(commond).read().strip())
        except:
            bident = 0
        return bident


class Subjects(list):
    """list of Subject"""
    def __init__(self, _in=list()):
        self.extend(_in)
    
    def _get_max_score(self):
        max_score = 0
        for subject in self:
            if subject.score > max_score:
                max_score = subject.score
        return max_score

    def best_subjects(self):
        max_score = self._get_max_score()
        return Subjects([subject for subject in self if subject.score == max_score])
    
    def get_item(self, item):
        for subject in self:
            if subject.acc == item.acc:
                return subject
        else:
            return None
    
    def __contains__(self, item):
        return item.acc in [i.acc for i in self]


class Query(Sequence):
    '''Sequence 子类 额外添加了一些方法'''
    def __init__(self, acc, title, sequence, start, end):
        super().__init__(acc, title, sequence, start, end)

    def blastn(self, db=BLASTDB, qcov_hsp_perc=70, perc_identity=70, negative_seqidlist='negative.id', task='blastn'):
        '''
        调用外部blastn, 返回结果, 可能为空
        '''
        self.subjects = Subjects()
        commond = "echo '%s' | blastn -db %s  -qcov_hsp_perc %s -max_hsps 1 -perc_identity %s \
                            -outfmt '6 qacc qstart qend sacc sstart send pident score sseq stitle' \
                            -negative_seqidlist %s -task %s 2> /dev/null" % (self.sequence, db, qcov_hsp_perc, perc_identity, negative_seqidlist, task)
        res = os.popen(commond).readlines()
        for x in res:
            _, _, _, sacc, sstart, send, pident, score, sseq, stitle = x.rstrip().split('\t')
            subject = Subject(sacc, stitle, sseq, int(sstart), int(send), float(pident), int(score))
            self.subjects.append(subject)

    def get_best_subjects(self): # list
        return self.subjects.best_subjects()
    
    def get_backbone(self, backbone_):
        commond = "echo '%s' | blastn -subject backbone.fasta\
                    -outfmt '6 sacc sstart send stitle sseq pident' -task blastn -max_hsps 1\
                    2> /dev/null" % (self.sequence)
        res = os.popen(commond).read().split()
        if res:
            sacc, sstart, send, _, sseq, _ = res
            backbone = Sequence(sacc, backbone_, sseq, int(sstart), int(send))
        else:
            backbone = None
        return backbone


class Triplet:
    def __init__(self, query, backbone, subject, align=True):
        self.query = query
        self.backbone = backbone
        self.subject = subject
        if align:
            self._align()
    
    def _to_fasta(self, _out):
        _query = self.query.fasta_()
        _backbone = self.backbone.fasta_()
        _subject = self.subject.fasta_()
        f = open(_out, 'w')
        f.write(''.join([_query, _backbone, _subject]))
        f.close()

    def _align(self):
        _tmp = '%s_%s_%s.fasta' % (self.query.acc, self.query.start, self.query.end)
        self._to_fasta(_tmp)
        ress = os.popen('mafft --auto --quiet %s' % (_tmp)).read().split('>')[1:]
        tasks = [self.query, self.backbone, self.subject]
        for i, res in enumerate(ress):
            tasks[i].sequence = ''.join(res.split('\n')[1:])
        os.remove(_tmp)
    
    def kaks(self): # if ks of subject > ks of backbone return True
        _tmp = '%s_%s_%s' % (self.query.acc, self.query.start, self.query.end)
        with open('%s.axt' % _tmp, 'w') as f:
            f.write('%s-%s\n' % (self.query.acc, self.backbone.acc))
            f.write('%s\n%s\n' % (self.query.sequence, self.backbone.sequence))
            f.write('\n')
            f.write('%s-%s\n' % (self.query.acc, self.subject.acc))
            f.write('%s\n%s\n' % (self.query.sequence, self.subject.sequence))
        commond = "KaKs_Calculator -i %s.axt -o %s.axt.kaks -m NG >> /dev/null 2>&1" % (_tmp, _tmp)
        os.system(commond)
        with open('%s.axt.kaks' % _tmp) as f:
            f.readline()
            qb_ks = f.readline().rstrip().split()[3]
            qs_ks = f.readline().rstrip().split()[3]
        # print('\tdS(backbone): %s, dS(subject): %s' % (f(b_ks), f(s_ks)))
        os.remove('%s.axt' % _tmp)
        os.remove('%s.axt.kaks' % _tmp)
        f = lambda x: float(x) if x != 'NA' else 0
        return (f(qb_ks), f(qs_ks))

    def _identity(self, x, y):
        count = 0
        for x_i, y_i in zip(x, y):
            if x_i == y_i:
                count += 1
        return count / len(x)

    def bootstrap_support(self, sample_size=1000):
        query_seq = np.array(list(self.query.sequence))
        backbone_seq = np.array(list(self.backbone.sequence))
        subject_seq = np.array(list(self.subject.sequence))
        _length = len(query_seq)
        np.random.seed(123)
        idx = np.random.randint(0, _length, size=(sample_size, _length))
        query_seqs, backbone_seqs, subject_seqs = query_seq[idx], backbone_seq[idx], subject_seq[idx]
        count = 0
        for query_seq_i, backbone_seq_i, subject_seq_i in zip(query_seqs, backbone_seqs, subject_seqs):
            qs_ident, qb_ident = self._identity(query_seq_i, subject_seq_i), self._identity(query_seq_i, backbone_seq_i)
            if  qs_ident > qb_ident:
                count += 1
        return count / sample_size


def one_window(start, *args, **kwargs):
    # status: I: backbone in subject; V: bootstrap and dS pass; E: empty
    end = start + window_size - 1
    query_i = query.extract_(start, end)
    query_i = Query(query_i.acc, query_i.title, query_i.sequence, start, end)
    query_i.blastn()
    subjects = query_i.get_best_subjects()

    backbone = query_i.get_backbone(backbone_title)

    winner = []
    if backbone and len(subjects)!=0:
        for subject in subjects:
            triplet = Triplet(query_i, backbone, subject, align=True)
            bootstrap_value = triplet.bootstrap_support()
            try:
                qb_ks, qs_ks = triplet.kaks()
            except:
                qb_ks, qs_ks = 1, 0
                print(str(start), end=' ,')
            if (bootstrap_value > bootstrap_threshod) and (qs_ks < qb_ks):
                winner.append((subject.title, bootstrap_value))
        if not winner:
            winner.append((backbone.title, bootstrap_value))
    elif backbone and len(subjects)==0:
        winner.append((backbone.title, 1))
    elif len(Subjects)!=0:
        for subject in subjects:
            winner.append((subject.title, 1))
    return start, winner


def main(start, stop, step):
    # 创建一个工作者线程池
    pool = Pool(10) 
    results = pool.map(one_window, range(start, stop, step))
    #close the pool and wait for the work to finish
    pool.close() 
    pool.join()
    results = {x: y for x, y in results}
    return results


if __name__ == '__main__':
    path = "/home/zeng/Desktop/recombination_task/task2/"
    tasks = ['SARS-CoV-2', 'RmYN02', 'RaTG13', 'PangolinGD', 'PangolinGX']
    cds = ['ORF1ab', 'S', 'M', 'N']
    backbones = ['RmYN02', 'RaTG13', 'PangolinGD', 'PangolinGX', 'BatSL']

    window_size = 501
    step = 3
    #b_length = 1000
    #b_ident_threshod = 80
    bootstrap_times = 500
    bootstrap_threshod = 0.8

    for i in range(1):
        print("\nTask: %s" % tasks[i], end=' ')
        os.chdir(path + tasks[i])
        backbone_title = backbones[i]
        for cds_i in cds:
            if cds_i:
                print("\n\tCDS: %s, " % cds_i, end=" ")
                query_file = tasks[i]+'_'+cds_i+'.fasta'
                query = Sequence()
                query.read_seqs(query_file)
                
                # one_window(1)

                # for i in range(1, 2500, 3):
                    # print(one_window(i))
                stop = len(query) - window_size + 2
                # results = {i: one_window(i) for i in range(21202, 21142, 3)}
                results = main(1, stop, step)

                with open(query_file.replace('.fasta', '.501.json'), 'w') as f:
                    f.write(json.dumps(results, indent='\t', sort_keys=True, separators=(',', ': ')))


