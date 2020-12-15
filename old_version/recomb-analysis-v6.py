# requirements: blastn, mafft, KaKs_calculator
# input: query, blast_database, bnegative_seqidlist, backbone

import os, sys 
import numpy as np
import re, copy, time, json
from collections import Counter, defaultdict
from multiprocessing import Pool, Process
from multiprocessing.dummy import Pool as ThreadPool

BLASTDB = "all-cov/all-cov"

class Sequence:
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

    def load_seqs(self, file):
        sequence = ''
        with open(file) as f:
            for line in f:
                if line.startswith(">"):
                    self.acc, *title = line[1:-1].split()
                    self.title = ' '.join(title)
                else:
                    sequence += line.strip()
                self.sequence = sequence
            else:
                self.file = file
                self.start = 1
                self.end = len(self.sequence)
    
    def extract_part(self, start, end):
        sequence = self.sequence[start-1: end]
        return Sequence(self.acc, self.title, sequence, start, end)

    def to_fasta(self):
        return ">%s %s_%s_%s\n%s\n" % (self.acc, self.start, self.end, self.title, self.sequence)


class Subject(Sequence):
    def __init__(self, acc, title, sequence, start, end, pident=0, score=0):
        super().__init__(acc, title, sequence, start, end)
        self.pident = pident
        self.score = score
    
    def get_bident(self, query, length=1000, db=BLASTDB):
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
    def __init__(self, acc, title, sequence, start, end):
        super().__init__(acc, title, sequence, start, end)

    def blastn(self, db=BLASTDB):
        self.subjects = Subjects()
        commond = "echo '%s' | blastn -db %s  -qcov_hsp_perc 70 -max_hsps 1 -perc_identity 70 \
                            -outfmt '6 qacc qstart qend sacc sstart send pident score sseq stitle' \
                            -negative_seqidlist negative.id -task blastn 2> /dev/null" % (self.sequence, db)
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
        sacc, sstart, send, _, sseq, _ = res
        backbone = Sequence(sacc, backbone_, sseq, int(sstart), int(send))
        return backbone


class Triplet:
    def __init__(self, query, backbone, subject, align=True):
        self.query = query
        self.backbone = backbone
        self.subject = subject
        if align:
            self._align()
    
    def _to_fasta(self, _out):
        _query = self.query.to_fasta()
        _backbone = self.backbone.to_fasta()
        _subject = self.subject.to_fasta()
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
        #try:
        #    cir = f(qs_ks) < f(qb_ks)
        #except:
        #    cir = True
        #return cir
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
            qs_ident = self._identity(query_seq_i, subject_seq_i)
            qb_ident = self._identity(query_seq_i, backbone_seq_i)
            if  qs_ident> qb_ident:
                count += 1
        return count / sample_size


def one_window(start, *args, **kwargs):
    # status: I: backbone in subject; V: bootstrap and dS pass; E: empty
    end = start + window_size - 1
    query_i = query.extract_part(start, end)
    query_i = Query(query_i.acc, query_i.title, query_i.sequence, start, end)
    query_i.blastn()
    subjects = query_i.get_best_subjects()
    #subjects = Subjects([subject for subject in subjects if subject.get_bident(query) >= b_ident_threshod])

    try:
        backbone = query_i.get_backbone(backbone_title)
    except:
        backbone = False

    if backbone:
        if backbone in subjects:
            winner = [('I', backbone.title, 1)]
            #print('Postion: %s-%s\n\tWinner: %s, Status: I' % (start, end, backbone.title))
        else:
            winner = []
            bootstrap_value_list = []
            #print('Position: %s-%s' % (start, end))
            for subject in subjects:
                triplet = Triplet(query_i, backbone, subject, align=True)
                bootstrap_value = triplet.bootstrap_support()
                bootstrap_value_list.append(bootstrap_value)
                try:
                    qb_ks, qs_ks = triplet.kaks()
                except:
                    qb_ks, qs_ks = 1, 0
                    print(str(start), end=' ,')                
                #print('\tSubject: %s, Pident: %s, Bootstrap value: %s, dS(QB): %s, dS(QS): %s' % (subject.title, subject.pident, bootstrap_value, qb_ks, qs_ks))
                if (bootstrap_value > bootstrap_threshod) and (qs_ks < qb_ks):
                    winner.append(('V', subject.title, bootstrap_value))
                    #print("\tWinner: %s, Status: V" % subject.title)
            if not winner:
                try:
                    bootstrap_value = sum(bootstrap_value_list)/len(bootstrap_value_list)
                except:
                    print("no subjects", start)
                    bootstrap_value = 0
                winner.append(('E', backbone.title, bootstrap_value))
                #print('\tWinner: %s, Status: E' % (backbone.title))
    else:
        winner = []
        for subject in subjects:
            winner.append(('X', subject.title, 0))
            #print("Position: %s-%s\n\tWinner: %s, Status: X" % (start, end, subject.title))
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
    path = "/home/zeng/Desktop/recombination_task/"
    tasks = ['SARS-CoV-2', 'RmYN02', 'RaTG13', 'PangolinGD', 'PangolinGX']
    cds = ['ORF1ab', 'S', 'M', 'N']
    backbones = ['RmYN02', 'RaTG13', 'PangolinGD', 'PangolinGX', 'BatSL']

    window_size = 501
    step = 3
    #b_length = 1000
    #b_ident_threshod = 80
    bootstrap_times = 500
    bootstrap_threshod = 0.8

    for i in range(5):
        print("\nTask: %s" % tasks[i], end=' ')
        os.chdir(path + tasks[i])
        backbone_title = backbones[i]
        for cds_i in cds:
            if cds_i:
                print("\n\tCDS: %s, " % cds_i, end=" ")
                query_file = tasks[i]+'_'+cds_i+'.fasta'
                query = Sequence()
                query.load_seqs(query_file)

                #one_window(1312)
                stop = len(query) - window_size + 2
                #results = {i: one_window(i) for i in range(21202, 21142, 3)}
                results = main(1, stop, step)

                with open(query_file.replace('.fasta', '.501.json'), 'w') as f:
                    f.write(json.dumps(results, indent='\t', sort_keys=True, separators=(',', ': ')))


