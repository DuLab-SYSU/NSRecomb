import os, sys 
import numpy as np
import re, copy, time, json
from collections import Counter, defaultdict
from multiprocessing import Pool, Process
from multiprocessing.dummy import Pool as ThreadPool


class Sequence:
    def __init__(self, acc=None, title=None, sequence=None, start=None, end=None):
        self.acc = acc
        self.title = title
        self.sequence = sequence
        self.start = start
        self.end = end
    
    def __len__(self):
        return len(self.sequence)
    
    def load_seqs(self, file):
        self.sequence = ''
        with open(file) as f:
            for line in f:
                if line.startswith(">"):
                    self.acc, *title = line[1:-1].split()
                    self.title = ' '.join(title)
                else:
                    self.sequence += line.strip()
            else:
                self.file = file
                self.start = 1
                self.end = len(self.sequence)
    
    def extract_part(self, start, end):
        start = start if start > 0 else 1
        sequence = self.sequence[start-1: end]
        return Sequence(self.acc, self.title, sequence, start, end)
    
    def update(self, seq):
        self.sequence = seq

    def fasta(self):
        return ">%s %s..%s..%s\n%s\n" % (self.acc, self.start, self.end, self.title, self.sequence)


class Subject(Sequence):
    def __init__(self, acc, title, sequence, start, end, pident=0, score=0):
        super().__init__(acc, title, sequence, start, end)
        self.pident = pident
        self.score = score
    
    def get_bident(self, to_compare, length=1000, db='all-cov/all-cov'):
        start = self.start-length if self.start-length > 0 else 1
        end = self.end+length
        commond ="blastdbcmd -db %s -entry %s -range %s-%s | \
                  blastn -subject %s -outfmt '6 pident' \
                  -max_hsps 1" % (db, self.acc, start, end, to_compare.file)
        try:
            bident = float(os.popen(commond).read().strip())
        except:
            bident = 0
        return bident

    def __eq__(self, item):
        return self.acc == item.acc


class Result(list):
    """list of Subject class"""
    def __init__(self, _in=list()):
        for i in _in:
            self.append(i)
        
    def top1(self):
        results = sorted(self, key=lambda x: x.score, reverse=True)
        return Result([subject for subject in results if subject.score == results[0].score])
    
    def get_item(self, item):
        try:
            return [value for value in self if value.acc == item.acc][0]
        except:
            raise KeyError('No such item!')
    
    def __contains__(self, item):
        return item.acc in [i.acc for i in self]

    def intersection(self, item):
        return Result([i for i in self if i in item])


class Query(Sequence):
    def __init__(self, acc, title, sequence, start, end):
        super().__init__(acc, title, sequence, start, end)
        self.results = Result()

    def blastn(self, mask, db='all-cov/all-cov'):
        commond = "echo '%s' | blastn -db %s -max_hsps 1 \
                        -outfmt '6 qacc qstart qend sacc sstart send pident score sseq stitle' \
                        -negative_seqidlist %s -task blastn \
                        2> /dev/null" % (self.sequence, db, mask)
        for x in os.popen(commond).readlines():
            _, _, _, sacc, sstart, send, pident, score, sseq, stitle = x.rstrip().split('\t')        
            subject = Subject(sacc, stitle, sseq, int(sstart), int(send), float(pident), int(score))
            self.results.append(subject)

    def get_top1(self):
        return self.results.top1()
    
    def get_item(self, item):
        return self.results.get_item(item)


class Triplet:
    def __init__(self, query, background, subject):
        self.query = query
        self.background = background
        self.subject = subject
    
    def _fasta(self, _out):
        query_ = self.query.fasta()
        background_ = self.background.fasta()
        subject_ = self.subject.fasta()
        f = open(_out, 'w')
        f.write(''.join([query_, background_, subject_]))
        f.close()
        # return query_, background_, subject_

    def align(self):
        _tmp = '%s_%s_%s.fasta' % (self.query.acc, self.query.start, self.query.end)
        self._fasta(_tmp)
        ress = os.popen('mafft --auto --quiet %s' % (_tmp)).read().split('>')[1:]
        tasks = [self.query, self.background, self.subject]
        for i, res in enumerate(ress):
            tasks[i].update(''.join(res.split('\n')[1:]))
        os.remove(_tmp)
    
    def kaks(self): # if ks of subject > ks of background return True
        _tmp = '%s_%s_%s' % (self.query.acc, self.query.start, self.query.end)
        # print(len(self.query.sequence), len(self.background.sequence), len(self.subject.sequence))
        # print(self.query.sequence, self.background.sequence, self.subject.sequence, sep='\n')
        with open('%s.axt' % _tmp, 'w') as f:
            f.write('%s-%s\n' % (self.query.title, self.background.title))
            f.write('%s\n%s\n' % (self.query.sequence, self.background.sequence))
            f.write('\n')
            f.write('%s-%s\n' % (self.query.title, self.subject.title))
            f.write('%s\n%s\n' % (self.query.sequence, self.subject.sequence))
        commond = "KaKs_Calculator -i %s.axt -o %s.axt.kaks -m LWL >> /dev/null 2>&1" % (_tmp, _tmp)
        os.system(commond)
        with open('%s.axt.kaks' % _tmp) as f:
            f.readline()
            b_ks = f.readline().rstrip().split()[3]
            s_ks = f.readline().rstrip().split()[3]
            f = lambda x: float(x) if x != 'NA' else 0
        print('\tdS(back): %s, dS(subject): %s' % (f(b_ks), f(s_ks)))
        os.remove('%s.axt' % _tmp)
        os.remove('%s.axt.kaks' % _tmp)
        return f(b_ks) > f(s_ks)


    def _identity(self, x, y):
        count = 0
        for x_i, y_i in zip(x.sequence, y.sequence):
            if x_i == y_i:
                count += 1
        return count / len(x)

    def is_subject(self):
        background_pident = self._identity(self.query, self.background)
        subject_pident = self._identity(self.query, self.subject)
        # print('back: ', background_pident, 'subject: ', subject_pident)
        return background_pident < subject_pident


class Bootstrap:
    def __init__(self, triplet):
        self.triplet = triplet
        self._query = np.array(list(triplet.query.sequence))
        self._background = np.array(list(triplet.background.sequence))
        self._subject = np.array(list(triplet.subject.sequence))
        self._length = len(triplet.query.sequence)

    def sampling(self, sample_size=1000):
        samples = []
        idx = np.random.randint(0, self._length, size=(sample_size, self._length))
        querys = self._query[idx]
        backgrounds = self._background[idx]
        subjects = self._subject[idx]
        for query, background, subject in zip(querys, backgrounds, subjects):
            tmp = copy.deepcopy(self.triplet)
            tmp.query.update(''.join(query))
            tmp.background.update(''.join(background))
            tmp.subject.update(''.join(subject))
            samples.append(tmp)
        return samples


def get_background(sequence, mask, start, end):
    sequence_background = sequence.extract_part(start, end)
    query_background = Query(sequence_background.acc, sequence_background.title, sequence_background.sequence, sequence_background.start, sequence_background.end)
    query_background.blastn(mask)
    back = query_background.results.top1()
    return back[0]


def one_window(start, *args, **kwargs):
    # status: I: subject coincidence with background; V: bootstrap verified; E: bootstarp reject
    end = start + window_size - 1
    b_start = start - b_length
    b_end = end + b_length

    sequence_part = sequence.extract_part(start, end)
    query = Query(sequence_part.acc, sequence_part.title, sequence_part.sequence, sequence_part.start, sequence_part.end)
    query.blastn(mask)
    background = query.get_item(get_background(sequence, mask, b_start, b_end))
    subjects = Result([subject for subject in query.get_top1() if subject.get_bident(sequence) > b_ident_threshod])
    if background in subjects or len(subjects)==0:
        winner = [('I', background.title, 1)]
        print('Postion: %s-%s, Background: %s, Rejected due to in subjects, winner: %s' % (start, end, background.title,background.title))
    else:
        winner = []
        supports = 0
        print('Postion: %s-%s, Background: %s, Pident: %s' % (start, end, background.title, background.pident))
        for subject in subjects:
            triplet = Triplet(query, background, subject)
            triplet.align()
            if True:
            # if triplet.kaks(): # ks(subject) < ks(background)
                bootstrap = Bootstrap(triplet)
                replicates = bootstrap.sampling(bootstrap_times)
                passed = [1 for replicate in replicates if replicate.is_subject()]
                bootstrap_support = sum(passed) / bootstrap_times
                supports += bootstrap_support
                print("\tSubject: %s, Pident: %s, Bootstrap support: %s" % (subject.title, subject.pident, bootstrap_support))
                if bootstrap_support >= bootstrap_threshod:
                    winner.append(('V', subject.title, bootstrap_support))
                    print("\tWinner: %s" % (subject.title))
            # else:
            #     print('\tSubject: %s, Reject due to dS(subject) > dS(background)' % subject.title)
        if not winner:
            winner.append(('E', background.title, supports/len(subjects)))
            print('\tWinner: %s' % (background.title))
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
    os.chdir("/home/zeng/Desktop/recombination-analysis-0328")
    # 'bat_SL_ZC45_cds.fasta', 'bat_SL_ZXC21_cds.fasta', 'idlist.2019-ncov-ratg-batlike-pangolin', 'idlist.2019-ncov-ratg-batlike-pangolin'
    files = ['ratg13-cds.fasta', 'pangolin-gd-cds.fasta', 'pangolin-gx-cds.fasta']
    masks = ['idlist.2019-ncov-ratg', 'idlist.2019-ncov-ratg-pangolin-gd', 'idlist.2019-ncov-ratg-pangolin']

    # files = ["ratg13-cds.fasta"]
    # masks = ['idlist.2019-ncov-ratg']

    window_size = 500
    step = 3
    b_length = 1000
    b_ident_threshod = 80
    bootstrap_times = 500
    bootstrap_threshod = 0.8

    for mask, file in zip(masks, files):
        sequence = Sequence()
        sequence.load_seqs(file)
        stop = len(sequence.sequence) - window_size
        # one_window(8399)
        # results = {i: one_window(i)[1] for i in range(8740, 8802, 3)}
        results = main(1, stop, step)
        with open(file.replace('.fasta', '.%s.json' % bootstrap_threshod), 'w') as f:
            f.write(json.dumps(results, indent='\t', sort_keys=True, separators=(',', ': ')))


