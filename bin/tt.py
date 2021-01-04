import os
import subprocess
from Bio import SeqIO
from io import StringIO


ENV = os.environ.copy()
# ENV['PATH']

# p = subprocess.run("KaKs_Calculator -h", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='unicode_escape', env={'PATH': '../package/KaKs_Calculator2.0/bin'})
# print(p.returncode)
# print(p.stdout)
# print(p.stderr)

# p2 = subprocess.run("mafft.bat -help", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='unicode_escape', cwd='../package/mafft-linux64')
# print(p2.returncode)
# print(p2.stdout)
# print(p2.stderr)

q = SeqIO.read('../data/SARS-CoV-2.ORF1ab.fasta', 'fasta')
# print(q.format('fasta'))

p4 = subprocess.run(
    'echo %s' % repr(q.format('fasta')),
    shell=True,
    encoding='unicode_escape',
    stdout=subprocess.PIPE
)

# print(p4.stdout)


p3 = subprocess.Popen(
    'blastn -db %s -outfmt "6 sseqid" -negative_seqidlist %s -max_target_seqs 1'
    % ('../db/all-cov', '../data/SARS-CoV-2.acc'),
    shell=True,
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    encoding='unicode_escape',
    env={'PATH': '../package/ncbi-blast-2.11.0+/bin'})

# print(p3.returncode)
# print(p3.stdout)
# print(p3.stderr)

stdout, stderr = p3.communicate(repr(q.format('fasta')))
print(stdout)
# print(stderr)
