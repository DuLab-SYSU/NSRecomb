import subprocess
from tempfile import NamedTemporaryFile


def string_to_file(string):
    file_like_obj = NamedTemporaryFile()
    file_like_obj.write(string.encode('utf-8'))
    file_like_obj.flush()
    file_like_obj.seek(0)
    return file_like_obj


def ms_client(nsam, rho, seq_length, f, lame):
    client = subprocess.Popen(
        'ms %s 1 -T -r %s %s -c %s %s | tail +4 | grep -v //'
        % (nsam, rho, seq_length, f, lame),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    stdout, stderr = client.communicate()
    print(stdout.decode('utf-8'))
    print(stderr.decode('utf-8'))
    return stdout.decode('utf-8')


def seq_gen_client(treefile, seq_length, scale, codon_rate):
    client = subprocess.Popen(
        'seq-gen -mHKY -l %s -s %s -p 50 -q %s -c %s -of'
        % (seq_length, scale, treefile, codon_rate),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    stdout, stderr = client.communicate()
    print(stdout.decode('utf-8'))
    print(stderr.decode('utf-8'))
    return stdout.decode('utf-8')


trees = ms_client(50, 3.0, 3000, 0.5, 500)
trees_handle = string_to_file(trees)
seqs = seq_gen_client(trees_handle.name, 3000, 3.0, '0.6 0.7 0.9')
with open('test_db', 'w') as w:
    w.write(seqs)
