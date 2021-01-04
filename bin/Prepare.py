import os
import subprocess
# from Bio.Blast.Applications import NcbimakeblastdbCommandline
import argparse


def build_db(db_fasta, out_db, input_type='fasta'):
    """
    Build local Blast v5 database
    """
    subprocess.run(
        'makeblastdb -dbtype nucl -in %s -input_type %s -parse_seqids -out %s'
        % (db_fasta, input_type, out_db),
        shell=True,
        env={'PATH': '../package/ncbi-blast-2.11.0+/bin'}
    )


def convert_to_binary(seqid_file_in, seqid_file_out):
    """
    Convert a accession list to binary format
    """
    subprocess.run(
        "blastdb_aliastool -seqid_file_in %s -seqid_file_out %s"
        % (seqid_file_in, seqid_file_out),
        shell=True,
        env={'PATH': '../package/ncbi-blast-2.11.0+/bin'}
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Required preparation work for recombintaion detection.')
    parser.add_argument('--seqid_file', '-s', help='File path for text mode seq id in file', required=False, default=None)
    parser.add_argument('--seqout_file', '-o', help='File path for binary mode seq id out file', required=False, default=None)
    parser.add_argument('--db_file', '-db', help='File path for database', required=False, default=None)
    parser.add_argument('--db_name', '-dbn', help='Alias for database', required=False, default=None)
    parser.add_argument('--task', '-t', help='Kind of task:\n\
                        1: build local blast database, required to specify --db_file --db_name\n\
                        2: convrt negative sequence id text file to binary file', required=True, type=int)
    args = parser.parse_args()

    seqid_in_file = args.seqid_file
    seqid_out_file = args.seqout_file
    db_file = args.db_file
    db_name = args.db_name
    task = args.task

    if task == 2:
        convert_to_binary(seqid_in_file, '../data/' + seqid_out_file)
    else:
        if not os.path.exists("../db/"):
            os.makedirs("../db/")
        build_db(db_file, '../db/' + db_name)
