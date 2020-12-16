import os
from Bio.Blast.Applications import NcbimakeblastdbCommandline
import argparse


def build_db(db_fasta, out_db, input_type='fasta'):
    """
    Build local Blast v5 database
    """
    client = NcbimakeblastdbCommandline(
        dbtype='nucl',
        input_file=db_fasta,
        input_type=input_type,
        parse_seqids=True, out=out_db
    )
    client()


def convert_to_binary(seqid_file_in, seqid_file_out):
    """
    Convert a accession list to binary format
    """
    command = "blastdb_aliastool -seqid_file_in %s -seqid_file_out %s" \
              % (seqid_file_in, seqid_file_out)
    os.system(command)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Trans negative sequence id file from text mode to binary mode.')
    parser.add_argument('--seqid_file', '-s', help='File path for text mode seq id in file', required=True)
    parser.add_argument('--out_file', '-o', help='File path for binary mode seq id out file', required=True)
    parser.add_argument('--database', '-db', help='File path for database', required=False)
    parser.add_argument('--db_name', '-dbn', help='Alias for database', required=False)
    args = parser.parse_args()

    seqid_in_file = args.seqid_file
    seqid_out_file = args.out_file
    db_file = args.database
    db_name = args.db_name

    convert_to_binary(seqid_in_file, seqid_out_file)
    if not os.path.exists("db/"):
        os.makedirs("db/")
    if db_file and db_name:
        build_db(db_file, db_name)
