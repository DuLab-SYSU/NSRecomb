import os
import subprocess


def build_db(db_fasta, out_db, input_type='fasta'):
    """
    Build local Blast v5 database
    """
    subprocess.run(
        'makeblastdb -dbtype nucl -in %s -input_type %s -parse_seqids -out %s'
        % (db_fasta, input_type, out_db),
        shell=True,
        env={'PATH': BLAST_PATH}
    )


def convert_to_binary(seqid_file_in, seqid_file_out):
    """
    Convert a accession list to binary format
    """
    subprocess.run(
        "blastdb_aliastool -seqid_file_in %s -seqid_file_out %s"
        % (seqid_file_in, seqid_file_out),
        shell=True,
        env={'PATH': BLAST_PATH}
    )


if __name__ == '__main__':
    BLAST_PATH = '../package/ncbi-blast-2.11.0+/bin'
