import os
from Bio.Blast.Applications import NcbimakeblastdbCommandline


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
    seqid_file = "../sample/q.id"
    convert_to_binary(seqid_file, '../sample/q.acc')
    if not os.path.exists("db/"):
        os.makedirs("db/")
    build_db('../sample/test_db', 'db/test')
