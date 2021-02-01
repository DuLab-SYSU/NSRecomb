from io import StringIO
import tempfile
import subprocess

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Blast import NCBIXML


def get_hits(query_i, db, negative_seqidlist, qcov_hsp_perc=70):
    '''
    For each window, blastn algorithm is used to search the hit with highest score from specified database.
    Parameters: 
        query_i            : SeqRecord , a specified window of query sequence
        db                 : string    , path for blast database
        negative_seqidlist : string    , path for negative sequence id list in binary format
        qcov_hsp_perc      : int       , least coverage for each hsp of a hit
    Return: list, a list of hits with highest score (SeqRecord)
    '''
    blastn_client = subprocess.run(
        'blastn -db %s -outfmt 5 -max_hsps 1 -qcov_hsp_perc %s -task blastn -negative_seqidlist %s'
        % (db, qcov_hsp_perc, negative_seqidlist),
        shell=True,
        input=repr(query_i.format('fasta')),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='unicode_escape',
        env={'PATH': BLAST_PATH}
    )

    best_hits = []
    try:
        blast_result = NCBIXML.read(StringIO(blastn_client.stdout))
        max_score = 0
        for alignment in blast_result.alignments:
            for hsp in alignment.hsps:
                if hsp.score >= max_score:
                    hit = SeqRecord(Seq(hsp.sbjct), id=alignment.hit_id, name=alignment.hit_id, description=alignment.title)
                    hit.annotations['start'] = hsp.sbjct_start
                    hit.annotations['end'] = hsp.sbjct_end
                    hit.annotations['score'] = hsp.score
                    hit.annotations['identity'] = hsp.identities / hsp.align_length
                    best_hits.append(hit)
                    max_score = hsp.score
    except:
        #! handle no hits
        best_hits.append(SeqRecord(Seq('-'*len(query_i)), id='NA', name='NA', description='NA'))
    return best_hits


def get_gene_backbone(query, db, negative_seqidlist):
    '''
    Get the backbone sequence for a given gene sequence.
    Parameters:
        query              : SeqRecord , full gene sequence
        db                 : string    , relative path for blast database
        negative_seqidlist : string    , relative path for negative sequence id list in binary format
    Return: SeqRecord, backbone sequence record
    '''
    blastn_client = subprocess.run(
        'blastn -db %s -outfmt "6 sseqid" -negative_seqidlist %s -max_target_seqs 1 -max_hsps 1'
        % (db, negative_seqidlist),
        shell=True,
        input=repr(query.format('fasta')),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='unicode_escape',
        env={'PATH': BLAST_PATH}
    )
    #! deal with no backbone
    backbone_id = blastn_client.stdout.strip('\n')
    if backbone_id:
        blastdbcmd_client = subprocess.run(
            'blastdbcmd -db %s -entry "%s"' % (db, backbone_id),
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            encoding='unicode_escape',
            env={'PATH': BLAST_PATH}
        )
        backbone = SeqIO.read(StringIO(blastdbcmd_client.stdout), 'fasta')
    else:
        raise ValueError("""
                         No backbone sequences in the database, which may due to the
                         long divergent time between your query sequence and the target
                         database.\n
                         Please change your database or specify the backbone
                         by your own with switch -b.
                         """)
    return backbone


def get_window_backbone(query_i, backbone, qcov_hsp_perc=70):
    '''
    Get the coresponding backbone sequence for a specified window query
    Parameters:
        query_i       : SeqRecord , window query
        backbone      : SeqRecord , specified backbone
        qcov_hsp_perc : int       , least coverage for each hsp of a hit
    '''
    named_backbone = tempfile.NamedTemporaryFile(mode='w+')
    SeqIO.write(backbone, named_backbone, 'fasta')
    named_backbone.flush()
    named_backbone.seek(0)
    
    blastn_client = subprocess.run(
        'blastn -subject %s -outfmt 5 -max_hsps 1 -qcov_hsp_perc %s -task blastn'
        % (named_backbone.name, qcov_hsp_perc),
        shell=True,
        input=repr(query_i.format('fasta')),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='unicode_escape',
        env={'PATH': BLAST_PATH}
    )

    try:
        blast_result = NCBIXML.read(StringIO(blastn_client.stdout))
        hit = blast_result.alignments[0]
        hsp = hit.hsps[0]
        backbone_i = SeqRecord(Seq(hsp.sbjct), id=hit.hit_id, name=hit.hit_id, description=hit.title)
        backbone_i.annotations['start'] = hsp.sbjct_start
        backbone_i.annotations['end'] = hsp.sbjct_end
    except:
        #! handle no backbone
        backbone_i = SeqRecord(Seq('-'*len(query_i)), id='NA', name='NA', description='NA')
    return backbone_i


if __name__ == '__main__':
    BLAST_PATH = '../package/ncbi-blast-2.11.0+/bin'

    query = SeqIO.read('../example/SARS-CoV-2.ORF1ab.fasta', 'fasta')
    query_i = query[0: 501]

    backbone = get_gene_backbone(query, db='../db/all_cov_unify', negative_seqidlist='../example/SARS-CoV-2.acc')
    backbone_i = get_window_backbone(query_i, backbone)
    hits = get_hits(query_i, db='../db/all-cov', negative_seqidlist='../example/SARS-CoV-2.acc', qcov_hsp_perc=70)
    print('Backbone: ')
    print(backbone, '\n')
    print("Backbone_i: ")
    print(backbone_i, '\n')
    print('# of best hits: ', len(hits))
    for hit in hits:
        print(hit)
