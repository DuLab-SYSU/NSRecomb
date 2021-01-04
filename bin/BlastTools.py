from io import StringIO
import tempfile
import subprocess

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Blast import NCBIXML

# blast_result.query
# alignment.hit_id, len(alignment.hsps), alignment.length
# hsp.score, hsp.expect, hsp.identities, hsp.align_length,
# hsp.identities/hsp.align_length,
# hsp.query_start, hsp.query_end, hsp.sbjct_start, hsp.sbjct_end


def get_hits(query_i, db, negative_seqidlist, qcov_hsp_perc=70):
    '''
    For each window, use blastn algorithm to search the hit with highest score form specified database.
    Parameters: 
        query_i            : SeqRecord , a specified window of query sequence
        db                 : string    , relative path for blast database
        negative_seqidlist : string    , relative path for negative sequence id list in binary format
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
        env={'PATH': '../package/ncbi-blast-2.11.0+/bin'}
    )

    blast_result = NCBIXML.read(StringIO(blastn_client.stdout))

    best_hits = []
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
    #TODO deal with no hits
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
        env={'PATH': '../package/ncbi-blast-2.11.0+/bin'}
    )
    # print(blastn_client.stdout)

    blastdbcmd_client = subprocess.run(
        'blastdbcmd -db %s -entry "%s"' % (db, blastn_client.stdout),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='unicode_escape',
        env={'PATH': '../package/ncbi-blast-2.11.0+/bin'}
    )
    # print(blastdbcmd_client.stdout)
    backbone = SeqIO.read(StringIO(blastdbcmd_client.stdout), 'fasta')
    #TODO deal with no backbone
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
        env={'PATH': '../package/ncbi-blast-2.11.0+/bin'}
    )

    blast_result = NCBIXML.read(StringIO(blastn_client.stdout))

    # print('backbone_i', len(blast_result.alignments), len(blast_result.alignments[0].hsps[0].sbjct))
    if len(blast_result.alignments) != 0:
        hit = blast_result.alignments[0]
        hsp = hit.hsps[0]
        backbone_i = SeqRecord(Seq(hsp.sbjct), id=hit.hit_id, name=hit.hit_id, description=hit.title)
        backbone_i.annotations['start'] = hsp.sbjct_start
        backbone_i.annotations['end'] = hsp.sbjct_end
    else:
        #TODO deal with no backbone
        backbone_i = SeqRecord(Seq('-'*len(query_i)), id='NA', name='NA', description='NA')
        # print(backbone_i)
    return backbone_i


if __name__ == '__main__':
    query = SeqIO.read('../data/SARS-CoV-2.ORF1ab.fasta', 'fasta')
    query_i = query[0: 501]

    backbone = get_gene_backbone(query, db='../db/all-cov', negative_seqidlist='../data/SARS-CoV-2.acc')
    backbone_i = get_window_backbone(query_i, backbone)
    hits = get_hits(query_i, db='../db/all-cov', negative_seqidlist='../data/SARS-CoV-2.acc', qcov_hsp_perc=70)
    print('Backbone: ')
    print(backbone, '\n')
    print("Backbone_i: ")
    print(backbone_i, '\n')
    print('# of best hits: ', len(hits))
    for hit in hits:
        print(hit)
