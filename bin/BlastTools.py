from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from io import StringIO
import tempfile
import subprocess

# blast_result.query
# alignment.hit_id, len(alignment.hsps), alignment.length
# hsp.score, hsp.expect, hsp.identities, hsp.align_length,
# hsp.identities/hsp.align_length,
# hsp.query_start, hsp.query_end, hsp.sbjct_start, hsp.sbjct_end


def get_hits(query_i,
             db='../db/test',
             negative_seqidlist='../sample/q.acc',
             outfmt=5,
             task='blastn',
             qcov_hsp_perc=70):
    blastn_client = NcbiblastnCommandline(
        db=db,
        negative_seqidlist=negative_seqidlist,
        outfmt=outfmt,
        task=task,
        qcov_hsp_perc=qcov_hsp_perc
    )
    stdout, _ = blastn_client(stdin=query_i.format('fasta'))
    blast_result = NCBIXML.read(StringIO(stdout))

    best_hits = []
    max_score = 0
    for alignment in blast_result.alignments:
        for hsp in alignment.hsps:
            if hsp.score >= max_score:
                hit = SeqRecord(Seq(hsp.sbjct), id=alignment.hit_id)
                hit.annotations['start'] = hsp.sbjct_start
                hit.annotations['end'] = hsp.sbjct_end
                hit.annotations['score'] = hsp.score
                hit.annotations['identity'] = hsp.identities/hsp.align_length
                best_hits.append(hit)
                max_score = hsp.score
    return best_hits


def _get_backbone(query, db='../db/test', negative_seqidlist='../sample/q.acc'):
    blastn_client = NcbiblastnCommandline(
        db=db,
        outfmt=5,
        negative_seqidlist=negative_seqidlist,
        max_target_seqs=1,
        max_hsps=1
    )
    stdout, _ = blastn_client(stdin=query.format('fasta'))
    blast_result = NCBIXML.read(StringIO(stdout))

    backbone_id = blast_result.alignments[0].hit_id

    blastdbcmd_client = subprocess.Popen(
        'blastdbcmd -db %s -entry "%s"' % (db, backbone_id),
        shell=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    stdout, _ = blastdbcmd_client.communicate()
    backbone = SeqIO.read(StringIO(stdout.decode('utf-8')), 'fasta')
    return backbone


def get_backbone(query_i, backbone, qcov_hsp_perc=70):
    named_backbone = tempfile.NamedTemporaryFile(mode='w+')
    SeqIO.write(backbone, named_backbone, 'fasta')
    named_backbone.flush()
    named_backbone.seek(0)
    blastn_client = NcbiblastnCommandline(
        subject=named_backbone.name,
        outfmt=5,
        max_hsps=1,
        qcov_hsp_perc=qcov_hsp_perc,
        task='blastn'
    )
    stdout, _ = blastn_client(stdin=query_i.format('fasta'))
    blast_result = NCBIXML.read(StringIO(stdout))

    # print('backbone_i', len(blast_result.alignments), len(blast_result.alignments[0].hsps[0].sbjct))
    if len(blast_result.alignments) != 0:
        hit = blast_result.alignments[0]
        hsp = hit.hsps[0]
        backbone_i = SeqRecord(Seq(hsp.sbjct), id=hit.hit_id)
        backbone_i.annotations['start'] = hsp.sbjct_start
        backbone_i.annotations['end'] = hsp.sbjct_end
    else:
        backbone_i = SeqRecord(Seq('-'*len(query_i)), id='NA')
        # print(backbone_i)
    return backbone_i


if __name__ == '__main__':
    query = SeqIO.read('../sample/query', 'fasta')
    query_i = query[0: 501]

    backbone = _get_backbone(query, db='../db/test', negative_seqidlist='../sample/q.acc')
    backbone_i = get_backbone(query_i, backbone)
    hits = get_hits(query_i,
                    db='../db/test',
                    negative_seqidlist='../sample/q.acc',
                    outfmt=5,
                    task='blastn',
                    qcov_hsp_perc=70)
    print('Backbone: ')
    print(backbone, '\n')
    print("Backbone_i: ")
    print(backbone_i, '\n')
    print('# of best hits: ', len(hits))
    for hit in hits:
        print(hit)
