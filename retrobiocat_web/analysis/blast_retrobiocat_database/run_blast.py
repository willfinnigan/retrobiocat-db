import os
import time
from decimal import Decimal
from io import StringIO

import numpy as np
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline

from retrobiocat_web.analysis.blast_retrobiocat_database.make_blast_db import make_blast_db
from retrobiocat_web.logging import add_logger
from retrobiocat_web.analysis.blast_retrobiocat_database.blast_dirs import BLAST_DB, BLAST_DIR

LOGGER = add_logger('Run_BioCatDB_Blast_DB', level='DEBUG')

def run_blast(query_seq):
    protein_name = 'query_seq'

    fasta_path = f"{BLAST_DIR}/{protein_name}.fasta"
    with open(fasta_path, 'w') as file:
        file.write(f'>{protein_name}\n')
        file.write(f"{query_seq}\n")

    output = NcbiblastpCommandline(query=fasta_path, db=BLAST_DB, evalue=0.5,
                                   num_threads=1, outfmt=5, num_alignments=1000)()[0]
    blast_record = NCBIXML.read(StringIO(output))

    os.remove(fasta_path)
    LOGGER.debug(f"sequence blasted")

    return blast_record

def process_blast_record(query_seq, blast_record):
    """ Process the blast record generating a list of alignments"""

    query_length = len(query_seq)

    alignment_names, alignment_scores, identities, coverages = [], [], [], []
    for alignment in blast_record.alignments:
        if len(alignment.hsps) == 1:
            hsp = alignment.hsps[0]
            subject_name = alignment.title.replace(f"{alignment.hit_id} ", "")
            coverage = round((hsp.align_length / query_length),3)
            identity = round((hsp.identities / hsp.align_length),3)

            score = _calc_alignment_score(hsp.bits, query_length, hsp.align_length)
            alignment_names.append(subject_name)
            alignment_scores.append(score)
            identities.append(identity)
            coverages.append(coverage)

    LOGGER.debug(f"processed {len(blast_record.alignments)} alignments")

    blast_result = {}
    for name, score, identity, coverage in zip(alignment_names, alignment_scores, identities, coverages):
        blast_result[name] = {'alignment_score': score,
                              'identity': identity,
                              'coverage': coverage}

    return blast_result

def _calc_alignment_score(bitscore, query_length, subject_length):
    two = Decimal(2)
    bitscore = Decimal(bitscore)
    x = np.power(two, -bitscore) * query_length * subject_length
    alignment_score = int(-np.log10(x))
    return alignment_score


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    query_seq = 'MSKHIGIFGLGAMGTALAAKYLEHGYKTSVWNRTTAKAIPLVEQGAKLASTISEGVNANDLIIICLLNNQVVEDALRDALQTLPSKTIVNLTNGTPNQARKLADFVTSHGARYIHGGIMAVPTMIGSPHAVLLYSGESLELFQSIESHLSLLGMSKYLGTDAGSASLHDLALLSGMYGLFSGFLHAVALIKSGQDTSTTATGLLPLLTPWLSAMTGYLSSIAKQIDDGDYATQGSNLGMQLAGVENIIRAGEEQRVSSQMILPIKALIEQAVGEGHGGEDLSALIEYFKVGKNVD'

    t0 = time.time()
    make_blast_db()
    blast_record = run_blast(query_seq)
    result = process_blast_record(query_seq, blast_record)
    t1 = time.time()

    print(result)
    print(f"Time = {round(t1-t0, 2)}")