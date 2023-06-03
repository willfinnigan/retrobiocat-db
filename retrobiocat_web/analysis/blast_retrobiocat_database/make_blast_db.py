import os
import shutil
import subprocess as sp

from retrobiocat_web.analysis.blast_retrobiocat_database.blast_dirs import BLAST_DB, BLAST_DIR
from retrobiocat_web.logging import add_logger
from retrobiocat_web.mongo.model_queries.sequence_queries import sequences_with_protein_sequence

LOGGER = add_logger('Make_BioCatDB_Blast_DB', level='WARNING')

def make_blast_db():
    """ Create a blast database for a single enzyme type """
    LOGGER.debug(f"Making blast database for all sequences")

    _recreate_blast_dir()
    _make_db_fasta()

    command = f"makeblastdb -in {BLAST_DB} -dbtype prot"
    sp.run(command, shell=True)

def _recreate_blast_dir():
    if os.path.exists(BLAST_DIR):
        shutil.rmtree(BLAST_DIR)
    os.mkdir(BLAST_DIR)

def _make_db_fasta():
    """ Create a fasta file containing all the sequences of an enzyme type """

    seqs = sequences_with_protein_sequence()

    with open(f"{BLAST_DB}", 'w') as file:
        for seq in list(seqs):
            name = seq.enzyme_name
            seq = seq.sequence.replace('\n', '')

            file.write(f'>{name}\n')
            file.write(f"{seq}\n")




if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    make_blast_db()
