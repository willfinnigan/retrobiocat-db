from pathlib import Path

"""
Paths for data folders used by all the modules in analysis
"""

all_by_all_blast_folder = str(Path(__file__).parents[0]) + '/analysis_data/all_by_all_blast'
ssn_folder = str(Path(__file__).parents[0]) + '/analysis_data/ssn'
biocatdb_blast_folder = str(Path(__file__).parents[0]) + '/analysis_data/all_seqs'

if __name__ == '__main__':
    print(all_by_all_blast_folder)