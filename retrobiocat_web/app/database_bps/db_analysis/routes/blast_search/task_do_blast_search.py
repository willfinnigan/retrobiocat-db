from retrobiocat_web.analysis.blast_retrobiocat_database.make_blast_db import make_blast_db
from retrobiocat_web.analysis.blast_retrobiocat_database.run_blast import run_blast, process_blast_record
from retrobiocat_web.app.database_bps.tables.get_table_data.get_sequence_table_data import \
    parse_sequence_data_to_tabulator
from retrobiocat_web.mongo.model_queries.sequence_queries import seqs_from_list_of_names


def task_do_blast_search(sequence):
    make_blast_db()
    blast_record = run_blast(sequence)
    blast_result = process_blast_record(sequence, blast_record)
    sequences = seqs_from_list_of_names(list(blast_result.keys()))
    table_data = parse_sequence_data_to_tabulator(sequences, add_activity_counts=True)
    table_data = _add_blast_result(table_data, blast_result)
    return table_data

def _add_blast_result(table_data, result_dict):
    for i, seq_data in enumerate(table_data):
        name = seq_data['enzyme_name']
        table_data[i]['alignment_score'] = result_dict[name]['alignment_score']
        table_data[i]['identity'] = result_dict[name]['identity']
        table_data[i]['coverage'] = result_dict[name]['coverage']
    return table_data

if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    query_seq = 'MSKHIGIFGLGAMGTALAAKYLEHGYKTSVWNRTTAKAIPLVEQGAKLASTISEGVNANDLIIICLLNNQVVEDALRDALQTLPSKTIVNLTNGTPNQARKLADFVTSHGARYIHGGIMAVPTMIGSPHAVLLYSGESLELFQSIESHLSLLGMSKYLGTDAGSASLHDLALLSGMYGLFSGFLHAVALIKSGQDTSTTATGLLPLLTPWLSAMTGYLSSIAKQIDDGDYATQGSNLGMQLAGVENIIRAGEEQRVSSQMILPIKALIEQAVGEGHGGEDLSALIEYFKVGKNVD'

    table_data = task_do_blast_search(query_seq)
    print(table_data)