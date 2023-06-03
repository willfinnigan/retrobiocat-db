from retrobiocat_web.mongo.model_queries import sequence_queries
from retrobiocat_web.mongo.model_queries.activity_queries import activity_for_sequence


def get_sequence_table_data_for_paper(paper):
    seqs = sequence_queries.seqs_from_paper(paper, select_related=True)
    enzyme_data = parse_sequence_data_to_tabulator(seqs)
    enzyme_alt_names_to_use = sequence_queries.get_alt_naming_for_paper(paper)
    enzyme_data = add_paper_naming(enzyme_data, enzyme_alt_names_to_use)
    return enzyme_data

def get_sequence_table_data_for_enzyme_type(enzyme_type_abbreviation):
    seqs = sequence_queries.seq_objs_of_enzyme_type(enzyme_type_abbreviation, select_related=True)
    enzyme_data = parse_sequence_data_to_tabulator(seqs)
    return enzyme_data

def get_sequence_table_data_for_all_enzymes():
    enzyme_data = sequence_queries.all_seqs(select_related=True)
    enzyme_data = parse_sequence_data_to_tabulator(enzyme_data)
    return enzyme_data

def parse_sequence_data_to_tabulator(sequences, add_activity_counts=False, trim_sequences=True):
    enzyme_data = []

    for seq in sequences:
        data = {'_id': str(seq.id),
                'enzyme_type': seq.enzyme_type,
                'enzyme_name': seq.enzyme_name,
                'sequence': seq.sequence,
                'sequence_unavailable': str(seq.sequence_unavailable).replace('False', ''),
                'accession': seq.accession,
                'pdb': seq.pdb,
                'mutant_of': seq.mutant_of,
                'notes': seq.notes,
                'papers': seq.papers,
                }

        if seq.owner is None:
            data['owner'] = ''
        else:
            data['owner'] = f"{seq.owner.first_name} {seq.owner.last_name}"

        if seq.reviewed == True:
            data['reviewed'] = 'True'
        else:
            data['reviewed'] = ''

        other_names_data = seq.other_names_data
        data['other_names'] = [other_data['name'] for other_data in other_names_data]

        data['papers'] = len(seq.papers)

        if trim_sequences == True and len(data['sequence']) > 50:
            data['sequence'] = data['sequence'][0:50] + '...'

        if add_activity_counts:
            data['activity'] = activity_for_sequence(seq, count_only=True)

        enzyme_data.append(data)

    return enzyme_data

def add_paper_naming(enzyme_data, enzyme_alt_names_to_use):
    for i, row in enumerate(enzyme_data):
        default_name = row['enzyme_name']
        if default_name in enzyme_alt_names_to_use:
            enzyme_data[i]['paper_name'] = enzyme_alt_names_to_use[default_name]
        else:
            enzyme_data[i]['paper_name'] = ''

    return enzyme_data

def enzyme_names_to_paper_names(enzyme_names, enzyme_alt_names_to_use):
    paper_names = []
    for name in enzyme_names:
        if name in enzyme_alt_names_to_use:
            paper_names.append(enzyme_alt_names_to_use[name])
        else:
            paper_names.append(name)
    return paper_names