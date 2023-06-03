from unidecode import unidecode
from retrobiocat_web.mongo.models.biocatdb_models import Sequence, EnzymeType
import mongoengine as db
from Bio.Seq import Seq

INVALID_NAME_CHARS = [".", "(", ")", "'", "/", ","]

def sanitise_string(text):
    for char in text:
        if char in INVALID_NAME_CHARS:
            text = text.replace(char, '')

    return unidecode(str(text))

def sanitise_sequence(seq_string):
    seq_string = seq_string.replace('\n', '')
    seq_string = seq_string.replace(' ', '')
    seq_string = seq_string.upper()
    return seq_string

def sequence_check(sequence):
    amino_acids_list = ['X', 'V', 'G', 'F', 'E', 'N', 'P', 'Q', 'M', 'K', 'T', 'S', 'W', 'A', 'R', 'D', 'L', 'Y', 'H',
                        'I', 'C', '*']

    bad_chars = []
    for letter in sequence:
        if letter.upper() not in amino_acids_list:
            bad_chars.append(letter)

    return bad_chars

def check_seq_dict(seq_dict):
    """Checks the sequence dictionary for any problems"""

    issues = []

    # check that enzyme_name and enzyme_type are present
    if not any(key in seq_dict for key in ['enzyme_name', 'enzyme_type']):
        issues.append('Must include a enzyme_name and enzyme_type')
        return False, issues

    # check name isn't blank
    if seq_dict['enzyme_name'] == '':
        issues.append(f"Sequence must have a name")

    # check name doesn't exist
    if Sequence.objects(enzyme_name=seq_dict['enzyme_name']).count() != 0:
        issues.append(f"Name already exists")

    # check enzyme_type is valid
    if EnzymeType.objects(enzyme_type=seq_dict['enzyme_type']).count() == 0:
        issues.append(f"Enzyme type {seq_dict['enzyme_type']} does not exist")

    # check for invalid chars in name
    if any(l in seq_dict['enzyme_name'] for l in INVALID_NAME_CHARS):
        issues.append(f"Invalid chars in name:  ({INVALID_NAME_CHARS})")

    # check if sequence uses invalid chars
    if sequence_check(seq_dict.get('sequence', '')) == False:
        issues.append(f"Amino acid sequence for {seq_dict['enzyme_name']} uses incorrect amino acid characters")

    return issues

def clean_seq_lists(seq_obj):
    def remove_empty_space(seq_list):
        if ' ' in seq_list:
            seq_list.remove(' ')
        if '' in seq_list:
            seq_list.remove('')
        return seq_list

    def remove_duplicates(seq_list):
        seq_list = list(set(seq_list))
        return seq_list

    seq_obj.other_names = remove_empty_space(seq_obj.other_names)
    seq_obj.other_names = remove_duplicates(seq_obj.other_names)
    seq_obj.other_identifiers = remove_empty_space(seq_obj.other_identifiers)
    seq_obj.other_identifiers = remove_duplicates(seq_obj.other_identifiers)

def check_seq_has_protein_seq_or_unavailable(seq):
    """Returns False if protein has no sequence and isn't marked unavailable"""

    if seq.sequence == '' or seq.sequence == " " or seq.sequence == None:
        if seq.sequence_unavailable != True:
            return False

    return True

def sequence_name_check(name):
    """ Return False if this name already exists, either as an enzyme_name or in other_names """

    # if blank name, don't allow
    if name == '':
        return False

    name_query = db.Q(enzyme_name=name)
    other_names_query = db.Q(other_names=name)
    other_names_data_query = db.Q(other_names_data__name=name)

    if Sequence.objects(name_query | other_names_query | other_names_data_query ).first() == None:
        return False

    return True


def is_sequence_DNA(sequence):
    return set(list(sequence.upper())) == {'A', 'T', 'G', 'C'}

def if_sequence_is_DNA_then_convert_to_protein(sequence):
    if is_sequence_DNA(sequence):
        dna_seq = Seq(sequence.upper())
        sequence = str(dna_seq.translate())
    return sequence


if __name__ == '__main__':
    print(is_sequence_DNA('atgcccatg'))
    print(if_sequence_is_DNA_then_convert_to_protein('atgcccatg'))
    print(is_sequence_DNA('MGATCYQ'))
    print(if_sequence_is_DNA_then_convert_to_protein('MGATCYQ'))