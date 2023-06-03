from retrobiocat_web.mongo.models.biocatdb_models import Sequence
import mongoengine as db


def does_sequence_match_check(enzyme_type, seq_to_check):
    seqs = Sequence.objects(db.Q(enzyme_type=enzyme_type) &
                            db.Q(sequence__ne=None) &
                            db.Q(sequence__ne=''))

    for seq in seqs:
        if seq.sequence == seq_to_check:
            return seq.enzyme_name

    return ''


def find_identical_sequences(enzyme_type):
    """Return lists of all sequences which match"""

    seqs = Sequence.objects(db.Q(enzyme_type=enzyme_type) &
                            db.Q(sequence__ne=None) &
                            db.Q(sequence__ne=''))


    # keep a dict of every protein sequence seen.  Add names to values of dict
    seq_dict = {}
    for seq in seqs:
        if seq.sequence not in seq_dict:
            seq_dict[seq.sequence] = [seq.enzyme_name]
        else:
            seq_dict[seq.sequence].append(seq.enzyme_name)

    # get the lists of seq_dict names that contain more than a single name
    matches = []
    for key, value in seq_dict.items():
        if len(value) > 1:
            matches.append(value)

    return matches


if __name__ == "__main__":
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection(host='138.68.135.53', database='test')

    et = 'AlOx'

    matches = find_identical_sequences(et)

    print(matches)