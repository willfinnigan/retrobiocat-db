import json


def get_seq_table_entry(seq):
    seq_dict = json.loads(seq.to_json())
    seq_dict['_id'] = str(seq_dict['_id'])
    seq_dict.pop("edits_by", None)
    seq_dict.pop("papers", None)
    seq_dict.pop("added_by", None)
    seq_dict['owner'] = "refresh to load"

    return seq_dict