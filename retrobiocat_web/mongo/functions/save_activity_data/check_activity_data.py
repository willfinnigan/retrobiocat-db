from rdkit import Chem
import numpy as np


def check_column(data_dict, col):
    issues = []
    if col in data_dict:
        if data_dict[col] == "" or data_dict[col] is None:
            return [f"Row {data_dict['n']}: Must specificy {col}"]
    else:
        return [f"Row {data_dict['n']}: Must specificy {col}"]
    return []

def check_all_required_are_numbers(data_dict):
    cols = ['conversion', 'specific_activity', 'km', 'kcat', 'mw']
    for col in cols:
        if col in data_dict:
            if not isinstance(data_dict[col], (int, float, complex)) or isinstance(data_dict[col], bool):
                if data_dict[col] != '':
                    return [f"Row {data_dict['n']}: {col} data is type string ({data_dict[col]})"]
    return []

def check_required_columns(data_dict):
    issues = []

    required_cols = ['reaction', 'enzyme_name']
    for col in required_cols:
        issues += check_column(data_dict, col)
    return issues


def check_all_have_binary(data_dict):
    if 'binary' not in data_dict:
        return [f"Row {data_dict['n']}: No activity data has been entered - no binary"]
    if data_dict['binary'] != 1 and data_dict['binary'] != 0:
        return [f"Row {data_dict['n']}: No activity data has been entered - no binary"]
    return []

def check_seqs_are_defined(data_dict, paper_seqs):
    issues = []
    names = []
    for seq in paper_seqs:
        names.append(seq.enzyme_name)
        for data in seq.other_names_data:
            names.append(data.name)

    enzyme_name = data_dict.get('enzyme_name', None)
    if enzyme_name not in names:
        issues.append(f"Enzymes must be defined in sequence tab first - {enzyme_name}")

    return issues





