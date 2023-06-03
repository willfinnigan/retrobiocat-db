import numpy as np

formulation_standardisation_dict = {'lysate': 'lysate',
                                    'crude lysate': 'lysate',
                                    'clarified cell extract': 'lysate',
                                    'crude': 'lysate',
                                    'crude extract': 'lysate',
                                    'cell-free extract': 'lysate',
                                    'crude cell free extract': 'lysate',
                                    'crude lyase': 'lysate',
                                    'lyate': 'lysate',

                                    'purified': 'purified',
                                    'pure enzyme': 'purified',
                                    'purified enzyme': 'purified',
                                    'pure protein': 'purified',
                                    'pure': 'purified',
                                    'purified protein': 'purified',
                                    'Purified enzyme': 'purified',
                                    'purified proteinpurified protein': 'purified',
                                    'purifed protein': 'purified',

                                    'whole cell': 'whole cell',
                                    'whole cells': 'whole cell',
                                    'whole-cell': 'whole cell',
                                    'e.coli cells': 'whole cell',
                                    'wet cells': 'whole cell',
                                    'whole cell (e. coli)': 'whole cell',

                                    'lyophilised cell extract': 'lyophilised cell extract',
                                    'lyophilized cell-free extract': 'lyophilised cell extract',
                                    'lyophilised cells': 'lyophilised cell extract',
                                    'cce': 'lyophilised cell extract',
                                    'crude powder': 'lyophilised cell extract',
                                    'cfe': 'lyophilised cell extract',
                                    'cell free extract': 'lyophilised cell extract',
                                    'lyophilised lysate': 'lyophilised cell extract',

                                    'lyophilized cells\xa0E. coli': 'lyophilised cells',
                                    'lyophisised cells': 'lyophilised cells',
                                    'lyophilized cells e. coli': 'lyophilised cells',
                                    ' lyophilised cells': 'lyophilised cells',

                                    }


def remove_trailing_space(string):
    if not isinstance(string, str):
        return string

    if len(string) == 0:
        return string

    while string[-1] == " ":
        string = string[:-1]
        if len(string) == 0:
            return string

    while string[1] == " ":
        string = string[1:]
        if len(string) == 0:
            return string

    return string

def make_lowercase(string):
    if not isinstance(string, str):
        return string
    return string.lower()

def empty_string_to_none(string):
    if not isinstance(string, str):
        return None
    if len(string) == 0:
        return None
    return string

def standardise_forumation(string):
    string = make_lowercase(string)
    string = remove_trailing_space(string)
    string = empty_string_to_none(string)
    if string in formulation_standardisation_dict:
        return formulation_standardisation_dict[string]
    else:
        print(string)
        return string




if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    from retrobiocat_web.analysis.data_query import get_data

    data_query = get_data.DataQuery(log_level=1, remove_negative=False)
    formulations = list(data_query.get_activity_df()['formulation'].unique())
    print(formulations)

    std_formulations = [standardise_forumation(s) for s in formulations]
    print(list(set(std_formulations)))