
def get_level(amount, levels):
    cat = 'None'
    current_cat_amount = 0

    for category, minimum in levels.items():
        if amount >= minimum:
            if minimum > current_cat_amount:
                cat = category
                current_cat_amount = minimum
    return cat

def binary_from_category(data_dict):
    if data_dict.get('binary', '') == '':
        if data_dict.get('categorical', '') != '':
            if data_dict['categorical'] == 'None':
                data_dict['binary'] = 0
            else:
                data_dict['binary'] = 1
    return data_dict

def category_from_conversion(data_dict, levels):
    if type(data_dict.get('conversion','')) != str:
        if data_dict.get('categorical', '') == '':
            data_dict['categorical'] = get_level(data_dict['conversion'], levels)
    return data_dict

def category_from_sa(data_dict, levels):

    if type(data_dict.get('specific_activity','')) != str:
        if data_dict.get('categorical', '') == '':
            data_dict['categorical'] = get_level(data_dict['specific_activity'], levels)
    return data_dict

def sa_from_kinetics(data_dict, default_conc):
    if type(data_dict.get('km', '')) != str and type(data_dict.get('kcat', '')) != str and type(data_dict.get('mw', '')) != str:
        if type(data_dict.get('specific_activity','')) == str:
            data_dict['specific_activity'] = calc_sa(data_dict['kcat'], data_dict['km'], data_dict['mw'], default_conc)
            data_dict['substrate_1_conc'] = f"{default_conc} mM - calculated from kinetics"
    return data_dict

def calc_sa(kcat, km, mw, s):
    if kcat == 0 or s == 0 or km == 0:
        umol_min_mg = 0
    else:
        umols_enz = 1000/mw
        vmax = kcat*umols_enz
        umol_min_mg = vmax * (s/(s+km))
        umol_min_mg = round(umol_min_mg, 2)
    return umol_min_mg

def remove_empty_columns(data_dict):
    for key in list(data_dict.keys()):
        if data_dict[key] == '' and key != '_id':
            data_dict.pop(key)
    return data_dict

def make_sure_binary_is_int(data_dict):
    try:
        if 'binary' in data_dict:
            data_dict['binary'] = int(data_dict['binary'])
    except:
        if 'binary' in data_dict:
            data_dict.pop('binary')
    return data_dict
