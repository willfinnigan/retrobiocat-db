from retrobiocat_web.mongo.model_queries.activity_queries import activity_in_paper


def data_dict_from_activity_obj(act, enzyme_alt_names_to_use):
    data = {'_id': str(act.id),
            'reaction': act.reaction,
            'enzyme_name': act.enzyme_name,
            'substrate_1_smiles': act.substrate_1_smiles,
            'substrate_2_smiles': act.substrate_2_smiles,
            'product_1_smiles': act.product_1_smiles,
            'temperature': act.temperature,
            'ph': act.ph,
            'solvent': act.solvent,
            'other_conditions': act.other_conditions,
            'notes': act.notes,
            'reaction_vol': act.reaction_vol,
            'formulation': act.formulation,
            'biocat_conc': act.biocat_conc,
            'kcat': act.kcat,
            'km': act.km,
            'mw': act.mw,
            'substrate_1_conc': act.substrate_1_conc,
            'substrate_2_conc': act.substrate_2_conc,
            'specific_activity': act.specific_activity,
            'conversion': act.conversion,
            'conversion_time': act.conversion_time,
            'selectivity': act.selectivity,
            'categorical': act.categorical}

    if act.binary:
        data['binary'] = 1
    else:
        data['binary'] = 0

    if act.enzyme_name in enzyme_alt_names_to_use:
        data['enzyme_name'] = enzyme_alt_names_to_use[act.enzyme_name]

    data = {k: v for k, v in data.items() if v is not None}

    return data

def get_activity_data(paper, enzyme_alt_names_to_use):

    activity_data = []
    for i, act_obj in enumerate(activity_in_paper(paper)):
        data = data_dict_from_activity_obj(act_obj, enzyme_alt_names_to_use)
        data['n'] = i+1
        activity_data.append(data)

    return activity_data