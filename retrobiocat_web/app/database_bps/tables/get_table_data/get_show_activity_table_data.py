import json
from retrobiocat_web.mongo.model_queries.get_activity_data import get_activity_data

def activity_data_from_args(args):
    enzyme_type = args.get('enzyme_type', None)

    if 'enzyme_name' in args:
        enzyme_names = [args['enzyme_name']]
    else:
        enzyme_names = args.get('enzyme_names', None)
        if enzyme_names is not None:
            enzyme_names = json.loads(enzyme_names)

    reaction = args.get('reaction', None)
    paper_id = args.get('paper_id', None)
    solvent = args.get('solvent', None)
    formulation = args.get('formulation', None)
    only_reviewed = args.get('reviewed', True)
    include_negative = args.get('include_negative', True)
    include_auto_generated = args.get('include_auto_generated', False)
    ph_left = args.get('ph_left', None)
    ph_right = args.get('ph_right', None)
    temp_left = args.get('temp_left', None)
    temp_right = args.get('temp_right', None)

    smiles = args.get('smiles', None)
    #if smiles is not None:
    #    smiles = json.loads(smiles)
    if 'smi_col' in args:
        smi_col = args.get('smi_col', 'product_1_smiles')
    else:
        smi_col = args.get('mol_type', 'product_1_smiles')

    print(args)

    return get_activity_data(enzyme_type=enzyme_type, reaction=reaction,
                              paper_id=paper_id, enzyme_names=enzyme_names,
                              formulation=formulation, solvent=solvent,
                              ph_left=ph_left, ph_right=ph_right,
                              temp_left=temp_left, temp_right=temp_right,
                              smi_col=smi_col, smiles=smiles,
                              only_reviewed=only_reviewed,
                              include_negative=include_negative,
                              include_auto_generated=include_auto_generated,
                              as_pymongo=True)

