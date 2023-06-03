from distutils.util import strtobool

from retrobiocat_web.analysis.data_query import get_data


def data_query_from_args(args):
    reaction = args.get('reaction', None)
    enzyme_type = args.get('enzyme_type', None)
    enzyme_names = args.get('enzyme_names', None)
    paper_id = args.get('paper_id', None)
    only_reviewed = bool(strtobool(args.get('only_reviewed', 'True')))
    mol_type = args.get('mol_type', 'substrate_1_smiles')
    smiles = args.get('smiles', None)

    dq = get_data.DataQuery(reaction=reaction,
                            enzyme_type=enzyme_type,
                            enzyme_names=enzyme_names,
                            paper_id=paper_id,
                            only_reviewed=only_reviewed,
                            smiles=smiles,
                            smi_col=mol_type)

    return dq

if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    args = {'enzyme_type': 'AmOx', 'mol_type': 'product_1_smiles'}
    dq = data_query_from_args(args)
    print(dq.unique_smiles())