from retrobiocat_web.analysis.data_query import get_data
import copy

# 1. unique substrates
# 2. order substrates // cluster
# 3. max number?

cat_num = {'High': 5,
           'Medium': 4,
           'Low': 3,
           'Active': 2,
           'Inactive': 1}

def get_smi_cats(activity_data, smi_col):

    smi_cats = {}

    for act_dict in activity_data:
        smi = act_dict.get(smi_col, None)
        if smi is not None and smi != '':
            level = act_dict.get('categorical', 'None')
            if level == 'None' or level == '':
                binary = act_dict.get('binary', 'False')
                if binary == 'False':
                    level = 'Inactive'
                else:
                    level = 'Active'

            if smi is None:
                pass
            elif cat_num.get(smi_cats.get(smi, None), 0) < cat_num[level]:
                smi_cats[smi] = level

    return smi_cats

def order_smiles(list_smis, num_enzymes):

    if len(list_smis) == 0:
        return []

    # group by number of enzymes, but group 2-3, 4-5 ect together
    grouped_by_number = {}
    for smi in list_smis:
        num = num_enzymes[smi]
        if num not in grouped_by_number:
            grouped_by_number[num] = []
        grouped_by_number[num].append(smi)

    # group 2-3, 4-5 ect together
    combined_groups = {}
    if 1 in grouped_by_number:
        combined_groups[1] = grouped_by_number[1]

    for i in range(2, max(grouped_by_number.keys())+1, 2):
        group_1 = grouped_by_number.get(i, [])
        group_2 = grouped_by_number.get(i+1, [])
        combined_groups[i] = group_1 + group_2

    new_order = []
    for num, group_smis in combined_groups.items():
        if len(group_smis) > 3:
            subs, _ = get_data.order_smis(group_smis)
            new_order += subs
        else:
            new_order += group_smis

    new_order.reverse()

    return new_order


def sort_smis_to_cats(smi_cats, num_enzymes):
    cats = {'High': [],
            'Medium': [],
            'Low': [],
            'Active': [],
            'Inactive': []}

    for smi, level in smi_cats.items():
        cats[level].append(smi)

    for level, list_smis in cats.items():
        preordered_smis = sorted(list_smis, key=lambda x: num_enzymes[x], reverse=True)
        ordered_smis = order_smiles(preordered_smis, num_enzymes)
        cats[level] = ordered_smis

    return cats

def get_svg_cats(cats, molSize=(250, 100)):
    svg_cats = {}

    for level, list_smis in cats.items():
        svg_map, svg_url_map = get_data.get_mol_svgs(list_smis, molSize=molSize)
        svg_cats[level] = svg_map

    return svg_cats


def get_substrate_info_dict(activity_data, smi_col, enz_col='enzyme_name', paper_col='paper'):
    smi_info = {}
    smi_stats = {}
    num_enzymes = {}

    for act_dict in activity_data:
        smi = act_dict.get(smi_col, None)

        if smi is not None and smi != '':
            if smi not in smi_info:
                smi_info[smi] = []
            smi_info[smi].append(act_dict)

            if smi not in smi_stats:
                smi_stats[smi] = {'enzymes': set(), 'papers': set()}

            smi_stats[smi]['enzymes'].add(act_dict.get(enz_col, None))
            smi_stats[smi]['papers'].add(act_dict.get(paper_col, None))

    for smi in smi_info:
        smi_stats[smi]['enzymes'] = list(smi_stats[smi]['enzymes'])
        smi_stats[smi]['num_enzymes'] = len(smi_stats[smi]['enzymes'])
        smi_stats[smi]['papers'] = list(smi_stats[smi]['papers'])
        smi_stats[smi]['num_papers'] = len(smi_stats[smi]['papers'])
        num_enzymes[smi] = len(smi_stats[smi]['enzymes'])

    return smi_info, smi_stats, num_enzymes


def get_substrate_table(activity_data, smi_col='substrate_1_smiles'):
    smi_info, smi_stats, num_enzymes = get_substrate_info_dict(activity_data, smi_col)
    smi_cats = get_smi_cats(activity_data, smi_col)
    cats = sort_smis_to_cats(smi_cats, num_enzymes)
    svg_cats = get_svg_cats(cats, molSize=(250, 100))


    return svg_cats, smi_stats






if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection(host='138.68.135.53', database='test')

    import time

    t0 = time.time()
    dq = get_data.DataQuery(enzyme_type='IRED',
                            reaction='Reductive amination',
                            smi_col='substrate_1_smiles')
    t1 = time.time()
    svg_map = get_substrate_table(dq.activity_data, smi_col='substrate_1_smiles')
    t2 = time.time()

    print(f'Time for query = {round(t1-t0,2)} seconds')
    print(f'Time for svgs = {round(t2 - t1, 2)} seconds')
    print(f'Time total = {round(t2 - t0, 2)} seconds')
    #dq.activity_data

    from rdkit import Chem
    import gzip

    #smis = dq.unique_smiles()
    #mols = [Chem.MolFromSmiles(smi) for smi in smis]
    #for mol in mols:
        #mol.SetProp('MyProp', 'foo')

    #writer = Chem.SDWriter('out.sdf')
    #writer.SetProps(['MyProp'])
    #for mol in mols:
     #   writer.write(mol)
    #writer.close()

