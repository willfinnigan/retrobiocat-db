from retrobiocat_web.app.database_bps.curation.functions.apply_reactions_during_activity_curation.apply_reaction import \
    retro_rxn, fwd_rxn
from retrobiocat_web.mongo.model_queries.reaction_queries import reactions_of_type
from retrobiocat_web.mongo.model_queries.sequence_queries import seq_obj_from_name


def get_reaction_name_from_product(product_smi, possible_reactions):
    product_smi = product_smi.replace('@', '')
    for reaction in possible_reactions:
        outcomes = retro_rxn(reaction.name, product_smi, True, True)
        if len(outcomes) != 0:
            return reaction.name
    return ''

def get_reaction_name_from_substrates(substrate_1, substrate_2, possible_reactions):
    substrate_1 = substrate_1.replace('@', '')
    substrate_2 = substrate_2.replace('@', '')

    for reaction in possible_reactions:
        outcomes = fwd_rxn(substrate_1, substrate_2, reaction.name, True, True)
        if len(outcomes) != 0:
            return reaction.name

    return ''


class ReactionName_AutoDetector():

    def __init__(self):
        self.enzyme_type_possible_reactions = {}

    def run(self, rows):
        rows = self.add_enzyme_types(rows)
        rows_by_product, rows_by_substrates = self.combine_rows_by_smis_and_enzyme_types(rows)
        product_rows = self.detect_name_by_product(rows_by_product)
        substrate_rows = self.detect_name_by_substrates(rows_by_substrates)
        new_rows = product_rows + substrate_rows
        return new_rows

    def get_possible_reactions(self, enzyme_type):
        if enzyme_type not in self.enzyme_type_possible_reactions:
            self.enzyme_type_possible_reactions[enzyme_type] = reactions_of_type([enzyme_type])
        return self.enzyme_type_possible_reactions[enzyme_type]

    def add_enzyme_types(self, rows):
        enzyme_names = []
        for i, row_data in enumerate(rows):
            name = row_data.get('enzyme_name', None)
            if name is not None:
                enzyme_names.append(name)

        enzyme_types = {}
        for name in enzyme_names:
            seq_obj = seq_obj_from_name(name, include_other_names=True)
            if seq_obj is not None:
                enzyme_types[name] = seq_obj.enzyme_type

        for i, row_data in enumerate(rows):
            name = row_data.get('enzyme_name', None)
            if name in enzyme_types:
                rows[i]['enzyme_type'] = enzyme_types[name]
        return rows

    def combine_rows_by_smis_and_enzyme_types(self, rows):
        rows_by_product = {}
        rows_by_substrates = {}

        for i, row_data in enumerate(rows):
            enzyme_type = row_data.get('enzyme_type', '')
            product = row_data.get('product_1_smiles', '')
            substrate_1 = row_data.get('substrate_1_smiles', '')
            substrate_2 = row_data.get('substrate_2_smiles', '')

            if product != '':
                product_enzyme_type = (product, enzyme_type)
                if product_enzyme_type not in rows_by_product:
                    rows_by_product[product_enzyme_type] = []
                rows_by_product[product_enzyme_type].append(row_data)

            elif substrate_1 != '':
                substrates_enzyme_type = (substrate_1, substrate_2, enzyme_type)
                if substrates_enzyme_type not in rows_by_substrates:
                    rows_by_substrates[substrates_enzyme_type] = []
                rows_by_substrates[substrates_enzyme_type].append(row_data)

        return rows_by_product, rows_by_substrates


    def detect_name_by_product(self, rows_by_product):
        product_rows = []
        for p_et_tuple, rows in rows_by_product.items():
            product = p_et_tuple[0]
            enzyme_type = p_et_tuple[1]
            possible_reactions = self.get_possible_reactions(enzyme_type)
            reaction_name = get_reaction_name_from_product(product, possible_reactions)

            for i, row_data in enumerate(rows):
                row_data['reaction'] = reaction_name
                product_rows.append(row_data)

        return product_rows

    def detect_name_by_substrates(self, rows_by_substrates):
        substrate_rows = []
        for s1_s2_et_tuple, rows in rows_by_substrates.items():
            substrate_1 = s1_s2_et_tuple[0]
            substrate_2 = s1_s2_et_tuple[1]
            enzyme_type = s1_s2_et_tuple[2]
            possible_reactions = self.get_possible_reactions(enzyme_type)

            reaction_name = get_reaction_name_from_substrates(substrate_1, substrate_2, possible_reactions)

            for i, row_data in enumerate(rows):
                row_data['reaction'] = reaction_name
                substrate_rows.append(row_data)

        return substrate_rows