from typing import List

from retrobiocat_web.app.database_bps.curation.functions.apply_reactions_during_activity_curation.apply_reaction import \
    retro_rxn, fwd_rxn


class ActivityTable_Rxn_Apply():

    def __init__(self):
        pass

    def run(self, rows):
        product_groups, substrate_groups = self.group(rows)
        product_results = self.get_product_rows(product_groups)
        substrate_results = self.get_substrate_rows(substrate_groups)
        new_rows = product_results + substrate_results
        return new_rows

    def group(self, rows):
        product_groups = {}
        substrate_groups = {}
        for row in rows:

            reaction_name = row.get('reaction', None)
            n = row.get('n', None)
            if reaction_name is None or n is None:
                continue

            product = row.get('product_1_smiles', '')
            substrate_1 = row.get('substrate_1_smiles', '')
            substrate_2 = row.get('substrate_2_smiles', '')

            if product != '' and substrate_1 == '':
                product_tuple = (product, reaction_name)
                if product_tuple not in product_groups:
                    product_groups[product_tuple] = []
                product_groups[product_tuple].append(row)

            if substrate_1 != '' and product == '':
                substrate_tuple = (substrate_1, substrate_2, reaction_name)
                if substrate_tuple not in substrate_groups:
                    substrate_groups[substrate_tuple] = []
                substrate_groups[substrate_tuple].append(row)

        return product_groups, substrate_groups

    def get_product_rows(self, product_groups) -> List[dict]:
        results = []
        for product_tuple, rows in product_groups.items():
            product = product_tuple[0]
            reaction_name = product_tuple[1]

            outcomes = retro_rxn(reaction_name, product, True, True)

            for row in rows:
                new_row = {'n': row['n']}
                if len(outcomes) == 1:
                    if len(outcomes[0]) == 1:
                        new_row['substrate_1_smiles'] = outcomes[0][0]
                    elif len(outcomes[0]) == 2:
                        new_row['substrate_1_smiles'] = outcomes[0][0]
                        new_row['substrate_2_smiles'] = outcomes[0][1]
                    results.append(new_row)

        return results

    def get_substrate_rows(self, substrate_groups) -> List[dict]:
        results = []
        for substrate_tuple, rows in substrate_groups.items():
            substrate_1 = substrate_tuple[0]
            substrate_2 = substrate_tuple[1]
            reaction_name = substrate_tuple[2]

            outcomes = fwd_rxn(substrate_1, substrate_2, reaction_name, True, True)

            for row in rows:
                new_row = {'n': row['n']}
                if len(outcomes) == 1:
                    if len(outcomes[0]) == 1:
                        new_row['product_1_smiles'] = outcomes[0][0]
                        results.append(new_row)

        return results


