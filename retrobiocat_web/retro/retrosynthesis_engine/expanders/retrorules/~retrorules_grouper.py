
class RetroRules_Grouper():

    def group(self, precursor_dict, rule_scores, rule_mx_ids):
        groups, product_map = self._sort_reactions_into_groups(precursor_dict)
        groups, product_map = self._groups_to_uuid(groups, product_map)
        new_precursor_dict = self._get_new_precursor_dict(product_map)
        metadata = self._metadata_uuid(groups, rule_scores, rule_mx_ids)
        return new_precursor_dict, metadata

    def no_group(self, precursor_dict, rule_scores, rule_mx_ids):
        new_names_precursor_dict = {}
        metadata = {}
        for rule_id in precursor_dict:
            rr = rule_id[-4:]
            uuid = rule_id[6:10]
            unique = f"RR-Biosyn-{rr}-{uuid}"

            new_names_precursor_dict[unique] = precursor_dict[rule_id]
            metadata[unique] = {rule_id: {'Score': round(rule_scores[rule_id], 2),
                                          'MNXR_ids': rule_mx_ids[rule_id]}}

    def _sort_reactions_into_groups(self, precursor_dict):
        groups = {}
        product_map = {}
        for reaction, multi_products in precursor_dict.items():
            for products in multi_products:
                productKey = str(sorted(products))
                if productKey not in groups:
                    groups[productKey] = []
                groups[productKey].append(reaction)
                product_map[productKey] = products
        return groups, product_map

    def _groups_to_uuid(self, groups, product_map):
        new_groups = {}
        new_product_map = {}
        for productKey in groups:
            rr = groups[productKey][0][-4:]
            uuid = groups[productKey][0][6:10]
            unique = f"RR-Biosyn-{rr}-{uuid}"
            new_groups[unique] = groups[productKey]
            new_product_map[unique] = product_map[productKey]
        return new_groups, new_product_map

    def _metadata_uuid(self, groups, rule_scores, rule_mx_ids):
        metadata = {}
        for uuid in groups:
            uuid_meta_dict = {}
            for rule_id in groups[uuid]:
                uuid_meta_dict[rule_id] = {'Score': round(rule_scores[rule_id],2),
                                           'MNXR_ids': rule_mx_ids[rule_id]}
            metadata[uuid] = uuid_meta_dict
        return metadata

    def _get_new_precursor_dict(self, product_map):
        new_precursor_dict = {}
        for uuid, products in product_map.items():
            if uuid not in new_precursor_dict:
                new_precursor_dict[uuid] = []
            new_precursor_dict[uuid].append(product_map[uuid])
        return new_precursor_dict


if __name__ == '__main__':

    test_ids = {'RR-02-7b86f0392ef772c6-02-F': ['1'],
                 'RR-02-ecb9868584935ce5-02-F': ['2'],
                 'RR-02-bf1b054a87db912e-02-F': ['3'],
                 'RR-02-3cd3056034952c90-02-F': ['5'],
                 'RR-02-8a61c601b46ddd02-02-F': ['4']}

    test_scores = {'RR-02-7b86f0392ef772c6-02-F': 1,
                'RR-02-ecb9868584935ce5-02-F': 0.9,
                'RR-02-bf1b054a87db912e-02-F': 0.5,
                'RR-02-3cd3056034952c90-02-F': 0.4,
                'RR-02-8a61c601b46ddd02-02-F': 0.4}

    test_dict = {'RR-02-7b86f0392ef772c6-02-F': [['CCCC=O']],
                 'RR-02-ecb9868584935ce5-02-F': [['CCCCOC(C)=O']],
                 'RR-02-bf1b054a87db912e-02-F': [['CCCC=O']],
                 'RR-02-3cd3056034952c90-02-F': [['CCCCOC(C)=O']],
                 'RR-02-8a61c601b46ddd02-02-F': [['CCCCCCCCCCCCCC(=O)OCC']]}

    grouper = RetroRules_Grouper()
    new_precursor_dict, new_metadata = grouper.group(test_dict, test_ids, test_scores)