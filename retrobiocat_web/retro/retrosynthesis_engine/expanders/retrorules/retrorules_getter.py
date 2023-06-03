import sqlite3
from pathlib import Path

import pandas as pd
from rdkit import DataStructs, Chem
from rdkit.Chem import rdChemReactions, rdFingerprintGenerator

from retrobiocat_web.logging import add_logger
from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig

data_folder = str(Path(__file__).parents[4]) + '/data/retrorules'


class RetroRules_Getter():
    fp_gen = rdFingerprintGenerator.GetRDKitFPGenerator()

    def __init__(self, config=None, log_level='WARNING'):
        self.logger = add_logger('RetroBioCat_Expander', level=log_level)

        self.config = config
        if self.config is None:
            self.config = RetrosynthesisConfig()

    def get_similar_smiles_ids(self, target_smi):
        threshold = self.config.rr_threshold

        df = pd.read_hdf(f'{data_folder}/rr_smi.h5')
        mol = Chem.MolFromSmiles(target_smi)
        fp = self.fp_gen.GetFingerprint(mol)
        df['sim'] = DataStructs.BulkTanimotoSimilarity(fp, list(df['fp']))
        df = df[df['sim'] >= threshold]

        smi_ids = list(df['Substrate_ID'])
        sims = list(df['sim'])

        self.logger.debug(f"Smi_Ids: {smi_ids}")
        self.logger.debug(f"Similarities: {sims}")
        smi_id_sim_dict = {smi_id: sim for smi_id, sim in zip(smi_ids, sims)}

        return smi_id_sim_dict

    def get_retrorules_rev(self, smi_id_sim_dict):

        smi_ids = str(list(smi_id_sim_dict.keys()))
        smi_ids = smi_ids.replace('[', '(')
        smi_ids = smi_ids.replace(']', ')')

        diameter = self.config.rr_diameter

        cmd = f'''SELECT *
                  FROM "RETRORULES"
                  WHERE Rule_usage IN ('both', 'reverse')
                    AND Diameter={str(diameter)}
                    AND Substrate_ID IN {smi_ids}
               '''

        conn = sqlite3.connect(f"{data_folder}/retrorules.db")
        cursor = conn.cursor()
        cursor.execute(cmd)
        result = cursor.fetchall()
        conn.close()

        self.logger.debug(f"Rules = {result}")

        return result

    def get_dicts(self, rules, smi_id_similarity_dict):
        rxn_smarts = {}
        scores = {}
        similarities = {}
        rxn_ids = {}
        ec_numbers = {}
        for rule in rules:
            rule_id = rule[1]
            if rule_id not in rxn_ids:
                rxn_ids[rule_id] = []

            rxn_ids[rule_id].append(rule[3])
            rxn_smarts[rule_id] = rule[6]
            scores[rule_id] = rule[14]
            similarities[rule_id] = smi_id_similarity_dict.get(rule[7], 0)

            if rule[15] is not None:
                ec_numbers[rule_id] = str(rule[15]).split(',')
            else:
                ec_numbers[rule_id] = []

            if rule[7] not in list(smi_id_similarity_dict.keys()):
                self.logger.warning(f"Warning - {rule[7]} not in similarity dict")

        return rxn_smarts, scores, rxn_ids, similarities, ec_numbers

    def create_final_dicts(self, rxn_smarts, scores, rxn_ids, similarities, ec_numbers):
        smarts_dict = {}
        metadata = {}
        for rule_id in rxn_smarts:
            smarts_dict[rule_id] = [rxn_smarts[rule_id]]
            metadata[rule_id] = {'score': scores[rule_id],
                                 'similarity': similarities[rule_id],
                                 'rxn_ids': rxn_ids[rule_id],
                                 'ec_numbers': ec_numbers[rule_id],
                                 'rule_id': rule_id}
        return smarts_dict, metadata

    def filter_by_biological_score(self, smarts_dict, metadata):
        score_cutoff = self.config.rr_score_threshold
        to_remove = []
        for name in metadata:
            if metadata[name]['score'] < score_cutoff:
                to_remove.append(name)

        for name in to_remove:
            smarts_dict.pop(name)
            metadata.pop(name)

        self.logger.debug(f'Removed {len(to_remove)} rules below the biological score of {score_cutoff}')

        return smarts_dict, metadata

    def filter_by_combined_score(self, smarts_dict, metadata):
        combined_score_cutoff = self.config.rr_combined_score_threshold
        to_remove = []
        for name in metadata:
            combined_score = metadata[name]['score'] * metadata[name]['similarity']
            if combined_score < combined_score_cutoff:
                to_remove.append(name)

        for name in to_remove:
            smarts_dict.pop(name)
            metadata.pop(name)

        self.logger.debug(f'Removed {len(to_remove)} rules below the combined biological score of {combined_score_cutoff}')

        return smarts_dict, metadata

    def get_rules(self, smi):
        self.logger.info(f'-- Get rules for {smi} --')
        smi_id_sim_dict = self.get_similar_smiles_ids(smi)
        rules = self.get_retrorules_rev(smi_id_sim_dict)
        rxn_smarts, scores, rxn_ids, similarities, ec_numbers = self.get_dicts(rules, smi_id_sim_dict)
        smarts, metadata = self.create_final_dicts(rxn_smarts, scores, rxn_ids, similarities, ec_numbers)
        smarts, metadata = self.filter_by_biological_score(smarts, metadata)
        smarts, metadata = self.filter_by_combined_score(smarts, metadata)
        self.logger.debug(smarts)
        self.logger.debug(metadata)
        return smarts, metadata

    def make_rdkit_rxns(self, smarts_dict):
        rxn_dict = {}
        for rule_id, smarts in smarts_dict.items():
            rxn = rdChemReactions.ReactionFromSmarts(smarts)
            rxn_dict[rule_id] = [rxn]
        return rxn_dict


if __name__ == "__main__":
    test_smi = 'NC[C@H](O)c1ccccc1'
    rr_get = RetroRules_Getter(log_level='DEBUG')
    rr_get.config.rr_diameter = 4
    rr_get.config.rr_threshold = 0.5
    rr_get.config.rr_combined_score_threshold = 0.5
    smarts, metadata = rr_get.get_rules(test_smi)

