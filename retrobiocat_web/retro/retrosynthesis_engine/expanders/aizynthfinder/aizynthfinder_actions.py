from rdkit import Chem
import numpy as np
import pandas as pd
from pathlib import Path

from retrobiocat_web.logging import add_logger
from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig
from retrobiocat_web.retro.retrosynthesis_engine.rule_application.rdchiral.initialization import \
    rdchiralReaction
from retrobiocat_web.retro.retrosynthesis_engine.rule_application.rule_applicator import RuleApplicator
from retrobiocat_web.retro.retrosynthesis_engine.expanders.aizynthfinder import load_keras_models, fingerprints

data_folder = str(Path(__file__).parents[4]) + '/data/aizynthfinder'


class ActionFilter():
    filter_model = None

    def __init__(self, config=None, log_level='WARNING'):
        self.logger = add_logger('AIZynthfinder_Filter', level=log_level)

        self.config = config
        if self.config is None:
            self.config = RetrosynthesisConfig()

    @classmethod
    def load_model(cls):
        if cls.filter_model == None:
            filter_path = data_folder + '/filter_policy_all.hdf5'
            cls.filter_model = load_keras_models.LocalKerasModel(filter_path)

    def feasibility(self, reaction):
        """
        Computes if a given reaction is feasible by giving the reaction fingerprint to a network model
        """

        cuttoff = self.config.aizynth_filter_cutoff

        if self.filter_model == None:
            self.load_model()

        reaction_as_mols = [[Chem.MolFromSmiles(smi) for smi in reaction[0]],
                            [Chem.MolFromSmiles(smi) for smi in reaction[1]]]

        prod_fp = fingerprints.get_mol_fingerprint(reaction_as_mols[0][0], len(self.filter_model))
        prod_fp = prod_fp.reshape([1, len(self.filter_model)])
        reaction_fp = fingerprints.get_reaction_fingerprint(reaction_as_mols, len(self.filter_model))
        reaction_fp = reaction_fp.reshape([1, len(self.filter_model)])
        kwargs = {'input_1': prod_fp, 'input_2': reaction_fp}

        prob = self.filter_model.predict(prod_fp, reaction_fp, **kwargs)[0][0]
        if prob >= cuttoff:
            feasible = True
        else:
            feasible = False
        return feasible, prob

    def filter(self, target, precursor_dict, metadata, log=False):
        if not self.config.aizynth_use_feasability_filter:
            return precursor_dict, metadata

        keep_metadata = {}
        keep_precursor_dict = {}
        count_filtered = 0
        for name, substrates in precursor_dict.items():
            for substrate_set in substrates:
                reaction = [[target], substrate_set]
                feasible, prob = self.feasibility(reaction)
                if feasible is True:
                    keep_metadata[name] = metadata[name]
                    keep_metadata[name]['filter_probability'] = round(float(prob), 5)
                    keep_precursor_dict[name] = substrates
                else:
                    count_filtered += 1
                    if log == True:
                        self.logger.debug(f"{feasible}, Prob. {prob}")

        self.logger.info(f"Filtered {count_filtered} reactions")
        self.logger.info(f"{len(keep_precursor_dict)} reaction remaining")

        return keep_precursor_dict, keep_metadata

class ActionGetter():
    policy_model = None
    templates = None

    def __init__(self, config=None, log_level='WARNING'):
        self.logger = add_logger('AIZynthfinder_Filter', level=log_level)

        self.config = config
        if self.config is None:
            self.config = RetrosynthesisConfig()

        self.actions = {}

    @classmethod
    def load_model(cls):
        if cls.policy_model == None:
            policy_path = data_folder + '/uspto_model.hdf5'
            cls.policy_model = load_keras_models.LocalKerasModel(policy_path)
        if cls.templates == None:
            templates_path = data_folder + '/uspto_templates.hdf5'
            cls.templates = pd.read_hdf(templates_path, "table")


    def _get_saved_actions_if_available(self, smi, first_only=False):
        if smi not in self.actions:
            return None

        if first_only:
            reactions = [self.actions[smi][0]]
        else:
            reactions = self.actions.pop(smi)[1:]
        return reactions

    def get_actions(self, smi, first_only=False):
        saved_actions = self._get_saved_actions_if_available(smi, first_only)
        if saved_actions is not None:
            return saved_actions

        reactions = []
        priors = []
        template_column = self.config.aizynth_template_column

        mol = Chem.MolFromSmiles(smi)

        all_transforms_prop = self._predict(mol)

        probable_transforms_idx = self._cutoff_predictions(all_transforms_prop)

        possible_moves = self.templates.iloc[probable_transforms_idx]
        probs = all_transforms_prop[probable_transforms_idx]

        priors.extend(probs)
        for idx, (move_index, move) in enumerate(possible_moves.iterrows()):
            metadata = dict(move)
            del metadata[template_column]
            metadata["policy_probability"] = round(float(probs[idx]), 5)
            metadata["template_code"] = move_index

            reaction = {'smarts': move[template_column],
                        'metadata': metadata,
                        'prior': priors[idx]}

            reactions.append(reaction)

        if first_only:
            self.actions[smi] = reactions
            reactions = [reactions[0]]

        return reactions

    def get_rxns(self, smile):
        if self.policy_model == None:
            self.load_model()

        reactions = self.get_actions(smile)
        rxns = {}
        metadata = {}

        for reaction in reactions:
            name = f"Chem_{reaction['metadata']['classification']}"
            num = 1
            extra_string = f"__{num}"
            while name+extra_string in rxns:
                extra_string = f"__{num}"
                num += 1
            name = name+extra_string
            rxns[name] = [reaction['smarts']]
            metadata[name] = reaction['metadata']
        return rxns, metadata

    def _predict(self, mol):
        fingerprint = fingerprints.get_mol_fingerprint(mol, 2, nbits=len(self.policy_model))
        fp_arr = fingerprint.reshape([1, len(self.policy_model)])
        return np.array(self.policy_model.predict(fp_arr)).flatten()

    def _cutoff_predictions(self, predictions):
        """
        Get the top transformations, by selecting those that have:
            * cumulative probability less than a threshold (cutoff_cumulative)
            * or at most N (cutoff_number)
        """
        cutoff_cumulative = self.config.aizynth_cutoff_cumulative
        cutoff_number = self.config.aizynth_cutoff_number

        sortidx = np.argsort(predictions)[::-1]
        cumsum = np.cumsum(predictions[sortidx])
        if any(cumsum >= cutoff_cumulative):
            maxidx = np.argmin(cumsum < cutoff_cumulative)
        else:
            maxidx = len(cumsum)

        maxidx = min(maxidx, cutoff_number) or 1
        return sortidx[:maxidx]


if __name__ == '__main__':
    aizynth_action_getter = ActionGetter()
    aizynth_action_filter = ActionFilter(log_level='INFO')
    smi = 'O=C(CC(=O)N1CCn2c(nnc2C(F)(F)F)C1)Cc1cc(F)c(F)cc1F'
    reactions, metadata = aizynth_action_getter.get_rxns(smi)

    rule_applier = RuleApplicator()
    rxn_dict = rule_applier.apply_rdchiral(smi, reactions)

    rxn_dict, metadata = aizynth_action_filter.filter(smi, rxn_dict, metadata)
