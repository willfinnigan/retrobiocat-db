from retrobiocat_web.retro.network_pathway.rdkit_utils import rdkit_smile
from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig
from retrobiocat_web.logging import add_logger
from retrobiocat_web.retro.retrosynthesis_engine.reaction_processing.reaction_parser import ReactionParser
from retrobiocat_web.retro.retrosynthesis_engine.reaction_processing.reaction_ranker import ReactionRanker
from retrobiocat_web.retro.retrosynthesis_engine.rule_application.rule_applicator import RuleApplicator


class CustomRxn_Expander():

    def __init__(self, network, config=None, log_level='WARNING'):

        self.logger = add_logger('CustomRxn_Expander', level=log_level)

        self.config = config
        if self.config is None:
            self.config = RetrosynthesisConfig()

        self.rule_applicator = RuleApplicator(config=config, log_level=log_level)
        self.reaction_parser = ReactionParser(network, config=config, log_level=log_level)
        self.reaction_ranker = ReactionRanker(config=config, log_level=log_level)


    def expand(self, target_smi, substrate_smis, rxn_name):

        target_smi = rdkit_smile(target_smi)
        substrate_smis = [rdkit_smile(smi) for smi in substrate_smis]
        reactions = self.reaction_parser.parse(target_smi, {rxn_name: [substrate_smis]}, rxn_type='custom')

        if len(reactions) > 0:
            self.logger.info(f"Created custom reaction {target_smi}<--{substrate_smis}")
        else:
            self.logger.warning(f"Custom reaction {target_smi}<--{substrate_smis} was removed during parsing")

        return reactions

    def expand_custom_rules(self, target_smi, smarts, multi_step_smarts=None):
        smi = rdkit_smile(target_smi)

        rxns = self._smarts_to_rxns(smarts)
        multi_step_rxns = self._multi_smarts_to_multi_rxns(multi_step_smarts)
        result = self._apply_custom_rules(smi, rxns, multi_step_rxns)

        reactions = self.reaction_parser.parse(smi, result, rxn_type='custom')
        return reactions

    def _apply_custom_rules(self, smi, rxns, multi_step_rxns):
        if self.config.custom_expander_rule_mode == 'rdchiral':
            result = self.rule_applicator.apply_rdchiral(smi, {'custom_rxn': rxns}, multistep_rxns={'custom_multi_rxn': multi_step_rxns})
        elif self.config.custom_expander_rule_mode == 'rdkit':
            result = self.rule_applicator.apply_rdkit(smi, {'custom_rxn': rxns}, multistep_rxns={'custom_multi_rxn': multi_step_rxns})
        else:
            self.logger.error(f"Rule mode for custom expander is invalid ({self.config.custom_expander_rule_mode})")
            Exception(f"Rule mode for custom expander is invalid ({self.config.custom_expander_rule_mode})")
            result = []
        return result

    def _smarts_to_rxns(self, smarts):
        if self.config.custom_expander_rule_mode == 'rdchiral':
            rxns = self.rule_applicator.smarts_to_rdchiral(smarts)
        elif self.config.custom_expander_rule_mode == 'rdkit':
            rxns = self.rule_applicator.smarts_to_rdkit(smarts)
        else:
            self.logger.error(f"Rule mode for custom expander is invalid ({self.config.custom_expander_rule_mode})")
            Exception(f"Rule mode for custom expander is invalid ({self.config.custom_expander_rule_mode})")
            rxns = []
        return rxns

    def _multi_smarts_to_multi_rxns(self, multi_smarts):
        multi_rxns = []

        for group_smas in multi_smarts:
            group_rxns = []
            for step_smas in group_smas:
                if self.config.custom_expander_rule_mode == 'rdchiral':
                    step_rxns = self.rule_applicator.smarts_to_rdchiral(step_smas)
                elif self.config.custom_expander_rule_mode == 'rdkit':
                    step_rxns = self.rule_applicator.smarts_to_rdkit(step_smas)
                else:
                    self.logger.error(f"Rule mode for custom expander is invalid ({self.config.custom_expander_rule_mode})")
                    Exception(f"Rule mode for custom expander is invalid ({self.config.custom_expander_rule_mode})")
                    step_rxns = []
                group_rxns.append(step_rxns)
            multi_rxns.append(group_rxns)
        return multi_rxns


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()
