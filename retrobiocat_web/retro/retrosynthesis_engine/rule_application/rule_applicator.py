from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig
from retrobiocat_web.retro.retrosynthesis_engine.rule_application.components.rdchiral_applicator import RdChiral_Applicator
from retrobiocat_web.retro.retrosynthesis_engine.rule_application.components.rdkit_applicator import RdKit_Applicator
from retrobiocat_web.retro.retrosynthesis_engine.rule_application.components.multi_step_applicator import Multi_step_applicator
from retrobiocat_web import logging

class RuleApplicator():

    def __init__(self, config=None, log_level='WARNING'):
        self.config = config
        if self.config is None:
            self.config = RetrosynthesisConfig()

        self.rdkit_applicator = RdKit_Applicator(config=self.config, log_level=log_level)
        self.rdchiral_applicator = RdChiral_Applicator(config=self.config, log_level=log_level)
        self.multi_step_applicator = Multi_step_applicator(config=self.config, log_level=log_level)
        self.logger = logging.add_logger('RuleApplicator', level=log_level)

    def apply_rdkit(self, smi, rxns, multistep_rxns=None):

        product_dict = self.rdkit_applicator.apply_rules(smi, rxns)

        if multistep_rxns is not None:
            multi_product_dict = self.multi_step_applicator.apply_multi_step_rules(self.rdkit_applicator, smi, multistep_rxns)
            product_dict.update(multi_product_dict)

        return product_dict

    def apply_rdchiral(self, smi, rxns, multistep_rxns=None):

        product_dict = self.rdchiral_applicator.apply_rules(smi, rxns)

        if multistep_rxns is not None:
            multi_product_dict = self.multi_step_applicator.apply_multi_step_rules(self.rdchiral_applicator, smi,
                                                                                   multistep_rxns)
            product_dict.update(multi_product_dict)

        return product_dict

    def smarts_to_rdchiral(self, smarts, remove_incorrect_atom_numbering=False):
        return self.rdchiral_applicator.get_rxns(smarts, remove_incorrect_atom_numbering=remove_incorrect_atom_numbering)

    def smarts_to_rdkit(self, smarts):
        return self.rdkit_applicator.get_rxns(smarts)

if __name__ == '__main__':

    applier = RuleApplicator(log_level='DEBUG')

    rxn1 = '[#6X4;z1:2]-[#7X3;z0:3]>>[#6X3;z1:2]=[#7X2;z0:3]'
    rxn2 = '[#6:1]-[#6@@H;X4:3](-[#6:2])-[#7X3;z0:4]>>[#6:1]-[#6H0;X3:3](-[#6:2])=[#7X2;z0:4]'
    rxn3 = '[#6:1]-[#6@H;X4:3](-[#6:2])-[#7X3;z0:4]>>[#6:1]-[#6H0;X3:3](-[#6:2])=[#7X2;z0:4]'
    rxns = {'Imine reduction': applier.smarts_to_rdchiral([rxn1, rxn2, rxn3])}

    products = applier.apply_rdchiral('[C@H]1(C2=CC=CC=C2)NCCCC1', rxns)

