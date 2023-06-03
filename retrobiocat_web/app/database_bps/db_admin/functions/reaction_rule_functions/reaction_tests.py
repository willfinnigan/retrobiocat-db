import yaml
from retrobiocat_web.retro.retrosynthesis_engine.rule_application.rule_applicator import RuleApplicator

from retrobiocat_web.retro.retrosynthesis_engine.rule_application.rdchiral.initialization import rdchiralReaction

from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType
from retrobiocat_web.mongo.models.reaction_models import Reaction
from retrobiocat_web.retro.network_pathway.rdkit_utils import rdkit_smile
from retrobiocat_web.app.retrobiocat.functions.get_images import smiles_rxn_to_svg
import copy


class ReactionTester(object):

    def __init__(self):
        self.state = 'danger'
        self.msgs = {'danger': 'Tests failed',
                     'success': 'Tests passed',
                     'warning': 'Tests passed with warnings'}
        self.issues = []

    def get_msg(self):
        return self.msgs[self.state]

    def run(self, selection, name, smarts, cofactors, positive_tests, negative_tests,
            type, experimental, two_step, requires_absence_of_water, example_rxn_smiles, steps_yaml,
            product_seeds, substrate_1_seeds, substrate_2_seeds):

        self.state = 'success'

        self._test_name(name, selection)
        self._test_cofactors_dict(cofactors)
        self._test_example_rxn_smiles(example_rxn_smiles)
        self._test_seeds_are_mols(product_seeds, substrate_1_seeds, substrate_2_seeds)
        list_rxns = self._test_smarts(smarts)

        multi_step_rxns = self._test_steps(steps_yaml)

        if list_rxns is not None:
            self._positive_tests(positive_tests, list_rxns, multi_step_rxns)
            self._negative_tests(negative_tests, list_rxns, multi_step_rxns)


    def _test_name(self, rxn_name, selection):
        if rxn_name == '' or rxn_name == ' ' or rxn_name == 'Empty template':
            self.issues.append('Reaction must have a name')
            self.state = 'danger'

        if rxn_name != selection:
            if len(Reaction.objects(name=rxn_name)) != 0:
                self.issues.append('Reaction name already exists')
                self.state = 'danger'

    def _test_cofactors_dict(self, cofactors):

        try:
            cofactors_dict = yaml.load(cofactors, Loader=yaml.FullLoader)
        except:
            self.state = 'danger'
            self.issues.append('Could not load cofactors yaml')
            return

        if cofactors_dict is None:
            self.state = 'danger'
            self.issues.append('At least one entry must be made for enzymes/cofactors')
            return

        for enz in list(cofactors_dict.keys()):
            if enz not in list(EnzymeType.objects().distinct('enzyme_type')):
                self.state = 'danger'
                self.issues.append(f'Enzyme type {enz} is not defined')

            elif type(cofactors_dict[enz]) != dict:
                self.state = 'danger'
                self.issues.append(f'Cofactors for {enz} are not structured as a dictionary')

            elif 'cofactors_plus' not in cofactors_dict[enz] or 'cofactors_minus' not in cofactors_dict[enz]:
                self.state = 'danger'
                self.issues.append(f'Cofactors must be defined')

            elif type(cofactors_dict[enz]['cofactors_plus']) != list or type(cofactors_dict[enz]['cofactors_minus']) != list:
                self.state = 'danger'
                self.issues.append(f'Cofactors must be given as lists')

    def _test_smarts(self, smarts):
        try:
            smarts_list = yaml.load(smarts, Loader=yaml.FullLoader)
        except:
            self.state = 'danger'
            self.issues.append('Could not load smarts yaml')
            return

        try:
            rxn_list = []
            for sma in smarts_list:
                rxn_list.append(rdchiralReaction(sma))
        except:
            self.state = 'danger'
            self.issues.append('Could not load reactions')
            return

        return rxn_list

    def _test_steps(self, steps):

        multi_step_rxns = []
        try:
            groups_steps_list = yaml.load(steps, Loader=yaml.FullLoader)
        except:
            self.state = 'danger'
            self.issues.append('Could not load steps yaml')
            return

        try:
            for steps_list in groups_steps_list:
                grouped_smas = []
                for name in steps_list:
                    rxn_q = Reaction.objects(name=name)
                    if len(rxn_q) == 0:
                        self.state = 'danger'
                        self.issues.append(f'Reaction: "{name}" does not exist')
                    elif len(rxn_q) > 1:
                        self.state = 'danger'
                        self.issues.append(f'Error - multiple reactions named "{name}"')

                    step_rules = [rdchiralReaction(sma) for sma in rxn_q[0].smarts]
                    grouped_smas.append(step_rules)
                multi_step_rxns.append(grouped_smas)

        except:
            self.state = 'danger'
            self.issues.append('Could not load reactions')
            return

        return multi_step_rxns

    def _positive_tests(self, positive_tests, list_rxns, list_multi_step_rxns):
        rule_applicator = RuleApplicator()
        rxns = {'tests': list_rxns}
        multi_step_rxns = {'tests': list_multi_step_rxns}

        try:
            positive_tests = yaml.load(positive_tests, Loader=yaml.FullLoader)
        except:
            self.state = 'danger'
            self.issues.append('Could not load positive tests yaml')
            return

        for test_product in positive_tests:
            try:
                rdkit_smile(test_product)
            except:
                self.state = 'danger'
                self.issues.append(f'Positive test SMILE: {test_product} not accepted by rdkit')
                return

        try:
            for test_product in positive_tests:
                reaction_outcomes = self._apply_reactions(rule_applicator, test_product, rxns, multi_step_rxns)
                if len(reaction_outcomes) == 0:
                    self.state = 'danger'
                    self.issues.append(f'Reaction not in outcomes for tested positive product: {test_product}')
        except Exception as e:
            self.state = 'danger'
            self.issues.append('Problem running positive tests')
            self.issues.append(str(e))
            return

        return True

    def _negative_tests(self, negative_tests, list_rxns, list_multi_step_rxns):
        rule_applicator = RuleApplicator()
        rxns = {'tests': list_rxns}
        multi_step_rxns = {'tests': list_multi_step_rxns}

        try:
            negative_tests = yaml.load(negative_tests, Loader=yaml.FullLoader)
        except:
            self.state = 'danger'
            self.issues.append('Could not load negative tests yaml')
            return

        for test_product in negative_tests:
            try:
                rdkit_smile(test_product)
            except:
                self.state = 'danger'
                self.issues.append(f'Negative test SMILE: {test_product} not accepted by rdkit')
                return

        for test_product in negative_tests:
            reaction_outcomes = self._apply_reactions(rule_applicator, test_product, rxns, multi_step_rxns)
            if len(reaction_outcomes) != 0:
                self.state = 'danger'
                self.issues.append(f'Reaction should not be outcomes for tested negative product: {test_product}')

        try:
            for test_product in negative_tests:
                reaction_outcomes = self._apply_reactions(rule_applicator, test_product, rxns, multi_step_rxns)
                if len(reaction_outcomes) != 0:
                    self.state = 'danger'
                    self.issues.append(f'Reaction should not be outcomes for tested negative product: {test_product}')
        except:
            self.state = 'danger'
            self.issues.append('Problem running negative tests')
            return

        return True

    def _apply_reactions(self, rule_applicator: RuleApplicator, test_product: str, rxns, multi_step_rxns):
        reaction_outcomes = []

        rule_applicator.config.combine_enantiomers = True
        reaction_outcomes += rule_applicator.apply_rdchiral(test_product, rxns, multistep_rxns=multi_step_rxns)

        rule_applicator.config.combine_enantiomers = False
        reaction_outcomes += rule_applicator.apply_rdchiral(test_product, rxns, multistep_rxns=multi_step_rxns)

        return reaction_outcomes

    def _test_example_rxn_smiles(self, example_rxn_smiles):
        if example_rxn_smiles == "":
            return

        try:
            smiles_rxn_to_svg(example_rxn_smiles)
        except:
            self.issues.append('Example rxn smiles incorrect')
            self.state = 'danger'

    def _test_seeds_are_mols(self, product_seeds, substrate_1_seeds, substrate_2_seeds):
        try:
            product_seeds = yaml.load(product_seeds, Loader=yaml.FullLoader)
            substrate_1_seeds = yaml.load(substrate_1_seeds, Loader=yaml.FullLoader)
            substrate_2_seeds = yaml.load(substrate_2_seeds, Loader=yaml.FullLoader)

        except:
            self.issues.append(f'Problem loading seed smiles')
            self.state = 'danger'
            return

        print(product_seeds)
        print(substrate_1_seeds)
        print(substrate_2_seeds)

        if product_seeds is not None:
            for smi in product_seeds:
                if rdkit_smile(smi) is None:
                    self.issues.append(f'Example product seed {smi} is not accepted by RdKit')
                    self.state = 'danger'

        if substrate_1_seeds is not None:
            for smi in substrate_1_seeds:
                if rdkit_smile(smi) is None:
                    self.issues.append(f'Example substrate 1 seed {smi} is not accepted by RdKit')
                    self.state = 'danger'

        if substrate_2_seeds is not None:
            for smi in substrate_2_seeds:
                if rdkit_smile(smi) is None:
                    self.issues.append(f'Example substrate 2 seed {smi} is not accepted by RdKit')
                    self.state = 'danger'



