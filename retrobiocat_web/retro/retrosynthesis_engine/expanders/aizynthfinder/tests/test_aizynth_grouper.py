from retrobiocat_web.retro.retrosynthesis_engine.expanders.aizynthfinder.aizynthfinder_grouper import AIZynth_Grouper
from retrobiocat_web.retro.retrosynthesis_engine.mols_and_reactions.molecules_and_reactions import Reaction, Molecule
import pytest

@pytest.fixture()
def reactions():
    target_mol = Molecule('CCO', 1)
    substrate_mol = Molecule('CC=O', 1)

    rxn_1 = Reaction('test_rxn__1', target_mol, [substrate_mol],
                     metadata={'policy_probability': 0.5, 'library_occurence': 100},
                     rxn_type='aizynthfinder')

    rxn_2 = Reaction('test_rxn__2', target_mol, [substrate_mol],
                     metadata={'policy_probability': 0.1, 'library_occurence': 50},
                     rxn_type='aizynthfinder')

    return [rxn_1, rxn_2]

class Test_AIZynth_Grouper():

    def test_reactions_can_be_grouped_by_name(self, reactions):
        grouper = AIZynth_Grouper()
        groups = grouper.group(reactions)
        assert list(groups.keys()) == ['test_rxn']

    def test_reactions_are_grouped_by_substrate(self, reactions):
        grouper = AIZynth_Grouper()
        groups = grouper.group(reactions)
        assert len(groups['test_rxn']) == 1
        assert len(groups['test_rxn'][0]) == 2

    def test_can_get_a_representative_reaction(self, reactions):
        grouper = AIZynth_Grouper()
        rep_reactions = grouper.group_and_select_representatives(reactions)
        assert rep_reactions[0].metadata['policy_probability'] == 0.6



