import pytest

from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig
from retrobiocat_web.retro.retrosynthesis_engine.mols_and_reactions.starting_material_evaluator import \
    StartingMaterialEvaluator


class Test_SourceMol_Lookup():

    def test_can_lookup_molecule(self):
        sme = StartingMaterialEvaluator()
        smi = 'CCC=O'
        available, info = sme.eval(smi)
        assert available == 1

    def test_random_molecule_not_available(self):
        sme = StartingMaterialEvaluator()
        smi = 'not_a_molecule'
        available, info = sme.eval(smi)
        assert available == 0

    def test_info_dict_is_returned_with_available(self):
        sme = StartingMaterialEvaluator()
        smi = 'CCC=O'
        available, info = sme.eval(smi)
        assert 'askcos' in info

    def test_can_query_metabolism(self):
        sme = StartingMaterialEvaluator()
        smi = 'CCC(=O)C(=O)O'
        available, info = sme.eval(smi, mode='metabolites')
        assert available == 1
        assert 'askcos' not in info


    @pytest.mark.parametrize('allowed,expected', [[True, 1], [False, 0]])
    def test_chiral_molecules_are_found_or_not_when_banned(self, allowed, expected):
        retro_config = RetrosynthesisConfig()
        retro_config.source_mols_can_be_chiral = allowed
        sme = StartingMaterialEvaluator(config=retro_config)
        smi = 'C[C@H](N)c1ccccc1'
        available, info = sme.eval(smi)
        assert available == expected



