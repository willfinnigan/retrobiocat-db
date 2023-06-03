from retrobiocat_web.retro.retrosynthesis_engine.mols_and_reactions.complexity_evaluator import \
    Complexity_Evaluator


class Test_Complexity_Evaluator():

    def test_can_get_complexity(self):
        evaluator = Complexity_Evaluator()
        smi = 'CCC=O'
        score = evaluator.get_complexity(smi)
        assert score == 1.2394