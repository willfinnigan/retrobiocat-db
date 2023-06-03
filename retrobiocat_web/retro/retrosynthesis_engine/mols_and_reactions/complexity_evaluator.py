from retrobiocat_web.retro.retrosynthesis_engine.mols_and_reactions.scscore.standalone_model_numpy import \
    sc_scorer

class Complexity_Evaluator():

    def get_complexity(self, smi, round_to=4):
        return get_complexity(smi, round_to)

def get_complexity(smi, round_to=4):
    (smile, score) = sc_scorer.get_score_from_smi(smi)
    return round(score, round_to)

if __name__ == '__main__':
    smis = ['CCCOCCC', 'CCCNc1ccccc1', 'CCC=O']
    eval = Complexity_Evaluator()
    for smi in smis:
        complexity = eval.get_complexity(smi)
        print(complexity)