from retrobiocat_web.retro.pathway_search.pathway_generation.group_pathways import group_pathways
from retrobiocat_web.retro.pathway_search.pathway_generation.evaluation import PathwayEvaluator


def run_evaluate_pathways(pathways, weights):
    pathways = group_pathways(pathways)
    pathway_evaluator = PathwayEvaluator(pathways)
    pathway_evaluator.weights = {'Normalised num enzyme steps': weights[0],
                                 'Normalised change in complexity': weights[1],
                                 'starting_material': weights[2],
                                 'postitive_enzyme_steps_score': weights[3]}
    pathway_evaluator.diversity_weight = weights[4]
    pathway_evaluator.evaluate()
    return pathway_evaluator