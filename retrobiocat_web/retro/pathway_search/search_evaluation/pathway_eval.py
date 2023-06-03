import numpy as np
import time
from retrobiocat_web.retro.network_pathway.network import Network
from retrobiocat_web.retro.pathway_search.pathway_generation.best_first_search import BFS
from retrobiocat_web.retro.pathway_search.pathway_generation.evaluation import PathwayEvaluator
from retrobiocat_web.mongo.models.retro_tests import TestCascade, TestCascadeRun, TestCascadeResult
from retrobiocat_web.mongo.models.reaction_models import Reaction
from retrobiocat_web import __version__
from retrobiocat_web.retro.pathway_search.mcts.mcts import MCTS

default_network_settings = {}

default_pathway_settings = {'max_pathways': 25000,
                            'min_weight': 1,
                            'mode': 'pathway_explorer'}

default_weights = [1, 1, 1, 1, 1]

class TrialSetEvaluator(object):

    def __init__(self, test_run_id,
                 log=False):

        self.test_run_obj = TestCascadeRun.objects(id=test_run_id).first()

        self.network_settings = self.test_run_obj.settings.get('network_settings')
        self.pathway_settings = self.test_run_obj.settings.get('pathway_settings')
        self.weights = self.test_run_obj.settings.get('weights')
        self.max_steps = self.test_run_obj.settings.get('max_steps')

        self.print_log = log

    def create_pathways_using_pathway_explorer(self, target_product):
        network = Network()
        network.update_settings(self.network_settings)
        network.generate(target_product, self.max_steps)
        bfs = BFS(network=network, max_pathways=self.pathway_settings['max_pathways'], min_weight=self.pathway_settings['min_weight'])
        bfs.run()
        pathways = bfs.get_pathways()
        if len(pathways) == self.pathway_settings['max_pathways']:
            self.log(f"  ~max pathways reached")

        #pathways = group_pathways(pathways)
        return pathways

    def create_pathways_using_mcts(self, target_product):
        mcts = MCTS(target_product,
                    max_length=self.max_steps,
                    max_search_time=self.pathway_settings['max_search_time'],
                    exploration=self.pathway_settings['exploration'],
                    max_pathways=self.pathway_settings['max_pathways'])
        mcts.run()
        return mcts.pathways

    def evaluate_pathways(self, pathways):
        pathway_evaluator = PathwayEvaluator(pathways)
        pathway_evaluator.weights = {'Normalised num enzyme steps': self.weights[0],
                                     'Normalised change in complexity': self.weights[1],
                                     'starting_material': self.weights[2],
                                     'postitive_enzyme_steps_score': self.weights[3]}
        pathway_evaluator.diversity_weight = self.weights[4]

        pathway_evaluator.evaluate()

        return pathway_evaluator

    def identify_locations_of_test_substrates(self, df, list_test_substrates, list_test_enzymes):
        list_indexes_abs = []
        list_rows_abs = []
        list_indexes_any = []
        list_rows_any = []

        for index, row in df.iterrows():

            # For absolute ranking
            if self._are_literature_end_nodes_present(list_test_substrates, row['End nodes']):
                enzymes = self._get_reaction_enzymes(row['Pathway'])
                if (all(elem in enzymes for elem in list_test_enzymes)
                    and len(list_test_enzymes) == row['Pathway'].scores.num_enzyme_steps) \
                        or len(list_test_enzymes) == 0:

                    if len(list_indexes_abs) == 0:
                        self.log(f"  ~ absolute ranking {index}, and enzymes match")
                    list_indexes_abs.append(index)
                    list_rows_abs.append(row)
                else:
                    if len(list_indexes_abs) == 0:
                        self.log(f"  ~ absolute ranking {index}, but enzymes do not match")

            # For any ranking
            if self._are_literature_end_nodes_present(list_test_substrates, row['Substrate nodes']):
                enzymes = self._get_reaction_enzymes(row['Pathway'])
                if all(elem in enzymes for elem in list_test_enzymes) or len(list_test_enzymes) == 0:
                    if len(list_indexes_any) == 0:
                        self.log(f"  ~ any ranking {index}, and enzymes match")
                    list_indexes_any.append(index)
                    list_rows_any.append(row)
                else:
                    if len(list_indexes_any) == 0:
                        self.log(f"  ~ any ranking {index}, but enzymes do not match")

        return list_indexes_any, list_rows_any, list_indexes_abs, list_rows_abs

    def run(self, test_id):
        test_obj = TestCascade.objects(id=test_id).first()
        self.log(f"Running test {test_obj.name}")

        t0 = time.time()

        self.log(f" - create pathways for {test_obj.target_smiles} using {self.pathway_settings['mode']}")
        if self.pathway_settings['mode'] == 'pathway_explorer':
            pathways = self.create_pathways_using_pathway_explorer(test_obj.target_smiles)
        elif self.pathway_settings['mode'] == 'mcts':
            pathways = self.create_pathways_using_mcts(test_obj.target_smiles)

        self.log(f" - evaluate pathways for {test_obj.target_smiles}")
        evaluator = self.evaluate_pathways(pathways)

        t1 = time.time()
        t = round(t1 - t0, 0)
        self.log(f'Complete - runtime of {t} seconds.')

        self.log(f" - find ranking of literature pathway")
        lit_end_nodes = test_obj.starting_smiles
        indexes_any, rows_any, indexes_abs, rows_abs = self.identify_locations_of_test_substrates(evaluator.df,
                                                                                                  lit_end_nodes,
                                                                                                  test_obj.enzymes)
        absolute_rank = None
        if len(indexes_abs) != 0:
            absolute_rank = indexes_abs[0] + 1

        any_rank = None
        if len(indexes_any) != 0:
            any_rank = indexes_any[0] + 1

        self.create_new_results_entry(test_id, t, any_rank, absolute_rank, len(pathways))
        self.update_run_scores()

    def log(self, msg):
        if self.print_log == True:
            print(msg)

    def create_new_results_entry(self, test_obj, t, ranking_any, ranking_specific, num_pathways):
        results_entry = TestCascadeResult(test_cascade=test_obj,
                                           time_taken=t,
                                           ranking_any=ranking_any,
                                           ranking_specific=ranking_specific,
                                           total_pathways=num_pathways)
        results_entry.save()

        self.test_run_obj.results.append(results_entry)
        self.test_run_obj.save()

        return results_entry

    def update_run_scores(self):
        test_result_objs = self.test_run_obj.results
        total_tests = len(test_result_objs)

        top_5 = 0
        top_25 = 0
        top_200 = 0
        for result in test_result_objs:
            ranking = result.ranking_any
            if ranking is None:
                ranking = np.inf
            if ranking <= 5:
                top_5 += 1
            if ranking <= 25:
                top_25 += 1
            if ranking <= 200:
                top_200 += 1

        self.test_run_obj.top_5_score = int((top_5/total_tests)*100)
        self.test_run_obj.top_25_score = int((top_25 / total_tests) * 100)
        self.test_run_obj.top_200_score = int((top_200 / total_tests) * 100)
        self.test_run_obj.save()

    def _get_reaction_enzymes(self, pathway):
        enzymes = []
        for reaction in pathway.reactions:
            potential_enzymes = pathway.network.graph.nodes[reaction]['attributes']['possible_enzymes']
            enzymes.extend(potential_enzymes)
        return enzymes

    def _are_literature_end_nodes_present(self, lit_end_nodes, pathway_nodes):
        for node in lit_end_nodes:
            if node not in pathway_nodes:
                return False
        return True

def make_new_test_run(run_name,
                      max_steps=4,
                      weights=None,
                      network_settings=None,
                      pathway_settings=None,
                      overwrite_ok=False):

    if TestCascadeRun.objects(run_name=run_name).count() != 0:
        if overwrite_ok == True:
            print('deleting old test run with matching run name')
            TestCascadeRun.objects(run_name=run_name).first().delete()
        else:
            print('ERROR - test run with this name already exists')
            return None

    if network_settings is None:
        network_settings = default_network_settings

    if pathway_settings is None:
        pathway_settings = default_pathway_settings

    if weights is None:
        weights = default_weights

    max_steps = max_steps
    reaction_rules_last_updated = Reaction.objects().order_by('+last_updated')[0].last_updated

    settings = {'network_settings': network_settings,
                'pathway_settings': pathway_settings,
                'weights': weights,
                'max_steps': max_steps}

    test_run = TestCascadeRun(run_name=run_name,
                              settings=settings,
                              reaction_rules_last_updated=reaction_rules_last_updated,
                              retrobiocat_version=__version__)
    test_run.save()

    return test_run

def load_test_ids(specific_type=None):
    if specific_type is None:
        test_cascades_ids = TestCascade.objects().distinct('id')
    else:
        test_cascades_ids = TestCascade.objects(test_type=specific_type).distinct('id')
    return test_cascades_ids


def run_test(test_id, test_run_id, log=True):
    test_set_evaluator = TestSetEvaluator(test_run_id, log=log)
    test_set_evaluator.run(test_id)


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

   # test_run = make_new_test_run('test_run_1_no_literature', weights=[1, 1, 1, 0, 1], overwrite_ok=True)
    pathway_settings = {'max_search_time': 60, 'exploration': 100, 'max_pathways': 25000, 'mode': 'mcts'}
    test_run = make_new_test_run('test_run_mcts_1_60s',
                                 pathway_settings=pathway_settings,
                                 max_steps=5,
                                 weights=[1, 1, 1, 1, 1],
                                 overwrite_ok=True)

    test_ids = load_test_ids(specific_type='biocatalysis')


    for test_id in test_ids:
        run_test(test_id, test_run.id)