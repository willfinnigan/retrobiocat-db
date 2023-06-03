import numpy as np


class PathwayScores():

    def __init__(self, pathway, calc_scores=True, multiple_end_substrates_policy='max'):
        self.pathway = pathway
        self.multiple_end_substrates_policy = multiple_end_substrates_policy

        # scores
        self.change_in_complexity = 0
        self.starting_material = 0  # 1 = all buyable
        self.num_enzyme_steps = 0
        self.positive_enzyme_steps_score = 0

        self.simple_score = 0

        if calc_scores==True:
            self.score()

    def scores_dict(self):
        scores_dict = {'change_in_complexity' : self.change_in_complexity,
                       'starting_material' : self.starting_material,
                       'num_enzyme_steps' : self.num_enzyme_steps,
                       'postitive_enzyme_steps_score': self.positive_enzyme_steps_score,
                       'simple_score': self.simple_score}

        return scores_dict

    def scores_list(self):
        scores_list = []
        scores_dict = self.scores_dict()
        for name in scores_dict:
            scores_list.append(scores_dict[name])
        return scores_list

    def scores_from_dict(self, scores_dict):
        self.change_in_complexity = scores_dict['change_in_complexity']
        self.starting_material = scores_dict['starting_material']
        self.num_enzyme_steps = scores_dict['num_enzyme_steps']
        self.positive_enzyme_steps_score = scores_dict['postitive_enzyme_steps_score']

    def score(self):
        self.starting_material = self.score_starting_material(self.pathway)
        self.num_intermediates = len(self.pathway.substrates) - 2
        self.num_enzyme_steps = self.score_num_enzyme_steps(self.pathway)

        self.change_in_complexity = self.score_complexity_change(self.pathway)
        if self.num_enzyme_steps != 0:
            self.complexity_per_enzyme = self.change_in_complexity / self.num_enzyme_steps
        if self.num_intermediates != 0:
            self.complexity_per_intermediate = self.change_in_complexity / self.num_intermediates

        self.positive_enzyme_steps_score = self.score_positive_enzyme_steps(self.pathway, self.num_enzyme_steps)

        self.simple_score = self.calc_simple_score()

    def calc_simple_score(self):
        complexity = self.change_in_complexity
        steps = self.num_enzyme_steps
        enz_steps = self.positive_enzyme_steps_score

        if steps == 0:
            steps = 1

        simple_score = ((complexity/steps) + (complexity*enz_steps)) / 2
        return simple_score

    def score_complexity_change(self, pathway):

        graph = pathway.network.graph
        target_smiles = pathway.network.target_smiles
        end_nodes = self.pathway.end_nodes
        end_complexity = []


        if type(pathway.end_nodes) != list:
            end_nodes = [end_nodes]
            print('Error end nodes are not a list - ' + str(end_nodes))

        for node in pathway.end_nodes:
            try:
                end_complexity.append(graph.nodes[node]['attributes']['complexity'])
            except:
                print('Error node not in graph - ' + str(node))
                end_complexity.append(0)

        try:
            if self.multiple_end_substrates_policy == 'max':
                end = np.max(end_complexity)
            elif self.multiple_end_substrates_policy == 'min':
                end = np.min(end_complexity)
            elif self.multiple_end_substrates_policy == 'mean':
                end = np.mean(end_complexity)
            else:
                end = np.mean(end_complexity)
        except:
            print('Error calculating avg complexity change from end nodes ' + str(end_nodes))
            end = 0

        start = graph.nodes[target_smiles]['attributes']['complexity']

        return start - end

    def score_starting_material(self, pathway):

        # average of scores
        scores = []
        for node in pathway.end_nodes:
            if pathway.network.graph.nodes[node]['attributes']['node_type'] == 'reaction':
                raise Exception('End node is a reaction')
            else:
                scores.append(pathway.network.graph.nodes[node]['attributes']['is_starting_material'])

        if len(scores) != 0:
            score = sum(scores) / len(scores)
        else:
            score = 0

        return score

    def score_num_enzyme_steps(self, pathway):

        enzyme_steps = 0
        for node in pathway.reactions:
            if pathway.network.graph.nodes[node]['attributes']['reaction_type'] == 'retrorules':
                enzyme_steps += 1
            elif pathway.network.graph.nodes[node]['attributes']['reaction_type'] == 'retrobiocat':
                enzyme_steps += pathway.network.graph.nodes[node]['attributes']['metadata']['is_enzyme']

        return enzyme_steps

    def score_positive_enzyme_steps(self, pathway, num_enzyme_steps):
        total_score = 0
        for node in pathway.reactions:
            enz = pathway.reaction_enzyme_map[node]
            if pathway.network.graph.nodes[node]['attributes']['reaction_type'] == 'retrobiocat':
                if pathway.network.graph.nodes[node]['attributes']['metadata']['is_enzyme'] == 1:
                    if 'specificity_scores' in pathway.network.graph.nodes[node]['attributes']['metadata']:
                        total_score += pathway.network.graph.nodes[node]['attributes']['metadata']['specificity_scores'][enz]

        if num_enzyme_steps == 0:
            positive_enzyme_steps_score = 0
        else:
            positive_enzyme_steps_score = total_score / num_enzyme_steps

        return round(positive_enzyme_steps_score,2)

