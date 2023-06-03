import pandas as pd

class PathwayEvaluator():

    def __init__(self, list_pathways):
        self.df = None
        self.pathways = list_pathways
        self.weights = {'Normalised num enzyme steps': 1,
                        'Normalised change in complexity': 1,
                        'starting_material': 1,
                        'postitive_enzyme_steps_score': 1}

        self.diversity_weight = 1


    def evaluate(self):
        self.make_evalation_df()
        self.calculate_normalised_scores()
        self.calc_weighted_scores()
        self.df.sort_values('Weighted Score', ascending=False, inplace=True)
        self.df = self.df.reset_index(drop=True)

    def make_evalation_df(self):

        to_df = {'Pathway': [],
                 'Substrate nodes': [],
                 'Reaction nodes': [],
                 'End nodes': []}

        for pathway in self.pathways:
            to_df['Pathway'].append(pathway)
            to_df['Substrate nodes'].append(pathway.substrates)
            to_df['Reaction nodes'].append(pathway.reaction_names)
            to_df['End nodes'].append(pathway.end_nodes)

            scores_dict = pathway.scores.scores_dict()
            for name in scores_dict:
                if name not in to_df:
                    to_df[name] = []
                to_df[name].append(scores_dict[name])

        self.df = pd.DataFrame(to_df)
        return self.df

    def take_best_pathways(self, number):
        self.df = self.df.head(number)

    def calculate_normalised_scores(self):
        def normalise_complexity(df):
            correct_neg_change_in_complexity = 0
            while (df['change_in_complexity'].max() + correct_neg_change_in_complexity) <= 0:
                correct_neg_change_in_complexity += 0.05
            df['Normalised change in complexity'] = df['change_in_complexity'] / (df['change_in_complexity'].max() + correct_neg_change_in_complexity)

            return df

        def normalise_num_steps(df):
            df['Rescaled num enzyme steps'] = df['num_enzyme_steps'].max() - df['num_enzyme_steps']
            df['Normalised num enzyme steps'] = df['Rescaled num enzyme steps'] / df['Rescaled num enzyme steps'].max()
            return df


        self.df = normalise_complexity(self.df)
        self.df = normalise_num_steps(self.df)

    def calc_weighted_scores(self):

        scores_to_df = []

        for index, row in self.df.iterrows():
            score = 0
            for name in self.weights:
                score += row[name] * self.weights[name]
            scores_to_df.append(score)

        self.df['Weighted Score'] = scores_to_df

        if self.diversity_weight != 0:
            self.rank_diversity()

        return self.df

    def apply_diversity(self, top=25, reaction_only=False):

        diversity = []
        div_dict = {}  # a dict of node names and how many times they have appeared

        for index, row in self.df.iterrows():  # loop through the df

            if index <= top:   # we're only going to adjust the diversity in the top 25 or so
                pathway = row['Pathway']
                pathway_diversity = 0   # a penalty to be applied based on how many identical reactions / substrates

                for node in pathway.list_nodes:  # loop through nodes looking for ones which have come up before
                    if (reaction_only == False) or pathway.network.graph.nodes[node]['attributes']['node_type'] == 'reaction':
                        node_name = pathway.network.graph.nodes[node]['attributes']['name']
                        if node_name not in div_dict:
                            div_dict[node_name] = 0  # store number of times a node appears
                        pathway_diversity += div_dict[node_name]  # if node has appeared before, add to diversity penalty
                        div_dict[node_name] += 1
                diversity.append(pathway_diversity)
            else:
                diversity.append(0)  # pathways beyond the top 25 or so dont have a penalty applied

        self.df['Diversity'] = diversity
        self.df['Diversity'] = self.df['Diversity'] / self.df['Diversity'].max()

        self.df['Diversity Weighted Score'] = (self.df['Diversity']*self.diversity_weight) + self.df['Weighted Score']


    def rank_diversity(self, n=10):

        for i in range(n):
            self.apply_diversity()
            self.df.sort_values('Diversity Weighted Score', ascending=False, inplace=True)