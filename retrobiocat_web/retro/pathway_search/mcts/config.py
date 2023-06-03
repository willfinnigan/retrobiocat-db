import numpy as np

class MCTS_Config():

    def __init__(self):

        self.biosynthesis_expansion = False
        self.chemistry_expansion = True
        self.biocatalysis_expansion = True
        self.maxLength = 4
        self.exploration = 2
        self.max_search_time = 120
        self.max_pathways = np.inf
        self.stop_expansion_if_nonbuyable_at_max_length = True

        # mcts scoring
        self.adjust_non_buyable_penalties_with_complexity = True
        self.min_penalty_for_non_buyable = -0.1
        self.complexity_for_minus_1 = 0  # complexity above this will always score -1
        self.complexity_for_zero = -1  # complexity below this will always score min penalty

        # pathways
        self.only_solved = True

        # clustering
        self.dont_cluster_single_step_routes = True
        self.distance_threshold = 2

        # costing
        self.in_stock_cost = 1
        self.not_in_stock_cost = 10
        self.reaction_cost = 1
        self.reaction_yield = 0.8


    def update_from_dict(self, attr_dict):
        current_dict = self.to_dict()
        for key, value in attr_dict.items():
            if key in current_dict:
                setattr(self, key, value)
        return self

    def to_dict(self):
        return self.__dict__

