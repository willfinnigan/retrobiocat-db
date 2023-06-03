import retrobiocat_web.retro.network_pathway.graph_functions

""" 
Possible scores are:
    'change_in_complexity' 
    'starting_material'
    'num_intermediates'
    'num_enzyme_steps'
    'complexity_per_enzyme'
    'complexity_per_intermediate'
    'postitive_enzyme_steps_score'
"""

default_scores_to_use = ['change_in_complexity',
                         'starting_material',
                         'postitive_enzyme_steps_score']

def _get_scores_list(pathway, scores_to_use):
    scores_list = []
    scores_dict = pathway.scores.scores_dict()
    for name in scores_dict:
        if name in scores_to_use:
            scores_list.append(scores_dict[name])
    return scores_list

def _generate_end_nodes_dict(list_pathways, scores_to_use):
    """
    Iterates over pathways, sorting pathways into groups of: [end_nodes][reactions][scores]

    Pathways which are identical for these three are grouped together.
    """

    end_nodes_dict = {}
    for pathway in list_pathways:
        end_nodes = str(sorted(pathway.end_nodes))
        if end_nodes not in end_nodes_dict:
            end_nodes_dict[end_nodes] = {}

        reactions = retrobiocat_web.retro.network_pathway.graph_functions.get_reaction_names(pathway.reactions, pathway.network.graph)
        reactions = str(sorted(reactions))
        if reactions not in end_nodes_dict[end_nodes]:
            end_nodes_dict[end_nodes][reactions] = {}

        scores = _get_scores_list(pathway, scores_to_use)
        scores = str(scores)
        if scores not in end_nodes_dict[end_nodes][reactions]:
            end_nodes_dict[end_nodes][reactions][scores] = []
        end_nodes_dict[end_nodes][reactions][scores].append(pathway)

    return end_nodes_dict

def _get_grouped_pathways(end_nodes_dict):
    """ From end_nodes_dict, makes a list of grouped_pathway lists"""
    grouped_pathways = []
    for end_nodes in end_nodes_dict:
        for reactions in end_nodes_dict[end_nodes]:
            for scores in end_nodes_dict[end_nodes][reactions]:
                grouped_pathways.append(end_nodes_dict[end_nodes][reactions][scores])

    return grouped_pathways

def _collapse_groups(groups):
    """
    From a group of pathways, choose one to represent the rest
    """

    def sort_group(group):
        scores = []
        pathway_indexes = []
        for i, pathway in enumerate(group):
            enzyme_score = float(pathway.scores.positive_enzyme_steps_score)
            scores.append(enzyme_score)
            pathway_indexes.append(i)

        zipped_pairs = zip(scores, pathway_indexes)

        sorted_indexes = [x for _, x in sorted(zipped_pairs)]

        sorted_group = []
        for index in sorted_indexes:
            sorted_group.append(group[index])

        return sorted_group

    pathways = []

    for group in groups:
        sorted_group = sort_group(group)
        pathway = sorted_group.pop(0)
        for other_pathway in sorted_group:
            pathway.other_varients.append(other_pathway)
            pathway.other_varients_as_nodes.append(other_pathway.list_nodes)
        pathways.append(pathway)

    return pathways

def group_pathways(pathways, scores_to_use=None):
    if scores_to_use == None:
        scores_to_use = default_scores_to_use

    end_nodes_dict = _generate_end_nodes_dict(pathways, scores_to_use)  # groups by end_nodes, reactions, scores
    grouped_pathways = _get_grouped_pathways(end_nodes_dict)  # converts dict to a list of lists
    new_pathways = _collapse_groups(grouped_pathways)
    return new_pathways


if __name__ == '__main__':
    from retrobiocat_web.retro.network_pathway.network import Network
    from retrobiocat_web.retro.pathway_search.pathway_generation.best_first_search import BFS

    network = Network(max_nodes=300)
    network.generate('[C@H]1(C2=CC=CC=C2)NCCCC1', 5)

    bfs = BFS(network=network, max_pathways=200)
    bfs.run()
    pathways = bfs.get_pathways()

    pathways = group_pathways(pathways)
    pathways = pathways[0:10]