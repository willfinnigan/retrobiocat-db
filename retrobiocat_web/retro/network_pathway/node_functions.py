def sort_by_score(list_nodes, list_scores, reverse=False):
    scores_nodes = list(zip(list_scores, list_nodes))
    scores_nodes.sort(reverse=reverse)
    nodes_sorted = [node for score, node in scores_nodes]
    return nodes_sorted

