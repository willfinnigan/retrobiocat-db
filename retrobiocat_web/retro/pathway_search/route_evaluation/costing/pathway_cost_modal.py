import retrobiocat_web.retro.network_pathway.graph_functions
from retrobiocat_web.logging import add_logger

class Pathway_Coster():

    def __init__(self, instock_cost=1, not_instock_cost=10, reaction_cost=1, reaction_yield=0.8, name='default', log_level='WARNING'):
        self.in_stock_cost = instock_cost
        self.not_in_stock_cost = not_instock_cost
        self.reaction_cost = reaction_cost
        self.reaction_yield = reaction_yield
        self.name = name

        self.logger = add_logger(f'Pathway coster ({name})', level=log_level)

    def cost(self, pathway):
        self._add_costs_to_end_nodes(pathway)
        self._iteratively_cost_other_nodes(pathway)
        cost = pathway.network.graph.nodes[pathway.target_smiles]['attributes']['cost']
        self.logger.debug(f"Cost of pathway is {cost}")
        pathway.costs[self.name] = cost
        return cost

    def order_by_cost(self, pathways):
        costs = [self.cost(p) for p in pathways]
        costs_pathways = list(zip(costs, pathways))
        sorted_costs_pathways = sorted(costs_pathways, key=lambda x: x[0])
        pathways_sorted = [p for c, p in sorted_costs_pathways]
        costs_sorted = [c for c, p in sorted_costs_pathways]
        return pathways_sorted, costs_sorted

    def order_clusters_by_cost(self, clusters):
        ordered_clusters = []
        best_costs = []
        for group in clusters:
            ordered, costs = self.order_by_cost(group)
            ordered_clusters.append(ordered)
            best_costs.append(costs[0])

        costed_clusters = list(zip(ordered_clusters, best_costs))
        ordered_costed_clusters = sorted(costed_clusters, key=lambda x: x[1])
        return [cluster for cluster, best_cost in ordered_costed_clusters]

    def _add_costs_to_end_nodes(self, pathway):
        self.logger.debug(f"Getting costs for end nodes: {pathway.end_nodes}")
        for node in pathway.end_nodes:
            if pathway.sub_graph.nodes[node]['attributes']['is_starting_material'] == 1:
                pathway.sub_graph.nodes[node]['attributes']['cost'] = self.in_stock_cost
            else:
                pathway.sub_graph.nodes[node]['attributes']['cost'] = self.not_in_stock_cost

            self.logger.debug(f"Cost for {node} is {pathway.sub_graph.nodes[node]['attributes']['cost']}")

    def _iteratively_cost_other_nodes(self, pathway):
        self.logger.debug("Costing intermediates through to target")
        nodes = pathway.end_nodes
        step = 1
        while pathway.target_smiles not in nodes:
            self.logger.debug(f"step {step}")
            reactions = retrobiocat_web.retro.network_pathway.graph_functions.get_reaction_nodes_of_list_substrates(pathway.sub_graph, nodes)
            nodes = []

            for reaction in reactions:
                precursors = retrobiocat_web.retro.network_pathway.graph_functions.get_reaction_substrates(pathway.sub_graph, reaction)
                if self._are_all_precursors_costed(pathway, precursors):
                    nodes += self._add_costs_to_reaction_substrates(pathway, reaction, precursors)
            step += 1

    def _add_costs_to_reaction_substrates(self, pathway, reaction, precursors):
        product = retrobiocat_web.retro.network_pathway.graph_functions.get_reaction_products(pathway.sub_graph, reaction)[0]  # move towards target, only ever 1
        cost_precursors = self._lookup_cost_of_nodes(pathway, precursors)
        cost = (cost_precursors/self.reaction_yield) + self.reaction_cost
        pathway.sub_graph.nodes[product]['attributes']['cost'] = cost

        self.logger.debug(f"{precursors}--{reaction}-->{product} cost is {cost}")

        return [product]

    def _are_all_precursors_costed(self, pathway, precursors):
        for node in precursors:
            if 'cost' not in pathway.network.graph.nodes[node]['attributes']:
                self.logger.debug("Not all precursors are costed, skip this reaction for now")
                return False
        return True

    def _lookup_cost_of_nodes(self, pathway, list_nodes):
        cost = 0
        for node in list_nodes:
            if 'cost' not in pathway.network.graph.nodes[node]['attributes']:
                self.logger.error(f"No cost attribute in {node}, in pathway {pathway.list_nodes}")
                raise Exception(f"No cost attribute in {node}, in pathway {pathway.list_nodes}")

            cost += pathway.network.graph.nodes[node]['attributes']['cost']
        return cost


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    from retrobiocat_web.retro.network_pathway.network import Network
    from retrobiocat_web.retro.network_pathway.pathway.pathway import Pathway
    from retrobiocat_web.retro.retrosynthesis_engine.retrosynthesis_engine import RetrosynthesisEngine

    network = Network(target_smiles='c1ccc([C@@H]2CCCCN2)cc1')
    retro_engine = RetrosynthesisEngine(network)
    rxns, smis = [], ['c1ccc([C@@H]2CCCCN2)cc1']
    product_smis, reaction_names = retro_engine.custom_reaction('c1ccc([C@@H]2CCCCN2)cc1', ['c1ccc(C2=NCCCC2)cc1'],
                                                                'Imine reduction')
    rxns += reaction_names
    smis += product_smis
    product_smis, reaction_names = retro_engine.custom_reaction('c1ccc(C2=NCCCC2)cc1', ['NCCCCC(=O)c1ccccc1'],
                                                                'Imine formation')
    rxns += reaction_names
    smis += product_smis

    pathway1 = Pathway(rxns + smis,
                      network, log_level='DEBUG')
    coster = Pathway_Coster(log_level='DEBUG')
    cost = coster.cost(pathway1)
    print(cost)

    product_smis, reaction_names = retro_engine.custom_reaction('NCCCCC(=O)c1ccccc1', ['O=CCCCC(=O)c1ccccc1'],
                                                                'Aldehyde amination')
    rxns += reaction_names
    smis += product_smis

    product_smis, reaction_names = retro_engine.custom_reaction('O=CCCCC(=O)c1ccccc1', ['O=C(O)CCCC(=O)c1ccccc1'],
                                                                'Carboxylic acid reduction')
    rxns += reaction_names
    smis += product_smis

    pathway2 = Pathway(rxns + smis,
                      network, log_level='DEBUG')
    coster = Pathway_Coster(log_level='DEBUG')
    cost = coster.cost(pathway2)
    print(cost)

    print()
    pathways = coster.order_by_cost([pathway1, pathway2])
    for p in pathways:
        print(coster.cost(p))

