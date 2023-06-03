from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

import retrobiocat_web.retro.network_pathway.graph_functions
from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig

from retrobiocat_web.retro.retrosynthesis_engine.reaction_processing.components.small_molecule_remover import SmallMolRemover
from retrobiocat_web.logging import add_logger

class ReactionProcessor():
    """
    ReactionProcessor takes reactions made by ReactionCreator,
    and processes them so:
    - removes any small molecules
    - delete any reactions which go backwards
    - delete any reactions which are duplicates of reactions already in the network graph
    - delete any cyclic reactions, or any reactions with specified disallowed products

    The main function to do this is the process method
    """

    def __init__(self, network, config=None, log_level='WARNING'):
        self.network = network
        self.small_mol_remover = SmallMolRemover(config=config, log_level=log_level)

        self.config = config
        if self.config is None:
            self.config = RetrosynthesisConfig()

        self.logger = add_logger('ReactionProcessor', level=log_level)

    def process(self, reactions, disallowed_products=None, duplicates_must_match_name=True, must_change_ring_num=False):
        num_start = len(reactions)

        if not self.config.allow_cyclic_reactions:
            reactions = self.remove_cyclic_reactions(reactions)

        if disallowed_products is not None:
            reactions = self.remove_disallowed_products(reactions, disallowed_products)

        if self.config.remove_small_mols:
            reactions = self.remove_small_mols(reactions)

        if not self.config.allow_backwards:
            reactions = self.drop_backwards_reactions(reactions)

        if not self.config.allow_duplicates:
            reactions = self.drop_duplicate_reactions(reactions, duplicates_must_match_name)

        if must_change_ring_num:
            reactions = self.drop_non_ring_change_reactions(reactions)

        self.logger.info(f"{num_start - len(reactions)} reactions removed during reaction processing")
        return reactions

    def remove_cyclic_reactions(self, reactions):
        self.logger.debug(f"Removing any reactions which are cyclic ")
        to_keep = []
        for reaction in reactions:
            if not self._is_reaction_cyclic(reaction):
                to_keep.append(reaction)
        self.logger.debug(f"...{len(reactions) - len(to_keep)} removed")
        return to_keep

    def remove_disallowed_products(self, reactions, disallowed_products):
        self.logger.debug(f"Removing any reactions with disallowed products ({disallowed_products})")
        to_keep = []
        for reaction in reactions:
            if not self._does_reaction_contain_a_disallowed_product(reaction, disallowed_products):
                to_keep.append(reaction)
        self.logger.debug(f"...{len(reactions) - len(to_keep)} removed")
        return to_keep

    def remove_small_mols(self, reactions):
        self.logger.debug(f"Processing removal of small molecules")
        for reaction in reactions:
            kept, removed = self.small_mol_remover.remove(reaction.substrates)
            reaction.substrates = kept
            reaction.metadata.update({'small_mols': removed})
        return reactions

    def drop_backwards_reactions(self, reactions):
        self.logger.debug(f"Removing any backwards reactions")
        to_keep = []
        for reaction in reactions:
            if not self._does_reaction_go_backwards(reaction):
                to_keep.append(reaction)
        self.logger.debug(f"...{len(reactions) - len(to_keep)} removed")
        return to_keep

    def drop_duplicate_reactions(self, reactions, duplicates_must_match_name):
        self.logger.debug(f"Removing any duplicates of existing reactions")
        to_keep = []
        for reaction in reactions:
            if not self._does_reaction_already_exist(reaction, duplicates_must_match_name):
                to_keep.append(reaction)
        self.logger.debug(f"...{len(reactions) - len(to_keep)} removed")
        return to_keep

    def drop_non_ring_change_reactions(self, reactions):
        self.logger.debug(f"Removing any reactions which don't change ring number")
        to_keep = []
        for reaction in reactions:
            if self._does_reaction_change_ring_num(reaction):
                to_keep.append(reaction)
        self.logger.debug(f"...{len(reactions) - len(to_keep)} removed")
        return to_keep

    def _does_reaction_already_exist(self, new_reaction, duplicates_must_match_name):
        """Duplicate reactions are reactions which have a matching name, and have the same substrates"""

        existing_reactions = retrobiocat_web.retro.network_pathway.graph_functions.get_reactions_with_smi_as_reactant(self.network.graph, new_reaction.reactant.smi)
        names = retrobiocat_web.retro.network_pathway.graph_functions.get_reaction_names(existing_reactions, self.network.graph)

        for reaction, name in zip(existing_reactions, names):
            if (new_reaction.name == name) or (not duplicates_must_match_name):
                substrates = retrobiocat_web.retro.network_pathway.graph_functions.get_reaction_substrates(self.network.graph, reaction)
                if sorted(substrates) == sorted([s.smi for s in new_reaction.substrates]):
                    return True

        return False

    def _does_reaction_go_backwards(self, reaction):
        """If the precursors of the reaction substrates (smi's generated by applying rules) match
        the reaction reactant (the smi which the rules were applied to),
        then the reaction is going backwards """

        for substrate in reaction.substrates:
            if substrate.smi in self.network.graph.nodes:
                precursors = retrobiocat_web.retro.network_pathway.graph_functions.substrate_precursors(self.network.graph, substrate.smi)
                self.logger.debug(f"{substrate.smi} precursors = {precursors}")
            else:
                precursors = []

            if reaction.reactant.smi in precursors:
                return True
        return False

    def _is_reaction_cyclic(self, reaction):
        """ If any of the reaction substrates matches the reaction reactant,
        then the reaction is cyclic"""
        for substrate in reaction.substrates:
            if substrate.smi == reaction.reactant.smi:
                return True
        return False

    def _does_reaction_contain_a_disallowed_product(self, reaction, disallowed_products: list):
        for substrate in reaction.substrates:
            if substrate.smi in disallowed_products:
                return True
        return False

    def _does_reaction_change_ring_num(self, reaction):
        mol_r = Chem.MolFromSmiles(reaction.substrates[0].smi)
        mol_p = Chem.MolFromSmiles(reaction.reactant.smi)
        change = rdMolDescriptors.CalcNumRings(mol_p) - rdMolDescriptors.CalcNumRings(mol_r)
        if change >= 1:
            return True
        return False



if __name__ == '__main__':
    from retrobiocat_web import Network
    from retrobiocat_web.retro.retrosynthesis_engine.reaction_processing.components.reaction_creator import \
        ReactionCreator
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    network = Network(target_smiles='[C@H]1(C2=CC=CC=C2)NCCCC1')
    output = {'Imine reduction': [['c1ccc(C2=NCCCC2)cc1', 'N'], ['C1=N[C@H](c2ccccc2)CCC1', 'N']]}

    reaction_creator = ReactionCreator(network, log_level='DEBUG')
    reactions = reaction_creator.create_new_reactions(network.target_smiles, output)

    reaction_processor = ReactionProcessor(network, log_level='DEBUG')
    reaction_processor.config.remove_small_mols = False
    reactions = reaction_processor.process(reactions)