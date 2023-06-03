
from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig
from retrobiocat_web.retro.retrosynthesis_engine.reaction_processing.components.reaction_creator import \
    ReactionCreator
from retrobiocat_web.retro.retrosynthesis_engine.reaction_processing.components.reaction_processor import \
    ReactionProcessor


class ReactionParser(object):

    def __init__(self, network, config=None, log_level='WARNING'):
        self.config = config
        if self.config is None:
            self.config = RetrosynthesisConfig()

        self.creator = ReactionCreator(network, config=self.config, log_level=log_level)
        self.processor = ReactionProcessor(network, config=self.config, log_level=log_level)

    def parse(self, reactant_smi, products_dict,
              metadata_dict=None, rxn_type='',
              disallowed_products=None):
        reactions = self.creator.create_new_reactions(reactant_smi, products_dict, metadata=metadata_dict, rxn_type=rxn_type)
        reactions = self.processor.process(reactions, disallowed_products=disallowed_products)
        return reactions


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection

    make_default_connection()
    from retrobiocat_web import Network

    network = Network(target_smiles='[C@H]1(C2=CC=CC=C2)NCCCC1')
    output = {'Imine reduction': [['c1ccc(C2=NCCCC2)cc1', 'N'], ['C1=N[C@H](c2ccccc2)CCC1', 'N']]}
    reaction_creator = ReactionCreator(network, log_level='INFO')
    reactions = reaction_creator.create_new_reactions(network.target_smiles, output)

    reaction_processor = ReactionProcessor(network, log_level='DEBUG')
    reaction_processor.config.remove_small_mols = True
    reactions = reaction_processor.process(reactions)

