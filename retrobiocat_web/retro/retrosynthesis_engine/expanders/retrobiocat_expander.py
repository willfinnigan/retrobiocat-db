from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig
from retrobiocat_web.logging import add_logger
from retrobiocat_web.retro.retrosynthesis_engine.expanders.retrobiocat_rules.rxn_class import \
    RetroBioCat_Reactions
from retrobiocat_web.retro.retrosynthesis_engine.reaction_processing.components.reaction_creator import ReactionCreator
from retrobiocat_web.retro.retrosynthesis_engine.reaction_processing.components.reaction_processor import \
    ReactionProcessor
from retrobiocat_web.retro.retrosynthesis_engine.reaction_processing.reaction_parser import ReactionParser
from retrobiocat_web.retro.retrosynthesis_engine.reaction_processing.reaction_ranker import ReactionRanker
from retrobiocat_web.retro.retrosynthesis_engine.rule_application.rule_applicator import RuleApplicator


class RetroBioCat_Expander():

    def __init__(self, network, rxn_class=None, config=None, log_level='WARNING'):

        self.logger = add_logger('RetroBioCat_Expander', level=log_level)

        self.config = config
        if self.config is None:
            self.config = RetrosynthesisConfig()

        self.rule_applicator = RuleApplicator(config=config, log_level=log_level)
        self.reaction_creator = ReactionCreator(network, config=config, log_level=log_level)
        self.reaction_processor = ReactionProcessor(network, config=config, log_level=log_level)
        self.reaction_ranker = ReactionRanker(config=config, log_level=log_level)

        self.rxn_class = rxn_class
        if self.rxn_class is None:
            self.rxn_class = RetroBioCat_Reactions(config=config, log_level=log_level)
            self.logger.debug("Loaded RetroBioCat reactions")

        self.options_map = {}  # smis as keys, list_options as values

    def expand(self, smi):
        reactions = self._get_reactions(smi)
        return reactions

    def get_expansion_options(self, smi):
        if smi in self.options_map:
            return self.options_map[smi]

        reactions = self._get_reactions(smi)
        options = self._reactions_to_options(smi, reactions)

        self.options_map[smi] = options

        return options

    def _get_reactions(self, smi):
        result = self.rule_applicator.apply_rdchiral(smi, self.rxn_class.rxns,
                                                     multistep_rxns=self.rxn_class.multi_step_rxns)
        reactions = self.reaction_creator.create_new_reactions(smi, result, rxn_type='retrobiocat')
        reactions = self.reaction_processor.process(reactions)
        reactions = self.reaction_ranker.rank_by_complexity(reactions, self.config.max_reactions)
        [self._add_metadata_to_reaction(reaction) for reaction in reactions]
        return reactions

    def _reactions_to_options(self, smi, reactions):
        options = []
        for i, reaction in enumerate(reactions):
            opt = ExpansionOption_RetroBioCat(self, smi, reaction, num=i+1)
            options.append(opt)
        return options

    def _add_metadata_to_reaction(self, reaction):
        reaction.metadata['possible_enzymes'] = self.rxn_class.reaction_enzyme_map.get(reaction.name, [''])
        reaction.metadata['selected_enzyme'] = reaction.metadata['possible_enzymes'][0]
        reaction.metadata['enzyme_cofactors'] = self.rxn_class.reactionEnzymeCofactorDict.get(reaction.name, {})

        reaction.metadata['is_enzyme'] = 1
        if reaction.metadata['selected_enzyme'] != 'Chemical':
            reaction.metadata['is_enzyme'] = 0

class ExpansionOption_RetroBioCat():
    def __init__(self, expander, smi, reaction, num=None):
        self.expander = expander
        self.smi = smi
        self.reaction = reaction
        self.processed_reactions = None
        self.score = reaction.complexity_change
        self.num = num

    def __repr__(self):
        return f"{self.num}. {self.reaction.name} for {self.smi}"

    def get_processed_reactions(self):
        """ Return the processed reaction objects for this reaction """
        if self.processed_reactions is None:
            reactions = self.get_reactions()
            self.processed_reactions = self.expander.reaction_processor.process([self.reaction])

        return self.processed_reactions

    def get_reactions(self):
        """ Return the unprocessed reaction objects for this reaction """
        return [self.reaction]




if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    from retrobiocat_web.retro.network_pathway.network import Network
    network = Network(target_smiles='CCCC=O')
    rbc_expander = RetroBioCat_Expander(network, log_level='INFO')
    reactions = rbc_expander.expand('CCCC=O')