from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig
from retrobiocat_web.logging import add_logger
from retrobiocat_web.retro.retrosynthesis_engine.expanders.retrorules.retrorules_getter import RetroRules_Getter
from retrobiocat_web.retro.retrosynthesis_engine.expanders.retrorules.retrorules_processor import RetroRules_Processor
from retrobiocat_web.retro.retrosynthesis_engine.reaction_processing.components.reaction_creator import ReactionCreator
from retrobiocat_web.retro.retrosynthesis_engine.reaction_processing.components.reaction_processor import \
    ReactionProcessor
from retrobiocat_web.retro.retrosynthesis_engine.reaction_processing.reaction_ranker import ReactionRanker
from retrobiocat_web.retro.retrosynthesis_engine.rule_application.rule_applicator import RuleApplicator


class RetroRules_Expander():

    def __init__(self, network, config=None, log_level='WARNING'):

        self.logger = add_logger('RetroRules_Expander', level=log_level)

        self.config = config
        if self.config is None:
            self.config = RetrosynthesisConfig()

        self.rule_applicator = RuleApplicator(config=self.config, log_level=log_level)
        self.reaction_creator = ReactionCreator(network, config=self.config, log_level=log_level)
        self.reaction_processor = ReactionProcessor(network, config=self.config, log_level=log_level)
        self.reaction_ranker = ReactionRanker(config=self.config, log_level=log_level)

        self.retrorules_getter = RetroRules_Getter(config=self.config, log_level=log_level)
        self.retrorules_processor = RetroRules_Processor(config=self.config, log_level=log_level)

        self.options_map = {}

    def expand(self, smi, disallowed_products=None):

        smarts, metadata = self.retrorules_getter.get_rules(smi)

        rxns = {}
        for rule_id, smarts in smarts.items():
            rxns[rule_id] = self.rule_applicator.smarts_to_rdkit(smarts)

        rxnsSubstratesDict = self.rule_applicator.apply_rdkit(smi, rxns)

        rxnsSubstratesDict = self.retrorules_processor.process(rxnsSubstratesDict)

        reactions = self.reaction_creator.create_new_reactions(smi, rxnsSubstratesDict,
                                                                metadata=metadata,
                                                                rxn_type='retrorules')

        self._convert_name_to_ec_where_possible(reactions)
        reactions = self.reaction_processor.process(reactions, duplicates_must_match_name=self.config.rr_duplicates_must_match_name)
        reactions = self.reaction_ranker.rank_by_complexity(reactions, self.config.max_reactions)
        return reactions

    def _convert_name_to_ec_where_possible(self, reactions):
        for reaction in reactions:
            if len(reaction.metadata['ec_numbers']) > 0:
                reaction.name = reaction.metadata['ec_numbers'][0]

    def get_expansion_options(self, smi):
        if smi in self.options_map:
            return self.options_map[smi]

        smarts, metadata = self.retrorules_getter.get_rules(smi)
        options = []

        for i, name in enumerate(smarts):
            opt = ExpansionOption_RetroRules(self, smi, name, smarts[name], metadata[name], num=i+1)
            options.append(opt)

        self.options_map[smi] = options

        return options



class ExpansionOption_RetroRules():
    def __init__(self, expander, smi, name, smarts, metadata, num=None):
        self.expander = expander
        self.smi = smi
        self.name = name
        self.smarts = smarts
        self.metadata = metadata
        self.score = metadata['score']
        self.similarity = metadata['similarity']
        self.reactions = None
        self.processed_reactions = None
        self.num = num

        if len(self.metadata['ec_numbers']) > 0:
            self.name = self.metadata['ec_numbers'][0]

    def __repr__(self):
        return f"Num. {self.num} RetroRules Expansion Option for {self.smi}"

    def get_processed_reactions(self):
        """ Return the processed reaction objects for this reaction """

        if self.processed_reactions is None:
            reactions = self.get_reactions()
            self.processed_reactions = self._process_reactions(reactions)

        return self.processed_reactions

    def get_reactions(self):
        """ Return the unprocessed reaction objects for this reaction """
        if self.reactions is None:
            self.reactions = self._execute_option()

        return self.reactions

    def _execute_option(self):
        rxns = {self.name: self.expander.rule_applicator.smarts_to_rdkit(self.smarts)}
        metadata = {self.name: self.metadata}
        rxnsSubstratesDict = self.expander.rule_applicator.apply_rdkit(self.smi, rxns)

        reactions = self.expander.reaction_creator.create_new_reactions(self.smi, rxnsSubstratesDict,
                                                                        metadata=metadata,
                                                                        rxn_type='retrorules')
        return reactions

    def _process_reactions(self, reactions):
        reactions = self.expander.reaction_processor.process(reactions, duplicates_must_match_name=self.expander.mcts_config.rr_duplicates_must_match_name)
        return reactions


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    from retrobiocat_web.retro.network_pathway.network import Network
    network = Network(target_smiles='CCCC=O')
    rr_expander = RetroRules_Expander(network, log_level='WARNING')

    reactions = rr_expander.expand('CCCC=O')


