from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig
from retrobiocat_web.logging import add_logger
from retrobiocat_web.retro.retrosynthesis_engine.expanders.aizynthfinder.aizynthfinder_actions import \
    ActionGetter, ActionFilter
from retrobiocat_web.retro.retrosynthesis_engine.expanders.aizynthfinder.aizynthfinder_grouper import AIZynth_Grouper
from retrobiocat_web.retro.retrosynthesis_engine.reaction_processing.components.reaction_creator import ReactionCreator
from retrobiocat_web.retro.retrosynthesis_engine.reaction_processing.components.reaction_processor import \
    ReactionProcessor
from retrobiocat_web.retro.retrosynthesis_engine.reaction_processing.reaction_parser import ReactionParser
from retrobiocat_web.retro.retrosynthesis_engine.reaction_processing.reaction_ranker import ReactionRanker
from retrobiocat_web.retro.retrosynthesis_engine.rule_application.rule_applicator import RuleApplicator


class AIZynthfinder_Expander():

    def __init__(self, network, config=None, log_level='WARNING'):

        self.logger = add_logger('AIZythfinder_Expander', level=log_level)

        self.config = config
        if self.config is None:
            self.config = RetrosynthesisConfig()

        self.rule_applicator = RuleApplicator(config=config, log_level=log_level)

        self.reaction_creator = ReactionCreator(network, config=self.config, log_level=log_level)
        self.reaction_processor = ReactionProcessor(network, config=self.config, log_level=log_level)
        self.reaction_ranker = ReactionRanker(config=config, log_level=log_level)

        self.aizynth_getter = ActionGetter(config=config, log_level=log_level)
        self.aizynth_filter = ActionFilter(config=config, log_level=log_level)
        self.aizynth_grouper = AIZynth_Grouper(config=config)

        self.options_map = {}  # smis as keys, list_options as values

    def expand(self, smi):

        # expand is generate options, create reactions.
        options = self.get_expansion_options(smi)
        reactions = []
        for opt in options:
            reactions.extend(opt.get_reactions())

        # process reactions as a batch rather than individually
        if self.config.aizynth_group_duplicates:
            reactions = self.aizynth_grouper.group_and_select_representatives(reactions)
        reactions = self.reaction_processor.process(reactions)

        # rank the reactions
        if self.config.aizynth_reaction_filter == 'complexity':
            reactions = self.reaction_ranker.rank_by_complexity(reactions, self.config.max_reactions)
        elif self.config.aizynth_reaction_filter == 'policy':
            reactions = self.reaction_ranker.rank_by_aizynthfinder_policy_score(reactions, self.config.max_reactions)
        else:
            self.logger.error(f"Error - aizynthfinder reaction filter mode {self.config.aizynth_reaction_filter} not recognised")

        return reactions

    def get_expansion_options(self, smi):
        if smi in self.options_map:
            return self.options_map[smi]

        smarts_dict, metadata = self.aizynth_getter.get_rxns(smi)
        options = []

        for i, name in enumerate(smarts_dict):
            opt = ExpansionOption_AIZynthfinder(self, smi, name, smarts_dict[name], metadata[name], num=i+1)
            options.append(opt)

        self.options_map[smi] = options

        return options

class ExpansionOption_AIZynthfinder():

    def __init__(self, expander, smi, name, smarts, metadata, num=None):
        self.expander = expander
        self.smi = smi
        self.name = name
        self.smarts = smarts
        self.metadata = metadata
        self.score = metadata['policy_probability']
        self.reactions = None
        self.processed_reactions = None
        self.num = num

    def __repr__(self):
        return f"Num. {self.num} AIZynthfinder Expansion Option for {self.smi}"

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
        rxns = {self.name: self.expander.rule_applicator.smarts_to_rdchiral(self.smarts)}
        metadata = {self.name: self.metadata}
        rxnsSubstratesDict = self.expander.rule_applicator.apply_rdchiral(self.smi, rxns)

        if self.expander.config.aizynth_use_feasability_filter:
            rxnsSubstratesDict, metadata = self.expander.aizynth_filter.filter(self.smi, rxnsSubstratesDict, metadata)

        reactions = self.expander.reaction_creator.create_new_reactions(self.smi, rxnsSubstratesDict,
                                                                        metadata=metadata,
                                                                        rxn_type='aizynthfinder')

        return reactions

    def _process_reactions(self, reactions):
        if self.expander.config.aizynth_group_duplicates:
            reactions = self.expander.aizynth_grouper.group_and_select_representatives(reactions)
        reactions = self.expander.reaction_processor.process(reactions, duplicates_must_match_name=False)
        return reactions






if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    from retrobiocat_web.retro.network_pathway.network import Network

    network = Network(target_smiles='CCO', log_level='DEBUG')

    expander = AIZynthfinder_Expander(network, log_level='DEBUG')
    options = expander.get_expansion_options("CCO")
    reaction = options[0].get_processed_reactions()
