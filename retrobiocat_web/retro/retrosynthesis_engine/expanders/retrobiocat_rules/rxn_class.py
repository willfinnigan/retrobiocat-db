from retrobiocat_web.logging import add_logger
from retrobiocat_web.mongo.model_queries import reaction_queries
from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig


class RetroBioCat_Reactions():

    def __init__(self, config=None, log_level='WARNING'):
        self.logger = add_logger('RetroBioCat_Rxns', level=log_level)

        self.config = config
        if self.config is None:
            self.config = RetrosynthesisConfig()

        self.rxns_strings = None
        self.rules_by_type = None

        query_result = reaction_queries.get_reactions(include_experimental=self.config.include_experimental,
                                                      include_two_step=self.config.include_two_step,
                                                      include_requires_absence_of_water=self.config.include_requires_absence_of_water)
        self.rxns = reaction_queries.load_rxns(query_result)
        self.multi_step_rxns = reaction_queries.load_multi_step(query_result)
        self.reactions, self.enzymes, self.reaction_enzyme_map = reaction_queries.load_reactions_and_enzymes(query_result)
        self.reactionEnzymeCofactorDict = reaction_queries.load_cofactors(query_result)

    def load_additional_info(self):
        query_result = reaction_queries.get_reactions()
        self.rxns_strings = reaction_queries.load_rxn_strings(query_result)
        self.rules_by_type = reaction_queries.load_rules_by_type(query_result)

if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    import time
    for i in range(10):
        t0 = time.time()
        rxns_class = RetroBioCat_Reactions()
        t1 = time.time()