class RetrosynthesisConfig():

    def __init__(self):

        self.target_always_not_buyable = True

        # rule application
        self.combine_enantiomers = True
        self.check_chiral_products = True  # best to leave this
        self.allow_chiral_symmetry = False  # best to leave this
        self.clean_brackets = True
        self.force_rdkit_smis = True

        # reaction parsing
        self.allow_backwards = False
        self.allow_duplicates = False
        self.allow_cyclic_reactions = False
        self.remove_small_mols = False
        self.small_precursors = ['N', 'O', 'O=O', 'H+', '[H+]']

        # evaluators
        self.source_mol_mode = 'building_blocks'
        self.source_mol_vendors = ['askcos', 'sigma', 'apollo', 'flurochem', 'alfa', 'lifechem']
        self.source_mols_can_be_chiral = True

        # graph control
        self.prune_on_substrates = False

        # expanders general
        self.max_reactions = None

        # retrobiocat
        self.include_experimental = False
        self.include_two_step = True
        self.include_requires_absence_of_water = False

        # retrorules
        self.rr_diameter = 6
        self.rr_threshold = 0.25
        self.rr_score_threshold = 0
        self.rr_combined_score_threshold = 0.25
        self.rr_remove_cofactors = True
        self.rr_group_duplicates = True
        self.rr_duplicates_must_match_name = True

        # aisynthfinder
        self.aizynth_reaction_filter = 'policy'  # or complexity or policy
        self.aizynth_filter_cutoff = 0.05
        self.aizynth_use_feasability_filter = False
        self.aizynth_template_column = 'retro_template'
        self.aizynth_cutoff_cumulative = 0.995
        self.aizynth_cutoff_number = 50
        self.aizynth_group_duplicates = True
        self.aizynth_group_policy = 'sum'  # or max

        # ringbreaker
        self.ringbreaker_cutoff = 10
        self.ringbreaker_cutoff_cumulative = 0.995


        # custom_rule_set_expander
        self.custom_expander_rule_mode = 'rdchiral'

    def update_from_dict(self, attr_dict):
        current_dict = self.to_dict()
        for key, value in attr_dict.items():
            if key in current_dict:
                setattr(self, key, value)
        return self

    def to_dict(self):
        return self.__dict__

if __name__ == '__main__':
    config = RetrosynthesisConfig()
    dict_attr = config.to_dict()
    new_config = RetrosynthesisConfig().update_from_dict(dict_attr)
