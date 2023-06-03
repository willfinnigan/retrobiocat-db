from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig


class AIZynth_Grouper():
    """ Group reactions with the same name and substrates together into a single reaction """

    def __init__(self, config=None):
        self.config = config
        if self.config is None:
            self.config = RetrosynthesisConfig()

    def group_and_select_representatives(self, reactions):
        groups = self.group(reactions)
        grouped_reactions = []
        for name, reaction_groups in groups.items():
            for reactions in reaction_groups:
                rep_rxn = self._get_representative_reaction(reactions)
                grouped_reactions.append(rep_rxn)
        return grouped_reactions

    def _get_representative_reaction(self, reactions):
        rep_rxn = reactions.pop(0)
        rep_rxn.metadata['num_rxns_grouped'] = 1

        for rxn in reactions:
            rep_rxn.metadata['num_rxns_grouped'] += 1
            if self.config.aizynth_group_policy == 'sum':
                rep_rxn.metadata['policy_probability'] += rxn.metadata['policy_probability']
                rep_rxn.metadata['library_occurence'] += rxn.metadata['library_occurence']
        return rep_rxn


    def group(self, reactions):
        groups = self._group_by_name(reactions)
        for name, reactions in groups.items():
            groups[name] = self._group_by_substrate(reactions)
        return groups

    def _group_by_substrate(self, reactions):
        same_substrates = {}
        for reaction in reactions:
            rxn_substrates = str(sorted([s.smi for s in reaction.substrates]))
            if rxn_substrates not in same_substrates:
                same_substrates[rxn_substrates] = []
            same_substrates[rxn_substrates].append(reaction)

        same_substrates_lists = []
        for sub_string, reacs in same_substrates.items():
            same_substrates_lists.append(reacs)

        return same_substrates_lists

    def _group_by_name(self, reactions):
        matching_names = {}
        for reaction in reactions:
            name = self._name_without_numbering(reaction)
            if name not in matching_names:
                matching_names[name] = []
            matching_names[name].append(reaction)
        return matching_names

    def _name_without_numbering(self, reaction):
        i = reaction.name.find("__")
        name = reaction.name[0:i]
        return name





