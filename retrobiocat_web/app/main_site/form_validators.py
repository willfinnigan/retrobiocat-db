from wtforms import ValidationError

import retrobiocat_web.retro.network_pathway.rdkit_utils
from retrobiocat_web.mongo.model_queries import specificity_data_query


def is_accepted_by_rdkit(form, field):
    if retrobiocat_web.retro.network_pathway.rdkit_utils.rdkit_smile(field.data) == None:
        if field.data != '':
            raise ValidationError('SMILES not accepted by rdkit')

def is_reaction(form, field):
    reaction_names = list(field.data.split(", "))
    for reaction in reaction_names:
        if reaction not in (specificity_data_query.get_reactions_in_db() + ['All']):
            if reaction != '':
                raise ValidationError('Reaction not defined in main_site')

def is_enzyme(form, field):
    enzyme_names = list(field.data.split(", "))
    for enzyme in enzyme_names:
        if enzyme not in (specificity_data_query.get_enzymes_in_db() + ['All']):
            if enzyme != ['']:
                raise ValidationError('Enzyme not defined in main_site')
