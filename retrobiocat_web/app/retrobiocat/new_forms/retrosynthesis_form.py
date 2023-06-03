from decimal import Decimal

from flask_wtf import FlaskForm
from wtforms import StringField, BooleanField, SubmitField, IntegerField, DecimalField, SelectField, SelectMultipleField
from wtforms.validators import NumberRange, AnyOf

from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig


class Retrosynthesis_Config_Form(FlaskForm):

    # rule application
    combine_enantiomers = BooleanField('Combine enantiomers for racemic starting material', default=True)

    # reaction parsing
    allow_backwards = BooleanField('Allow backwards steps')
    remove_small_mols = BooleanField('Remove small molecules (eg NH3)', default=False)

    # source molecules
    source_mols_can_be_chiral = BooleanField('Source mols can be chiral', default=True)
    source_mol_mode = SelectField('Mode', default='building_blocks',
                                  choices={'building_blocks': 'building_blocks',
                                   'metabolites': 'metabolites'})

    # building blocks
    zinc = BooleanField("zinc", default=False)
    mcule = BooleanField("mcule", default=False)
    molport = BooleanField("molport", default=False)
    askcos = BooleanField("askcos", default=True)
    sigma = BooleanField("sigma", default=True)
    apollo = BooleanField("apollo", default=True)
    flurochem = BooleanField("flurochem", default=True)
    alfa = BooleanField("alfa", default=True)
    lifechem = BooleanField("lifechem", default=True)

    # metabolites
    ecoli = BooleanField("ecoli", default=True)

    # max reactions before applying filtering
    max_reactions = IntegerField('Max reactions', default=50, validators=[NumberRange(min=1, max=100)])

    # retrobiocat
    include_experimental = BooleanField('Include experimental reaction rules', default=False)
    include_two_step = BooleanField(
        'Include reaction rules for multi-step reactions, which are also included as single step rules', default=True)
    include_requires_absence_of_water = BooleanField('Include reactions which require an absence of water', default=False)
    only_reviewed = BooleanField('Use only reviewed substrate specificity data', default=True)

    # retrorules
    rr_diameter = IntegerField('RetroRules diameter', default=8, validators=[AnyOf([2,4,6,8,10,12,14,16])])
    rr_threshold = DecimalField('RetroRules similarity threshold', default=0.25, validators=[NumberRange(min=0, max=1)])
    rr_score_threshold = DecimalField('RetroRules biological score threshold', default=0, validators=[NumberRange(min=0, max=1)])
    rr_combined_score_threshold = DecimalField('RetroRules combine score threshold', default=0.25, validators=[NumberRange(min=0, max=1)])
    rr_remove_cofactors = BooleanField('Remove cofactors from RetroRules output', default=True)
    rr_duplicates_must_match_name = BooleanField('Duplicate reactions must have matching reaction name', default=True)


    # aisynthfinder
    aizynth_reaction_filter = SelectField('Reaction filter', default='policy',
                                            choices={'policy': 'policy',
                                                     'complexity': 'complexity'})
    aizynth_use_feasability_filter = BooleanField('Use the AIZynthfinder feasability filter', default=False)
    aizynth_filter_cutoff = DecimalField('AIZynthfinder feasability filter cutoff', default=0.05, validators=[NumberRange(min=0, max=1)])
    aizynth_template_column = 'retro_template'
    aizynth_cutoff_cumulative = DecimalField('AIZynthfinder cumulative cutoff', default=0.995, places=3, validators=[NumberRange(min=0, max=1)])
    aizynth_cutoff_number = IntegerField('AIZynthfinder max reactions', default=50, validators=[NumberRange(min=1, max=100)])

    def form_to_config_dict(self):
        config = RetrosynthesisConfig()
        attr_dict = {}
        for fieldname, value in self.data.items():
            if type(value) is Decimal:
                value = float(value)
            attr_dict[fieldname] = value
        attr_dict['source_mol_vendors'] = self.get_vendors()
        config.update_from_dict(attr_dict)
        return config.to_dict()

    def get_vendors(self):
        if self.source_mol_mode.data == 'building_blocks':
            return self._get_building_block_vendors()
        elif self.source_mol_mode.data == 'metabolites':
            return self._get_metabolites_vendors()
        else:
            Exception(f'Source mol mode is invalid ({self.mode.data})')
            return []

    def _get_building_block_vendors(self):
        vendors = []
        if self.askcos.data:
            vendors.append('askcos')
        if self.zinc.data:
            vendors.append('zinc')
        if self.mcule.data:
            vendors.append('mcule')
        if self.molport.data:
            vendors.append('molport')
        if self.sigma.data:
            vendors.append('sigma')
        if self.apollo.data:
            vendors.append('apollo')
        if self.flurochem.data:
            vendors.append('flurochem')
        if self.alfa.data:
            vendors.append('alfa')
        if self.lifechem.data:
            vendors.append('lifechem')

        return vendors

    def _get_metabolites_vendors(self):
        vendors = []
        if self.ecoli.data:
            vendors.append('ecmdb')

        return vendors




if __name__ == '__main__':
    from retrobiocat_web.app.app import create_app
    app = create_app()
    with app.test_request_context('/'):
        f = Retrosynthesis_Config_Form()





