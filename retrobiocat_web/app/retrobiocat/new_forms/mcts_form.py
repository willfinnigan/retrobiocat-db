from decimal import Decimal

from flask_wtf import FlaskForm
from wtforms import StringField, BooleanField, SubmitField, IntegerField, DecimalField, SelectField
from wtforms.validators import DataRequired, NumberRange

from retrobiocat_web.app.retrobiocat.new_forms.custom_validators import is_accepted_by_rdkit


class MCTSExploreForm(FlaskForm):
    target_smiles = StringField('Target SMILES', validators=[DataRequired(), is_accepted_by_rdkit])

    exploration = DecimalField("Exploration", default=2, validators=[NumberRange(min=0.1, max=1000)])
    max_search_time = IntegerField("Search time", default=120, validators=[NumberRange(min=10, max=240)])
    maxLength = IntegerField('Max steps', default=4, validators=[NumberRange(min=1, max=6)])

    # expanders
    chemistry_expansion = BooleanField('Chemistry Expansion (AIZynthfinder)', default=True)
    biocatalysis_expansion = BooleanField('Biocatalysis Expansion (RetroBioCat)', default=True)
    biosynthesis_expansion = BooleanField('Biosynthesis Expansion (RetroRules)', default=False)

    # expansion optioons
    stop_expansion_if_nonbuyable_at_max_length = BooleanField('Stop expansion if a max length node is not available', default=True)

    # ucb
    adjust_non_buyable_penalties_with_complexity = BooleanField('Use complexity to adjust penalties for mol not in stock', default=True)
    min_penalty_for_non_buyable = DecimalField("Minimum penalty for a mol not in stock", default=-0.1, validators=[NumberRange(min=-1, max=0)])
    complexity_for_minus_1 = DecimalField("Complexity limit for max penalty", default=0, validators=[NumberRange(min=-0.5, max=0.5)])
    complexity_for_zero = DecimalField("Complexity limit for no penalty", default=-1, validators=[NumberRange(min=-2, max=0)])
    only_solved = BooleanField('Only return solved pathways', default=True)

    # clustering
    dont_cluster_single_step_routes = BooleanField("Don't cluster single step routes", default=True)
    distance_threshold = DecimalField("Clustering distance threshold", default=2, validators=[NumberRange(min=0.1, max=10)])

    # costing
    in_stock_cost = IntegerField("In stock cost", default=1, validators=[NumberRange(min=1, max=10)])
    not_in_stock_cost = IntegerField("In stock cost", default=10, validators=[NumberRange(min=1, max=100)])
    reaction_cost = IntegerField("Reaction cost", default=1, validators=[NumberRange(min=1, max=10)])
    reaction_yield = DecimalField("Reaction yield", default=0.8, validators=[NumberRange(min=0.1, max=1)])

    submit = SubmitField('Start')

    def form_to_config_dict(self):
        attr_dict = {}
        for fieldname, value in self.data.items():
            if type(value) is Decimal:
                value = float(value)
            attr_dict[fieldname] = value
        return attr_dict

