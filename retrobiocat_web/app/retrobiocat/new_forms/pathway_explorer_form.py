from decimal import Decimal

from flask_wtf import FlaskForm
from wtforms import StringField, BooleanField, SubmitField, IntegerField, DecimalField, SelectField
from wtforms.validators import DataRequired, NumberRange, ValidationError

from retrobiocat_web.app.retrobiocat.new_forms.custom_validators import is_accepted_by_rdkit


class PathwayExploreForm(FlaskForm):
    target_smiles = StringField('Target SMILES', validators=[DataRequired(), is_accepted_by_rdkit])
    number_steps = IntegerField('Max steps', default=4, validators=[NumberRange(min=1, max=5)])
    max_nodes = IntegerField('Max nodes', default=400, validators=[NumberRange(min=100, max=800)])

    weight_complexity = IntegerField('Weight Complexity Change', default=1, validators=[NumberRange(min=-1, max=10)])
    weight_num_enzymes = IntegerField('Weight Number of Enzymes', default=1, validators=[NumberRange(min=-1, max=10)])
    weight_starting = IntegerField('Weight Starting Material', default=1, validators=[NumberRange(min=-1, max=10)])
    weight_known_enzymes = IntegerField('Weight Known Enzyme Steps', default=1, validators=[NumberRange(min=-1, max=10)])
    weight_diversity = IntegerField('Weight Diversity', default=1, validators=[NumberRange(min=-1, max=10)])

    max_pathways = IntegerField('Max pathways to initially make', default=40000, validators=[NumberRange(min=100, max=80000)])
    keep_best_pathways = IntegerField('Keep only best pathways after scoring', default=2500, validators=[NumberRange(min=50, max=100000)])
    min_weight = DecimalField('Min complexity weight', default=1, validators=[NumberRange(min=0.1, max=2)])

    submit = SubmitField('Start')

    def form_to_config_dict(self):
        attr_dict = {}
        for fieldname, value in self.data.items():
            if type(value) is Decimal:
                value = float(value)
            attr_dict[fieldname] = value
        return attr_dict