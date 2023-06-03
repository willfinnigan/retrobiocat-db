from decimal import Decimal
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, IntegerField
from wtforms.validators import DataRequired, NumberRange

from retrobiocat_web.app.retrobiocat.new_forms.custom_validators import is_accepted_by_rdkit

class NetworkExploreForm(FlaskForm):
    target_smiles = StringField('Target SMILES', validators=[DataRequired(), is_accepted_by_rdkit])
    number_steps = IntegerField('Initial steps', default=0, validators=[NumberRange(min=0, max=10)])
    max_initial_nodes = IntegerField('Max initial nodes', default=20, validators=[NumberRange(min=1, max=80)])
    submit = SubmitField('Start')


    def form_to_config_dict(self):
        attr_dict = {}
        for fieldname, value in self.data.items():
            if type(value) is Decimal:
                value = float(value)
            attr_dict[fieldname] = value
        return attr_dict