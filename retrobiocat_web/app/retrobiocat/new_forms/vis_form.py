from decimal import Decimal
from flask_wtf import FlaskForm
from wtforms import BooleanField, DecimalField, SelectField
from wtforms.validators import NumberRange
from retrobiocat_web.retro.visualisation.config import Visualiser_Config


class Visualiser_Config_Form(FlaskForm):

    colour_arrows = SelectField('Colour reactions', default='None',
                                                    choices={'None': 'None',
                                                             'Complexity change': 'Complexity change'})

    def form_to_config_dict(self):
        config = Visualiser_Config()
        attr_dict = {}
        for fieldname, value in self.data.items():
            if type(value) is Decimal:
                value = float(value)
            attr_dict[fieldname] = value
        config.update_from_dict(attr_dict)
        return config.to_dict()