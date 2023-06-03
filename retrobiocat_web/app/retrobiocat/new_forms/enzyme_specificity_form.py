from decimal import Decimal
from flask_wtf import FlaskForm
from wtforms import BooleanField, DecimalField
from wtforms.validators import NumberRange

from retrobiocat_web.retro.enzyme_identification.config import Specificity_Scorer_Config


class Specificity_Scorer_Config_Form(FlaskForm):
    run_similarity_search = BooleanField('Calculate substrate similarities', default=True)
    similarity_threshold = DecimalField('Substrate similarity threshold', default=0.55, validators=[NumberRange(min=0, max=1)])
    products_only = BooleanField('Only use products for calculating similarity (faster)', default=True)
    only_active = BooleanField('Only show active enzyme activity', default=True)
    only_reviewed = BooleanField('Use only reviewed substrate specificity data', default=True)

    def form_to_config_dict(self):
        config = Specificity_Scorer_Config()
        attr_dict = {}
        for fieldname, value in self.data.items():
            if type(value) is Decimal:
                value = float(value)
            attr_dict[fieldname] = value
        config.update_from_dict(attr_dict)
        return config.to_dict()