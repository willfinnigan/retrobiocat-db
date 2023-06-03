from flask import url_for, redirect, render_template
from flask_wtf import FlaskForm
from wtforms import SelectField, IntegerField, StringField, BooleanField, SubmitField
from wtforms.fields.html5 import DecimalField
from wtforms.validators import DataRequired, NumberRange

from retrobiocat_web.app.database_bps.db_analysis import bp
from retrobiocat_web.app.main_site.form_validators import is_reaction, is_enzyme, is_accepted_by_rdkit
from retrobiocat_web.mongo.model_queries import activity_queries


# smi_col_choices = [('product_1_smiles', 'Product'),
#                    ('substrate_1_smiles', 'Substrate 1'),
#                    ('substrate_2_smiles', 'Substrate 2')]


class SubstrateForm(FlaskForm):
    enzyme_type = SelectField('Enzyme type', validators=[DataRequired(), is_enzyme])
    reaction = SelectField('Reaction', validators=[is_reaction])
    num_choices = IntegerField('Max enzymes per product', default=10, validators=[NumberRange(min=1, max=100)])
    target_smiles = StringField('Product SMILES', validators=[is_accepted_by_rdkit])
    similarity = DecimalField('Similarity cutoff', default=0.5, validators=[NumberRange(min=0.05, max=1)])
    submit = SubmitField('Submit')

    def set_choices(self):
        self.enzyme_type.choices = [(c, c) for c in ['All'] + activity_queries.enzyme_types_for_which_there_are_activities()]
        self.reaction.choices = [(c, c) for c in ['All'] + activity_queries.reactions_for_which_there_are_activities()]

    def validate(self):
        if not FlaskForm.validate(self):
            return False
        result = True
        seen = set()

        if self.enzyme_type.data == 'All' and self.reaction.data == 'All' and self.target_smiles.data == "":
            self.target_smiles.errors.append("Can't both be All with no product")
            self.enzyme_type.errors.append("Can't both be All with no product")
            self.reaction.errors.append("Can't both be All with no product")
            result = False
        else:
            seen.add(self.enzyme_type.data)
            seen.add(self.reaction.data)
        return result

def reaction_enzyme_and_target_smiles_are_all_empty(form_data):
    if form_data['enzyme_type'] == '' and form_data['reaction'] == '' and form_data['target_smiles'] == '':
        return True
    return False

@bp.route('/substrate_specificity_form',  methods=['GET', 'POST'])
def substrate_specificity_form():
    form = SubstrateForm()
    form.set_choices()

    if form.validate_on_submit() == True:
        if reaction_enzyme_and_target_smiles_are_all_empty(form.data):
            return redirect(url_for('.substrate_specificity_form'))
        else:
            keys = ['enzyme_type', 'reaction', 'target_smiles', 'similarity', 'num_choices']
            kwargs = {k: form.data[k] for k in keys}
            return redirect(url_for('.substrate_specificity', **kwargs))

    return render_template('substrate_similarity_search/substrate_similarity_form_page.html', form=form)
