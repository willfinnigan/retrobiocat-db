from flask_wtf import FlaskForm
from wtforms import StringField, BooleanField, SubmitField, IntegerField, DecimalField, SelectField
from wtforms.validators import DataRequired, NumberRange, ValidationError
from retrobiocat_web.mongo.models.biocatdb_models import SSN_record, EnzymeType, Sequence, Activity
from retrobiocat_web.mongo.models.reaction_models import Reaction

def enzyme_type_selected(form, field):
    if field.data == '-':
        raise ValidationError('Must choose an enzyme type')

class SSN_Form(FlaskForm):
    enzyme_type = SelectField('Enzyme name')
    alignment_score = IntegerField('Score', validators=[NumberRange(min=40, max=300)])
    hide_mutants = BooleanField('Hide mutants', default=False)
    only_biocatdb = BooleanField('Only RetroBioCat database sequences', default=False)
    submit = SubmitField('Submit')

    def set_choices(self):
        ssn_records = SSN_record.objects().distinct('enzyme_type')
        enzyme_types = EnzymeType.objects()

        list_enzyme_types = []
        enzyme_descriptions = {}
        for enz_type in enzyme_types:
            if enz_type in ssn_records:
                enzyme_descriptions[enz_type.enzyme_type] = f"{enz_type.enzyme_type} - {enz_type.full_name}"
                list_enzyme_types.append(enz_type.enzyme_type)

        list_enzyme_types = sorted(list_enzyme_types)

        self.enzyme_type.choices = []
        for key in list_enzyme_types:
            self.enzyme_type.choices.append((key, enzyme_descriptions[key]))


class HeatmapForm(FlaskForm):
    enzyme_type = SelectField('Enzyme type')
    reaction = SelectField('Reaction')
    only_reviewed = BooleanField('Only reviewed data', default=True)
    molecule = SelectField('Reaction', choices=[('product_1_smiles', 'Product'),
                                                ('substrate_1_smiles', 'Substrate 1'),
                                                ('substrate_2_smiles', 'Substrate 2')])

    max_molecules = IntegerField(default=1000)
    max_enzymes = IntegerField(default=200)

    remove_negative = SelectField('Remove negative data', choices=[('False', 'Do not remove'),
                                                                   ('only_negative',
                                                                    'Remove molecules where only negative'),
                                                                   ('all', 'Remove all negative data')])
    submit = SubmitField('Generate heatmap')

    def set_choices(self):

        self.reaction.choices = [(c, c) for c in ['-'] + (list(Activity.objects().distinct('reaction')))]

        enzyme_types_with_data = list(Activity.objects().distinct('enzyme_type'))
        enzyme_types = EnzymeType.objects(enzyme_type__in=enzyme_types_with_data)

        list_enzyme_types = []
        enzyme_descriptions = {}
        for enz_type in enzyme_types:
            enzyme_descriptions[enz_type.enzyme_type] = f"{enz_type.enzyme_type} - {enz_type.full_name}"
            list_enzyme_types.append(enz_type.enzyme_type)

        list_enzyme_types = sorted(list_enzyme_types)

        self.enzyme_type.choices = [('-', '-')]
        for key in list_enzyme_types:
            self.enzyme_type.choices.append((key, enzyme_descriptions[key]))

    def validate(self):
        if not FlaskForm.validate(self):
            return False
        result = True
        seen = set()

        if self.enzyme_type.data == '-' and self.reaction.data == '-':
            self.enzyme_type.errors.append("Can't both be -")
            self.reaction.errors.append("Can't both be -")
            result = False
        else:
            seen.add(self.enzyme_type.data)
            seen.add(self.reaction.data)
        return result

class ScopeForm(FlaskForm):
    enzyme_type = SelectField('Enzyme type')
    reaction = SelectField('Reaction')
    only_reviewed = BooleanField('Only reviewed data', default=True)
    molecule = SelectField('Reaction', choices=[('product_1_smiles', 'Product'),
                                                ('substrate_1_smiles', 'Substrate 1'),
                                                ('substrate_2_smiles', 'Substrate 2')])

    remove_negative = SelectField('Negative data', choices=[('all', 'Remove all'),
                                                            ('False', 'Do not remove'),
                                                            ('only_negative',
                                                            'Remove molecules where only negative')])
    submit = SubmitField('Generate scope view')

    def set_choices(self):

        self.reaction.choices = [(c, c) for c in ['-'] + (list(Activity.objects().distinct('reaction')))]

        enzyme_types_with_data = list(Activity.objects().distinct('enzyme_type'))
        enzyme_types = EnzymeType.objects(enzyme_type__in=enzyme_types_with_data)

        list_enzyme_types = []
        enzyme_descriptions = {}
        for enz_type in enzyme_types:
            enzyme_descriptions[enz_type.enzyme_type] = f"{enz_type.enzyme_type} - {enz_type.full_name}"
            list_enzyme_types.append(enz_type.enzyme_type)

        list_enzyme_types = sorted(list_enzyme_types)

        self.enzyme_type.choices = [('-', '-')]
        for key in list_enzyme_types:
            self.enzyme_type.choices.append((key, enzyme_descriptions[key]))

    def validate(self):
        if not FlaskForm.validate(self):
            return False
        result = True
        seen = set()

        if self.enzyme_type.data == '-' and self.reaction.data == '-':
            self.enzyme_type.errors.append("Can't both be -")
            self.reaction.errors.append("Can't both be -")
            result = False
        else:
            seen.add(self.enzyme_type.data)
            seen.add(self.reaction.data)
        return result


if __name__ == '__main__':
    import time
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection(host='138.68.135.53', database='test')

    t0 = time.time()
    r_query = Activity.objects().distinct('reaction')
    t1 = time.time()
    print(f"Activity query took {round(t1-t0, 2)} seconds")

    t0 = time.time()
    r_query = Reaction.objects().distinct('name')
    t1 = time.time()
    print(f"Name query took {round(t1-t0, 2)} seconds")