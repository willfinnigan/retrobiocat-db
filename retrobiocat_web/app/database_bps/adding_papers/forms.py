from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, BooleanField, DateField
from wtforms.validators import DataRequired,  ValidationError
from flask_wtf import FlaskForm
from wtforms import StringField, BooleanField, SubmitField, SelectField
from wtforms.validators import DataRequired

from retrobiocat_web.mongo.model_queries.paper_queries import check_if_paper_exists
from retrobiocat_web.mongo.models.biocatdb_models import Sequence
from retrobiocat_web.mongo.models.reaction_models import Reaction

def is_doi_taken(form, field):
    doi = field.data
    if check_if_paper_exists(doi):
        raise ValidationError(f'{field.data} is already in the database')

class PaperInfo(FlaskForm):
    short_cit = StringField(validators=[DataRequired()])
    doi = StringField(validators=[DataRequired(), is_doi_taken])
    journal = StringField()
    date = DateField(validators=[DataRequired()])
    title = StringField()
    authors = StringField()
    self_assign = BooleanField(default=False)
    high_importance = BooleanField(default=False)
    tags = StringField()
    submit = SubmitField('Save')

class PapersSearch(FlaskForm):
    enzyme_type = SelectField('Enzyme type')
    enzyme_name = SelectField('Enzyme name')
    reaction = SelectField('Reaction')
    only_reviewed = BooleanField('Only reviewed data', default=True)
    submit = SubmitField('Submit')

    def set_choices(self):
        self.enzyme_type.choices = [(c, c) for c in ['All'] + (list(Sequence.objects().distinct('enzyme_type')))]
        self.enzyme_name.choices = [(c, c) for c in ['All'] + (list(Sequence.objects().distinct('enzyme_name')))]
        self.reaction.choices = [(c, c) for c in ['All'] + (list(Reaction.objects().distinct('name')))]

class PaperSearch(FlaskForm):
    keyword = StringField('Keywords', validators=[DataRequired()])
    enzyme_type = SelectField('Enzyme name')
    submit = SubmitField('Search')

    def set_choices(self):
        self.enzyme_type.choices = [(c, c) for c in (list(Sequence.objects().distinct('enzyme_type')))]