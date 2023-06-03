from flask import render_template, flash, redirect, url_for, current_app
from flask_security import roles_required

from retrobiocat_web.app.database_bps.db_admin import bp
from retrobiocat_web.mongo.models.reaction_models import Reaction
from retrobiocat_web.analysis.enzyme_type_info_scores import enzyme_scores

from flask_wtf import FlaskForm
from wtforms import StringField,  TextAreaField, SubmitField
from wtforms.validators import DataRequired, ValidationError, Length
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType


def is_type_taken(form, field):
    for obj in EnzymeType.objects():
        if field.data == obj.enzyme_type:
            raise ValidationError(f'{field.data} is already an enzyme type in the database')

class EnzymeTypeForm(FlaskForm):
    enzyme_type = StringField(validators=[DataRequired(), Length(max=120), is_type_taken])
    full_name = StringField()
    other_abbreviations = StringField()
    description = TextAreaField(validators=[Length(max=255)])
    rep_reaction = StringField()
    submit = SubmitField()


@bp.route('/add_enzyme_type', methods=['GET', 'POST'])
@roles_required('enzyme_types_admin')
def add_enzyme_type():
    form = EnzymeTypeForm()

    if form.validate_on_submit() == True:
        if form.description.data == '':
            form.description.data = None
        if form.other_abbreviations.data == '':
            other_abbreviations = None
        else:
            other_abbreviations = form.other_abbreviations.data.split(', ')

        rep_reaction_name = form.rep_reaction.data
        rep_reaction_query = Reaction.objects(name=rep_reaction_name)
        if len(rep_reaction_query) == 1:
            rep_reaction = rep_reaction_query[0]
        else:
            rep_reaction = None

        enz_type = EnzymeType(enzyme_type=form.enzyme_type.data,
                              full_name=form.full_name.data,
                              other_abbreviations=other_abbreviations,
                              description=form.description.data,
                              rep_reaction=rep_reaction)
        enz_type.save()

        job_id = f"update_enzyme_scores_for_new_enzyme_type_{form.enzyme_type.data}_form.full_name.data"
        current_app.task_queue.enqueue(enzyme_scores.update_scores, job_id=job_id)

        flash("Data added successfully", 'success')
        return redirect(url_for('db_admin.add_enzyme_type'))

    return render_template('enzyme_type/add_enzyme_type.html', form=form)

@bp.route('/edit_enzyme_types', methods=['GET', 'POST'])
@roles_required('enzyme_types_admin')
def edit_enzyme_types():
    headings = ['Name', 'Full name', 'Other abbreviations', 'Description', 'Num rules', 'Rep reaction']
    enzyme_types = EnzymeType.objects().distinct("enzyme_type")
    enzyme_types.sort()

    enz_type_dict_list = EnzymeType.objects().order_by('enzyme_type').as_pymongo()
    renamed_enz_type_dict_list = []

    for enz_type_dict in enz_type_dict_list:
        new_enz_type_dict = {}
        enz_type = enz_type_dict.get('enzyme_type')
        new_enz_type_dict['Name'] = enz_type
        new_enz_type_dict['Full name'] = enz_type_dict.get('full_name', '')
        new_enz_type_dict['Other abbreviations'] = enz_type_dict.get('other_abbreviations', '')
        new_enz_type_dict['Description'] = enz_type_dict.get('description', '')
        new_enz_type_dict['Num rules'] = len(Reaction.objects(enzyme_types=enz_type))
        new_enz_type_dict['Rep reaction'] = enz_type_dict.get('rep_reaction', '')
        renamed_enz_type_dict_list.append(new_enz_type_dict)

    return render_template('enzyme_type/edit_enzyme_types.html',
                           headings=headings, rows=renamed_enz_type_dict_list, enzyme_types=enzyme_types)