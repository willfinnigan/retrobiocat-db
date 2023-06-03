from wtforms import StringField
from wtforms.validators import DataRequired

from flask import render_template, flash
from flask_security import roles_required
from werkzeug.utils import secure_filename
from flask_wtf import FlaskForm
from wtforms import SubmitField
from flask_wtf.file import FileField, FileRequired

from retrobiocat_web.app.retrobiocat import bp


class UploadCascadeTestsForm(FlaskForm):
    test_set_name = StringField("Test set name", validators=[DataRequired()])
    excel_file = FileField("Excel file", validators=[FileRequired()])
    submit = SubmitField('Submit')

@roles_required('admin')
@bp.route('/create_cascade_tests', methods=['GET', 'POST'])
def create_cascade_tests():
    form = UploadCascadeTestsForm()
    if form.validate_on_submit():
        filename = secure_filename(form.excel_file.data.filename)
        form.excel_file.data.save(filename)
        if '.xlsx' not in filename:
            flash("Test file must be .xlsx", "danger")
        else:
            issues = load_test_set(form.test_set_name, filename)
            for issue in issues:
                flash(issue, "danger")
            if len(issues) == 0:
                flash('Test set loaded with no problems', 'success')

    return render_template('cascade_testing/upload_cascade_tests.html', form=form)


def load_test_set(name, filename):
    issues = []
    print('filename')


    return issues


