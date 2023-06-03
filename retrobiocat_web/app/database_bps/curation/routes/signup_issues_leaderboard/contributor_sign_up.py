from flask import render_template,  redirect, url_for
from flask_security import current_user

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.mongo.models.user_models import Role
from retrobiocat_web.app.app import user_datastore
from flask_wtf import FlaskForm
from wtforms import  SubmitField

class ContributorSignup(FlaskForm):
    submit = SubmitField('Become a contributor')

@bp.route('/contributor_sign_up', methods=['GET', 'POST'])
def contributor_sign_up():

    form = ContributorSignup()

    if form.validate_on_submit() == True:
        if current_user.is_authenticated:
            user = user_datastore.get_user(current_user.id)
            contributor_role = Role.objects(name='contributor')[0]
            user_datastore.add_role_to_user(user, contributor_role)
            return redirect(url_for("curation.contributor_is_signed_up"))

    return render_template('contributor_sign_up/contributor_signup.html', form=form)

@bp.route('/contributor_is_signed_up', methods=['GET'])
def contributor_is_signed_up():
    return render_template('contributor_sign_up/contributor_is_signed_up.html')