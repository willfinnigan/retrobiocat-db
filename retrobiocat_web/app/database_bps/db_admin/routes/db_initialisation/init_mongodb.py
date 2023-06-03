
from flask import render_template, flash, redirect, url_for, current_app, jsonify
from flask_security import roles_required
from werkzeug.utils import secure_filename
from flask_wtf import FlaskForm
from wtforms import SubmitField
from flask_wtf.file import FileField

from retrobiocat_web.app.database_bps.db_admin import bp
from retrobiocat_web.mongo.functions import mongo_dump
from retrobiocat_web.app.database_bps.db_admin.functions import download_starting_files, create_retrorules_database
from mongoengine.connection import get_db

class InitDB(FlaskForm):
    mongo_dump = FileField("mongo dump")
    submit = SubmitField('Submit')

@bp.route('/init_db', methods=['GET', 'POST'])
@roles_required('admin')
def init_db():
    form = InitDB()

    db_name = get_db()

    if form.validate_on_submit() == True:
        if form.mongo_dump.data != None:
            filename = secure_filename(form.mongo_dump.data.filename)
            save_path = f"{mongo_dump.MONGO_FOLDER}/{filename}"
            form.mongo_dump.data.save(save_path)
            try:
                if 'mongo_dump.gz' in filename:
                    current_app.db_queue.enqueue(mongo_dump.execute_mongo_restore, save_path)
                    flash("mongo dump job started ok - please check terminal for progress", "success")
                else:
                    print("ERROR File must be mongo_dump.gz")
                    Exception('File must be mongo_dump.gz')
                return redirect(url_for('main_site.home'))

            except Exception as e:
                flash(f"Problem loading mongo dump - {e}", "fail")

    return render_template('init_db/init_db.html', form=form, db_name=db_name)

@bp.route('/_load_aizynthfinder_files', methods=['GET', 'POST'])
@roles_required('admin')
def load_aizynthfinder_files():
    current_app.db_queue.enqueue(download_starting_files.download_aizynthfinder_files)
    result = {'status': 'success',
              'msg': f'Load AIZynthfinder files job queued',
              'issues': []}

    return jsonify(result=result)

@bp.route('/_load_ringbreaker_files', methods=['GET', 'POST'])
@roles_required('admin')
def load_ringbreaker_files():
    current_app.db_queue.enqueue(download_starting_files.download_ringbreaker_files)
    result = {'status': 'success',
              'msg': f'Load AIZynthfinder files job queued',
              'issues': []}

    return jsonify(result=result)

@bp.route('/_load_building_blocks', methods=['GET', 'POST'])
@roles_required('admin')
def load_building_blocks():

    current_app.db_queue.enqueue(download_starting_files.download_building_blocks)
    result = {'status': 'success',
              'msg': f'Load buyable molecules job queued',
              'issues': []}

    return jsonify(result=result)

@bp.route('/_load_retrorules', methods=['GET', 'POST'])
@roles_required('admin')
def load_retrorules():

    current_app.db_queue.enqueue(create_retrorules_database.run)
    result = {'status': 'success',
              'msg': f'Load RetroRules files job queued',
              'issues': []}

    return jsonify(result=result)

@bp.route('/_create_analysis_folders', methods=['GET', 'POST'])
@roles_required('admin')
def create_analysis_folders():
    download_starting_files.create_other_folders()

    result = {'status': 'success',
              'msg': f'Folders created',
              'issues': []}

    return jsonify(result=result)


