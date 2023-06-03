from flask import render_template, current_app, redirect, url_for
from flask_wtf import FlaskForm
from wtforms import ValidationError, StringField, SubmitField
from wtforms.validators import DataRequired
from wtforms.widgets import TextArea

from retrobiocat_web.app.database_bps.db_analysis import bp
from retrobiocat_web.app.database_bps.db_analysis.functions.task_queues import get_queue, get_job_registeries
from retrobiocat_web.app.database_bps.db_analysis.routes.blast_search.task_do_blast_search import task_do_blast_search
from retrobiocat_web.app.main_site.functions.get_queue_task_details import queue_task_details
from retrobiocat_web.mongo.functions.save_sequence_data.save_sequence_functions import sequence_check, sanitise_sequence
from retrobiocat_web.mongo.model_queries.enzyme_type_queries import all_enzyme_type_strings

QUEUE_NAME = 'db'

def is_protein_seq(form, field):
    protein_seq = sanitise_sequence(field.data)
    bad_chars = sequence_check(protein_seq)

    if len(bad_chars) != 0:
            raise ValidationError(f'non protein characters not allowes - ({bad_chars})')

class BlastForm(FlaskForm):
    sequence = StringField('Target SMILES', widget=TextArea(), validators=[DataRequired(), is_protein_seq])
    submit = SubmitField('Go')


@bp.route('/blast_search', methods=['GET', 'POST'])
def blast_search():
    form = BlastForm()

    if form.validate_on_submit():
        sequence = sanitise_sequence(form.data['sequence'])
        queue = current_app.redis_queues_dict[QUEUE_NAME]
        task = queue.enqueue(task_do_blast_search, sequence)
        task_id = task.get_id()
        return redirect(url_for('db_analysis.blast_search_result', task_id=task_id))

    return render_template('blast_search/blast_search.html', form=form)


@bp.route('/blast_search_result/<task_id>/', methods=['GET'])
def blast_search_result(task_id):
    queue = get_queue(current_app, QUEUE_NAME)
    finished, started = get_job_registeries(queue)

    print(started.get_job_ids())

    # if the job is finished, render the page
    if task_id in list(finished.get_job_ids()):
        task = queue.fetch_job(task_id)
        if task != None:
            table_data = task.result

            seq_table_options = {'table_height': '80vh',
                                 'show_header_filters': True,
                                 'lock_enzyme_types': True}

            enzyme_types = sorted(list(all_enzyme_type_strings()))

            title = 'Results for BLAST search'

            return render_template('show_sequences/show_blast_result_sequences.html',
                                   seq_data=table_data,
                                   seq_table_options=seq_table_options,
                                   enzyme_types=enzyme_types,
                                   title=title)

    # render the loading page
    queue_details, task_details = queue_task_details(task_id, QUEUE_NAME)
    return render_template('queue_loading.html', task_queue=QUEUE_NAME, task_id=task_id,
                           queue_details=queue_details, task_details='',
                           title='Finding homologous sequences', ajax_timer=1500, refresh_timer=30000)