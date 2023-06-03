from flask import render_template, flash, redirect, url_for
from flask_security import roles_required, current_user

from retrobiocat_web.app.database_bps.tables import bp
from retrobiocat_web.app.database_bps.tables.get_table_data.get_sequence_table_data import \
    get_sequence_table_data_for_enzyme_type, get_sequence_table_data_for_all_enzymes
from retrobiocat_web.mongo.model_queries.enzyme_type_queries import all_enzyme_type_strings, enzyme_type_from_name
from retrobiocat_web.app.app import user_datastore

@bp.route('/edit_sequences', methods=['GET', 'POST'])
@roles_required('super_contributor')
def edit_sequences():
    enzyme_data = get_sequence_table_data_for_all_enzymes()
    enzyme_types = sorted(all_enzyme_type_strings())

    seq_table_options = {'table_height': '80vh',
                         'show_header_filters': True,
                         'lock_enzyme_types': False}

    return render_template('edit_sequences/edit_sequences.html',
                           seq_data=enzyme_data,
                           seq_table_options=seq_table_options,
                           enzyme_types=enzyme_types,
                           title=f"Super contributor access to all sequences")


@bp.route('/enz_champ_seqs/<enzyme_type>', methods=['GET'])
@roles_required('enzyme_champion')
def enzyme_champion_seq(enzyme_type):
    user = user_datastore.get_user(current_user.id)
    enzyme_type_obj = enzyme_type_from_name(enzyme_type)
    if enzyme_type_obj not in user.enzyme_champion:
        flash('No access', 'danger')
        return redirect(url_for('main_site.home'))

    enzyme_data = get_sequence_table_data_for_enzyme_type(enzyme_type)
    enzyme_types = sorted(all_enzyme_type_strings())

    seq_table_options = {'table_height': '90vh',
                         'show_header_filters': True,
                         'lock_enzyme_types': True}

    return render_template('edit_sequences/edit_sequences.html',
                           seq_data=enzyme_data,
                           seq_table_options=seq_table_options,
                           enzyme_types=enzyme_types,
                           title=f"Enzyme champion for {enzyme_type} sequences")