import json

from retrobiocat_web.analysis.mol_name_lookup import get_smi_from_name
from retrobiocat_web.app.database_bps.curation import bp
from flask import request, jsonify
from flask_security import roles_required, current_user

from werkzeug.utils import secure_filename
import os
import pandas as pd
import numpy as np

from retrobiocat_web.app.database_bps.curation.functions.clean_str_of_spaces import remove_end_spaces_if_str
from retrobiocat_web.app.database_bps.curation.functions.mol_obj_to_json import get_new_mol_json
from retrobiocat_web.app.database_bps.curation.routes.ajax.error_results import not_a_post_request_error, \
    not_an_excel_file, too_many_excel_rows
from retrobiocat_web.mongo.modal_updates.molecules_CUD import make_new_activity_molecule
from retrobiocat_web.mongo.model_queries.paper_queries import paper_from_id




@bp.route('/_upload_molecules_excel', methods=['GET', 'POST'])
@roles_required('contributor')
def upload_molecules_excel():
    """Route for uploading excel file"""

    if request.method != 'POST':
        return not_a_post_request_error()

    excel_file = request.files['file_mol']
    paper = paper_from_id(request.form['paper_id_field'])
    filename = secure_filename(excel_file.filename)

    if filename[-5:] != '.xlsx':
        return not_an_excel_file()

    excel_file.save(filename)
    df = pd.read_excel(filename)

    # limit number of rows to 200 if user is not super contributor
    max_rows = 200
    if not current_user.has_role('super_contributor') and len(df.index) > max_rows:
        return too_many_excel_rows(max_rows)

    data_list = process_uploaded_moleclues_excel(df)
    os.remove(filename)

    list_new_mols, issues = add_molecules(data_list, paper)

    new_mols_data = []
    for new_mol in list_new_mols:
        new_mols_data.append(get_new_mol_json(new_mol))

    if len(issues) == 0:
        result = {'status': 'success',
                  'msg': 'Molecules saved and added to paper',
                  'issues': [],
                  'new_mols_data': new_mols_data}
    else:
        result = {'status': 'warning',
                  'msg': 'Molecules saved with some issues:',
                  'issues': issues,
                  'new_mols_data': new_mols_data}

    return jsonify(result=result)


def process_uploaded_moleclues_excel(df):
    """Process the uploaded excel into a data list for adding molecules"""

    cols = ['name', 'chem_name', 'smiles']

    cols_to_keep = [c for c in cols if c in list(df.columns)]
    df = df[cols_to_keep]
    df.replace(np.nan, '', inplace=True)

    data_list = df.to_dict(orient='records')

    return data_list


def add_molecules(data_list, paper):
    list_new_mols, issues = [], []
    for i, data in enumerate(data_list):
        smi = None
        chem_name = data.get('chem_name', None)
        chem_name = remove_end_spaces_if_str(chem_name)
        if chem_name == '':
            chem_name = None
        if chem_name is not None:
            chem_name = str(chem_name)
            smi = get_smi_from_name(chem_name)
            if smi is None:
                issues.append(f"Couldn't find smi for {chem_name}")

        if smi is None:
            smi = data.get('smiles', None)
        if smi is not None:
            smi = str(smi)
            smi = remove_end_spaces_if_str(smi)

        name = data.get('name', None)
        if name is not None:
            name = str(name)
            name = remove_end_spaces_if_str(name)

        new_mol = make_new_activity_molecule(smi, paper, name=name, chem_name=chem_name)
        if new_mol is None:
            issues.append(f"Issue adding molecule: {name}, {chem_name}, {smi}")
        else:
            list_new_mols.append(new_mol)

    return list_new_mols, issues


