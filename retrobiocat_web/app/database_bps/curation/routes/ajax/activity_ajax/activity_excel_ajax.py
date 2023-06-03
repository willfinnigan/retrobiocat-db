from flask import request, jsonify
from flask_security import roles_required, current_user
import numpy as np
from werkzeug.utils import secure_filename
import os
import pandas as pd

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.routes.ajax.activity_ajax.errors import not_post_error, \
    error_loading_excel_file, too_long_error
from retrobiocat_web.mongo.model_queries.molecule_queries import get_molecules_in_paper
from retrobiocat_web.mongo.model_queries.paper_queries import paper_from_id


@bp.route('/_upload_activity_excel',methods=['GET', 'POST'])
@roles_required('contributor')
def upload_activity_excel():
    if request.method != 'POST':
        return not_post_error()

    excel_file = request.files['file']
    filename = secure_filename(excel_file.filename)
    paper = paper_from_id(request.form['paper_id_field'])

    if filename[-5:] != '.xlsx':
        return error_loading_excel_file(['Error loading excel file into a dataframe', (str(e))])

    try:
        excel_file.save(filename)
        df = pd.read_excel(filename)
    except Exception as e:
        return error_loading_excel_file(['Error loading excel file into a dataframe', (str(e))])

    # limit number of rows to 2000 if user is not super contributor
    if not current_user.has_role('super_contributor') and len(df.index) > 2000:
        return too_long_error()

    data_list = process_uploaded_excel(df)
    data_list = clear_empty_rows(data_list)
    data_list = remove_enzyme_name_end_spaces(data_list)
    data_list = get_smiles_from_mol_names(data_list, paper)
    os.remove(filename)

    result = {'status': 'success',
              'msg': 'Data loaded from excel - not yet saved..',
              'issues': [],
              'data_list': list(data_list)}
    return jsonify(result=result)


def process_uploaded_excel(df):
    col_rename = {"Reaction": "reaction",
                  "Enzyme name": "enzyme_name",
                  "Substrate 1 SMILES": "substrate_1_smiles",
                  "Substrate 2 SMILES": "substrate_2_smiles",
                  "Product 1 SMILES": "product_1_smiles",
                  "Substrate 1 Name": "substrate_1_name",
                  "Substrate 2 Name": "substrate_2_name",
                  "Product 1 Name": "product_1_name",
                  "Temperature": "temperature",
                  "pH": "ph",
                  "Solvent": "solvent",
                  "Other conditions": "other_conditions",
                  "Notes": "notes",
                  "Reaction volume (ml)": "reaction_vol",
                  "Biocatalyst Formulation": "formulation",
                  "Biocatalyst Concentration (mg/ml)": "biocat_conc",
                  "kcat (min-1)": "kcat",
                  "KM (mM)": "km",
                  "Enz MW (Da)": "mw",
                  "Substrate 1 conc (mM)": "substrate_1_conc",
                  "Substrate 2 conc (mM)": "substrate_2_conc",
                  "Specific activity (U/mg)": "specific_activity",
                  "Conversion (%)": "conversion",
                  "Conversion time (hrs)": "conversion_time",
                  "Selectivity": "selectivity",
                  "Categorical": "categorical",
                  "Binary": "binary"
                }

    df.rename(columns=col_rename, inplace=True)
    cols = [col for col in list(col_rename.values()) if col in list(df.columns)]
    df = df[cols]
    df.replace(np.nan, '', inplace=True)

    data_list = df.to_dict(orient='records')

    return data_list

def does_row_have_data(row):
    for key, value in row.items():
        if value != '' and value is not None:
            return True
    return False

def clear_empty_rows(data_list):
    new_data_list = []
    for row in data_list:
        if does_row_have_data(row):
            new_data_list.append(row)
    return new_data_list

def remove_enzyme_name_end_spaces(data_list):
    for i, data in enumerate(data_list):
        if data.get('enzyme_name', '') != '':
            while data_list[i]['enzyme_name'][-1] == ' ':
                data_list[i]['enzyme_name'] = data_list[i]['enzyme_name'][:-1]
                if data.get('enzyme_name', '') == '':
                    break
    return data_list

def get_smiles_from_mol_names(data_list, paper):

    mols = get_molecules_in_paper(paper)
    mol_dict = {mol.name: mol.smi for mol in mols}

    columns = {'substrate_1_name': 'substrate_1_smiles',
               'substrate_2_name': 'substrate_2_smiles',
               'product_1_name': 'product_1_smiles'}

    for i, data in enumerate(data_list):
        for col in columns:
            name = data.get(col, None)
            if name is None or name not in mol_dict:
                pass
            else:
                smi_col = columns[col]
                data_list[i][smi_col] = mol_dict[name]

    return data_list

