from flask import jsonify


def no_permission():
    result = {'status': 'danger',
              'msg': 'No permission to create new molecules for this paper',
              'issues': []}
    return jsonify(result=result)

def no_name_found():
    result = {'status': 'danger',
              'msg': 'Could not revolve name -> smiles',
              'issues': []}
    return jsonify(result=result)

def not_post_request():
    result = {'status': 'danger',
              'msg': 'Request is not a post request',
              'issues': []}
    return jsonify(result=result)

def error_processing_mol_images(issues):
    result = {'status': 'danger',
             'msg': 'Problem uploading / processing files',
             'issues': issues}
    return jsonify(result=result)

def error_creating_mol():
    result = {'status': 'danger',
              'msg': 'Error creating new molecule',
              'issues': ['New molecule creation failed'],
              'update_dict': {}}
    return jsonify(result=result)

def mol_not_found():
    result = {'status': 'danger',
              'msg': 'Molecule not found',
              'issues': ['Could not alter molecule as it was not found using specified id'],
              'update_dict': {}}
    return jsonify(result=result)

def general_error(issues):
    result = {'status': 'danger',
              'msg': 'Error',
              'issues': issues,
              'update_dict': {}}
    return jsonify(result=result)

def molecule_name_already_exists(mol_name):
    result = {'status': 'danger',
              'msg': 'Molecule name already exists',
              'issues': [f'{mol_name} is already assigned to a molecule, please pick a different name'],
              'update_dict': {}}
    return jsonify(result=result)

def molecule_name_cannot_be_empty():
    result = {'status': 'danger',
              'msg': 'Molecule name can not be empty',
              'issues': [],
              'update_dict': {}}
    return jsonify(result=result)