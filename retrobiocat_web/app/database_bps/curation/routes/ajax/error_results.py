from flask import jsonify


def no_access_to_edit_paper():
    result = {'status': 'danger',
              'msg': 'Can not edit paper',
              'issues': ['Users has no access to edit this paper']}
    return jsonify(result=result)


def not_a_post_request_error():
    result = {'status': 'danger',
              'msg': 'Error processing file',
              'issues': ['Method is not POST']}
    return jsonify(result=result)

def not_an_excel_file():
    result = {'status': 'danger',
              'msg': 'Error processing file',
              'issues': ['File does not end in .xlsx']}
    return jsonify(result=result)

def too_many_excel_rows(max_rows):
    result = {'status': 'danger',
              'msg': f'Can not load more than {max_rows} rows',
              'issues': ['Please email your excel to an admin for addition']}
    return jsonify(result=result)
