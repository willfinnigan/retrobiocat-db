from flask import jsonify


def too_long_error():
    result = {'status': 'danger',
              'msg': 'Can not load more than 2000 rows at a time',
              'issues': ['Please email your excel to an admin for addition'],
              'data_list': []}
    return jsonify(result=result)

def not_post_error():
    result = {'status': 'danger',
              'msg': 'Error processing file',
              'issues': ['Method is not POST']}
    return jsonify(result=result)

def error_loading_excel_file(issues):
    result = {'status': 'danger',
              'msg': 'Error loading excel file',
              'issues': issues}
    return jsonify(result=result)