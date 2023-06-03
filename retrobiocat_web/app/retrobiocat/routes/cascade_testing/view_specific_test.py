from retrobiocat_web.app.retrobiocat import bp
from flask import render_template, request
from retrobiocat_web.mongo.models.retro_tests import TestCascade, TestCascadeRun, TestCascadeResult
from flask_security import roles_required, current_user, auth_required


@roles_required('admin')
@bp.route('/specific_cascade_test/', methods=['GET'])
def specific_cascade_test():
    test_id = request.args.get("test_id")
    test_cascade = TestCascade.objects(id=test_id).first().select_related()

    results_table = []
    test_data = {'test_name': test_cascade.name}

    return render_template('cascade_testing/specific_test.html',
                           results_table=results_table,
                           test_data=test_data)

