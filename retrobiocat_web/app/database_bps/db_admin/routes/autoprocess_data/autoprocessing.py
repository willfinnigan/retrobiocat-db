from retrobiocat_web.app.database_bps.db_admin import bp
from flask import render_template,  request, jsonify, current_app
from retrobiocat_web.mongo.models.reaction_misc_models import AutoProcessingRule
from flask_security import roles_required
from retrobiocat_web.app.database_bps.db_admin.routes.autoprocess_data.task_autoprocess import task_autoprocess


@bp.route('/autoprocessing_page')
@roles_required('admin')
def autoprocessing_page():
    rules = AutoProcessingRule.objects().as_pymongo()

    # This makes sure everything is a string (ie '_id')
    str_rules = []
    for i, rule in enumerate(rules):
        new_rule = {}
        for key in rule:
            new_rule[key] = str(rule[key])
        str_rules.append(new_rule)

    return render_template('autoprocessing/autoprocessing.html', rules=str_rules)

@bp.route('/_update_autoprocessing_rule', methods=['GET', 'POST'])
@roles_required('admin')
def update_autoprocessing_rule():
    rule_id = request.form['rule_id']

    rule = AutoProcessingRule.objects(id=rule_id)[0]

    ignore_substrate_two = request.form['ignore_substrate_two']
    if ignore_substrate_two == 'false':
        ignore_substrate_two = False
    elif ignore_substrate_two == 'true':
        ignore_substrate_two = True

    reactions = request.form['reactions']
    reactions = reactions.replace("[", "")
    reactions = reactions.replace("]", "")
    reactions = reactions.replace('''"''', "")
    reactions = reactions.replace("'", "")
    reactions = reactions.split(', ')

    rule.multi_step_reaction = request.form['multi_step_reaction']
    rule.reactions = reactions
    rule.min_steps = int(request.form['min_steps'])
    rule.max_steps = int(request.form['max_steps'])
    rule.ignore_substrate_two = ignore_substrate_two

    rule.save()

    result = {'status': 'success',
              'msg': 'Autoprocessing rule updated',
              'issues': []}

    return jsonify(result=result)


@bp.route('/_delete_autoprocessing_rule', methods=['GET', 'POST'])
@roles_required('admin')
def delete_autoprocessing_rule():
    rule_id = request.form['rule_id']
    rule = AutoProcessingRule.objects(id=rule_id)[0]
    rule.delete()

    result = {'status': 'success',
              'msg': 'Autoprocessing rule deleted',
              'issues': []}

    return jsonify(result=result)

@bp.route('/_run_autoprocessing', methods=['GET', 'POST'])
@roles_required('admin')
def run_autoprocessing():
    current_app.db_queue.enqueue(task_autoprocess)

    result = {'status': 'success',
              'msg': 'Autoprocessing job added to database queue. (Do not click again)',
              'issues': []}

    return jsonify(result=result)

@bp.route('/_create_new_autoprocessing_rule', methods=['GET', 'POST'])
@roles_required('admin')
def create_new_autoprocessing_rule():
    new_rule = AutoProcessingRule(multi_step_reaction='New rule',
                                  reactions=[],
                                  min_steps=1, max_steps=1,
                                  ignore_substrate_two=True)
    new_rule.save()

    result = {'status': 'success',
              'msg': 'Autoprocessing job added to database queue. (Do not click again)',
              'issues': []}

    return jsonify(result=result)
