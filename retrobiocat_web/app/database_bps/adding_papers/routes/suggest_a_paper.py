from flask import render_template, jsonify,  request, flash
from flask_security import roles_required, current_user, auth_required
from retrobiocat_web.app.database_bps.adding_papers import bp
from retrobiocat_web.mongo.models.other_models import PaperSuggestion
from retrobiocat_web.mongo.model_queries import paper_queries
from retrobiocat_web.app.app import user_datastore

@bp.route('/suggest_new_paper/', methods=["GET"])
def suggest_new_paper():
    suggestion_id = ""
    doi = ""
    tags = ""
    notes = ""
    comments = []
    owner = ""
    status = ''
    can_save = False
    can_delete = False
    enable_comments = False
    if current_user.is_authenticated:
        can_save = True

    return render_template('paper_suggestions/paper_suggestion.html',
                           suggestion_id=suggestion_id,
                           doi=doi,
                           notes=notes,
                           tags=tags,
                           owner=owner,
                           comments=comments,
                           status=status,
                           can_save=can_save,
                           can_delete=can_delete,
                           enable_comments=enable_comments)


@bp.route('/paper_suggestion/<suggestion_id>')
def paper_suggestion(suggestion_id):
    suggestion = PaperSuggestion.objects(id=suggestion_id).first()

    doi = suggestion.doi
    notes = suggestion.notes
    tags = suggestion.tags
    status = suggestion.status
    owner = ""
    if suggestion.owner is not None:
        owner = f"{suggestion.owner.first_name} {suggestion.owner.last_name}, {suggestion.owner.affiliation}"

    can_save = False
    can_delete = False
    enable_comments = True
    if current_user.is_authenticated:
        user = user_datastore.get_user(current_user.id)
        if suggestion.owner == user or user.has_role('rxn_rules_admin'):
            can_save = True
            can_delete = True
    else:
        user = None

    comments = []
    for comment in suggestion.comments:
        comment_can_edit = False
        comment_can_delete = False
        if current_user.has_role('paper_adder') or comment.owner == user:
            comment_can_edit = True
            comment_can_delete = True

        new_comment = {'user': f"{comment.owner.first_name} {comment.owner.last_name}, {comment.owner.affiliation}",
                       'date': comment.date.strftime("%d/%m/%Y, %H:%M:%S"),
                       'comment': comment.text,
                       'comment_id': str(comment.id),
                       'can_edit': comment_can_edit,
                       'can_delete': comment_can_delete
                       }
        comments.append(new_comment)

    exists = paper_queries.check_if_paper_exists(doi)

    return render_template('paper_suggestions/paper_suggestion.html',
                           suggestion_id=suggestion_id,
                           doi=doi,
                           notes=notes,
                           tags=tags,
                           owner=owner,
                           exists=exists,
                           comments=comments,
                           status=status,
                           can_save=can_save,
                           can_delete=can_delete,
                           enable_comments=enable_comments)



@bp.route('/_save_paper_suggestion', methods=['GET', 'POST'])
@auth_required()
def save_paper_suggestion():
    user = user_datastore.get_user(current_user.id)
    suggestion_id = request.form['suggestion_id']
    doi = request.form['doi'].lower().replace(' ', '')
    notes = str(request.form['notes'])
    tags = str(request.form['tags'])

    if suggestion_id == '':

        if paper_queries.check_if_paper_exists(doi):
            result = {'status': 'danger',
                      'msg': 'DOI is already entered into the database',
                      'issues': ['Suggestion not saved because a paper with this DOI is already in the database'],
                      'suggestion_id': ''}
            return jsonify(result=result)


        suggestion = PaperSuggestion(doi=doi,
                                     tags=tags,
                                     notes=notes,
                                     owner=user)
        suggestion.save()

        result = {'status': 'success',
                  'msg': 'Paper suggestion saved',
                  'issues': [],
                  'suggestion_id': str(suggestion.id)}
        return jsonify(result=result)

    else:
        suggestion = PaperSuggestion.objects(id=suggestion_id).first()
        suggestion.doi = doi
        suggestion.notes = notes
        suggestion.tags = tags
        suggestion.save()

    result = {'status': 'success',
              'msg': 'Paper suggestion saved',
              'issues': [],
              'suggestion_id': str(suggestion.id)}
    return jsonify(result=result)


@bp.route('/_open_close_paper_suggestion', methods=['GET', 'POST'])
@roles_required('paper_adder')
def open_close_paper_suggestion():
    id = request.form['suggestion_id']
    open_close = request.form['open_close']

    suggestion = PaperSuggestion.objects(id=id).first()
    if open_close == 'Open':
        suggestion.status = 'Open'
    elif open_close == 'Close':
        suggestion.status = 'Closed'
    suggestion.save()

    flash(f'{open_close} paper suggestion complete', 'success')
    result = {'status': 'success',
              'msg': f'{open_close} paper suggestion complete',
              'issues': []}
    return jsonify(result=result)


@bp.route('/paper_suggestions_table')
def paper_suggestions_table():
    suggestions = PaperSuggestion.objects().only('id', 'date', 'doi', 'tags', 'notes', 'owner', 'status').select_related()

    paper_suggestions_data = []
    for sug in suggestions:
        data = {}
        data['_id'] = str(sug.id)
        data['date'] = sug.date.isoformat()
        data['doi'] = sug.doi
        data['notes'] = sug.notes
        data['tags'] = sug.tags
        if sug.owner is not None:
            data['owner'] = f"{sug.owner.first_name} {sug.owner.last_name}, {sug.owner.affiliation}"
        else:
            data['owner'] = ""
        data['status'] = sug.status
        paper_suggestions_data.append(data)

    return render_template('paper_suggestions/paper_suggestions_table.html',
                           paper_suggestions_data=paper_suggestions_data)

@bp.route('/_delete_paper_suggestion', methods=['GET', 'POST'])
def delete_paper_suggestion():

    suggestion_id = request.form['suggestion_id']
    sug = PaperSuggestion.objects(id=suggestion_id).first()

    if current_user.is_authenticated:
        user = user_datastore.get_user(current_user.id)
        if sug.owner == user or user.has_role('paper_adder'):
            sug.delete()
            for comment in sug.comments:
                comment.delete()

            flash('Paper suggestion deleted', 'success')
            result = {'status': 'success',
                      'msg': 'Suggestion deleted',
                      'issues': []}
            return jsonify(result=result)

    result = {'status': 'danger',
              'msg': 'No access to delete this suggestion',
              'issues': []}
    return jsonify(result=result)