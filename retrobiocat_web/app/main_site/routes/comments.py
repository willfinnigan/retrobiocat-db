
from retrobiocat_web.app.retrobiocat import bp
from flask import jsonify, request,  flash

from flask_security import current_user, auth_required
from retrobiocat_web.app.app import user_datastore
from retrobiocat_web.mongo.models.reaction_misc_models import Issue, ReactionSuggestion
from retrobiocat_web.mongo.models.comments import Comment
from retrobiocat_web.mongo.models.biocatdb_models import Paper, ActivityIssue
from retrobiocat_web.mongo.models.other_models import PaperSuggestion

@bp.route('/_submit_comment', methods=['GET', 'POST'])
@auth_required()
def submit_comment():
    if current_user.is_authenticated:
        user = user_datastore.get_user(current_user.id)
    else:
        user = None

    parent_type = request.form['parent_type']
    parent_id = request.form['parent_id']
    comment_text = request.form['comment']
    comment_id = request.form['comment_id']


    if comment_id != '':
        comment_obj = Comment.objects(id=comment_id)[0]
        if current_user.has_role('rxn_rules_admin') or comment_obj.owner == user:
            comment_obj.text = comment_text
            comment_obj.save()
        else:
            result = {'status': 'danger',
                      'msg': 'Could not edit comment - no access',
                      'issues': []}
            return jsonify(result=result)

    else:
        comment_obj = Comment(owner=user,
                              text=comment_text)
        comment_obj.save()

        if parent_type == 'issue':
            issue = Issue.objects(id=parent_id)[0]
            issue.comments.append(comment_obj)
            issue.save()
        elif parent_type == 'activity_issue':
            issue = ActivityIssue.objects(id=parent_id)[0]
            issue.comments.append(comment_obj)
            issue.save()
        elif parent_type == 'reaction_suggestion':
            r_sug = ReactionSuggestion.objects(id=parent_id)[0]
            r_sug.comments.append(comment_obj)
            r_sug.save()
        elif parent_type == 'paper':
            paper = Paper.objects(id=parent_id)[0]
            paper.comments.append(comment_obj)
            paper.save()
        elif parent_type == 'paper_suggestion':
            suggestion = PaperSuggestion.objects(id=parent_id)[0]
            suggestion .comments.append(comment_obj)
            suggestion .save()

    result = {'status': 'success',
              'msg': 'Comment submitted',
              'issues': []}
    return jsonify(result=result)


@bp.route('/_delete_comment', methods=['GET', 'POST'])
@auth_required()
def delete_comment():
    if current_user.is_authenticated:
        user = user_datastore.get_user(current_user.id)
    else:
        user = None

    comment_id = request.form['comment_id']
    comment_obj = Comment.objects(id=comment_id)[0]
    if comment_obj.owner == user or current_user.has_role('rxn_rules_admin'):
        comment_obj.delete()

        result = {'status': 'success',
                  'msg': 'Comment deleted',
                  'issues': []}
        flash('Comment deleted', 'success')
        return jsonify(result=result)

    else:

        result = {'status': 'danger',
                  'msg': 'Could not delete - You are not the owner of this comment and do not have access',
                  'issues': []}
        flash('You are not the owner of this comment and do not have access', 'danger')
        return jsonify(result=result)


@bp.route('/_load_comment_info', methods=['GET', 'POST'])
def load_comment_info():
    comment_id = request.form['comment_id']

    if comment_id != '':
        comment_obj = Comment.objects(id=comment_id)[0]
        comment_txt = comment_obj.text
    else:
        comment_txt = ''

    result = {'comment_text': comment_txt}
    return jsonify(result=result)

