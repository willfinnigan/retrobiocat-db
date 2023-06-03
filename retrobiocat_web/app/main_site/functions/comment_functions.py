from flask_login import current_user

def get_comments_for_a_paper(paper, user):
    """Get comments for for a paper"""
    comments = []
    for comment in paper.comments:
        try:
            owner = comment.owner
        except:
            owner = None

        comment_can_edit = False
        comment_can_delete = False
        if current_user.has_role('admin') or owner == user:
            comment_can_edit = True
            comment_can_delete = True

        if owner is not None:
            comment_by = f"{comment.owner.first_name} {comment.owner.last_name}, {comment.owner.affiliation}"
        else:
            comment_by = f"Unknown user"

        new_comment = {'user': comment_by,
                       'date': comment.date.strftime("%d/%m/%Y, %H:%M:%S"),
                       'comment': comment.text,
                       'comment_id': str(comment.id),
                       'can_edit': comment_can_edit,
                       'can_delete': comment_can_delete
                       }
        comments.append(new_comment)

    return comments