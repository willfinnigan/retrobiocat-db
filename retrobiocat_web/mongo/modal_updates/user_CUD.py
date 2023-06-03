
def paper_was_unassigned(user):
    try:
        user.total_unassigns += 1
        user.save()
        return []
    except Exception as e:
        return ["Error updating user's unassigns", str(e)]
