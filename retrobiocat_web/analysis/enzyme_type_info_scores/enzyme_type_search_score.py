from retrobiocat_web.mongo.models.biocatdb_models import EnzymeTypeSearchHistory, EnzymeType
from datetime import datetime
import json

score_dict = {'10': 'Thorough',
              '8': 'Good',
              '5': 'Modest',
              '2': 'Basic'}

time_dict = {'10': "in the last month",
             '9': "in the last 6 months",
             '8': "in the last year",
             '6': "in the last 18 months",
             '4': 'in the last 2 years',
             '1': 'more than 2 years ago'}

def get_time_score(days):
    if int(days/30) < 1:
        time_score = 10
    elif int(days/30) < 6:
        time_score = 9
    elif int(days/30) < 12:
        time_score = 8
    elif int(days/30) < 18:
        time_score = 6
    elif int(days/30) < 24:
        time_score = 4
    elif int(days/30) < 48:
        time_score = 1
    else:
        time_score = 0
    return time_score

def calc_score(enzyme_type_obj):
    search_history = EnzymeTypeSearchHistory.objects(enzyme_type=enzyme_type_obj).order_by("date").select_related()

    return_score = 0
    score_string = "No history"
    now = datetime.utcnow().date()
    for history in search_history:
        days_since = (now - history.date.date()).days
        time_score = get_time_score(days_since)
        score = time_score * history.score

        if score > return_score:
            return_score = score
            try:
                score_string = f"'{score_dict[str(history.score)]}' search carried out {time_dict[str(time_score)]} by {history.user.first_name} {history.user.last_name}"
            except:
                score_string = f"'{score_dict[str(history.score)]}' search carried out {time_dict[str(time_score)]} by user"

    return return_score, score_string

def get_history_data(enzyme_type):
    enz_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]
    paper_searches = EnzymeTypeSearchHistory.objects(enzyme_type=enz_type_obj).select_related()

    return_data = []
    for history in paper_searches:
        user = {}
        if history.user is not None:
            user = json.loads(history.user.to_json())
        return_data.append({'user': f"{user.get('first_name', '')} {user.get('last_name', '')}, {user.get('affiliation', '')}",
                            'date': history.date.strftime("%m/%d/%Y"),
                            'score': score_dict[str(history.score)]})

    return return_data


if __name__ == "__main__":
    now = datetime.utcnow().date()
    then = datetime.utcnow().date()
    time_since = (then-now).days
    print(time_since)

