import mongoengine as db
from retrobiocat_web.mongo.models.reaction_models import Reaction

default_enzyme_score_dict = {"search_score": 0, "search_score_string": "",
                             "papers_complete_string": "0 out of 0", "papers_complete": 0,
                             "priority_complete": 0, "num_priority": 0,
                             "num_priority_complete": 0, "priority_complete_string": "0 out of 0",
                             "%_proteins": 0, "num_unique_products": 0, "%_unique_products": 0,
                             "% all activity": 0, "num_enzymes": 0, "num_enzymes_w_protein": 0,
                             "num_uniref50": 0, "num_activity": 0, "pc_active": 0, "pc_inactive": 0,
                             "total_score": 0}


class EnzymeType(db.Document):
    enzyme_type = db.StringField(max_length=120, unique=True, required=True)
    full_name = db.StringField(default='')
    description = db.StringField(default='')
    other_abbreviations = db.ListField(db.StringField())
    bioinformatics_status = db.StringField(default='Idle')
    database_score = db.FloatField(default=0)
    database_score_dict = db.DictField(default=default_enzyme_score_dict)
    rep_reaction = db.ReferenceField(Reaction, reverse_delete_rule=1)

    def __unicode__(self):
        return self.enzyme_type

    def __str__(self):
        return self.enzyme_type

    meta = {'indexes': ['enzyme_type']}