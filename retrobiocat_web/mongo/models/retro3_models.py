from datetime import datetime

import mongoengine as db

from retrobiocat_web.mongo.models.user_models import User

class SearchResult(db.Document):
    datetime = db.DateTimeField(default=datetime.utcnow)
    datetime_completed = db.DateTimeField()
    mcts_config_dict = db.DictField()
    source_mol_config_dict = db.DictField()
    status = db.StringField(default='Queued')
    network_size = db.IntField()
    all_pathways = db.ListField(db.DictField())   # {"reaction_ids": []}
    num_all_pathways = db.IntField()
    solved_pathways = db.ListField(db.DictField())
    num_solved_pathways = db.IntField()
    solved_pathways_with_enzymes = db.ListField(db.DictField())
    num_solved_pathways_solved_pathways_with_enzymes = db.IntField()

class Project(db.Document):
    target_smi = db.StringField()
    project_name = db.StringField(default="")
    expansion_config_dict = db.DictField(default={})
    source_mol_config_dict = db.DictField(default={})
    search_results = db.ListField(db.ReferenceField(SearchResult), default=[], reverse_delete_rule=4)
    last_search_completed = db.DateTimeField()
    num_reactions = db.IntField(default=0)
    saved_pathways = db.ListField(db.ListField(db.StringField()), default=[])
    blocked_reactions = db.DictField(default={})
    owner = db.ReferenceField(User, reverse_delete_rule=4)

    meta = {'indexes': ['owner']}
