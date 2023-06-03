import mongoengine as db
from datetime import datetime

class TestCascade(db.Document):
    name = db.StringField()
    doi = db.StringField()
    target_smiles = db.StringField()
    starting_smiles = db.ListField(db.StringField())
    enzymes = db.ListField(db.StringField())
    test_type = db.StringField()
    notes = db.StringField()

class TestCascadeResult(db.Document):
    test_cascade = db.ReferenceField(TestCascade, reverse_delete_rule=db.CASCADE)
    ranking_any = db.IntField()
    ranking_specific = db.IntField()
    time_taken = db.IntField()
    total_pathways = db.IntField()

class TestCascadeRun(db.Document):
    run_name = db.StringField()
    date = db.DateTimeField(default=datetime.utcnow)
    retrobiocat_version = db.StringField()
    reaction_rules_last_updated = db.DateTimeField()
    settings = db.DictField()
    results = db.ListField(db.ReferenceField(TestCascadeResult, reverse_delete_rule=db.CASCADE))
    top_5_score = db.IntField()
    top_25_score = db.IntField()
    top_200_score = db.IntField()




