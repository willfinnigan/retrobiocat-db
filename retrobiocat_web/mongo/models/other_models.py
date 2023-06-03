import mongoengine as db
from datetime import datetime
from retrobiocat_web.mongo.models.comments import Comment
from retrobiocat_web.mongo.models.user_models import User

class AutoJobStatus(db.Document):
    job_type = db.StringField()
    running = db.BooleanField(default=False)
    pause_time = db.IntField(default=30)
    current_status = db.StringField(default='Paused')
    last_ran = db.DateTimeField(default=None)
    next_unassign_update = db.DateTimeField(default=datetime.utcnow)

class Chebi(db.Document):
    code = db.StringField(required=True, unique=True)
    smiles = db.StringField()

class Rhea(db.Document):
    code = db.StringField(required=True, unique=True)
    chebis = db.ListField(db.StringField())
    ecs = db.ListField(db.StringField())
    equation = db.StringField()


class PaperSuggestion(db.Document):
    """Suggestions for papers to be added"""
    doi = db.StringField()
    notes = db.StringField()
    tags = db.StringField()
    owner = db.ReferenceField(User)
    date = db.DateTimeField(default=datetime.utcnow)
    comments = db.ListField(db.ReferenceField(Comment, reverse_delete_rule=4))
    status = db.StringField(default='Open')