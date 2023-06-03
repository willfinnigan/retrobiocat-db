from retrobiocat_web.mongo.models.user_models import User
import mongoengine as db
from retrobiocat_web.mongo.models.reaction_models import Reaction
from datetime import datetime
from retrobiocat_web.mongo.models.comments import Comment


class Issue(db.Document):
    reaction = db.ReferenceField(Reaction)
    issue_reaction_smiles = db.StringField()
    issue_reaction_svg = db.StringField()
    raised_by = db.ReferenceField(User, reverse_delete_rule=2)
    status = db.StringField(default='Open')
    comments = db.ListField(db.ReferenceField(Comment, reverse_delete_rule=4))
    public = db.BooleanField(default=False)
    date = db.DateTimeField(default=datetime.utcnow)

class ReactionSuggestion(db.Document):
    name = db.StringField()
    smarts = db.ListField(db.StringField())
    details = db.StringField()
    owner = db.ReferenceField(User, reverse_delete_rule=2)
    status = db.StringField(default='Open')
    comments = db.ListField(db.ReferenceField(Comment, reverse_delete_rule=4))
    date = db.DateTimeField(default=datetime.utcnow)

class AutoProcessingRule(db.Document):
    multi_step_reaction = db.StringField()
    reactions = db.ListField(db.StringField())
    min_steps = db.IntField()
    max_steps = db.IntField()
    ignore_substrate_two = db.BooleanField()
