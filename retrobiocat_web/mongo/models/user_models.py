import mongoengine as db
from flask_security import UserMixin, RoleMixin, hash_password, MongoEngineUserDatastore
import datetime
import uuid

from retrobiocat_web.mongo.models.enzyme_type import EnzymeType


class Role(db.Document, RoleMixin):
    name = db.StringField(max_length=80, unique=True)
    description = db.StringField(max_length=255)

    def __unicode__(self):
        return self.name

    def __str__(self):
        return self.name

class User(db.Document, UserMixin):
    email = db.StringField(max_length=255, unique=True)
    password = db.StringField(max_length=255)
    active = db.BooleanField(default=True)
    fs_uniquifier = db.StringField(max_length=255)
    confirmed_at = db.DateTimeField()
    roles = db.ListField(db.ReferenceField(Role), default=[])
    first_name = db.StringField(max_length=255)
    last_name = db.StringField(max_length=255)
    affiliation = db.StringField(max_length=255)
    email_opt_in = db.BooleanField()
    enzyme_champion = db.ListField(db.ReferenceField(EnzymeType))
    enzyme_teams = db.ListField(db.ReferenceField(EnzymeType))

    unassigns = db.IntField(default=1)
    total_unassigns = db.IntField(default=0)

    api_key = db.StringField()

    current_login_at = db.DateTimeField()
    current_login_ip = db.StringField()
    last_login_at = db.DateTimeField()
    last_login_ip = db.StringField()
    login_count = db.IntField()

    def generate_new_api_key(self):
        self.api_key = uuid.uuid4()

    def __unicode__(self):
        return f"{self.first_name} {self.last_name}"

    def __str__(self):
        return f"{self.first_name} {self.last_name}"

    @db.queryset_manager
    def contributors(doc_cls, queryset):
        c_role = Role.objects(name='contributor')[0]
        return queryset.filter(roles=c_role)


user_datastore = MongoEngineUserDatastore(db, User, Role)


def delete_initial_user(app):
    user = User.objects(email=app.config['ADMIN_EMAIL']).first()
    if user is not None:
        user.delete()

def create_initial_user(app):
    admin = user_datastore.find_or_create_role('admin', description='admin role')
    enzyme_types_admin = user_datastore.find_or_create_role('enzyme_types_admin', description='enzyme_types_admin')
    contributor = user_datastore.find_or_create_role('contributor', description='contributor')
    super_contributor = user_datastore.find_or_create_role('super_contributor', description='contributor')
    rxn_rules_admin = user_datastore.find_or_create_role('rxn_rules_admin', description='rxn_rules_admin')
    paper_adder = user_datastore.find_or_create_role('paper_adder', description='paper_finder')
    experimental = user_datastore.find_or_create_role('experimental', description='experimental')
    enzyme_champion = user_datastore.find_or_create_role('enzyme_champion', description='enzyme_champion')
    enzyme_teams = user_datastore.find_or_create_role('enzyme_teams', description='enzyme_teams')
    full_data_access = user_datastore.find_or_create_role('full_data_access', description='full_data_access')

    if not user_datastore.get_user(app.config['ADMIN_EMAIL']):
        user = user_datastore.create_user(email=app.config['ADMIN_EMAIL'],
                                          password=hash_password(app.config['ADMIN_PASSWORD']),
                                          api_key=app.config['ADMIN_API_KEY'],
                                          first_name='RetroBioCat',
                                          last_name='Admin',
                                          affiliation='RetroBioCat',
                                          confirmed_at=datetime.datetime.now())
        user_datastore.add_role_to_user(user, admin)
        user_datastore.add_role_to_user(user, enzyme_types_admin)
        user_datastore.add_role_to_user(user, contributor)
        user_datastore.add_role_to_user(user, super_contributor)
        user_datastore.add_role_to_user(user, rxn_rules_admin)
        user_datastore.add_role_to_user(user, paper_adder)
        user_datastore.add_role_to_user(user, experimental)
        user_datastore.add_role_to_user(user, enzyme_champion)
        user_datastore.add_role_to_user(user, enzyme_teams)
        user_datastore.add_role_to_user(user, full_data_access)

