from flask import Flask
from redis import Redis
import rq
from retrobiocat_web.config import Config
from flask_talisman import Talisman
from flask_wtf.csrf import CSRFProtect
from flask_jsglue import JSGlue
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from flask_security import Security, MongoEngineUserDatastore, hash_password, current_user
from flask_mail import Mail
from flask_admin import Admin
from flask_session import Session
from flask_mongoengine import MongoEngine
import datetime
from mongoengine import disconnect
from mongoengine import Q
from datetime import timedelta
import rdkit

csrf = CSRFProtect()
jsglue = JSGlue()
limiter = Limiter(key_func=get_remote_address, default_limits=["100/10second"])
talisman = Talisman(content_security_policy=False)
mail = Mail()
admin_ext = Admin()
db = MongoEngine()
session = Session()

from retrobiocat_web.mongo.models.user_models import User, Role, create_initial_user, user_datastore
from retrobiocat_web.app.user_model_forms import ExtendedConfirmRegisterForm, ExtendedRegisterForm

from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Sequence, Paper, Molecule, Activity, Tag, ActivityIssue, ActivityMol
from retrobiocat_web.mongo.models.reaction_models import Reaction
from retrobiocat_web.mongo.models.reaction_misc_models import Issue, ReactionSuggestion
from retrobiocat_web.mongo.models.comments import Comment
from retrobiocat_web.app.admin import MyAdminIndexView, MyModelView

from retrobiocat_web.app import main_site, retrobiocat
from retrobiocat_web.app.database_bps import curation, tables, adding_papers, db_admin, db_analysis


def create_app(config_class=Config, use_talisman=True):
    print("Create app...")
    app = Flask(__name__, static_folder='static')

    app.config.from_object(config_class)

    print("Init task queues...")
    app.redis = Redis.from_url(app.config['REDIS_URL'])
    app.osra_queue = rq.Queue('osra', connection=app.redis, job_timeout=300, default_timeout=300)
    app.task_queue = rq.Queue('tasks', connection=app.redis, job_timeout=600, default_timeout=600)
    app.network_queue = rq.Queue('network', connection=app.redis, job_timeout=300, default_timeout=300)
    app.pathway_queue = rq.Queue('pathway', connection=app.redis, job_timeout=600, default_timeout=600)
    app.mcts_queue = rq.Queue('mcts', connection=app.redis, job_timeout=1000, default_timeout=1000)
    app.retrorules_queue = rq.Queue('retrorules', connection=app.redis, job_timeout=300, default_timeout=300)
    app.db_queue = rq.Queue('db', connection=app.redis, job_timeout=900, default_timeout=900)
    app.blast_queue = rq.Queue('blast', connection=app.redis, job_timeout=3600, default_timeout=3600)
    app.process_blasts_queue = rq.Queue('process_blasts', connection=app.redis, job_timeout=1200, default_timeout=1200)
    app.alignment_queue = rq.Queue('alignment', connection=app.redis, job_timeout=1800, default_timeout=1800)
    app.preprocess_queue = rq.Queue('preprocess', connection=app.redis, job_timeout=1200, default_timeout=1200)
    app.auto_jobs = rq.Queue('auto_jobs', connection=app.redis, job_timeout=300, default_timeout=600)
    app.embeddings_queue = rq.Queue('embeddings', connection=app.redis, job_timeout=300, default_timeout=600)
    app.heatmap_queue = rq.Queue('heatmap', connection=app.redis, job_timeout=45, result_ttl=1800, default_timeout=45)
    app.scope_queue = rq.Queue('scope', connection=app.redis, job_timeout=100, default_timeout=100)
    app.substrate_summary_queue = rq.Queue('substrate_summary', connection=app.redis, job_timeout=300, default_timeout=300)
    app.redis_queues = [app.task_queue, app.network_queue, app.pathway_queue, app.retrorules_queue,
                        app.db_queue, app.blast_queue, app.alignment_queue, app.process_blasts_queue,
                        app.preprocess_queue, app.auto_jobs, app.osra_queue, app.embeddings_queue,
                        app.heatmap_queue, app.scope_queue, app.substrate_summary_queue]
    app.redis_queues_dict = {'osra': app.osra_queue,
                             'tasks': app.task_queue,
                             'network': app.network_queue,
                             'pathway': app.pathway_queue,
                             'mcts': app.mcts_queue,
                             'retrorules': app.retrorules_queue,
                             'db': app.db_queue,
                             'blast': app.blast_queue,
                             'process_blasts': app.process_blasts_queue,
                             'alignment': app.alignment_queue,
                             'preprocess': app.preprocess_queue,
                             'auto_jobs': app.auto_jobs,
                             'embeddings': app.embeddings_queue,
                             'heatmap': app.heatmap_queue,
                             'scope': app.scope_queue,
                             'substrate_summary': app.substrate_summary_queue
                             }

    print("Init addons...")
    print(rdkit.__version__)
    if use_talisman == True:
        talisman.init_app(app, content_security_policy=False)

    csrf.init_app(app)
    disconnect()
    db.init_app(app)
    jsglue.init_app(app)
    session.init_app(app)
    limiter.init_app(app)
    mail.init_app(app)
    security = Security(app, user_datastore,
                        confirm_register_form=ExtendedConfirmRegisterForm,
                        register_form=ExtendedRegisterForm)

    print("Prepare admin views..")
    admin_ext.init_app(app, index_view=MyAdminIndexView())
    admin_ext.add_view(MyModelView(User))
    admin_ext.add_view(MyModelView(Role))
    admin_ext.add_view(MyModelView(Tag))
    admin_ext.add_view(MyModelView(EnzymeType))
    admin_ext.add_view(MyModelView(Sequence))
    admin_ext.add_view(MyModelView(Paper))
    admin_ext.add_view(MyModelView(Molecule))
    admin_ext.add_view(MyModelView(Activity))
    admin_ext.add_view(MyModelView(Reaction))
    admin_ext.add_view(MyModelView(Issue))
    admin_ext.add_view(MyModelView(Comment))
    admin_ext.add_view(MyModelView(ReactionSuggestion))
    admin_ext.add_view(MyModelView(ActivityIssue))
    admin_ext.add_view(MyModelView(ActivityMol))

    # Create a user to test with
    @app.before_first_request
    def create_user():
        create_initial_user(app)

    @app.context_processor
    def inject_login_mode():
        inject_dict = {}
        inject_dict['login_mode'] = app.config['USE_EMAIL_CONFIRMATION']

        if current_user.is_authenticated:
            user = User.objects(id=current_user.id).select_related()[0]
            if user.has_role('enzyme_teams') and user.enzyme_teams is not None:
                inject_dict['enzyme_teams'] = [enz_type_obj.enzyme_type for enz_type_obj in user.enzyme_teams]
            if user.has_role('enzyme_champion') and user.enzyme_champion is not None:
                inject_dict['enzyme_champion'] = [enz_type_obj.enzyme_type for enz_type_obj in user.enzyme_champion]
            if user.has_role('contributor'):
                inject_dict['user_papers_need_data'] = Paper.objects(Q(owner=user) & (Q(status__icontains='required') | Q(status__icontains='curation'))).count()

            inject_dict['total_team_notifications'] = 0
            inject_dict['team_notifications'] = {}
            inject_dict['champ_seq_notifications'] = {}
            inject_dict['champ_notifications'] = {}

            if 'enzyme_teams' in inject_dict:
                for enz_type in inject_dict['enzyme_teams']:
                    num_papers = Paper.objects(Q(tags=enz_type) & Q(owner=None) & (Q(status__icontains='required') | Q(status__icontains='curation'))).count()
                    inject_dict['team_notifications'][enz_type] = num_papers
                    inject_dict['total_team_notifications'] += num_papers
            if 'enzyme_champion' in inject_dict:
                for enz_type in inject_dict['enzyme_champion']:
                    num_papers = Paper.objects(Q(tags=enz_type) & Q(status__icontains='review')).count()
                    num_seqs = Sequence.objects(Q(enzyme_type=enz_type) & ((Q(sequence=None)|Q(sequence='')) & (Q(sequence_unavailable__ne=True)))).count()
                    inject_dict['champ_notifications'][enz_type] = num_papers
                    inject_dict['champ_seq_notifications'][enz_type] = num_seqs
                    inject_dict['total_team_notifications'] += num_papers + num_seqs

        return inject_dict

    print("Register blueprints...")
    with app.app_context():
        app.register_blueprint(main_site.bp)
        app.register_blueprint(retrobiocat.bp)
        app.register_blueprint(db_analysis.bp)
        app.register_blueprint(curation.bp)
        app.register_blueprint(tables.bp)
        app.register_blueprint(adding_papers.bp)
        app.register_blueprint(db_admin.bp)

        return app









