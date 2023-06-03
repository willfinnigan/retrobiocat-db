from flask import Blueprint

bp = Blueprint('db_admin',
               __name__,
               template_folder='templates',
               static_folder='static',
               static_url_path='/db_admin/static'
               )

from retrobiocat_web.app.database_bps.db_admin.routes.enzyme_type import (add_or_edit_enzyme_type,
                                                                          enzyme_types_ajax)

from retrobiocat_web.app.database_bps.db_admin.routes.db_initialisation import (init_mongodb,
                                                                                misc_admin_functions)

from retrobiocat_web.app.database_bps.db_admin.routes.reaction_rules import edit_rxn_rules


from retrobiocat_web.app.database_bps.db_admin.routes import (auto_jobs_admin,
                                                              download_data,
                                                              merge_identical_enzymes)
from retrobiocat_web.app.database_bps.db_admin.routes.autoprocess_data import autoprocessing


