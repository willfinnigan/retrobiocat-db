from flask import Blueprint

bp = Blueprint('tables',
               __name__,
               template_folder='templates',
               static_folder='static',
               static_url_path='/tables/static'
               )


# edit tables
from retrobiocat_web.app.database_bps.tables.routes.edit_tables import (paper_edit_tables,
                                                                        paper_signup_tables,
                                                                        sequence_edit_tables)

# query tables
from retrobiocat_web.app.database_bps.tables.routes.query_tables import (show_papers,
                                                                         show_papers_for_molecule,
                                                                         show_sequences,
                                                                         show_activity)

# ajax
from retrobiocat_web.app.database_bps.tables.routes.ajax import (delete_paper_ajax,
                                                                 paper_load_ajax,
                                                                 self_assign_ajax,
                                                                 sequence_load_ajax,
                                                                 substrate_specificity_ajax)
