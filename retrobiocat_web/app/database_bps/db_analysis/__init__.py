from flask import Blueprint

bp = Blueprint('db_analysis',
               __name__,
               template_folder='templates',
               static_folder='static',
               static_url_path='/db_analysis/static'
               )

from retrobiocat_web.app.database_bps.db_analysis.routes.ssn import ssn, ssn_ajax, download_ajax
from retrobiocat_web.app.database_bps.db_analysis.routes.substrate_summary import substrate_summary_routes, substrate_summary_ajax
from retrobiocat_web.app.database_bps.db_analysis.routes.substrate_similarity_search import substrate_similarity_form
from retrobiocat_web.app.database_bps.db_analysis.routes.substrate_similarity_search import substrate_similarity_route
from retrobiocat_web.app.database_bps.db_analysis.routes.blast_search import blast_search
from retrobiocat_web.app.database_bps.db_analysis.routes.substrate_grid import substrate_grid
from retrobiocat_web.app.database_bps.db_analysis.routes import (bioinformatics,
                                                                 data_api,
                                                                 db_stats,
                                                                 enzyme_homepage,
                                                                 enzyme_toolbox_launch_page,
                                                                 heatmap_route,
                                                                 scope_route,
                                                                 search,
                                                                 substrate_table,
                                                                 summary_page)