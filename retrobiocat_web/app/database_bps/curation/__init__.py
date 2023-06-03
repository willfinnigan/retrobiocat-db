from flask import Blueprint

bp = Blueprint('curation',
               __name__,
               template_folder='templates',
               static_folder='static',
               static_url_path='/curation/static'
               )

from retrobiocat_web.app.database_bps.curation.routes.main_tab_routes import (admin_tab,
                                                                              overview_tab,
                                                                              paper_metadata_tab,
                                                                              sequences_tab,
                                                                              molecules_tab,
                                                                              molecules_osr_tab,
                                                                              activity_main_table_tab)

# overview ajax
from retrobiocat_web.app.database_bps.curation.routes.ajax.overview_tab_ajax import (paper_assignment_ajax,
                                                                                     paper_issues_ajax)

# paper metadata ajax
from retrobiocat_web.app.database_bps.curation.routes.ajax.paper_tab_ajax import (crossref_pubmed_lookup_ajax,
                                                                                  update_paper_metadata_ajax)

# sequence curation ajax
from retrobiocat_web.app.database_bps.curation.routes.ajax.sequence_ajax import (sequence_update_ajax,
                                                                                 sequence_excel_upload)

# molecules ajax
from retrobiocat_web.app.database_bps.curation.routes.ajax.molecules_ajax import (new_molecule_ajax,
                                                                                  delete_molecules_ajax,
                                                                                  update_molecule_ajax,
                                                                                  osra_ajax,
                                                                                  molecules_excel_upload_ajax)

# activity ajax
from retrobiocat_web.app.database_bps.curation.routes.ajax.activity_ajax import (activity_excel_ajax,
                                                                                 update_activity_ajax,
                                                                                 apply_reaction_rules_ajax,
                                                                                 apply_reaction_rules_to_table_ajax,
                                                                                 reaction_name_from_smiles_ajax)

# admin ajax
from retrobiocat_web.app.database_bps.curation.routes.ajax.admin_tab_ajax import curation_admin_ajax

# reviewing for all paper, sequences and activity
from retrobiocat_web.app.database_bps.curation.routes.ajax import reviewing_ajax

# misc
from retrobiocat_web.app.database_bps.curation.routes.signup_issues_leaderboard import (leaderboard,
                                                                                        data_issues,
                                                                                        contributor_sign_up,
                                                                                        join_leave_team)
from retrobiocat_web.app.database_bps.curation.routes import mutant_generator
