from flask import Blueprint

bp = Blueprint('adding_papers',
               __name__,
               template_folder='templates',
               static_folder='static',
               static_url_path='/adding_papers/static'
               )


from retrobiocat_web.app.database_bps.adding_papers.routes import (paper_searching_adding,
                                                                   add_new_paper,
                                                                   suggest_a_paper)