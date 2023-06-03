import sys
from rq import Connection, Worker

from retrobiocat_web.app.app import create_app
from environs import Env

env = Env()

if __name__ == '__main__':
    production_mode = env.bool('PRODUCTION', False)

    app = create_app(use_talisman=production_mode)
    app.app_context().push()

    with Connection(app.redis):
        qs = sys.argv[1:] or ['tasks', 'network', 'pathway', 'mcts', 'db', 'process_blasts',
                              'alignment', 'blast', 'preprocess', 'osra', 'embeddings',
                              'heatmap', 'scope', 'substrate_summary', 'auto_jobs']
        w = Worker(qs, log_job_description=False)
        w.work()

