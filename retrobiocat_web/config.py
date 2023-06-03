import os
basedir = os.path.abspath(os.path.dirname(__file__))
import redis
from environs import Env

env = Env()

class Config(object):
    SECRET_KEY = env.str('SECRET_KEY', 'testing_key')
    SECURITY_PASSWORD_SALT = env.str("SECURITY_PASSWORD_SALT", '8439842')

    PRODUCTION = env.bool('PRODUCTION', False)

    AUTO_JOBS = env.bool('AUTO_JOBS', False)
    START_AUTO_JOBS_RUNNING = env.bool('START_AUTO_JOBS_RUNNING', False)

    REDIS_URL = env.str('REDIS_URL', 'redis://')
    MONGODB_HOST = env.str('MONGO_HOST', 'mongodb://')
    MONGODB_DB = 'retrobiocat_database'
    MONGODB_PORT = env.str('MONGODB_PORT', '27017')
    MONGODB_CONNECT = False
    MONGODB_SETTINGS = {'db': MONGODB_DB,
                        'host': f'{MONGODB_HOST}:{MONGODB_PORT}/{MONGODB_DB}',
                        'connect': MONGODB_CONNECT}

    ADMIN_EMAIL = env.str('ADMIN_EMAIL', 'admin1@email.com')
    ADMIN_PASSWORD = env.str('ADMIN_PASSWORD', 'password1')
    ADMIN_API_KEY = env.str('ADMIN_API_KEY', 'wifheoifh23jb34pgbp3i4ijbg3434urh34fyb32iuewiubp4g4g3bbovcnwo43fu4t43ibg43ijb34gibj')

    USE_EMAIL_CONFIRMATION = env.bool('USE_EMAIL_CONFIRMATION', False)
    SECURITY_REGISTERABLE = True
    SECURITY_SEND_REGISTER_EMAIL = USE_EMAIL_CONFIRMATION
    SECURITY_CONFIRMABLE = USE_EMAIL_CONFIRMATION
    SECURITY_CHANGEABLE = USE_EMAIL_CONFIRMATION
    SECURITY_RECOVERABLE = USE_EMAIL_CONFIRMATION
    SECURITY_LOGIN_WITHOUT_CONFIRMATION = not USE_EMAIL_CONFIRMATION

    if SECURITY_CONFIRMABLE == True:
        SECURITY_POST_REGISTER_VIEW = '/confirm'
        SECURITY_POST_CHANGE_VIEW = '/change'
    else:
        SECURITY_POST_REGISTER_VIEW = '/'
        SECURITY_POST_CHANGE_VIEW = '/'

    MAIL_SERVER = env.str('MAIL_SERVER', '')
    MAIL_PORT = env.str('MAIL_PORT', '')
    MAIL_USE_SSL = True
    MAIL_USERNAME = env.str('EMAIL_ADDRESS', '')
    MAIL_PASSWORD = env.str('EMAIL_PASSWORD', '')

    CACHE_TYPE = "redis"
    CACHE_DEFAULT_TIMEOUT = 300
    CACHE_REDIS_HOST = REDIS_URL

    SESSION_TYPE = "redis"
    SESSION_REDIS = redis.from_url(REDIS_URL)
    SESSION_USE_SIGNER = True

    # micro service for osra, because its a nightmare to install
    OSRA_API_HOST = env.str('OSRA_API_HOST', 'http://0.0.0.0:8080')

    # optionally used by data entry portal - turns of elements in templates
    OSRA_ENABLED = env.bool('OSRA_ENABLED', True)
    MARVINJS_ENABLED = env.bool('MARVINJS_ENABLED', True)

    # for internal versions, hide the webpages promoting data curation
    HIDE_CURATION = env.bool('HIDE_CURATION', False)

    ''' These variables control whether retrobiocat can access external resources '''
    # Used by network/pathway explorer
    ALLOW_PUBCHEM_LOOKUP = env.bool('ALLOW_PUBCHEM_LOOKUP', True)
    ALLOW_BIOCATHUB_INTEGRATION = env.bool('ALLOW_BIOCATHUB_INTEGRATION', True)

    # used to get identifiers of enzymes if identifier not present in the database
    # used to generate SSN - if this is disabled, SSN's can not be made
    ALLOW_UNIPROT = env.bool('ALLOW_UNIPROT', True)
    ALLOW_SSN_CREATION = env.bool('ALLOW_SSN_CREATION', True)

