import mongoengine as db
import os

def make_default_connection(host='localhost', database='retrobiocat_database'):
    MONGODB_HOST = host
    MONGODB_DB = database
    MONGODB_PORT = os.environ.get('MONGODB_PORT') or 27017
    MONGO_USERNAME = os.environ.get('MONGO_USERNAME') or ''
    MONGO_PASSWORD = os.environ.get('MONGO_PASSWORD') or ''

    db.connect(MONGODB_DB,
               host=MONGODB_HOST,
               port=MONGODB_PORT,
               username=MONGO_USERNAME,
               password=MONGO_PASSWORD,
               alias='default')

    return db

if __name__ == '__main__':
    make_default_connection()

    from retrobiocat_web.mongo.model_queries.specificity_data_query import query_specificity_data
    df = query_specificity_data([], ['CAR'])

