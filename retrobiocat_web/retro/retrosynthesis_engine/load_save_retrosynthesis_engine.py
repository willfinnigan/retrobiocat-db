import pickle as pkl
from uuid import uuid4

from pymongo import MongoClient

def save_retro_engine(retro_engine):

    picked_data = pkl.dumps(retro_engine)
    uid = uuid4()

    client = MongoClient()  # add DB url in the constructor if needed
    db = client.test

    # insertion
    db.data.insert_one({
        'uuid': uid,
        'data': picked_data
    })

    # retrieval
    result = db.data.find_one({'uuid': uid})
    assert pkl.loads(result['data']) == data

    return True