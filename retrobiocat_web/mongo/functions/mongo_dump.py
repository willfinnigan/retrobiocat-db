import zipfile
from pathlib import Path
import subprocess as sp
from flask import current_app
from retrobiocat_web.app.app import db
from retrobiocat_web.mongo.models.user_models import create_initial_user, delete_initial_user

MONGO_FOLDER = str(Path(__file__).parents[1]) + '/mongo_dump'

def execute_mongo_dump():
    print("Executing mongo dump..")
    output_path = f"{MONGO_FOLDER}/mongo_dump.gz"

    command = f'''mongodump --uri="{current_app.config['MONGODB_HOST']}:{current_app.config['MONGODB_PORT']}" --archive={output_path} --verbose'''

    print(f"CMD = {command}")
    sp.run(command, shell=True)
    print('command complete')
    zip_path = f"{MONGO_FOLDER}/mongo.zip"
    zipfile.ZipFile(zip_path, mode='w').write(output_path, arcname='mongo_dump.gz')
    print('zip complete')
    return zip_path

def execute_mongo_restore(filename, nsFrom='retrobiocat_database'):
    nsTo = current_app.config['MONGODB_DB']
    db.connection.drop_database(nsTo)
    db.connection.drop_database('test')
    print("Executing mongo restore..")
    uri = f"--uri='{current_app.config['MONGODB_HOST']}:{current_app.config['MONGODB_PORT']}'"
    include = f"--nsInclude={nsFrom}.*"
    from_and_to = f"--nsFrom='{nsFrom}.*' --nsTo='{nsTo}.*'"
    other_params = f"--archive={filename} --verbose"

    command = f'''mongorestore {uri} {include} {from_and_to} {other_params}'''

    print(f"CMD = {command}")
    sp.run(command, shell=True)

    delete_initial_user(current_app)
    create_initial_user(current_app)

