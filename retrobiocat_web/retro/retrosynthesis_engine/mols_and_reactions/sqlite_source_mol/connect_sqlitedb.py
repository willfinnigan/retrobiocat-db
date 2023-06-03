import sqlite3
import sqlalchemy
from pathlib import Path

DEFAULT_DB_PATH = str(Path(__file__).parents[4]) + '/data/buyability/source_mols.db'

class SQLite_Database():

    def __init__(self, path):
        if path is None:
            path = DEFAULT_DB_PATH

        self.engine = None
        self.path = path
        self.conn = self.connect()

    def connect(self):
        return sqlite3.connect(self.path)

    def sqlalchemy_engine(self):
        if self.engine is None:
            conn_string = f"sqlite:////{self.path}"
            self.engine = sqlalchemy.create_engine(conn_string)
        return self.engine

