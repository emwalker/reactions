import os
import json


basepath = os.path.dirname(__file__)
DB_PATH = os.path.abspath(os.path.join(basepath, "../db/isotopes.json"))


class Results(object):

    @property
    def json(self):
        return {}


class Studies(object):

    _database = None

    @classmethod
    def db(cls):
        if cls._database is None:
            with open(DB_PATH) as fh:
                config = json.load(fh)
            cls._database = cls(config)
        return cls._database

    def __init__(self, config):
        self._config = config

    def isotopes(self, array):
        return Results()
