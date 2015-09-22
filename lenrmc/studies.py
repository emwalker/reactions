import os
import json
from collections import defaultdict


basepath = os.path.dirname(__file__)
DB_PATH = os.path.abspath(os.path.join(basepath, "../db/isotopes.json"))


class Result(object):

    def __init__(self, label, desc, row):
        self.label = label
        self.desc = desc
        self._row = row

    @property
    def json(self):
        row = self._row.copy()
        row['shortDescription'] = self.desc
        return row


class Results(object):

    def __init__(self, results):
        self._results = results

    @property
    def json(self):
        return [r.json for r in self._results]


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
        self._isotopes = defaultdict(list)
        for study in config.get('studies', []):
            desc = study['shortDescription']
            for row in study['isotopes']:
                label = row['label']
                result = Result(label, desc, row)
                self._isotopes[label].append(result)

    def isotopes(self, array):
        results = []
        for label in array:
            for result in self._isotopes[label]:
                results.append(result)
        return Results(results)
