import os
import json
from collections import defaultdict


basepath = os.path.dirname(__file__)
DB_PATH = os.path.abspath(os.path.join(basepath, "../db/isotopes.json"))


class Result(object):

    def __init__(self, label, citation, row):
        self.label = label
        self.citation = citation
        self._row = row

    @property
    def json(self):
        row = self._row.copy()
        row['citation'] = self.citation
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
            citation = study['citation']
            for row in study['isotopes']:
                label = row['label']
                result = Result(label, citation, row)
                self._isotopes[label].append(result)

    def isotopes(self, array):
        results = []
        for label in array:
            for result in self._isotopes[label]:
                results.append(result)
        return Results(results)
