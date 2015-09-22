import os
import json
from collections import defaultdict


basepath = os.path.dirname(__file__)
DB_PATH = os.path.abspath(os.path.join(basepath, "../db/isotopes.json"))


class Result(object):

    def __init__(self, study, row):
        self._study = study
        self._row = row
        self.label = study['label']
        self.nuclide = row['label']
        self.description = study['shortDescription']
        self.reference_line = '[{}] {}'.format(self.label, self.description)
        self._row = row

    @property
    def json(self):
        row = self._row.copy()
        row['shortDescription'] = self.description
        return row

    @property
    def reference_mark(self):
        change = self._row['change']
        if 'increase' == change:
            glyph = '↑'
        elif 'decrease' == change:
            glyph = '↓'
        else:
            raise BadResult('do not recognize change: {}'.format(change))
        return '{} {} [{}]'.format(glyph, self.nuclide, self.label)


class Results(object):

    def __init__(self, results):
        self._results = results

    @property
    def json(self):
        return [r.json for r in self._results]

    def __iter__(self):
        return iter(self._results)


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
            for row in study['isotopes']:
                nuclide = row['label']
                result = Result(study, row)
                self._isotopes[nuclide].append(result)

    def isotopes(self, it):
        results = []
        for label in it:
            for result in self._isotopes[label]:
                results.append(result)
        return Results(results)
