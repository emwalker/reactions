"""
Include references to academic studies that are relevant to a given reaction.
"""
import os
import json
from collections import defaultdict


BASEPATH = os.path.dirname(__file__)
DB_PATH = os.path.abspath(os.path.join(BASEPATH, "../db/isotopes.json"))


class Result:
    """Capture one possible result of a given reaction."""

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
        """Convert a result to a JSON-serializable dict."""
        row = self._row.copy()
        row['shortDescription'] = self.description
        return row

    def reference_mark(self, expected):
        """Add an academic reference to a given result."""
        agreement = self._row['change'] == expected
        glyph = '✓' if agreement else '✗'
        return agreement, '{} {} [{}]'.format(glyph, self.nuclide, self.label)


class Results:
    """Hold a collection of possible results for a given reaction."""

    def __init__(self, results):
        self._results = results

    @property
    def json(self):
        """Convert a set of results to a JSON-serializable list."""
        return [r.json for r in self._results]

    def __iter__(self):
        return iter(self._results)


class Studies:
    """Collect together studies that are relevant to a given reaction."""

    _database = None

    @classmethod
    def data(cls):
        """Factory method that returns an initialized database of studies."""
        if cls._database is None:
            with open(DB_PATH) as file:
                config = json.load(file)
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

    def isotopes(self, labels):
        """Return a set of results for a given set of isotope labels."""
        results = []
        for label in labels:
            for result in self._isotopes[label]:
                results.append(result)
        return Results(results)
