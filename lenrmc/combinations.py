import os
import json
import sqlite3
import hashlib
import logging
import pickle
import operator
from os.path import expanduser


LENRMC_DIR = os.path.join(expanduser('~'), '.lenrmc')


class RejectCombination(RuntimeError):
    pass


def vectors3(integer):
    for i in range(integer):
        j = integer - i
        for k in range(j):
            yield (j - k, k, i)


def make_connection():
    try:
        os.mkdir(LENRMC_DIR, 0o755)
    except OSError:
        pass
    conn = sqlite3.connect('{}/lenrmc.db'.format(LENRMC_DIR))
    try:
        conn.execute("""
        create table reactions (parents text, reaction text, q_value_kev real)
        """)
    except sqlite3.OperationalError:
        pass
    return conn


def regular_combinations(totals):
    mass_number, atomic_number = totals
    seen = set()
    for masses in vectors3(mass_number):
        for protons in vectors3(atomic_number):
            pairs = []
            try:
                for i, m in enumerate(masses):
                    p = protons[i]
                    if m < p:
                        raise RejectCombination
                    pair = (m, p)
                    if (0, 0) == pair:
                        continue
                    pairs.append(pair)
            except RejectCombination:
                continue
            pairs = tuple(sorted(pairs))
            if pairs in seen:
                continue
            seen.add(pairs)
            yield pairs


def regular_outcomes(reactants):
    numbers = [num * n.numbers for num, n in reactants]
    mass_number, atomic_number = tuple(map(operator.add, *numbers))
    return regular_combinations((mass_number, atomic_number))




class Combinations(object):

    _connection = None

    @classmethod
    def load(cls, **kwargs):
        reactants = kwargs['reactants']
        del kwargs['reactants']
        outcomes = kwargs.get('outcomes', regular_outcomes)
        return cls(reactants, outcomes, **kwargs)

    def __init__(self, reactants, outcomes, **kwargs):
        self._reactants = list(reactants)
        self._outcomes = outcomes
        self._kwargs = kwargs
        self._lower_bound = float(kwargs['lower_bound']) if 'lower_bound' in kwargs else None
        self.cache_key = self._cache_key()

    @property
    def connection(self):
        if self._connection is None:
            self._connection = make_connection()
        return self._connection

    def _cached_results(self):
        cursor = self.connection.execute(
            "select reaction from reactions where parents = ?",
            (self.cache_key,))
        array = list(cursor)
        if array:
            logging.info('reading previously computed values from cache')
            return (pickle.loads(r[0]) for r in array)
        return None

    def _cache_results(self, results):
        self.connection.executemany("""
        insert into reactions (parents, reaction) values (?, ?)
        """, ((self.cache_key, pickle.dumps(r)) for r in results))
        self.connection.commit()

    def _cache_key(self):
        parents = [(num, n.signature) for num, n in sorted(self._reactants, key=self._sort_key)]
        signature = {'parents': parents}
        for field in ('lower_bound',):
            signature['lower_bound'] = self._kwargs.get('lower_bound')
        string = json.dumps(signature, sort_keys=True).encode('utf-8')
        key = hashlib.sha1(string).hexdigest()
        return key

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self._reactants)

    def _sort_key(self, pair):
        num, nuclide = pair
        return nuclide.signature

    def _daughters(self):
        nuclides = Nuclides.db()
        pairs = []
        for _pairs in self._outcomes():
            daughters = [nuclides.isomers[pair] for pair in _pairs]
            if not all(daughters):
                continue
            yield from itertools.product(*daughters)

    def reactions(self):
        results = self._cached_results()
        if results:
            yield from results
        else:
            results = []
            for daughters in self._daughters():
                rvalues = ((1, d) for d in daughters)
                r = Reaction(self._reactants, rvalues)
                if not self._allowed(r):
                    continue
                yield r
                results.append(r)
            self._cache_results(results)

    def _allowed(self, r):
        return self._lower_bound is None or r.q_value_kev > self._lower_bound
