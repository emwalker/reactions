import os
import json
import sqlite3
import hashlib
import logging
import pickle
import operator
import itertools
from os.path import expanduser

from .nubase import Nuclides


LENRMC_DIR = os.path.join(expanduser('~'), '.lenrmc')


class RejectCombination(RuntimeError):
    pass


class GammaPhoton(object):

    def __init__(self):
        self.mass_number = 0
        self.full_label = self.label = 'ɣ'
        self.is_stable = False
        self.spin_and_parity = '1-'
        self.numbers = (0, 0)
        self.notes = {'ɣ'}


class Reaction(object):

    @classmethod
    def load(cls, **kwargs):
        nuclides = Nuclides.db()
        reactants = ((num, nuclides[s]) for num, s in kwargs['reactants'])
        daughters = ((num, nuclides[s]) for num, s in kwargs['daughters'])
        return cls(reactants, daughters)

    _noteworthy = {
        '4He': 'α',
        'n':   'n',
        't':   't',
    }

    def __init__(self, lvalues, rvalues):
        self._lvalues = list(lvalues)
        self._rvalues = list(rvalues)
        self.q_value_kev = self._q_value_kev()
        self.is_stable = self._is_stable()

    @property
    def notes(self):
        notes = set()
        for _, d in self._rvalues:
            note = self._noteworthy.get(d.label)
            if note:
                notes.add(note)
            for _, p in self._lvalues:
                if self._neutron_transfer(d, p):
                    notes.add('n-transfer')
        if self.is_stable:
            notes.add('stable')
        if self.has_gamma:
            notes.add('ɣ')
            self._rvalues.append((1, GammaPhoton()))
        for num, d in self._rvalues:
            notes |= d.notes
        return notes

    def _neutron_transfer(self, d, p):
        return d.numbers == tuple(map(operator.add, p.numbers, (1, 0)))

    def _is_stable(self):
        return all(d.is_stable for num, d in self._rvalues)

    @property
    def has_gamma(self):
        if 1 < len(self._rvalues):
            return False
        return all(num == 1 for num, c in self._rvalues)

    def _q_value_kev(self):
        lvalues = sum(num * i.mass_excess_kev for num, i in self._lvalues)
        rvalues = sum(num * i.mass_excess_kev for num, i in self._rvalues)
        return lvalues - rvalues


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
            daughters = []
            try:
                for i, m in enumerate(masses):
                    p = protons[i]
                    if m < p:
                        raise RejectCombination
                    pair = (m, p)
                    if (0, 0) == pair:
                        continue
                    daughters.append(pair)
            except RejectCombination:
                continue
            daughters = tuple(sorted(daughters))
            if daughters in seen:
                continue
            seen.add(daughters)
            yield daughters


def add_numbers(*numbers):
    return tuple(map(operator.add, *numbers))


def regular_outcomes(reactants):
    numbers = [num * n.numbers for num, n in reactants]
    mass_number, atomic_number = add_numbers(*numbers)
    return regular_combinations((mass_number, atomic_number))


def normalize(pair):
    o_mass, o_atomic = pair
    if 0 == o_atomic:
        daughters = [(1, 0) for i in range(o_mass)]
    elif o_mass == o_atomic:
        daughters = [(1, 1) for i in range(o_mass)]
    else:
        daughters = [(o_mass, o_atomic)]
    return daughters


class PionExchangeModel(object):

    def __call__(self, reactants):
        num0, smaller = min(reactants, key=lambda t: t[1].mass_number)
        num1, larger  = max(reactants, key=lambda t: t[1].mass_number)
        assert 1 == num0 == num1
        s_mass, s_atomic = smaller.numbers
        if s_mass > 1:
            seen = set()
            p_add = add_numbers(smaller.numbers, (0, 1))
            p_sub = add_numbers(larger.numbers,  (1, 0))
            yield from self._combinations(seen, p_add, p_sub)


class PionExchangeAndDecayModel(PionExchangeModel):

    def _combinations(self, seen, p_add, p_sub):
        outcomes = regular_combinations(p_add)
        for pairs in outcomes:
            daughters = normalize(p_sub)
            for pair in pairs:
                daughters.extend(normalize(pair))
            daughters = tuple(sorted(daughters))
            if daughters in seen:
                continue
            seen.add(daughters)
            yield daughters


class StrictPionExchangeModel(PionExchangeModel):

    def _combinations(self, seen, p_add, p_sub):
        daughters = tuple(sorted(normalize(p_add) + normalize(p_sub)))
        if daughters not in seen:
            yield daughters
        seen.add(daughters)


MODELS = {
    'regular':              regular_outcomes,
    'pion-exchange':        PionExchangeAndDecayModel(),
    'strict-pion-exchange': StrictPionExchangeModel(),
}


class Combinations(object):

    _connection = None

    @classmethod
    def load(cls, **kwargs):
        reactants = kwargs['reactants']
        del kwargs['reactants']
        return cls(reactants, **kwargs)

    def __init__(self, reactants, **kwargs):
        self._reactants = list(reactants)
        self.model_name = kwargs.get('model') or 'regular'
        self._model = MODELS[self.model_name]
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
        signature['model'] = self.model_name
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
        for _daughters in self._model(self._reactants):
            daughters = [nuclides.isomers[pair] for pair in _daughters]
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
