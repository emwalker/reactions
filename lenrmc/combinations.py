import os
import json
import math
import sqlite3
import hashlib
import logging
import pickle
import operator
import itertools
from os.path import expanduser

from .nubase import Energy, Nuclides, Electron, ElectronNeutrino
from .calculations import (
    IsotopicDecay,
    GamowSuppressionFactor,
    GeigerNuttal,
    Gamow2,
    ReactionEnergy,
)


LENRMC_DIR = os.path.join(expanduser('~'), '.lenrmc')


class RejectCombination(RuntimeError):
    pass


class GammaPhoton(object):

    def __init__(self):
        self.mass_number = 0
        self.is_baryon = False
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
        return cls(reactants, daughters, **kwargs)

    _noteworthy = {
        '4He': 'α',
        'n':   'n',
        't':   't',
    }

    def __init__(self, lvalues, rvalues, **kwargs):
        self._lvalues = list(lvalues)
        self.rvalues = list(rvalues)
        self.q_value = ReactionEnergy(self).value
        self._model = kwargs.get('model')
        self.is_stable = self._is_stable()
        self.any_excited = self._any_excited()
        self.lvalue_delim = self.rvalue_delim = '+'
        # Photons are not counted as a daughter.
        self.daughter_count = sum(n for (n, d) in self.rvalues)
        if self.is_single_body and 1 < len(self._lvalues):
            self.rvalues.append((1, GammaPhoton()))
        g = self.gamow(**kwargs)
        self._gamow = g.value() if g else -1

    @property
    def lvalues(self):
        if 'stimulated-decay' == self._model:
            return [(1, Electron())] + self._lvalues
        return self._lvalues

    @property
    def notes(self):
        notes = set()
        for _, d in self.rvalues:
            note = self._noteworthy.get(d.label)
            if note:
                notes.add(note)
            for _, p in self._lvalues:
                if self._neutron_transfer(d, p):
                    notes.add('n-transfer')
        if self.is_stable:
            notes.add('in nature')
        for num, d in self.rvalues:
            notes |= d.notes
        return notes

    def _neutron_transfer(self, d, p):
        return d.numbers == tuple(map(operator.add, p.numbers, (1, 0)))

    def _is_stable(self):
        return all(d.is_stable for num, d in self.rvalues)

    def _any_excited(self):
        combined = self.rvalues + self._lvalues
        return any(n.is_excited for num, n in combined)

    def _decay_components(self):
        values = [p for num, p in self.rvalues if p.is_baryon]
        if 2 != len(values):
            return None
        n0, n1 = values
        smaller = min(values, key=lambda v: v.mass_number)
        larger = n0 if n1.mass_number == smaller.mass_number else n1
        num, parent = max(self._lvalues, key=lambda t: t[1].atomic_number)
        return parent, (smaller, larger)

    @property
    def has_electron_parent(self):
        return any(isinstance(p, Electron) for num, p in self._lvalues)

    @property
    def is_single_body(self):
        if 1 < len(self.rvalues):
            return False
        return all(num == 1 for num, c in self.rvalues)

    def geiger_nuttal(self):
        return GeigerNuttal.load(self._decay_components(), self.q_value)

    def gamow(self, **kwargs):
        return GamowSuppressionFactor.load(
            self._decay_components(),
            self.q_value,
        )

    def gamow2(self):
        return Gamow2.load(self._decay_components(), self.q_value)

    def decay(self, **kwargs):
        return IsotopicDecay.load(
            self._decay_components(),
            self.q_value,
            **kwargs
        )


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


def normalize(pair):
    o_mass, o_atomic = pair
    if o_mass < 0:
        return [(o_mass, o_atomic)]
    if 0 == o_atomic:
        daughters = [(1, 0) for i in range(o_mass)]
    elif o_mass == o_atomic:
        daughters = [(1, 1) for i in range(o_mass)]
    else:
        daughters = [(o_mass, o_atomic)]
    return daughters


class Model(object):

    def sort_key(self, reaction):
        kev = reaction.q_value.kev
        return kev > 0, kev

    def parents(self, parents, daughters):
        yield parents

    def _smaller_and_larger(self, reactants):
        "This helper method assumes a two-body reaction."
        num0, smaller = min(reactants, key=lambda t: t[1].mass_number)
        num1, larger  = max(reactants, key=lambda t: t[1].mass_number)
        assert 1 == num0 == num1
        return [(num0, smaller), (num1, larger)]


class StandardModel(Model):

    def __call__(self, reactants):
        numbers = [num * n.numbers for num, n in reactants]
        mass_number, atomic_number = add_numbers(*numbers)
        return regular_combinations((mass_number, atomic_number))


class PionExchangeModel(Model):

    # Careful!
    # 7Li  -> 7Be
    # 51Ni -> 52Co
    _transformations = [
        [(0,  1), (0, -1)],
        [(0, -1), (0,  1)],
    ]

    def __call__(self, reactants):
        (num0, smaller), (num1, larger) = self._smaller_and_larger(reactants)
        s_mass, s_atomic = smaller.numbers
        seen = set()
        for adj_left, adj_right in self._transformations:
            p_left = add_numbers(smaller.numbers, adj_left)
            p_right = add_numbers(larger.numbers, adj_right)
            yield from self._combinations(seen, p_left, p_right)


class PionExchangeAndDecayModel(PionExchangeModel):

    def _combinations(self, seen, p_left, p_right):
        outcomes = regular_combinations(p_left)
        for pairs in outcomes:
            daughters = normalize(p_right)
            for pair in pairs:
                daughters.extend(normalize(pair))
            daughters = tuple(sorted(daughters))
            if daughters in seen:
                continue
            seen.add(daughters)
            yield daughters


class StrictPionExchangeModel(PionExchangeModel):

    def _combinations(self, seen, p_left, p_right):
        daughters = tuple(sorted(normalize(p_left) + normalize(p_right)))
        if daughters not in seen:
            yield daughters
        seen.add(daughters)


class ElectronMediatedDecayModel(Model):

    _transformations = [
        # electron capture
        [( 0, -1),  [(0, 0), (0, -1)]],
        # β- decay
        [( 0,  1),  [(0, 0), (0, -1)]],
        # β-β- decay
        [( 0,  2),  [(0, 0), (0, 0), (0, -1), (0, -1)]],
        # ejection of a proton
        [(-1, -1),  [(1, 1)]],
        # ejection of an alpha particle
        [(-4, -2),  [(4, 2), (0, -1)]],
        # [(-8, -4),  [(4, 2), (4, 2), (0, -1)]],
        # [(-12, -6), [(4, 2), (4, 2), (4, 2), (0, -1)]],
    ]

    def parents(self, parents, daughters):
        yield list(parents) + [(1, Electron())]

    def __call__(self, reactants):
        assert 1 == len(reactants)
        num, n = reactants[0]
        assert 1 == num
        for transformation, rvalues2 in self._transformations:
            rvalue = add_numbers(n.numbers, transformation)
            yield tuple(normalize(rvalue) + rvalues2)


class SeparatedNuclideModel(Model):
    """Find those parents who, if separated from a specific
    nuclide, would fission in an exothermic reaction, if fission
    were feasible.
    """
    nuclides = Nuclides.db()

    def parents(self, fission_daughters, daughters):
        assert 1 == len(fission_daughters)
        num, d = fission_daughters[0]
        assert 2 == len(daughters)
        parent_numbers = add_numbers(daughters[0].numbers, daughters[1].numbers)
        for parent in self.nuclides.isomers[parent_numbers]:
            if not parent.in_nature:
                continue
            if parent.mass_number == d.mass_number:
                continue
            yield [(1, parent)] + [(1, Electron())]

    def __call__(self, reactants):
        assert 1 == len(reactants)
        num, n0 = reactants[0]
        assert 1 == num
        for n1 in Nuclides.db():
            yield tuple([n0.numbers, n1.numbers])


class ElectronMediatedFissionModel(Model):

    def sort_key(self, reaction):
        g = reaction._gamow
        kev = reaction.q_value.kev
        return -g, kev > 0, kev

    def __call__(self, reactants):
        assert 1 == len(reactants)
        num, n0 = reactants[0]
        assert 1 == num
        return regular_combinations(n0.numbers)


MODELS = {
    'standard':              StandardModel(),
    'pion-exchange':         PionExchangeAndDecayModel(),
    'strict-pion-exchange':  StrictPionExchangeModel(),
    'induced-decay':         ElectronMediatedDecayModel(),
    'separated-nuclide':     SeparatedNuclideModel(),
    'induced-fission':       ElectronMediatedFissionModel(),
}


class Combinations(object):

    _connection = None

    @classmethod
    def load(cls, **kwargs):
        parents = kwargs['reactants']
        del kwargs['reactants']
        return cls(parents, **kwargs)

    @classmethod
    def connection(cls):
        if cls._connection is None:
            cls._connection = make_connection()
        return cls._connection

    def __init__(self, parents, **kwargs):
        self.model_name = kwargs.get('model') or 'standard'
        self._model = MODELS[self.model_name]
        self._parents = list(parents)
        self._kwargs = kwargs
        self._lower_bound = float(kwargs.get('lower_bound', 0))
        self._upper_bound = float(kwargs.get('upper_bound', 500000))
        self._excited = kwargs.get('excited')
        self.daughter_count = {int(c) for c in kwargs.get('daughter_count', '').split(',') if c}

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self._parents)

    def sort_key(self, reactions):
        return self._model.sort_key(reactions)

    def _reactions(self):
        nuclides = Nuclides.db()
        pairs = []
        for _daughters in self._model(self._parents):
            daughters = [nuclides.isomers[pair] for pair in _daughters]
            if not all(daughters):
                continue
            yield from itertools.product(*daughters)

    def reactions(self):
        for daughters in self._reactions():
            all_parents = self._model.parents(self._parents, daughters)
            for parents in all_parents:
                rvalues = ((1, d) for d in daughters)
                r = Reaction(parents, rvalues, **self._kwargs)
                if not self._allowed(r):
                    continue
                yield r

    def _allowed(self, r):
        conditions = [
            r.q_value.kev >  self._lower_bound,
            r.q_value.kev <= self._upper_bound,
        ]
        if self.daughter_count:
            conditions.append(
                any(c == r.daughter_count for c in self.daughter_count)
            )
        if not self._excited:
            conditions.append(not r.any_excited)
        return all(conditions)
