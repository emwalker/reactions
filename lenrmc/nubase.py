from __future__ import absolute_import
import os
import sys
import re
import operator
import itertools
from collections import defaultdict


basepath = os.path.dirname(__file__)
DB_PATH = os.path.abspath(os.path.join(basepath, "../db/nubtab12.asc"))


ALTERNATE_LABELS = {
    '1 n': 'n',
    '1H':  'p',
    '2H':  'd',
    '3H':  't',
}


class RejectCombination(RuntimeError):
    pass


class BadNubaseRow(RuntimeError):
    pass


class HalfLife(object):

    def __init__(self, value, unit):
        self.value = value
        self.unit = unit

    @property
    def seconds(self):
        if 's' == self.unit:
            return float(self.value)
        raise ValueError('do not know how to convert unit: {}'.format(self.unit))

    def __str__(self):
        return '{} {}'.format(self.value, self.unit)


def first_match(pattern, string):
    match = re.search(pattern, string)
    if not match:
        return None
    return match.group()


class Nuclide(object):

    _columns = (
        (  4, 'id'                      ),
        (  7, 'atomicNumber'            ),
        (  9, 'atomicNumberExtra'       ),
        ( 19, 'nuclide'                 ),
        ( 39, 'massExcess'              ),
        ( 61, 'excitationEnergy'        ),
        ( 69, 'halfLife'                ),
        ( 71, 'halfLifeUnit'            ),
        ( 79, 'unknown'                 ),
        ( 93, 'spinAndParity'           ),
        ( 96, 'ensdfArchiveFileYear'    ),
        (105, 'reference'               ),
        (110, 'yearOfDiscovery'         ),
        ( -1, 'decayModesAndIntensities'),
    )

    _not_excited = {
        '1 n',
        '3Li',
        '4Li',
        '4H',
        '5H',
        '5He',
        '5Li',
        '6He',
        '8Be',
    }

    @classmethod
    def load(cls, **kwargs):
        line = kwargs['line']
        row = {}
        endcol_prev = 0
        for endcol, field in cls._columns:
            text = line[endcol_prev:endcol]
            row[field] = text.strip()
            endcol_prev = endcol
        return cls(row)

    def __init__(self, row):
        self._row = row
        self._label = row['nuclide']
        self.atomic_number = int(first_match(r'\d+', self._row['atomicNumber']))
        self.mass_number = int(first_match(r'\d+', self._row['nuclide']))
        g = re.search(r'IS=([\d\.]+)', self._row['decayModesAndIntensities'])
        self.isotopic_abundance = float(g.group(1)) if g else 0.
        self.numbers = (self.mass_number, self.atomic_number)
        if self.is_excited:
            label, self._excitation_level = self._label[:-1], self._label[-1]
            self.label = ALTERNATE_LABELS.get(label, label)
            self.full_label = '{} ({})'.format(self.label, self._excitation_level)
        else:
            label, self._excitation_level = self._label, '0'
            self.label = ALTERNATE_LABELS.get(label, label)
            self.full_label = self.label
        self.signature = (self.label, self._excitation_level)
        kev = first_match(r'[\d\.]+', self._row['massExcess'])
        if not kev:
            raise BadNubaseRow('no mass excess: {}'.format(row))
        self.mass_excess_kev = float(kev)

    @property
    def is_excited(self):
        if self.isotopic_abundance:
            return False
        if self._label in self._not_excited:
            return False
        return any(self._label.endswith(s) for s in 'ijmnpqrx')

    @property
    def half_life(self):
        return HalfLife(self._row['halfLife'], self._row['halfLifeUnit'])

    def json(self):
        return {
            'halfLife':     self.half_life.seconds,
            'atomicNumber': self.atomic_number,
            'massNumber':   self.mass_number,
        }

    def __iter__(self):
        return self.json().iteritems()

    def __eq__(self, o):
        return self.signature == o.signature

    def __hash__(self):
        return hash(self.signature)


class Nuclides(object):

    _nuclides = None

    @classmethod
    def db(cls):
        if cls._nuclides is None:
            cls._nuclides = cls.load(path=DB_PATH)
        return cls._nuclides

    @classmethod
    def load(cls, **kwargs):
        path = kwargs['path']
        nuclides = []
        with open(path) as fh:
            for line in fh:
                try:
                    n = Nuclide.load(line=line)
                    nuclides.append(n)
                except BadNubaseRow:
                    continue
        return cls(nuclides)

    def __init__(self, nuclides):
        self._nuclides = list(nuclides)
        self._by_label = {}
        self._by_signature = {}
        self.isomers = defaultdict(list)
        for n in self._nuclides:
            self._by_label[n._label] = n
            self._by_signature[n.signature] = n
            self.isomers[n.numbers].append(n)

    def get(self, signature):
        return self._by_signature.get(signature)


def vectors3(integer):
    for i in range(integer):
        j = integer - i
        for k in range(j):
            yield (j - k, k, i)


def possible_daughters(totals):
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


class GammaPhoton(object):

    def __init__(self):
        self.mass_number = 0
        self.full_label = self.label = 'ɣ'


class Reaction(object):

    @classmethod
    def load(cls, **kwargs):
        nuclides = Nuclides.db()
        reactants = ((num, nuclides.get(s)) for num, s in kwargs['reactants'])
        daughters = ((num, nuclides.get(s)) for num, s in kwargs['daughters'])
        return cls(reactants, daughters)

    def __init__(self, lvalues, rvalues):
        self._lvalues = list(lvalues)
        self._rvalues = list(rvalues)

    def _sort_key(self, a):
        return a[0].mass_number, a[0].label

    def _fancy_side(self, side):
        isotopes = defaultdict(lambda: 0)
        for c, i in side:
            isotopes[i] += c
        if self._gamma(isotopes):
            isotopes[GammaPhoton()] += 1
        values = [
            '{}×{}'.format(c, i.full_label) if c > 1 else i.full_label
            for i, c
            in sorted(isotopes.items(), key=self._sort_key)
        ]
        return ' + '.join(values)

    def _gamma(self, isotopes):
        if 1 < len(isotopes):
            return False
        return all(c == 1 for c in isotopes.values())

    @property
    def q_kev(self):
        lvalues = sum(num * i.mass_excess_kev for num, i in self._lvalues)
        rvalues = sum(num * i.mass_excess_kev for num, i in self._rvalues)
        return lvalues - rvalues

    @property
    def fancy(self):
        kev = self.q_kev
        sign = '+' if kev >= 0 else '-'
        return '{} → {} {} {:.0f} keV'.format(
            self._fancy_side(self._lvalues),
            self._fancy_side(self._rvalues),
            sign,
            abs(kev),
        )


class Combinations(object):

    @classmethod
    def load(cls, **kwargs):
        nuclides = Nuclides.db()
        reactants = [
            (num, (s, '0')) if isinstance(s, str) else (num, s)
            for num, s in kwargs['reactants']
        ]
        return cls((num, nuclides.get(s)) for num, s in reactants)

    def __init__(self, reactants):
        self._reactants = list(reactants)

    def _outcomes(self):
        numbers = [num * n.numbers for num, n in self._reactants]
        mass_number, atomic_number = tuple(map(operator.add, *numbers))
        return possible_daughters((mass_number, atomic_number))

    def _daughters(self):
        nuclides = Nuclides.db()
        pairs = []
        for _pairs in self._outcomes():
            daughters = [nuclides.isomers[pair] for pair in _pairs]
            if not all(daughters):
                continue
            yield from itertools.product(*daughters)

    def _reactions(self):
        for daughters in self._daughters():
            rvalues = ((1, d) for d in daughters)
            yield Reaction(self._reactants, rvalues)

    def json(self):
        return [r.fancy for r in self._reactions()]
