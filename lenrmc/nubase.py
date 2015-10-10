from __future__ import absolute_import
import os
import sys
import re
import operator
import itertools
from collections import defaultdict


basepath = os.path.dirname(__file__)
NUBASE_PATH = os.path.abspath(os.path.join(basepath, "../db/nubtab12.asc"))


ALTERNATE_LABELS = {
    '1 n':  'n',
    '1H':   'p',
    '2H':   'd',
    '3H':   't',
    '12Cx': '12C',
    '10Bx': '10B',
    '30Px': '30P',
}

ELEMENTS = {
    'e-':  -1,
    'n':    0,
    'H':    1,
    'He':   2,
    'Be':   4,
    'Li':   3,
    'C':    6,
    'Al':  13,
    'Mn':  25,
    'Fe':  26,
    'Co':  27,
    'Ni':  28,
    'Cu':  29,
    'Zn':  30,
    'Pd':  46,
}


class Electron(object):

    def __init__(self):
        self.mass_number = 0
        self.full_label = self._label = self.label = 'e-'
        self.is_stable = True
        self.spin_and_parity = '1/2+'
        self.atomic_number = 0
        self.numbers = (0, -1)
        self.signature = (self.label, '0')
        self.notes = {'e-'}
        self.mass_excess_kev = 510.998928


class ElectronNeutrino(object):

    def __init__(self):
        self.mass_number = 0
        self.full_label = self._label = self.label = 'νe'
        self.is_stable = False
        self.spin_and_parity = '1/2+'
        self.atomic_number = 0
        self.signature = ('νe', '0')
        self.numbers = (0, 0)
        self.notes = {'νe'}
        self.mass_excess_kev = 0.00023


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
        (  4, 'massNumber'              ),
        (  7, 'atomicNumber'            ),
        (  9, 'atomicNumberExtra'       ),
        ( 18, 'nuclide'                 ),
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
        '48Ni',
        '49Ni',
        '50Ni',
        '51Ni',
        '52Ni',
        '53Ni',
        '54Ni',
        '55Ni',
        '56Ni',
        '59Ni',
        '60Ni',
        '61Ni',
        '62Ni',
        '63Ni',
        '64Ni',
        '65Ni',
        '66Ni',
        '67Ni',
        '68Ni',
        '69Ni',
        '70Ni',
        '71Ni',
        '72Ni',
        '74Ni',
        '73Ni',
        '75Ni',
        '76Ni',
        '77Ni',
        '78Ni',
        '79Ni',
    }

    @classmethod
    def load(cls, **kwargs):
        line = kwargs['line']
        row = {}
        endcol_prev = 0
        for endcol, field in cls._columns:
            text = line[endcol_prev:endcol].strip()
            if text:
                row[field] = text
            endcol_prev = endcol
        return cls(row)

    def __init__(self, row):
        if not 'massExcess' in row:
            raise BadNubaseRow('no mass excess: {}'.format(row))
        self._row = row
        self._label = row['nuclide']
        self.atomic_number = int(first_match(r'\d+', self._row['atomicNumber']))
        self.mass_number = int(self._row['massNumber'])
        decays = self._row.get('decayModesAndIntensities', '')
        g = re.search(r'IS=([\d\.]+)', decays)
        self.isotopic_abundance = float(g.group(1)) if g else 0.
        self.is_stable = g is not None
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
        kev = first_match(r'[\d\.\-]+', self._row['massExcess'])
        self.mass_excess_kev = float(kev)
        self.spin_and_parity = None
        if 'spinAndParity' in self._row:
            self.spin_and_parity = ' '.join(self._row['spinAndParity'].split())

    _noteworthy = {
        'A':    '→α',
        'B-':   '→β-',
        'B+':   '→β+',
        'B+p':  '→β+p',
        'B+A':  '→β+α',
        'B-n':  '→β-n',
        'B-2n': '→β-2n',
        'B-3n': '→β-3n',
        'B+SF': '→β+SF',
        'B-SF': '→β-SF',
        'B-A':  '→β-α',
        'B-d':  '→β-d',
        'n':    '→n',
        '2n':   '→2n',
        'p':    '→p',
        '2p':   '→2p',
        'EC':   '→ε',
        'IT':   '→IT',
        'SF':   '→SF',
    }

    @property
    def notes(self):
        it = re.split(r'[;=~<]', self._row.get('decayModesAndIntensities', ''))
        return {self._noteworthy.get(token) for token in filter(None, it)} - {None}

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

    def __repr__(self):
        return 'Nuclide({})'.format(self.full_label)


class Nuclides(object):

    _nuclides = None

    @classmethod
    def db(cls):
        if cls._nuclides is None:
            cls._nuclides = cls.load(path=NUBASE_PATH)
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
        nuclides.extend([Electron(), ElectronNeutrino()])
        return cls(nuclides)

    def __init__(self, nuclides):
        self._nuclides = list(nuclides)
        self._by_label = {}
        self._by_signature = {}
        self._by_atomic_number = defaultdict(list)
        self.isomers = defaultdict(list)
        for n in self._nuclides:
            self._by_label[n._label] = n
            self._by_signature[n.signature] = n
            self._by_atomic_number[n.atomic_number].append(n)
            self.isomers[n.numbers].append(n)

    def atomic_number(self, number):
        return self._by_atomic_number[number]

    def get(self, signature):
        return self._by_signature.get(signature)

    def __getitem__(self, signature):
        return self._by_signature[signature]


def parse_spec(spec, **kwargs):
    nuclides = Nuclides.db()
    reactants = []
    for label in (l.strip() for l in spec.split('+')):
        n = nuclides.get((label, '0'))
        if n:
            reactants.append([(1, n)])
        else:
            number = ELEMENTS[label]
            ns = nuclides.atomic_number(number)
            if kwargs.get('unstable'):
                it = ((1, n) for n in ns)
            else:
                it = ((1, n) for n in ns if n.is_stable)
            reactants.append(it)
    return itertools.product(*reactants)
