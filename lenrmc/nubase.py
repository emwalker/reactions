from __future__ import absolute_import
import os
import sys
import re
import operator
import itertools
from collections import defaultdict


basepath = os.path.dirname(__file__)
NUBASE_PATH = os.path.abspath(os.path.join(basepath, "../db/nubtab12.asc"))

# Fine structure constant
# e^2/(4 * pi * epsilon_0), in MeV fm
FSC = 1.439976

ALTERNATE_LABELS = {
    '1 n':    'n',
    '1H':     'p',
    '2H':     'd',
    '3H':     't',
    '12Cx':   '12C',
    '8Bx':    '8B',
    '10Bx':   '10B',
    '28Px':   '28P',
    '30Px':   '30P',
    '45Vx':   '45V',
    '180Tam': '180Ta',
}

ELEMENTS = {
    'e-':  -1,
    'n':    0,
    'H':    1,
    'He':   2,
    'Li':   3,
    'Be':   4,
    'B':    5,
    'C':    6,
    'N':    7,
    'O':    8,
    'F':    9,
    'Ne':  10,
    'Na':  11,
    'Mg':  12,
    'Al':  13,
    'Si':  14,
    'P':   15,
    'S':   16,
    'Cl':  17,
    'Ar':  18,
    'K':   19,
    'Ca':  20,
    'Sc':  21,
    'Ti':  22,
    'V':   23,
    'Cr':  24,
    'Mn':  25,
    'Fe':  26,
    'Co':  27,
    'Ni':  28,
    'Cu':  29,
    'Zn':  30,
    'Ga':  31,
    'Ge':  32,
    'As':  33,
    'Se':  34,
    'Br':  35,
    'Kr':  36,
    'Rb':  37,
    'Sr':  38,
    'Y':   39,
    'Zr':  40,
    'Nb':  41,
    'Mo':  42,
    'Tc':  43,
    'Ru':  44,
    'Rh':  45,
    'Pd':  46,
    'Ag':  47,
    'Cd':  48,
    'In':  49,
    'Sn':  50,
    'Sb':  51,
    'Te':  52,
    'I':   53,
    'Xe':  54,
    'Cs':  55,
    'Ba':  56,
    'La':  57,
    'Ce':  58,
    'Pr':  59,
    'Nd':  60,
    'Pm':  61,
    'Sm':  62,
    'Eu':  63,
    'Gd':  64,
    'Tb':  65,
    'Dy':  66,
    'Ho':  67,
    'Er':  68,
    'Tm':  69,
    'Yb':  70,
    'Lu':  71,
    'Hf':  72,
    'Ta':  73,
    'W':   74,
    'Re':  75,
    'Os':  76,
    'Ir':  77,
    'Pt':  78,
    'Au':  79,
    'Hg':  80,
    'Tl':  81,
    'Pb':  82,
    'Bi':  83,
    'Po':  84,
    'At':  85,
    'Rn':  86,
    'Fr':  87,
    'Ra':  88,
    'Ac':  89,
    'Th':  90,
    'Pa':  91,
    'U':   92,
    'Np':  93,
    'Pu':  94,
    'Am':  95,
    'Cm':  96,
    'Bk':  97,
    'Cf':  98,
    'Es':  99,
    'Fm':  100,
    'Md':  101,
    'No':  102,
    'Lr':  103,
    'Rf':  104,
    'Db':  105,
    'Sg':  106,
    'Bh':  107,
    'Hs':  108,
}


class Energy(object):

    @classmethod
    def load(cls, **kwargs):
        if 'kev' in kwargs:
            return cls(kwargs['kev'])
        if 'mev' in kwargs:
            return cls(1e3 * kwargs['mev'])
        raise ValueError('do not know how to load energy')

    def __init__(self, energy_kev):
        self.kev = energy_kev

    @property
    def mev(self):
        return 1e-3 * self.kev


class Distance(object):

    @classmethod
    def load(cls, **kwargs):
        return cls(kwargs['fermis'])

    def __init__(self, fermis):
        self.fermis = fermis


class HalfLife(object):

    def __init__(self, value, unit):
        self.value = value
        self.unit = unit

    @property
    def seconds(self):
        if 's' == self.unit:
            return float(self.value)
        if 'd' == self.unit:
            return 86400 * float(self.value)
        if 'y' == self.unit:
            return 3.154e+7 * float(self.value)
        if 'ky' == self.unit:
            return 3.154e+10 * float(self.value)
        if 'Py' == self.unit:
            return 3.154e+23 * float(self.value)
        if 'ms' == self.unit:
            return 0.001 * float(self.value)
        raise ValueError('do not know how to convert unit: {}'.format(self.unit))

    def __str__(self):
        return '{} {}'.format(self.value, self.unit)


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
        self.mass_excess_kev = 0
        self.is_excited = False

    def __repr__(self):
        return 'Electron'


class ElectronNeutrino(object):

    def __init__(self):
        self.mass_number = 0
        self.full_label = self._label = self.label = 'νe'
        self.is_stable = False
        self.spin_and_parity = '1/2+'
        self.atomic_number = 0
        self.signature = ('νe', '0')
        self.numbers = (0, 0)
        self.notes = set() #{'νe'}
        self.mass_excess_kev = 0.00023
        self.is_excited = False

    def __repr__(self):
        return 'ElectronNeutrino'


class BadNubaseRow(RuntimeError):
    pass


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

    def coulomb_barrier(self, o, radius):
        """V(r) for a given radius r in fermis."""
        B = FSC * self.atomic_number * o.atomic_number / float(radius)
        return Energy.load(mev=B)

    def coulomb_barrier_width(self, o, q_value):
        """The width of the Coulomb barrier through which an alpha particle would
        need to escape.
        """
        b = FSC * self.atomic_number * o.atomic_number / q_value.mev
        return Distance.load(fermis=b)


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


def stable_nuclides(nuclides, unstable):
    if unstable:
        return ((1, n) for n in nuclides)
    return ((1, n) for n in nuclides if n.is_stable and not n.is_excited)


def parse_spec(spec, **kwargs):
    unstable = kwargs.get('unstable')
    nuclides = Nuclides.db()
    reactants = []
    for label in (l.strip() for l in spec.split('+')):
        n = nuclides.get((label, '0'))
        if n:
            reactants.append([(1, n)])
        elif 'all' == label:
            row = []
            parent_ub = kwargs.get('parent_ub', 1000)
            for number in ELEMENTS.values():
                if number > parent_ub:
                    continue
                ns = nuclides.atomic_number(number)
                row.extend(stable_nuclides(ns, unstable))
            reactants.append(row)
        else:
            number = ELEMENTS[label]
            ns = nuclides.atomic_number(number)
            it = stable_nuclides(ns, unstable)
            reactants.append(it)
    return itertools.product(*reactants)
