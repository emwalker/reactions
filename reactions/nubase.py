"""
Provide a Python interface into the Nubase isotope data file.
"""
# pylint: disable=too-many-instance-attributes, too-few-public-methods
from __future__ import absolute_import
import os
import re
import itertools
from collections import defaultdict

from .units import Energy, HalfLife
from .constants import DALTON_KEV


BASEPATH = os.path.dirname(__file__)
NUBASE_PATH = os.path.abspath(os.path.join(BASEPATH, "../db/nubtab12.asc"))


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

TRACE_ISOTOPES = {
    '3H':    ['https://en.wikipedia.org/wiki/Hydrogen'],
    '7Be':   ['https://en.wikipedia.org/wiki/Beryllium'],
    '14C':   ['https://en.wikipedia.org/wiki/Carbon'],
    '18F':   ['https://en.wikipedia.org/wiki/Fluorine'],
    '22Na':  ['https://en.wikipedia.org/wiki/Sodium'],
    '26Al':  ['https://en.wikipedia.org/wiki/Aluminium'],
    '32Si':  ['https://en.wikipedia.org/wiki/Silicon'],
    '36Cl':  ['https://en.wikipedia.org/wiki/Chlorine'],
    '39Ar':  ['https://en.wikipedia.org/wiki/Argon'],
    '41Ca':  ['https://en.wikipedia.org/wiki/Calcium'],
    '53Mn':  ['https://en.wikipedia.org/wiki/Manganese'],
    '59Ni':  ['https://en.wikipedia.org/wiki/Nickel'],
    '79Se':  ['https://en.wikipedia.org/wiki/Selenium'],
    '81Kr':  ['https://en.wikipedia.org/wiki/Krypton'],
    '90Sr':  ['https://en.wikipedia.org/wiki/Strontium'],
    '93Zr':  ['https://en.wikipedia.org/wiki/Zirconium'],
    '99Tc':  ['https://en.wikipedia.org/wiki/Technetium'],
    '107Pd': ['https://en.wikipedia.org/wiki/Palladium'],
    '126Sn': ['https://en.wikipedia.org/wiki/Tin'],
    '129I':  ['https://en.wikipedia.org/wiki/Iodine'],
    '135Cs': ['https://en.wikipedia.org/wiki/Caesium'],
    '137Cs': ['https://en.wikipedia.org/wiki/Caesium'],
    '146Sm': ['https://en.wikipedia.org/wiki/Samarium'],
    '147Pm': ['https://en.wikipedia.org/wiki/Promethium'],
    '182Hf': ['https://en.wikipedia.org/wiki/Caesium'],
    '210Bi': ['https://en.wikipedia.org/wiki/Bismuth'],
    '210Pb': ['https://en.wikipedia.org/wiki/Lead'],
    '210Po': ['https://en.wikipedia.org/wiki/Polonium'],
    '221Fr': ['https://en.wikipedia.org/wiki/Francium'],
    '222Rn': ['https://en.wikipedia.org/wiki/Radon'],
    '223Fr': ['https://en.wikipedia.org/wiki/Francium'],
    '223Ra': ['https://en.wikipedia.org/wiki/Radium'],
    '224Ra': ['https://en.wikipedia.org/wiki/Radium'],
    '225Ac': ['https://en.wikipedia.org/wiki/Actinium'],
    '225Ra': ['https://en.wikipedia.org/wiki/Radium'],
    '226Ra': ['https://en.wikipedia.org/wiki/Radium'],
    '227Th': ['https://en.wikipedia.org/wiki/Thorium'],
    '228Ra': ['https://en.wikipedia.org/wiki/Radium'],
    '228Th': ['https://en.wikipedia.org/wiki/Thorium'],
    '229Th': ['https://en.wikipedia.org/wiki/Thorium'],
    '230Th': ['https://en.wikipedia.org/wiki/Thorium'],
    '231Th': ['https://en.wikipedia.org/wiki/Thorium'],
    '233Pa': ['https://en.wikipedia.org/wiki/Protactinium'],
    '233U':  ['https://en.wikipedia.org/wiki/Uranium'],
    '234Pa': ['https://en.wikipedia.org/wiki/Protactinium'],
    '234Th': ['https://en.wikipedia.org/wiki/Thorium'],
    '236U':  ['https://en.wikipedia.org/wiki/Uranium'],
    '237Np': ['https://en.wikipedia.org/wiki/Neptunium'],
    '238Pu': ['https://en.wikipedia.org/wiki/Plutonium'],
    '239Np': ['https://en.wikipedia.org/wiki/Neptunium'],
    '239Pu': ['https://en.wikipedia.org/wiki/Plutonium'],
    '240Pu': ['https://en.wikipedia.org/wiki/Plutonium'],
    '244Pu': ['https://en.wikipedia.org/wiki/Plutonium'],
}


class Electron:
    """Model an electron."""

    def __init__(self):
        self.mass_number = 0
        self.full_label = self.label = self.initial_label = 'e-'
        self.is_stable = True
        self.spin_and_parity = '1/2+'
        self.atomic_number = 0
        self.numbers = (0, -1)
        self.signature = (self.label, '0')
        self.notes = set()
        self.mass_excess_kev = 0
        self.is_excited = False
        self.in_nature = True
        self.is_trace = False
        self.is_baryon = False

    def __repr__(self):
        return 'Electron'


class ElectronNeutrino:
    """Model an electon neutrino."""

    def __init__(self):
        self.mass_number = 0
        self.full_label = self.label = self.initial_label = 'ν'
        self.is_stable = False
        self.spin_and_parity = '1/2+'
        self.atomic_number = 0
        self.signature = ('ν', '0')
        self.numbers = (0, 0)
        self.notes = set() #{'ν'}
        self.mass_excess_kev = 0.00023
        self.is_excited = False
        self.in_nature = True
        self.is_trace = False
        self.is_baryon = False

    def __repr__(self):
        return 'ElectronNeutrino'


class BadNubaseRow(RuntimeError):
    """Sentinel exception raised when a row in the Nubase file cannot be parsed."""
    pass


def first_match(pattern, string):
    """Return the first match of <pattern>"""
    match = re.search(pattern, string)
    if not match:
        return None
    return match.group()


class Nuclide:
    """Model an individual nuclide that will participate in reactions."""

    _columns = (
        (4, 'massNumber'),
        (7, 'atomicNumber'),
        (9, 'atomicNumberExtra'),
        (18, 'nuclide'),
        (39, 'massExcess'),
        (61, 'excitationEnergy'),
        (69, 'halfLife'),
        (71, 'halfLifeUnit'),
        (79, 'unknown'),
        (93, 'spinAndParity'),
        (96, 'ensdfArchiveFileYear'),
        (105, 'reference'),
        (110, 'yearOfDiscovery'),
        (-1, 'decayModesAndIntensities'),
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
        """Load the nuclide from a line in the Nubase file."""
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
        self.row = row
        self._initialize_basic_fields()
        self._initialize_isotope_fields()
        self._initialize_isomer_fields()
        self._initialize_calculated_fields()

    def _initialize_basic_fields(self):
        self.initial_label = self.row['nuclide']
        self.atomic_number = int(first_match(r'\d+', self.row['atomicNumber']))
        self.is_baryon = True
        self.mass_number = int(self.row['massNumber'])
        self.neutron_number = self.mass_number - self.atomic_number
        self.numbers = (self.mass_number, self.atomic_number)
        self.spin_and_parity = None
        if 'spinAndParity' in self.row:
            self.spin_and_parity = ' '.join(self.row['spinAndParity'].split())

    def _initialize_isotope_fields(self):
        decays = self.row.get('decayModesAndIntensities', '')
        matches = re.search(r'IS=([\d\.]+)', decays)
        self.isotopic_abundance = float(matches.group(1)) if matches else 0.
        self.is_stable = matches is not None
        self.is_trace = self.initial_label in TRACE_ISOTOPES
        self.in_nature = self.is_stable or self.is_trace

    def _initialize_isomer_fields(self):
        if self.is_excited:
            label, self._excitation_level = self.initial_label[:-1], self.initial_label[-1]
            self.label = ALTERNATE_LABELS.get(label, label)
            self.full_label = '{} ({})'.format(self.label, self._excitation_level)
        else:
            label, self._excitation_level = self.initial_label, '0'
            self.label = ALTERNATE_LABELS.get(label, label)
            self.full_label = self.label
        self.signature = (self.label, self._excitation_level)

    def _initialize_calculated_fields(self):
        kev = first_match(r'[\d\.\-]+', self.row['massExcess'])
        self.mass_excess_kev = float(kev)
        self.mass = Energy.load(kev=self.mass_number * DALTON_KEV + self.mass_excess_kev)

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
        """Are there any notes to include with this isotope?"""
        matches = re.split(r'[;=~<]', self.row.get('decayModesAndIntensities', ''))
        notes = {self._noteworthy.get(token) for token in filter(None, matches)} - {None}
        if self.is_trace:
            notes.add('trace')
        return notes

    @property
    def is_excited(self):
        """Is the nuclide an isomer in an excited state?"""
        if self.isotopic_abundance:
            return False
        if self.initial_label in self._not_excited:
            return False
        matches = re.findall('[a-z]', self.initial_label)
        if len(matches) > 1:
            return True
        return False

    @property
    def half_life(self):
        """What is the half-life of this nuclide?"""
        return HalfLife(self.row['halfLife'], self.row['halfLifeUnit'])

    def json(self):
        """Return a JSON-serializable dict."""
        return {
            'halfLife':     self.half_life.seconds,
            'atomicNumber': self.atomic_number,
            'massNumber':   self.mass_number,
        }

    def __iter__(self):
        return iter(self.json().items())

    def __eq__(self, o):
        if isinstance(o, self.__class__):
            return self.signature == o.signature
        return False

    def __hash__(self):
        return hash(self.signature)

    def __repr__(self):
        return 'Nuclide({})'.format(self.full_label)


class Nuclides:
    """Provide a database of nuclides by atomic number, mass number, etc."""

    _nuclides = None

    @classmethod
    def data(cls):
        """Return a memoized singleton of the database of nuclides."""
        if cls._nuclides is None:
            cls._nuclides = cls.load(path=NUBASE_PATH)
        return cls._nuclides

    @classmethod
    def load(cls, **kwargs):
        """Load the database of nuclides from a file."""
        path = kwargs['path']
        nuclides = []
        with open(path) as file:
            for line in file:
                try:
                    nuclide = Nuclide.load(line=line)
                    nuclides.append(nuclide)
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
        self.trace = {}
        self._index_nuclides()

    def _index_nuclides(self):
        for nuclide in self._nuclides:
            self._by_label[nuclide.initial_label] = nuclide
            self._by_signature[nuclide.signature] = nuclide
            self._by_atomic_number[nuclide.atomic_number].append(nuclide)
            self.isomers[nuclide.numbers].append(nuclide)

    def atomic_number(self, atomic_number):
        """What is the nuclide for this number?"""
        return self._by_atomic_number[atomic_number]

    def get(self, signature):
        """Return a nuclide for a given signature."""
        return self._by_signature.get(signature)

    def __iter__(self):
        return iter(self._nuclides)

    def __getitem__(self, signature):
        return self._by_signature[signature]


def stable_nuclides(nuclides, unstable_parents):
    """Return an interateor of (1, nuclide) tuples."""
    if unstable_parents:
        return ((1, n) for n in nuclides)
    return ((1, n) for n in nuclides if n.in_nature and not n.is_excited)


def parse_spec(spec, **kwargs):
    """Parse the reaction spec provided from the command line."""

    unstable_parents = kwargs.get('unstable_parents')
    database = Nuclides.data()
    reactants = []

    for label in (l.strip() for l in spec.split('+')):
        nuclide = database.get((label, '0'))
        if nuclide:
            reactants.append([(1, nuclide)])
            continue

        if label == 'all':
            row = []
            parent_ub = kwargs.get('parent_ub', 1000)
            for number in ELEMENTS.values():
                if number > parent_ub:
                    continue
                nuclides = database.atomic_number(number)
                row.extend(stable_nuclides(nuclides, unstable_parents))
            reactants.append(row)
            continue

        number = ELEMENTS[label]
        nuclides = database.atomic_number(number)
        iterator = stable_nuclides(nuclides, unstable_parents)
        reactants.append(iterator)

    return itertools.product(*reactants)
