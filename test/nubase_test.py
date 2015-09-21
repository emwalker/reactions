import unittest
import os.path
import sys

from lenrmc.nubase import (
    DB_PATH,
    Nuclide,
    Nuclides,
    Combinations,
    vectors3,
    possible_daughters,
)


class NuclideTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(DB_PATH) as fh:
            cls.lines = list(fh)

    def test_basic_fields(self):
        n = Nuclide.load(line=self.lines[0])
        self.assertEqual(0, n.atomic_number)
        self.assertEqual('613.9 s', str(n.half_life))

    def test_column_widths(self):
        n = Nuclide.load(line=self.lines[16])
        self.assertEqual({
            'atomicNumber':             '003',
            'atomicNumberExtra':        '8W',
            'decayModesAndIntensities': 'IT=100',
            'ensdfArchiveFileYear':     '02',
            'excitationEnergy':         '3562.88    0.10',
            'halfLife':                 '56',
            'halfLifeUnit':             'as',
            'id':                       '006',
            'massExcess':               '17649.76     0.10',
            'nuclide':                  '6Lii',
            'reference':                '81Ro02e',
            'spinAndParity':            '0+      T=1',
            'unknown':                  '14',
            'yearOfDiscovery':          '1981',
        }, n._row)

    def test_mass_number(self):
        n = Nuclide.load(line=self.lines[0])
        self.assertEqual(1, n.mass_number)

    def test_json(self):
        n = Nuclide.load(line=self.lines[0])
        self.assertEqual({
            'atomicNumber': 0,
            'massNumber':   1,
            'halfLife':     613.9,
        }, n.json())

    def test_isotopic_abundance_1(self):
        "There should be an abundance for a stable isotope."
        n = Nuclide.load(line=self.lines[1])
        self.assertEqual(99.9885, n.isotopic_abundance)
        n = Nuclide.load(line=self.lines[2])
        self.assertEqual(0.0115, n.isotopic_abundance)

    def test_isotopic_abundance_2(self):
        "The abundance should be zero for unstable isotopes."
        n = Nuclide.load(line=self.lines[3])
        self.assertEqual(0., n.isotopic_abundance)


class NuclidesTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.nuclides = Nuclides.db()

    def test_nuclide(self):
        self.assertEqual('7Li', self.nuclides.get('7Li').label)

    def test_isomers(self):
        n = self.nuclides.get('7Li')
        self.assertEqual((7, 3), n.numbers)
        ns = self.nuclides.isomers[n.numbers]
        self.assertEqual(['7Li', '7Lii'], [n.label for n in ns])


class ReactionCombinationsTest(unittest.TestCase):

    def test_simple_case(self):
        self.assertEqual([
            (3, 0, 0),
            (2, 1, 0),
            (1, 2, 0),
            (2, 0, 1),
            (1, 1, 1),
            (1, 0, 2),
        ], list(vectors3(3)))

    def test_sums(self):
        sums = [sum(t) for t in vectors3(5)]
        self.assertEqual(15, len(sums))
        self.assertTrue(all(v == 5 for v in sums))

    def test_triples(self):
        it = possible_daughters((2, 1))
        self.assertEqual([
            ((2, 1),),
            ((1, 0), (1, 1)),
        ], list(it))

    def test_possible_daughters(self):
        ts = list(possible_daughters((6, 3)))
        self.assertEqual(19, len(ts))
        self.assertTrue(all(6 == sum(m for m,a in t) for t in ts))
        self.assertTrue(all(3 == sum(a for m,a in t) for t in ts))


class CombinationsTest(unittest.TestCase):

    def test_outcomes(self):
        c = Combinations.load(reactants=[(1, '7Li'), (1, 'p')])
        outcomes = list(c._outcomes())
        self.assertEqual(42, len(outcomes))
        self.assertEqual(((8, 4),), outcomes[0])

    def test_daughters(self):
        c = Combinations.load(reactants=[(1, '7Li'), (1, 'p')])
        self.assertEqual([
            ['8Be'],
            ['8Bei'],
            ['8Bej'],
            ['n', '7Be'],
            ['n', '7Bei'],
            ['p', '7Li'],
            ['p', '7Lii'],
            ['d', '6Li'],
            ['d', '6Lii'],
            ['t', '5Li'],
            ['3He', '5He'],
            ['3Li', '5H'],
            ['4H', '4Li'],
            ['4He', '4He'],
            ['n', 'n', '6Be'],
            ['n', 'p', '6Li'],
            ['n', 'p', '6Lii'],
            ['p', 'p', '6He'],
            ['n', 'd', '5Li'],
            ['p', 'd', '5He'],
            ['n', 't', '4Li'],
            ['n', '3He', '4He'],
            ['n', '3Li', '4H'],
            ['p', 't', '4He'],
            ['p', '3He', '4H'],
            ['d', 'd', '4He'],
            ['d', 't', '3He']]
        , c.json())
