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

    def test_signature(self):
        n = Nuclide.load(line=self.lines[15])
        self.assertEqual(('6Li', '0'), n.signature)
        n = Nuclide.load(line=self.lines[16])
        self.assertEqual(('6Li', 'i'), n.signature)



class NuclidesTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.nuclides = Nuclides.db()

    def test_nuclide(self):
        self.assertEqual('7Li', self.nuclides.get(('7Li', '0')).label)

    def test_isomers(self):
        n = self.nuclides.get(('7Li', '0'))
        self.assertEqual((7, 3), n.numbers)
        ns = self.nuclides.isomers[n.numbers]
        self.assertEqual([('7Li', '0'), ('7Li', 'i')], [n.signature for n in ns])


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

    @classmethod
    def setUpClass(cls):
        cls.maxDiff = None

    def test_outcomes(self):
        c = Combinations.load(reactants=[(1, '7Li'), (1, 'p')])
        outcomes = list(c._outcomes())
        self.assertEqual(42, len(outcomes))
        self.assertEqual(((8, 4),), outcomes[0])

    def test_reactions(self):
        c = Combinations.load(reactants=[(1, '7Li'), (1, 'p')])
        self.assertEqual(
            ['p + 7Li → 8Be',
             'p + 7Li → 8Be (i)',
             'p + 7Li → 8Be (j)',
             'p + 7Li → n + 7Be',
             'p + 7Li → n + 7Be (i)',
             'p + 7Li → p + 7Li',
             'p + 7Li → p + 7Li (i)',
             'p + 7Li → d + 6Li',
             'p + 7Li → d + 6Li (i)',
             'p + 7Li → t + 5Li',
             'p + 7Li → 3He + 5He',
             'p + 7Li → 3Li + 5H',
             'p + 7Li → 4H + 4Li',
             'p + 7Li → 2×4He',
             'p + 7Li → 2×n + 6Be',
             'p + 7Li → n + p + 6Li',
             'p + 7Li → n + p + 6Li (i)',
             'p + 7Li → 2×p + 6He',
             'p + 7Li → n + d + 5Li',
             'p + 7Li → p + d + 5He',
             'p + 7Li → n + t + 4Li',
             'p + 7Li → n + 3He + 4He',
             'p + 7Li → n + 3Li + 4H',
             'p + 7Li → p + t + 4He',
             'p + 7Li → p + 3He + 4H',
             'p + 7Li → 2×d + 4He',
             'p + 7Li → d + 3He + t']
        , c.json())

    def test_reactions_2(self):
        c = Combinations.load(reactants=[(1, '6Li'), (1, '6Li')])
        self.assertIn('2×6Li → 3×4He', c.json())
