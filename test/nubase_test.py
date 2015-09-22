import unittest
import os.path
import sys

from lenrmc.nubase import (
    Combinations,
    DB_PATH,
    Nuclide,
    Nuclides,
    possible_daughters,
    Reaction,
    System,
    vectors3,
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
            'massNumber':               '006',
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

    def test_60Ni(self):
        n = Nuclide.load(line=self.lines[781])
        self.assertEqual('60Ni', n.label)
        self.assertEqual(-64472.5, n.mass_excess_kev)

    def test_59Ni(self):
        n = Nuclide.load(line=self.lines[763])
        self.assertEqual('59Ni', n.label)
        self.assertEqual(-61156.1, n.mass_excess_kev)


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


class PossibleDaughtersTest(unittest.TestCase):

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


class ReactionsTest(unittest.TestCase):

    def test_q_value(self):
        r = Reaction.load(
            reactants=[(1, ('p', '0')), (1, ('7Li', '0'))],
            daughters=[(2, ('4He', '0'))],
        )
        self.assertEqual(17346.2443, r.q_value_kev)
        line, refs = r.terminal
        self.assertIn('17346 keV', line)

    def test_4He_note(self):
        r = Reaction.load(
            reactants=[(1, ('6Li', '0')), (1, ('6Li', '0'))],
            daughters=[(3, ('4He', '0'))],
        )
        self.assertIn('4He', r.notes)

    def test_n_note(self):
        r = Reaction.load(
            reactants=[(1, ('p', '0')), (1, ('7Li', '0'))],
            daughters=[(1, ('n', '0')), (1, ('7Be', '0'))],
        )
        self.assertIn('n', r.notes)

    def test_n_transfer_note(self):
        r = Reaction.load(
            reactants=[(1, ('7Li', '0')), (1, ('60Ni', '0'))],
            daughters=[(1, ('6Li', '0')), (1, ('61Ni', '0'))],
        )
        self.assertIn('n-transfer', r.notes)

    def test_stable_note(self):
        r = Reaction.load(
            reactants=[(1, ('7Li', '0')), (1, ('60Ni', '0'))],
            daughters=[(1, ('6Li', '0')), (1, ('61Ni', '0'))],
        )
        self.assertIn('stable', r.notes)

    def test_no_stable_note(self):
        r = Reaction.load(
            reactants=[(1, ('7Li', '0')), (1, ('60Ni', '0'))],
            daughters=[(1, ('8Be', '0')), (1, ('59Co', '0'))],
        )
        self.assertEqual(set(), r.notes)


class SystemTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.maxDiff = None

    def test_outcomes(self):
        c = Combinations.load(reactants=[(1, '7Li'), (1, 'p')])
        outcomes = list(c._outcomes())
        self.assertEqual(42, len(outcomes))
        self.assertEqual(((8, 4),), outcomes[0])

    def test_spins(self):
        s = System.parse('p+7Li', spins=True)
        self.assertEqual(
            ['p + 7Li → 2·4He + 17346 keV                   4He, stable               1/2+, 3/2-           0+, 0+',
             'p + 7Li → ɣ + 8Be + 17254 keV                 ɣ                         1/2+, 3/2-           0+, 1-',
             'p + 7Li → ɣ + 8Be (i) + 628 keV               ɣ                         1/2+, 3/2-           1-, 2+ frg T=1',
             'p + 7Li → 2·d + 4He + -6500 keV               4He, n-transfer, stable   1/2+, 3/2-           0+, 1+, 1+',
             'p + 7Li → d + 6Li + -5027 keV                 n-transfer, stable        1/2+, 3/2-           1+, 1+',
             'p + 7Li → p + 7Li + 0 keV                     stable                    1/2+, 3/2-           1/2+, 3/2-',
             'p + 7Li → d + 6Li (i) + -8589 keV             n-transfer                1/2+, 3/2-           0+ T=1, 1+',
             'p + 7Li → p + d + 5He + -9460 keV             n-transfer                1/2+, 3/2-           1+, 1/2+, 3/2-',
             'p + 7Li → d + 3He + t + -20821 keV            n-transfer                1/2+, 3/2-           1+, 1/2+, 1/2+',
             'p + 7Li → p + t + 4He + -2468 keV             4He                       1/2+, 3/2-           0+, 1/2+, 1/2+',
             'p + 7Li → 3He + 5He + -3966 keV                                         1/2+, 3/2-           1/2+, 3/2-',
             'p + 7Li → t + 5Li + -4434 keV                                           1/2+, 3/2-           1/2+, 3/2-',
             'p + 7Li → 2·p + 6He + -9974 keV                                         1/2+, 3/2-           0+, 1/2+, 1/2+',
             'p + 7Li → p + 7Li (i) + -11243 keV                                      1/2+, 3/2-           1/2+, 3/2- T=3/2',
             'p + 7Li → p + 3He + 4H + -24644 keV                                     1/2+, 3/2-           1/2+, 1/2+, 2-',
             'p + 7Li → 4H + 4Li + -27744 keV                                         1/2+, 3/2-           2-, 2-',
             'p + 7Li → 3Li + 5H + -39364 keV                                         1/2+, 3/2-           (1/2+)',
             'p + 7Li → n + d + 5Li + -10691 keV            n, n-transfer             1/2+, 3/2-           1+, 1/2+, 3/2-',
             'p + 7Li → n + 3He + 4He + -3231 keV           4He, n                    1/2+, 3/2-           0+, 1/2+, 1/2+',
             'p + 7Li → ɣ + 8Be (j) + -10240 keV            ɣ                         1/2+, 3/2-           0+ T=2, 1-',
             'p + 7Li → n + 7Be + -1644 keV                 n                         1/2+, 3/2-           1/2+, 3/2-',
             'p + 7Li → n + p + 6Li + -7251 keV             n                         1/2+, 3/2-           1+, 1/2+, 1/2+',
             'p + 7Li → n + p + 6Li (i) + -10814 keV        n                         1/2+, 3/2-           0+ T=1, 1/2+, 1/2+',
             'p + 7Li → 2·n + 6Be + -12322 keV              n                         1/2+, 3/2-           0+, 1/2+, 1/2+',
             'p + 7Li → n + 7Be (i) + -12625 keV            n                         1/2+, 3/2-           1/2+, 3/2- T=3/2',
             'p + 7Li → n + t + 4Li + -26145 keV            n                         1/2+, 3/2-           1/2+, 1/2+, 2-',
             'p + 7Li → n + 3Li + 4H + -39165 keV           n                         1/2+, 3/2-           1/2+, 2-']
        , s.terminal())

    def test_references(self):
        s = System.parse('p+7Li', references=True)
        self.assertEqual(
            ['p + 7Li → 2·4He + 17346 keV                   4He, stable',
             'p + 7Li → ɣ + 8Be + 17254 keV                 ɣ',
             'p + 7Li → ɣ + 8Be (i) + 628 keV               ɣ',
             'p + 7Li → 2·d + 4He + -6500 keV               4He, n-transfer, stable',
             'p + 7Li → d + 6Li + -5027 keV                 n-transfer, stable          ↑ 6Li [15L]',
             'p + 7Li → p + 7Li + 0 keV                     stable                      ↓ 7Li [15L]',
             'p + 7Li → d + 6Li (i) + -8589 keV             n-transfer                  ↑ 6Li [15L]',
             'p + 7Li → p + d + 5He + -9460 keV             n-transfer',
             'p + 7Li → d + 3He + t + -20821 keV            n-transfer',
             'p + 7Li → p + t + 4He + -2468 keV             4He',
             'p + 7Li → 3He + 5He + -3966 keV',
             'p + 7Li → t + 5Li + -4434 keV',
             'p + 7Li → 2·p + 6He + -9974 keV',
             'p + 7Li → p + 7Li (i) + -11243 keV                                        ↓ 7Li [15L]',
             'p + 7Li → p + 3He + 4H + -24644 keV',
             'p + 7Li → 4H + 4Li + -27744 keV',
             'p + 7Li → 3Li + 5H + -39364 keV',
             'p + 7Li → n + d + 5Li + -10691 keV            n, n-transfer',
             'p + 7Li → n + 3He + 4He + -3231 keV           4He, n',
             'p + 7Li → ɣ + 8Be (j) + -10240 keV            ɣ',
             'p + 7Li → n + 7Be + -1644 keV                 n',
             'p + 7Li → n + p + 6Li + -7251 keV             n                           ↑ 6Li [15L]',
             'p + 7Li → n + p + 6Li (i) + -10814 keV        n                           ↑ 6Li [15L]',
             'p + 7Li → 2·n + 6Be + -12322 keV              n',
             'p + 7Li → n + 7Be (i) + -12625 keV            n',
             'p + 7Li → n + t + 4Li + -26145 keV            n',
             'p + 7Li → n + 3Li + 4H + -39165 keV           n',
             '',
             '[15L] 2015 Lugano E-Cat test by Levi et al.']
        , s.terminal())

    def test_reactions_2(self):
        c = Combinations.load(reactants=[(1, '6Li'), (1, '6Li')])
        self.assertTrue(any('2·6Li → 3·4He + 20899 keV' in r.terminal[0] for r in c.reactions()))
