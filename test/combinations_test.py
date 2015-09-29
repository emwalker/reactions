import unittest

from lenrmc.nubase import parse_spec
from lenrmc.combinations import vectors3, regular_combinations, regular_outcomes


def pion_exchange_combinations(total):
    pass


class PossibleDaughtersTest(unittest.TestCase):

    def test_outcomes(self):
        reactants = list(parse_spec('p+7Li'))[0]
        outcomes = list(regular_outcomes(reactants))
        self.assertEqual(42, len(outcomes))
        self.assertEqual(((8, 4),), outcomes[0])

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
        it = regular_combinations((2, 1))
        self.assertEqual([
            ((2, 1),),
            ((1, 0), (1, 1)),
        ], list(it))

    def test_possible_daughters(self):
        ts = list(regular_combinations((6, 3)))
        self.assertEqual(19, len(ts))
        self.assertTrue(all(6 == sum(m for m,a in t) for t in ts))
        self.assertTrue(all(3 == sum(a for m,a in t) for t in ts))


class PionExchangeAndSimultaneousDecayTest(unittest.TestCase):

    def test_simple_case(self):
        possibilities = parse_spec('p+d')
        # c = pion_exchange_and_simultaneous_decay
        # self.assertEqual([],
        #     list(c._outcomes()))
