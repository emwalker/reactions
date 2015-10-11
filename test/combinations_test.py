import unittest

from lenrmc.nubase import parse_spec
from lenrmc.combinations import (
    ElectronStimulatedDecayModel,
    PionExchangeAndDecayModel,
    Reaction,
    regular_combinations,
    regular_outcomes,
    StrictPionExchangeModel,
    vectors3,
)


class ReactionsTest(unittest.TestCase):

    def test_q_value(self):
        r = Reaction.load(
            reactants=[(1, ('p', '0')), (1, ('7Li', '0'))],
            daughters=[(2, ('4He', '0'))],
        )
        self.assertEqual(17346.2443, r.q_value_kev)

    def test_electron_capture_q_value(self):
        r = Reaction.load(
            reactants=[(1, ('e-', '0')), (1, ('63Cu', '0'))],
            daughters=[(1, ('νe', '0')), (1, ('63Ni', '0'))],
        )
        # TODO - fix until it goes to zero?
        self.assertEqual(-67, int(r.q_value_kev))

    def test_4He_note(self):
        r = Reaction.load(
            reactants=[(1, ('6Li', '0')), (1, ('6Li', '0'))],
            daughters=[(3, ('4He', '0'))],
        )
        self.assertIn('α', r.notes)

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
        self.assertIn('in nature', r.notes)

    def test_no_stable_note(self):
        r = Reaction.load(
            reactants=[(1, ('7Li', '0')), (1, ('60Ni', '0'))],
            daughters=[(1, ('8Be', '0')), (1, ('59Co', '0'))],
        )
        self.assertNotIn('in nature', r.notes)

    def test_beta_decay_note(self):
        # Reaction is fictional
        r = Reaction.load(
            reactants=[(1, ('7Li', '0')), (1, ('60Ni', '0'))],
            daughters=[(1, ('t', '0')), (1, ('t', '0'))],
        )
        self.assertEqual({'→β-', 't'}, r.notes)


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

    @classmethod
    def setUpClass(cls):
        cls.model = PionExchangeAndDecayModel()

    def test_p_d(self):
        reactants = list(parse_spec('p+d'))[0]
        self.assertEqual([], list(self.model(reactants)))

    def test_d_d(self):
        reactants = list(parse_spec('d+d'))[0]
        self.assertEqual([
            ((1, 0), (1, 0), (1, 1), (1, 1)),
        ], list(self.model(reactants)))

    def test_d_3He(self):
        reactants = list(parse_spec('d+3He'))[0]
        self.assertEqual([
            ((1, 1), (1, 1), (3, 1)),
        ], list(self.model(reactants)))


class StrictPionExchangeTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.model = StrictPionExchangeModel()

    def test_p_d(self):
        reactants = list(parse_spec('p+d'))[0]
        self.assertEqual(
            [((1, 0), (1, 0), (1, 2)),
             ((1, 0), (1, 1), (1, 1)),
        ], list(self.model(reactants)))

    def test_d_d(self):
        reactants = list(parse_spec('d+d'))[0]
        self.assertEqual([
            ((1, 0), (1, 0), (1, 1), (1, 1)),
        ], list(self.model(reactants)))

    def test_d_3He(self):
        reactants = list(parse_spec('d+3He'))[0]
        self.assertEqual([
            ((1, 1), (1, 1), (3, 1)),
            ((1, 0), (1, 0), (1, 1), (1, 1), (1, 1)),
        ], list(self.model(reactants)))


class ElectronStimulatedDecayModelTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.model = ElectronStimulatedDecayModel()

    def test_p(self):
        results = [list(self.model(r)) for r in parse_spec('p')]
        self.assertEqual([
            [((1, 0), (0, 0), (0, -1)),
             ((-3, -1), (4, 2), (0, -1)),
             ((-7, -3), (4, 2), (4, 2)),
             ((-11, -5), (4, 2), (4, 2), (4, 2))]]
        , results)

    def test_n(self):
        # Will fail at a later stage
        results = [list(self.model(r)) for r in parse_spec('n')]
        self.assertEqual([
            [((1, -1), (0, 0), (0, -1)),
             ((-3, -2), (4, 2), (0, -1)),
             ((-7, -4), (4, 2), (4, 2)),
             ((-11, -6), (4, 2), (4, 2), (4, 2))]]
        , results)

    def test_t(self):
        results = [list(self.model(r)) for r in parse_spec('t')]
        self.assertEqual([
            [((1, 0), (1, 0), (1, 0), (0, 0), (0, -1)),
             ((-1, -1), (4, 2), (0, -1)),
             ((-5, -3), (4, 2), (4, 2)),
             ((-9, -5), (4, 2), (4, 2), (4, 2))]]
        , results)

    def test_H(self):
        results = [list(self.model(r)) for r in parse_spec('H')]
        self.assertEqual([
            [((1, 0), (0, 0), (0, -1)),
             ((-3, -1), (4, 2), (0, -1)),
             ((-7, -3), (4, 2), (4, 2)),
             ((-11, -5), (4, 2), (4, 2), (4, 2))],
            [((1, 0), (1, 0), (0, 0), (0, -1)),
             ((-2, -1), (4, 2), (0, -1)),
             ((-6, -3), (4, 2), (4, 2)),
             ((-10, -5), (4, 2), (4, 2), (4, 2))]]
        , results)

    def test_Li(self):
        results = [list(self.model(r)) for r in parse_spec('Li')]
        self.assertEqual([
            [((6, 2), (0, 0), (0, -1)),
             ((2, 1), (4, 2), (0, -1)),
             ((-2, -1), (4, 2), (4, 2)),
             ((-6, -3), (4, 2), (4, 2), (4, 2))],
            [((7, 2), (0, 0), (0, -1)),
             ((3, 1), (4, 2), (0, -1)),
             ((-1, -1), (4, 2), (4, 2)),
             ((-5, -3), (4, 2), (4, 2), (4, 2))]]
        , results)
