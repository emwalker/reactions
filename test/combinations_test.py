import unittest
import math

from lenrmc.nubase import parse_spec
from lenrmc.system import System
from lenrmc.combinations import (
    ElectronMediatedDecayModel,
    PionExchangeAndDecayModel,
    Reaction,
    regular_combinations,
    StandardModel,
    StrictPionExchangeModel,
    vectors3,
)


class ReactionsTest(unittest.TestCase):

    def test_daughter_count(self):
        r = Reaction.load(
            reactants=[(1, ('p', '0')), (1, ('7Li', '0'))],
            daughters=[(2, ('4He', '0'))],
        )
        self.assertEqual(2, r.daughter_count)

    def test_q_value(self):
        r = Reaction.load(
            reactants=[(1, ('p', '0')), (1, ('7Li', '0'))],
            daughters=[(2, ('4He', '0'))],
        )
        self.assertEqual(17346.2443, r.q_value.kev)

    def test_electron_capture_q_value(self):
        r = Reaction.load(
            reactants=[(1, ('e-', '0')), (1, ('63Cu', '0'))],
            daughters=[(1, ('ν', '0')), (1, ('63Ni', '0'))],
        )
        # TODO - fix until it goes to zero?
        self.assertEqual(-67, int(r.q_value.kev))

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
        self.assertEqual({'→β-', 't', 'trace'}, r.notes)

    def test_geiger_nuttal_law_1(self):
        r = Reaction.load(
            reactants=[(1, ('185Re', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('181Ta', '0'))],
        )
        self.assertEqual(24, int(r.geiger_nuttal_law()))

    def test_geiger_nuttal_law_2(self):
        r = Reaction.load(
            reactants=[(1, ('144Nd', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('140Ce', '0'))],
        )
        self.assertEqual(14, int(r.geiger_nuttal_law()))

    def test_geiger_nuttal_law_3(self):
        r = Reaction.load(
            reactants=[(1, ('212Po', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('208Pb', '0'))],
        )
        self.assertEqual(-6, int(r.geiger_nuttal_law()))

    def test_gamow_factor_1(self):
        r = Reaction.load(
            reactants=[(1, ('212Po', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('208Pb', '0'))],
        )
        self.assertEqual('4.92e+12', '{:.2e}'.format(r.gamow_factor()))

    def test_gamow_factor_2(self):
        r = Reaction.load(
            reactants=[(1, ('185Re', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('181Ta', '0'))],
        )
        self.assertEqual('2.36e+13', '{:.2e}'.format(r.gamow_factor()))

    def test_gamow_supression_factor(self):
        r = Reaction.load(
            reactants=[(1, ('185Re', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('181Ta', '0'))],
        )
        self.assertEqual(52, int(r.gamow_supression_factor()))

    def test_gamow_supression_factor_2(self):
        r = Reaction.load(
            reactants=[(1, ('58Fe', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('54Cr', '0'))],
        )
        self.assertTrue(math.isnan(r.gamow_supression_factor()))

    def test_gamow_supression_factor_3(self):
        r = Reaction.load(
            reactants=[(1, ('190Pt', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('186Os', '0'))],
        )
        self.assertEqual(40, int(r.gamow_supression_factor()))

    def test_gamow_supression_factor_4(self):
        r = Reaction.load(
            reactants=[(1, ('241Am', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('237Np', '0'))],
        )
        self.assertEqual(31, int(r.gamow_supression_factor()))


class PossibleDaughtersTest(unittest.TestCase):

    def test_outcomes(self):
        reactants = list(parse_spec('p+7Li'))[0]
        outcomes = list(StandardModel()(reactants))
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


class ElectronMediatedDecayModelTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.model = ElectronMediatedDecayModel()

    def test_p(self):
        results = [list(self.model(r)) for r in parse_spec('p')]
        self.assertEqual([
            [((1, 0), (0, 0), (0, -1)),
             ((1, 2), (0, 0), (0, -1)),
             ((1, 3), (0, 0), (0, 0), (0, -1), (0, -1)),
             ((1, 1),),
             ((-3, -1), (4, 2), (0, -1))]]
        , results)

    def test_n(self):
        # Will fail at a later stage
        results = [list(self.model(r)) for r in parse_spec('n')]
        self.assertEqual([
            [((1, -1), (0, 0), (0, -1)),
             ((1, 1), (0, 0), (0, -1)),
             ((1, 2), (0, 0), (0, 0), (0, -1), (0, -1)),
             ((0, -1), (1, 1)),
             ((-3, -2), (4, 2), (0, -1))]]
        , results)

    def test_t(self):
        results = [list(self.model(r)) for r in parse_spec('t')]
        self.assertEqual([
            [((1, 0), (1, 0), (1, 0), (0, 0), (0, -1)),
             ((3, 2), (0, 0), (0, -1)),
             ((1, 1), (1, 1), (1, 1), (0, 0), (0, 0), (0, -1), (0, -1)),
             ((1, 0), (1, 0), (1, 1)),
             ((-1, -1), (4, 2), (0, -1))]]
        , results)

    def test_H(self):
        results = [list(self.model(r)) for r in parse_spec('H')]
        self.assertEqual([
            [((1, 0), (0, 0), (0, -1)),
             ((1, 2), (0, 0), (0, -1)),
             ((1, 3), (0, 0), (0, 0), (0, -1), (0, -1)),
             ((1, 1),),
             ((-3, -1), (4, 2), (0, -1))],
            [((1, 0), (1, 0), (0, 0), (0, -1)),
             ((1, 1), (1, 1), (0, 0), (0, -1)),
             ((2, 3), (0, 0), (0, 0), (0, -1), (0, -1)),
             ((1, 0), (1, 1)),
             ((-2, -1), (4, 2), (0, -1))],
            [((1, 0), (1, 0), (1, 0), (0, 0), (0, -1)),
             ((3, 2), (0, 0), (0, -1)),
             ((1, 1), (1, 1), (1, 1), (0, 0), (0, 0), (0, -1), (0, -1)),
             ((1, 0), (1, 0), (1, 1)),
             ((-1, -1), (4, 2), (0, -1))]]
        , results)

    def test_Li(self):
        results = [list(self.model(r)) for r in parse_spec('Li')]
        self.assertEqual([
            [((6, 2), (0, 0), (0, -1)),
             ((6, 4), (0, 0), (0, -1)),
             ((6, 5), (0, 0), (0, 0), (0, -1), (0, -1)),
             ((5, 2), (1, 1)),
             ((2, 1), (4, 2), (0, -1))],
            [((7, 2), (0, 0), (0, -1)),
             ((7, 4), (0, 0), (0, -1)),
             ((7, 5), (0, 0), (0, 0), (0, -1), (0, -1)),
             ((6, 2), (1, 1)),
             ((3, 1), (4, 2), (0, -1))]]
        , results)

    def test_90Sr(self):
        results = [list(self.model(r)) for r in parse_spec('90Sr')]
        self.assertEqual([
            [((90, 37), (0, 0), (0, -1)),
             ((90, 39), (0, 0), (0, -1)),
             ((90, 40), (0, 0), (0, 0), (0, -1), (0, -1)),
             ((89, 37), (1, 1)),
             ((86, 36), (4, 2), (0, -1))]]
        , results)

    def test_90Sr_2(self):
        s = System.parse('90Sr', model='induced-decay', lb=-1000)
        reactions = list(r for c in s._combinations for r in c.reactions())
        self.assertEqual(2, len(reactions))
        reaction = reactions[0]
        self.assertEqual(('90Y', '0'), reaction.rvalues[0][1].signature)
        self.assertEqual(545.9997699999949, reaction.q_value.kev)
