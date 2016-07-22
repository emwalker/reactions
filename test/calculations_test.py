import unittest
import math

import numpy as np

from lenrmc.units import Energy, Distance
from lenrmc.nubase import Nuclides
from lenrmc.combinations import Reaction
from lenrmc.calculations import DecayPower, CoulombBarrier


class CoulombBarrierTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        nuclides = Nuclides.db()
        cls.c = CoulombBarrier(
            nuclides.get(('4He', '0')),
            nuclides.get(('208Pb', '0'))
        )

    def test_coulomb_barrier(self):
        radius = Distance(fermis=1)
        np.testing.assert_approx_equal(236156, self.c.height(radius).kev)

    def test_coulomb_barrier_width(self):
        q_value = Energy.load(mev=6)
        np.testing.assert_approx_equal(39.359344, self.c.width(q_value).fermis)


class GamowSuppressionFactorTest(unittest.TestCase):

    def test_gamow_1(self):
        c = Reaction.load(
            reactants=[(1, ('185Re', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('181Ta', '0'))],
        ).gamow()
        self.assertEqual(52, int(c.value()))

    def test_gamow_2(self):
        c = Reaction.load(
            reactants=[(1, ('58Fe', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('54Cr', '0'))],
        ).gamow()
        self.assertTrue(math.isnan(c.value()))

    def test_gamow_3(self):
        c = Reaction.load(
            reactants=[(1, ('190Pt', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('186Os', '0'))],
        ).gamow()
        self.assertEqual(40, int(c.value()))

    def test_gamow_4(self):
        c = Reaction.load(
            reactants=[(1, ('241Am', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('237Np', '0'))],
        ).gamow()
        self.assertEqual(31, int(c.value()))


class Gamow2Test(unittest.TestCase):

    def test_gamow_factor_1(self):
        c = Reaction.load(
            reactants=[(1, ('212Po', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('208Pb', '0'))],
        ).gamow2()
        self.assertEqual('4.92e+12', '{:.2e}'.format(c.value()))

    def test_gamow_factor_2(self):
        c = Reaction.load(
            reactants=[(1, ('185Re', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('181Ta', '0'))],
        ).gamow2()
        self.assertEqual('2.36e+13', '{:.2e}'.format(c.value()))


class GeigerNuttalLawTest(unittest.TestCase):

    def test_geiger_nuttal_law_1(self):
        c = Reaction.load(
            reactants=[(1, ('185Re', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('181Ta', '0'))],
        ).geiger_nuttal()
        self.assertEqual(24, int(c.value()))

    def test_geiger_nuttal_law_2(self):
        c = Reaction.load(
            reactants=[(1, ('144Nd', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('140Ce', '0'))],
        ).geiger_nuttal()
        self.assertEqual(14, int(c.value()))

    def test_geiger_nuttal_law_3(self):
        c = Reaction.load(
            reactants=[(1, ('212Po', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('208Pb', '0'))],
        ).geiger_nuttal()
        self.assertEqual(-6, int(c.value()))


class AlphaDecayTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.c = Reaction.load(
            reactants=[(1, ('212Po', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('208Pb', '0'))],
        ).alpha_decay()

    def test_nuclear_separation(self):
        np.testing.assert_almost_equal(9.01, self.c.nuclear_separation, 2)

    def test_barrier_height(self):
        np.testing.assert_almost_equal(26.196, self.c.barrier_height, 2)

    def test_alpha_energy(self):
        np.testing.assert_almost_equal(8.785, self.c.alpha_energy, 2)

    def test_alpha_velocity(self):
        np.testing.assert_approx_equal(2.06e7, self.c.alpha_velocity, significant=3)

    def test_barrier_assault_frequency(self):
        np.testing.assert_approx_equal(1.14e21, self.c.barrier_assault_frequency, significant=3)

    def test_gamow_factor(self):
        np.testing.assert_approx_equal(16.77, self.c.gamow_factor, significant=3)

    def test_tunneling_probability(self):
        np.testing.assert_approx_equal(2.697e-15, self.c.tunneling_probability, significant=4)

    def test_decay_constant(self):
        np.testing.assert_approx_equal(3080639, self.c.decay_constant, significant=4)

    def test_half_life(self):
        np.testing.assert_approx_equal(2.2495e-7, self.c.half_life.seconds, significant=4)


class DecayPowerTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.c = Reaction.load(
            reactants=[(1, ('190Pt', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('186Os', '0'))]
        ).alpha_decay().power(moles=1)

    def test_remaining(self):
        np.testing.assert_approx_equal(6.022e23, self.c.remaining(seconds=1), significant=4)
        np.testing.assert_approx_equal(6.022e23, self.c.remaining(seconds=100), significant=4)
        np.testing.assert_approx_equal(6.022e23, self.c.remaining(seconds=3.154e7), significant=4)
        np.testing.assert_approx_equal(1.378e23, self.c.remaining(seconds=1e20), significant=4)

    def test_activity(self):
        np.testing.assert_approx_equal(8880, self.c.activity(seconds=1), significant=4)
        np.testing.assert_approx_equal(2032, self.c.activity(seconds=1e20), significant=4)

    def test_power(self):
        np.testing.assert_approx_equal(4.627712259789981e-09, self.c.power(seconds=1).watts, significant=4)
        np.testing.assert_approx_equal(1.05908972475364e-09, self.c.power(seconds=1e20).watts, significant=4)
