import unittest
import math

import numpy as np

from lenrmc.units import Energy, Distance
from lenrmc.nubase import Nuclides
from lenrmc.combinations import Reaction
from lenrmc.calculations import DecayPower, CoulombBarrier


class CoulombBarrierTest(unittest.TestCase):

    nuclides = Nuclides.db()
    c = CoulombBarrier(
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

    c = Reaction.load(
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
        np.testing.assert_approx_equal(20594372.446746927, self.c.alpha_velocity)

    def test_barrier_assault_frequency(self):
        np.testing.assert_approx_equal(1.142244329315791e+21, self.c.barrier_assault_frequency)

    def test_gamow_factor(self):
        np.testing.assert_approx_equal(16.77331705003231, self.c.gamow_factor)

    def test_tunneling_probability(self):
        np.testing.assert_approx_equal(2.697006070113634e-15, self.c.tunneling_probability)

    def test_decay_constant(self):
        np.testing.assert_approx_equal(3080639, self.c.decay_constant)

    def test_half_life(self):
        np.testing.assert_approx_equal(2.2495326451918876e-07, self.c.half_life.seconds)


class DecayPowerTest(unittest.TestCase):

    pt190 = Reaction.load(
        reactants=[(1, ('190Pt', '0'))],
        daughters=[(1, ('4He', '0')), (1, ('186Os', '0'))]
    ).alpha_decay().power(moles=1)

    am241 = Reaction.load(
        reactants=[(1, ('241Am', '0'))],
        daughters=[(1, ('4He', '0')), (1, ('237Np', '0'))]
    ).alpha_decay().power(moles=1)

    def test_remaining_190Pt(self):
        np.testing.assert_approx_equal(6.02214129e+23, self.pt190.remaining(seconds=1))
        np.testing.assert_approx_equal(6.02214129e+23, self.pt190.remaining(seconds=100))
        np.testing.assert_approx_equal(6.02214129e+23, self.pt190.remaining(seconds=3.154e7))
        np.testing.assert_approx_equal(1.378216190464504e+23, self.pt190.remaining(seconds=1e20))

    def test_activity_190Pt(self):
        np.testing.assert_approx_equal(8880.567784325336, self.pt190.activity(seconds=1))
        np.testing.assert_approx_equal(2032.390425843806, self.pt190.activity(seconds=1e20))

    def test_power_190Pt(self):
        np.testing.assert_approx_equal(4.627712259789981e-09, self.pt190.power(seconds=1).watts)
        np.testing.assert_approx_equal(1.05908972475364e-09, self.pt190.power(seconds=1e20).watts)

    def test_remaining_241Am(self):
        np.testing.assert_approx_equal(6.022141289646228e+23, self.am241.remaining(seconds=1))
        np.testing.assert_approx_equal(6.02214125462277e+23, self.am241.remaining(seconds=100))
        np.testing.assert_approx_equal(6.01099364222362e+23, self.am241.remaining(seconds=3.154e7))
        np.testing.assert_approx_equal(0.0, self.am241.remaining(seconds=1e20))

    def test_activity_241Am(self):
        np.testing.assert_approx_equal(35377229832344.93, self.am241.activity(seconds=1))
        np.testing.assert_approx_equal(0.0, self.am241.activity(seconds=1e20))

    def test_power_241Am(self):
        np.testing.assert_approx_equal(31.955283773279884, self.am241.power(seconds=1).watts)
        np.testing.assert_approx_equal(0.0, self.am241.power(seconds=1e20).watts)
