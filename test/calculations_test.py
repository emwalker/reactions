import unittest
import math

import numpy as np
import pandas as pd
import pandas.util.testing as pdt

from lenrmc.units import Energy, Distance
from lenrmc.nubase import Nuclides
from lenrmc.system import System
from lenrmc.combinations import Reaction
from lenrmc.calculations import CoulombBarrier


nuclides = Nuclides.db()


class CoulombBarrierTest(unittest.TestCase):

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


class IsotopicAlphaDecayTest(unittest.TestCase):

    c = Reaction.load(
        reactants=[(1, ('212Po', '0'))],
        daughters=[(1, ('4He', '0')), (1, ('208Pb', '0'))],
    )._alpha_decay_calculation(moles=1)

    def test_parent(self):
        self.assertEqual(nuclides.get(('212Po', '0')), self.c.parent)

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


class DecayTest(unittest.TestCase):

    pt190 = Reaction.load(
        reactants=[(1, ('190Pt', '0'))],
        daughters=[(1, ('4He', '0')), (1, ('186Os', '0'))]
    ).alpha_decay(moles=1)

    def test_atomic_number(self):
        self.assertEqual(78, self.pt190.parent.atomic_number)

    def test_abundance(self):
        self.assertEqual(0.012, self.pt190.parent.isotopic_abundance)


class AlphaDecayTest(unittest.TestCase):

    pt190 = System.load('190Pt', model='induced-decay').alpha_decay(moles=1)

    def test_remaining_190Pt(self):
        kwargs = dict(moles=1, isotopic_fraction=1)
        np.testing.assert_approx_equal(6.02214129e+23, self.pt190.remaining_active_atoms(seconds=1, **kwargs))
        np.testing.assert_approx_equal(6.02214129e+23, self.pt190.remaining_active_atoms(seconds=100, **kwargs))
        np.testing.assert_approx_equal(6.02214129e+23, self.pt190.remaining_active_atoms(seconds=3.154e7, **kwargs))
        np.testing.assert_approx_equal(1.378216190464504e+23, self.pt190.remaining_active_atoms(seconds=1e20, **kwargs))

    def test_activity_190Pt(self):
        kwargs = dict(moles=1, isotopic_fraction=1)
        np.testing.assert_approx_equal(8880.567784325336, self.pt190.activity(seconds=1, **kwargs))
        np.testing.assert_approx_equal(2032.390425843806, self.pt190.activity(seconds=1e20, **kwargs))

    def test_power_190Pt(self):
        kwargs = dict(moles=1, isotopic_fraction=1)
        np.testing.assert_approx_equal(4.627712259789981e-09, self.pt190.power(seconds=1, **kwargs).watts)
        np.testing.assert_approx_equal(1.05908972475364e-09, self.pt190.power(seconds=1e20, **kwargs).watts)

    def test_241Am(self):
        kwargs = dict(isotopic_fraction=1, moles=1)
        am241 = System.load('241Am', model='induced-decay').alpha_decay(moles=1)
        # Remaining
        np.testing.assert_approx_equal(6.022141289646228e+23, am241.remaining_active_atoms(seconds=1, **kwargs))
        np.testing.assert_approx_equal(6.02214125462277e+23, am241.remaining_active_atoms(seconds=100, **kwargs))
        np.testing.assert_approx_equal(6.01099364222362e+23, am241.remaining_active_atoms(seconds=3.154e7, **kwargs))
        np.testing.assert_approx_equal(0.0, am241.remaining_active_atoms(seconds=1e20, **kwargs))
        # Activity
        np.testing.assert_approx_equal(35377229832344.93, am241.activity(seconds=1, **kwargs))
        np.testing.assert_approx_equal(0.0, am241.activity(seconds=1e20, **kwargs))
        # Power
        np.testing.assert_approx_equal(31.955283773279884, am241.power(seconds=1, **kwargs).watts)
        np.testing.assert_approx_equal(0.0, am241.power(seconds=1e20, **kwargs).watts)

    def test_241Am_power(self):
        # Sanity check against a number found in Wikipedia
        # https://en.wikipedia.org/wiki/Isotopes_of_americium
        # (1 kg * 1000 g/kg) / (243 g/mole)
        moles = 1e3 / 243
        kwargs = dict(moles=moles, isotopic_fraction=1)
        am241 = System.load('241Am', model='induced-decay').alpha_decay(**kwargs)
        # Should be 114 watts/kg
        np.testing.assert_approx_equal(131.50322540444398, am241.power(seconds=1, **kwargs).watts)
        np.testing.assert_approx_equal(0.0, am241.power(seconds=1e20, **kwargs).watts)

    def test_screened_190Pt(self):
        kwargs = dict(screening=11, moles=1, isotopic_fraction=1)
        pt190 = System.load('190Pt', model='induced-decay').alpha_decay(**kwargs)
        # Remaining
        np.testing.assert_approx_equal(6.02214129e+23, pt190.remaining_active_atoms(seconds=1, **kwargs))
        np.testing.assert_approx_equal(6.02214129e+23, pt190.remaining_active_atoms(seconds=100, **kwargs))
        np.testing.assert_approx_equal(6.021836528709286e+23, pt190.remaining_active_atoms(seconds=3.154e7, **kwargs))
        np.testing.assert_approx_equal(0.0, pt190.remaining_active_atoms(seconds=1e20, **kwargs))
        # Activity
        np.testing.assert_approx_equal(966293603266.9395, pt190.activity(seconds=1, **kwargs))
        np.testing.assert_approx_equal(0.0, pt190.activity(seconds=1e20, **kwargs))
        # Power
        np.testing.assert_approx_equal(0.5035408616876824, pt190.power(seconds=1, **kwargs).watts)
        np.testing.assert_approx_equal(0.0, pt190.power(seconds=1e20, **kwargs).watts)


class AlphaDecayColumnTest(unittest.TestCase):

    ##
    # Model production of of 22,522,522,523 4He/s from 0.02193926719 mol pt
    # over a period of 4440 seconds
    # http://lenr-canr.org/acrobat/MilesMcorrelatio.pdf
    #

    screening = 32.045
    moles = 0.02193926719
    active_fraction = 1e-6

    pt190 = Reaction.load(
        reactants=[(1, ('190Pt', '0'))],
        daughters=[(1, ('4He', '0')), (1, ('186Os', '0'))]
    ).alpha_decay(screening=screening, moles=moles * active_fraction)

    miles = System.load('Pt', model='induced-decay') \
        .alpha_decay(screening=screening, moles=moles)
    df = miles.dataframe(seconds=1, moles=moles, active_fraction=active_fraction)

    def test_elemental_Pt(self):
        pt = System.load('Pt', model='induced-decay').alpha_decay(moles=1)
        np.testing.assert_approx_equal(1.0656681343649925, pt.activity(seconds=1, moles=1, active_fraction=1))

    def test_screened_Pt(self):
        pt = System.load('Pt', model='induced-decay').alpha_decay(screening=11, moles=1)
        np.testing.assert_approx_equal(115955233.71509394, pt.activity(seconds=1, moles=1, active_fraction=1))
        np.testing.assert_approx_equal(6.042490391601354e-05, pt.power(seconds=1, moles=1, active_fraction=1).watts)

    def test_miles_4He_study(self):
        activity = self.miles.activity(seconds=1, moles=self.moles, active_fraction=self.active_fraction)
        np.testing.assert_approx_equal(22522522523, activity, significant=4)
        power = self.miles.power(seconds=1, moles=self.moles, active_fraction=self.active_fraction)
        np.testing.assert_approx_equal(0.008740219935185282, power.watts)

    def test_atomic_number(self):
        np.testing.assert_allclose([
            78,
            78,
            78,
            78,
            78,
            78,
        ], self.df.atomic_number)

    def test_isotope(self):
        np.testing.assert_equal([
            '190Pt',
            '192Pt',
            '194Pt',
            '195Pt',
            '196Pt',
            '198Pt',
        ], self.df.isotope.values)

    def test_gamow_factor(self):
        np.testing.assert_allclose([
            20.737672,
            28.145014,
            42.571041,
            52.247275,
            68.578169,
            240.011295,
        ], self.df.gamow_factor)

    def test_decay_constant(self):
        np.testing.assert_allclose([
            6.879694e+002,
            2.180137e-004,
            5.082914e-017,
            1.757218e-025,
            9.529659e-040,
            4.286183e-189
        ], self.df.decay_constant, rtol=1e-6)

    def test_isotopic_abundance(self):
        np.testing.assert_allclose([
            1.200000e-02,
            7.820000e-01,
            3.286000e+01,
            3.378000e+01,
            2.521000e+01,
            7.360000e+00,
        ], self.df.isotopic_abundance)

    def test_isotopic_fraction(self):
        np.testing.assert_allclose([
            1.200000e-04,
            7.820000e-03,
            3.286000e-01,
            3.378000e-01,
            2.521000e-01,
            7.360000e-02,
        ], self.df.isotopic_fraction)

    def test_starting_moles(self):
        np.testing.assert_allclose([
            2.632712e-06,
            1.715651e-04,
            7.209243e-03,
            7.411084e-03,
            5.530889e-03,
            1.614730e-03,
        ], self.df.starting_moles, rtol=1e-6)

    def test_active_fraction(self):
        np.testing.assert_allclose([
            1.000000e-06,
            1.000000e-06,
            1.000000e-06,
            1.000000e-06,
            1.000000e-06,
            1.000000e-06,
        ], self.df.active_fraction, rtol=1e-6)

    def test_starting_active_moles(self):
        np.testing.assert_allclose([
            2.632712e-12,
            1.715651e-10,
            7.209243e-09,
            7.411084e-09,
            5.530889e-09,
            1.614730e-09,
        ], self.df.starting_active_moles, rtol=1e-6)

    def test_starting_active_atoms(self):
        np.testing.assert_allclose([
            1.585456e+12,
            1.033189e+14,
            4.341508e+15,
            4.463060e+15,
            3.330780e+15,
            9.724133e+14,
        ], self.df.starting_active_atoms, rtol=1e-6)

    def test_remaining_active_atoms(self):
        np.testing.assert_allclose([
            2.623280e-287,
            1.032964e+014,
            4.341508e+015,
            4.463060e+015,
            3.330780e+015,
            9.724133e+014,
        ], self.df.remaining_active_atoms, rtol=1e-6)

    def test_activity(self):
        np.testing.assert_allclose([
            1.804736e-284,
            2.252003e+010,
            2.206751e-001,
            7.842568e-010,
            3.174119e-024,
            4.167941e-174,
        ], self.df.activity, rtol=1e-6)

    def test_watts(self):
        np.testing.assert_allclose([
            9.404578e-297,
            8.740220e-003,
            5.380783e-014,
            1.477268e-022,
            4.129855e-037,
            7.117470e-188,
        ], self.df.watts, rtol=1e-6)

    def test_total_activity(self):
        df = self.df
        np.testing.assert_approx_equal(22522522523, df.activity.sum(), significant=4)
        np.testing.assert_approx_equal(1.8047361656075522e-284, df.activity[df.isotope == '190Pt'][0])

    def test_total_power(self):
        np.testing.assert_approx_equal(0.008740219935185282, self.df.watts.sum())
