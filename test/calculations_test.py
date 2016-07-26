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

    pt190 = System.load('190Pt', model='induced-decay') \
        .alpha_decay(moles=1, seconds=1, isotopic_fraction=1)

    def test_remaining_190Pt(self):
        np.testing.assert_approx_equal(6.02214129e+23, self.pt190.remaining_active_atoms())
        np.testing.assert_approx_equal(6.02214129e+23, self.pt190.remaining_active_atoms(seconds=100))
        np.testing.assert_approx_equal(6.02214129e+23, self.pt190.remaining_active_atoms(seconds=3.154e7))
        np.testing.assert_approx_equal(1.3760309564333265e+23, self.pt190.remaining_active_atoms(seconds=1e20))

    def test_activity_190Pt(self):
        np.testing.assert_approx_equal(8890.123782859177, self.pt190.activity(seconds=1))
        np.testing.assert_approx_equal(2031.3514649767333, self.pt190.activity(seconds=1e20))

    def test_power_190Pt(self):
        np.testing.assert_approx_equal(4.63269194269355e-09, self.pt190.power(seconds=1).watts)
        np.testing.assert_approx_equal(1.0585483165848424e-09, self.pt190.power(seconds=1e20).watts)

    def test_241Am(self):
        scenario = System.load('241Am', model='induced-decay') \
            .alpha_decay(isotopic_fraction=1, moles=1, seconds=1)
        # Remaining
        np.testing.assert_approx_equal(6.022141289646228e+23, scenario.remaining_active_atoms())
        np.testing.assert_approx_equal(6.02214125462277e+23, scenario.remaining_active_atoms(seconds=100))
        np.testing.assert_approx_equal(6.0109829989678545e+23, scenario.remaining_active_atoms(seconds=3.154e7))
        np.testing.assert_approx_equal(0.0, scenario.remaining_active_atoms(seconds=1e20))
        # Activity
        np.testing.assert_approx_equal(35411037706734.586, scenario.activity(seconds=1))
        np.testing.assert_approx_equal(0.0, scenario.activity(seconds=1e20))
        # Power
        np.testing.assert_approx_equal(31.985821501219938, scenario.power(seconds=1).watts)
        np.testing.assert_approx_equal(0.0, scenario.power(seconds=1e20).watts)

    def test_241Am_power(self):
        # Sanity check against a number found in Wikipedia
        # https://en.wikipedia.org/wiki/Isotopes_of_americium
        # (1 kg * 1000 g/kg) / (243 g/mole)
        moles = 1e3 / 243
        scenario = System.load('241Am', model='induced-decay') \
            .alpha_decay(seconds=1, isotopic_fraction=1, moles=moles)
        # Should be 114 watts/kg
        np.testing.assert_approx_equal(131.62889506674873, scenario.power().watts)
        np.testing.assert_approx_equal(0.0, scenario.power(seconds=1e20).watts)

    def test_screened_190Pt(self):
        pt190 = System.load('190Pt', model='induced-decay') \
            .alpha_decay(seconds=1, screening=11, moles=1, isotopic_fraction=1)
        # Remaining
        np.testing.assert_approx_equal(6.02214129e+23, pt190.remaining_active_atoms())
        np.testing.assert_approx_equal(6.02214129e+23, pt190.remaining_active_atoms(seconds=100))
        np.testing.assert_approx_equal(6.021836528709286e+23, pt190.remaining_active_atoms(seconds=3.154e7))
        np.testing.assert_approx_equal(0.0, pt190.remaining_active_atoms(seconds=1e20))
        # Activity
        np.testing.assert_approx_equal(967090273179.6211, pt190.activity(seconds=1))
        np.testing.assert_approx_equal(0.0, pt190.activity(seconds=1e20))
        # Power
        np.testing.assert_approx_equal(0.5039560107199808, pt190.power().watts)
        np.testing.assert_approx_equal(0.0, pt190.power(seconds=1e20).watts)


class PlatinumAlphaDecayTest(unittest.TestCase):

    ##
    # Model production of of 22,522,522,523 4He/s from 0.02193926719 mol pt
    # over a period of 4440 seconds
    # http://lenr-canr.org/acrobat/MilesMcorrelatio.pdf
    #

    screening = 32.046
    moles = 0.02193926719
    active_fraction = 1e-6

    scenario = System.load('Pt', model='induced-decay') \
        .alpha_decay(
            seconds=1,
            moles=moles,
            active_fraction=active_fraction,
            screening=screening,
        )

    def test_miles_4He_study(self):
        np.testing.assert_approx_equal(22522522523, self.scenario.activity(), significant=3)
        np.testing.assert_approx_equal(0.00876108327676111, self.scenario.power().watts)

    def test_elemental_Pt(self):
        scenario = System.load('Pt', model='induced-decay').alpha_decay(seconds=1, moles=1, active_fraction=1)
        np.testing.assert_approx_equal(1.0668148541893832, scenario.activity())

    def test_screened_Pt(self):
        scenario = System.load('Pt', model='induced-decay').alpha_decay(screening=11, seconds=1, moles=1, active_fraction=1)
        np.testing.assert_approx_equal(116050834.10601069, scenario.activity())
        np.testing.assert_approx_equal(6.0474721800430736e-05, scenario.power().watts)

    def test_parent_z(self):
        np.testing.assert_allclose([
            78, 78, 78, 78, 78, 78,
        ], self.scenario.df.parent_z)

    def test_heavier_daughter_z(self):
        np.testing.assert_allclose([
            76, 76, 76, 76, 76, 76
        ], self.scenario.df.heavier_daughter_z)

    def test_lighter_daughter_a(self):
        np.testing.assert_allclose([
            4, 4, 4, 4, 4, 4
        ], self.scenario.df.lighter_daughter_a)

    def test_heavier_daughter_a(self):
        np.testing.assert_allclose([
            186, 188, 190, 191, 192, 194
        ], self.scenario.df.heavier_daughter_a)

    def test_screened_heavier_z(self):
        np.testing.assert_allclose(43.954, self.scenario.df.screened_heavier_daughter_z[0])

    def test_isotope(self):
        np.testing.assert_equal([
            '190Pt',
            '192Pt',
            '194Pt',
            '195Pt',
            '196Pt',
            '198Pt',
        ], self.scenario.df.isotope.values)

    def test_barrier_height_mev(self):
        np.testing.assert_allclose([
            14.458853,
            14.418562,
            14.378778,
            14.359072,
            14.339488,
            14.300682
        ], self.scenario.df.barrier_height_mev)

    def test_alpha_mass_mev(self):
        np.testing.assert_allclose([
            3728.40116,
            3728.40116,
            3728.40116,
            3728.40116,
            3728.40116,
            3728.40116,
        ], self.scenario.df.alpha_mass_mev)

    def test_heavier_daughter_mass_mev(self):
        np.testing.assert_allclose([
            173214.892946,
            175079.744168,
            176945.16219,
            177878.968751,
            178810.975812,
            180677.410634,
        ], self.scenario.df.heavier_daughter_mass_mev)

    def test_alpha_ke_mev(self):
        np.testing.assert_allclose([
            3.183951,
            2.371874,
            1.490479,
            1.151548,
            0.795497,
            0.104429,
        ], self.scenario.df.alpha_ke_mev, rtol=1e-5)

    def test_alpha_velocity_m_per_s(self):
        np.testing.assert_allclose([
            12389607.706447,
            10693508.438471,
            8476903.862868,
            7451013.125387,
            6192895.791641,
            2243808.726547
        ], self.scenario.df.alpha_velocity_m_per_s)

    def test_nuclear_separation_fm(self):
        np.testing.assert_allclose([
            8.754802,
            8.779266,
            8.803558,
            8.81564,
            8.827679,
            8.851634,
        ], self.scenario.df.nuclear_separation_fm)

    def test_barrier_assault_frequency(self):
        np.testing.assert_allclose([
            7.075892e+20,
            6.090206e+20,
            4.814476e+20,
            4.226020e+20,
            3.507658e+20,
            1.267455e+20
        ], self.scenario.df.barrier_assault_frequency, rtol=1e-6)

    def test_gamow_factor(self):
        np.testing.assert_allclose([
            20.736394,
            28.143476,
            42.569015,
            52.244928,
            68.575286,
            240.002875
        ], self.scenario.df.gamow_factor)

    def test_tunneling_probability(self):
        np.testing.assert_allclose([
            9.740854e-019,
            3.588289e-025,
            1.059311e-037,
            4.174765e-046,
            2.730631e-060,
            3.436779e-209
        ], self.scenario.df.tunneling_probability, rtol=1e-6)

    def test_decay_constant(self):
        np.testing.assert_allclose([
            6.892524e+002,
            2.185342e-004,
            5.100026e-017,
            1.764264e-025,
            9.578119e-040,
            4.355961e-189
        ], self.scenario.df.decay_constant, rtol=1e-6)

    def test_half_life(self):
        np.testing.assert_allclose([
            1.005651e-003,
            3.171802e+003,
            1.359105e+016,
            3.928818e+024,
            7.236778e+038,
            1.591261e+188
        ], self.scenario.df.half_life, rtol=1e-6)

    def test_isotopic_abundance(self):
        np.testing.assert_allclose([
            1.200000e-02,
            7.820000e-01,
            3.286000e+01,
            3.378000e+01,
            2.521000e+01,
            7.360000e+00,
        ], self.scenario.df.isotopic_abundance)

    def test_isotopic_fraction(self):
        np.testing.assert_allclose([
            1.200000e-04,
            7.820000e-03,
            3.286000e-01,
            3.378000e-01,
            2.521000e-01,
            7.360000e-02,
        ], self.scenario.df.isotopic_fraction)

    def test_starting_moles(self):
        np.testing.assert_allclose([
            2.632712e-06,
            1.715651e-04,
            7.209243e-03,
            7.411084e-03,
            5.530889e-03,
            1.614730e-03,
        ], self.scenario.df.starting_moles, rtol=1e-6)

    def test_active_fraction(self):
        np.testing.assert_allclose([
            1.000000e-06,
            1.000000e-06,
            1.000000e-06,
            1.000000e-06,
            1.000000e-06,
            1.000000e-06,
        ], self.scenario.df.active_fraction, rtol=1e-6)

    def test_starting_active_moles(self):
        np.testing.assert_allclose([
            2.632712e-12,
            1.715651e-10,
            7.209243e-09,
            7.411084e-09,
            5.530889e-09,
            1.614730e-09,
        ], self.scenario.df.starting_active_moles, rtol=1e-6)

    def test_starting_active_atoms(self):
        np.testing.assert_allclose([
            1.585456e+12,
            1.033189e+14,
            4.341508e+15,
            4.463060e+15,
            3.330780e+15,
            9.724133e+14,
        ], self.scenario.df.starting_active_atoms, rtol=1e-6)

    def test_remaining_active_atoms(self):
        np.testing.assert_allclose([
            7.272133e-288,
            1.032963e+014,
            4.341508e+015,
            4.463060e+015,
            3.330780e+015,
            9.724133e+014
        ], self.scenario.df.remaining_active_atoms, rtol=1e-6)

    def test_activity(self):
        np.testing.assert_allclose([
            5.012335e-285,
            2.257378e+010,
            2.214180e-001,
            7.874016e-010,
            3.190260e-024,
            4.235794e-174
        ], self.scenario.df.activity, rtol=1e-6)

    def test_watts(self):
        np.testing.assert_allclose([
            2.611955e-297,
            8.761083e-003,
            5.398897e-014,
            1.483192e-022,
            4.150856e-037,
            7.233340e-188
        ], self.scenario.df.watts, rtol=1e-6)

    def test_total_activity(self):
        df = self.scenario.df
        np.testing.assert_approx_equal(22522522523, df.activity.sum(), significant=3)
        np.testing.assert_approx_equal(5.01233468945662e-285, df.activity[df.isotope == '190Pt'][0])

    def test_total_power(self):
        np.testing.assert_approx_equal(0.00876108327676111, self.scenario.df.watts.sum())


class PoloniumAlphaDecayTest(unittest.TestCase):

    scenario = System.load('212Po', model='induced-decay') \
        .alpha_decay(seconds=1, moles=1, active_fraction=1, isotopic_fraction=1)

    def test_nuclear_separation(self):
        np.testing.assert_allclose([9.014871826539528], self.scenario.df.nuclear_separation_fm)

    def test_barrier_height(self):
        np.testing.assert_allclose([26.1967118938676], self.scenario.df.barrier_height_mev, rtol=1e-4)

    def test_alpha_ke(self):
        np.testing.assert_allclose([8.78], self.scenario.df.alpha_ke_mev, rtol=1e-3)

    def test_radius_for_alpha_ke_fm(self):
        np.testing.assert_allclose([26.89749430523918], self.scenario.df.radius_for_alpha_ke_fm, rtol=1e-3)

    def test_barrier_width_fm(self):
        np.testing.assert_allclose([17.882622478699652], self.scenario.df.barrier_width_fm, rtol=1e-3)

    def test_barrier_assault_frequency(self):
        np.testing.assert_allclose([1.142126641655716e21], self.scenario.df.barrier_assault_frequency, rtol=1e-3)

    def test_alpha_v_over_c_m_per_s(self):
        np.testing.assert_allclose([0.068648], self.scenario.df.alpha_v_over_c_m_per_s, rtol=1e-5)

    def test_alpha_velocity_m_per_s(self):
        np.testing.assert_allclose([2.06e7], self.scenario.df.alpha_velocity_m_per_s, rtol=1e-3)

    def test_tunneling_probability(self):
        np.testing.assert_allclose([2.636693524272448e-15], self.scenario.df.tunneling_probability, rtol=1e-1)
