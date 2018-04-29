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

    @classmethod
    def setUpClass(cls):
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
        np.testing.assert_approx_equal(51.6243422304808, c.value())

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
        np.testing.assert_approx_equal(39.086112497598315, c.value())

    def test_gamow_4(self):
        c = Reaction.load(
            reactants=[(1, ('241Am', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('237Np', '0'))],
        ).gamow()
        np.testing.assert_approx_equal(30.05677680959287, c.value())

    def test_gamow_5(self):
        c = Reaction.load(
            reactants=[(1, ('8Be', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('4He', '0'))],
        ).gamow()
        np.testing.assert_approx_equal(5.474197090566792, c.value())

    def test_gamow_6(self):
        c = Reaction.load(
            reactants=[(1, ('106Pd', '0'))],
            daughters=[(1, ('46Ca', '0')), (1, ('60Fe', '0'))],
        ).gamow()
        np.testing.assert_approx_equal(142.1826806102932, c.value())

    def test_gamow_7(self):
        c = Reaction.load(
            reactants=[(1, ('106Pd', '0'))],
            daughters=[(1, ('90Sr', '0')), (1, ('16O', '0'))],
        ).gamow()
        np.testing.assert_approx_equal(456.3426542907934, c.value())


class Gamow2Test(unittest.TestCase):

    def test_gamow_factor_1(self):
        c = Reaction.load(
            reactants=[(1, ('212Po', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('208Pb', '0'))],
        ).gamow2()
        np.testing.assert_approx_equal(4920433262606.988, c.value())

    def test_gamow_factor_2(self):
        c = Reaction.load(
            reactants=[(1, ('185Re', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('181Ta', '0'))],
        ).gamow2()
        np.testing.assert_approx_equal(23635102096516.79, c.value())

    def test_gamow_factor_3(self):
        c = Reaction.load(
            reactants=[(1, ('8Be', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('4He', '0'))],
        ).gamow2()
        np.testing.assert_approx_equal(3725992198095.755, c.value())


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

    @classmethod
    def setUpClass(cls):
        cls.pt190 = Reaction.load(
            reactants=[(1, ('190Pt', '0'))],
            daughters=[(1, ('4He', '0')), (1, ('186Os', '0'))]
        ).decay(moles=1)

    def test_atomic_number(self):
        self.assertEqual(78, self.pt190.parent.atomic_number)

    def test_abundance(self):
        self.assertEqual(0.012, self.pt190.parent.isotopic_abundance)


class AlphaDecayTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.pt190 = System.load('190Pt', model='induced-decay') \
            .hyperphysics(moles=1, seconds=1, isotopic_fraction=1)

    def test_remaining_190Pt(self):
        np.testing.assert_approx_equal(6.02214129e+23, self.pt190.remaining_active_atoms())
        np.testing.assert_approx_equal(6.02214129e+23, self.pt190.remaining_active_atoms(seconds=100))
        np.testing.assert_approx_equal(6.02214129e+23, self.pt190.remaining_active_atoms(seconds=3.154e7))
        np.testing.assert_approx_equal(1.3788626498604035e+23, self.pt190.remaining_active_atoms(seconds=1e20))

    def test_activity_190Pt(self):
        np.testing.assert_approx_equal(8877.743730232462, self.pt190.activity(seconds=1))
        np.testing.assert_approx_equal(2032.697118046182, self.pt190.activity(seconds=1e20))

    def test_power_190Pt(self):
        np.testing.assert_approx_equal(4.626240629814815e-09, self.pt190.power(seconds=1).watts)
        np.testing.assert_approx_equal(1.0592495437313659e-09, self.pt190.power(seconds=1e20).watts)

    def test_241Am(self):
        scenario = System.load('241Am', model='induced-decay') \
            .hyperphysics(isotopic_fraction=1, moles=1, seconds=1)
        # Remaining
        np.testing.assert_approx_equal(6.022141289646228e+23, scenario.remaining_active_atoms())
        np.testing.assert_approx_equal(6.02214125462277e+23, scenario.remaining_active_atoms(seconds=100))
        np.testing.assert_approx_equal(6.010996320511956e+23, scenario.remaining_active_atoms(seconds=3.154e7))
        np.testing.assert_approx_equal(0.0, scenario.remaining_active_atoms(seconds=1e20))
        # Activity
        np.testing.assert_approx_equal(35368722366460.305, scenario.activity(seconds=1))
        np.testing.assert_approx_equal(0.0, scenario.activity(seconds=1e20))
        # Power
        np.testing.assert_approx_equal(31.94759921211376, scenario.power(seconds=1).watts)
        np.testing.assert_approx_equal(0.0, scenario.power(seconds=1e20).watts)

    def test_241Am_power(self):
        # Sanity check against a number found in Wikipedia
        # https://en.wikipedia.org/wiki/Isotopes_of_americium
        # (1 kg * 1000 g/kg) / (243 g/mole)
        moles = 1e3 / 243
        scenario = System.load('241Am', model='induced-decay') \
            .hyperphysics(seconds=1, isotopic_fraction=1, moles=moles)
        # Should be 114 watts/kg
        np.testing.assert_approx_equal(131.4716016959414, scenario.power().watts)
        np.testing.assert_approx_equal(0.0, scenario.power(seconds=1e20).watts)

    def test_screened_190Pt(self):
        pt190 = System.load('190Pt', model='induced-decay') \
            .hyperphysics(seconds=1, screening=11, moles=1, isotopic_fraction=1)
        # Remaining
        np.testing.assert_approx_equal(6.02214129e+23, pt190.remaining_active_atoms())
        np.testing.assert_approx_equal(6.02214129e+23, pt190.remaining_active_atoms(seconds=100))
        np.testing.assert_approx_equal(6.021836528709286e+23, pt190.remaining_active_atoms(seconds=3.154e7))
        np.testing.assert_approx_equal(0.0, pt190.remaining_active_atoms(seconds=1e20))
        # Activity
        np.testing.assert_approx_equal(965968183447.2651, pt190.activity(seconds=1))
        np.testing.assert_approx_equal(0.0, pt190.activity(seconds=1e20))
        # Power
        np.testing.assert_approx_equal(0.5033712836465414, pt190.power().watts)
        np.testing.assert_approx_equal(0.0, pt190.power(seconds=1e20).watts)


##
# Model production of of 22,522,522,523 4He/s from 0.02193926719 mol pt
# over a period of 4440 seconds
# http://lenr-canr.org/acrobat/MilesMcorrelatio.pdf
#
screening = 32.046
moles = 0.02193926719
active_fraction = 1e-6


class HyperphysicsPlatinumAlphaDecayTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.scenario = System.load('Pt', model='induced-decay') \
            .hyperphysics(
                seconds=1,
                moles=moles,
                active_fraction=active_fraction,
                screening=screening,
            )

    def test_miles_4He_study(self):
        np.testing.assert_approx_equal(22522522523, self.scenario.activity(), significant=3)
        np.testing.assert_approx_equal(0.008753265223103896, self.scenario.power().watts)

    def test_elemental_Pt(self):
        scenario = System.load('Pt', model='induced-decay').hyperphysics(seconds=1, moles=1, active_fraction=1)
        np.testing.assert_approx_equal(1.0668148541893832, scenario.activity(), significant=3)

    def test_screened_Pt(self):
        scenario = System.load('Pt', model='induced-decay').hyperphysics(screening=11, seconds=1, moles=1, active_fraction=1)
        np.testing.assert_approx_equal(116050834.10601069, scenario.activity(), significant=3)
        np.testing.assert_approx_equal(6.0474721800430736e-05, scenario.power().watts, significant=3)

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

    def test_parent(self):
        np.testing.assert_equal([
            '190Pt',
            '192Pt',
            '194Pt',
            '195Pt',
            '196Pt',
            '198Pt',
        ], self.scenario.df.parent.values)

    def test_gamow_factor(self):
        np.testing.assert_allclose([
            20.736756,
            28.143922,
            42.569622,
            52.245641,
            68.576178,
            240.005615
        ], self.scenario.df.gamow_factor)

    def test_tunneling_probability(self):
        np.testing.assert_allclose([
            9.733810e-019,
            3.585087e-025,
            1.058025e-037,
            4.168814e-046,
            2.725768e-060,
            3.417993e-209
        ], self.scenario.df.tunneling_probability, rtol=1e-6)

    def test_partial_decay_constant(self):
        np.testing.assert_allclose([
            6.887539e+002,
            2.183392e-004,
            5.093837e-017,
            1.761749e-025,
            9.561062e-040,
            4.332151e-189
        ], self.scenario.df.partial_decay_constant, rtol=1e-6)

    def test_isotope_decay_constant(self):
        np.testing.assert_allclose([
            6.887539e+002,
            2.183392e-004,
            5.093837e-017,
            1.761749e-025,
            9.561062e-040,
            4.332151e-189
        ], self.scenario.df.isotope_decay_constant, rtol=1e-6)

    def test_partial_half_life(self):
        np.testing.assert_allclose([
            1.006379e-003,
            3.174635e+003,
            1.360756e+016,
            3.934426e+024,
            7.249688e+038,
            1.600007e+188
        ], self.scenario.df.partial_half_life, rtol=1e-6)

    def test_isotopic_abundance(self):
        np.testing.assert_allclose([
            1.200000e-02,
            7.820000e-01,
            3.286000e+01,
            3.378000e+01,
            2.521000e+01,
            7.360000e+00,
        ], self.scenario.df.isotopic_abundance)

    def test_parent_fraction(self):
        np.testing.assert_allclose([
            1.200000e-04,
            7.820000e-03,
            3.286000e-01,
            3.378000e-01,
            2.521000e-01,
            7.360000e-02,
        ], self.scenario.df.parent_fraction)

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
            1.197061e-287,
            1.032964e+014,
            4.341508e+015,
            4.463060e+015,
            3.330780e+015,
            9.724133e+014
        ], self.scenario.df.remaining_active_atoms, rtol=1e-6)

    def test_partial_activity(self):
        np.testing.assert_allclose([
            8.244806e-285,
            2.255364e+010,
            2.211493e-001,
            7.862791e-010,
            3.184579e-024,
            4.212641e-174
        ], self.scenario.df.partial_activity, rtol=1e-6)

    def test_watts(self):
        np.testing.assert_allclose([
            4.296413e-297,
            8.753265e-003,
            5.392346e-014,
            1.481078e-022,
            4.143464e-037,
            7.193802e-188
        ], self.scenario.df.watts, rtol=1e-6)

    def test_total_activity(self):
        df = self.scenario.df
        np.testing.assert_approx_equal(22522522523, df.partial_activity.sum(), significant=3)
        np.testing.assert_approx_equal(8.244805816629383e-285, df.partial_activity[df.parent == '190Pt'][0])

    def test_total_power(self):
        np.testing.assert_approx_equal(0.008753265223103896, self.scenario.df.watts.sum())


class HyperphysicsPoloniumAlphaDecayTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.scenario = System.load('212Po', model='induced-decay') \
            .hyperphysics(seconds=1, moles=1, active_fraction=1, isotopic_fraction=1)

    def test_nuclear_separation_radius_fm(self):
        np.testing.assert_allclose([9.014871826539528], self.scenario.df.nuclear_separation_radius_fm)

    def test_barrier_height(self):
        np.testing.assert_allclose([26.1967118938676], self.scenario.df.barrier_height_mev, rtol=1e-4)

    def test_lighter_ke(self):
        np.testing.assert_allclose([8.78], self.scenario.df.lighter_ke_mev, rtol=1e-3)

    def test_barrier_assault_frequency(self):
        np.testing.assert_allclose([1.142126641655716e21], self.scenario.df.barrier_assault_frequency, rtol=1e-3)

    def test_lighter_velocity_m_per_s(self):
        np.testing.assert_allclose([2.06e7], self.scenario.df.lighter_velocity_m_per_s, rtol=1e-3)

    def test_tunneling_probability(self):
        np.testing.assert_allclose([2.636693524272448e-15], self.scenario.df.tunneling_probability, rtol=1e-1)


class InducedFission106PdTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.scenario = System.load('106Pd',
            model='induced-fission',
            daughters=[(1, ('90Sr', '0')), (1, ('16O', '0'))],
        ).hyperphysics(
            seconds=1,
            moles=1,
            active_fraction=1,
            isotopic_fraction=1,
        )

    def test_daughters(self):
        np.testing.assert_equal(['90Sr, 16O'], self.scenario.df.daughters.values)

    def test_gamow_factor(self):
        np.testing.assert_allclose([529.519671], self.scenario.df.gamow_factor)


class HermesDecayTest(unittest.TestCase):

    ##
    # Others:
    # 218Rn, 219Rn, 220Rn, 222Rn, 226Rn
    # 221Fr
    # 223Ra, 224Ra, 226Ra
    # 225Ac, 227Ac
    # 227Th, 228Th, 229Th, 230Th, 232Th
    # 231Pa
    # 233U, 234U, 235U, 236U, 238U
    # 237Np
    # 238Pu, 239Pu, 240Pu, 244Pu
    # 244Cm, 245Cm, 248Cm
    # 249Cf, 252Cf
    #

    @classmethod
    def setUpClass(cls):
        cls.parameters = dict(seconds=1, moles=1, active_fraction=1, isotopic_fraction=1)

    def compare(self, isotope):
        system = System.load(isotope, model='induced-decay')
        hp_scenario = system.hyperphysics(**self.parameters)
        hm_scenario = system.hermes(**self.parameters)
        return hp_scenario, hm_scenario

    def test_209Bi(self):
        dfe, dfa = self.compare('209Bi')
        np.testing.assert_allclose(dfe.df.gamow_factor, dfa.df.gamow_factor + 8.5, rtol=1e-2)

    def test_211Bi(self):
        dfe, dfa = self.compare('211Bi')
        np.testing.assert_allclose(dfe.df.gamow_factor, dfa.df.gamow_factor + 4.0, rtol=1e-2)

    def test_212Bi(self):
        dfe, dfa = self.compare('212Bi')
        np.testing.assert_allclose(dfe.df.gamow_factor, dfa.df.gamow_factor + 4.4, rtol=1e-2)

    def test_213Bi(self):
        dfe, dfa = self.compare('213Bi')
        np.testing.assert_allclose(dfe.df.gamow_factor, dfa.df.gamow_factor + 4.6, rtol=1e-2)

    def test_210Po(self):
        dfe, dfa = self.compare('210Po')
        np.testing.assert_allclose(dfe.df.gamow_factor, dfa.df.gamow_factor + 5.2, rtol=1e-2)

    def test_211Po(self):
        dfe, dfa = self.compare('211Po')
        np.testing.assert_allclose(dfe.df.gamow_factor, dfa.df.gamow_factor + 3.6, rtol=1e-2)

    def test_212Po(self):
        dfe, dfa = self.compare('212Po')
        np.testing.assert_allclose(dfe.df.gamow_factor, dfa.df.gamow_factor + 2.9, rtol=1e-2)

    def test_214Po(self):
        dfe, dfa = self.compare('214Po')
        np.testing.assert_allclose(dfe.df.gamow_factor, dfa.df.gamow_factor + 3.4, rtol=1e-2)

    def test_215Po(self):
        dfe, dfa = self.compare('215Po')
        np.testing.assert_allclose(dfe.df.gamow_factor, dfa.df.gamow_factor + 3.6, rtol=1e-2)

    def test_216Po(self):
        dfe, dfa = self.compare('216Po')
        np.testing.assert_allclose(dfe.df.gamow_factor, dfa.df.gamow_factor + 3.9, rtol=1e-2)

    def test_218Po(self):
        dfe, dfa = self.compare('218Po')
        np.testing.assert_allclose(dfe.df.gamow_factor, dfa.df.gamow_factor + 4.6, rtol=1e-2)

    def test_215At(self):
        dfe, dfa = self.compare('215At')
        np.testing.assert_allclose(dfe.df.gamow_factor, dfa.df.gamow_factor + 3.3, rtol=1e-2)

    def test_217At(self):
        dfe, dfa = self.compare('217At')
        np.testing.assert_allclose(dfe.df.gamow_factor, dfa.df.gamow_factor + 3.9, rtol=1e-2)

    def test_218At(self):
        dfe, dfa = self.compare('218At')
        np.testing.assert_allclose(dfe.df.gamow_factor, dfa.df.gamow_factor + 4.0, rtol=1e-2)

    def test_241Am(self):
        dfe, dfa = self.compare('241Am')
        np.testing.assert_allclose(dfe.df.gamow_factor, dfa.df.gamow_factor + 5.8, rtol=1e-2)


class DependenceOfDecayConstantOnScreening(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.system = System.load(
            '190Pt',
            model='induced-fission',
            daughters=[(1, ('4He', '0')), (1, ('186Os', '0'))],
        )

    def test_relationship(self):
        decay_constants = []
        half_lives = []
        for screening in range(0, 10):
            scenario = self.system.hyperphysics(
                seconds=1,
                moles=1,
                active_fraction=1,
                isotopic_fraction=1,
                screening=screening,
            )
            constant = scenario.df.partial_decay_constant.values[0]
            decay_constants.append(constant)
            half_lives_seconds = 0.693 / constant
            half_lives_years = half_lives_seconds / 3.154e+7
            half_lives.append(half_lives_years)

        np.testing.assert_allclose([
            1.47e-20,
            8.09e-20,
            4.42e-19,
            2.40e-18,
            1.30e-17,
            7.04e-17,
            3.79e-16,
            2.03e-15,
            1.08e-14,
            5.76e-14,
        ], decay_constants, rtol=1e-1)

        np.testing.assert_allclose([
            1.48e+12,
            2.71e+11,
            4.96e+10,
            9.12e+09,
            1.68e+09,
            3.11e+08,
            5.79e+07,
            1.08e+07,
            2.02e+06,
            3.81e+05,
        ], half_lives, rtol=1e-1)
