import math
from collections import defaultdict

import numpy as np
import pandas as pd
import scipy.constants as cs

from .constants import FINE_STRUCTURE_CONSTANT_MEV_FM, HBAR_MEV_S
from .units import Energy, Power, HalfLife, Distance


##
# Choices of constant
# r0 -- Fermi model nuclear radius
#
#   1.1 fm
#   https://goo.gl/1OVppU
#
#   1.14 fm
#   http://www.iaea.org/inis/collection/NCLCollectionStore/_Public/26/031/26031671.pdf
#
#   1.2 fm
#   http://hyperphysics.phy-astr.gsu.edu/hbase/nuclear/nucuni.html#c4
#   http://inside.mines.edu/~fsarazin/phgn310/PDFs/NukeLect1.pdf
#   http://www.physics.umd.edu/courses/Phys741/xji/chapter7.pdf
#   http://www.jpoffline.com/physics_docs/y3s6/nuclearphysics_ln.pdf
#   https://www.southampton.ac.uk/~ab1u06//teaching/phys3002/course/03_diffraction.pdf
#
#   1.25 fm
#   https://en.wikipedia.org/wiki/Atomic_nucleus#Nuclear_models
#
#   1.57 fm
#   https://goo.gl/nNFgSA
#


class FragmentCalculationMixin(object):

    @classmethod
    def load(cls, components, q_value, **kwargs):
        if components is None:
            return None
        parent, daughters = components
        if 'daughters' in kwargs:
            del kwargs['daughters']
        return cls(parent, daughters, q_value, **kwargs)

    def __init__(self, parent, daughters, q_value, **kwargs):
        self.parent = parent
        self.smaller, self.larger = daughters
        self.q_value = q_value
        self.kwargs = kwargs


class CoulombBarrier(object):

    def __init__(self, n0, n1):
        self.n0 = n0
        self.n1 = n1
        self._base = FINE_STRUCTURE_CONSTANT_MEV_FM * n0.atomic_number * n1.atomic_number

    def height(self, radius):
        height = self._base / float(radius.fermis)
        return Energy.load(mev=height)

    def width(self, q_value):
        width = self._base / q_value.mev
        return Distance.load(fermis=width)


class ReactionEnergy(object):

    def __init__(self, reaction):
        self.reaction = reaction
        self.value = Energy.load(kev=self._kev())

    def _kev(self):
        lvalues = sum(num * i.mass_excess_kev for num, i in self.reaction._lvalues)
        rvalues = sum(num * i.mass_excess_kev for num, i in self.reaction.rvalues)
        return lvalues - rvalues


class GeigerNuttal(FragmentCalculationMixin):

    def value(self):
        return -46.83 + 1.454 * self.larger.atomic_number / math.sqrt(self.q_value.mev)


class Gamow2(FragmentCalculationMixin):
    """Gamow factor for alpha particle tunneling.

    Assumes one of the daughters is a heavy nucleus and the other an alpha particle.
    """

    def value(self):
        Q = self.q_value.mev
        x = Q / CoulombBarrier(self.smaller, self.larger).width(self.q_value).fermis
        t0 = math.sqrt((2 * self.smaller.mass.mev)/(HBAR_MEV_S**2 * Q))
        t1 = self.smaller.atomic_number * self.larger.atomic_number * FINE_STRUCTURE_CONSTANT_MEV_FM
        t2 = math.acos(math.sqrt(x)) - math.sqrt(x * (1 - x))
        return t0 * t1 * t2


class GamowSuppressionFactor(FragmentCalculationMixin):
    """Gamow suppression factor in log10 units

    From Hermes: https://www.lenr-forum.com/forum/index.php/Thread/3434-Document-Isotopic-Composition
    -of-Rossi-Fuel-Sample-Unverified/?postID=29085#post29085.

    Estimate distance between 2 spherical nuclei when they touch
    Estimate distance when supplied energy has overcome Coulomb barrier
    if r >= 1 there is no Coulomb barrier to cross!

    - from Hermes's comments
    """

    def value(self):
        screening = self.kwargs.get('screening') or 0
        A  = self.larger.mass_number
        Z  = self.larger.atomic_number - screening
        A4 = self.smaller.mass_number
        Z4 = self.smaller.atomic_number
        Q  = self.q_value.mev
        if Q <= 0:
            return math.nan
        # Distances in fm
        # Changed from 1.1 to 1.2, and from .333333 to 1/3
        rs = 1.2 * (pow(A, 1/3) + pow(A4, 1/3))
        rc = float(Z) * Z4 * 1.43998 / Q
        r  = 1 if rc <= 0 else rs / rc
        G  = 0 if r >= 1 else math.acos(math.sqrt(r)) - math.sqrt(r * (1. - r))
        m  = (float(A) * A4) / (A + A4)
        return 0.2708122 * Z * Z4 * G * math.sqrt(m / Q)


class IsotopicDecay(FragmentCalculationMixin):
    """From http://hyperphysics.phy-astr.gsu.edu/hbase/nuclear/alpdec.html
    """

    def __init__(self, parent, daughters, q_value, **kwargs):
        self.parent = parent
        self.parent_z = parent.atomic_number
        self.parent_a = parent.mass_number
        self.daughters = daughters
        self.q_value = q_value
        self.kwargs = kwargs
        self.smaller, self.larger = daughters
        self.data = self._data()

    def _data(self):
        heavier_daughter_z = self.larger.atomic_number
        return {
            'parent': self.parent.label,
            'daughters': ', '.join([self.larger.label, self.smaller.label]),
            'parent_z': self.parent_z,
            'parent_a': self.parent_a,
            'lighter_daughter_z': self.smaller.atomic_number,
            'lighter_daughter_a': self.smaller.mass_number,
            'heavier_daughter_z': heavier_daughter_z,
            'heavier_daughter_a': self.larger.mass_number,
            'lighter_mass_mev': self.smaller.mass.mev,
            'heavier_daughter_mass_mev': self.larger.mass.mev,
            'q_value_mev': self.q_value.mev,
            'isotopic_abundance': self.parent.isotopic_abundance,
            'deposited_q_value_joules': self.q_value.joules,
        }

    def __getitem__(self, key):
        return self.data[key]


class DecayScenario(object):

    speed_of_light, _, _ = cs.physical_constants['speed of light in vacuum']
    # hbar * c, in units of MeV.fm
    hbarc = 197.327
    # e^2 / 4pi, in units of MeV.fm
    # Google: (electron charge)^2 / (4 * pi * epsilon_0) -> 2.30707751e-28 m^3 kg / s^2
    # Wolfram Alpha: 2.30707751e-28 m^3 kg / s^2 in MeV fm -> 1.4399645 MeV fm
    # Changed to 1.43998 to match Hermes's value
    e2_4pi = 1.43998
    avogadros_number, _, _ = cs.physical_constants['Avogadro constant']

    def __init__(self, base_df, reactions, **kwargs):
        self.base_df = base_df.copy()
        self.reactions = reactions
        self.kwargs = kwargs
        self.df = self.calculate(base_df, kwargs)

    def to_csv(self, io):
        self.df.to_csv(io, index=False)

    def to_string(self):
        return self.df.to_string()

    def calculate(self, df, kwargs):
        df = df.copy()
        df = self.calculate_preliminaries(df, kwargs)
        df = self.calculate_gamow_factor(df, kwargs)
        df = self.calculate_decay_constant(df, kwargs)
        df = self.calculate_products(df, kwargs)
        return df

    def calculate_gamow_factor(self, df, kwargs):
        raise NotImplementedError

    def recalculate(self, **kwargs):
        merged = {**self.kwargs, **kwargs}
        if merged == self.kwargs:
            return self
        return self.__class__(self.base_df, self.reactions, **merged)

    def activity(self, **kwargs):
        return self.recalculate(**kwargs).df.partial_activity.sum()

    def power(self, **kwargs):
        watts = self.recalculate(**kwargs).df.watts.sum()
        return Power.load(watts=watts)

    def remaining_active_atoms(self, **kwargs):
        return self.recalculate(**kwargs).df.remaining_active_atoms.sum()

    def calculate_preliminaries(self, df, kwargs):
        df['screening'] = kwargs.get('screening') or 0
        df['screened_heavier_daughter_z'] = df.heavier_daughter_z - df.screening
        df['lighter_ke_mev'] = df.q_value_mev / (1 + df.lighter_mass_mev / df.heavier_daughter_mass_mev)
        df['nuclear_separation_radius_fm'] = 1.2 * (np.power(df.lighter_daughter_a, 1./3) + np.power(df.heavier_daughter_a, 1./3))
        df['lighter_velocity_m_per_s'] = np.sqrt(2 * df.lighter_ke_mev / df.lighter_mass_mev) * self.speed_of_light
        df['barrier_assault_frequency'] = df.lighter_velocity_m_per_s * math.pow(10, 15) / (2 * df.nuclear_separation_radius_fm)
        # rc = float(Z) * Z4 * 1.43998 / Q
        df['coulomb_barrier_radius_fm'] = df.screened_heavier_daughter_z * df.lighter_daughter_z * self.e2_4pi / df.q_value_mev
        # r  = rs / rc
        df['radius_ratio'] = df.nuclear_separation_radius_fm / df.coulomb_barrier_radius_fm
        return df

    def calculate_decay_constant(self, df, kwargs):
        df['tunneling_probability'] = np.exp(-2 * df.gamow_factor)
        df['partial_decay_constant'] = df.tunneling_probability * df.barrier_assault_frequency
        df['isotope_decay_constant'] = df.groupby(['parent_a', 'parent_z']).partial_decay_constant.transform(np.sum)
        df['partial_half_life'] = np.where(df.partial_decay_constant > 0, math.log(2) / df.partial_decay_constant, math.inf)
        return df

    def calculate_products(self, df, kwargs):
        elapsed = kwargs['seconds']
        df['starting_moles'] = kwargs['moles'] * (kwargs.get('isotopic_fraction') or df.parent_fraction)
        df['active_fraction'] = kwargs.get('active_fraction') or 1
        df['starting_active_moles'] = df.starting_moles * df.active_fraction
        df['starting_active_atoms'] = df.starting_active_moles * self.avogadros_number
        df['remaining_active_atoms'] = df.starting_active_atoms * np.exp(-df.isotope_decay_constant * elapsed)
        df['partial_activity'] = df.partial_decay_constant * df.remaining_active_atoms
        df['watts'] = df.partial_activity * df.deposited_q_value_joules
        return df


class HyperphysicsDecayScenario(DecayScenario):

    ##
    # From http://hyperphysics.phy-astr.gsu.edu/hbase/nuclear/alpdec.html.
    #

    def calculate_gamow_factor(self, df, kwargs):
        df['barrier_height_mev'] = 2 * df.screened_heavier_daughter_z * self.e2_4pi / df.nuclear_separation_radius_fm
        r = df.lighter_ke_mev / df.barrier_height_mev
        ph = math.sqrt(2) * np.sqrt(df.lighter_mass_mev / df.lighter_ke_mev)
        G = np.where(r >= 1, 0, np.arccos(np.sqrt(r)) - np.sqrt(r * (1 - r)))
        df['gamow_factor'] = self.e2_4pi / self.hbarc * df.screened_heavier_daughter_z * df.lighter_daughter_z * G * ph
        return df


class HermesDecayScenario(DecayScenario):

    def calculate_gamow_factor(self, df, kwargs):
        # m  = (float(A) * A4) / (A + A4)
        m = (df.heavier_daughter_a * df.lighter_daughter_a) / (df.heavier_daughter_a + df.lighter_daughter_a)
        # G  = 0 if r >= 1 else math.acos(math.sqrt(r)) - math.sqrt(r * (1. - r))
        G = np.where(df.radius_ratio >= 1, 0, np.arccos(np.sqrt(df.radius_ratio)) - np.sqrt(df.radius_ratio * (1 - df.radius_ratio)))
        # return 0.2708122 * Z * Z4 * G * math.sqrt(m / Q)
        df['gamow_factor'] = 0.2708122 * df.screened_heavier_daughter_z * df.lighter_daughter_z * G * np.sqrt(m / df.q_value_mev)
        return df


class Decay(object):

    initial_column_names = [
        'parent_z',
        'parent_a',
        'parent',
        'daughters',
        'heavier_daughter_z',
        'lighter_daughter_a',
        'heavier_daughter_a',
        'lighter_daughter_z',
        'lighter_mass_mev',
        'heavier_daughter_mass_mev',
        'q_value_mev',
        'isotopic_abundance',
        'deposited_q_value_joules',
    ]

    @classmethod
    def load(cls, **kwargs):
        copy = kwargs.copy()
        reactions = copy['reactions']
        del copy['reactions']
        return cls(reactions, **copy)

    def __init__(self, reactions, **kwargs):
        self.reactions = list(reactions)
        self.decays = []
        for d in (r.decay(**kwargs) for c, r in self.reactions):
            if d is None:
                continue
            self.decays.append(d)
        self.kwargs = kwargs
        self.df = self._initial_dataframe()

    def _initial_dataframe(self):
        rows = []
        for d in self.decays:
            rows.append([d[c] for c in self.initial_column_names])
        if len(rows) < 1:
            df = pd.DataFrame(columns=self.initial_column_names)
        else:
            df = pd.DataFrame(rows, columns=self.initial_column_names)
        df['parent_fraction'] = df.isotopic_abundance / 100.
        return df

    def hyperphysics(self, **kwargs):
        merged = {**self.kwargs, **kwargs}
        return HyperphysicsDecayScenario(self.df, self.reactions, **merged)

    def hermes(self, **kwargs):
        merged = {**self.kwargs, **kwargs}
        return HermesDecayScenario(self.df, self.reactions, **merged)
