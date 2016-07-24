import math
from collections import defaultdict

import numpy as np
import pandas as pd
import scipy.constants as cs

from .constants import FINE_STRUCTURE_CONSTANT_MEV_FM, HBAR_MEV_S
from .units import Energy, Power, HalfLife, Distance


class AlphaCalculationMixin(object):

    @classmethod
    def load(cls, components, q_value, **kwargs):
        if components is None:
            return None
        parent, daughters = components
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


class GeigerNuttal(AlphaCalculationMixin):

    def value(self):
        return -46.83 + 1.454 * self.larger.atomic_number / math.sqrt(self.q_value.mev)


class Gamow2(AlphaCalculationMixin):
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


class GamowSuppressionFactor(AlphaCalculationMixin):
    """Gamow suppression factor in log10 units

    From Hermes: https://www.lenr-forum.com/forum/index.php/Thread/3434-Document-Isotopic-Composition
    -of-Rossi-Fuel-Sample-Unverified/?postID=29085#post29085.

    Estimate distance between 2 spherical nuclei when they touch
    Estimate distance when supplied energy has overcome Coulomb barrier
    if r >= 1 there is no Coulomb barrier to cross!

    - from Hermes's comments
    """

    def value(self):
        A  = self.larger.mass_number
        Z  = self.larger.atomic_number
        A4 = self.smaller.mass_number
        Z4 = self.smaller.atomic_number
        Q  = self.q_value.mev
        if Q < 0:
            return math.nan
        # Distances in fm
        rs = 1.1 * (pow(A, .333333) + pow(A4, .333333))
        rc = float(Z) * Z4 * 1.43998 / Q
        r  = rs / rc
        G  = 0 if r >= 1 else math.acos(math.sqrt(r)) - math.sqrt(r * (1. - r))
        m  = (float(A) * A4) / (A + A4)
        return 0.2708122 * Z * Z4 * G * math.sqrt(m / Q)


class IsotopicAlphaDecayCalculation(AlphaCalculationMixin):
    """From http://hyperphysics.phy-astr.gsu.edu/hbase/nuclear/alpdec.html
    """

    speed_of_light = 3 * math.pow(10, 8)
    hbarc = 197.33

    def __init__(self, parent, daughters, q_value, **kwargs):
        self.parent = parent
        self.daughters = daughters
        self.q_value = q_value
        self.kwargs = kwargs
        self.parent_isotopic_abundance = parent.isotopic_abundance
        self.screening = kwargs.get('screening') or 0
        self.smaller, self.larger = daughters
        self.A4, self.Z4 = self.smaller.mass_number, self.smaller.atomic_number
        self.A,  self.Z  = self.larger.mass_number,  self.larger.atomic_number
        self.screened_Z = self.Z - self.screening
        self.alpha_mass = self.smaller.mass.mev
        # Ea = Q / (1 + m/M)
        self.alpha_energy = self.q_value.mev / (1 + self.alpha_mass / self.larger.mass.mev)
        # Units in fm
        self.nuclear_separation = 1.2 * (math.pow(self.A4, 1./3) - (-1) * math.pow(self.A, 1./3))

    @property
    def barrier_height(self):
        "Units in MeV"
        return 2 * self.screened_Z * 1.44 / self.nuclear_separation

    @property
    def alpha_velocity(self):
        "Units in m/s"
        return math.sqrt(2 * self.alpha_energy / self.alpha_mass) * self.speed_of_light

    @property
    def barrier_assault_frequency(self):
        "Units in s^-1"
        return self.alpha_velocity * math.pow(10, 15) / (2 * self.nuclear_separation)

    @property
    def gamow_factor(self):
        x = self.alpha_energy / self.barrier_height
        ph = math.sqrt(2 * self.alpha_mass / ((self.hbarc ** 2) * self.alpha_energy))
        return ph * 2 * self.screened_Z * 1.44 * (math.acos(math.sqrt(x)) - math.sqrt(x * (1 - x)))

    @property
    def tunneling_probability(self):
        return math.exp(-2 * self.gamow_factor)

    @property
    def decay_constant(self):
        return self.tunneling_probability * self.barrier_assault_frequency

    @property
    def half_life(self):
        seconds = 0.693 / self.decay_constant
        return HalfLife(seconds, 's')


class DecayScenario(object):

    avogadros_number, _, _ = cs.physical_constants['Avogadro constant']

    def __init__(self, base_df, **kwargs):
        self.base_df = base_df.copy()
        self.kwargs = kwargs
        self.df = self._calculate(base_df, kwargs)

    def _calculate(self, df, kwargs):
        elapsed = kwargs['seconds']
        df = df.copy()
        df = self._gamow_factor(df)
        df['starting_moles'] = kwargs['moles'] * (kwargs.get('isotopic_fraction') or df.isotopic_fraction)
        df['active_fraction'] = kwargs.get('active_fraction') or 1
        df['starting_active_moles'] = df.starting_moles * df.active_fraction
        df['starting_active_atoms'] = df.starting_active_moles * self.avogadros_number
        df['remaining_active_atoms'] = df.starting_active_atoms * np.exp(-df.decay_constant * elapsed)
        df['activity'] = df.decay_constant * df.remaining_active_atoms
        df['watts'] = df.activity * df.deposited_q_value_joules
        return df

    def recalculate(self, **kwargs):
        merged = {**self.kwargs, **kwargs}
        if merged == self.kwargs:
            return self
        return DecayScenario(self.base_df, **merged)

    def activity(self, **kwargs):
        return self.recalculate(**kwargs).df.activity.sum()

    def power(self, **kwargs):
        watts = self.recalculate(**kwargs).df.watts.sum()
        return Power.load(watts=watts)

    def remaining_active_atoms(self, **kwargs):
        return self.recalculate(**kwargs).df.remaining_active_atoms.sum()

    def _gamow_factor(self, df):
        return df


class AlphaDecay(object):

    initial_column_names = [
        'atomic_number',
        'isotope',
        'screening',
        'gamow_factor',
        'decay_constant',
        'isotopic_abundance',
        'deposited_q_value_kev',
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
        self.decays = defaultdict(list)
        for d in (r.alpha_decay(**kwargs) for r in self.reactions):
            if d is None:
                continue
            self.decays[d.parent.atomic_number].append(d)
        self.kwargs = kwargs
        self.df = self._initial_dataframe()

    def _initial_dataframe(self):
        rows = []
        for atomic_number, decays in self.decays.items():
            for d in decays:
                rows.append([
                    atomic_number,
                    d.parent.label,
                    d.screening,
                    d.gamow_factor,
                    d.decay_constant,
                    d.parent.isotopic_abundance,
                    d.q_value.kev,
                    d.q_value.joules,
                ])
        df = pd.DataFrame(rows, columns=self.initial_column_names)
        df['isotopic_fraction'] = df.isotopic_abundance / 100.
        return df

    def scenario(self, **kwargs):
        merged = {**self.kwargs, **kwargs}
        return DecayScenario(self.df, **merged)
