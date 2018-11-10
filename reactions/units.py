"""
Helper classes to facilitate converting to desired units.
"""
# pylint: disable=too-many-return-statements, too-few-public-methods
import math


class Energy:
    """Model a quantity of energy."""

    @classmethod
    def load(cls, **kwargs):
        """Factory method that produces an instance for a given set of arguments."""
        if 'kev' in kwargs:
            return cls(kwargs['kev'])
        if 'mev' in kwargs:
            return cls(1e3 * kwargs['mev'])
        raise ValueError('do not know how to load energy')

    def __init__(self, energy_kev):
        self.kev = energy_kev

    @property
    def mev(self):
        """The energy expressed in megaelectronvolts."""
        return 1e-3 * self.kev

    @property
    def joules(self):
        """The energy expressed in Joules."""
        return self.kev / 6.2415091e15


class Power:
    """Model a quantity of power."""

    @classmethod
    def load(cls, **kwargs):
        """No-op factory method that returns an instance."""
        return cls(kwargs['watts'])

    def __init__(self, watts):
        self.watts = watts


class Distance:
    """Model a distance."""

    @classmethod
    def load(cls, **kwargs):
        """No-op factory method that returns an instance."""
        return cls(kwargs['fermis'])

    def __init__(self, fermis):
        self.fermis = fermis


class HalfLife:
    """Model the half-life of a radionuclide."""

    def __init__(self, value, unit):
        self.value = value
        self.unit = unit

    @property
    def seconds(self):
        """Convert the half-life to seconds."""
        if math.inf == self.value:
            return self.value
        if self.unit == 'ms':
            return 0.001 * float(self.value)
        if self.unit == 's':
            return float(self.value)
        if self.unit == 'd':
            return 86400 * float(self.value)
        if self.unit == 'y':
            return 3.154e+7 * float(self.value)
        if self.unit == 'ky':
            return 3.154e+10 * float(self.value)
        if self.unit == 'Py':
            return 3.154e+23 * float(self.value)
        raise ValueError('do not know how to convert unit: {}'.format(self.unit))

    def __str__(self):
        return '{} {}'.format(self.value, self.unit)
