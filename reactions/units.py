import math


class Energy(object):

    @classmethod
    def load(cls, **kwargs):
        if 'kev' in kwargs:
            return cls(kwargs['kev'])
        if 'mev' in kwargs:
            return cls(1e3 * kwargs['mev'])
        raise ValueError('do not know how to load energy')

    def __init__(self, energy_kev):
        self.kev = energy_kev

    @property
    def mev(self):
        return 1e-3 * self.kev

    @property
    def joules(self):
        return self.kev / 6.2415091e15


class Power(object):

    @classmethod
    def load(cls, **kwargs):
        return cls(kwargs['watts'])

    def __init__(self, watts):
        self.watts = watts


class Distance(object):

    @classmethod
    def load(cls, **kwargs):
        return cls(kwargs['fermis'])

    def __init__(self, fermis):
        self.fermis = fermis


class HalfLife(object):

    def __init__(self, value, unit):
        self.value = value
        self.unit = unit

    @property
    def seconds(self):
        if math.inf == self.value:
            return self.value
        if 's' == self.unit:
            return float(self.value)
        if 'd' == self.unit:
            return 86400 * float(self.value)
        if 'y' == self.unit:
            return 3.154e+7 * float(self.value)
        if 'ky' == self.unit:
            return 3.154e+10 * float(self.value)
        if 'Py' == self.unit:
            return 3.154e+23 * float(self.value)
        if 'ms' == self.unit:
            return 0.001 * float(self.value)
        raise ValueError('do not know how to convert unit: {}'.format(self.unit))

    def __str__(self):
        return '{} {}'.format(self.value, self.unit)
