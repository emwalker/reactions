"""
Carry out some calculations of interest that were demonstrated by Koonin.
"""
# pylint: disable=invalid-name, too-few-public-methods, no-self-use
import math
import numpy as np


class LogLambda:
    """Calculate the (log) rates for a given set of reactions."""

    _rate_constants = {
        'p+p': 8e-40,
        'd+p': 5.2e-22,
        'p+t': 4.8e-21,
        'd+d': 1.5e-16,
        'd+t': 1.3e-14,
    }

    # Hartrees
    _nucleon_masses = {
        'p': 3.449964e7,
        'd': 6.89463e7,
        't': 1.032447e8,
    }

    # Bohr radius
    bohr_radius = 0.53e-8

    def __init__(self, reactants, scale_factors):
        self.reactants = reactants
        self.scale_factors = scale_factors

    def rates(self):
        """What are the reaction rates for the given reactions?"""
        reaction = self._reaction(self.reactants)
        masses = [self._nucleon_masses[r] for r in self.reactants]
        A = self._rate_constants[reaction]
        nucleon_mass = self._nucleon_masses['p']
        reduced_mass = np.prod(masses) / sum(masses)
        mu_ratio = reduced_mass / nucleon_mass

        for scale in self.scale_factors:
            term1 = math.log10(A/self.bohr_radius**3)
            term2 = 3 * math.log10(mu_ratio)
            term3 = -79 * math.sqrt(mu_ratio) * math.sqrt(1./scale)
            rate = 6.5 + term1 + term2 + term3
            yield rate

    def _reaction(self, pair):
        return '+'.join(sorted(pair))
