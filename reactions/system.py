"""
Model a system of reactions specified by a reaction spec on the
command line.
"""
# pylint: disable=too-few-public-methods
from .nubase import parse_spec
from .combinations import Combinations
from .calculations import Decay


class Options:
    """Holds command-line options."""

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self.simple = kwargs.get('simple', False)
        self.references = not self.simple and kwargs.get('references')
        self.spins = not self.simple and kwargs.get('spins')
        self.gamow = kwargs.get('gamow')
        self.ascii = kwargs.get('ascii', False)
        self.notes = not self.simple


class System:
    """Models a system of reactions."""

    @classmethod
    def load(cls, string, **kwargs):
        """Factory method that returns an instance for a given set of options."""
        combinations = []
        system = filter(None, (rs.strip() for rs in string.split(',')))
        for spec in system:
            for reactants in parse_spec(spec, **kwargs):
                combination = Combinations.load(reactants=reactants, **kwargs)
                combinations.append(combination)
        return cls(combinations, **kwargs)

    def __init__(self, combinations, **kwargs):
        self.combinations = list(combinations)
        self._kwargs = kwargs

    def reactions(self):
        """Returns the various nuclear reactions that can result from the
        given parent nuclides, or that satisfy the input arguments.
        """
        for combination in self.combinations:
            for reaction in combination.reactions():
                yield combination, reaction

    def hyperphysics(self, **kwargs):
        """Carry out a set of decay calculations described in a Hyperphysics
        model.
        """
        return self._decay().hyperphysics(**kwargs)

    def hermes(self, **kwargs):
        """Cary out a set of decay calculations described by Hermes."""
        return self._decay().hermes(**kwargs)

    def _decay(self):
        return Decay.load(reactions=self.reactions())
