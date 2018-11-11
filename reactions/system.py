"""
Model a system of reactions specified by a reaction spec on the
command line.
"""
# pylint: disable=too-few-public-methods, invalid-name
from .nubase import parse_spec
from .combinations import Combinations
from .calculations import Decay
from .views import SystemTerminalView


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

    def to_terminal(self, io, **kwargs):
        """Print the system to the provided io object."""
        SystemTerminalView(self, io, **kwargs).call()

    def _decay(self):
        return Decay.load(reactions=self.reactions())
