from .nubase import parse_spec
from .combinations import Combinations


class System(object):

    @classmethod
    def parse(cls, string, **kwargs):
        system = filter(None, (rs.strip() for rs in string.split(',')))
        combinations = []
        for spec in system:
            for reactants in parse_spec(spec):
                c = Combinations.load(reactants=reactants, **kwargs)
                combinations.append(c)
        return cls(combinations, **kwargs)

    def __init__(self, combinations, **kwargs):
        self._combinations = list(combinations)
        self._kwargs = kwargs

    def reactions(self):
        for c in self._combinations:
            yield from c.reactions()
