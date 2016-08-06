from .nubase import parse_spec
from .combinations import Combinations
from .calculations import Decay


class Options(object):

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self.simple = kwargs.get('simple', False)
        self.references = not self.simple and kwargs.get('references')
        self.spins = not self.simple and kwargs.get('spins')
        self.gamow = kwargs.get('gamow')
        self.ascii = kwargs.get('ascii', False)
        self.notes = not self.simple


class System(object):

    @classmethod
    def load(cls, string, **kwargs):
        combinations = []
        system = filter(None, (rs.strip() for rs in string.split(',')))
        for spec in system:
            for reactants in parse_spec(spec, **kwargs):
                c = Combinations.load(reactants=reactants, **kwargs)
                combinations.append(c)
        return cls(combinations, **kwargs)

    def __init__(self, combinations, **kwargs):
        self.combinations = list(combinations)
        self._kwargs = kwargs

    def reactions(self):
        for c in self.combinations:
            for r in c.reactions():
                yield c, r

    def decay(self, **kwargs):
        return Decay.load(reactions=self.reactions()).scenario(**kwargs)
