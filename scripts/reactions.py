import argparse


from lenrmc.nubase import Combinations


class App(object):

    def __init__(self, **kwargs):
        self.kwargs = kwargs

    def run(self):
        kwargs = dict(self.kwargs)
        if 'reactants' in kwargs:
            del kwargs['reactants']
        c = Combinations.load(reactants=self._reactants, **kwargs)
        for reaction in c.terminal():
            print(reaction)

    @property
    def _reactants(self):
        return [
            (1, r)
            for r in filter(None, self.kwargs.get('reactants', '').split())
        ]


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reactants', type=str)
    parser.add_argument('--lb', dest='lower_bound', help='lower bound in keV')
    parser.set_defaults(lower_bound=-1000)
    return parser.parse_args()


if '__main__' == __name__:
    opts = parse_arguments()
    App(**vars(opts)).run()
