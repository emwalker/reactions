import argparse

from lenrmc.system import System, Options
from lenrmc.terminal import TerminalView, StudiesTerminalView


class App(object):

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        if 'studies' == self.kwargs.get('view') or kwargs.get('studies'):
            self.view_cls = StudiesTerminalView
        else:
            self.view_cls = TerminalView

    def run(self):
        s = System.load(self.kwargs['system'], **self.kwargs)
        options = Options(**self.kwargs)
        for line in self.view_cls(s).lines(options):
            print(line)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('system', type=str)
    parser.add_argument('--lb', dest='lower_bound')
    parser.add_argument('--ub', dest='upper_bound')
    parser.add_argument('--spins', dest='spins', action='store_true')
    parser.add_argument('--references', dest='references', action='store_true')
    parser.add_argument('--view', type=str, dest='view')
    parser.add_argument('--studies', dest='studies', action='store_true')
    parser.add_argument('--model', dest='model')
    parser.add_argument('--unstable', dest='unstable', action='store_true')
    parser.add_argument('--excited', dest='excited', action='store_true')
    parser.add_argument('--ascii', dest='ascii', action='store_true')
    parser.add_argument('--parent-ub', dest='parent_ub', type=int)
    parser.add_argument('--simple', dest='simple', action='store_true')
    parser.add_argument('--gamow', dest='gamow', action='store_true')
    parser.add_argument('--daughter-count', dest='daughter_count')
    parser.set_defaults(
        lower_bound    = 0,
        upper_bound    = 500000,
        spins          = False,
        references     = False,
        view           = 'default',
        model          = 'standard',
        unstable       = False,
        ascii          = False,
        excited        = False,
        simple         = False,
        gamow          = False,
        parent_ub      = 1000,
        daughter_count = '',
    )
    return parser.parse_args()


if '__main__' == __name__:
    opts = parse_arguments()
    App(**vars(opts)).run()
