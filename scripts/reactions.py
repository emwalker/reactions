import argparse

from lenrmc.system import System
from lenrmc.terminal import TerminalView, StudiesTerminalView


class App(object):

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        if 'studies' == self.kwargs.get('view') or kwargs.get('studies'):
            self.view_cls = StudiesTerminalView
        else:
            self.view_cls = TerminalView

    def run(self):
        s = System.parse(self.kwargs['system'], **self.kwargs)
        for line in self.view_cls(s).lines(**self.kwargs):
            print(line)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('system', type=str)
    parser.add_argument('--lb', dest='lower_bound')
    parser.add_argument('--spins', dest='spins', action='store_true')
    parser.add_argument('--references', dest='references', action='store_true')
    parser.add_argument('--view', type=str, dest='view')
    parser.add_argument('--studies', dest='studies', action='store_true')
    parser.add_argument('--model', dest='model')
    parser.add_argument('--unstable', dest='unstable', action='store_true')
    parser.set_defaults(
        lower_bound = 0,
        spins       = False,
        references  = True,
        view        = 'default',
        model       = 'regular',
        unstable    = False,
    )
    return parser.parse_args()


if '__main__' == __name__:
    opts = parse_arguments()
    App(**vars(opts)).run()
