# pylint: disable=missing-docstring, wrong-import-position
import argparse
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from reactions.system import System


class App:
    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self.system = System.load(self.kwargs['system_spec'], **self.kwargs)

    def call(self):
        if self.kwargs.get('decay_power'):
            self.print_decay_power()
            return
        self.print_possible_reactions()

    def print_decay_power(self):
        scenario = self.system.hyperphysics(**self.kwargs)
        if self.kwargs.get('format') == 'csv':
            scenario.to_csv(sys.stdout)
        else:
            scenario.to_terminal(sys.stdout)

    def print_possible_reactions(self):
        self.system.to_terminal(sys.stdout, **self.kwargs)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('system_spec', type=str)
    parser.add_argument('--lb', dest='lower_bound')
    parser.add_argument('--ub', dest='upper_bound')
    parser.add_argument('--spins', dest='spins', action='store_true')
    parser.add_argument('--references', dest='references', action='store_true')
    parser.add_argument('--view', type=str, dest='view')
    parser.add_argument('--studies', dest='studies', action='store_true')
    parser.add_argument('--model', dest='model')
    parser.add_argument('--unstable-parents', dest='unstable_parents', action='store_true')
    parser.add_argument('--excited', dest='excited', action='store_true')
    parser.add_argument('--ascii', dest='ascii', action='store_true')
    parser.add_argument('--parent-ub', dest='parent_ub', type=int)
    parser.add_argument('--simple', dest='simple', action='store_true')
    parser.add_argument('--gamow', dest='gamow', action='store_true')
    parser.add_argument('--decay-power', dest='decay_power', action='store_true')
    parser.add_argument('--screening', dest='screening', type=float)
    parser.add_argument('--moles', dest='moles', type=float)
    parser.add_argument('--seconds', dest='seconds', type=float)
    parser.add_argument('--active-fraction', dest='active_fraction', type=float)
    parser.add_argument('--format', dest='format')
    parser.add_argument('--daughter-count', dest='daughter_count')
    parser.set_defaults(
        active_fraction=1,
        ascii=False,
        daughter_count='',
        decay_power=False,
        excited=False,
        format=None,
        gamow=False,
        lower_bound=0,
        model='standard',
        moles=1,
        parent_ub=1000,
        references=False,
        screening=0,
        seconds=1,
        simple=False,
        spins=False,
        unstable_parents=False,
        upper_bound=500000,
        view='default',
    )
    return parser.parse_args()


if __name__ == '__main__':
    ARGS = parse_arguments()
    try:
        App(**vars(ARGS)).call()
    except KeyboardInterrupt:
        print('command canceled.')
