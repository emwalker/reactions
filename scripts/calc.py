# pylint: disable=missing-docstring, wrong-import-position
import argparse
import sys
import os
import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from reactions.system import System, Options
from reactions.terminal import TerminalView, StudiesTerminalView


class App:
    def __init__(self, **kwargs):
        self.kwargs = kwargs
        if self.kwargs.get('view') == 'studies' or kwargs.get('studies'):
            self.view_cls = StudiesTerminalView
        else:
            self.view_cls = TerminalView

    def run(self):
        system = System.load(self.kwargs['system'], **self.kwargs)
        if self.kwargs.get('decay_power'):
            scenario = system.hyperphysics(**self.kwargs)
            if self.kwargs.get('format') == 'csv':
                scenario.to_csv(sys.stdout)
            else:
                self.print_scenario(scenario)
        else:
            options = Options(**self.kwargs)
            for line in self.view_cls(system).lines(options):
                print(line)

    def print_scenario(self, scenario):
        print()
        print('At second:      {}'.format(self.kwargs.get('seconds')))
        print('Starting moles: {}'.format(self.kwargs.get('moles')))
        print('Activity:       {:.2e}'.format(scenario.activity()))
        print('Watts:          {:.2e}'.format(scenario.power().watts))
        print()

        frame = scenario.df[[
            'parent',
            'daughters',
            'parent_fraction',
            'q_value_mev',
            'starting_moles',
            'gamow_factor',
            'partial_half_life',
            'partial_activity',
            'watts',
        ]]

        if frame.empty:
            print('No active isotopes.')
        else:
            with pd.option_context('display.max_rows', 999, 'display.max_columns', 10):
                frame = frame.dropna().sort(['watts', 'gamow_factor'], ascending=[0, 1])
                print(frame)
        print()


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
        lower_bound=0,
        upper_bound=500000,
        spins=False,
        references=False,
        view='default',
        model='standard',
        unstable_parents=False,
        ascii=False,
        excited=False,
        simple=False,
        gamow=False,
        parent_ub=1000,
        daughter_count='',
        decay_power=False,
        screening=0,
        seconds=1,
        moles=1,
        active_fraction=1,
        format=None,
    )
    return parser.parse_args()


if __name__ == '__main__':
    ARGS = parse_arguments()
    try:
        App(**vars(ARGS)).run()
    except KeyboardInterrupt:
        print('command canceled.')
