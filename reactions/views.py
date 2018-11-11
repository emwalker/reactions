"""
Print various objects to an io stream, e.g., a terminal, formatted in a human-
readable way.
"""
#pylint: disable=invalid-name, too-few-public-methods

import pandas as pd

from .terminal import TerminalView, StudiesTerminalView, Options


class DecayTerminalView:
    """Serialize a decay scenario to something that can be printed to a
    terminal.
    """

    def __init__(self, scenario, io, **kwargs):
        self.scenario = scenario
        self.io = io
        self.kwargs = kwargs

    def call(self):
        """Print to the io object."""
        self.io.write('')
        self.io.write('At second:      {}'.format(self.kwargs.get('seconds')))
        self.io.write('Starting moles: {}'.format(self.kwargs.get('moles')))
        self.io.write('Activity:       {:.2e}'.format(self.scenario.activity()))
        self.io.write('Watts:          {:.2e}'.format(self.scenario.power().watts))
        self.io.write('')

        df = self.scenario.df[[
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

        if df.empty:
            self.io.write('No active isotopes.')
        else:
            with pd.option_context('display.max_rows', 999, 'display.max_columns', 10):
                df = df.dropna().sort_values(['watts', 'gamow_factor'], ascending=[0, 1])
                self.io.write(df.to_string() + '\n')
        self.io.write('')


class SystemTerminalView:
    """Print out a system of reactions to a terminal."""

    def __init__(self, system, io, **kwargs):
        self.system = system
        self.io = io
        self.kwargs = kwargs

    def call(self):
        """Print to the io object."""
        options = Options(**self.kwargs)
        for line in self._view_cls(self.system).lines(options):
            self.io.write(line + '\n')

    @property
    def _view_cls(self):
        if self.kwargs.get('view') == 'studies' or self.kwargs.get('studies'):
            return StudiesTerminalView
        return TerminalView
