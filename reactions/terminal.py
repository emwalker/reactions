"""
Various helper classes for serializing various reactions in a form that can
be printed to a terminal.
"""
# pylint: disable=no-self-use, too-few-public-methods, too-many-instance-attributes
# pylint: disable=no-member
import re
from collections import defaultdict

from .studies import Studies


STUDIES = Studies.data()


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


class TerminalView:
    """Base class responsible for formatting reactions for printing out
    to the console.
    """

    def __init__(self, system):
        self._system = system

    _kwargs = {'selective': True}

    def reactions(self, cls):
        """Return a sorted list consisting of one line per reaction in the system."""
        reactions = (cls(c, r, **self._kwargs) for c, r in self._system.reactions())
        return sorted(self._filter(reactions), key=lambda l: l.sort_key, reverse=True)

    def _filter(self, reactions):
        return reactions

    def _line_class(self, options):
        return AsciiTerminalLine if options.ascii else UnicodeTerminalLine

    def lines(self, options):
        """Return all of the lies to be printed out to the terminal, together with references
        if any.
        """
        refs = set()
        lines = []
        for reaction in self.reactions(self._line_class(options)):
            line, _refs = reaction.terminal(options)
            lines.append(line)
            refs |= set(_refs)
        if refs and options.references:
            lines.extend([''] + sorted(refs))
        return lines


class StudiesTerminalView(TerminalView):
    """Terminal view that includes references to relevant studies."""

    _not_observed = {'ɣ', 'n'}
    _kwargs = {}

    def _sort_key(self, reaction):
        length = len(reaction.agreements)
        sign = 1 if reaction.agreement > 0 else -1
        return reaction.agreement, sign * length

    def _filter(self, reactions):
        for reaction in reactions:
            if reaction.agreement is None:
                continue
            if any(n in self._not_observed for n in reaction.notes):
                continue
            yield reaction


class TerminalLine:
    """Abstract base class representing a single line of terminal output."""

    _notes_template = '{:<55} {:<25}'

    def __init__(self, combinations, reaction, **kwargs):
        self.combinations = combinations
        self.reaction = reaction
        self.sort_key = combinations.sort_key(reaction)
        self.q_value_kev = reaction.q_value.kev
        self.notes = [self.format(s) for s in reaction.notes]
        self._lvalues = reaction.lvalues
        self._rvalues = reaction.rvalues
        self.references = []
        self.marks = []
        self.agreements = []
        self._add_references(self._lvalues, 'decrease', **kwargs)
        self._add_references(self._rvalues, 'increase')
        # Cases where a daughter was found in a study
        self.agreement = sum(self.agreements) if self.agreements else None

    def _spin_and_parity(self, string, values):
        spins_and_parities = (n.spin_and_parity for num, n in sorted(values, key=self._sort_key))
        spins_and_parities = filter(None, spins_and_parities)
        string = '{} {:<20}'.format(string, ', '.join(sorted(spins_and_parities)))
        return string

    def _add_references(self, values, expected, **kwargs):
        selective = kwargs.get('selective')
        for result in STUDIES.isotopes(n.label for num, n in values):
            agreement, mark = result.reference_mark(expected)
            self.agreements.append(1 if agreement else -1)
            if selective and agreement:
                continue
            self.marks.append(mark)
            self.references.append(result.reference_line)

    def _add_marks(self, string):
        string += '   {}'.format(', '.join(self._format_mark(m) for m in self.marks))
        return string

    def _format_mark(self, mark):
        return '{:>13}'.format(self.format(mark))

    def _sort_key(self, pair):
        _, nuclide = pair
        return nuclide.mass_number, nuclide.label

    def _fancy_side(self, delim, side):
        isotopes = defaultdict(lambda: 0)
        for num, nuclide in side:
            isotopes[nuclide] += num
        values = []
        nuclides = sorted(((num, n) for n, num in isotopes.items()), key=self._sort_key)
        for num, nuclide in nuclides:
            label = self.format(nuclide.full_label)
            string = self._multi_daughter_template.format(num, label) if num > 1 else label
            values.append(string)
        return ' {} '.format(delim).join(values)

    def _add_gamow(self, string):
        gamow = self.reaction.gamow_value
        if gamow is None:
            return string
        return '{} [{:.0f}]'.format(string, gamow)

    def terminal(self, options):
        """Returns a single line of terminal output together with any relevant
        academic references.
        """
        kev = self.q_value_kev
        string = self._reaction_template.format(
            self._fancy_side(self.reaction.lvalue_delim, self._lvalues),
            self._fancy_side(self.reaction.rvalue_delim, self._rvalues),
            kev,
        )
        if options.gamow:
            string = self._add_gamow(string)
        if options.notes:
            string = self._notes_template.format(string, ', '.join(sorted(self.notes)))
        if options.spins:
            string = self._spin_and_parity(string, self._lvalues)
            string = self._spin_and_parity(string, self._rvalues)
        if options.references:
            string = self._add_marks(string)
        return string.strip(), self.references


class UnicodeTerminalLine(TerminalLine):
    """A terminal line for a terminal that supports unicode characters."""

    _multi_daughter_template = '{}·{}'
    _reaction_template = '{} → {} + {:.0f} keV'

    def format(self, string):
        """No-op method that returns the unmodified string."""
        return string


class AsciiTerminalLine(TerminalLine):
    """A terminal line that includes only ASCII characters."""

    _multi_daughter_template = '{}*{}'
    _reaction_template = '{} => {} + {:.0f} keV'

    _translated_patterns = [
        ('→', '->'),
        ('β', 'B'),
        ('ε', 'EC'),
        ('α', 'A'),
        ('ɣ', 'gamma'),
        ('ν', 'neutrino'),
        ('✗', 'x'),
        ('✓' 'a'),
    ]

    def format(self, string):
        """Replaces unicode characters with ASCII representations."""
        for before, after in self._translated_patterns:
            string = re.sub(before, after, string)
        return string
