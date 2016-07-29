import re
import math
from collections import defaultdict

from .studies import Studies


_studies = Studies.db()


class TerminalView(object):

    def __init__(self, system):
        self._system = system

    _kwargs = {'selective': True}

    def reactions(self, cls):
        reactions = (cls(c, r, **self._kwargs) for c, r in self._system.reactions())
        return sorted(self._filter(reactions), key=lambda l: l.sort_key, reverse=True)

    def _filter(self, reactions):
        return reactions

    def lines(self, options):
        refs = set()
        lines = []
        cls = AsciiTerminalLine if options.ascii else UnicodeTerminalLine
        for r in self.reactions(cls):
            line, _refs = r.terminal(options)
            lines.append(line)
            refs |= set(_refs)
        if refs and options.references:
            lines.extend([''] + sorted(refs))
        return lines


class StudiesTerminalView(TerminalView):

    _not_observed = {'ɣ', 'n'}
    _kwargs = {}

    def _sort_key(self, reaction):
        length = len(reaction._agreements)
        sign = 1 if reaction.agreement > 0 else -1
        return reaction.agreement, sign * length

    def _filter(self, reactions):
        for r in reactions:
            if r.agreement is None:
                continue
            if any(n in self._not_observed for n in r.notes):
                continue
            yield r


class TerminalLine(object):

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
        self._agreements = []
        self._add_references(self._lvalues, 'decrease', **kwargs)
        self._add_references(self._rvalues, 'increase')
        # Cases where a daughter was found in a study
        self.agreement = sum(self._agreements) if self._agreements else None

    def _spin_and_parity(self, string, values):
        spins_and_parities = (n.spin_and_parity for num, n in sorted(values, key=self._sort_key))
        spins_and_parities = filter(None, spins_and_parities)
        string = '{} {:<20}'.format(string, ', '.join(sorted(spins_and_parities)))
        return string

    def _add_references(self, values, expected, **kwargs):
        selective = kwargs.get('selective')
        for result in _studies.isotopes(n.label for num, n in values):
            agreement, mark = result.reference_mark(expected)
            self._agreements.append(1 if agreement else -1)
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
        num, n = pair
        return n.mass_number, n.label

    def _fancy_side(self, delim, side):
        isotopes = defaultdict(lambda: 0)
        for num, n in side:
            isotopes[n] += num
        values = []
        nuclides = sorted(((num, n) for n, num in isotopes.items()), key=self._sort_key)
        for num, n in nuclides:
            label = self.format(n.full_label)
            string = self._multi_daughter_template.format(num, label) if num > 1 else label
            values.append(string)
        return ' {} '.format(delim).join(values)

    def _add_gamow(self, string):
        G = self.reaction._gamow
        if G is None:
            return string
        return '{} [{:.0f}]'.format(string, G)

    def terminal(self, options):
        kev = self.q_value_kev
        sign = '+' if kev >= 0 else '-'
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

    _multi_daughter_template = '{}·{}'
    _reaction_template = '{} → {} + {:.0f} keV'

    def format(self, string):
        return string


class AsciiTerminalLine(TerminalLine):

    _multi_daughter_template = '{}*{}'
    _reaction_template = '{} => {} + {:.0f} keV'

    _translated_patterns = [
        ('→',   '->'),
        ('β',   'B'),
        ('ε',   'EC'),
        ('α',   'A'),
        ('ɣ',   'gamma'),
        ('ν',   'neutrino'),
        ('✗',   'x'),
        ('✓'    'a'),
    ]

    def format(self, string):
        for before, after in self._translated_patterns:
            string = re.sub(before, after, string)
        return string
