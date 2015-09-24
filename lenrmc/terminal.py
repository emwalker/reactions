from collections import defaultdict

from .studies import Studies


_studies = Studies.db()


class TerminalView(object):

    def __init__(self, system):
        self._system = system

    _desirable = {
        'n-transfer': 2,
        'stable':     3,
        '4He':        1,
        'ɣ':         -4,
        'n':         -5,
    }

    def _sort_key(self, reaction):
        desirable = sum(self._desirable.get(n, 0) for n in reaction.notes)
        return reaction.q_value_kev > 0, desirable, reaction.q_value_kev

    def _reactions(self):
        reactions = (TerminalLine(r) for r in self._system.reactions())
        return sorted(self._filter(reactions), key=self._sort_key, reverse=True)

    def _filter(self, reactions):
        return reactions

    def lines(self, **kwargs):
        refs = set()
        lines = []
        for r in self._reactions():
            line, _refs = r.terminal(**kwargs)
            lines.append(line)
            refs |= set(_refs)
        if refs and kwargs.get('references'):
            lines.extend([''] + sorted(refs))
        return lines


class StudiesTerminalView(TerminalView):

    def _sort_key(self, reaction):
        desirable = sum(self._desirable.get(n, 0) for n in reaction.notes)
        return reaction.agreement

    def _filter(self, reactions):
        return (r for r in reactions if r.agreement is not None)


class TerminalLine(object):

    def __init__(self, reaction, **kwargs):
        self._reaction = reaction
        self.q_value_kev = reaction.q_value_kev
        self.notes = reaction.notes
        self._lvalues = reaction._lvalues
        self._rvalues = reaction._rvalues
        self.references = []
        self.marks = []
        self._agreements = []
        self._add_references(self._lvalues, 'decrease', selective=True)
        self._add_references(self._rvalues, 'increase')
        self.agreement = sum(self._agreements) if self._agreements else None

    def terminal(self, **kwargs):
        kev = self.q_value_kev
        sign = '+' if kev >= 0 else '-'
        string = '{} → {} + {:.0f} keV'.format(
            self._fancy_side(self._lvalues),
            self._fancy_side(self._rvalues),
            kev,
        )
        string = '{:<55} {:<25}'.format(string, ', '.join(sorted(self.notes)))
        if kwargs.get('spins'):
            string = self._spin_and_parity(string, self._lvalues)
            string = self._spin_and_parity(string, self._rvalues)
        if kwargs.get('references'):
            string = self._add_marks(string)
        return string.strip(), self.references

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
        string += '   {}'.format(', '.join(self.marks))
        return string

    def _sort_key(self, pair):
        num, n = pair
        return n.mass_number, n.label

    def _fancy_side(self, side):
        isotopes = defaultdict(lambda: 0)
        for num, n in side:
            isotopes[n] += num
        values = [
            '{}·{}'.format(num, n.full_label) if num > 1 else n.full_label
            for num, n
            in sorted(((num, n) for n, num in isotopes.items()), key=self._sort_key)
        ]
        return ' + '.join(values)
