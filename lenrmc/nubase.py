import re


class HalfLife(object):

    def __init__(self, value, unit):
        self.value = value
        self.unit = unit

    @property
    def seconds(self):
        if 's' == self.unit:
            return float(self.value)
        raise ValueError('do not know how to convert unit: {}'.format(self.unit))

    def __str__(self):
        return '{} {}'.format(self.value, self.unit)


class Nuclide(object):

    COLUMNS = (
        (  4, '_id'                      ),
        (  9, '_atomicNumber'            ),
        ( 19, '_nuclide'                 ),
        ( 39, '_massExcess'              ),
        ( 61, '_excitationEnergy'        ),
        ( 69, '_halfLife'                ),
        ( 71, '_halfLifeUnit'            ),
        ( 79, '_unknown'                 ),
        ( 93, '_spinAndParity'           ),
        ( 96, '_ensdfArchiveFileYear'    ),
        (105, '_reference'               ),
        (110, '_yearOfDiscovery'         ),
        ( -1, '_decayModesAndIntensities'),
    )

    def __init__(self, line):
        self.line = line
        self._raw = {}
        endcol_prev = 0
        for endcol, field in self.COLUMNS:
            text = line[endcol_prev:endcol]
            self._raw[field] = text.strip()
            endcol_prev = endcol
        self.atomic_number = int(re.search(r'\d+', self._raw['_atomicNumber']).group())
        self.mass_number = int(re.search(r'\d+', self._raw['_nuclide']).group())
        g = re.search(r'IS=([\d\.]+)', self._raw['_decayModesAndIntensities'])
        self.isotopic_abundance = float(g.group(1)) if g else 0.

    @property
    def half_life(self):
        return HalfLife(self._raw['_halfLife'], self._raw['_halfLifeUnit'])

    def json(self):
        return {
            'halfLife':     self.half_life.seconds,
            'atomicNumber': self.atomic_number,
            'massNumber':   self.mass_number,
        }

    def __iter__(self):
        return self.json().iteritems()
