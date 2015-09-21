import re


class HalfLife(object):

    def __init__(self, value, unit):
        self.value = value
        self.unit = unit

    def __str__(self):
        return '{} {}'.format(self.value, self.unit)


class NubaseRow(object):

    COLUMNS = (
        (  4, '_id'                      ),
        (  8, '_atomicNumber'            ),
        ( 20, '_nuclide'                 ),
        ( 40, '_massExcess'              ),
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

    @property
    def atomic_number(self):
        return int(re.search(r'\d+', self._raw['_atomicNumber']).group())

    @property
    def half_life(self):
        return HalfLife(self._raw['_halfLife'], self._raw['_halfLifeUnit'])
