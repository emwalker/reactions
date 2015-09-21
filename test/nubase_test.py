import unittest
import os.path
import sys

from lenrmc.nubase import NubaseRow


class NubaseTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        basepath = os.path.dirname(__file__)
        path = os.path.abspath(os.path.join(basepath, "../db/nubtab12.asc"))
        with open(path) as fh:
            cls.lines = list(fh)

    def test_basic_fieds(self):
        r = NubaseRow(self.lines[0])
        self.assertEqual(0, r.atomic_number)
        self.assertEqual('613.9 s', str(r.half_life))

    def test_column_widths(self):
        r = NubaseRow(self.lines[16])
        self.assertEqual({
            '_atomicNumber': '0038W',
            '_decayModesAndIntensities': 'IT=100',
            '_ensdfArchiveFileYear': '02',
            '_excitationEnergy': '3562.88    0.10',
            '_halfLife': '56',
            '_halfLifeUnit': 'as',
            '_id': '006',
            '_massExcess': '17649.76     0.10',
            '_nuclide': '6Lii',
            '_reference': '81Ro02e',
            '_spinAndParity': '0+      T=1',
            '_unknown': '14',
            '_yearOfDiscovery': '1981',
        }, r._raw)

    def test_mass_number(self):
        r = NubaseRow(self.lines[0])
        self.assertEqual(1, r.mass_number)

    def test_json(self):
        r = NubaseRow(self.lines[0])
        self.assertEqual({
            'atomicNumber': 0,
            'massNumber':   1,
            'halfLife':     613.9,
        }, r.json())
