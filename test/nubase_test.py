import unittest
import os.path
import sys

from lenrmc.nubase import Nuclide


class NuclideTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        basepath = os.path.dirname(__file__)
        path = os.path.abspath(os.path.join(basepath, "../db/nubtab12.asc"))
        with open(path) as fh:
            cls.lines = list(fh)

    def test_basic_fields(self):
        n = Nuclide.load(line=self.lines[0])
        self.assertEqual(0, n.atomic_number)
        self.assertEqual('613.9 s', str(n.half_life))

    def test_column_widths(self):
        n = Nuclide.load(line=self.lines[16])
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
        }, n._row)

    def test_mass_number(self):
        n = Nuclide.load(line=self.lines[0])
        self.assertEqual(1, n.mass_number)

    def test_json(self):
        n = Nuclide.load(line=self.lines[0])
        self.assertEqual({
            'atomicNumber': 0,
            'massNumber':   1,
            'halfLife':     613.9,
        }, n.json())

    def test_isotopic_abundance_1(self):
        "There should be an abundance for a stable isotope."
        n = Nuclide.load(line=self.lines[1])
        self.assertEqual(99.9885, n.isotopic_abundance)
        n = Nuclide.load(line=self.lines[2])
        self.assertEqual(0.0115, n.isotopic_abundance)

    def test_isotopic_abundance_2(self):
        "The abundance should be zero for unstable isotopes."
        n = Nuclide.load(line=self.lines[3])
        self.assertEqual(0., n.isotopic_abundance)
