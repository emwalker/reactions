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

    def test_case_1(self):
        r = NubaseRow(self.lines[0])
        self.assertEqual(0, r.atomic_number)
        self.assertEqual('613.9 s', str(r.half_life))
