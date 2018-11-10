# pylint: disable=missing-docstring
import unittest

import numpy as np

from reactions.units import Energy


class EnergyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.q = Energy.load(mev=3)

    def test_kev(self):
        np.testing.assert_approx_equal(3000.0, self.q.kev)

    def test_mev(self):
        np.testing.assert_approx_equal(3, self.q.mev)

    def test_joules(self):
        np.testing.assert_approx_equal(4.806529882332463e-13, self.q.joules)
