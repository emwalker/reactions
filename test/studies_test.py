# pylint: disable=missing-docstring
import unittest

from reactions.studies import Studies


class StudiesTest(unittest.TestCase):
    def test_simple_case(self):
        results = Studies.data().isotopes(['6Li'])
        self.assertEqual([
            {'shortDescription': '2015 Lugano E-Cat test by Levi et al.',
             'label': '6Li',
             'change': 'increase'}
        ], results.json)
