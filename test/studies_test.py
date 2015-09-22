import unittest

from lenrmc.studies import Studies


class StudiesTest(unittest.TestCase):

    def test_simple_case(self):
        results = Studies.db().isotopes(["6Li"])
        self.assertEqual({
            "6Li": [{"citation": "Lugano", "direction": "up"}]
        }, results.json)
