import unittest

from reactions.system import System, Options
from reactions.terminal import TerminalView, TerminalLine, StudiesTerminalView


class SystemTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.maxDiff = None

    def test_reactions_2(self):
        s = System.load('6Li+6Li')
        t = TerminalView(s)
        options = Options()
        self.assertTrue(any('2·6Li → 3·4He + 20899 keV' in l for l in t.lines(options)))

    def test_element_shorthand(self):
        s = System.load('H+Li')
        self.assertEqual(6, len(s.combinations))


class TestAscii(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.maxDiff = None

    def test_output(self):
        s = System.load('p+d', lower_bound=-3000)
        v = TerminalView(s)
        self.assertEqual(
            ['p + d => gamma + 3He + 5493 keV                         gamma, in nature',
             'p + d => p + d + 0 keV                                  in nature, n-transfer',
             'p + d => n + 2*p + -2225 keV                            ->B-, n']
        , v.lines(Options(ascii=True)))
