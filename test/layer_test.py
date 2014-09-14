# coding: UTF-8
import unittest

from nose.tools import nottest

from lenrmc.base import *
from lenrmc.layer import DetectorLayer


class LayerTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        model = Model()
        nickel_1cm = model.layer_for('Nickel', '1cm')
        lead_1cm = model.layer_for('Lead', '1cm')
        air_20cm = model.layer_for('Air, Dry', '20cm')
        cls.detector = DetectorLayer('6150AD-b', diameter='10cm', distance='20cm', efficiency=0.26)
        path1 = model.layer_of([nickel_1cm, lead_1cm, air_20cm, cls.detector], name='E-Cat (6150AD-b)')
        cls.layers = [
            model.layer_for('Pyrex', '1cm'),
            nickel_1cm,
            air_20cm,
            path1,
        ]
        cls.model = model

    @nottest
    def test_photons_escaping_through_layer(self):
        layer = self.layers[0]
        counts = [
            ['58Ni(d,p)59Ni', 1],
            ['58Ni(d,p)59Ni', 1],
            ['60Ni(d,p)61Ni', 1],
            ['61Ni(d,p)62Ni', 1],
            ['62Ni(d,p)63Ni', 1],
            ['64Ni(d,p)65Ni', 1],
            ['58Ni(p,ɣ)59Cu', 1],
            ['60Ni(p,ɣ)61Cu', 1],
            ['60Ni(p,ɣ)61Cu', 1],
            ['61Ni(p,α)58Co', 1],
            ['61Ni(p,α)58Co', 1],
            ['61Ni(p,ɣ)62Cu', 1],
            ['62Ni(p,ɣ)63Cu', 1],
            ['64Ni(p,ɣ)65Cu', 1],
            ['d(p,ɣ)3He',     1],
            ['d(d,ɣ)4He',     1],
        ]
        df = pd.DataFrame(counts, columns=['transition', 'photon_count'])
        values = layer.escaping_photons(self.model, df).values.tolist()
        self.assertEqual(values, [
            ['Pyrex',            '1cm',            566884433791.0281],
            ['Nickel',           '1cm',            7.30679905354258e-258],
            ['Lead',             '1cm',            0.0],
            ['Lead',             '2cm',            0.0],
            ['Lead',             '5cm',            0.0],
            ['Air, Dry',         '20cm',           6.153816211202043e+16],
            ['E-Cat',            '1cm/1cm/20cm',   0.0],
            ['E-Cat (6150AD-b)', '1cm/1cm/20cm/D', 0.0],
        ])


if __name__ == '__main__':
    unittest.main()
