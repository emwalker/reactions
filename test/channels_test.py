# coding: UTF-8
import unittest
from lenrmc.base import *


class ChannelsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        model = Model.from_csv('fixtures/materials.csv', 'fixtures/channels.csv')
        cls.materials = [
            model.material_for('Pyrex', '1cm'),
            model.material_for('Nickel', '1cm'),
            model.material_for('Lead', '1cm'),
            model.material_for('Lead', '2cm'),
            model.material_for('Lead', '5cm'),
            model.material_for('Air, Dry', '100cm'),
        ]
        cls.model = model

    def test_escaping_15_kev_xrays(self):
        channel = self.model.generations[-1].channel_for('58Ni(d,p)59Ni', 'Ni K x-rays')
        values = channel.escaping_photons(self.materials, 6.3979e+16).values.tolist()
        self.assertEqual(values, [
            ['Pyrex',    '1cm',   566884433791.0281],
            ['Nickel',   '1cm',   7.30679905354258e-258],
            ['Lead',     '1cm',   0.0],
            ['Lead',     '2cm',   0.0],
            ['Lead',     '5cm',   0.0],
            ['Air, Dry', '100cm', 5.267115195919887e+16],
        ])

    def test_escaping_ep_annihilation_photons(self):
        "Test for slowly β+ decaying 59Ni → 59Co, a decay that has a half-life of ~ 7e4 years."
        channel = self.model.generations[-1].channel_for('58Ni(d,p)59Ni', 'β-β+ annihilation photons')
        values = channel.escaping_photons(self.materials, 41).values.tolist()
        self.assertEqual(values, [
            ['Pyrex',    '1cm',   34.27415483196742],
            ['Nickel',   '1cm',   20.204807006274784],
            ['Lead',     '1cm',   9.957639710984768],
            ['Lead',     '2cm',   2.418404600336118],
            ['Lead',     '5cm',   0.03464548262924468],
            ['Air, Dry', '100cm', 40.60396783279807],
        ])


if __name__ == '__main__':
    unittest.main()
