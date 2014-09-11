# coding: UTF-8
import unittest
import math

import numpy as np
from lenrmc.base import *
from lenrmc.layer import DetectorLayer


class ChannelsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        model = Model.from_csv('fixtures/materials.csv', 'fixtures/channels.csv')
        nickel_1cm = model.layer_for('Nickel', '1cm')
        lead_1cm = model.layer_for('Lead', '1cm')
        air_20cm = model.layer_for('Air, Dry', '20cm')
        cls.detector = DetectorLayer('6150AD-b', diameter='10cm', distance='20cm', efficiency=0.26)
        system1 = model.layer_of([nickel_1cm, lead_1cm, air_20cm], name='E-Cat')
        system2 = model.layer_of([nickel_1cm, lead_1cm, air_20cm, cls.detector], name='E-Cat (6150AD-b)')
        cls.materials = [
            model.layer_for('Pyrex', '1cm'),
            nickel_1cm,
            lead_1cm,
            model.layer_for('Lead', '2cm'),
            model.layer_for('Lead', '5cm'),
            air_20cm,
            system1,
            system2,
        ]
        cls.model = model

    def test_detector_transmitted_fraction(self):
        self.assertEqual(self.detector.solid_angle, math.pi * 0.05971499970933625)
        self.assertEqual(self.detector.relative_coverage, 0.014928749927334062)
        self.assertEqual(self.detector.transmitted_fraction, 0.003881474981106856)

    def test_escaping_15_kev_xrays(self):
        channel = self.model.generations[-1].channel_for('58Ni(d,p)59Ni', 'bremsstrahlung')
        values = channel.escaping_photons(self.materials, 6.3979e+16).values.tolist()
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

    def test_escaping_ep_annihilation_photons(self):
        "Test for slowly β+ decaying 59Ni → 59Co, a decay that has a half-life of ~ 7e4 years."
        channel = self.model.generations[-1].channel_for('58Ni(d,p)59Ni', 'β-β+ annihilation photons')
        values = channel.escaping_photons(self.materials, 41).values.tolist()
        self.assertEqual(values, [
            ['Pyrex',            '1cm',            34.27415483196742],
            ['Nickel',           '1cm',            20.204807006274784],
            ['Lead',             '1cm',            9.957639710984768],
            ['Lead',             '2cm',            2.418404600336118],
            ['Lead',             '5cm',            0.03464548262924468],
            ['Air, Dry',         '20cm',           40.92048574868956],
            ['E-Cat',            '1cm/1cm/20cm',   4.897609807430786],
            ['E-Cat (6150AD-b)', '1cm/1cm/20cm/D', 0.019009949934766165],
        ])

    def test_compound_layer(self):
        "A compound layer should keep track of the photon transmission of each successive layer."
        channel = self.model.generations[-1].channel_for('58Ni(d,p)59Ni', 'β-β+ annihilation photons')
        compound = self.materials[-1]
        multi_df = channel.escaping_photons([compound], 41)
        self.assertEqual(multi_df.last.values.tolist(), [
            ['E-Cat (6150AD-b)', '1cm',            20.204807006274784],
            ['E-Cat (6150AD-b)', '1cm/1cm',        4.907126551182075],
            ['E-Cat (6150AD-b)', '1cm/1cm/20cm',   4.897609807430786],
            ['E-Cat (6150AD-b)', '1cm/1cm/20cm/D', 0.019009949934766165]
        ])


if __name__ == '__main__':
    unittest.main()
