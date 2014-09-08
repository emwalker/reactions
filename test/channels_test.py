import unittest
from lenrmc.base import *


class ChannelsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Model.from_csv('fixtures/materials.csv', 'fixtures/channels.csv')
        cls.channel = cls.model.generations[-1].channel_for('58Ni(d,p)59Ni', 'Ni K x-rays')

    def test_escaping_photons(self):
        values = self.channel.escaping_photons(6.3979e+16).values.tolist()
        self.assertEqual(values, [
            ['Pyrex',  '1cm', 566884433791.0281],
            ['Nickel', '1cm', 7.30679905354258e-258],
            ['Lead',   '1cm', 0.0],
            ['Lead',   '2cm', 0.0],
            ['Lead',   '5cm', 0.0],
        ])


if __name__ == '__main__':
    unittest.main()
