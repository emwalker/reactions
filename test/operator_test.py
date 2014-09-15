import unittest

import pandas as pd
from nose.tools import nottest

from lenrmc import operator


layers = operator.table(
    'Layers',
    values=[['Lead', '1cm']],
    columns=['material', 'material_thickness']
)

layer_details = layers.merge(operator.materials)

system1 = operator.table(
    'Photon energies x Photon counts',
    values=[[511e-3, 2], [67e-3, 2]],
    columns=['photon_energy', 'photons_per_second'],
)

system2 = operator.table(
    'Photon energies x Photon counts',
    values=[[8e-3, 2], [20e-3, 2]],
    columns=['photon_energy', 'photons_per_second'],
)


class TableTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.maxDiff = None

    def test_columns(self):
        self.assertEqual(layer_details.columns.tolist(), [
            'material',
            'material_thickness',
            'photon_energy',
            'mu_over_rho',
            'mu_en_over_rho',
            'density',
        ])

    def test_length(self):
        self.assertEqual(len(layer_details), 58)

    def test_values(self):
        self.assertEqual(layer_details.values[-5:].tolist(), [
            ['Lead', '1cm', 6.0, 0.043910000000000005, 0.02744, 11.34],
            ['Lead', '1cm', 8.0, 0.04675, 0.029889999999999996, 11.34],
            ['Lead', '1cm', 10.0, 0.04972, 0.031810000000000005, 11.34],
            ['Lead', '1cm', 15.0, 0.05658, 0.03478, 11.34],
            ['Lead', '1cm', 20.0, 0.062060000000000004, 0.03595, 11.34],
        ])

    def test_name(self):
        self.assertEqual(layer_details.name, 'Layers x Materials')

    def test_combination_with_operator(self):
        "An operator applied to table results in a table."
        result = operator.transmitted_fraction() * layer_details
        self.assertEqual(type(result), operator.Table)

    def test_does_not_apply_operator_from_right(self):
        "Currently a table cannot be applied to an operator."
        self.assertRaises(AttributeError, lambda: layer_details * operator.transmitted_fraction())


class EscapingPhotonsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.maxDiff = None
        cls.table = layer_details.merge(system1, on=[operator.closest('photon_energy')])

    def test_initial_columns(self):
        self.assertEqual(self.table.columns.tolist(), [
            'material',
            'material_thickness',
            'photon_energy',
            'mu_over_rho',
            'mu_en_over_rho',
            'density',
            'photons_per_second',
        ])

    def test_application_of_an_operation(self):
        self.assertTrue('transmitted_fraction' not in self.table.columns)
        table = operator.transmitted_fraction() * self.table
        self.assertTrue('transmitted_fraction' in table.columns)
        self.assertEqual(table.values.tolist(),
            [['Lead', '1cm', 0.08, 2.419, 1.916, 11.34, 2, 1.2208668215604055e-12],
             ['Lead', '1cm', 0.6, 0.1248, 0.06819, 11.34, 2, 0.24286926124353095]]
        )
        self.assertEqual(table.name, 'Transmitted fraction x Layers x Materials x Photon energies x Photon counts')

    def test_length(self):
        "There is only one result because there were only two rows in the system operator."
        self.assertEqual(len(self.table), 2)

    def test_escaping_photons(self):
        table = operator.escaping_photons() * operator.transmitted_fraction() * self.table
        view = table.df[['material', 'photon_energy', 'escaping_photons']]
        self.assertEqual(view.values.tolist(), [
            ['Lead', 0.08, 2.441733643120811e-12],
            ['Lead', 0.6, 0.4857385224870619]
        ])

    def test_name(self):
        self.assertEqual(self.table.name, 'Layers x Materials x Photon energies x Photon counts')


class ArrayTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ta0 = operator.Array('Trial 1', [system1, system2])
        cls.ta1 = cls.ta0.merge(layer_details, on=[operator.closest('photon_energy')])

    def test_count(self):
        self.assertEqual(len(self.ta1), 2)

    def test_columns(self):
        for table in self.ta1:
            self.assertEqual(table.columns.tolist(), [
                'material',
                'material_thickness',
                'photon_energy',
                'mu_over_rho',
                'mu_en_over_rho',
                'density',
                'photons_per_second',
            ])

    def test_initial_values(self):
        "Array.tolist() returns a multi-dimensional array with values for each operator."
        self.assertEqual(self.ta1.tolist(),
            [[['Lead', '1cm', 0.08, 2.419, 1.916, 11.34, 2],
              ['Lead', '1cm', 0.6, 0.1248, 0.06819, 11.34, 2]],
             [['Lead', '1cm', 0.008, 228.7, 220.7, 11.34, 2],
              ['Lead', '1cm', 0.02, 86.36, 68.99, 11.34, 2]]]
        )


if __name__ == '__main__':
  unittest.main()
