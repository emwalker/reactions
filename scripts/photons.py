import argparse
import os
import sys

from tabulate import tabulate

sys.path.append('lib')
from lenrmc import operator


class Report(object):
    def __init__(self, array):
        self.array = array


class TransitionReport(Report):
    def write_to(self, io):
        tables = self.array.partition_by(['channel', 'transition'])
        for table in tables:
            self._write_channel_to(io, table)

    def _write_channel_to(self, io, table):
        energy_mev, channel, transition = table.unique_values_at(['photon_energy', 'channel', 'transition'])
        energy_kev = energy_mev * 1000
        df_view = table.df[['material', 'material_thickness', 'escaping_photons']]
        tbl = tabulate(df_view, headers=df_view.columns.tolist())
        io.write("""
Escaping {:.0f} keV {} photons from {} transitions:

{}


""".format(energy_kev, channel, transition, tbl))



class App(object):
    def run(self):
        channels_by_layer = self._attenuators() * self._photon_counts() * self._channel_data()
        photons = operator.escaping_photons() * operator.transmitted_fraction() * channels_by_layer
        TransitionReport(photons).write_to(sys.stdout)

    def _channel_data(self):
        columns = ['transition', 'channel', 'photon_energy']
        table = operator.table(
            'System 1',
            values=[['58Ni(d,p)59Ni', 'bremsstrahlung', 511e-3]],
            columns=columns
        )
        system = table.merge(operator.materials, on=[operator.closest('photon_energy')])
        return system

    def _photon_counts(self):
        counts = operator.table(
            'Photon counts',
            values=[['58Ni(d,p)59Ni', 'bremsstrahlung', 5e10]],
            columns=['transition', 'channel', 'photon_count']
        )
        return counts

    def _attenuators(self):
        nickel_1cm = operator.attenuation('Nickel (1cm)', 'Nickel', '1cm')
        lead_1cm = operator.attenuation('Lead (1cm)', 'Lead', '1cm')
        air_20cm = operator.attenuation('Air', 'Air, Dry', '20cm')
        detector = operator.detector_attenuation('6150AD-b', diameter='10cm', distance='20cm', efficiency=0.26)
        #path1 = operator.attenuation('E-Cat (6150AD-b)', [nickel_1cm, lead_1cm, air_20cm, detector])
        array = operator.Array('First system', [
            operator.attenuation('Pyrex (1cm)', 'Pyrex', '1cm'),
            nickel_1cm,
            lead_1cm,
            operator.attenuation('Lead (2cm)', 'Lead', '2cm'),
            operator.attenuation('Lead (5cm)', 'Lead', '5cm'),
            air_20cm,
            #path1,
        ])
        return array


if __name__ == '__main__':
    App().run()
