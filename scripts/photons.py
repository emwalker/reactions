# coding: utf-8
import argparse
import os
import sys

import pandas as pd

from lenrmc import operator


pd.set_option('display.float_format', lambda f: '{:.2e}'.format(f))


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
        df_view = table.df[['material', 'material_thickness', 'photon_count', 'escaping_photons']]
        io.write("""
Escaping {:.0f} keV {} photons from {} transitions:

{}


""".format(energy_kev, channel, transition, repr(df_view)))


class LayerReport(Report):
    def write_to(self, io):
        tables = self.array.partition_by(['material', 'material_thickness'])
        for table in tables:
            self._write_channel_to(io, table)

    def _write_channel_to(self, io, table):
        material, thickness = table.unique_values_at(['material', 'material_thickness'])
        df_view = table.df[['transition', 'channel', 'photon_energy', 'photons_per_event', 'photon_count', 'escaping_photons']]
        io.write("""
Photons escaping through {} of {}:

{}


""".format(thickness, material, repr(df_view)))


class App(object):
    def run(self):
        channels_by_layer = self._layers() * operator.transitions * self._channel_data()
        fractions = operator.transmitted_fraction() * channels_by_layer
        escaping = operator.escaping_photons() * operator.photon_counts(7e14) * fractions
        for cls in (TransitionReport, LayerReport):
            cls(escaping).write_to(sys.stdout)

    def _channel_data(self):
        system = operator.channels.merge(operator.materials, on=[operator.closest('photon_energy')])
        return system

    def _layers(self):
        nickel_1cm = operator.layer('Nickel (1cm)', 'Nickel', '1cm')
        lead_1cm = operator.layer('Lead (1cm)', 'Lead', '1cm')
        air_20cm = operator.layer('Air', 'Air, Dry', '20cm')
        detector = operator.detector('6150AD-b', diameter='10cm', distance='20cm', efficiency=0.26)
        #path1 = operator.attenuation('E-Cat (6150AD-b)', [nickel_1cm, lead_1cm, air_20cm, detector])
        array = operator.Array('First system', [
            operator.layer('Pyrex (1cm)', 'Pyrex', '1cm'),
            nickel_1cm,
            lead_1cm,
            operator.layer('Lead (2cm)', 'Lead', '2cm'),
            operator.layer('Lead (5cm)', 'Lead', '5cm'),
            air_20cm,
            #path1,
        ])
        return array


if __name__ == '__main__':
    App().run()