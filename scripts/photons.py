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
        df_view = table.df[['material', 'material_thickness', 'photon_per_second', 'escaping_photons']]
        io.write("""
Escaping {:.0f} keV {} photons from {} transitions:

{}


""".format(energy_kev, channel, transition, repr(df_view)))


class LayerReport(Report):
    report_columns = [
        'transition',
        'channel',
        'relative_cross_section',
        'photon_energy',
        'photons_per_event',
        'events_per_transition_per_second',
        'photons_per_second',
        'escaping_photons',
        '6150ADb_detected_photons',
    ]

    scalar_columns = [
        'material',
        'material_thickness',
        'events_per_second',
        '6150ADb_detected_fraction'
    ]

    def write_to(self, io):
        tables = self.array.partition_by(['material', 'material_thickness'])
        for table in tables:
            self._write_channel_to(io, table)

    def _write_channel_to(self, io, table):
        material, thickness, events, detected = table.unique_values_at(self.scalar_columns)
        df_view = table.df[self.report_columns]
        io.write("""
Photons from a total of {} transitions per second, escaping through {} of {}:

{}


""".format(events, thickness, material, repr(df_view)))


class App(object):
    def run(self):
        channels_by_layer = self._layers() * operator.transitions * self._channel_data()
        fractions = operator.transmitted_fraction() * channels_by_layer
        escaping = operator.escaping_photons() * operator.photon_counts(7e14) * fractions
        detected = self._detected() * escaping
        for cls in (LayerReport,):
            cls(detected).write_to(sys.stdout)

    def _channel_data(self):
        system = operator.channels.merge(operator.materials, on=[operator.closest('photon_energy')])
        return system

    def _detected(self):
        return operator.detector('6150ADb', diameter='10cm', distance='20cm', efficiency=0.26)

    def _layers(self):
        nickel_1cm = operator.layer('Nickel (1cm)', 'Nickel', '1cm')
        lead_1cm = operator.layer('Lead (1cm)', 'Lead', '1cm')
        air_20cm = operator.layer('Air', 'Air, Dry', '20cm')
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
