import argparse
import os
import sys

from tabulate import tabulate

sys.path.append('lib')
from lenrmc.base import Model
from lenrmc.layer import DetectorLayer


class Report(object):
    def __init__(self, model, layers):
        self.model = model
        self.layers = layers
        self.channels = self.model.generations[-1].channels


class TransitionReport(Report):
    def write_to(self, io):
        for channel in self.channels:
            self._write_channel_to(io, channel)

    def _write_channel_to(self, io, channel):
        df_full = channel.escaping_photons(self.layers, 6.3979e+16).df
        df_view = df_full[['material', 'material_thickness', 'escaping_photons']]
        energy = channel.photon_energy.tolist()[0] * 1e3
        table = tabulate(df_view, headers=df_view.columns.tolist())
        io.write("""
Escaping {:.0f} keV {} photons from {} transitions:

{}


""".format(energy, channel.name, channel.transition, table))


class LayerReport(Report):
    def write_to(self, io):
        for layer in self.layers:
            self._write_layer_to(io, layer)

    def _write_layer_to(self, io, layer):
        df_full = layer.escaping_photons(self.model, 6.3979e+16).df
        df_view = df_full[['channel', 'photon_energy', 'escaping_photons']]
        table = tabulate(df_view, headers=df_view.columns.tolist())
        io.write("""
Escaping {:.0f} keV {} photons for {} of {}:

{}


""".format(energy, channel.name, layer.thickness, layer.name, table))


class App(object):
    def run(self):
        model = Model()
        report = self._report_for(model)
        report.write_to(sys.stdout)

    def _report_for(self, model):
        nickel_1cm = model.layer_for('Nickel', '1cm')
        lead_1cm = model.layer_for('Lead', '1cm')
        air_20cm = model.layer_for('Air, Dry', '20cm')
        detector = DetectorLayer('6150AD-b', diameter='10cm', distance='20cm', efficiency=0.26)
        path1 = model.layer_of([nickel_1cm, lead_1cm, air_20cm, detector], name='E-Cat (6150AD-b)')
        layers = [
            model.layer_for('Pyrex', '1cm'),
            nickel_1cm,
            lead_1cm,
            model.layer_for('Lead', '2cm'),
            model.layer_for('Lead', '5cm'),
            air_20cm,
            path1,
        ]
        return TransitionReport(model, layers)


if __name__ == '__main__':
    App().run()
