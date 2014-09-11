import re

import pandas as pd
import numpy as np

from .layer import Layer


class Channel(object):
    def __init__(self, model, df):
        self.model = model
        self.df = df
        self.photon_energy = df.photon_energy

    def escaping_photons(self, materials, photon_count):
        df = self.df.copy()
        df = self._transmitted_fraction(materials, df)
        df['escaping_photons'] = df.transmitted_fraction * photon_count
        return df[['material', 'material_thickness', 'escaping_photons']]

    def _transmitted_fraction(self, materials, df):
        array = map(lambda m: m.transmitted_fraction_for(df.photon_energy), materials)
        df_fractions = pd.concat(array)
        df_fractions['join_key'] = 1
        df['join_key'] = 1
        df = df.merge(df_fractions, on='join_key')
        return df


class Generation(object):
    def __init__(self, model, df):
        self.model = model
        self.df = df

    def channel_for(self, transition, name):
        df = self.df[(self.df.transition == transition) & (self.df.channel == name)]
        return Channel(self, df)


class Model(object):
    @classmethod
    def from_csv(cls, materials_path, initpath):
        df_materials = pd.read_csv(materials_path, header=0)
        df_generation = pd.read_csv(initpath, header=0)
        return cls(df_materials, df_generation)

    def __init__(self, df_materials, df_generation):
        self.df_materials = df_materials
        self.generations = [Generation(self, df_generation)]

    def layer_for(self, name, thickness):
        df = self.df_materials[self.df_materials.material == name]
        if len(df) < 1:
            raise Exception('no material found: {}'.format(name))
        return Layer(name, thickness, df)

    def layer_of(self, *args, **kwargs):
        return Layer.of(*args, **kwargs)
