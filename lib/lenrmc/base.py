import re
import os

import pandas as pd
import numpy as np

from .layer import Layer
from .helpers import MultiDataFrame


class Channel(object):
    def __init__(self, model, df):
        self.model = model
        self.df = df
        self.photon_energy = df.photon_energy

    def escaping_photons(self, materials, photon_count):
        df_original = self.df.copy()
        dfs = self._transmitted_fraction(materials, df_original)
        views = []
        for df in dfs:
            df['escaping_photons'] = df.transmitted_fraction * photon_count
            view = df[['material', 'material_thickness', 'escaping_photons']]
            views.append(view)
        return MultiDataFrame(views)

    def _transmitted_fraction(self, materials, df_original):
        array = []
        for material in materials:
            multi_df = material.transmitted_fraction_for(df_original.photon_energy)
            array.append(multi_df.dataframes)
        combineds = []
        for dfs in array:
            df = df_original.copy()
            df_fractions = pd.concat(dfs)
            df_fractions['join_key'] = 1
            df['join_key'] = 1
            df = df.merge(df_fractions, on='join_key')
            combineds.append(df)
        return combineds


class Generation(object):
    def __init__(self, model, df):
        self.model = model
        self.df = df

    def channel_for(self, transition, name):
        df = self.df[(self.df.transition == transition) & (self.df.channel == name)]
        return Channel(self, df)


class Model(object):
    materials_path = os.path.join(os.path.dirname(__file__), '../../db/materials.csv')
    channels_path  = os.path.join(os.path.dirname(__file__), '../../db/channels.csv')

    @classmethod
    def from_csv(cls, materials_path, channels_path):
        return cls(materials_path, channels_path)

    def __init__(self, materials_path=None, channels_path=None):
        self.df_materials = pd.read_csv(materials_path or self.materials_path, header=0)
        df_generation = pd.read_csv(channels_path or self.channels_path, header=0)
        self.generations = [Generation(self, df_generation)]

    def layer_for(self, name, thickness):
        df = self.df_materials[self.df_materials.material == name]
        if len(df) < 1:
            raise Exception('no material found: {}'.format(name))
        return Layer(name, thickness, df)

    def layer_of(self, *args, **kwargs):
        return Layer.of(*args, **kwargs)
