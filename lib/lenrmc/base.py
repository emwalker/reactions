import re
import pandas as pd
import numpy as np


class Channel(object):
    def __init__(self, df):
        self.df = df

    @property
    def escaped_photons(self):
        return self.df


class Generation(object):
    def __init__(self, engine, df):
        self.engine = engine
        df = self._escaping_photons(df)
        df = self._nist_energies(df)
        self.df = df

    def channel_for(self, transition, name):
        df = self.df[(self.df.transition == transition) & (self.df.channel == name)]
        return Channel(df)

    def _nist_energies(self, df):
        return df

    def _escaping_photons(self, df):
        df = self.engine.add_materials(df)
        def thickness(value):
            match = re.match(r'(\d+)\s*cm', value)
            if match:
                return float(match.group(1))
            raise Exception('do not recognize thickness: {}'.format(value))
        thickness_cm = df.material_thickness.apply(thickness)
        df['transmitted_fraction'] = np.exp(-df.mass_attenuation_coefficient * df.density * thickness_cm)
        df['escaping_photons'] = df.photon_count * df.transmitted_fraction
        return df


class Engine(object):
    @classmethod
    def from_csv(cls, materials_path, initpath):
        df_materials = pd.read_csv(materials_path, header=0)
        df_generation = pd.read_csv(initpath, header=0)
        return cls(df_materials, df_generation)

    def __init__(self, df_materials, df_generation):
        self.df_materials = df_materials
        self.generations = [Generation(self, df_generation)]

    def add_materials(self, df):
        df_after = df.merge(self.df_materials, on=['material'])
        if len(df) != len(df_after):
            raise Exception('one or more materials not found')
        return df_after
