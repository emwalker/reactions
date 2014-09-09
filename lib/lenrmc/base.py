import re
import pandas as pd
import numpy as np


class Channel(object):
    def __init__(self, model, df):
        self.model = model
        self.df = df

    def escaping_photons(self, materials, photon_count):
        df = self.df.copy()
        df = self._transmitted_fraction(materials, df)
        df['escaping_photons'] = df.transmitted_fraction * photon_count
        return df[['material', 'material_thickness', 'escaping_photons']]

    def _transmitted_fraction(self, materials, df):
        combined = pd.concat(map(lambda m: m.df, materials))
        array = map(lambda m: m.transmitted_fraction_for(df.photon_energy), materials)
        df_fractions = pd.concat(array)
        df_fractions['join_key'] = 1
        df['join_key'] = 1
        df = df.merge(df_fractions, on='join_key')
        return df


class Material(object):
    def __init__(self, name, thickness, df):
        self.df = df.copy()
        self.density = df.density.iloc[0]
        self.df['name'] = self.name = name
        self.df['material_thickness'] = self.thickness = thickness
        self.df = df.sort('photon_energy')
        self.photon_energy = df.photon_energy.values
        self.mu_over_rho = df.mu_over_rho.values

    def transmitted_fraction_for(self, photon_energy):
        energies = photon_energy.copy()
        energies.sort()
        idx = np.searchsorted(self.photon_energy, energies)
        mu_over_rho = self.mu_over_rho[idx]
        if np.isnan(mu_over_rho).any():
            raise Exception('no suitable energy found for {}'.format(name))
        df = pd.DataFrame({
            'material': self.name,
            'mu_over_rho': mu_over_rho,
            'density': self.density,
            'material_thickness': self.thickness,
        })
        df = self._transmitted_fraction(df)
        return df

    def _transmitted_fraction(self, df):
        """Calculate the fraction of photons that escape through various materials.

        mu_over_row: mass attenuation coefficient, in cm^2/g.  Values taken from tables
        such as this one: http://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z82.html.
        For more information, see: http://physics.nist.gov/PhysRefData/XrayMassCoef/chap2.html.
        """
        def thickness(value):
            match = re.match(r'(\d+)\s*cm', value)
            if match:
                return float(match.group(1))
            raise Exception('do not recognize thickness: {}'.format(value))
        thickness_cm = df.material_thickness.apply(thickness)
        df['transmitted_fraction'] = np.exp(-df.mu_over_rho * df.density * thickness_cm)
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

    def material_for(self, name, thickness):
        df = self.df_materials[self.df_materials.material == name]
        if len(df) < 1:
            raise Exception('no material found: {}'.format(name))
        return Material(name, thickness, df)
