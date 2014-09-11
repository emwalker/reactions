import math

import pandas as pd
import numpy as np

from .helpers import centimeters_to_float, MultiDataFrame


class CompositeLayer(object):
    def __init__(self, layers, **kwargs):
        self.name = kwargs.get('name', '/'.join(map(lambda l: l.name, layers)))
        self.layers = layers

    def transmitted_fraction_for(self, photon_energies):
        if self.layers:
            dfs = []
            fractions = []
            thicknesses = []
            result = None
            for material in self.layers:
                df = self._transmitted_fraction(material, photon_energies, fractions, thicknesses)
                dfs.append(df)
            multi_df = MultiDataFrame(dfs)
        else:
            multi_df = None
        return multi_df

    def _transmitted_fraction(self, material, photon_energies, fractions, thicknesses):
        thicknesses.append(material.thickness)
        result = material.transmitted_fraction_for(photon_energies).last
        fractions.append(result.transmitted_fraction)
        df = result.copy()
        df['transmitted_fraction'] = reduce(lambda l,r: l * r, fractions)
        df['material'] = self.name
        df['material_thickness'] = '/'.join(thicknesses)
        return df


class Layer(object):
    @classmethod
    def of(cls, layers, **kwargs):
        return CompositeLayer(layers, **kwargs)

    def __init__(self, name, thickness, df):
        self.df = df.copy()
        self.density = df.density.iloc[0]
        self.df['name'] = self.name = name
        self.df['material_thickness'] = self.thickness = thickness
        self.df = df.sort('photon_energy')
        self.photon_energy = df.photon_energy.values
        self.mu_over_rho = df.mu_over_rho.values

    def transmitted_fraction_for(self, photon_energies):
        energies = photon_energies.copy()
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
        return MultiDataFrame([df])

    def _transmitted_fraction(self, df):
        """Calculate the fraction of photons that escape through various materials.

        mu_over_row: mass attenuation coefficient, in cm^2/g.  Values taken from tables
        such as this one: http://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z82.html.
        For more information, see: http://physics.nist.gov/PhysRefData/XrayMassCoef/chap2.html.
        """
        thickness_cm = df.material_thickness.apply(centimeters_to_float)
        df['transmitted_fraction'] = np.exp(-df.mu_over_rho * df.density * thickness_cm)
        return df


class DetectorLayer(Layer):
    def __init__(self, name, **kwargs):
        self.name = name
        self.diameter = centimeters_to_float(kwargs['diameter'])
        self.distance = centimeters_to_float(kwargs['distance'])
        self.efficiency = float(kwargs['efficiency'])
        self.thickness = 'D'
        ratio = (1/2.) * self.diameter / self.distance
        self.solid_angle = 2 * math.pi * (1 - math.cos(math.atan(ratio)))
        self.relative_coverage = self.solid_angle / (4 * math.pi)
        self.transmitted_fraction = self.relative_coverage * self.efficiency

    def transmitted_fraction_for(self, photon_energies):
        df = pd.DataFrame({
            'material': [self.name],
            'material_thickness': [self.thickness],
            'transmitted_fraction': [self.transmitted_fraction],
        })
        return MultiDataFrame([df])
