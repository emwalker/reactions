import pandas as pd


class Channel(object):
    def __init__(self, df):
        self.df = df

    @property
    def escaped_photons(self):
        return self.df


class Generation(object):
    def __init__(self, df):
        self.df = df
        self._calculate_escaping_photons()

    def channel_for(self, transition, name):
        df = self.df[(self.df.transition == transition) & (self.df.channel == name)]
        return Channel(df)

    def _calculate_escaping_photons(self):
        pass


class Engine(object):
    @classmethod
    def from_csv(cls, materials_path, initpath):
        df_materials = pd.read_csv(materials_path, header=0)
        df_generation = pd.read_csv(initpath, header=0)
        return cls(df_materials, df_generation)

    def __init__(self, df_materials, df_generation):
        self.df_materials = df_materials
        self.generations = [Generation(df_generation)]
