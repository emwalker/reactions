import os
import math

import pandas as pd
import numpy as np

from .helpers import centimeters_to_float


class Table(object):
    def __init__(self, name, df, **kwargs):
        self.name = name
        self.df = df
        self.is_canonical = bool(kwargs.get('canonical'))

    def __mul__(self, o):
        return o.apply_table(self)

    def merge(self, o, **kwargs):
        merge_keys, merge_functions = self._partition_keys(kwargs.get('on'))
        df_left, df_right = self.df, o.df
        for prepare_function in merge_functions:
            df_left, df_right = prepare_function(df_left, df_right, self.is_canonical, o.is_canonical)
        if len(merge_keys):
            df = df_left.merge(df_right, on=merge_keys)
        else:
            try:
                df = df_left.merge(df_right)
            except pd.tools.merge.MergeError:
                df = self._cross_product(df_left, df_right)
        df = df.drop_duplicates(take_last=True)
        table = Table('{} x {}'.format(self.name, o.name), df)
        return table

    def apply_table(self, table):
        return table.merge(self)

    def apply_array(self, array):
        objects = []
        for previous in reversed(array):
            obj = previous * self
            objects.append(obj)
        return Array("{}'".format(array.name), objects)

    @property
    def objects(self):
        return [self]

    def __len__(self):
        return len(self.df)

    def __eq__(self, o):
        df1, df2 = self.df, o.df
        return np.all(self.df == o.df)

    def apply_operator(self, operator):
        name = '{} x {}'.format(operator.name, self.name)
        result = Table(name, self.df.copy())
        # Operators are applied from right to left
        for function in reversed(operator.transformations):
            result = function(result)
        return result

    @property
    def columns(self):
        return self.df.columns

    @property
    def values(self):
        return self.df.values

    def unique_value_at(self, column):
        values = self.df[column].unique()
        if len(values) > 1:
            raise ValueError('more than one value for {}'.format(column))
        return values[0]

    def unique_values_at(self, columns):
        return (self.unique_value_at(column) for column in columns)

    def _cross_product(self, df_left, df_right):
        df_left, df_right = df_left.copy(), df_right.copy()
        for df in (df_left, df_right):
            if '__join_key' not in df:
                continue
            raise MergeError('the column "__join_key" already exists')
        df_left['__join_key'] = 1
        df_right['__join_key'] = 1
        product = df_left.merge(df_right, on='__join_key')
        del product['__join_key']
        return product

    def _partition_keys(self, keys):
        if not keys:
            return [], []
        merge_keys, merge_functions = [], []
        for item in keys:
            if MergeFunction == type(item):
                merge_keys.append(item.column)
                merge_functions.append(item.prepare)
            else:
                merge_keys.append(item)
        return merge_keys, merge_functions


def read_csv(name, relpath, **kwargs):
    basepath = os.path.dirname(__file__)
    abspath = os.path.join(basepath, relpath)
    return Table(name, pd.read_csv(abspath, header=0), **kwargs)


class Array(object):
    def __init__(self, name, objects):
        self.name = name
        self.objects = objects

    def merge(self, o, **kwargs):
        objects = []
        for left in self.objects:
            for right in o.objects:
                table = left.merge(right, **kwargs)
                objects.append(table)
        return Array("{}'".format(self.name), objects)

    def apply_operator(self, operator):
        name = '{} x {}'.format(self.name, operator.name)
        results = map(lambda t: t.apply_operator(operator), self.objects)
        return Array(name, results)

    def __mul__(self, o):
        return o.apply_array(self)

    def __getitem__(self, idx):
        return self.objects[idx]

    def __len__(self):
        return len(self.objects)

    def __iter__(self):
        return iter(self.objects)

    def tolist(self):
        return map(lambda o: o.values.tolist(), self.objects)

    def __eq__(self, o):
        for i, left in enumerate(self):
            right = o[i]
            if ~(left == right):
                return False
        return True

    @property
    def combined(self):
        df = pd.concat(map(lambda o: o.df, self.objects))
        return Table(self.name, df)

    def partition_by(self, columns):
        grouped = self.combined.df.groupby(columns)
        tables = []
        for compound_key, group in grouped:
            table = Table(', '.join(compound_key), group)
            tables.append(table)
        return tables


class MergeFunction(object):
    def __init__(self, column, prepare):
        self.column = column
        self.prepare = prepare


def table(name, **kwargs):
    df = pd.DataFrame(kwargs['values'], columns=kwargs['columns'])
    return Table(name, df)


def closest(column):
    def prepare(left, right, left_canonical, right_canonical):
        if bool(left_canonical) != bool(right_canonical):
            left, right = [left, right] if left_canonical else [right, left]
        else:
            left, right = [left, right] if len(left) >= len(right) else [right, left]
        left = left.sort(column)
        right = right.sort(column)
        left_indexes = np.searchsorted(left[column].values, right[column])
        right[column] = left[column].values[left_indexes]
        return left, right
    return MergeFunction(column, prepare)


class Operator(object):
    def __init__(self, name, transformations):
        self.name = name
        self.transformations = transformations

    def __mul__(self, o):
        return o.apply_operator(self)

    def apply_array(self, array):
        results = map(lambda o: o.apply_operator(self), array.objects)
        return Array("{}'".format(array.name), results)

    def apply_operator(self, operator):
        name = '{} x {}'.format(self.name, operator.name)
        return Operator(name, operator.transformations + self.transformations)


##
#  Helpers
#

def layer(name, *args):
    if len(args) == 1:
        layer = reduce(lambda l,r: l * r, reversed(args[0]))
    elif len(args) == 2:
        df = pd.DataFrame([args], columns=['material', 'material_thickness'])
        layer = Table('{} (layer)'.format(name), df)
    else:
        raise ValueError('do not know how to handle arguments: {}'.format(args))
    return layer


##
#  Operators
#

def escaping_photons():
    def photons(table):
        df = table.df
        df['escaping_photons'] = df.transmitted_fraction * df.photons_per_second
        return table
    return Operator('Escaping photons', [photons])


def transmitted_fraction():
    def transmitted(table):
        df = table.df
        thickness_cm = df.material_thickness.apply(centimeters_to_float)
        df['transmitted_fraction'] = np.exp(-df.mu_over_rho * df.density * thickness_cm)
        return table
    return Operator('Transmitted fraction', [transmitted])


def detector(name, **kwargs):
    diameter = centimeters_to_float(kwargs['diameter'])
    distance = centimeters_to_float(kwargs['distance'])
    efficiency = kwargs['efficiency']
    ratio = (1/2.) * diameter / float(distance)
    solid_angle = 2 * math.pi * (1 - math.cos(math.atan(ratio)))
    def transmitted(table):
        df = table.df
        relative_coverage = solid_angle / (4 * math.pi)
        df['transmitted_fraction'] = df.transmitted_fraction * relative_coverage * efficiency
        return table
    return Operator(name, [transmitted])


def photon_counts(events):
    def counts(table):
        df = table.df
        df['photons_per_event'] = np.where(
            df.photons_per_transition.isnull(),
            df.q_max_mev / df.photon_energy,
            df.photons_per_transition
        )
        df['events_per_transition_per_second'] = np.where(
            df.half_life_seconds.isnull(),
            1,
            0.693/df.half_life_seconds
        )
        df['combined_cross_section'] = df.cross_section_barns.sum()
        df['fractional_cross_section'] = df.cross_section_barns / df.combined_cross_section
        rate = df.fractional_cross_section * \
            df.isotopic_abundance * df.events_per_transition_per_second
        df['events_per_second'] = events
        df['photons_per_second'] = rate * df.photons_per_event * df.events_per_second
        return table
    return Operator('Photon counts', [counts])

##
#  Useful data sets
#

materials   = read_csv('Materials',   '../db/materials.csv',   canonical=True)
channels    = read_csv('Channels',    '../db/channels.csv',    canonical=True)
transitions = read_csv('Transitions', '../db/transitions.csv', canonical=True)
