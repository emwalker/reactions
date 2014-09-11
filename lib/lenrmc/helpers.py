import re

import pandas as pd


def centimeters_to_float(string):
    match = re.match(r'(\d+)\s*cm', string)
    if match:
        return float(match.group(1))
    raise Exception('do not recognize value: {}'.format(string))


class MultiDataFrame(object):
    def __init__(self, dataframes):
        self.dataframes = dataframes

    @property
    def last(self):
        return self.dataframes[-1]

    @property
    def df(self):
        values = map(lambda d: d.values, self.dataframes)
        values = map(lambda v: v[-1].tolist(), values)
        return pd.DataFrame(values, columns=self.last.columns)

    @property
    def values(self):
        return self.df.values
