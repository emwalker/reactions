import re

import pandas as pd


def centimeters_to_float(string):
    match = re.match(r'([\d\.]+)\s*cm', string)
    if match:
        return float(match.group(1))
    raise ValueError('do not recognize value: {}'.format(string))
