import re


def centimeters_to_float(string):
    match = re.match(r'(\d+)\s*cm', string)
    if match:
        return float(match.group(1))
    raise Exception('do not recognize value: {}'.format(string))
