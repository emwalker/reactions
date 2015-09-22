import argparse


from lenrmc.nubase import System


class App(object):

    def __init__(self, **kwargs):
        self.kwargs = kwargs

    def run(self):
        s = System.parse(self.kwargs['system'], **self.kwargs)
        for reaction in s.terminal():
            print(reaction)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('system', type=str)
    parser.add_argument('--lb', dest='lower_bound', help='lower bound in keV')
    parser.set_defaults(lower_bound=-1000)
    return parser.parse_args()


if '__main__' == __name__:
    opts = parse_arguments()
    App(**vars(opts)).run()
