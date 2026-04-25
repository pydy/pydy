#!/usr/bin/env python

import os
import argparse

from packaging.version import parse as parse_version
import IPython

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--notebook", action='store_true',
                    help="Launches the example in a Jupyter notebook.")
args = parser.parse_args()

if __name__ == "__main__":
    if args.notebook:
        if parse_version(IPython.__version__) < parse_version('3.0'):
            os.system('ipython notebook')
        else:
            os.system('jupyter notebook')
    else:
        os.system('python chaos_pendulum.py')
