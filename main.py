"""
Main file for running entropy estimators
"""

import argparse
from support import *

parser = argparse.ArgumentParser()
parser.add_argument("-pmin", type=float,
                    help="Minimum mass")
parser.add_argument("-L", type=int,
                    help="Polynoimal degree. Default c0*log(1/pmin)")
parser.add_argument("-M", type=int,
                    help="M/n is the right-end of approximation interval. Default c1*log(1/pmin)")
parser.add_argument("-fin", type=str,
                    help="fingerprint data file")
parser.add_argument("-hist", type=str,
                    help="histogram data file")
args = parser.parse_args()


if args.pmin is None:
    raise Exception('Please input pmin!')


support = Support(pmin=args.pmin)
if args.L is not None:
    support.degree = args.L
if args.M is not None:
    support.ratio = args.M


if args.fin is not None:
    fin = np.loadtxt(args.fin, dtype=int)
elif args.hist is not None:
    hist = np.loadtxt(args.hist, dtype=int)
    fin = hist_to_fin(hist)
else:
    raise Exception('Please input fingerprint or histogram!')


print("")
print("Parameters:")
print("Preset pmin value\t=%d" % support.pmin)
print("Polynoimal degree\t=%d" % support.degree)
print("Approximation interval\t=[%.2e,%.2f/n]" % (support.pmin, support.ratio))
print("")
print("Results:")
print("Sample size\t=%d" % get_sample_size(fin))
print("Polynomial\t=%d" % int(support.estimate(fin)))
print("Plug-in\t\t=%d" % int(support.estimate_plug(fin)))
print("")
