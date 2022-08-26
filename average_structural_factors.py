
from mdtools.utils import *
import argparse

# parameters
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="Input file for the crystal system")
# parser.add_argument("-o", "--output", type=str, help="Prefix for the output trajectory and state files")
parser.add_argument("-n", "--n_chains", type=int, help="Number of protein chains in the system", default=4)
parser.add_argument("-N", "--n_epochs", type=int, help="Number of epochs to average over")
args = parser.parse_args()

for phase in ['pos', 'neg', 'zero']:
    for k in range(args.n_chains):
        fname = args.input+f'_epoch_{i}_chainwise_{phase}_subtraj_{k}_avg'
        average_structure_factors(fname, max_frame=args.n_epochs)
