from mdtools.utils import compute_net_dipole_moment
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input file name (to be expanded with _i)")
parser.add_argument("-o", "--output", type=str, help="output file name")
parser.add_argument("-p", "--partial_charges_fname", type=str, default='partial_charges.npy')
parser.add_argument("-n", "--n_chains", type=int, help="number of chains to process")
args = parser.parse_args()

partial_charges = np.load(args.partial_charges_fname)
input_lst = [f'{args.input}_{i}' for i in range(args.n_chains)]
if args.output is None:
    args.output = args.input + '_net_dipole_moment'
output_lst = [f'{args.output}_{i}' for i in range(args.n_chains)]
out = compute_net_dipole_moment(partial_charges, input=input_lst, output=output_lst)
