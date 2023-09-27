#!/usr/bin/env python3

import sys
from argparse import ArgumentParser
import numpy as np


# Parse args
parser = ArgumentParser(description='n choose 2')
parser.add_argument("-n", "--size", type=int, help="Number of elements")
parser.add_argument("-o", "--out", type=str, help="Output file")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

# Define functions
def numpy_combinations(x):
    idx = np.stack(np.triu_indices(len(x), k=1), axis=-1)
    return x[idx]

# Run
c = numpy_combinations(np.arange(args.size))

out_fh = open(args.out, 'w')
np.savetxt(out_fh, c, fmt = "%d")
out_fh.close()
