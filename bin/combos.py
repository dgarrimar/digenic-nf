#!/usr/bin/env python3

import sys
from argparse import ArgumentParser
import itertools

# Parse args
parser = ArgumentParser(description='Find possible pairs given a list of variants')
parser.add_argument("-i", "--input", type=str, help="Input variants file")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

# Run
with open(args.input, 'r') as f:
    lines = f.read().splitlines()

for pair in itertools.combinations(lines, 2):
    print(str(pair[0])+' '+str(pair[1]))
