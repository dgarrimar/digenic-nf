#!/usr/bin/env python3

import sys
from argparse import ArgumentParser
import re
import itertools

# Parse args
parser = ArgumentParser(description='Find possible pairs (different chromosomes) given a list of variants')
parser.add_argument("-i", "--input", type=str, help="Input variants file")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

# Run
with open(args.input, 'r') as f:
    lines = f.read().splitlines()

chrs = [str(i+1) for i in range(22)]
chrs.append('X')

for chr in chrs:
    In = list(filter(lambda x: re.match(chr + ":", x), lines))
    Out = list(filter(lambda y: not re.match(chr + ":", y), lines))
    prod = [item for item in itertools.product(In,Out)]
    lines = [i for i in lines if i not in In]
    for pair in prod:
        print(pair[0]+' '+pair[1])
