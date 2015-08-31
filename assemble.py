#!/usr/bin/env python

import argparse

import util

parser = argparse.ArgumentParser(description='Assemble MPI data files. The output will say "mpi" instead of a node number')
parser.add_argument('filename_pattern', help='Node number is given by {}, e.g. out_1.000_{}.pn')
parser.add_argument('nodes_x', type=int)
parser.add_argument('nodes_y', type=int)
args = parser.parse_args()

output_filename = args.filename_pattern.format('mpi')

if args.filename_pattern.endswith(".sn"):
    util.formats.Composite(util.formats.Kinetic, args.filename_pattern, args.nodes_x, args.nodes_y).output(output_filename)
elif args.filename_pattern.endswith(".pn"):
    util.formats.Composite(util.formats.Moment, args.filename_pattern, args.nodes_x, args.nodes_y).output(output_filename)
else:
    print("Unsupported file")
