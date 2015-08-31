#!/usr/bin/env python

import sys
import util
import argparse

parser = argparse.ArgumentParser(description='Show a logarithmic plot of the given output file')
parser.add_argument('filename')
args = parser.parse_args()

if args.filename.endswith(".sn"):
    util.formats.Kinetic(args.filename).logplot()
elif args.filename.endswith(".pn"):
    util.formats.Moment(args.filename).logplot()
else:
    print("Unsupported file")
    sys.exit(-1)
