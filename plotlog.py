#!/usr/bin/env python

import sys
import util

filename = sys.argv[1]

if filename.endswith(".sn"):
    util.formats.Kinetic(filename).logplot()
elif filename.endswith(".pn"):
    util.formats.Moment(filename).logplot()
else:
    print("Unsupported file")
