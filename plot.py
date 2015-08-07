#!/usr/bin/env python

import sys
import util

filename = sys.argv[1]

if filename.endswith(".sn"):
    util.formats.Kinetic(filename).plot()
elif filename.endswith(".pn"):
    util.formats.Moment(filename).plot()
else:
    print("Unsupported file")
