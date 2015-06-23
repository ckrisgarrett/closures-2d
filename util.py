#!/usr/bin/env python

import os
import sys
import time
import tempfile

def nameprint(s):
    sys.stdout.write("%-53s" % s)
    sys.stdout.flush()

def resultprint(s):
    sys.stdout.write("%15s\n" % s)
    sys.stdout.flush()

def bulletprint(s, indent=1):
    nameprint(" %s %s" % ("-" * indent, s))


class Timed(object):
    def __enter__(self):
        self.start = time.time()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end = time.time()
        sys.stdout.write("%9.4fs" % (self.end - self.start))
        sys.stdout.flush()

        return False

class ResetFile(object):
    def __init__(self, *filenames):
        self.filenames = set(filenames)
        self.tmpfiles = dict()

    def __enter__(self):
        for s in self.filenames:
            try:
                with open(s, "rb") as f:
                    self.tmpfiles[s] = tempfile.NamedTemporaryFile()
                    self.tmpfiles[s].write(f.read())
                    self.tmpfiles[s].flush()
            except IOError:
                pass

    def __exit__(self, exc_type, exc_val, exc_tb):
        for s in self.filenames:
            try:
                self.tmpfiles[s].seek(0)
                with open(s, "wb") as f:
                    f.write(self.tmpfiles[s].read())
                self.tmpfiles[s].close()
            except KeyError:
                try:
                    os.unlink(s)
                except OSError:
                    pass
        return False
