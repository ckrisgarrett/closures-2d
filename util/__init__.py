import os
import sys
import time
import tempfile

from . import formats
from . import decks


def tee(filename, s, wait=False):
    if wait:
        sys.stdout.write(s)
        sys.stdout.flush()
        tail = ""
    else:
        print(s)
        tail = "\n"
    if filename:
        with open(filename, "a") as f:
            f.write(s + tail)

def result(filename, s, wait=False):
    tee(filename, "%20s" % s, wait)

def lresult(filename, s, wait=False):
    tee(filename, "%-20s" % s, wait)

def name(filename, s, wait=True):
    tee(filename, "%-50s" % s, wait)

def bullet(filename, s, indent=1, wait=True):
    name(filename, " %s %s" % ("-" * indent, s), wait)

def right(filename, s, wait=True):
    tee(filename, "%50s" % s, wait)

class Timed(object):
    def __init__(self, filename=None):
        self.filename = filename

    def __enter__(self):
        self.start = time.time()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end = time.time()
        tee(self.filename, "%-10s" % ("%.3fs" % (self.end - self.start)), wait=True)

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
