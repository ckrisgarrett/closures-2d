#!/usr/bin/env python

import os
import math
import shutil
import subprocess

from util import *
import formats


DEVNULL = open(os.devnull, "w+b")

class Test(object):
    def __init__(self):
        self.setup_files = dict()

    def run(self, commands):
        with ResetFile("out_0.000_0.pn", "out_0.000_0.sn", *self.setup_files.keys()):
            for dst,txt in self.setup_files.items():
                with open(dst, "w") as f:
                    f.write(txt)
            with Timed():
                for e in commands:
                    subprocess.call(e, stdout=DEVNULL)

    def kinetic_deck(self):
        return "QUAD_ORDER 12\nCFL_FACTOR 0.9"

    def moment_deck(self):
        filterno = {"none": 0, "hauck": 1, "sspline": 2, "lanczos": 3}
        return (
"""MOMENT_ORDER 3
QUAD_ORDER 30
FILTER_TYPE %d
FILTER_TUNE 100
CFL_FACTOR 0.9
""" % filterno[self.filter_type])

class RegressionTest(Test):
    def __init__(self):
        super(RegressionTest, self).__init__()
        self._relative_error = None

    @property
    def reference_path(self):
        return os.path.join("tests", "regression", self.solver, "%s.%s" % (self.initcond, self.ext))

    @property
    def current_path(self):
        return "out_1.030_0.%s" % self.ext

    @property
    def relative_error(self):
        if self._relative_error:
            return self._relative_error
        else:
            return abs(self.current - self.reference) / abs(self.reference)

    def show_result(self):
        error = self.relative_error
        if error == 0.0:
            resultprint("OK")
        else:
            resultprint("%0e" % self.relative_error)

    def run(self, commands):
        self.setup_files["input.deck"] = self.input_deck()
        with ResetFile(self.current_path):
            super(RegressionTest, self).run(commands)
            self.reference = self.parser(self.reference_path)
            self.current = self.parser(self.current_path)
        self.show_result()
        return self.current

    def input_deck(self):
        if self.initcond == "delta":
            sigma = 0.0
        else:
            sigma = 1.7
        condno = {"delta": 0, "gaussian": 0, "lattice": 1, "smooth": 2}
        return (
"""SOLVER %s
NUM_CELLS_X 36
NUM_CELLS_Y 71
NUM_MPI_PARTITIONS_X 1
NUM_MPI_PARTITIONS_Y 1
A_X -1.55
B_X 1.45
A_Y -1.62
B_Y 1.38
T_FINAL 1.03
OUT_DELTA_T 1000.1
GAUSSIAN_SIGMA %f
FLOOR 1e-4
INIT_COND %d
SIGMA 1.111""" % (self.solver, sigma, condno[self.initcond]))

class KineticRegressionTest(RegressionTest):
    def __init__(self, initcond):
        super(KineticRegressionTest, self).__init__()
        self.parser = formats.Kinetic
        self.ext = "sn"
        self.solver = "kinetic"
        self.initcond = initcond

    def run(self,commands):
        bulletprint("%s,%s,%s" % (os.path.split(commands[0])[1], self.solver, self.initcond))
        self.setup_files["kinetic.deck"] = self.kinetic_deck()
        return super(KineticRegressionTest, self).run(commands)

class MomentRegressionTest(RegressionTest):
    def __init__(self, filter_type, initcond):
        super(MomentRegressionTest, self).__init__()
        self.parser = formats.Moment
        self.ext = "pn"
        self.solver = "moment"
        self.filter_type = filter_type
        self.initcond = initcond

    def run(self, commands):
        bulletprint("%s,%s,%s(%s)" % (os.path.split(commands[0])[1],
            self.solver, self.initcond, self.filter_type))
        self.setup_files["moment.deck"] = self.moment_deck()
        self.message = "(%s)" % self.filter_type
        return super(MomentRegressionTest, self).run(commands)

    @property
    def reference_path(self):
        return os.path.join("tests", "regression", "%s-%s" % (self.solver, self.filter_type),
            "%s.%s" % (self.initcond, self.ext))

class ConvergenceTest(Test):
    def __init__(self):
        super(ConvergenceTest, self).__init__()

    def run(self, iterations):
        samples = []
        errors = []
        for i in reversed(range(iterations)):
            self.setup_files["input.deck"] = self.input_deck(i)
            bulletprint("dx / %d" % 2**i, indent=2)
            with ResetFile(self.current_path):
                super(ConvergenceTest, self).run([os.path.realpath("solver_omp.x")])
                samples.append(self.parser(self.current_path))
            if len(samples) > 1:
                errors.append(samples[-1] - samples[0])
            if len(errors) > 1:
                resultprint("%f" % math.log(abs(errors[-1]) / abs(errors[-2]), 2))
            else:
                resultprint("")
        return samples

    @property
    def current_path(self):
        return "out_0.330_0.%s" % self.ext

    def input_deck(self, scale):
        return ("""SOLVER %s
NUM_CELLS_X %d
NUM_CELLS_Y %d
NUM_MPI_PARTITIONS_X 1
NUM_MPI_PARTITIONS_Y 1
A_X -1.5
B_X 1.5
A_Y -1.5
B_Y 1.5
T_FINAL 0.33
OUT_DELTA_T 1000.1
GAUSSIAN_SIGMA 0.4
FLOOR 1e-4
INIT_COND 2
SIGMA 1.0""" % (self.solver, 32 * 2**(scale), 32 * 2**(scale)))

class KineticConvergenceTest(ConvergenceTest):
    def __init__(self):
        super(KineticConvergenceTest, self).__init__()
        self.parser = formats.Kinetic
        self.ext = "sn"
        self.solver = "kinetic"

    def run(self, iterations):
        bulletprint("kinetic")
        resultprint("")
        self.setup_files["kinetic.deck"] = self.kinetic_deck()
        return super(KineticConvergenceTest, self).run(iterations)

class MomentConvergenceTest(ConvergenceTest):
    def __init__(self, filter_type):
        super(MomentConvergenceTest, self).__init__()
        self.parser = formats.Moment
        self.ext = "pn"
        self.solver = "moment"
        self.filter_type = filter_type

    def run(self, iterations):
        bulletprint("moment(%s)" % self.filter_type)
        resultprint("")
        self.setup_files["moment.deck"] = self.moment_deck()
        return super(MomentConvergenceTest, self).run(iterations)

def _convergence(Parser, ext, solver, n):
    bulletprint(solver)
    resultprint("")
    samples = []
    errors = []

    for i in reversed(range(n)):
        bulletprint("dx / %d" % 2**i, indent=2)
        with ResetFile("input.deck", "%s.deck" % solver, "out_0.000_0.%s" % ext, "out_0.330_0.%s" % ext):
            shutil.copyfile(os.path.join("tests", "convergence", "%s%d.deck" % (solver, 2**i)), "input.deck")
            shutil.copyfile(os.path.join("tests", "%s.deck" % solver), "%s.deck" % solver)
            with Timed():
                subprocess.call([os.path.realpath("solver_omp.x")], stdout=DEVNULL)
            s = Parser("out_0.330_0.%s" % ext)
            samples.append(s)
        if len(samples) > 1:
            errors.append(samples[-1] - samples[0])
        if len(errors) > 1:
            resultprint("%f" % math.log(abs(errors[-1]) / abs(errors[-2]), 2))
        else:
            resultprint("")
    return samples


def convergence(solver, n):
    if solver == "kinetic":
        return _convergence(formats.Kinetic, "sn", solver, n)
    else:
        return _convergence(formats.Moment, "pn", solver, n)

def _regression(Parser, ext, program, initcond, solver, variant=""):
    bulletprint("%s, %s, %s" % (os.path.split(program[0])[1], solver, initcond))
    with ResetFile("input.deck", "%s.deck" % solver, "out_0.000_0.%s" % ext, "out_1.030_0.%s" % ext):
        shutil.copyfile(os.path.join("tests", "regression", solver, "%s.deck" % initcond), "input.deck")
        shutil.copyfile(os.path.join("tests", "%s.deck" % solver), "%s.deck" % solver)
        with Timed():
            subprocess.call(program, stdout=DEVNULL)
        current = Parser("out_1.030_0.%s" % ext)
        reference = Parser(os.path.join("tests", "regression", solver, "%s.%s" % (initcond, ext)))
        relativeError = abs(current - reference) / abs(reference)
        if relativeError == 0.0:
            resultprint("OK")
        else:
            resultprint("%0e" % relativeError)
        return current

def regression(program, initcond, solver, variant=""):
    if solver == "kinetic":
        return _regression(formats.Kinetic, "sn", program, initcond, solver)
    else:
        return _regression(formats.Moment, "pn", program, initcond, solver)
