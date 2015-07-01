import os
import sys
import math
import shutil
import subprocess

import util

TOLERANCE = 1e-10

REGENERATE = False
# To regenerate regression test answers (due to changes in the code), delete
# everything under tests/regression, flip this flag, and run regression tests
# in serial mode. The script will write out new answers.
# Remember to flip the flag back!

INITCOND = ["delta", "gaussian", "lattice", "smooth"]
SOLVER = ["kinetic", "moment"]
MOMENT_FILTER = ["none", "sspline", "hauck", "lanczos"]

DEVNULL = open(os.devnull, "w+b")


def regression_test(target, source, env):
    try:
        util.tee(str(target[0]), "regression tests:")
        for p in [str(x) for x in source]:
            if p.startswith("solver_mpi"):
                commands = [["mpirun", "-np", "4", p]]
                composite = True
            else:
                commands = [os.path.realpath(p)]
                composite = False
            for s in SOLVER:
                for i in INITCOND:
                    if s == "kinetic":
                        util.bullet(str(target[0]), "%s,%s,%s" % (p, s, i), wait=False)
                        if composite:
                            KineticCompositeRegressionTest(str(target[0]), i).run(commands)
                        else:
                            KineticRegressionTest(str(target[0]), i).run(commands)
                    elif s == "moment":
                        for f in MOMENT_FILTER:
                            util.bullet(str(target[0]), "%s,%s,%s(%s)" % (p, s, i, f), wait=False)
                            if composite:
                                MomentCompositeRegressionTest(str(target[0]), f, i).run(commands)
                            else:
                                MomentRegressionTest(str(target[0]), f, i).run(commands)
    except:
        try:
            os.unlink(str(target[0]))
        except:
            pass
        raise

def convergence_test(target, source, env):
    try:
        util.tee(str(target[0]), "convergence tests:")
        for s in SOLVER:
            if s == "kinetic":
                util.bullet(str(target[0]), "%s" % s, wait=False)
                KineticConvergenceTest(str(target[0])).run([os.path.realpath(str(source[0]))], 5)
            elif s == "moment":
                for f in MOMENT_FILTER:
                    util.bullet(str(target[0]), "%s(%s)" % (s, f), wait=False)
                    MomentConvergenceTest(str(target[0]), f).run([os.path.realpath(str(source[0]))], 5)
    except:
        try:
            os.unlink(str(target[0]))
        except:
            pass
        raise

def test_everything(target, source, env):
    with open(str(target[0]), "w") as f:
        f.write("tests complete\n")


class Test(object):
    def __init__(self):
        self.setup_files = dict()

    @property
    def initial_path(self):
        return "out_0.000_0.%s" % self.ext

    def run(self, commands):
        with util.ResetFile(*self.setup_files.keys()):
            for dst,txt in self.setup_files.items():
                with open(dst, "w") as f:
                    f.write(txt)
            with util.Timed(self.log_path):
                for e in commands:
                    subprocess.call(e, stdout=DEVNULL)


class RegressionTest(Test):
    def __init__(self):
        super(RegressionTest, self).__init__()
        self._relative_error = None
        self._mass_change = None

    @property
    def current_path(self):
        return "out_1.030_0.%s" % self.ext

    @property
    def gaussian_sigma(self):
        if self.initcond == "delta":
            return 0.0
        else:
            return 1.7

    @property
    def relative_error(self):
        if not self._relative_error:
            self._relative_error = abs(self.current -
                self.reference) / abs(self.reference)
        return self._relative_error


    @property
    def mass_change(self):
        if not self._mass_change:
            self._mass_change = self.current.mass - self.initial.mass
        return self._mass_change


    def run(self, commands):
        self.setup_files["input.deck"] = util.decks.input_deck(solver=self.solver,
            num_cells_x=36, num_cells_y=74, a_x=-1.55, b_x=1.45, a_y=-1.62, b_y=1.38,
            t_final=1.03, sigma=1.111, init_cond = self.initcond,
            gaussian_sigma=self.gaussian_sigma)
        with util.ResetFile(self.current_path, self.initial_path):
            util.right(self.log_path, "execution time:")
            util.tee(self.log_path, " " * 10, wait=True)
            super(RegressionTest, self).run(commands)
            util.tee(self.log_path, "")
            self.initial = self.parser(self.initial_path)
            self.current = self.parser(self.current_path)
            if REGENERATE:
                try:
                    os.makedirs(os.path.dirname(self.reference_path))
                except OSError:
                    pass
                shutil.copyfile(self.current_path, self.reference_path)
            self.reference = self.parser(self.reference_path)

        self.show_results()

        return self.current

    def show_results(self):
        util.right(self.log_path, "relative error:")
        if abs(self.relative_error) < TOLERANCE:
            util.result(self.log_path, "OK")
        else:
            util.result(self.log_path, "%0e" % self.relative_error)

        util.right(self.log_path, "mass change:")
        if self.initcond == "delta" or self.initcond == "smooth":
            if abs(self.mass_change) < TOLERANCE:
                util.result(self.log_path, "OK")
            else:
                util.result(self.log_path, "%0e" % self.mass_change)
        else:
            util.result(self.log_path, "--")

class CompositeRegressionTest(RegressionTest):
    def __init__(self, xNodes, yNodes):
        RegressionTest.__init__(self)
        self.xNodes = xNodes
        self.yNodes = yNodes

    @property
    def initial_path(self):
        return "out_0.000_{}.%s" % self.ext

    @property
    def current_path(self):
        return "out_1.030_{}.%s" % self.ext

    @property
    def initial_paths(self):
        return [self.initial_path.format(x) for x in range(self.xNodes * self.yNodes)]

    @property
    def current_paths(self):
        return [self.current_path.format(x) for x in range(self.xNodes * self.yNodes)]

    def run(self, commands):
        self.setup_files["input.deck"] = util.decks.input_deck(solver=self.solver,
            num_cells_x=36, num_cells_y=74, a_x=-1.55, b_x=1.45, a_y=-1.62, b_y=1.38,
            t_final=1.03, sigma=1.111, init_cond = self.initcond,
            gaussian_sigma=self.gaussian_sigma, num_mpi_partitions_x=self.xNodes,
            num_mpi_partitions_y=self.yNodes)
        with util.ResetFile(*(self.initial_paths + self.current_paths)):
            util.right(self.log_path, "execution time:")
            util.tee(self.log_path, " " * 10, wait=True)
            super(RegressionTest, self).run(commands)
            util.tee(self.log_path, "")
            self.initial = util.formats.Composite(self.parser,
                self.initial_path, self.xNodes, self.yNodes)
            self.current = util.formats.Composite(self.parser,
                self.current_path, self.xNodes, self.yNodes)
            self.reference = self.parser(self.reference_path)

        self.show_results()

        return self.current

class KineticRegressionTest(RegressionTest):
    def __init__(self, log_path, initcond):
        super(KineticRegressionTest, self).__init__()
        self.parser = util.formats.Kinetic
        self.ext = "sn"
        self.solver = "kinetic"
        self.initcond = initcond
        self.log_path = log_path

    @property
    def reference_path(self):
        return os.path.join("tests", "regression", self.solver, "%s.sn" % self.initcond)

    def run(self, commands):
        self.setup_files["kinetic.deck"] = util.decks.kinetic_deck()
        return super(KineticRegressionTest, self).run(commands)

class KineticCompositeRegressionTest(CompositeRegressionTest, KineticRegressionTest):
    def __init__(self, log_path, initcond):
        CompositeRegressionTest.__init__(self, 2, 2)
        KineticRegressionTest.__init__(self, log_path, initcond)

    def run(self, commands):
        self.setup_files["kinetic.deck"] = util.decks.kinetic_deck()
        return super(KineticCompositeRegressionTest, self).run(commands)

class MomentRegressionTest(RegressionTest):
    def __init__(self,log_path,  filter_type, initcond):
        super(MomentRegressionTest, self).__init__()
        self.parser = util.formats.Moment
        self.ext = "pn"
        self.solver = "moment"
        self.filter_type = filter_type
        self.initcond = initcond
        self.log_path = log_path

    @property
    def reference_path(self):
        return os.path.join("tests", "regression", self.solver, self.initcond, "%s.pn" % self.filter_type)

    def run(self, commands):
        self.setup_files["moment.deck"] = util.decks.moment_deck()
        return super(MomentRegressionTest, self).run(commands)

class MomentCompositeRegressionTest(CompositeRegressionTest, MomentRegressionTest):
    def __init__(self, log_path, filter_type, initcond):
        CompositeRegressionTest.__init__(self, 2, 2)
        MomentRegressionTest.__init__(self, log_path, filter_type, initcond)

    def run(self, commands):
        self.setup_files["moment.deck"] = util.decks.moment_deck()
        return super(MomentCompositeRegressionTest, self).run(commands)

class ConvergenceTest(Test):
    def __init__(self):
        super(ConvergenceTest, self).__init__()

    @property
    def current_path(self):
        return "out_0.330_0.%s" % self.ext

    def run(self, commands, iterations):
        samples = []
        errors = []
        for i in reversed(range(iterations)):
            self.setup_files["input.deck"] = util.decks.input_deck(solver=self.solver,
                num_cells_x=32*2**i, num_cells_y=32*2**i, t_final=0.33, init_cond="smooth")
            util.bullet(self.log_path, "dx / %d" % 2**i, indent=2)
            with util.ResetFile(self.current_path, self.initial_path):
                super(ConvergenceTest, self).run(commands)
                samples.append(self.parser(self.current_path))
            if len(samples) > 1:
                errors.append(samples[-1] - samples[0])
            if len(errors) > 1:
                util.result(self.log_path,
                    "%f" % math.log(abs(errors[-1]) / abs(errors[-2]), 2))
            else:
                util.result(self.log_path, "")
        return samples

class KineticConvergenceTest(ConvergenceTest):
    def __init__(self, log_path):
        super(KineticConvergenceTest, self).__init__()
        self.parser = util.formats.Kinetic
        self.ext = "sn"
        self.solver = "kinetic"
        self.log_path = log_path

    def run(self, commands, iterations):
        self.setup_files["kinetic.deck"] = util.decks.kinetic_deck()
        return super(KineticConvergenceTest, self).run(commands, iterations)

class MomentConvergenceTest(ConvergenceTest):
    def __init__(self, log_path, filter_type):
        super(MomentConvergenceTest, self).__init__()
        self.parser = util.formats.Moment
        self.ext = "pn"
        self.solver = "moment"
        self.filter_type = filter_type
        self.log_path = log_path

    def run(self, commands, iterations):
        self.setup_files["moment.deck"] = util.decks.moment_deck()
        return super(MomentConvergenceTest, self).run(commands, iterations)
