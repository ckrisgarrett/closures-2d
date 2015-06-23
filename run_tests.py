#!/usr/bin/env python

import os
import sys
import argparse
import subprocess

import tests
from util import *

INITCOND = ["delta", "gaussian", "lattice", "smooth"]
SOLVER = ["kinetic", "moment"]
MOMENT_FILTER = ["none", "sspline", "hauck", "lanczos"]
OPTION = ["mpi", "omp", "cuda", "serial"]

DEVNULL = open(os.devnull, "w+b")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test convergence rates and compare against " +
        "test vectors. Called without arguments, this script just checks that everything compiles.")
    parser.add_argument("-a", "--all", action="store_true", help="run all tests")
    parser.add_argument("-c", "--convergence", action="store_true", help="test convergence rate")
    parser.add_argument("-r", "--regression", action="store_true", help="check against test vectors")
    parser.add_argument("-i", "--initcond", action="append", choices=INITCOND + ["all"], default=[])
    parser.add_argument("-s", "--solver", action="append", choices=SOLVER + ["all"], default=[])
    parser.add_argument("-o", "--option", action="append", choices=OPTION + ["all"], default=[])
    parser.add_argument("--moment-filter", action="append", choices=MOMENT_FILTER + ["all"], default=["all"])
    parser.add_argument("--convergence-iterations", type=int, default=5)

    args = parser.parse_args()

    if args.all or "all" in args.initcond:
        args.initcond = INITCOND
    if args.all or "all" in args.solver:
        args.solver = SOLVER
    if args.all or "all" in args.option:
        args.option = OPTION
    if args.all or "all" in args.moment_filter:
        args.moment_filter = MOMENT_FILTER

    if args.all or ("mpi" in args.option and "omp" in args.option):
        args.option.append("mpiomp")
    if args.all or ("mpi" in args.option and "cuda" in args.option):
        args.option.append("mpicuda")
    if args.all or ("omp" in args.option and "cuda" in args.option):
        args.option.append("ompcuda")
    if args.all or ("mpi" in args.option and "omp" in args.option and "cuda" in args.option):
        args.option.append("mpiompcuda")

    programs = [[os.path.realpath("solver_%s.x" % s)] for s in args.option]

    #TODO add conditional compilation
    with ResetFile("solver_serial.x", "solver_omp.x", "solver_mpi.x", "solver_cuda.x",
            "solver_mpiomp.x", "solver_mpicuda.x", "solver_ompcuda.x",
            "solver_mpiompcuda.x", "fobj_cuda.o", "solve_flux_cuda.o"):
        nameprint("compiling...")
        with Timed():
            rc = subprocess.call(["make", "all"], stderr=subprocess.STDOUT, stdout=DEVNULL)
        if rc is 0:
            resultprint("OK")
        else:
            resultprint("FAIL")
            sys.exit(1)

        if args.regression or args.all:
            print("regression tests:")
            for p in programs:
                for s in args.solver:
                    for i in args.initcond:
                        if s == "kinetic":
                            tests.KineticRegressionTest(i).run(p)
                        elif s == "moment":
                            for f in args.moment_filter:
                                tests.MomentRegressionTest(f, i).run(p)

        if args.convergence or args.all:
            print("convergence tests:")
            for s in args.solver:
                if s == "kinetic":
                    tests.KineticConvergenceTest().run(args.convergence_iterations)
                elif s == "moment":
                    for f in args.moment_filter:
                        tests.MomentConvergenceTest(f).run(args.convergence_iterations)
