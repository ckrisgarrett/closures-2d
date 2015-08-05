Dependencies
============

This software depends on

* GSL - GNU Scientific Library
* BLAS
* LAPACK
* OpenMP (optional, provided by many compilers)
* Open MPI (optional)
* PAPI (optional)

If the required libraries are installed at nonstandard locations, or
additional compiler customization is necessary, users should edit `config.py`.
This Python file declares several common build options, such as library
paths, include paths, compiler flags, etc.

Building
========

This software is built using [SCons](http://scons.org/). In case SCons is not
installed locally, a copy is included in the `build/` directory. To use it,
ensure that Python 2.x is installed and
replace `scons` in the following build commands with `python build/scons.py`.

To build the basic configuration, simply run

    scons

This will build `solver_serial`, the serial version without any acceleration
or other extras. To clean up build files, run

    scons -c

To use MPI and/or OMP acceleration, pass SCons the appropriate flag on the
command line (`--mpi` and `--omp`). Note that both can be combined in the
same build.

    scons --omp --mpi

In this case, SCons will build `solver_serial`, `solver_omp`, `solver_mpi`, and
`solver_mpiomp` by default. To build only a particular target, specify it on
the command line. Note that build targets are only available with the
corresponding flag. For example, running

    scons solver_mpi

fails with

    scons: Reading SConscript files ...
    scons: done reading SConscript files.
    scons: Building targets ...
    scons: *** Do not know how to make File target `mpi' (closures-2d/mpi).  Stop.
    scons: building terminated because of errors.

To build only the MPI and OMP accelerated version, run

    scons --omp --mpi solver_mpiomp

Non-default variants (`solver_profile` and `solver_debug`) can also be built in
this way.


Issues
======

* Where is input.deck?
* Does this work on mac?
* Name papi build something different.
* How to build debug?
* Need to give every parameter available.
* Where are plot.py and other python tools?

