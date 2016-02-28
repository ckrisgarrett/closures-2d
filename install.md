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

Copy the files in the `examples` directory to the root directory.
This should contain `config.py` and `input.deck`.
In `config.py`, you will need to supply the appropriate paths for the GSL, 
BLAS, and Lapack libraries.
You may also need to change the name of the libraries.
For instance, if you use openblas, you can get rid of `lapack` in cpu libs
and replace `blas` with `openblas`.

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
`solver_mpiomp` by default. Running `scons` builds only the default targets.
Specifying build options (e.g. `--omp`) only makes
the build targets available, it does *not* build them.
To build a particular target, specify it on
the command line. For example, running

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

To add PAPI profiling, add `--papi` to the command line,
e.g. `scons --omp --papi` builds `solver_serial` and `solver_omp`, both
with PAPI compiled in.

Example Build Commands
======================

* `scons` builds `solver_serial`
* `scons --omp solver_omp` builds `solver_omp`
* `scons --mpi solver_mpi` builds `solver_mpi`
* `scons --omp --mpi solver_mpiomp` builds `solver_mpiomp`
* `scons solver_debug` builds `solver_debug`
* `scons solver_profile` builds `solver_profile`




