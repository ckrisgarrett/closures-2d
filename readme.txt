****************************************
    Description of Files
****************************************

- Makefile.gnu Makefile.titan Makefile.yona
These are makefiles for: a general gnu c++ compiler setup for one machine, titan which has
MPI and OpenMP, and Yona which has MPI, OpenMP, and CUDA capabilities.  I always do the
command "cp Makefile.gnu Makefile" if I want to use the GNU compiler so all I have to do is
type make to make the entire project.

Look at the makefile to decide which version you want to make.  The easiest from Makefile.gnu
is to use just use make which will create solver_omp.x which is an openmp solver.

Also note that you may have to change the parameters for the location of the GSL library
for your makefile.


- input.deck: General parameters for the program
- kinetic.deck: Parameters for Sn
- moment.deck: Parameters for Pn/FPn
- momopt.deck: Parameters for Mn/PPn


- plot_pn.py, plot_sn.py
Python scripts to plot the results of pn/sn runs.  plot_pn works for Mn/PPn too.


- put_data_together.py
This python script puts the data output of an MPI implementation of the entropy
program into 1 file.  In general, the data is output at several times and with
the following names out_<time>_<node>.dat.  The output of the python script is
the .dat files associated with one time put together with an output name of
the type outall_<time>.dat.  To use the script type:

python put_data_together.py <time1> <time2> ...


- src
This directory contains all the source files.


- yona.run, titan.run
Example run scripts for running a program on Yona or Titan using MPI/OpenMP.  Usage is:

qsub yona.run    or     qsub titan.run

Output is sent to output.txt and you can view program progress by running
tail -F output.txt
which shows the last ten lines of output.txt in real time.





****************************************
    Description of Deck Files
****************************************

First, you may not change the order or put comments in the deck files.  The code cannot
handle this right now!!!


--- input.deck ---
SOLVER                  Can be "moment, momopt, kinetic"
NUM_CELLS_X             Size of grid in x-direction
NUM_CELLS_Y             Size of grid in y-direction
NUM_MPI_PARTITIONS_X    Number of MPI tasks in x-direction
NUM_MPI_PARTITIONS_Y    Number of MPI tasks in y-direction
A_X                     Domain is [A_X, B_X] x [A_Y, B_Y]
B_X
A_Y
B_Y
T_FINAL                 Final time
OUT_DELTA_T             Times to output data
GAUSSIAN_SIGMA          Sigma for gaussian initial condition (0.0 gives linesource init cond)
FLOOR                   Smallest value for initial condition
INIT_COND               Should be 0.  All other initial conditions are experimental.
SIGMA                   Total/Scattering cross-section


--- kinetic.deck, moment.deck ---
QUAD_ORDER              Order of quadrature (number of points on z-axis)
CFL_FACTOR              Multiply cfl by this for max time step size
MOMENT_ORDER            The n for Pn
FILTER_TYPE             0 gives Pn, 1,2,3 give different FPn
FILTER_TUNE             Only used for FPn


--- momopt.deck (parameters not in above) ---
TOL                     Tolerance for Mn/PPn convergence
COND_H_MAX              Condition number maximum for Hessian
COND_H_MAX_BFGS         Condition number maximum for Hessian in BFGS sub iteration
MAX_ITER                Maximum number of iterations before regularization
MAX_BFGS_ITER           Max iterations for BFGS sub iterations
USE_CLEBSCH_GORDAN      1 Use Clebsch-Gordan, 0 Don't.  1 can only be used for Optimization Type 0
THETA                   Range 1.0 - 2.0.  For slope minmod.
DELTA_PPN               Delta for PPn under the sqrt
MOMENT_TYPE             0 = Mn, 1 = PPn
OPTIMIZATION_TYPE       0 = normal, 1 = change of basis, 2 = BFGS w/ change of basis
FILTER                  Always set this to 0 (This was an experiment)
NUM_CUDA_CARDS          Can ignore this for non-GPU runs
NUM_THREADS_PER_CUDA_CARD Can ignore this for non-GPU runs


