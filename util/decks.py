def kinetic_deck(quad_order=9,
                 cfl_factor=0.9,
                 use_lebedev=1):
    return (
"""QUAD_ORDER %d
CFL_FACTOR %.17f
USE_LEBEDEV %d""" % (quad_order,
                    cfl_factor,
                    use_lebedev))

def moment_deck(moment_order=3,
                quad_order=30,
                filter_type="none",
                filter_tune=100,
                cfl_factor=0.9):
    filterno = {"none": 0, "hauck": 1, "sspline": 2, "lanczos": 3}
    return ("""
MOMENT_ORDER %d
QUAD_ORDER %d
FILTER_TYPE %d
FILTER_TUNE %d
CFL_FACTOR %.17f
""" % (moment_order,
       quad_order,
       filterno[filter_type],
       filter_tune,
       cfl_factor))

def momopt_deck(moment_order=3, quad_order=30, cfl_factor=0.9, tol=1e-4,
        cond_h_max=1e10, cond_h_max_bfgs=20, max_iter=100, max_bgfs_iter=5,
        use_cg=1, theta=2.0, delta_ppn=1e-10, moment_type="mn", opt_type=0,
        num_cuda_cards=1, num_threads_pre_cuda_card=2):
    momtypenum = {"mn": 0, "ppn": 1}
    return ("""
MOMENT_ORDER %d
QUAD_ORDER %d
CFL_FACTOR %.17f
TOL %.17f
COND_H_MAX %.17f
COND_H_MAX_BFGS %d
MAX_ITER %d
MAX_BFGS_ITER %d
USE_CLEBSCH_GORDAN %d
THETA %.17f
DELTA_PPN %.17f
MOMENT_TYPE %d
OPTIMIZATION_TYPE %d
NUM_CUDA_CARDS %d
NUM_THREADS_PER_CUDA_CARD %d
""" % (moment_order,
       quad_order,
       cfl_factor,
       tol,
       cond_h_max,
       cond_h_max_bfgs,
       max_iter,
       max_bgfs_iter,
       use_cg,
       theta,
       delta_ppn,
       momtypenum[moment_type],
       opt_type,
       num_cuda_cards,
       num_threads_pre_cuda_card))

def input_deck(solver="kinetic",
               num_cells_x=50,
               num_cells_y=50,
               num_mpi_partitions_x=1,
               num_mpi_partitions_y=1,
               a_x=-1.5,
               b_x=1.5,
               a_y=-1.5,
               b_y=1.5,
               t_final=1.0,
               out_delta_t=1000.1,
               gaussian_sigma=0.4,
               floor=1e-4,
               init_cond="delta",
               sigma=1.0):
    condno = {"delta": 0, "gaussian": 0, "lattice": 1, "smooth": 2}
    return ("""
SOLVER %s
NUM_CELLS_X %d
NUM_CELLS_Y %d
NUM_MPI_PARTITIONS_X %d
NUM_MPI_PARTITIONS_Y %d
A_X %.17f
B_X %.17f
A_Y %.17f
B_Y %.17f
T_FINAL %.17f
OUT_DELTA_T %.17f
GAUSSIAN_SIGMA %.17f
FLOOR %.17f
INIT_COND %d
SIGMA %.17f
""" % (solver,
       num_cells_x,
       num_cells_y,
       num_mpi_partitions_x,
       num_mpi_partitions_y,
       a_x,
       b_x,
       a_y,
       b_y,
       t_final,
       out_delta_t,
       gaussian_sigma,
       floor,
       condno[init_cond],
       sigma))
