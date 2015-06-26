def kinetic_deck(quad_order=12,
                 cfl_factor=0.9):
    return (
"""QUAD_ORDER %d
CFL_FACTOR %f""" % (quad_order,
                    cfl_factor))

def moment_deck(moment_order=3,
                quad_order=30,
                filter_type="none",
                filter_tune=100,
                cfl_factor=0.9):
    filterno = {"none": 0, "hauck": 1, "sspline": 2, "lanczos": 3}
    return (
"""MOMENT_ORDER %d
QUAD_ORDER %d
FILTER_TYPE %d
FILTER_TUNE %d
CFL_FACTOR %f""" % (moment_order,
                    quad_order,
                    filterno[filter_type],
                    filter_tune,
                    cfl_factor))

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
    return (
"""SOLVER %s
NUM_CELLS_X %d
NUM_CELLS_Y %d
NUM_MPI_PARTITIONS_X %d
NUM_MPI_PARTITIONS_Y %d
A_X %f
B_X %f
A_Y %f
B_Y %f
T_FINAL %f
OUT_DELTA_T %f
GAUSSIAN_SIGMA %f
FLOOR %f
INIT_COND %d
SIGMA %f""" % (solver,
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
