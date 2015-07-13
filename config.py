###### required libraries ######
cpu_libs = ["gsl", "gslcblas", "blas", "lapack", "m"]
cuda_libs = ["cudart"]

###### additional search paths for CPU libraries ######
library_paths = []
include_paths = []

###### C++ options ######
cxx_compiler = "g++"
cxx_flags = ["-Wall", "-Wno-unknown-pragmas", "-O3"]
cxx_debug_flags = ["-g", "-Wall"]
cxx_profile_flags = ["-pg"]

###### OpenMP Options ######
omp_compiler_flags = ["-fopenmp"]
omp_linker_flags = ["-fopenmp"]

###### MPI Options ######
mpi_compiler = "mpic++"

###### CUDA Options ######
nvcc_flags = ["-O3", "-arch=sm_20"]
cuda_toolkit_path = "/opt/cuda"
cuda_sdk_path = "/opt/cuda"
