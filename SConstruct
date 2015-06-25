### additional search paths for nonstandard installations, etc. ###
library_paths = []
include_paths = []

### names of required libraries ###
# DO NOT prefix with "l"
cpu_libs = ["gsl", "blas", "cblas", "lapack", "rt"]
cuda_libs = ["cudart"]

# compilation options
AddOption("--omp", action="store_true")
AddOption("--mpi", action="store_true")
AddOption("--cuda", action="store_true")
AddOption("--all", action="store_true")

### Compiler definitions ###
# base compilers
cxx = Environment(CXX="g++", CCFLAGS=["-Wall", "-Wno-unknown-pragmas", "-O3"],
    LIBS=cpu_libs, LIBPATH=library_paths, CPPPATH=include_paths)

if GetOption("mpi") or GetOption("all"):
    mpic = cxx.Clone(CXX="mpic++")

if GetOption("cuda") or GetOption("all"):
    cudac = Environment()
    cudac['CUDA_TOOLKIT_PATH'] = "/opt/cuda"
    cudac['CUDA_SDK_PATH'] = "/opt/cuda"
    cudac['NVCCFLAGS'] = "-O3 -arch=sm_30"
    cudac.Tool("cuda")
    cudac.AppendUnique(LIBS=cuda_libs)

# now customize for various build configurations

serial = cxx.Clone()

debug = cxx.Clone(CCFLAGS=["-g", "-Wall"])

profile = cxx.Clone()
profile.AppendUnique(CCFLAGS=["-pg"])

if GetOption("omp") or GetOption("all"):
    omp = cxx.Clone()
    omp.AppendUnique(CCFLAGS=["-fopenmp"], LINKFLAGS=["-fopenmp"])

if GetOption("mpi") or GetOption("all"):
    mpi = mpic.Clone()

if GetOption("cuda") or GetOption("all"):
    cuda = cxx.Clone()
    cuda.AppendUnique(LIBS=cuda_libs)

if (GetOption("mpi") and GetOption("omp")) or GetOption("all"):
    mpiomp = mpic.Clone()
    mpiomp.AppendUnique(CCFLAGS=["-fopenmp"], LINKFLAGS=["-fopenmp"])

if (GetOption("mpi") and GetOption("cuda")) or GetOption("all"):
    mpicuda = mpic.Clone()
    mpicuda.AppendUnique(LIBS=cuda_libs)

if (GetOption("omp") and GetOption("cuda")) or GetOption("all"):
    ompcuda = omp.Clone()
    ompcuda.AppendUnique(LIBS=cuda_libs)

if (GetOption("mpi") and GetOption("omp") and GetOption("cuda")) or GetOption("all"):
    mpiompcuda = mpicuda.Clone()
    mpiompcuda.AppendUnique(CCFLAGS=["-fopenmp"], LINKFLAGS=["-fopenmp"])

###### end of user-configurable build commands ######

import os

src_files = ["main.cpp",
             "utils.cpp",
             "comm.cpp",
             "timer.cpp",
             "moment/moment_boundaries.cpp",
             "moment/moment_init.cpp",
             "moment/moment_output.cpp",
             "moment/moment_update.cpp",
             "kinetic/kinetic_boundaries.cpp",
             "kinetic/kinetic_init.cpp",
             "kinetic/kinetic_output.cpp",
             "kinetic/kinetic_update.cpp",
             "momopt/momopt_boundaries.cpp",
             "momopt/momopt_init.cpp",
             "momopt/momopt_output.cpp",
             "momopt/momopt_update.cpp",
             "momopt/opt/fobj.cpp",
             "momopt/opt/opt.cpp",
             "momopt/opt/linesearch.cpp",
             "momopt/opt/optcb.cpp",
             "momopt/opt/optbfgs.cpp",
             "dn/dn_boundaries.cpp",
             "dn/dn_init.cpp",
             "dn/dn_output.cpp",
             "dn/dn_update.cpp"]

# build targets
serial.VariantDir("build/serial", "src")
Default(serial.Program("solver_serial",
    [os.path.join("build/serial", x) for x in src_files]))

profile.VariantDir("build/profile", "src")
profile.Program("solver_profile",
    [os.path.join("build/profile", x) for x in src_files])

debug.VariantDir("build/debug", "src")
debug.Program("solver_debug",
    [os.path.join("build/debug", x) for x in src_files])

if GetOption("omp") or GetOption("all"):
    omp.VariantDir("build/omp", "src")
    omp.AppendUnique(CPPDEFINES=["USE_OPENMP"])
    Default(omp.Program("solver_omp",
        [os.path.join("build/omp", x) for x in src_files]))

if GetOption("mpi") or GetOption("all"):
    mpi.VariantDir("build/mpi", "src")
    mpi.AppendUnique(CPPDEFINES=["USE_MPI"])
    Default(mpi.Program("solver_mpi",
        [os.path.join("build/mpi", x) for x in src_files]))

if GetOption("cuda") or GetOption("all"):
    cudac.VariantDir("build/cudac", "src")
    cudac.AppendUnique(CPPDEFINES=["USE_CUDA_FLUX", "USE_CUDA_FOBJ"])
    fobj_cuda = cudac.Object("build/cudac/momopt/opt/fobj_cuda.cu")
    flux_cuda = cudac.Object("build/cudac/momopt/solve_flux_cuda.cu")

    cuda.VariantDir("build/cuda", "src")
    cuda.AppendUnique(CPPDEFINES=["USE_CUDA_FLUX", "USE_CUDA_FOBJ"])
    Default(cuda.Program("solver_cuda",
        [os.path.join("build/cuda", x) for x in src_files] + fobj_cuda + flux_cuda))

if (GetOption("mpi") and GetOption("omp")) or GetOption("all"):
    mpiomp.VariantDir("build/mpiomp", "src")
    mpiomp.AppendUnique(CPPDEFINES=["USE_OPENMP", "USE_MPI"])
    mpiomp.Program("solver_mpiomp",
        [os.path.join("build/mpiomp", x) for x in src_files])

if (GetOption("omp") and GetOption("cuda")) or GetOption("all"):
    ompcuda.VariantDir("build/ompcuda", "src")
    ompcuda.AppendUnique(CPPDEFINES=["USE_OPENMP", "USE_CUDA_FLUX", "USE_CUDA_FOBJ"])
    ompcuda.Program("solver_ompcuda",
        [os.path.join("build/ompcuda", x) for x in src_files] + fobj_cuda + flux_cuda)

if (GetOption("mpi") and GetOption("cuda")) or GetOption("all"):
    mpicuda.VariantDir("build/mpicuda", "src")
    mpicuda.AppendUnique(CPPDEFINES=["USE_MPI", "USE_CUDA_FLUX", "USE_CUDA_FOBJ"])
    mpicuda.Program("solver_mpicuda",
        [os.path.join("build/mpicuda", x) for x in src_files] + fobj_cuda + flux_cuda)

if (GetOption("mpi") and GetOption("omp") and GetOption("cuda")) or GetOption("all"):
    mpiompcuda.VariantDir("build/mpiompcuda", "src")
    mpiompcuda.AppendUnique(CPPDEFINES=["USE_MPI", "USE_OMP", "USE_CUDA_FLUX", "USE_CUDA_FOBJ"])
    mpiompcuda.Program("solver_mpiompcuda",
        [os.path.join("build/mpiompcuda", x) for x in src_files] + fobj_cuda + flux_cuda)
