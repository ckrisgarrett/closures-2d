import os
import config
import testing


binaries = []

src_files = ["main.cpp",
             "utils.cpp",
             "comm.cpp",
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


AddOption("--omp", action="store_true")
AddOption("--mpi", action="store_true")
AddOption("--cuda", action="store_true")
AddOption("--all", action="store_true")


cxx = Environment(CXX=config.cxx_compiler, CCFLAGS=config.cxx_flags,
    LIBS=config.cpu_libs, LIBPATH=config.library_paths, CPPPATH=config.include_paths)

if GetOption("mpi") or GetOption("all"):
    mpic = cxx.Clone(CXX=config.mpi_compiler)

if GetOption("cuda") or GetOption("all"):
    cudac = Environment()
    cudac['CUDA_TOOLKIT_PATH'] = config.cuda_toolkit_path
    cudac['CUDA_SDK_PATH'] = config.cuda_sdk_path
    cudac['NVCCFLAGS'] = config.nvcc_flags
    cudac.Tool("cuda")
    cudac.AppendUnique(LIBS=config.cuda_libs)


serial = cxx.Clone()

debug = cxx.Clone(CCFLAGS=config.cxx_debug_flags)

profile = cxx.Clone()
profile.AppendUnique(CCFLAGS=config.cxx_profile_flags)

if GetOption("omp") or GetOption("all"):
    omp = cxx.Clone()
    omp.AppendUnique(CCFLAGS=config.omp_compiler_flags,
        LINKFLAGS=config.omp_linker_flags)

if GetOption("mpi") or GetOption("all"):
    mpi = mpic.Clone()

if GetOption("cuda") or GetOption("all"):
    cuda = cxx.Clone()
    cuda.AppendUnique(LIBS=config.cuda_libs)

if (GetOption("mpi") and GetOption("omp")) or GetOption("all"):
    mpiomp = mpic.Clone()
    mpiomp.AppendUnique(CCFLAGS=config.omp_compiler_flags,
        LINKFLAGS=config.omp_linker_flags)

if (GetOption("mpi") and GetOption("cuda")) or GetOption("all"):
    mpicuda = mpic.Clone()
    mpicuda.AppendUnique(LIBS=config.cuda_libs)

if (GetOption("omp") and GetOption("cuda")) or GetOption("all"):
    ompcuda = omp.Clone()
    ompcuda.AppendUnique(LIBS=config.cuda_libs)

if (GetOption("mpi") and GetOption("omp") and GetOption("cuda")) or GetOption("all"):
    mpiompcuda = mpicuda.Clone()
    mpiompcuda.AppendUnique(CCFLAGS=config.omp_compiler_flags,
        LINKFLAGS=config.omp_linker_flags)


serial.VariantDir("build/serial", "src")
binaries.append(serial.Program("solver_serial",
    [os.path.join("build/serial", x) for x in src_files]))
convergence_binary = binaries[-1]
Default(binaries[-1])

profile.VariantDir("build/profile", "src")
profile.Program("solver_profile",
    [os.path.join("build/profile", x) for x in src_files])

debug.VariantDir("build/debug", "src")
debug.Program("solver_debug",
    [os.path.join("build/debug", x) for x in src_files])

if GetOption("omp") or GetOption("all"):
    omp.VariantDir("build/omp", "src")
    omp.AppendUnique(CPPDEFINES=["USE_OPENMP"])
    binaries.append(omp.Program("solver_omp",
        [os.path.join("build/omp", x) for x in src_files]))
    convergence_binary = binaries[-1]
    Default(binaries[-1])

if GetOption("mpi") or GetOption("all"):
    mpi.VariantDir("build/mpi", "src")
    mpi.AppendUnique(CPPDEFINES=["USE_MPI"])
    binaries.append(mpi.Program("solver_mpi",
        [os.path.join("build/mpi", x) for x in src_files]))
    Default(binaries[-1])

if GetOption("cuda") or GetOption("all"):
    cudac.VariantDir("build/cudac", "src")
    cudac.AppendUnique(CPPDEFINES=["USE_CUDA_FLUX", "USE_CUDA_FOBJ"])
    fobj_cuda = cudac.Object("build/cudac/momopt/opt/fobj_cuda.cu")
    flux_cuda = cudac.Object("build/cudac/momopt/solve_flux_cuda.cu")

    cuda.VariantDir("build/cuda", "src")
    cuda.AppendUnique(CPPDEFINES=["USE_CUDA_FLUX", "USE_CUDA_FOBJ"])
    binaries.append(cuda.Program("solver_cuda",
        [os.path.join("build/cuda", x) for x in src_files] + fobj_cuda + flux_cuda))
    Default(binaries[-1])

if (GetOption("mpi") and GetOption("omp")) or GetOption("all"):
    mpiomp.VariantDir("build/mpiomp", "src")
    mpiomp.AppendUnique(CPPDEFINES=["USE_OPENMP", "USE_MPI"])
    binaries.append(mpiomp.Program("solver_mpiomp",
        [os.path.join("build/mpiomp", x) for x in src_files]))

if (GetOption("omp") and GetOption("cuda")) or GetOption("all"):
    ompcuda.VariantDir("build/ompcuda", "src")
    ompcuda.AppendUnique(CPPDEFINES=["USE_OPENMP", "USE_CUDA_FLUX", "USE_CUDA_FOBJ"])
    binaries.append(ompcuda.Program("solver_ompcuda",
        [os.path.join("build/ompcuda", x) for x in src_files] + fobj_cuda + flux_cuda))

if (GetOption("mpi") and GetOption("cuda")) or GetOption("all"):
    mpicuda.VariantDir("build/mpicuda", "src")
    mpicuda.AppendUnique(CPPDEFINES=["USE_MPI", "USE_CUDA_FLUX", "USE_CUDA_FOBJ"])
    binaries.append(mpicuda.Program("solver_mpicuda",
        [os.path.join("build/mpicuda", x) for x in src_files] + fobj_cuda + flux_cuda))

if (GetOption("mpi") and GetOption("omp") and GetOption("cuda")) or GetOption("all"):
    mpiompcuda.VariantDir("build/mpiompcuda", "src")
    mpiompcuda.AppendUnique(CPPDEFINES=["USE_MPI", "USE_OMP", "USE_CUDA_FLUX", "USE_CUDA_FOBJ"])
    binaries.append(mpiompcuda.Program("solver_mpiompcuda",
        [os.path.join("build/mpiompcuda", x) for x in src_files] + fobj_cuda + flux_cuda))


tests = Environment()
tests.Append(BUILDERS={"Regression": Builder(action=testing.regression_test)})
tests.Append(BUILDERS={"Convergence": Builder(action=testing.convergence_test)})
tests.Append(BUILDERS={"Everything": Builder(action=testing.test_everything)})
tests.Alias("regression_tests", "regression.txt")
tests.Alias("convergence_tests", "convergence.txt")
tests.Alias("tests", "tests.txt")

tests.Regression("regression.txt", binaries)
tests.Convergence("convergence.txt", convergence_binary)
tests.Everything("tests.txt", ["regression.txt", "convergence.txt"])
