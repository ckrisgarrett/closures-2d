# The compile commands
CC_SERIAL  = g++ -Wall -O3
CC_DEBUG   = g++ -Wall -g
CC_PROFILE = g++ -Wall -pg -O3
CC_CUDA    = g++ -Wall -O3 -DUSE_CUDA_FLUX -DUSE_CUDA_FOBJ
CC_OMP     = g++ -Wall -O3 -fopenmp -DUSE_OPENMP
#CC_OMPCUDA = g++ -Wall -O3 -fopenmp -DUSE_OPENMP -DUSE_CUDA_FLUX -DUSE_CUDA_FOBJ
CC_OMPCUDA = g++ -Wall -O3 -fopenmp -DUSE_OPENMP -DUSE_CUDA_FLUX
CC_CUDA    = nvcc -O3 -arch=sm_30 -DUSE_CUDA_FLUX -DUSE_CUDA_FOBJ


# Source files
MOMENT_SRC_FILES = src/moment/moment_boundaries.cpp src/moment/moment_init.cpp \
                   src/moment/moment_output.cpp src/moment/moment_update.cpp
KINETIC_SRC_FILES = src/kinetic/kinetic_boundaries.cpp src/kinetic/kinetic_init.cpp \
                    src/kinetic/kinetic_output.cpp src/kinetic/kinetic_update.cpp
MOMOPT_SRC_FILES = src/momopt/momopt_boundaries.cpp src/momopt/momopt_init.cpp \
                   src/momopt/momopt_output.cpp src/momopt/momopt_update.cpp \
                   src/momopt/opt/fobj.cpp src/momopt/opt/opt.cpp \
                   src/momopt/opt/linesearch.cpp src/momopt/opt/optcb.cpp \
                   src/momopt/opt/optbfgs.cpp
DN_SRC_FILES = src/dn/dn_boundaries.cpp src/dn/dn_init.cpp \
                   src/dn/dn_output.cpp src/dn/dn_update.cpp
SRC_FILES = src/main.cpp src/utils.cpp src/comm.cpp src/timer.cpp \
            $(MOMENT_SRC_FILES) $(KINETIC_SRC_FILES) \
            $(MOMOPT_SRC_FILES) $(DN_SRC_FILES)

SRC_FILE_FOBJ_CUDA = src/momopt/opt/fobj_cuda.cu
OBJ_FILE_FOBJ_CUDA = fobj_cuda.o
SRC_FILE_FLUX_CUDA = src/momopt/solve_flux_cuda.cu
OBJ_FILE_FLUX_CUDA = solve_flux_cuda.o
OBJ_FILES_CUDA = $(OBJ_FILE_FOBJ_CUDA) $(OBJ_FILE_FLUX_CUDA)


# Libraries
CPU_LIB  = -lgsl -lblas -llapack -lrt
CUDA_LIB = -lcudart


# Build commands
.PHONY: default all
default:
	$(MAKE) omp
all:
	$(MAKE) serial
	$(MAKE) profile
	$(MAKE) debug
	$(MAKE) omp
	$(MAKE) ompcuda


.PHONY: serial
serial: 
	$(CC_SERIAL) -o solver_serial.x $(SRC_FILES) $(CPU_LIB)

.PHONY: profile
profile: 
	$(CC_PROFILE) -o solver_profile.x $(SRC_FILES) $(CPU_LIB)

.PHONY: debug
debug: 
	$(CC_DEBUG) -o solver_debug.x $(SRC_FILES) $(CPU_LIB)

.PHONY: omp
omp: 
	$(CC_OMP) -o solver_omp.x $(SRC_FILES) $(CPU_LIB)

.PHONY: cuda
cuda: 
	$(CC_CUDA) -c $(SRC_FILE_FOBJ_CUDA) -o $(OBJ_FILE_FOBJ_CUDA) $(INCLUDE)
	$(CC_CUDA) -c $(SRC_FILE_FLUX_CUDA) -o $(OBJ_FILE_FLUX_CUDA) $(INCLUDE)
	$(CC_CUDA) -o solver_cuda.x $(SRC_FILES) $(OBJ_FILES_CUDA) $(CPU_LIB) $(CUDA_LIB)

.PHONY: ompcuda
ompcuda: 
	$(CC_CUDA) -c $(SRC_FILE_FOBJ_CUDA) -o $(OBJ_FILE_FOBJ_CUDA) $(INCLUDE)
	$(CC_CUDA) -c $(SRC_FILE_FLUX_CUDA) -o $(OBJ_FILE_FLUX_CUDA) $(INCLUDE)
	$(CC_OMPCUDA) -o solver_ompcuda.x $(SRC_FILES) $(OBJ_FILES_CUDA) $(CPU_LIB) $(CUDA_LIB)


.PHONY: clean
clean:
	rm -f *.o

