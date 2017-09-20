OPS = $(OPS_INSTALL_PATH)
F_OPS = $(OPS)/fortran
F_OPS_LIB = $(F_OPS)/lib
F_OPS_MOD = $(F_OPS)/mod/$(OPS_COMPILER)
F_OPS_INC = $(F_OPS)/include

# flags for nvcc
# set NV_ARCH to select the correct one
ifndef NV_ARCH
  MESSAGE=select an NVIDA device to compile in CUDA, e.g. make NV_ARCH=KEPLER
  NV_ARCH=Kepler
endif
ifeq ($(NV_ARCH),Fermi)
  CODE_GEN_CUDA=-gencode arch=compute_20,code=sm_21
  PGI_CUDA_FORT_FLAGS=-Mcuda=cuda7.5,cc20
else
ifeq ($(NV_ARCH),Kepler)
  CODE_GEN_CUDA=-gencode arch=compute_35,code=sm_35
  PGI_CUDA_FORT_FLAGS=-Mcuda=cuda7.5,cc35
endif
endif

ifdef HDF5_INSTALL_PATH
  HDF5_INC              := -I$(HDF5_INSTALL_PATH)/include
  HDF5_LIB              := -L$(HDF5_INSTALL_PATH)/lib -lhdf5_hl -lhdf5 -lz
ifdef CHECKPOINTING
  HDF5_INC              := $(HDF5_INC) -DCHECKPOINTING
endif
endif


ifeq ($(OPS_COMPILER),intel)
  CC 		:= icc -std=c99
ifdef DEBUG
  #OPT	 	:= -O3 -g -xAVX -vec-report2 #-xAVX
  OPT	 	:= -O0 -g -traceback -check bounds
else
  OPT	 	:= -O3 -qopt-report #-xAVX
endif
  OMPFLAGS  	= -qopenmp #-openmp-report2
  CCFLAGS       := $(OPT) #-Wall -Werror -wd981
  CXX 		:= icpc
  CXXFLAGS 	:= $(CCFLAGS)
  CXXLINK	:= -lstdc++
  FC		:= ifort
  FFLAGS	:= $(OPT)
  MPICC 	:= $(MPI_INSTALL_PATH)/bin/mpicc -std=c99
  MPICXX 	:= $(MPI_INSTALL_PATH)/bin/mpicxx
  MPIF90 	:= $(MPI_INSTALL_PATH)/bin/mpif90
  MPIFLAGS 	:= $(CXXFLAGS) -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
  TARGETS	:= clean cfd4t_seq cfd4t_mpi cfd4t_openmp cfd4t_mpi_openmp cfd4t_cuda cfd4t_mpi_cuda #cfd4t_openacc cfd4t_mpi_openacc
else
ifeq ($(OPS_COMPILER),pgi)
ifdef DEBUG
  OPT           := -O0 -g -Minline -Kieee #-Minform=inform -Minfo=all
else
  OPT           := -O3 -Kieee -Minline #-Minform=inform -Minfo=all
endif
  OMPFLAGS      := -mp
  OpenACCFLAGS	:= -acc -ta=tesla:cc35 -ta=tesla:maxregcount:64 #-Minline
  CCFLAGS       := $(OPT) #-Wall -Werror -wd981
  CXX           := pgc++
  CXXFLAGS      := $(CCFLAGS)
  CUDA_ALIGN_FLAG := -D__x86_64 -D__align__\(n\)=__attribute__\(\(aligned\(n\)\)\) -D__location__\(a\)=__annotate__\(a\) -DCUDARTAPI=
  CXXLINK       := -lstdc++ -lgcc_s -pgcpplibs
  FC            := pgfortran $(PGI_CUDA_FORT_FLAGS)
  FFLAGS        := $(OPT) -ldl
  MPICC         := $(MPI_INSTALL_PATH)/bin/mpicc -std=c99
  MPICXX        := $(MPI_INSTALL_PATH)/bin/mpicxx
  MPIF90 	:= $(MPI_INSTALL_PATH)/bin/mpif90 $(PGI_CUDA_FORT_FLAGS)
  MPIFLAGS      := $(FFLAGS)
  CUDAFOR       = -DOPS_WITH_CUDAFOR
  TARGETS       := clean cfd4t_seq cfd4t_mpi cfd4t_openmp cfd4t_mpi_openmp cfd4t_cuda cfd4t_mpi_cuda #cfd4t_openacc cfd4t_mpi_openacc
else
print:
	@echo "unrecognised value for OPS_COMPILER"
endif
endif

FLINK = -L$(F_OPS_LIB) $(CXXLINK)
FMODS = -I$(F_OPS_MOD) -I$(F_OPS_INC)
FMODS_OPENACC = -I$(F_OPS_MOD)/openacc -I$(F_OPS_INC)

all: $(TARGETS)

cfd4t_seq: main_ops.F90 constants.F90
	$(FC) $(OPT) $(FMODS) -c constants.F90 ./MPI/initialise_kernel_seq_kernel.F90                         \
	./MPI/conres_kernel_seq_kernel.F90 ./MPI/momres_kernel_seq_kernel.F90                                 \
	./MPI/energyres_kernel_seq_kernel.F90 ./MPI/save_kernel_seq_kernel.F90                                \
	./MPI/calc_dt_kernel_seq_kernel.F90                                                                   \
	main_ops.F90
	$(MPIF90) $(OPT) $(FFLAGS) $(FLINK) main_ops.o constants.o                                            \
	initialise_kernel_seq_kernel.o conres_kernel_seq_kernel.o momres_kernel_seq_kernel.o                  \
	energyres_kernel_seq_kernel.o save_kernel_seq_kernel.o                                                \
	calc_dt_kernel_seq_kernel.o                                                                           \
	 -o ../bin/cfd4t_seq -lops_for_seq -lops_hdf5 $(HDF5_LIB)

cfd4t_openmp: main_ops.F90 constants.F90
	$(FC) $(OPT) $(FMODS) $(OMPFLAGS) -c constants.F90 ./MPI_OpenMP/initialise_kernel_omp_kernel.F90      \
	./MPI_OpenMP/dsdx_kernel_omp_kernel.F90 ./MPI_OpenMP/dvdx_kernel_omp_kernel.F90                       \
	main_ops.F90
	$(MPIF90) $(OPT) $(FLINK) $(FFLAGS) $(OMPFLAGS) main_ops.o constants.o                                \
	initialise_kernel_omp_kernel.o dsdx_kernel_omp_kernel.o dvdx_kernel_omp_kernel.o                      \
    -o ../bin/cfd4t_openmp -lops_for_seq -lops_hdf5 $(HDF5_LIB)

cfd4t_cuda: main_ops.F90 constants.F90
ifeq ($(OPS_COMPILER),pgi)
	$(FC) $(OPT) $(FMODS) $(CUDAFOR) -c constants.F90 ./CUDA/initialise_kernel_cuda_kernel.CUF            \
	./CUDA/dsdx_kernel_cuda_kernel.CUF ./CUDA/dvdx_kernel_cuda_kernel.CUF                                 \
	main_ops.F90
	$(MPIF90) $(OPT) $(FLINK) $(FFLAGS) $(CUDAFOR) main_ops.o constants.o initialise_kernel_cuda_kernel.o \
	dsdx_kernel_cuda_kernel.o dvdx_kernel_cuda_kernel.o                                                   \
    -o cfd4t_cuda -lops_for_cuda -lops_hdf5 $(HDF5_LIB)
else
	@echo "cuda cannot be built with $(OPS_COMPILER) compilers .. requires a CUDA FORTRAN compiler (e.g. pgi)"
endif

cfd4t_mpi: main_ops.F90 constants.F90
	$(MPIF90) $(OPT) $(FMODS) -DOPS_MPI -c constants.F90 ./MPI/initialise_kernel_seq_kernel.F90           \
	./MPI/conres_kernel_seq_kernel.F90 ./MPI/momres_kernel_seq_kernel.F90                                 \
	./MPI/energyres_kernel_seq_kernel.F90 ./MPI/save_kernel_seq_kernel.F90                                \
	./MPI/calc_dt_kernel_seq_kernel.F90                                                                   \
	main_ops.F90
	$(MPIF90) $(OPT) $(MPIFLAGS) $(FLINK) -DOPS_MPI main_ops.o constants.o                                \
	initialise_kernel_seq_kernel.o conres_kernel_seq_kernel.o momres_kernel_seq_kernel.o                  \
	energyres_kernel_seq_kernel.o save_kernel_seq_kernel.o                                                \
	calc_dt_kernel_seq_kernel.o                                                                           \
	-o ../bin/cfd4t_mpi -lops_for_mpi $(HDF5_LIB)

cfd4t_mpi_openmp: main_ops.F90 constants.F90
	$(MPIF90) $(OPT) $(FMODS) $(OMPFLAGS) -DOPS_MPI -c constants.F90                                      \
	./MPI_OpenMP/initialise_kernel_omp_kernel.F90 ./MPI_OpenMP/dsdx_kernel_omp_kernel.F90                 \
	./MPI_OpenMP/dvdx_kernel_omp_kernel.F90                                                               \
	main_ops.F90
	$(MPIF90) $(OPT) $(MPIFLAGS) $(FLINK) $(OMPFLAGS) -DOPS_MPI main_ops.o constants.o                    \
	initialise_kernel_omp_kernel.o dsdx_kernel_omp_kernel.o dvdx_kernel_omp_kernel.o                      \
    -o ../bin/cfd4t_mpi_openmp -lops_for_mpi $(HDF5_LIB)

cfd4t_mpi_cuda: main_ops.F90 constants.F90
ifeq ($(OPS_COMPILER),pgi)
	$(MPIF90) $(OPT) $(FMODS) -DOPS_MPI $(CUDAFOR) -c constants.F90 ./CUDA/initialise_kernel_cuda_kernel.CUF \
	main_ops.F90
	./CUDA/dsdx_kernel_cuda_kernel.CUF ./CUDA/dvdx_kernel_cuda_kernel.CUF                                    \
	$(MPIF90) $(OPT) $(MPIFLAGS) $(FLINK) -DOPS_MPI $(CUDAFOR) main_ops.o constants.o                        \
	initialise_kernel_cuda_kernel.o                                                                          \
	dsdx_kernel_cuda_kernel.o dvdx_kernel_cuda_kernel.o                                                      \
    -o cfd4t_mpi_cuda -lops_for_mpi_cuda $(HDF5_LIB)
else
	@echo "cuda cannot be built with $(OPS_COMPILER) compilers .. requires a CUDA FORTRAN compiler (e.g. pgi)"
endif

cfd4t_openacc: main_ops.F90 constants.F90
	$(FC) $(OPT) $(FMODS_OPENACC) $(OpenACCFLAGS) -c constants.F90                                          \
	./MPI_OpenACC/initialise_kernel_openacc_kernel.F90                                                      \
	./MPI_OpenACC/dsdx_kernel_openacc_kernel.F90 ./MPI_OpenACC/dvdx_kernel_openacc_kernel.F90               \
	main_ops.F90
	$(MPIF90) $(OPT) $(FLINK) $(OpenACCFLAGS) main_ops.o constants.o initialise_kernel_openacc_kernel.o \
    -o ../bin/cfd4t_openacc -lops_for_cuda -lops_hdf5 $(HDF5_LIB)

cfd4t_mpi_openacc: main_ops.F90 constants.F90
	$(MPIF90) $(OPT) $(FMODS_OPENACC) $(OpenACCFLAGS) -DOPS_MPI -c constants.F90                            \
	./MPI_OpenACC/initialise_kernel_openacc_kernel.F90                                                      \
	./MPI_OpenACC/dsdx_kernel_openacc_kernel.F90 ./MPI_OpenACC/dvdx_kernel_openacc_kernel.F90               \
	main_ops.F90
	$(MPIF90) $(OPT) $(FLINK) $(OpenACCFLAGS) -DOPS_MPI main_ops.o constants.o initialise_kernel_openacc_kernel.o \
	dsdx_kernel_cuda_kernel.o dvdx_kernel_cuda_kernel.o                                                      \
    -o ../bin/cfd4t_mpi_openacc -lops_for_mpi_cuda $(HDF5_LIB)

clean:
	rm -f *.o
	rm -f *.mod
	rm -f $(EXEC)
	rm -f *~
	rm -f ../bin/*

