#
#  Compiler and flags
#


# Intel Fortran compiler
# ----------------------
FC = ifort   # in case -D_QMC_MPI or ifort in case not MPI mode
FC_FLAGS = -m64 -O3 -unroll  -warn all
NOOPT_FLAGS = -O0 -warn all

# MAGMA library path
# ------------------
#MAGMAPATH = ./magma-1.4.0
#MAGMALIB = -L$(MAGMAPATH)/lib -lmagma -lmagmablas -lmagma 
#MAGMAINC = -I$(MAGMAPATH)/include

# CUDA installation path
# ----------------------
#CUDAPATH = /usr/local/cuda
#CUDAPATH = /usr/lib/nvidia-cuda-toolkit/

# CUDA compiler and libraries
# ---------------------------
#NVCC = $(CUDAPATH)/bin/nvcc
#CU_FLAGS = -O3 -Xptxas -v -m 64 -arch sm_20 $(MAGMAINC)
#CUDALIB = -L$(CUDAPATH)/lib64 -lcublas -lcudart -lcuda
#CUDAINC = -I$(CUDAPATH)/include

# Intel C++ compiler
# ------------------
CXX = icpc # mpiicpc
CXX_FLAGS = -m64 -O3 -Wall -unroll# $(CUDAINC) $(MAGMAINC)

# C++ libraries
# -------------
CXXLIB = -lstdc++ #-lrt

# Path to Intel MKL library
# -------------------------
MKLPATH = $(MKLROOT)/lib/intel64
#MKLPATH = /opt/intel/mkl/lib/intel64

LAPACKLIB = -L$(MKLPATH) -mkl:sequential -lpthread 

# DQMC library
# ------------
DQMCLIB    = libdqmc.a

# SPRNG library
# ------------
SPRNGLIB   = libs/libsprng.a 

MPILIB = libs/libmpichf90.la

# Library linking
# ---------------
LIB        = $(CXXLIB) $(LAPACKLIB) #$(SPRNGLIB) $MPILIB # $(CUDALIB) $(MAGMALIB) 

# Archiver and flags
# ------------------
ARCH       = xiar
ARFLAG     = cr
RANLIB     = ranlib


# Program dependent flags
# -----------------------
#    -D_CKB    : use checkboard method for matrix B
#    -D_PREFIX : using prefix matrix multiplication
#    -D_OpenBC : use open boundary condition
#

PRG_FLAGS = #-D_CKB #-D_QMC_MPI #-DDQMC_ASQRD  # -DDQMC_PROFILE -DDQMC_CUDA -D_CKBP -D_SXX

FLAGS = $(FC_FLAGS) $(PRG_FLAGS)
