# Source and library selection derived from the public configuration variables.

GPU_PROVIDER_SRC :=
GPU_DEFINES :=
GPU_INC :=
GPU_LIBDIR :=
GPU_LIBS :=
GPU_LAPACK_INC :=
GPU_LAPACK_LIBDIR :=
GPU_LAPACK_LIBS :=
GPU_MODULE_FLAGS :=

ifeq ($(GPU_MODE),ON)
  ifeq ($(GPU_BACKEND),CUDA)
    # Prefer an explicit CUDA_DIR, then facility environment variables, then a
    # conventional default. XNet links CUDA libraries but compiles no .cu code.
    ifndef CUDA_DIR
      ifdef CUDATOOLKIT_HOME
        CUDA_DIR := $(CUDATOOLKIT_HOME)
      else ifdef CRAY_CUDATOOLKIT_DIR
        CUDA_DIR := $(CRAY_CUDATOOLKIT_DIR)
      else ifdef OLCF_CUDA_ROOT
        CUDA_DIR := $(OLCF_CUDA_ROOT)
      else ifneq ($(filter $(PE_ENV),CCE CRAY),)
        CUDA_DIR := /opt/nvidia/cudatoolkit/default
      else
        CUDA_DIR := /usr/local/cuda
      endif
    endif

    GPU_DEFINES += -DXNET_GPU -DXNET_CUDA
    GPU_INC += -I$(CUDA_DIR)/include
    GPU_LIBDIR += -L$(CUDA_DIR)/lib64
    GPU_LIBS += -lcudart -lcuda
    GPU_PROVIDER_SRC += $(SOURCE_DIR)/cudaf.F90

    ifeq ($(GPU_LAPACK_VER),CUBLAS)
      GPU_DEFINES += -DXNET_LA_CUBLAS
      GPU_LAPACK_INC += -I$(CUDA_DIR)/include
      GPU_LAPACK_LIBDIR += -L$(CUDA_DIR)/lib64
      GPU_LAPACK_LIBS += -lcublas -lcusparse -lcusolver
      GPU_PROVIDER_SRC += $(addprefix $(SOURCE_DIR)/,cublasf.F90 cusparsef.F90 cusolverf.F90)
    else ifeq ($(GPU_LAPACK_VER),MAGMA)
      MAGMA_DIR ?= $(OLCF_MAGMA_ROOT)
      GPU_DEFINES += -DXNET_LA_MAGMA
      GPU_LAPACK_INC += -I$(MAGMA_DIR)/include -I$(CUDA_DIR)/include
      GPU_LAPACK_LIBDIR += -L$(MAGMA_DIR)/lib -L$(CUDA_DIR)/lib64
      GPU_LAPACK_LIBS += -lmagma -lcublas -lcusparse -lcusolver
      GPU_PROVIDER_SRC += $(addprefix $(SOURCE_DIR)/,magmaf.F90 cublasf.F90 cusparsef.F90 cusolverf.F90)
    endif

  else ifeq ($(GPU_BACKEND),HIP)
    HIPFORT_DIR ?= $(OLCF_HIPFORT_ROOT)
    ROCM_DIR ?= $(ROCM_PATH)
    HIPFORT_MODULE_DIR := $(HIPFORT_DIR)/include/hipfort/amdgcn
    ifeq ($(wildcard $(HIPFORT_MODULE_DIR)),)
      $(error HIPFort module directory does not exist: $(HIPFORT_MODULE_DIR))
    endif
    GPU_MODULE_FLAGS := -I$(HIPFORT_MODULE_DIR)
    GPU_DEFINES += -DXNET_GPU -DXNET_HIP -DXNET_LA_ROCM
    GPU_INC += -I$(HIPFORT_MODULE_DIR)
    GPU_LIBDIR += -L$(HIPFORT_DIR)/lib
    GPU_LIBS += -lhipfort-amdgcn
    GPU_LAPACK_INC += -I$(HIPFORT_MODULE_DIR) -I$(ROCM_DIR)/include
    GPU_LAPACK_LIBDIR += -L$(ROCM_DIR)/lib
    GPU_LAPACK_LIBS += -lrocsparse -lrocsolver -lrocblas -lhipblas -lhipsparse -lamdhip64
    GPU_PROVIDER_SRC += $(addprefix $(SOURCE_DIR)/,hipf.F90 hipblasf.F90 rocblasf.F90 rocsparsef.F90 rocsolverf.F90)
  endif

  ifeq ($(OPENACC_MODE),ON)
    GPU_DEFINES += -DXNET_OACC
    GPU_PROVIDER_SRC += $(SOURCE_DIR)/openaccf.F90
    GPU_DIRECTIVE_FLAGS := $(OPENACC)
  else ifeq ($(OPENMP_OL_MODE),ON)
    GPU_DEFINES += -DXNET_OMP_OL
    GPU_PROVIDER_SRC += $(SOURCE_DIR)/openmpf.F90
    GPU_DIRECTIVE_FLAGS := $(OPENMP_OL)
  endif

  FFLAGS += $(GPU_DIRECTIVE_FLAGS)
  LDFLAGS += $(GPU_DIRECTIVE_FLAGS)
  LAPACK_INC += $(GPU_INC) $(GPU_LAPACK_INC)
  LAPACK_LIBDIR += $(GPU_LIBDIR) $(GPU_LAPACK_LIBDIR)
  LAPACK_LIBS += $(GPU_LIBS) $(GPU_LAPACK_LIBS)
endif

ifeq ($(MPI_MODE),ON)
  MPI_SRC := $(SOURCE_DIR)/xnet_parallel.F90
else
  MPI_SRC := $(SOURCE_DIR)/xnet_parallel_stubs.F90
endif

ifeq ($(EOS),STARKILLER)
  STARKILLER_HELMHOLTZ_PATH ?= $(ROOT_DIR)/tools/starkiller-helmholtz
  EOS_SRC := $(abspath $(STARKILLER_HELMHOLTZ_PATH))/actual_eos.F90 \
             $(abspath $(STARKILLER_HELMHOLTZ_PATH))/eos_type.F90 \
             $(SOURCE_DIR)/xnet_eos_starkiller.F90
else ifeq ($(EOS),BAHCALL)
  EOS_SRC := $(SOURCE_DIR)/xnet_eos_bahcall.F90
else ifeq ($(EOS),HELMHOLTZ)
  EOS_SRC := $(SOURCE_DIR)/xnet_eos_helm.F90 $(abspath $(HELMHOLTZ_PATH))/helmholtz.F90
endif

# Jacobian filenames follow MATRIX_SOLVER directly, including PARDISO_MKL.
JAC_SRC := $(SOURCE_DIR)/xnet_jacobian_$(MATRIX_SOLVER).F90

# MA48 and oneMKL PARDISO share the persisted sparse-index reader.
SPARSE_SRC :=
ifneq ($(filter MA48 PARDISO_MKL,$(MATRIX_SOLVER)),)
  SPARSE_SRC := $(SOURCE_DIR)/xnet_sparse.F90
endif

SOLVER_SRC :=
ifeq ($(MATRIX_SOLVER),MA41)
  SOLVER_SRC := $(MA41_SRC)
else ifeq ($(MATRIX_SOLVER),MA48)
  SOLVER_SRC := $(MA48_SRC)
endif
