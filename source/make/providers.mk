# MPI, EOS, solver, accelerator, and numerical-library provider selection.
# The accelerator implementation is selected explicitly for the production build.
override GPU_PROVIDER_SRC :=
override GPU_DEFINES :=
override GPU_INC :=
override GPU_LIBDIR :=
override GPU_LIBS :=
override GPU_LAPACK_INC :=
override GPU_LAPACK_LIBDIR :=
override GPU_LAPACK_LIBS :=
SELECTED_EXTERNAL_MODULE_DIRS := $(XNET_EXTERNAL_MODULE_DIRS)
ifeq ($(GPU_MODE),ON)
  ifeq ($(GPU_BACKEND),CUDA)
    CUDA_DIR ?= /usr/local/cuda
    override GPU_DEFINES += -DXNET_GPU -DXNET_CUDA
    override GPU_INC += -I$(CUDA_DIR)/include
    override GPU_LIBDIR += -L$(CUDA_DIR)/lib64
    override GPU_LIBS += -lcudart -lcuda
    override GPU_PROVIDER_SRC += $(XNET_DIR)/cudaf.F90
    ifeq ($(GPU_LAPACK_VER),CUBLAS)
      override GPU_DEFINES += -DXNET_LA_CUBLAS
      override GPU_LAPACK_INC += -I$(CUDA_DIR)/include
      override GPU_LAPACK_LIBDIR += -L$(CUDA_DIR)/lib64
      override GPU_LAPACK_LIBS += -lcublas -lcusparse -lcusolver
      override GPU_PROVIDER_SRC += $(addprefix $(XNET_DIR)/,cublasf.F90 cusparsef.F90 cusolverf.F90)
    else
      MAGMA_DIR ?= $(OLCF_MAGMA_ROOT)
      override GPU_DEFINES += -DXNET_LA_MAGMA
      override GPU_LAPACK_INC += -I$(MAGMA_DIR)/include -I$(CUDA_DIR)/include
      override GPU_LAPACK_LIBDIR += -L$(MAGMA_DIR)/lib -L$(CUDA_DIR)/lib64
      override GPU_LAPACK_LIBS += -lmagma -lcublas -lcusparse -lcusolver
      override GPU_PROVIDER_SRC += $(addprefix $(XNET_DIR)/,magmaf.F90 cublasf.F90 cusparsef.F90 cusolverf.F90)
    endif
  else
    HIPFORT_DIR ?= $(OLCF_HIPFORT_ROOT)
    ROCM_DIR ?= $(ROCM_PATH)
    SELECTED_EXTERNAL_MODULE_DIRS += $(HIPFORT_DIR)/include/hipfort/amdgcn
    override GPU_DEFINES += -DXNET_GPU -DXNET_HIP -DXNET_LA_ROCM
    override GPU_INC += -I$(HIPFORT_DIR)/include/hipfort/amdgcn
    override GPU_LIBDIR += -L$(HIPFORT_DIR)/lib
    override GPU_LIBS += -lhipfort-amdgcn
    override GPU_LAPACK_INC += -I$(HIPFORT_DIR)/include/hipfort/amdgcn -I$(ROCM_DIR)/include
    override GPU_LAPACK_LIBDIR += -L$(ROCM_DIR)/lib
    override GPU_LAPACK_LIBS += -lrocsparse -lrocsolver -lrocblas -lhipblas -lhipsparse -lamdhip64
    override GPU_PROVIDER_SRC += $(addprefix $(XNET_DIR)/,hipf.F90 hipblasf.F90 rocblasf.F90 rocsparsef.F90 rocsolverf.F90)
  endif
  ifeq ($(OPENACC_MODE),ON)
    override GPU_DEFINES += -DXNET_OACC
    override GPU_PROVIDER_SRC += $(XNET_DIR)/openaccf.F90
    GPU_DIRECTIVE_FLAGS := $(OPENACC)
  else
    override GPU_DEFINES += -DXNET_OMP_OL
    override GPU_PROVIDER_SRC += $(XNET_DIR)/openmpf.F90
    GPU_DIRECTIVE_FLAGS := $(OPENMP_OL)
  endif
  FFLAGS += $(GPU_DIRECTIVE_FLAGS)
  LDFLAGS += $(GPU_DIRECTIVE_FLAGS)
  LAPACK_INC += $(GPU_INC) $(GPU_LAPACK_INC)
  LAPACK_LIBDIR += $(GPU_LIBDIR) $(GPU_LAPACK_LIBDIR)
  LAPACK_LIBS += $(GPU_LIBS) $(GPU_LAPACK_LIBS)
endif

MISSING_SELECTED_EXTERNAL_MODULE_DIRS := $(filter-out $(wildcard $(SELECTED_EXTERNAL_MODULE_DIRS)),$(SELECTED_EXTERNAL_MODULE_DIRS))
ifneq ($(strip $(MISSING_SELECTED_EXTERNAL_MODULE_DIRS)),)
  $(error selected external module directory does not exist: $(MISSING_SELECTED_EXTERNAL_MODULE_DIRS))
endif

ifeq ($(MPI_MODE),ON)
  override MPI_SRC := $(XNET_DIR)/xnet_parallel.F90
else
  override MPI_SRC := $(XNET_DIR)/xnet_parallel_stubs.F90
endif
ifeq ($(EOS),STARKILLER)
  STARKILLER_HELMHOLTZ_PATH ?= $(ROOT_DIR)/tools/starkiller-helmholtz
  override EOS_SRC := $(abspath $(STARKILLER_HELMHOLTZ_PATH))/actual_eos.F90 \
             $(abspath $(STARKILLER_HELMHOLTZ_PATH))/eos_type.F90 $(XNET_DIR)/xnet_eos_starkiller.F90
else ifeq ($(EOS),BAHCALL)
  override EOS_SRC := $(XNET_DIR)/xnet_eos_bahcall.F90
else
  override EOS_SRC := $(XNET_DIR)/xnet_eos_helm.F90 $(abspath $(HELMHOLTZ_PATH))/helmholtz.F90
endif
ifeq ($(MATRIX_SOLVER),dense)
  override JAC_SRC := $(XNET_DIR)/xnet_jacobian_dense.F90
else ifeq ($(MATRIX_SOLVER),PARDISO)
  override JAC_SRC := $(XNET_DIR)/xnet_jacobian_PARDISO_MKL.F90
else
  override JAC_SRC := $(XNET_DIR)/xnet_jacobian_$(MATRIX_SOLVER).F90
endif
override SOLVER_SRC := $(if $(filter MA41,$(MATRIX_SOLVER)),$(MA41_DIR)/MA41.f,\
  $(if $(filter MA48,$(MATRIX_SOLVER)),$(MA48_DIR)/MA48.f))
XNET_NETLIB_FILES := dcopy.f ddot.f dlamch.f dnrm2.f dscal.f dswap.f idamax.f lsame.f xerbla.f \
  daxpy.f dger.f dgemv.f dtrmv.f $(if $(filter MA48,$(MATRIX_SOLVER)),dtrsv.f) \
  dgemm.f dtrmm.f dtrsm.f ieeeck.f iparmq.f ilaenv.f iladlc.f iladlr.f \
  dgetrf2.f dlaswp.f dlapy2.f dlaisnan.f dlarf.f dlarfg.f dgelq2.f dgeqr2.f disnan.f \
  dlarfb.f dlarft.f dlassq.f dorm2r.f dorml2.f dgelqf.f dgeqrf.f dlabad.f dlange.f \
  dlascl.f dlaset.f dormlq.f dormqr.f dtrtrs.f dgels.f dgetrf.f dgetrs.f dgesv.f
override LAPACK_SRC := $(if $(filter NETLIB,$(LAPACK_VER)),$(addprefix $(NETLIB_DIR)/,$(XNET_NETLIB_FILES)))
LAPACK_PROVIDER_ID = $(if $(LAPACK_SRC),sources:$(abspath $(LAPACK_SRC)),\
  $(if $(filter MKL,$(LAPACK_VER)),MKLROOT:$(abspath $(MKLROOT)),\
  $(if $(filter ATLAS,$(LAPACK_VER)),ATLAS_DIR:$(abspath $(ATLAS_DIR)),\
  $(if $(filter ESSL,$(LAPACK_VER)),ESSL_DIR:$(abspath $(ESSL_DIR)),\
  $(if $(filter ACCEL,$(LAPACK_VER)),system-framework:Accelerate,compiler-driver:$(COMPILE_FC))))))
SOLVER_PROVIDER_ID = $(if $(filter PARDISO,$(MATRIX_SOLVER)),MKLROOT:$(abspath $(MKLROOT)),\
  $(if $(SOLVER_SRC),sources:$(abspath $(SOLVER_SRC)),in-tree-jacobian:$(JAC_SRC)))
