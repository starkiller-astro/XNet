# Compiler and platform selections, defaults, and selector validation.
include $(XNET_DIR)/Makefile.opt
MPI_MODE ?= OFF
OPENMP_MODE ?= OFF
GPU_MODE ?= OFF
OPENACC_MODE ?= OFF
OPENMP_OL_MODE ?= OFF
GPU_BACKEND ?= CUDA
GPU_LAPACK_VER ?= CUBLAS
GPU_TARGET ?= Volta
XNET_EXTERNAL_MODULE_DIRS ?=
XNET_CPP ?= cpp
NVCC ?= nvcc
NVCCFLAGS ?= -O3 -m64
MISSING_EXTERNAL_MODULE_DIRS := $(filter-out $(wildcard $(XNET_EXTERNAL_MODULE_DIRS)),$(XNET_EXTERNAL_MODULE_DIRS))
ifneq ($(strip $(MISSING_EXTERNAL_MODULE_DIRS)),)
  $(error external module directory does not exist: $(MISSING_EXTERNAL_MODULE_DIRS))
endif
# Compatibility solver names select the same build configuration without
# recursive goal replay. Conflicting selectors fail while Make is still parsing.
SOLVER_GOALS := $(filter xnet_dense xnet_MA41 xnet_MA48 xnet_PARDISO,$(REQUESTED_GOALS))
ifneq ($(SOLVER_GOALS),)
  ifneq ($(words $(SOLVER_GOALS)),1)
    $(error select at most one solver compatibility target)
  endif
  IMPLIED_SOLVER := $(patsubst xnet_%,%,$(SOLVER_GOALS))
  ifeq ($(origin MATRIX_SOLVER),command line)
    ifneq ($(MATRIX_SOLVER),$(IMPLIED_SOLVER))
      $(error $(SOLVER_GOALS) conflicts with MATRIX_SOLVER=$(MATRIX_SOLVER))
    endif
  endif
  override MATRIX_SOLVER := $(IMPLIED_SOLVER)
endif

VALID_PE := GNU INTEL PGI NVIDIA NVHPC LLVM CCE CRAY XL
VALID_EOS := STARKILLER BAHCALL HELMHOLTZ
VALID_SOLVER := dense MA41 MA48 PARDISO
VALID_LAPACK := NETLIB MKL LIBSCI ATLAS ACCEL PGIBLAS ESSL
ifneq ($(filter $(PE_ENV),$(VALID_PE)),$(PE_ENV))
  $(error unsupported PE_ENV '$(PE_ENV)')
endif
ifneq ($(filter $(CMODE),OPT DEBUG),$(CMODE))
  $(error unsupported CMODE '$(CMODE)')
endif
ifneq ($(filter $(MPI_MODE),ON OFF),$(MPI_MODE))
  $(error MPI_MODE must be ON or OFF)
endif
ifneq ($(filter $(OPENMP_MODE),ON OFF),$(OPENMP_MODE))
  $(error OPENMP_MODE must be ON or OFF)
endif
ifneq ($(filter $(GPU_MODE),ON OFF),$(GPU_MODE))
  $(error GPU_MODE must be ON or OFF)
endif
ifneq ($(filter $(OPENACC_MODE),ON OFF),$(OPENACC_MODE))
  $(error OPENACC_MODE must be ON or OFF)
endif
ifneq ($(filter $(OPENMP_OL_MODE),ON OFF),$(OPENMP_OL_MODE))
  $(error OPENMP_OL_MODE must be ON or OFF)
endif
ifneq ($(filter $(EOS),$(VALID_EOS)),$(EOS))
  $(error unsupported EOS '$(EOS)')
endif
ifneq ($(filter $(MATRIX_SOLVER),$(VALID_SOLVER)),$(MATRIX_SOLVER))
  $(error unsupported MATRIX_SOLVER '$(MATRIX_SOLVER)')
endif

ifeq ($(GPU_MODE),OFF)
  ifeq ($(origin GPU_BACKEND),command line)
    ifneq ($(GPU_BACKEND),inactive)
      $(error GPU_BACKEND is active while GPU_MODE=OFF)
    endif
  endif
  ifeq ($(origin GPU_LAPACK_VER),command line)
    ifneq ($(GPU_LAPACK_VER),inactive)
      $(error GPU_LAPACK_VER is active while GPU_MODE=OFF)
    endif
  endif
  ifeq ($(origin GPU_TARGET),command line)
    ifneq ($(GPU_TARGET),inactive)
      $(error GPU_TARGET is active while GPU_MODE=OFF)
    endif
  endif
  ifneq ($(OPENACC_MODE)$(OPENMP_OL_MODE),OFFOFF)
    $(error accelerator directives require GPU_MODE=ON)
  endif
  override GPU_BACKEND := inactive
  override GPU_LAPACK_VER := inactive
  override GPU_TARGET := inactive
else
  ifeq ($(filter $(GPU_BACKEND),CUDA HIP),)
    $(error unsupported GPU_BACKEND '$(GPU_BACKEND)')
  endif
  ifneq ($(words $(filter ON,$(OPENACC_MODE) $(OPENMP_OL_MODE))),1)
    $(error GPU_MODE=ON requires exactly one accelerator directive mode)
  endif
  ifeq ($(GPU_BACKEND),CUDA)
    ifeq ($(filter $(GPU_LAPACK_VER),CUBLAS MAGMA),)
      $(error CUDA requires GPU_LAPACK_VER=CUBLAS or MAGMA)
    endif
    ifneq ($(filter Tesla,$(GPU_TARGET)),)
      override GPU_TARGET := $(filter-out Tesla,$(GPU_TARGET)) sm10 sm13
    endif
    ifneq ($(filter Fermi,$(GPU_TARGET)),)
      override GPU_TARGET := $(filter-out Fermi,$(GPU_TARGET)) sm20
    endif
    ifneq ($(filter Kepler,$(GPU_TARGET)),)
      override GPU_TARGET := $(filter-out Kepler,$(GPU_TARGET)) sm30 sm35
    endif
    ifneq ($(filter Maxwell,$(GPU_TARGET)),)
      override GPU_TARGET := $(filter-out Maxwell,$(GPU_TARGET)) sm50
    endif
    ifneq ($(filter Pascal,$(GPU_TARGET)),)
      override GPU_TARGET := $(filter-out Pascal,$(GPU_TARGET)) sm60
    endif
    ifneq ($(filter Volta,$(GPU_TARGET)),)
      override GPU_TARGET := $(filter-out Volta,$(GPU_TARGET)) sm70
    endif
    ifneq ($(filter Ampere,$(GPU_TARGET)),)
      override GPU_TARGET := $(filter-out Ampere,$(GPU_TARGET)) sm80
    endif
    override GPU_TARGET := $(sort $(GPU_TARGET))
    ifneq ($(filter-out sm10 sm13 sm20 sm30 sm35 sm50 sm60 sm70 sm80,$(GPU_TARGET)),)
      $(error unsupported CUDA GPU_TARGET '$(GPU_TARGET)')
    endif
    override NVCCFLAGS += $(foreach target,$(GPU_TARGET),-gencode arch=compute_$(patsubst sm%,%,$(target)),code=sm_$(patsubst sm%,%,$(target)))
  else
    ifneq ($(GPU_LAPACK_VER),ROCM)
      $(error HIP requires GPU_LAPACK_VER=ROCM)
    endif
    ifeq ($(origin GPU_TARGET),command line)
      ifneq ($(GPU_TARGET),inactive)
        $(error GPU_TARGET is valid only with GPU_BACKEND=CUDA)
      endif
    endif
    override GPU_TARGET := inactive
  endif
  ifneq ($(MATRIX_SOLVER),dense)
    $(error GPU builds require MATRIX_SOLVER=dense)
  endif
endif
ifeq ($(EOS),HELMHOLTZ)
  ifndef HELMHOLTZ_PATH
    $(error HELMHOLTZ_PATH must be set for EOS=HELMHOLTZ)
  endif
endif

include $(XNET_DIR)/Makefile.internal

PRIVATE_PROVIDER_VARIABLES := MPI_SRC EOS_SRC JAC_SRC GPU_PROVIDER_SRC LAPACK_SRC SOLVER_SRC \
  CORE_SRC DRIVER_SRC XNSE_MAIN_SRC PROBE_SRC SPARSE_SRC XNSE_CORE_SRC SETUP_CORE_SRC
PROVIDER_OVERRIDES := $(strip $(foreach variable,$(PRIVATE_PROVIDER_VARIABLES),\
  $(if $(filter command line,$(origin $(variable))),$(variable))))
ifneq ($(PROVIDER_OVERRIDES),)
  $(error provider variables are internal and may not be overridden: $(PROVIDER_OVERRIDES))
endif

ifneq ($(filter $(LAPACK_VER),$(VALID_LAPACK)),$(LAPACK_VER))
  $(error unsupported or unwired LAPACK_VER '$(LAPACK_VER)')
endif
ifeq ($(MATRIX_SOLVER),PARDISO)
  ifneq ($(LAPACK_VER),MKL)
    $(error unsupported standalone PARDISO: use LAPACK_VER=MKL)
  endif
endif
ifeq ($(LAPACK_VER),MKL)
  ifeq ($(strip $(LAPACK_LIBS)),)
    $(error LAPACK_VER=MKL requires resolved MKL_LIBS or a working MKL link tool)
  endif
endif
