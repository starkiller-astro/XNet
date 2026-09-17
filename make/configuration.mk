# Public build-choice validation. Set these variables on the Make command
# line. make/build.mk reads the tracked defaults from Makefile.opt before
# choosing BUILD_NAME:
#   PE_ENV: GNU, INTEL, PGI, NVIDIA, NVHPC, LLVM, CCE, CRAY, XL
#   CMODE: OPT, DEBUG
#   MPI_MODE, OPENMP_MODE, GPU_MODE: ON, OFF
#   GPU_BACKEND: CUDA, HIP
#   GPU_LAPACK_VER: CUBLAS, MAGMA, ROCM
#   OPENACC_MODE, OPENMP_OL_MODE: ON, OFF
#   EOS: STARKILLER, BAHCALL, HELMHOLTZ
#   MATRIX_SOLVER: dense, MA41, MA48, PARDISO_MKL
#   LAPACK_VER: NETLIB, MKL, LIBSCI, ATLAS, ACCEL, PGIBLAS, ESSL

# XNet sources use C-preprocessor macros before Fortran compilation. This may
# be set to another conventional C preprocessor when required by a toolchain.
XNET_CPP ?= cpp

VALID_PE := GNU INTEL PGI NVIDIA NVHPC LLVM CCE CRAY XL
VALID_EOS := STARKILLER BAHCALL HELMHOLTZ
VALID_SOLVER := dense MA41 MA48 PARDISO_MKL
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

# Accelerator libraries are inactive unless GPU_MODE=ON. A GPU build selects
# exactly one backend and one directive model.
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
  ifneq ($(OPENACC_MODE)$(OPENMP_OL_MODE),OFFOFF)
    $(error accelerator directives require GPU_MODE=ON)
  endif
  GPU_BACKEND := inactive
  GPU_LAPACK_VER := inactive
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
  else ifeq ($(GPU_BACKEND),HIP)
    ifneq ($(GPU_LAPACK_VER),ROCM)
      $(error HIP requires GPU_LAPACK_VER=ROCM)
    endif
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

# The maintained PARDISO implementation is the oneMKL interface. Selecting it
# chooses MKL unless the caller supplied an incompatible LAPACK_VER, which is
# rejected below.
ifeq ($(MATRIX_SOLVER),PARDISO_MKL)
  LAPACK_VER := MKL
endif

include Makefile.internal

ifneq ($(filter $(LAPACK_VER),$(VALID_LAPACK)),$(LAPACK_VER))
  $(error unsupported or unwired LAPACK_VER '$(LAPACK_VER)')
endif
ifeq ($(MATRIX_SOLVER),PARDISO_MKL)
  ifneq ($(LAPACK_VER),MKL)
    $(error MATRIX_SOLVER=PARDISO_MKL requires LAPACK_VER=MKL)
  endif
endif
ifeq ($(LAPACK_VER),MKL)
  ifeq ($(strip $(MKL_LIBS)),)
    $(error LAPACK_VER=MKL requires resolved MKL_LIBS or a working MKL link tool)
  endif
endif
