# Production source lists, object paths, and BUILD_DIR configuration checking.
CORE_SRC := $(addprefix $(SOURCE_DIR)/,xnet_controls.F90 xnet_data.F90 xnet_output.F90 \
  xnet_abundances.F90 xnet_conditions.F90 xnet_constants.F90 xnet_evolve.F90 xnet_fd.F90 \
  xnet_ffn.F90 xnet_flux.F90 xnet_gpu.F90 xnet_integrate.F90 xnet_integrate_bdf.F90 \
  xnet_integrate_be.F90 xnet_linalg.F90 xnet_match.F90 xnet_nnu.F90 xnet_nse.F90 \
  xnet_preprocess.F90 xnet_screening.F90 xnet_timers.F90 xnet_types.F90 xnet_util.F90)
DRIVER_SRC := $(addprefix $(SOURCE_DIR)/,model_input_ascii.F90 net.F90)
XNSE_MAIN_SRC := $(SOURCE_DIR)/nse_slice.F90
XNSE_CORE_SRC := $(addprefix $(SOURCE_DIR)/,xnet_abundances.F90 xnet_conditions.F90 xnet_controls.F90 \
  xnet_data.F90 xnet_constants.F90 xnet_fd.F90 xnet_ffn.F90 xnet_match.F90 xnet_nse.F90 \
  xnet_preprocess.F90 xnet_util.F90 xnet_nnu.F90 xnet_timers.F90 xnet_types.F90)
SETUP_CORE_SRC := $(addprefix $(SOURCE_DIR)/,xnet_conditions.F90 xnet_controls.F90 xnet_data.F90 \
  net_setup.F90 xnet_constants.F90 xnet_fd.F90 xnet_ffn.F90 xnet_nnu.F90 \
  xnet_preprocess.F90 xnet_util.F90 xnet_types.F90)

# XNet, EOS, Jacobian, and accelerator bindings are free-form Fortran. The
# optional MA41/MA48 and bundled NETLIB sources are fixed-form Fortran.
SOURCE_FREE_SRC := $(sort $(CORE_SRC) $(DRIVER_SRC) $(XNSE_CORE_SRC) \
                   $(SETUP_CORE_SRC) $(XNSE_MAIN_SRC) $(MPI_SRC))
EOS_FREE_SRC := $(sort $(EOS_SRC))
SOLVER_FREE_SRC := $(JAC_SRC)
GPU_FREE_SRC := $(sort $(GPU_PROVIDER_SRC))
SOLVER_FIXED_SRC := $(sort $(filter %.F %.f,$(SOLVER_SRC)))
LAPACK_FIXED_SRC := $(sort $(filter %.F %.f,$(LAPACK_SRC)))

ALL_SELECTED_SRC := $(SOURCE_FREE_SRC) $(EOS_FREE_SRC) $(SOLVER_FREE_SRC) \
                    $(GPU_FREE_SRC) $(SOLVER_FIXED_SRC) $(LAPACK_FIXED_SRC)
MISSING_SELECTED_SRC := $(filter-out $(wildcard $(ALL_SELECTED_SRC)),$(ALL_SELECTED_SRC))
ifneq ($(strip $(MISSING_SELECTED_SRC)),)
  $(error selected source does not exist: $(MISSING_SELECTED_SRC))
endif

# XNet's variadic accelerator macros require a conventional C preprocessor;
# native Fortran preprocessors are not sufficient for every supported compiler.
# The retained .f90 files also give every compiler the same macro expansion.
strip_cpp_flags = $(filter-out -cpp -fpp -Mpreprocess -eZ -qpreprocess -D% -U%,$(1))
COMPILE_FC := $(FC)
COMPILE_FFLAGS := $(call strip_cpp_flags,$(FFLAGS))
COMPILE_F90FLAGS := $(call strip_cpp_flags,$(F90FLAGS))

ifeq ($(filter INTEL PGI NVIDIA NVHPC,$(PE_ENV)),)
  ifeq ($(PE_ENV),XL)
    MODULE_FLAGS := -qmoddir=$(MOD_DIR) -I$(MOD_DIR)
  else
    MODULE_FLAGS := -J$(MOD_DIR) -I$(MOD_DIR)
  endif
else
  MODULE_FLAGS := -module $(MOD_DIR) -I$(MOD_DIR)
endif

# External cpp must be able to resolve mpif.h before the retained source is
# compiled from PP_DIR.
ifeq ($(MPI_MODE),ON)
  MPI_WRAPPER_SHOW := $(shell $(COMPILE_FC) --showme:compile 2>/dev/null || $(COMPILE_FC) -show 2>/dev/null || $(COMPILE_FC) --cray-print-opts=all 2>/dev/null || true)
  MPI_HEADER_DIR := $(shell header=`$(COMPILE_FC) -print-file-name=mpif.h 2>/dev/null`; test -f "$$header" && dirname "$$header" || true)
  MPI_INCLUDE_FLAGS := $(filter -I%,$(MPI_WRAPPER_SHOW)) $(if $(MPI_HEADER_DIR),-I$(MPI_HEADER_DIR))
endif

CPP_INCLUDE_FLAGS := -I$(SOURCE_DIR) $(sort $(addprefix -I,$(dir $(SOURCE_FREE_SRC) \
  $(EOS_FREE_SRC) $(SOLVER_FREE_SRC) $(GPU_FREE_SRC)))) \
  $(filter -I%,$(FFLAGS) $(F90FLAGS) $(LAPACK_INC) $(SOLVER_INC)) $(MPI_INCLUDE_FLAGS)
CPP_DEFINITION_FLAGS := $(filter -D% -U%,$(FFLAGS) $(F90FLAGS)) $(GPU_DEFINES)
CPP_EFFECTIVE_FLAGS := -P -C -nostdinc $(CPP_INCLUDE_FLAGS) $(CPP_DEFINITION_FLAGS)
ifneq ($(filter CCE CRAY,$(PE_ENV)),)
  ifeq ($(OPENMP_OL_MODE),ON)
    CRAY_OMP_PREPROCESS := yes
  endif
endif
ifeq ($(CRAY_OMP_PREPROCESS),)
  CRAY_OMP_PREPROCESS := no
endif

ifeq ($(CRAY_OMP_PREPROCESS),yes)
  CRAY_PREPROCESS_INPUT := $(ROOT_DIR)/make/crayftn_cpp.sh
else
  CRAY_PREPROCESS_INPUT :=
endif

SOURCE_OBJ_DIR := $(OBJ_DIR)/source
EOS_OBJ_DIR := $(OBJ_DIR)/eos
SOLVER_OBJ_DIR := $(OBJ_DIR)/solver
LAPACK_OBJ_DIR := $(OBJ_DIR)/lapack
GPU_OBJ_DIR := $(OBJ_DIR)/gpu
SOURCE_PP_DIR := $(PP_DIR)/source
EOS_PP_DIR := $(PP_DIR)/eos
SOLVER_PP_DIR := $(PP_DIR)/solver
GPU_PP_DIR := $(PP_DIR)/gpu

stem = $(basename $(notdir $(1)))
source_obj = $(SOURCE_OBJ_DIR)/$(call stem,$(1)).o
eos_obj = $(EOS_OBJ_DIR)/$(call stem,$(1)).o
solver_obj = $(SOLVER_OBJ_DIR)/$(call stem,$(1)).o
lapack_obj = $(LAPACK_OBJ_DIR)/$(call stem,$(1)).o
gpu_obj = $(GPU_OBJ_DIR)/$(call stem,$(1)).o
source_pp = $(SOURCE_PP_DIR)/$(call stem,$(1)).f90
eos_pp = $(EOS_PP_DIR)/$(call stem,$(1)).f90
solver_pp = $(SOLVER_PP_DIR)/$(call stem,$(1)).f90
gpu_pp = $(GPU_PP_DIR)/$(call stem,$(1)).f90
parent_dir = $(patsubst %/,%,$(dir $(1)))

EOS_OBJ := $(foreach source,$(EOS_FREE_SRC),$(call eos_obj,$(source)))
SOLVER_FREE_OBJ := $(foreach source,$(SOLVER_FREE_SRC),$(call solver_obj,$(source)))
GPU_OBJ := $(foreach source,$(GPU_FREE_SRC),$(call gpu_obj,$(source)))
SOLVER_FIXED_OBJ := $(foreach source,$(SOLVER_FIXED_SRC),$(call solver_obj,$(source)))
LAPACK_OBJ := $(foreach source,$(LAPACK_FIXED_SRC),$(call lapack_obj,$(source)))
SOLVER_OBJ := $(SOLVER_FREE_OBJ) $(SOLVER_FIXED_OBJ)

XNET_SOURCE_OBJ := $(foreach source,$(sort $(CORE_SRC) $(DRIVER_SRC) $(MPI_SRC)),$(call source_obj,$(source)))
XNSE_SOURCE_OBJ := $(foreach source,$(sort $(XNSE_CORE_SRC) $(XNSE_MAIN_SRC) $(MPI_SRC)),$(call source_obj,$(source)))
SETUP_SOURCE_OBJ := $(foreach source,$(sort $(SETUP_CORE_SRC) $(MPI_SRC)),$(call source_obj,$(source)))
XNET_OBJ := $(XNET_SOURCE_OBJ) $(EOS_OBJ) $(SOLVER_OBJ) $(LAPACK_OBJ) $(GPU_OBJ)
XNSE_OBJ := $(XNSE_SOURCE_OBJ) $(EOS_OBJ) $(LAPACK_OBJ)
SETUP_OBJ := $(SETUP_SOURCE_OBJ)
XNET_EXE := $(BIN_DIR)/xnet
XNSE_EXE := $(BIN_DIR)/xnse
NET_SETUP_EXE := $(BIN_DIR)/net_setup

# Record only paths relevant to the selected configuration. The complete
# source and flag records below provide practical stale-object protection.
CONFIG_NETLIB_DIR := $(if $(filter NETLIB,$(LAPACK_VER)),$(abspath $(NETLIB_DIR)))
CONFIG_MKLROOT := $(if $(filter MKL,$(LAPACK_VER)),$(abspath $(MKLROOT)))
CONFIG_MA41_DIR := $(if $(filter MA41,$(MATRIX_SOLVER)),$(abspath $(MA41_DIR)))
CONFIG_MA48_DIR := $(if $(filter MA48,$(MATRIX_SOLVER)),$(abspath $(MA48_DIR)))
CONFIG_CUDA_DIR := $(if $(filter CUDA,$(GPU_BACKEND)),$(abspath $(CUDA_DIR)))
CONFIG_MAGMA_DIR := $(if $(filter MAGMA,$(GPU_LAPACK_VER)),$(abspath $(MAGMA_DIR)))
CONFIG_HIPFORT_DIR := $(if $(filter HIP,$(GPU_BACKEND)),$(abspath $(HIPFORT_DIR)))
CONFIG_ROCM_DIR := $(if $(filter HIP,$(GPU_BACKEND)),$(abspath $(ROCM_DIR)))
CONFIG_HELMHOLTZ_PATH := $(if $(filter HELMHOLTZ,$(EOS)),$(abspath $(HELMHOLTZ_PATH)))

define XNET_CONFIG_NEWLINE


endef
XNET_CONFIG_CR := $(shell printf '\r')
CONFIG_RECORD_VALUES := $(ROOT_DIR) $(PE_ENV) $(CMODE) $(MACHINE) $(MPI_MODE) $(OPENMP_MODE) \
  $(GPU_MODE) $(GPU_BACKEND) $(GPU_LAPACK_VER) $(OPENACC_MODE) $(OPENMP_OL_MODE) \
  $(EOS) $(MATRIX_SOLVER) $(LAPACK_VER) $(COMPILE_FC) $(XNET_CPP) $(LDR) \
  $(COMPILE_FFLAGS) $(COMPILE_F90FLAGS) $(F77FLAGS) $(CPP_EFFECTIVE_FLAGS) $(LDFLAGS) \
  $(LAPACK_INC) $(LAPACK_LIBDIR) $(LAPACK_LIBS) $(SOLVER_INC) $(SOLVER_LIBDIR) \
  $(SOLVER_LIBS) $(MPI_SRC) $(EOS_SRC) $(JAC_SRC) $(SOLVER_SRC) $(LAPACK_SRC) \
  $(GPU_PROVIDER_SRC) $(GPU_MODULE_FLAGS) $(CONFIG_NETLIB_DIR) $(CONFIG_MKLROOT) \
  $(CONFIG_MA41_DIR) $(CONFIG_MA48_DIR) $(CONFIG_CUDA_DIR) $(CONFIG_MAGMA_DIR) \
  $(CONFIG_HIPFORT_DIR) $(CONFIG_ROCM_DIR) $(CONFIG_HELMHOLTZ_PATH)
CONFIG_UNSUPPORTED_QUOTE := $(findstring ',$(CONFIG_RECORD_VALUES))
CONFIG_UNSUPPORTED_LF := $(if $(findstring $(XNET_CONFIG_NEWLINE),$(CONFIG_RECORD_VALUES)),yes)
CONFIG_UNSUPPORTED_CR := $(if $(findstring $(XNET_CONFIG_CR),$(CONFIG_RECORD_VALUES)),yes)
ifneq ($(strip $(CONFIG_UNSUPPORTED_QUOTE)$(CONFIG_UNSUPPORTED_LF)$(CONFIG_UNSUPPORTED_CR)),)
  $(error config.txt cannot represent single quotes or line breaks in recorded settings; choose a different spelling or BUILD_DIR)
endif

CONFIG_RECORD_ARGUMENTS = \
  'SOURCE_ROOT=$(ROOT_DIR)' \
  'PE_ENV=$(PE_ENV)' 'CMODE=$(CMODE)' 'MACHINE=$(MACHINE)' \
  'MPI_MODE=$(MPI_MODE)' 'OPENMP_MODE=$(OPENMP_MODE)' \
  'GPU_MODE=$(GPU_MODE)' 'GPU_BACKEND=$(GPU_BACKEND)' \
  'GPU_LAPACK_VER=$(GPU_LAPACK_VER)' \
  'OPENACC_MODE=$(OPENACC_MODE)' 'OPENMP_OL_MODE=$(OPENMP_OL_MODE)' \
  'EOS=$(EOS)' 'MATRIX_SOLVER=$(MATRIX_SOLVER)' 'LAPACK_VER=$(LAPACK_VER)' \
  'FC=$(COMPILE_FC)' 'CPP=$(XNET_CPP)' 'LDR=$(LDR)' \
  'FFLAGS=$(COMPILE_FFLAGS)' 'F90FLAGS=$(COMPILE_F90FLAGS)' 'F77FLAGS=$(F77FLAGS)' \
  'CPPFLAGS=$(CPP_EFFECTIVE_FLAGS)' 'LDFLAGS=$(LDFLAGS)' \
  'LAPACK_INC=$(LAPACK_INC)' 'LAPACK_LIBDIR=$(LAPACK_LIBDIR)' 'LAPACK_LIBS=$(LAPACK_LIBS)' \
  'SOLVER_INC=$(SOLVER_INC)' 'SOLVER_LIBDIR=$(SOLVER_LIBDIR)' 'SOLVER_LIBS=$(SOLVER_LIBS)' \
  'PARALLEL_SOURCE=$(MPI_SRC)' 'EOS_SOURCES=$(EOS_SRC)' 'JACOBIAN_SOURCE=$(JAC_SRC)' \
  'SOLVER_SOURCES=$(SOLVER_SRC)' 'LAPACK_SOURCES=$(LAPACK_SRC)' \
  'GPU_SOURCES=$(GPU_PROVIDER_SRC)' 'GPU_MODULE_FLAGS=$(GPU_MODULE_FLAGS)' \
  'NETLIB_DIR=$(CONFIG_NETLIB_DIR)' 'MKLROOT=$(CONFIG_MKLROOT)' \
  'MA41_DIR=$(CONFIG_MA41_DIR)' 'MA48_DIR=$(CONFIG_MA48_DIR)' \
  'CUDA_DIR=$(CONFIG_CUDA_DIR)' 'MAGMA_DIR=$(CONFIG_MAGMA_DIR)' \
  'HIPFORT_DIR=$(CONFIG_HIPFORT_DIR)' 'ROCM_DIR=$(CONFIG_ROCM_DIR)' \
  'HELMHOLTZ_PATH=$(CONFIG_HELMHOLTZ_PATH)'

.PHONY: xnet xnse net_setup all xinab xnet_gpu FORCE
FORCE:
$(CONFIG): FORCE
	@$(ROOT_DIR)/make/update-config.sh '$@' $(CONFIG_RECORD_ARGUMENTS)

$(SOURCE_OBJ_DIR) $(EOS_OBJ_DIR) $(SOLVER_OBJ_DIR) $(LAPACK_OBJ_DIR) $(GPU_OBJ_DIR) \
$(SOURCE_PP_DIR) $(EOS_PP_DIR) $(SOLVER_PP_DIR) $(GPU_PP_DIR) $(MOD_DIR) $(BIN_DIR): | $(CONFIG)
	@mkdir -p '$@'

BUILD_LOGIC_INPUTS := $(ROOT_DIR)/Makefile $(ROOT_DIR)/Makefile.opt \
  $(ROOT_DIR)/Makefile.internal $(ROOT_DIR)/make/build.mk \
  $(ROOT_DIR)/make/configuration.mk $(ROOT_DIR)/make/providers.mk \
  $(ROOT_DIR)/make/sources.mk $(ROOT_DIR)/make/dependencies.mk \
  $(ROOT_DIR)/make/rules.mk $(ROOT_DIR)/make/update-config.sh \
  $(ROOT_DIR)/make/machines.mk $(ROOT_DIR)/make/machines/generic.mk \
  $(ROOT_DIR)/make/machines/cray-pe.mk
PP_STATIC_INPUTS := $(SOURCE_DIR)/xnet_macros.fh $(BUILD_LOGIC_INPUTS) $(CRAY_PREPROCESS_INPUT)
