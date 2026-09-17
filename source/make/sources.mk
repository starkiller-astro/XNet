# Source inventories, object mapping, artifact paths, and configuration reuse.
CORE_SRC := $(addprefix $(XNET_DIR)/,xnet_controls.F90 xnet_data.F90 xnet_output.F90 \
  xnet_abundances.F90 xnet_conditions.F90 xnet_constants.F90 xnet_evolve.F90 xnet_fd.F90 \
  xnet_ffn.F90 xnet_flux.F90 xnet_gpu.F90 xnet_integrate.F90 xnet_integrate_bdf.F90 \
  xnet_integrate_be.F90 xnet_linalg.F90 xnet_match.F90 xnet_nnu.F90 xnet_nse.F90 \
  xnet_preprocess.F90 xnet_screening.F90 xnet_timers.F90 xnet_types.F90 xnet_util.F90)
DRIVER_SRC := $(addprefix $(XNET_DIR)/,model_input_ascii.F90 net.F90)
XNSE_MAIN_SRC := $(XNET_DIR)/nse_slice.F90
PROBE_SRC := $(ROOT_DIR)/test/qualification/frontier/gpu_linalg_probe.F90
SPARSE_SRC := $(if $(filter MA48 PARDISO,$(MATRIX_SOLVER)),$(XNET_DIR)/xnet_sparse.F90)
XNSE_CORE_SRC := $(addprefix $(XNET_DIR)/,xnet_abundances.F90 xnet_conditions.F90 xnet_controls.F90 \
  xnet_data.F90 xnet_constants.F90 xnet_fd.F90 xnet_ffn.F90 xnet_match.F90 xnet_nse.F90 \
  xnet_preprocess.F90 xnet_util.F90 xnet_nnu.F90 xnet_timers.F90 xnet_types.F90)
SETUP_CORE_SRC := $(addprefix $(XNET_DIR)/,xnet_conditions.F90 xnet_controls.F90 xnet_data.F90 \
  net_setup.F90 xnet_constants.F90 xnet_fd.F90 xnet_ffn.F90 xnet_nnu.F90 \
  xnet_preprocess.F90 xnet_util.F90 xnet_types.F90)
SOURCE_FREE_SRC := $(sort $(CORE_SRC) $(DRIVER_SRC) $(XNSE_CORE_SRC) $(SETUP_CORE_SRC) \
                   $(XNSE_MAIN_SRC) $(MPI_SRC) $(PROBE_SRC))
EOS_FREE_SRC := $(sort $(EOS_SRC))
SOLVER_FREE_SRC := $(sort $(JAC_SRC) $(SPARSE_SRC) $(filter %.F90 %.f90,$(SOLVER_SRC)))
LAPACK_FREE_SRC := $(sort $(filter %.F90 %.f90,$(LAPACK_SRC)))
GPU_FREE_SRC := $(sort $(GPU_PROVIDER_SRC))
SOLVER_FIXED_SRC := $(sort $(filter %.F %.f,$(SOLVER_SRC)))
LAPACK_FIXED_SRC := $(sort $(filter %.F %.f,$(LAPACK_SRC)))
SOLVER_C_SRC := $(sort $(filter %.c,$(SOLVER_SRC)))
LAPACK_C_SRC := $(sort $(filter %.c,$(LAPACK_SRC)))
SOLVER_CXX_SRC := $(sort $(filter %.cpp %.cxx,$(SOLVER_SRC)))
LAPACK_CXX_SRC := $(sort $(filter %.cpp %.cxx,$(LAPACK_SRC)))
SOLVER_CUDA_SRC := $(sort $(filter %.cu,$(SOLVER_SRC)))
LAPACK_CUDA_SRC := $(sort $(filter %.cu,$(LAPACK_SRC)))
ALL_SELECTED_SRC := $(SOURCE_FREE_SRC) $(EOS_FREE_SRC) $(SOLVER_FREE_SRC) $(LAPACK_FREE_SRC) \
                    $(GPU_FREE_SRC) $(SOLVER_FIXED_SRC) $(LAPACK_FIXED_SRC) \
                    $(SOLVER_C_SRC) $(LAPACK_C_SRC) $(SOLVER_CXX_SRC) $(LAPACK_CXX_SRC) \
                    $(SOLVER_CUDA_SRC) $(LAPACK_CUDA_SRC)
MISSING_SELECTED_SRC := $(filter-out $(wildcard $(ALL_SELECTED_SRC)),$(ALL_SELECTED_SRC))
ifneq ($(strip $(MISSING_SELECTED_SRC)),)
  $(error selected provider source does not exist: $(MISSING_SELECTED_SRC))
endif

FRONTIER_REQUESTED := $(filter frontier_gpu_linalg_probe,$(REQUESTED_GOALS))
ifneq ($(FRONTIER_REQUESTED),)
  ifneq ($(PE_ENV):$(MACHINE):$(GPU_MODE):$(GPU_BACKEND):$(GPU_LAPACK_VER):$(OPENMP_OL_MODE):$(OPENACC_MODE),CRAY:frontier:ON:HIP:ROCM:ON:OFF)
    $(error frontier_gpu_linalg_probe requires exact Frontier CRAY/HIP/ROCM/OpenMP-offload selectors)
  endif
  ifneq ($(CMODE):$(MPI_MODE):$(OPENMP_MODE):$(EOS):$(MATRIX_SOLVER):$(LAPACK_VER),OPT:OFF:OFF:STARKILLER:dense:LIBSCI)
    $(error frontier_gpu_linalg_probe is gated to the exact retained Frontier selector row)
  endif
endif

strip_cpp_flags = $(filter-out -cpp -fpp -Mpreprocess -eZ -qpreprocess -D% -U%,$(1))
COMPILE_FC := $(if $(and $(filter CCE CRAY,$(PE_ENV)),$(XNET_CRAY_FTN)),$(XNET_CRAY_FTN),$(FC))
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
ifeq ($(MPI_MODE),ON)
  MPI_WRAPPER_SHOW := $(shell $(COMPILE_FC) --showme:compile 2>/dev/null || $(COMPILE_FC) -show 2>/dev/null || $(COMPILE_FC) --cray-print-opts=all 2>/dev/null || true)
  MPI_HEADER_DIR := $(shell header=`$(COMPILE_FC) -print-file-name=mpif.h 2>/dev/null`; test -f "$$header" && dirname "$$header" || true)
  MPI_INCLUDE_FLAGS := $(filter -I%,$(MPI_WRAPPER_SHOW)) $(if $(MPI_HEADER_DIR),-I$(MPI_HEADER_DIR))
endif
EXTERNAL_MODULE_FLAGS := $(addprefix -I,$(SELECTED_EXTERNAL_MODULE_DIRS))
CPP_INCLUDE_FLAGS := -I$(XNET_DIR) $(sort $(addprefix -I,$(dir $(SOURCE_FREE_SRC) $(EOS_FREE_SRC) \
  $(SOLVER_FREE_SRC) $(LAPACK_FREE_SRC) $(GPU_FREE_SRC)))) \
  $(filter -I%,$(FFLAGS) $(F90FLAGS) $(LAPACK_INC) $(SOLVER_INC)) $(MPI_INCLUDE_FLAGS)
CPP_DEFINITION_FLAGS := $(filter -D% -U%,$(FFLAGS) $(F90FLAGS)) $(GPU_DEFINES)
CPP_EFFECTIVE_FLAGS := -P -C -nostdinc $(CPP_INCLUDE_FLAGS) $(CPP_DEFINITION_FLAGS)
ifneq ($(filter CCE CRAY,$(PE_ENV)),)
  CRAY_RETAINED_PP := yes
else
  CRAY_RETAINED_PP := no
endif

SOURCE_OBJ_DIR := $(OBJ_DIR)/source
EOS_OBJ_DIR := $(OBJ_DIR)/eos
SOLVER_OBJ_DIR := $(OBJ_DIR)/solver
LAPACK_OBJ_DIR := $(OBJ_DIR)/lapack
GPU_OBJ_DIR := $(OBJ_DIR)/gpu
SOURCE_PP_DIR := $(PP_DIR)/source
EOS_PP_DIR := $(PP_DIR)/eos
SOLVER_PP_DIR := $(PP_DIR)/solver
LAPACK_PP_DIR := $(PP_DIR)/lapack
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
lapack_pp = $(LAPACK_PP_DIR)/$(call stem,$(1)).f90
gpu_pp = $(GPU_PP_DIR)/$(call stem,$(1)).f90
parent_dir = $(patsubst %/,%,$(dir $(1)))

EOS_OBJ := $(foreach source,$(EOS_FREE_SRC),$(call eos_obj,$(source)))
SOLVER_FREE_OBJ := $(foreach source,$(SOLVER_FREE_SRC),$(call solver_obj,$(source)))
LAPACK_FREE_OBJ := $(foreach source,$(LAPACK_FREE_SRC),$(call lapack_obj,$(source)))
GPU_OBJ := $(foreach source,$(GPU_FREE_SRC),$(call gpu_obj,$(source)))
SOLVER_FIXED_OBJ := $(foreach source,$(SOLVER_FIXED_SRC),$(call solver_obj,$(source)))
LAPACK_FIXED_OBJ := $(foreach source,$(LAPACK_FIXED_SRC),$(call lapack_obj,$(source)))
SOLVER_C_OBJ := $(foreach source,$(SOLVER_C_SRC),$(call solver_obj,$(source)))
LAPACK_C_OBJ := $(foreach source,$(LAPACK_C_SRC),$(call lapack_obj,$(source)))
SOLVER_CXX_OBJ := $(foreach source,$(SOLVER_CXX_SRC),$(call solver_obj,$(source)))
LAPACK_CXX_OBJ := $(foreach source,$(LAPACK_CXX_SRC),$(call lapack_obj,$(source)))
SOLVER_CUDA_OBJ := $(foreach source,$(SOLVER_CUDA_SRC),$(call solver_obj,$(source)))
LAPACK_CUDA_OBJ := $(foreach source,$(LAPACK_CUDA_SRC),$(call lapack_obj,$(source)))
LAPACK_OBJ_ALL := $(LAPACK_FREE_OBJ) $(LAPACK_FIXED_OBJ) $(LAPACK_C_OBJ) $(LAPACK_CXX_OBJ) $(LAPACK_CUDA_OBJ)
SOLVER_OBJ_ALL := $(SOLVER_FREE_OBJ) $(SOLVER_FIXED_OBJ) $(SOLVER_C_OBJ) $(SOLVER_CXX_OBJ) $(SOLVER_CUDA_OBJ)
XNET_SOURCE_OBJ := $(foreach source,$(sort $(CORE_SRC) $(DRIVER_SRC) $(MPI_SRC)),$(call source_obj,$(source)))
XNSE_SOURCE_OBJ := $(foreach source,$(sort $(XNSE_CORE_SRC) $(XNSE_MAIN_SRC) $(MPI_SRC)),$(call source_obj,$(source)))
SETUP_SOURCE_OBJ := $(foreach source,$(sort $(SETUP_CORE_SRC) $(MPI_SRC)),$(call source_obj,$(source)))
PROBE_SOURCE_OBJ := $(foreach source,$(sort $(CORE_SRC) $(MPI_SRC) $(PROBE_SRC)),$(call source_obj,$(source)))
XNET_OBJ := $(XNET_SOURCE_OBJ) $(EOS_OBJ) $(SOLVER_OBJ_ALL) $(LAPACK_OBJ_ALL) $(GPU_OBJ)
XNSE_OBJ := $(XNSE_SOURCE_OBJ) $(EOS_OBJ) $(LAPACK_OBJ_ALL)
SETUP_OBJ := $(SETUP_SOURCE_OBJ)
PROBE_OBJ := $(PROBE_SOURCE_OBJ) $(EOS_OBJ) $(SOLVER_OBJ_ALL) $(LAPACK_OBJ_ALL) $(GPU_OBJ)
XNET_EXE := $(BIN_DIR)/xnet
XNSE_EXE := $(BIN_DIR)/xnse
NET_SETUP_EXE := $(BIN_DIR)/net_setup
PROBE_EXE := $(BIN_DIR)/frontier_gpu_linalg_probe

XNET_CFG_SOURCE_ROOT := $(ROOT_DIR)
XNET_CFG_PE_ENV := $(PE_ENV)
XNET_CFG_CMODE := $(CMODE)
XNET_CFG_MACHINE := $(MACHINE)
XNET_CFG_MPI_MODE := $(MPI_MODE)
XNET_CFG_OPENMP_MODE := $(OPENMP_MODE)
XNET_CFG_GPU_MODE := $(GPU_MODE)
XNET_CFG_GPU_BACKEND := $(GPU_BACKEND)
XNET_CFG_GPU_LAPACK_VER := $(GPU_LAPACK_VER)
XNET_CFG_GPU_TARGET := $(GPU_TARGET)
XNET_CFG_OPENACC_MODE := $(OPENACC_MODE)
XNET_CFG_OPENMP_OL_MODE := $(OPENMP_OL_MODE)
XNET_CFG_EOS := $(EOS)
XNET_CFG_MATRIX_SOLVER := $(MATRIX_SOLVER)
XNET_CFG_LAPACK_VER := $(LAPACK_VER)
define XNET_CONFIG_NEWLINE


endef
XNET_CONFIG_CR := $(shell printf '\r')
CONFIG_RECORD_VALUES := $(XNET_CFG_SOURCE_ROOT) $(XNET_CFG_PE_ENV) $(XNET_CFG_CMODE) $(XNET_CFG_MACHINE) \
  $(XNET_CFG_MPI_MODE) $(XNET_CFG_OPENMP_MODE) $(XNET_CFG_GPU_MODE) $(XNET_CFG_GPU_BACKEND) \
  $(XNET_CFG_GPU_LAPACK_VER) $(XNET_CFG_GPU_TARGET) $(XNET_CFG_OPENACC_MODE) $(XNET_CFG_OPENMP_OL_MODE) \
  $(XNET_CFG_EOS) $(XNET_CFG_MATRIX_SOLVER) $(XNET_CFG_LAPACK_VER) $(COMPILE_FC) $(CC) $(CXX) $(NVCC) $(XNET_CPP) $(LDR) \
  $(COMPILE_FFLAGS) $(COMPILE_F90FLAGS) $(F77FLAGS) $(CFLAGS) $(CXXFLAGS) $(NVCCFLAGS) \
  $(CPP_EFFECTIVE_FLAGS) $(LDFLAGS) $(LAPACK_INC) $(LAPACK_LIBDIR) $(LAPACK_LIBS) \
  $(SOLVER_INC) $(SOLVER_LIBDIR) $(SOLVER_LIBS) $(GPU_PROVIDER_SRC) $(SELECTED_EXTERNAL_MODULE_DIRS) $(MPI_SRC) $(EOS_SRC) $(JAC_SRC) \
  $(SOLVER_PROVIDER_ID) $(LAPACK_PROVIDER_ID)
CONFIG_UNSUPPORTED_QUOTE := $(findstring ',$(CONFIG_RECORD_VALUES))
CONFIG_UNSUPPORTED_LF := $(if $(findstring $(XNET_CONFIG_NEWLINE),$(CONFIG_RECORD_VALUES)),yes)
CONFIG_UNSUPPORTED_CR := $(if $(findstring $(XNET_CONFIG_CR),$(CONFIG_RECORD_VALUES)),yes)
ifneq ($(strip $(CONFIG_UNSUPPORTED_QUOTE)$(CONFIG_UNSUPPORTED_LF)$(CONFIG_UNSUPPORTED_CR)),)
  $(error config.txt cannot represent single quotes or line breaks in recorded settings; choose a different spelling or BUILD_DIR)
endif

.PHONY: $(PUBLIC_TARGETS) FORCE
FORCE:
$(CONFIG): FORCE
	@test ! -L '$(@D)' || { echo 'refusing symlink BUILD_DIR' >&2; exit 2; }
	@mkdir -p '$(@D)'; { \
	 printf '%s\n' 'XNET_CONFIG_SCHEMA=1' 'SOURCE_ROOT=$(XNET_CFG_SOURCE_ROOT)' 'PE_ENV=$(XNET_CFG_PE_ENV)' 'CMODE=$(XNET_CFG_CMODE)' 'MACHINE=$(XNET_CFG_MACHINE)' 'MPI_MODE=$(XNET_CFG_MPI_MODE)' 'OPENMP_MODE=$(XNET_CFG_OPENMP_MODE)' 'GPU_MODE=$(XNET_CFG_GPU_MODE)' 'GPU_BACKEND=$(XNET_CFG_GPU_BACKEND)' 'GPU_LAPACK_VER=$(XNET_CFG_GPU_LAPACK_VER)' 'GPU_TARGET=$(XNET_CFG_GPU_TARGET)' 'OPENACC_MODE=$(XNET_CFG_OPENACC_MODE)' 'OPENMP_OL_MODE=$(XNET_CFG_OPENMP_OL_MODE)' 'EOS=$(XNET_CFG_EOS)' 'MATRIX_SOLVER=$(XNET_CFG_MATRIX_SOLVER)' 'LAPACK_VER=$(XNET_CFG_LAPACK_VER)' 'FC=$(COMPILE_FC)' 'CC=$(CC)' 'CXX=$(CXX)' 'NVCC=$(NVCC)' 'CPP=$(XNET_CPP)' 'LDR=$(LDR)' 'FFLAGS=$(COMPILE_FFLAGS)' 'F90FLAGS=$(COMPILE_F90FLAGS)' 'F77FLAGS=$(F77FLAGS)' 'CFLAGS=$(CFLAGS)' 'CXXFLAGS=$(CXXFLAGS)' 'NVCCFLAGS=$(NVCCFLAGS)' 'CPPFLAGS=$(CPP_EFFECTIVE_FLAGS)' 'LDFLAGS=$(LDFLAGS)' 'LAPACK_INC=$(LAPACK_INC)' 'LAPACK_LIBDIR=$(LAPACK_LIBDIR)' 'LAPACK_LIBS=$(LAPACK_LIBS)' 'SOLVER_INC=$(SOLVER_INC)' 'SOLVER_LIBDIR=$(SOLVER_LIBDIR)' 'SOLVER_LIBS=$(SOLVER_LIBS)' 'EXTERNAL_MODULE_DIRS=$(SELECTED_EXTERNAL_MODULE_DIRS)' 'PARALLEL_PROVIDER=$(MPI_SRC)' 'EOS_PROVIDERS=$(EOS_SRC)' 'JACOBIAN_PROVIDER=$(JAC_SRC)' 'GPU_PROVIDERS=$(GPU_PROVIDER_SRC)' 'SOLVER_PROVIDER=$(SOLVER_PROVIDER_ID)' 'LAPACK_PROVIDER=$(LAPACK_PROVIDER_ID)'; } > '$@.tmp'; \
	 if test -f '$@'; then cmp -s '$@' '$@.tmp' || { echo 'incompatible BUILD_DIR configuration; clean this BUILD_DIR or select another BUILD_DIR' >&2; rm -f '$@.tmp'; exit 2; }; rm -f '$@.tmp'; else mv '$@.tmp' '$@'; fi

$(SOURCE_OBJ_DIR) $(EOS_OBJ_DIR) $(SOLVER_OBJ_DIR) $(LAPACK_OBJ_DIR) $(GPU_OBJ_DIR) \
$(SOURCE_PP_DIR) $(EOS_PP_DIR) $(SOLVER_PP_DIR) $(LAPACK_PP_DIR) $(GPU_PP_DIR) \
$(MOD_DIR) $(BIN_DIR): | $(CONFIG)
	@mkdir -p '$@'

BUILD_LOGIC_INPUTS := $(XNET_DIR)/Makefile $(XNET_DIR)/Makefile.opt \
  $(XNET_DIR)/Makefile.internal $(XNET_DIR)/make/build.mk \
  $(XNET_DIR)/make/configuration.mk $(XNET_DIR)/make/providers.mk \
  $(XNET_DIR)/make/sources.mk $(XNET_DIR)/make/dependencies.mk \
  $(XNET_DIR)/make/rules.mk $(XNET_DIR)/make/machines.mk \
  $(XNET_DIR)/make/machines/generic.mk $(XNET_DIR)/make/machines/cray-pe.mk
PP_STATIC_INPUTS := $(XNET_DIR)/xnet_macros.fh $(BUILD_LOGIC_INPUTS) $(XNET_DIR)/crayftn_cpp.sh
