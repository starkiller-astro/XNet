# Frontier qualification build additions.
#
# The production build owns XNet and its reusable objects. This file adds only
# the facility probe used by this qualification; it is not a production XNet
# program.
include make/build.mk

PROBE_SRC := $(ROOT_DIR)/test/qualification/frontier/gpu_linalg_probe.F90
PROBE_PP := $(PP_DIR)/qualification/frontier/gpu_linalg_probe.f90
PROBE_OBJ := $(OBJ_DIR)/qualification/frontier/gpu_linalg_probe.o
PROBE_EXE := $(BIN_DIR)/frontier_gpu_linalg_probe
PROBE_SUPPORT_OBJ := $(filter-out $(call SOBJ,model_input_ascii) $(call SOBJ,net),$(XNET_OBJ))

.PHONY: frontier_gpu_linalg_probe
frontier_gpu_linalg_probe: $(PROBE_EXE)

$(dir $(PROBE_PP)) $(dir $(PROBE_OBJ)): | $(CONFIG)
	mkdir -p '$@'

$(PROBE_PP): $(PROBE_SRC) $(PP_STATIC_INPUTS) $(CONFIG) | $(dir $(PROBE_PP))
ifeq ($(CRAY_OMP_PREPROCESS),yes)
	@XNET_CPP_OUTPUT='$@.$$$$.tmp' XNET_CPP='$(XNET_CPP)' \
	  $(ROOT_DIR)/make/crayftn_cpp.sh $(CPP_EFFECTIVE_FLAGS) '$<' && mv '$@.$$$$.tmp' '$@'
else
	@$(XNET_CPP) $(CPP_EFFECTIVE_FLAGS) '$<' > '$@.$$$$.tmp' && mv '$@.$$$$.tmp' '$@'
endif

$(PROBE_OBJ): $(PROBE_PP) $(call SOBJ,xnet_controls) $(call SOBJ,xnet_gpu) \
    $(call SOBJ,xnet_linalg) $(call SOBJ,xnet_types) | $(dir $(PROBE_OBJ)) $(MOD_DIR)
	$(COMPILE_FC) $(COMPILE_FFLAGS) $(COMPILE_F90FLAGS) $(MODULE_FLAGS) \
	  $(LAPACK_INC) $(SOLVER_INC) $(CPP_INCLUDE_FLAGS) $(GPU_MODULE_FLAGS) \
	  -c '$<' -o '$@'

$(PROBE_EXE): $(PROBE_SUPPORT_OBJ) $(PROBE_OBJ) | $(BIN_DIR)
	$(LDR) $(LDFLAGS) -o '$@.$$$$.tmp' $(PROBE_SUPPORT_OBJ) $(PROBE_OBJ) \
	  $(SOLVER_LIBDIR) $(SOLVER_LIBS) $(LAPACK_LIBDIR) $(LAPACK_LIBS) && mv '$@.$$$$.tmp' '$@'
