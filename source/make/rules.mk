# Preprocessing, compilation, link, and public target rules.
define free_rule
$(3): $(1) $(PP_STATIC_INPUTS) $(CONFIG) | $(call parent_dir,$(3))
ifeq ($(CRAY_RETAINED_PP),yes)
	@XNET_CPP_OUTPUT='$$@.$$$$.tmp' XNET_CPP='$(XNET_CPP)' XNET_CRAY_FTN=: \
	  $(XNET_DIR)/crayftn_cpp.sh $(CPP_EFFECTIVE_FLAGS) '$(1)' && mv '$$@.$$$$.tmp' '$$@'
else
	@$(XNET_CPP) $(CPP_EFFECTIVE_FLAGS) '$(1)' > '$$@.$$$$.tmp' && mv '$$@.$$$$.tmp' '$$@'
endif
$(2): $(3) | $(call parent_dir,$(2)) $(MOD_DIR)
	$(COMPILE_FC) $(COMPILE_FFLAGS) $(COMPILE_F90FLAGS) $(MODULE_FLAGS) $(LAPACK_INC) $(SOLVER_INC) \
	  $(CPP_INCLUDE_FLAGS) $(EXTERNAL_MODULE_FLAGS) -c '$$<' -o '$$@'
endef
$(foreach source,$(SOURCE_FREE_SRC),$(eval $(call free_rule,$(source),$(call source_obj,$(source)),$(call source_pp,$(source)))))
$(foreach source,$(EOS_FREE_SRC),$(eval $(call free_rule,$(source),$(call eos_obj,$(source)),$(call eos_pp,$(source)))))
$(foreach source,$(SOLVER_FREE_SRC),$(eval $(call free_rule,$(source),$(call solver_obj,$(source)),$(call solver_pp,$(source)))))
$(foreach source,$(LAPACK_FREE_SRC),$(eval $(call free_rule,$(source),$(call lapack_obj,$(source)),$(call lapack_pp,$(source)))))
$(foreach source,$(GPU_FREE_SRC),$(eval $(call free_rule,$(source),$(call gpu_obj,$(source)),$(call gpu_pp,$(source)))))
define fixed_rule
$(2): $(1) $(CONFIG) $(BUILD_LOGIC_INPUTS) | $(call parent_dir,$(2)) $(MOD_DIR)
	$(COMPILE_FC) $(FFLAGS) $(F77FLAGS) $(MODULE_FLAGS) $(LAPACK_INC) $(SOLVER_INC) \
	  $(EXTERNAL_MODULE_FLAGS) -c '$$<' -o '$$@'
endef
$(foreach source,$(SOLVER_FIXED_SRC),$(eval $(call fixed_rule,$(source),$(call solver_obj,$(source)))))
$(foreach source,$(LAPACK_FIXED_SRC),$(eval $(call fixed_rule,$(source),$(call lapack_obj,$(source)))))

define c_rule
$(2): $(1) $(CONFIG) $(BUILD_LOGIC_INPUTS) | $(call parent_dir,$(2))
	$(CC) $(CFLAGS) $(LAPACK_INC) $(SOLVER_INC) -c '$$<' -o '$$@'
endef
$(foreach source,$(SOLVER_C_SRC),$(eval $(call c_rule,$(source),$(call solver_obj,$(source)))))
$(foreach source,$(LAPACK_C_SRC),$(eval $(call c_rule,$(source),$(call lapack_obj,$(source)))))

define cxx_rule
$(2): $(1) $(CONFIG) $(BUILD_LOGIC_INPUTS) | $(call parent_dir,$(2))
	$(CXX) $(CXXFLAGS) $(LAPACK_INC) $(SOLVER_INC) -c '$$<' -o '$$@'
endef
$(foreach source,$(SOLVER_CXX_SRC),$(eval $(call cxx_rule,$(source),$(call solver_obj,$(source)))))
$(foreach source,$(LAPACK_CXX_SRC),$(eval $(call cxx_rule,$(source),$(call lapack_obj,$(source)))))

define cuda_rule
$(2): $(1) $(CONFIG) $(BUILD_LOGIC_INPUTS) | $(call parent_dir,$(2))
	$(NVCC) $(NVCCFLAGS) $(LAPACK_INC) $(SOLVER_INC) -c '$$<' -o '$$@'
endef
$(foreach source,$(SOLVER_CUDA_SRC),$(eval $(call cuda_rule,$(source),$(call solver_obj,$(source)))))
$(foreach source,$(LAPACK_CUDA_SRC),$(eval $(call cuda_rule,$(source),$(call lapack_obj,$(source)))))

xnet xnet_dense xnet_MA41 xnet_MA48 xnet_PARDISO: $(XNET_EXE)
xnse: $(XNSE_EXE)
net_setup: $(NET_SETUP_EXE)
all: xnet xnse net_setup
frontier_gpu_linalg_probe: $(PROBE_EXE)
xinab: ; @echo 'xinab is unsupported: init_abund.F90 is not tracked' >&2; exit 2
xnet_gpu: ; @echo 'xnet_gpu is unsupported; use xnet with explicit GPU selectors' >&2; exit 2
$(XNET_EXE): $(XNET_OBJ) | $(BIN_DIR)
	$(LDR) $(LDFLAGS) -o '$@.$$$$.tmp' $(XNET_OBJ) $(SOLVER_LIBDIR) $(SOLVER_LIBS) \
	  $(LAPACK_LIBDIR) $(LAPACK_LIBS) && mv '$@.$$$$.tmp' '$@'
$(XNSE_EXE): $(XNSE_OBJ) | $(BIN_DIR)
	$(LDR) $(LDFLAGS) -o '$@.$$$$.tmp' $(XNSE_OBJ) $(LAPACK_LIBDIR) $(LAPACK_LIBS) && mv '$@.$$$$.tmp' '$@'
$(NET_SETUP_EXE): $(SETUP_OBJ) | $(BIN_DIR)
	$(LDR) $(LDFLAGS) -o '$@.$$$$.tmp' $(SETUP_OBJ) && mv '$@.$$$$.tmp' '$@'
$(PROBE_EXE): $(PROBE_OBJ) | $(BIN_DIR)
	$(LDR) $(LDFLAGS) -o '$@.$$$$.tmp' $(PROBE_OBJ) $(SOLVER_LIBDIR) $(SOLVER_LIBS) \
	  $(LAPACK_LIBDIR) $(LAPACK_LIBS) && mv '$@.$$$$.tmp' '$@'
print-XNET_EXE: $(CONFIG) ; @echo 'XNET_EXE = $(XNET_EXE)'
print-XNSE_EXE: $(CONFIG) ; @echo 'XNSE_EXE = $(XNSE_EXE)'
print-NET_SETUP_EXE: $(CONFIG) ; @echo 'NET_SETUP_EXE = $(NET_SETUP_EXE)'
print-PROBE_EXE: $(CONFIG) ; @echo 'PROBE_EXE = $(PROBE_EXE)'
print-%: $(CONFIG) ; @echo '$* = $($*)'
