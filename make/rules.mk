# Preprocessing, compilation, link, and public target rules.

# Free-form sources can live in several directories while their objects are
# grouped under obj/source, obj/eos, obj/solver, or obj/gpu. Each invocation
# below expands to two ordinary rules. For example:
#
#   pp/source/xnet_types.f90: source/xnet_types.F90
#   obj/source/xnet_types.o:  pp/source/xnet_types.f90
define compile_free_form
$(3): $(1) $(PP_STATIC_INPUTS) $(CONFIG) | $(call parent_dir,$(3))
ifeq ($(CRAY_OMP_PREPROCESS),yes)
	@XNET_CPP_OUTPUT='$$@.$$$$.tmp' XNET_CPP='$(XNET_CPP)' \
	  $(ROOT_DIR)/make/crayftn_cpp.sh $(CPP_EFFECTIVE_FLAGS) '$(1)' && mv '$$@.$$$$.tmp' '$$@'
else
	@$(XNET_CPP) $(CPP_EFFECTIVE_FLAGS) '$(1)' > '$$@.$$$$.tmp' && mv '$$@.$$$$.tmp' '$$@'
endif
$(2): $(3) | $(call parent_dir,$(2)) $(MOD_DIR)
	$(COMPILE_FC) $(COMPILE_FFLAGS) $(COMPILE_F90FLAGS) $(MODULE_FLAGS) \
	  $(LAPACK_INC) $(SOLVER_INC) $(CPP_INCLUDE_FLAGS) $(GPU_MODULE_FLAGS) \
	  -c '$$<' -o '$$@'
endef

$(foreach source,$(SOURCE_FREE_SRC),$(eval $(call compile_free_form,$(source),$(call source_obj,$(source)),$(call source_pp,$(source)))))
$(foreach source,$(EOS_FREE_SRC),$(eval $(call compile_free_form,$(source),$(call eos_obj,$(source)),$(call eos_pp,$(source)))))
$(foreach source,$(SOLVER_FREE_SRC),$(eval $(call compile_free_form,$(source),$(call solver_obj,$(source)),$(call solver_pp,$(source)))))
$(foreach source,$(GPU_FREE_SRC),$(eval $(call compile_free_form,$(source),$(call gpu_obj,$(source)),$(call gpu_pp,$(source)))))

# Fixed-form MA41/MA48 and bundled NETLIB sources do not require a retained
# preprocessing step.
define compile_fixed_form
$(2): $(1) $(CONFIG) $(BUILD_LOGIC_INPUTS) | $(call parent_dir,$(2)) $(MOD_DIR)
	$(COMPILE_FC) $(FFLAGS) $(F77FLAGS) $(MODULE_FLAGS) $(LAPACK_INC) \
	  $(SOLVER_INC) $(GPU_MODULE_FLAGS) -c '$$<' -o '$$@'
endef

$(foreach source,$(SOLVER_FIXED_SRC),$(eval $(call compile_fixed_form,$(source),$(call solver_obj,$(source)))))
$(foreach source,$(LAPACK_FIXED_SRC),$(eval $(call compile_fixed_form,$(source),$(call lapack_obj,$(source)))))

xnet: $(XNET_EXE)
xnse: $(XNSE_EXE)
net_setup: $(NET_SETUP_EXE)
all: xnet xnse net_setup
xinab: ; @echo 'xinab is unsupported: init_abund.F90 is not tracked' >&2; exit 2
xnet_gpu: ; @echo 'xnet_gpu is unsupported; use xnet with explicit GPU selectors' >&2; exit 2

$(XNET_EXE): $(XNET_OBJ) | $(BIN_DIR)
	$(LDR) $(LDFLAGS) -o '$@.$$$$.tmp' $(XNET_OBJ) $(SOLVER_LIBDIR) $(SOLVER_LIBS) \
	  $(LAPACK_LIBDIR) $(LAPACK_LIBS) && mv '$@.$$$$.tmp' '$@'
$(XNSE_EXE): $(XNSE_OBJ) | $(BIN_DIR)
	$(LDR) $(LDFLAGS) -o '$@.$$$$.tmp' $(XNSE_OBJ) $(LAPACK_LIBDIR) $(LAPACK_LIBS) && mv '$@.$$$$.tmp' '$@'
$(NET_SETUP_EXE): $(SETUP_OBJ) | $(BIN_DIR)
	$(LDR) $(LDFLAGS) -o '$@.$$$$.tmp' $(SETUP_OBJ) && mv '$@.$$$$.tmp' '$@'

print-%: $(CONFIG)
	@echo '$* = $($*)'
