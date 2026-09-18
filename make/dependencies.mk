# Manually maintained Fortran module build order for ordinary `make -j`.
# SOBJ and EOBJ map a module name to its object file. For example,
#
#   $(call SOBJ,xnet_constants): $(call SOBJ,xnet_types)
#
# means approximately:
#
#   build/.../obj/source/xnet_constants.o: build/.../obj/source/xnet_types.o
#
# because xnet_constants USEs the xnet_types module.
SOBJ = $(call source_obj,$(SOURCE_DIR)/$(1).F90)
EOBJ = $(call eos_obj,$(STARKILLER_HELMHOLTZ_PATH)/$(1).F90)
$(call SOBJ,xnet_constants): $(call SOBJ,xnet_types)
$(call SOBJ,xnet_parallel_stubs) $(call SOBJ,xnet_parallel): $(call SOBJ,xnet_types)
$(call SOBJ,xnet_util): $(call SOBJ,xnet_constants) $(call SOBJ,xnet_types) $(call source_obj,$(MPI_SRC))
$(call SOBJ,xnet_conditions): $(call SOBJ,xnet_controls) $(call SOBJ,xnet_types) $(call SOBJ,xnet_util)
$(call SOBJ,xnet_controls): $(call SOBJ,xnet_constants) $(call SOBJ,xnet_types) $(call SOBJ,xnet_util) $(call source_obj,$(MPI_SRC))
$(call SOBJ,xnet_data): $(call SOBJ,xnet_constants) $(call SOBJ,xnet_controls) $(call SOBJ,xnet_ffn) $(call SOBJ,xnet_nnu) $(call SOBJ,xnet_types) $(call SOBJ,xnet_util) $(call source_obj,$(MPI_SRC))
$(call SOBJ,xnet_fd) $(call SOBJ,xnet_timers): $(call SOBJ,xnet_types)
$(call SOBJ,xnet_ffn): $(call SOBJ,xnet_constants) $(call SOBJ,xnet_controls) $(call SOBJ,xnet_types) $(call SOBJ,xnet_conditions) $(call SOBJ,xnet_fd)
$(call SOBJ,xnet_abundances): $(call SOBJ,xnet_controls) $(call SOBJ,xnet_data) $(call SOBJ,xnet_types)
$(call SOBJ,xnet_nnu): $(call SOBJ,xnet_conditions) $(call SOBJ,xnet_controls) $(call SOBJ,xnet_types) $(call SOBJ,xnet_util)
$(call SOBJ,xnet_match): $(call SOBJ,xnet_controls) $(call SOBJ,xnet_data) $(call SOBJ,xnet_types) $(call SOBJ,xnet_util)
$(call SOBJ,xnet_preprocess): $(call SOBJ,xnet_constants) $(call SOBJ,xnet_data) $(call SOBJ,xnet_types) $(call SOBJ,xnet_util) $(call source_obj,$(MPI_SRC))
$(call SOBJ,xnet_gpu): $(call SOBJ,xnet_controls) $(call SOBJ,xnet_types) $(call SOBJ,xnet_util) $(GPU_OBJ)
$(call SOBJ,xnet_linalg): $(call SOBJ,xnet_constants) $(call SOBJ,xnet_controls) $(call SOBJ,xnet_gpu) $(call SOBJ,xnet_types) $(LAPACK_OBJ)
$(call SOBJ,xnet_screening): $(call SOBJ,xnet_abundances) $(call SOBJ,xnet_conditions) $(call SOBJ,xnet_constants) $(call SOBJ,xnet_controls) $(call SOBJ,xnet_data) $(call SOBJ,xnet_timers) $(call SOBJ,xnet_types) $(EOS_OBJ)
$(call EOBJ,eos_type): $(call SOBJ,xnet_types) $(call SOBJ,xnet_util)
$(call EOBJ,actual_eos): $(call EOBJ,eos_type) $(call SOBJ,xnet_util) $(call source_obj,$(MPI_SRC))
$(call eos_obj,$(SOURCE_DIR)/xnet_eos_starkiller.F90): $(call EOBJ,actual_eos) $(call EOBJ,eos_type) $(call SOBJ,xnet_fd) $(call SOBJ,xnet_constants) $(call SOBJ,xnet_controls) $(call SOBJ,xnet_data) $(call SOBJ,xnet_types) $(call SOBJ,xnet_util) $(call SOBJ,xnet_abundances)
$(call eos_obj,$(SOURCE_DIR)/xnet_eos_bahcall.F90): $(call SOBJ,xnet_abundances) $(call SOBJ,xnet_constants) $(call SOBJ,xnet_controls) $(call SOBJ,xnet_data) $(call SOBJ,xnet_fd) $(call SOBJ,xnet_types) $(call SOBJ,xnet_util)
$(call eos_obj,$(SOURCE_DIR)/xnet_eos_helm.F90): $(call SOBJ,xnet_abundances) $(call SOBJ,xnet_constants) $(call SOBJ,xnet_controls) $(call SOBJ,xnet_data) $(call SOBJ,xnet_fd) $(call SOBJ,xnet_types) $(call SOBJ,xnet_util)
$(call solver_obj,$(JAC_SRC)): $(call SOBJ,xnet_abundances) $(call SOBJ,xnet_conditions) $(call SOBJ,xnet_controls) $(call SOBJ,xnet_data) $(call SOBJ,xnet_linalg) $(if $(filter MA48 PARDISO_MKL,$(MATRIX_SOLVER)),$(call source_obj,$(MPI_SRC))) $(call SOBJ,xnet_timers) $(call SOBJ,xnet_types)
$(call SOBJ,net_setup): $(call SOBJ,xnet_controls) $(call SOBJ,xnet_preprocess)
$(call SOBJ,nse_slice): $(call SOBJ,xnet_controls) $(call SOBJ,xnet_match) $(call SOBJ,xnet_nse) $(call SOBJ,xnet_preprocess) $(call SOBJ,xnet_types) $(call SOBJ,xnet_data) $(call SOBJ,xnet_timers) $(EOS_OBJ)
$(call SOBJ,xnet_nse): $(call SOBJ,xnet_constants) $(call SOBJ,xnet_controls) $(call SOBJ,xnet_data) $(call SOBJ,xnet_timers) $(call SOBJ,xnet_types) $(call SOBJ,xnet_util) $(EOS_OBJ) $(LAPACK_OBJ)
$(call SOBJ,xnet_flux): $(call SOBJ,xnet_abundances) $(call SOBJ,xnet_conditions) $(call SOBJ,xnet_controls) $(call SOBJ,xnet_data) $(call SOBJ,xnet_match) $(call SOBJ,xnet_timers) $(call SOBJ,xnet_types)
$(call SOBJ,xnet_output): $(call SOBJ,xnet_abundances) $(call SOBJ,xnet_conditions) $(call SOBJ,xnet_controls) $(call SOBJ,xnet_data) $(call SOBJ,xnet_flux) $(call SOBJ,xnet_match) $(call source_obj,$(MPI_SRC)) $(call SOBJ,xnet_timers) $(call SOBJ,xnet_types)
$(call SOBJ,xnet_integrate): $(call SOBJ,xnet_abundances) $(call SOBJ,xnet_conditions) $(call SOBJ,xnet_constants) $(call SOBJ,xnet_controls) $(call SOBJ,xnet_data) $(call SOBJ,xnet_ffn) $(call SOBJ,xnet_linalg) $(call SOBJ,xnet_nnu) $(call SOBJ,xnet_screening) $(call SOBJ,xnet_timers) $(call SOBJ,xnet_types) $(call SOBJ,xnet_util) $(EOS_OBJ) $(LAPACK_OBJ)
$(call SOBJ,xnet_integrate_bdf) $(call SOBJ,xnet_integrate_be): $(call SOBJ,xnet_abundances) $(call SOBJ,xnet_conditions) $(call SOBJ,xnet_controls) $(call SOBJ,xnet_data) $(call SOBJ,xnet_integrate) $(call SOBJ,xnet_timers) $(call SOBJ,xnet_types) $(SOLVER_FREE_OBJ)
$(call SOBJ,xnet_evolve): $(call SOBJ,xnet_abundances) $(call SOBJ,xnet_conditions) $(call SOBJ,xnet_controls) $(call SOBJ,xnet_data) $(call SOBJ,xnet_integrate) $(call SOBJ,xnet_integrate_bdf) $(call SOBJ,xnet_integrate_be) $(call SOBJ,xnet_output) $(call SOBJ,xnet_timers) $(call SOBJ,xnet_types) $(call SOBJ,xnet_util)
$(call SOBJ,model_input_ascii): $(call SOBJ,xnet_abundances) $(call SOBJ,xnet_conditions) $(call SOBJ,xnet_controls) $(call SOBJ,xnet_data) $(call SOBJ,xnet_nnu) $(call SOBJ,xnet_nse) $(call SOBJ,xnet_types) $(call SOBJ,xnet_util)
$(call SOBJ,net): $(call SOBJ,xnet_abundances) $(call SOBJ,xnet_conditions) $(call SOBJ,xnet_controls) $(call SOBJ,xnet_data) $(call SOBJ,xnet_evolve) $(call SOBJ,xnet_flux) $(call SOBJ,xnet_gpu) $(call SOBJ,xnet_match) $(call SOBJ,xnet_nnu) $(call SOBJ,xnet_nse) $(call SOBJ,xnet_preprocess) $(call SOBJ,xnet_screening) $(call SOBJ,xnet_timers) $(call SOBJ,xnet_types) $(call SOBJ,xnet_util) $(EOS_OBJ) $(SOLVER_FREE_OBJ) $(call source_obj,$(MPI_SRC))
