SHELL = /bin/bash

ROOT_DIR := $(abspath ../..)
XNET_DIR := $(ROOT_DIR)/source
NETLIB_DIR := $(ROOT_DIR)/tools/LAPACK
BUILD_NET_SOURCE_DIR := $(ROOT_DIR)/test/build_net
HELM_TABLE := $(ROOT_DIR)/tools/starkiller-helmholtz/helm_table.dat
STARKILLER_DIR := $(ROOT_DIR)/tools/starkiller-helmholtz

include $(ROOT_DIR)/Makefile.opt
include $(ROOT_DIR)/Makefile.internal

BUILD_DIR := build/$(shell echo $(CMODE) | tr A-Z a-z)
TEST_EXE := $(BUILD_DIR)/xnet_contract_tests
NSE_TEST_SOURCE := $(BUILD_DIR)/xnet_nse.F90
MODULE_FLAGS := -J$(BUILD_DIR) -I$(BUILD_DIR) -I$(XNET_DIR)
LAPACK_OBJECTS := $(patsubst $(NETLIB_DIR)/%.f,$(BUILD_DIR)/lapack/%.o,$(LAPACK_SRC))
PREPROCESS_BUILD_DIR := $(BUILD_DIR)/preprocess
PREPROCESS_MODULE_FLAGS := -J$(PREPROCESS_BUILD_DIR) -I$(PREPROCESS_BUILD_DIR) -I$(XNET_DIR)
PREPROCESS_DIRECT_EXE := $(PREPROCESS_BUILD_DIR)/preprocess_direct
PREPROCESS_NET_SETUP_EXE := $(PREPROCESS_BUILD_DIR)/net_setup
PREPROCESS_VERIFY_EXE := $(PREPROCESS_BUILD_DIR)/verify_preprocess
BUILD_NET_READER_EXE := $(PREPROCESS_BUILD_DIR)/verify_build_net_readers
NUCLEAR_INPUT_MUTATOR_EXE := $(PREPROCESS_BUILD_DIR)/mutate_nuclear_input
PREPROCESS_WORK_DIR := $(abspath $(BUILD_DIR)/preprocess-work)
PREPROCESS_FIXTURE_DIR := $(abspath fixtures/preprocess)
BUILD_NET_BUILD_DIR := $(BUILD_DIR)/build-net
BUILD_NET_MODULE_FLAGS := -J$(BUILD_NET_BUILD_DIR) -I$(BUILD_NET_BUILD_DIR)
BUILD_NET_EXE := $(BUILD_NET_BUILD_DIR)/build_net
BUILD_NET_WORK_DIR := $(abspath $(BUILD_DIR)/build-net-work)
EOS_BUILD_DIR := $(BUILD_DIR)/eos
EOS_WORK_DIR := $(abspath $(BUILD_DIR)/eos-work)
STARKILLER_EOS_EXE := $(EOS_BUILD_DIR)/starkiller/eos_component_tests
BAHCALL_EOS_EXE := $(EOS_BUILD_DIR)/bahcall/eos_component_tests
EOS_TEST_EXES := $(STARKILLER_EOS_EXE) $(BAHCALL_EOS_EXE)
PRODUCTION_BUILD_NAME := unit-$(shell echo $(CMODE) | tr A-Z a-z)
PRODUCTION_XNET_EXE := $(ROOT_DIR)/build/$(PRODUCTION_BUILD_NAME)/bin/xnet
MODEL_INPUT_BUILD_DIR := $(BUILD_DIR)/model-input
MODEL_INPUT_MODULE_FLAGS := -J$(MODEL_INPUT_BUILD_DIR) -I$(MODEL_INPUT_BUILD_DIR) -I$(XNET_DIR)
MODEL_INPUT_TEST_EXE := $(MODEL_INPUT_BUILD_DIR)/model_input_component_tests

OBJECTS := \
  $(BUILD_DIR)/testdrive.o \
  $(BUILD_DIR)/xnet_types.o \
  $(BUILD_DIR)/xnet_constants.o \
  $(BUILD_DIR)/xnet_test_stubs.o \
  $(BUILD_DIR)/xnet_timers.o \
  $(BUILD_DIR)/xnet_util.o \
  $(BUILD_DIR)/xnet_conditions.o \
  $(BUILD_DIR)/xnet_abundances.o \
  $(BUILD_DIR)/xnet_nnu.o \
  $(BUILD_DIR)/xnet_nse.o \
  $(LAPACK_OBJECTS) \
  $(BUILD_DIR)/test_xnet_contracts.o

PREPROCESS_CORE_OBJECTS := \
  $(PREPROCESS_BUILD_DIR)/xnet_types.o \
  $(PREPROCESS_BUILD_DIR)/xnet_constants.o \
  $(PREPROCESS_BUILD_DIR)/xnet_parallel_stubs.o \
  $(PREPROCESS_BUILD_DIR)/xnet_util.o \
  $(PREPROCESS_BUILD_DIR)/xnet_controls.o \
  $(PREPROCESS_BUILD_DIR)/xnet_conditions.o \
  $(PREPROCESS_BUILD_DIR)/xnet_fd.o \
  $(PREPROCESS_BUILD_DIR)/xnet_ffn.o \
  $(PREPROCESS_BUILD_DIR)/xnet_nnu.o \
  $(PREPROCESS_BUILD_DIR)/xnet_data.o \
  $(PREPROCESS_BUILD_DIR)/xnet_preprocess.o

PREPROCESS_VERIFY_OBJECTS := \
  $(PREPROCESS_CORE_OBJECTS) \
  $(PREPROCESS_BUILD_DIR)/xnet_match.o \
  $(PREPROCESS_BUILD_DIR)/xnet_abundances.o \
  $(PREPROCESS_BUILD_DIR)/xnet_timers.o \
  $(PREPROCESS_BUILD_DIR)/xnet_jacobian_PARDISO.o \
  $(PREPROCESS_BUILD_DIR)/pardiso_test_stubs.o \
  $(PREPROCESS_BUILD_DIR)/verify_preprocess.o

BUILD_NET_OBJECTS := \
  $(BUILD_NET_BUILD_DIR)/net_module.o \
  $(BUILD_NET_BUILD_DIR)/partf_module.o \
  $(BUILD_NET_BUILD_DIR)/ffn_module.o \
  $(BUILD_NET_BUILD_DIR)/nnu_module.o \
  $(BUILD_NET_BUILD_DIR)/file_module.o \
  $(BUILD_NET_BUILD_DIR)/reaclib_reader.o

.DEFAULT_GOAL := test

.PHONY: test eos-test model-input-test clean production-xnet

test: eos-test model-input-test $(TEST_EXE) $(PREPROCESS_DIRECT_EXE) \
    $(PREPROCESS_NET_SETUP_EXE) $(PREPROCESS_VERIFY_EXE) $(BUILD_NET_READER_EXE) \
    $(NUCLEAR_INPUT_MUTATOR_EXE) $(BUILD_NET_EXE) production-xnet
	@$(TEST_EXE)
	@bash run_preprocess_contracts.sh $(PREPROCESS_FIXTURE_DIR) $(PREPROCESS_WORK_DIR) \
	  $(abspath $(PREPROCESS_DIRECT_EXE)) $(abspath $(PREPROCESS_NET_SETUP_EXE)) \
	  $(abspath $(PREPROCESS_VERIFY_EXE))
	@python3 run_build_net_contracts.py $(BUILD_NET_WORK_DIR) $(abspath $(BUILD_NET_EXE)) \
	  $(abspath $(PREPROCESS_NET_SETUP_EXE)) $(abspath $(BUILD_NET_READER_EXE)) \
	  $(abspath $(NUCLEAR_INPUT_MUTATOR_EXE)) $(abspath $(PRODUCTION_XNET_EXE)) \
	  $(abspath $(HELM_TABLE))

eos-test: $(EOS_TEST_EXES)
	@bash run_eos_contracts.sh $(EOS_WORK_DIR) $(abspath $(STARKILLER_EOS_EXE)) \
	  $(abspath $(BAHCALL_EOS_EXE)) $(abspath $(HELM_TABLE))

model-input-test: $(MODEL_INPUT_TEST_EXE)
	@bash run_model_input_contracts.sh $(MODEL_INPUT_BUILD_DIR) $(MODEL_INPUT_TEST_EXE)

$(TEST_EXE): $(OBJECTS)
	$(LDR) -o $@ $(OBJECTS) $(LDFLAGS) $(LAPACK_LIBDIR) $(LAPACK_LIBS)

$(BUILD_DIR):
	mkdir -p $@

$(MODEL_INPUT_BUILD_DIR):
	mkdir -p $@

$(MODEL_INPUT_BUILD_DIR)/testdrive.o: vendor/test-drive/testdrive.F90 | $(MODEL_INPUT_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(MODEL_INPUT_MODULE_FLAGS) -c $< -o $@

$(MODEL_INPUT_BUILD_DIR)/xnet_types.o: $(XNET_DIR)/xnet_types.F90 | $(MODEL_INPUT_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(MODEL_INPUT_MODULE_FLAGS) -c $< -o $@

$(MODEL_INPUT_BUILD_DIR)/xnet_constants.o: $(XNET_DIR)/xnet_constants.F90 \
    $(MODEL_INPUT_BUILD_DIR)/xnet_types.o | $(MODEL_INPUT_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(MODEL_INPUT_MODULE_FLAGS) -c $< -o $@

$(MODEL_INPUT_BUILD_DIR)/xnet_parallel_stubs.o: $(XNET_DIR)/xnet_parallel_stubs.F90 \
    $(MODEL_INPUT_BUILD_DIR)/xnet_types.o | $(MODEL_INPUT_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(MODEL_INPUT_MODULE_FLAGS) -c $< -o $@

$(MODEL_INPUT_BUILD_DIR)/model_input_component_stubs.o: support/model_input_component_stubs.F90 \
    $(MODEL_INPUT_BUILD_DIR)/xnet_types.o | $(MODEL_INPUT_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(MODEL_INPUT_MODULE_FLAGS) -c $< -o $@

$(MODEL_INPUT_BUILD_DIR)/xnet_util.o: $(XNET_DIR)/xnet_util.F90 \
    $(MODEL_INPUT_BUILD_DIR)/xnet_constants.o $(MODEL_INPUT_BUILD_DIR)/xnet_parallel_stubs.o \
    | $(MODEL_INPUT_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(MODEL_INPUT_MODULE_FLAGS) -c $< -o $@

$(MODEL_INPUT_BUILD_DIR)/model_input_ascii.o: $(XNET_DIR)/model_input_ascii.F90 \
    $(MODEL_INPUT_BUILD_DIR)/model_input_component_stubs.o $(MODEL_INPUT_BUILD_DIR)/xnet_util.o \
    | $(MODEL_INPUT_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(MODEL_INPUT_MODULE_FLAGS) -c $< -o $@

$(MODEL_INPUT_BUILD_DIR)/test_model_input_contracts.o: test_model_input_contracts.F90 \
    $(MODEL_INPUT_BUILD_DIR)/testdrive.o $(MODEL_INPUT_BUILD_DIR)/model_input_ascii.o \
    | $(MODEL_INPUT_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(MODEL_INPUT_MODULE_FLAGS) -c $< -o $@

$(MODEL_INPUT_TEST_EXE): $(MODEL_INPUT_BUILD_DIR)/testdrive.o \
    $(MODEL_INPUT_BUILD_DIR)/xnet_types.o $(MODEL_INPUT_BUILD_DIR)/xnet_constants.o \
    $(MODEL_INPUT_BUILD_DIR)/xnet_parallel_stubs.o \
    $(MODEL_INPUT_BUILD_DIR)/model_input_component_stubs.o \
    $(MODEL_INPUT_BUILD_DIR)/xnet_util.o $(MODEL_INPUT_BUILD_DIR)/model_input_ascii.o \
    $(MODEL_INPUT_BUILD_DIR)/test_model_input_contracts.o
	$(LDR) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/lapack:
	mkdir -p $@

$(PREPROCESS_BUILD_DIR):
	mkdir -p $@

$(BUILD_NET_BUILD_DIR):
	mkdir -p $@

production-xnet: | $(BUILD_DIR)
	$(MAKE) -C $(ROOT_DIR) BUILD_NAME=$(PRODUCTION_BUILD_NAME) clean
	$(MAKE) -C $(ROOT_DIR) BUILD_NAME=$(PRODUCTION_BUILD_NAME) -j CMODE=$(CMODE) xnet

$(BUILD_DIR)/testdrive.o: vendor/test-drive/testdrive.F90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(MODULE_FLAGS) -c $< -o $@

$(BUILD_DIR)/xnet_types.o: $(XNET_DIR)/xnet_types.F90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(MODULE_FLAGS) -c $< -o $@

$(BUILD_DIR)/xnet_constants.o: $(XNET_DIR)/xnet_constants.F90 \
    $(BUILD_DIR)/xnet_types.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(MODULE_FLAGS) -c $< -o $@

$(BUILD_DIR)/xnet_test_stubs.o: support/xnet_test_stubs.F90 \
    $(BUILD_DIR)/xnet_types.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(MODULE_FLAGS) -c $< -o $@

$(BUILD_DIR)/xnet_timers.o: $(XNET_DIR)/xnet_timers.F90 \
    $(BUILD_DIR)/xnet_types.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(MODULE_FLAGS) -c $< -o $@

$(BUILD_DIR)/xnet_util.o: $(XNET_DIR)/xnet_util.F90 \
    $(BUILD_DIR)/xnet_constants.o $(BUILD_DIR)/xnet_test_stubs.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(MODULE_FLAGS) -c $< -o $@

$(BUILD_DIR)/xnet_conditions.o: $(XNET_DIR)/xnet_conditions.F90 \
    $(BUILD_DIR)/xnet_types.o $(BUILD_DIR)/xnet_test_stubs.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(MODULE_FLAGS) -c $< -o $@

$(BUILD_DIR)/xnet_abundances.o: $(XNET_DIR)/xnet_abundances.F90 \
    $(BUILD_DIR)/xnet_constants.o $(BUILD_DIR)/xnet_test_stubs.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(MODULE_FLAGS) -c $< -o $@

$(BUILD_DIR)/xnet_nnu.o: $(XNET_DIR)/xnet_nnu.F90 \
    $(BUILD_DIR)/xnet_conditions.o $(BUILD_DIR)/xnet_util.o \
    $(BUILD_DIR)/xnet_test_stubs.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(MODULE_FLAGS) -c $< -o $@

$(NSE_TEST_SOURCE): $(XNET_DIR)/xnet_nse.F90 | $(BUILD_DIR)
	cp $< $@

$(BUILD_DIR)/xnet_nse.o: $(NSE_TEST_SOURCE) \
    $(BUILD_DIR)/xnet_constants.o $(BUILD_DIR)/xnet_test_stubs.o \
    $(BUILD_DIR)/xnet_timers.o $(BUILD_DIR)/xnet_util.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(MODULE_FLAGS) -c $< -o $@

$(BUILD_DIR)/lapack/%.o: $(NETLIB_DIR)/%.f | $(BUILD_DIR)/lapack
	$(FC) $(FFLAGS) $(F77FLAGS) $(LAPACK_INC) -c $< -o $@

$(BUILD_DIR)/test_xnet_contracts.o: test_xnet_contracts.F90 $(OBJECTS:$(BUILD_DIR)/test_xnet_contracts.o=) \
    | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_DIRECT_EXE): $(PREPROCESS_CORE_OBJECTS) $(PREPROCESS_BUILD_DIR)/preprocess_direct.o
	$(LDR) -o $@ $^ $(LDFLAGS)

$(PREPROCESS_NET_SETUP_EXE): $(PREPROCESS_CORE_OBJECTS) $(PREPROCESS_BUILD_DIR)/net_setup.o
	$(LDR) -o $@ $^ $(LDFLAGS)

$(PREPROCESS_VERIFY_EXE): $(PREPROCESS_VERIFY_OBJECTS)
	$(LDR) -o $@ $^ $(LDFLAGS)

$(BUILD_NET_READER_EXE): $(PREPROCESS_VERIFY_OBJECTS:$(PREPROCESS_BUILD_DIR)/verify_preprocess.o=) \
    $(PREPROCESS_BUILD_DIR)/verify_build_net_readers.o
	$(LDR) -o $@ $^ $(LDFLAGS)

$(NUCLEAR_INPUT_MUTATOR_EXE): $(PREPROCESS_BUILD_DIR)/mutate_nuclear_input.o
	$(LDR) -o $@ $^ $(LDFLAGS)

$(PREPROCESS_BUILD_DIR)/xnet_types.o: $(XNET_DIR)/xnet_types.F90 | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_BUILD_DIR)/xnet_constants.o: $(XNET_DIR)/xnet_constants.F90 \
    $(PREPROCESS_BUILD_DIR)/xnet_types.o | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_BUILD_DIR)/xnet_parallel_stubs.o: $(XNET_DIR)/xnet_parallel_stubs.F90 \
    $(PREPROCESS_BUILD_DIR)/xnet_types.o | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_BUILD_DIR)/xnet_util.o: $(XNET_DIR)/xnet_util.F90 \
    $(PREPROCESS_BUILD_DIR)/xnet_constants.o $(PREPROCESS_BUILD_DIR)/xnet_parallel_stubs.o \
    | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_BUILD_DIR)/xnet_controls.o: $(XNET_DIR)/xnet_controls.F90 \
    $(PREPROCESS_BUILD_DIR)/xnet_constants.o $(PREPROCESS_BUILD_DIR)/xnet_util.o \
    | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_BUILD_DIR)/xnet_conditions.o: $(XNET_DIR)/xnet_conditions.F90 \
    $(PREPROCESS_BUILD_DIR)/xnet_controls.o $(PREPROCESS_BUILD_DIR)/xnet_util.o \
    | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_BUILD_DIR)/xnet_fd.o: $(XNET_DIR)/xnet_fd.F90 \
    $(PREPROCESS_BUILD_DIR)/xnet_types.o | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_BUILD_DIR)/xnet_ffn.o: $(XNET_DIR)/xnet_ffn.F90 \
    $(PREPROCESS_BUILD_DIR)/xnet_conditions.o $(PREPROCESS_BUILD_DIR)/xnet_fd.o \
    | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_BUILD_DIR)/xnet_nnu.o: $(XNET_DIR)/xnet_nnu.F90 \
    $(PREPROCESS_BUILD_DIR)/xnet_conditions.o $(PREPROCESS_BUILD_DIR)/xnet_util.o \
    | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_BUILD_DIR)/xnet_data.o: $(XNET_DIR)/xnet_data.F90 \
    $(PREPROCESS_BUILD_DIR)/xnet_controls.o $(PREPROCESS_BUILD_DIR)/xnet_ffn.o \
    $(PREPROCESS_BUILD_DIR)/xnet_nnu.o | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_BUILD_DIR)/xnet_preprocess.o: $(XNET_DIR)/xnet_preprocess.F90 \
    $(PREPROCESS_BUILD_DIR)/xnet_data.o $(PREPROCESS_BUILD_DIR)/xnet_util.o \
    | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_BUILD_DIR)/xnet_match.o: $(XNET_DIR)/xnet_match.F90 \
    $(PREPROCESS_BUILD_DIR)/xnet_data.o $(PREPROCESS_BUILD_DIR)/xnet_util.o \
    | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_BUILD_DIR)/xnet_abundances.o: $(XNET_DIR)/xnet_abundances.F90 \
    $(PREPROCESS_BUILD_DIR)/xnet_data.o | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_BUILD_DIR)/xnet_timers.o: $(XNET_DIR)/xnet_timers.F90 \
    $(PREPROCESS_BUILD_DIR)/xnet_types.o | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_BUILD_DIR)/xnet_jacobian_PARDISO.o: $(XNET_DIR)/xnet_jacobian_PARDISO.F90 \
    $(PREPROCESS_BUILD_DIR)/xnet_abundances.o $(PREPROCESS_BUILD_DIR)/xnet_conditions.o \
    $(PREPROCESS_BUILD_DIR)/xnet_data.o $(PREPROCESS_BUILD_DIR)/xnet_timers.o \
    | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_BUILD_DIR)/preprocess_direct.o: support/preprocess_direct.F90 \
    $(PREPROCESS_BUILD_DIR)/xnet_preprocess.o | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_BUILD_DIR)/net_setup.o: $(XNET_DIR)/net_setup.F90 \
    $(PREPROCESS_BUILD_DIR)/xnet_preprocess.o | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_BUILD_DIR)/mutate_nuclear_input.o: support/mutate_nuclear_input.F90 \
    | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_BUILD_DIR)/pardiso_test_stubs.o: support/pardiso_test_stubs.F90 \
    $(PREPROCESS_BUILD_DIR)/xnet_types.o | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_BUILD_DIR)/verify_preprocess.o: support/verify_preprocess.F90 \
    $(PREPROCESS_BUILD_DIR)/xnet_jacobian_PARDISO.o $(PREPROCESS_BUILD_DIR)/xnet_match.o \
    | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(PREPROCESS_BUILD_DIR)/verify_build_net_readers.o: support/verify_build_net_readers.F90 \
    $(PREPROCESS_BUILD_DIR)/xnet_jacobian_PARDISO.o $(PREPROCESS_BUILD_DIR)/xnet_match.o \
    | $(PREPROCESS_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(PREPROCESS_MODULE_FLAGS) -c $< -o $@

$(BUILD_NET_EXE): $(BUILD_NET_OBJECTS)
	$(LDR) -o $@ $^ $(LDFLAGS)

$(BUILD_NET_BUILD_DIR)/net_module.o: $(BUILD_NET_SOURCE_DIR)/net_module.f90 | $(BUILD_NET_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(BUILD_NET_MODULE_FLAGS) -c $< -o $@

$(BUILD_NET_BUILD_DIR)/partf_module.o: $(BUILD_NET_SOURCE_DIR)/partf_module.f90 \
    $(BUILD_NET_BUILD_DIR)/net_module.o | $(BUILD_NET_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(BUILD_NET_MODULE_FLAGS) -c $< -o $@

$(BUILD_NET_BUILD_DIR)/ffn_module.o: $(BUILD_NET_SOURCE_DIR)/ffn_module.f90 \
    $(BUILD_NET_BUILD_DIR)/net_module.o $(BUILD_NET_BUILD_DIR)/partf_module.o \
    | $(BUILD_NET_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(BUILD_NET_MODULE_FLAGS) -c $< -o $@

$(BUILD_NET_BUILD_DIR)/nnu_module.o: $(BUILD_NET_SOURCE_DIR)/nnu_module.f90 \
    $(BUILD_NET_BUILD_DIR)/net_module.o | $(BUILD_NET_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(BUILD_NET_MODULE_FLAGS) -c $< -o $@

$(BUILD_NET_BUILD_DIR)/file_module.o: $(BUILD_NET_SOURCE_DIR)/file_module.f90 \
    $(BUILD_NET_BUILD_DIR)/net_module.o $(BUILD_NET_BUILD_DIR)/partf_module.o \
    $(BUILD_NET_BUILD_DIR)/ffn_module.o $(BUILD_NET_BUILD_DIR)/nnu_module.o \
    | $(BUILD_NET_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(BUILD_NET_MODULE_FLAGS) -c $< -o $@

$(BUILD_NET_BUILD_DIR)/reaclib_reader.o: $(BUILD_NET_SOURCE_DIR)/reaclib_reader.f90 \
    $(BUILD_NET_BUILD_DIR)/file_module.o | $(BUILD_NET_BUILD_DIR)
	$(FC) $(FFLAGS) $(F90FLAGS) $(BUILD_NET_MODULE_FLAGS) -c $< -o $@

# Compile each EOS implementation against the same component checks. The
# arguments are the build subdirectory, source, test macro, and extra objects.
define eos_rules
$(EOS_BUILD_DIR)/$(1):
	mkdir -p $$@

$(EOS_BUILD_DIR)/$(1)/testdrive.o: vendor/test-drive/testdrive.F90 | $(EOS_BUILD_DIR)/$(1)
	$$(FC) $$(FFLAGS) $$(F90FLAGS) -J$(EOS_BUILD_DIR)/$(1) -I$(EOS_BUILD_DIR)/$(1) \
	  -c $$< -o $$@

$(EOS_BUILD_DIR)/$(1)/xnet_types.o: $$(XNET_DIR)/xnet_types.F90 | $(EOS_BUILD_DIR)/$(1)
	$$(FC) $$(FFLAGS) $$(F90FLAGS) -J$(EOS_BUILD_DIR)/$(1) -I$(EOS_BUILD_DIR)/$(1) \
	  -c $$< -o $$@

$(EOS_BUILD_DIR)/$(1)/xnet_constants.o: $$(XNET_DIR)/xnet_constants.F90 \
    $(EOS_BUILD_DIR)/$(1)/xnet_types.o | $(EOS_BUILD_DIR)/$(1)
	$$(FC) $$(FFLAGS) $$(F90FLAGS) -J$(EOS_BUILD_DIR)/$(1) -I$(EOS_BUILD_DIR)/$(1) \
	  -c $$< -o $$@

$(EOS_BUILD_DIR)/$(1)/xnet_parallel_stubs.o: $$(XNET_DIR)/xnet_parallel_stubs.F90 \
    $(EOS_BUILD_DIR)/$(1)/xnet_types.o | $(EOS_BUILD_DIR)/$(1)
	$$(FC) $$(FFLAGS) $$(F90FLAGS) -J$(EOS_BUILD_DIR)/$(1) -I$(EOS_BUILD_DIR)/$(1) \
	  -c $$< -o $$@

$(EOS_BUILD_DIR)/$(1)/eos_component_stubs.o: support/eos_component_stubs.F90 \
    $(EOS_BUILD_DIR)/$(1)/xnet_types.o | $(EOS_BUILD_DIR)/$(1)
	$$(FC) $$(FFLAGS) $$(F90FLAGS) -J$(EOS_BUILD_DIR)/$(1) -I$(EOS_BUILD_DIR)/$(1) \
	  -c $$< -o $$@

$(EOS_BUILD_DIR)/$(1)/xnet_util.o: $$(XNET_DIR)/xnet_util.F90 \
    $(EOS_BUILD_DIR)/$(1)/xnet_constants.o $(EOS_BUILD_DIR)/$(1)/xnet_parallel_stubs.o \
    | $(EOS_BUILD_DIR)/$(1)
	$$(FC) $$(FFLAGS) $$(F90FLAGS) -J$(EOS_BUILD_DIR)/$(1) -I$(EOS_BUILD_DIR)/$(1) \
	  -I$$(XNET_DIR) -c $$< -o $$@

$(EOS_BUILD_DIR)/$(1)/xnet_abundances.o: $$(XNET_DIR)/xnet_abundances.F90 \
    $(EOS_BUILD_DIR)/$(1)/xnet_constants.o $(EOS_BUILD_DIR)/$(1)/eos_component_stubs.o \
    | $(EOS_BUILD_DIR)/$(1)
	$$(FC) $$(FFLAGS) $$(F90FLAGS) -J$(EOS_BUILD_DIR)/$(1) -I$(EOS_BUILD_DIR)/$(1) \
	  -I$$(XNET_DIR) -c $$< -o $$@

$(EOS_BUILD_DIR)/$(1)/xnet_fd.o: $$(XNET_DIR)/xnet_fd.F90 \
    $(EOS_BUILD_DIR)/$(1)/xnet_types.o | $(EOS_BUILD_DIR)/$(1)
	$$(FC) $$(FFLAGS) $$(F90FLAGS) -J$(EOS_BUILD_DIR)/$(1) -I$(EOS_BUILD_DIR)/$(1) \
	  -c $$< -o $$@

$(EOS_BUILD_DIR)/$(1)/xnet_eos.o: $(2) $(4) \
    $(EOS_BUILD_DIR)/$(1)/xnet_abundances.o $(EOS_BUILD_DIR)/$(1)/xnet_fd.o \
    $(EOS_BUILD_DIR)/$(1)/xnet_util.o | $(EOS_BUILD_DIR)/$(1)
	$$(FC) $$(FFLAGS) $$(F90FLAGS) -J$(EOS_BUILD_DIR)/$(1) -I$(EOS_BUILD_DIR)/$(1) \
	  -I$$(XNET_DIR) -c $$< -o $$@

$(EOS_BUILD_DIR)/$(1)/test_eos_contracts.o: test_eos_contracts.F90 \
    $(EOS_BUILD_DIR)/$(1)/testdrive.o $(EOS_BUILD_DIR)/$(1)/xnet_eos.o \
    | $(EOS_BUILD_DIR)/$(1)
	$$(FC) $$(FFLAGS) $$(F90FLAGS) -D$(3) -J$(EOS_BUILD_DIR)/$(1) \
	  -I$(EOS_BUILD_DIR)/$(1) -c $$< -o $$@

$(EOS_BUILD_DIR)/$(1)/eos_component_tests: \
    $(EOS_BUILD_DIR)/$(1)/testdrive.o \
    $(EOS_BUILD_DIR)/$(1)/xnet_types.o \
    $(EOS_BUILD_DIR)/$(1)/xnet_constants.o \
    $(EOS_BUILD_DIR)/$(1)/xnet_parallel_stubs.o \
    $(EOS_BUILD_DIR)/$(1)/eos_component_stubs.o \
    $(EOS_BUILD_DIR)/$(1)/xnet_util.o \
    $(EOS_BUILD_DIR)/$(1)/xnet_abundances.o \
    $(EOS_BUILD_DIR)/$(1)/xnet_fd.o \
    $(4) \
    $(EOS_BUILD_DIR)/$(1)/xnet_eos.o \
    $(EOS_BUILD_DIR)/$(1)/test_eos_contracts.o
	$$(LDR) -o $$@ $$^ $$(LDFLAGS)
endef

$(EOS_BUILD_DIR)/starkiller/eos_type.o: $(STARKILLER_DIR)/eos_type.F90 \
    $(EOS_BUILD_DIR)/starkiller/xnet_types.o $(EOS_BUILD_DIR)/starkiller/xnet_util.o \
    | $(EOS_BUILD_DIR)/starkiller
	$(FC) $(FFLAGS) $(F90FLAGS) -J$(EOS_BUILD_DIR)/starkiller -I$(EOS_BUILD_DIR)/starkiller \
	  -I$(XNET_DIR) -c $< -o $@

$(EOS_BUILD_DIR)/starkiller/actual_eos.o: $(STARKILLER_DIR)/actual_eos.F90 \
    $(EOS_BUILD_DIR)/starkiller/eos_type.o $(EOS_BUILD_DIR)/starkiller/xnet_util.o \
    $(EOS_BUILD_DIR)/starkiller/xnet_parallel_stubs.o | $(EOS_BUILD_DIR)/starkiller
	$(FC) $(FFLAGS) $(F90FLAGS) -J$(EOS_BUILD_DIR)/starkiller -I$(EOS_BUILD_DIR)/starkiller \
	  -I$(XNET_DIR) -c $< -o $@

$(eval $(call eos_rules,starkiller,$(XNET_DIR)/xnet_eos_starkiller.F90,TEST_STARKILLER,$(EOS_BUILD_DIR)/starkiller/eos_type.o $(EOS_BUILD_DIR)/starkiller/actual_eos.o))
$(eval $(call eos_rules,bahcall,$(XNET_DIR)/xnet_eos_bahcall.F90,TEST_BAHCALL,))

clean:
	rm -rf build
