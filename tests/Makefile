SHELL := /bin/bash

FC = gfortran

BUILD_ROOT := build
BUILD_DIR ?= $(BUILD_ROOT)
UNIT_BUILD_DIR := $(BUILD_ROOT)/unit
SELFCHECK_BUILD_DIR := $(BUILD_ROOT)/selfcheck
STRICT_KIND_BUILD_DIR := $(BUILD_ROOT)/strict-kind
MOD_DIR := $(BUILD_DIR)/mod
UNIT_EXE := $(BUILD_DIR)/unit_tests
FAILURE_EXE := $(BUILD_DIR)/failure_fixture
STRICT_KIND_EXE := $(BUILD_DIR)/strict_kind_tests

SOURCE_DIR := ../source
VENDOR_DIR := vendor/test-drive-v0.6.1

F90FLAGS = -g -Og -cpp -fimplicit-none -fallow-argument-mismatch \
           -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none
STRICT_KIND_F90FLAGS := -g -Og -cpp -fimplicit-none \
                        -fallow-argument-mismatch -ffree-line-length-none
COMPILE_FLAGS := $(F90FLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -I$(SOURCE_DIR)

PRODUCTION_OBJECTS := \
	$(BUILD_DIR)/xnet_types.o \
	$(BUILD_DIR)/xnet_timers.o \
	$(BUILD_DIR)/xnet_fd.o \
	$(BUILD_DIR)/xnet_parallel_stubs.o \
	$(BUILD_DIR)/xnet_constants.o \
	$(BUILD_DIR)/xnet_util.o

TESTDRIVE_OBJECT := $(BUILD_DIR)/testdrive.o
UNIT_OBJECTS := \
	$(BUILD_DIR)/test_types.o \
	$(BUILD_DIR)/test_timers.o \
	$(BUILD_DIR)/test_fd0h.o \
	$(BUILD_DIR)/test_safe_exp.o \
	$(BUILD_DIR)/test_input_parsing.o \
	$(BUILD_DIR)/test_name_ordered.o \
	$(BUILD_DIR)/unit_test_main.o
FAILURE_OBJECT := $(BUILD_DIR)/failure_fixture.o
STRICT_KIND_OBJECTS := \
	$(TESTDRIVE_OBJECT) \
	$(BUILD_DIR)/xnet_types.o \
	$(BUILD_DIR)/xnet_constants.o \
	$(BUILD_DIR)/test_constants.o

.PHONY: test run_unit strict_kind run_strict_kind selfcheck run_selfcheck failure_fixture_status \
	clean clean_build

test:
	@$(MAKE) --no-print-directory clean_build BUILD_DIR="$(UNIT_BUILD_DIR)"
	@$(MAKE) --no-print-directory run_unit BUILD_DIR="$(UNIT_BUILD_DIR)"
	@$(MAKE) --no-print-directory strict_kind

run_unit: $(UNIT_EXE)
	$(UNIT_EXE)

strict_kind:
	@$(MAKE) --no-print-directory clean_build BUILD_DIR="$(STRICT_KIND_BUILD_DIR)"
	@$(MAKE) --no-print-directory run_strict_kind \
	  BUILD_DIR="$(STRICT_KIND_BUILD_DIR)" F90FLAGS="$(STRICT_KIND_F90FLAGS)"

run_strict_kind: $(STRICT_KIND_EXE)
	$(STRICT_KIND_EXE)

selfcheck:
	@$(MAKE) --no-print-directory clean_build BUILD_DIR="$(SELFCHECK_BUILD_DIR)"
	@$(MAKE) --no-print-directory run_selfcheck BUILD_DIR="$(SELFCHECK_BUILD_DIR)"

run_selfcheck: $(UNIT_EXE) $(FAILURE_EXE)
	@$(UNIT_EXE)
	@if $(FAILURE_EXE); then \
	  echo "ERROR: the deliberate assertion failure returned zero"; \
	  exit 1; \
	else \
	  echo "PASS: the deliberate assertion failure returned nonzero"; \
	fi
	@if $(MAKE) --no-print-directory failure_fixture_status; then \
	  echo "ERROR: Make did not propagate the deliberate assertion failure"; \
	  exit 1; \
	else \
	  echo "PASS: Make propagated the deliberate assertion failure"; \
	fi

failure_fixture_status: $(FAILURE_EXE)
	$(FAILURE_EXE)

$(UNIT_EXE): $(TESTDRIVE_OBJECT) $(PRODUCTION_OBJECTS) $(UNIT_OBJECTS)
	$(FC) $(F90FLAGS) -o $@ $^

$(FAILURE_EXE): $(TESTDRIVE_OBJECT) $(FAILURE_OBJECT)
	$(FC) $(F90FLAGS) -o $@ $^

$(STRICT_KIND_EXE): $(STRICT_KIND_OBJECTS)
	$(FC) $(F90FLAGS) -o $@ $^

$(BUILD_DIR) $(MOD_DIR):
	mkdir -p $@

$(BUILD_DIR)/testdrive.o: $(VENDOR_DIR)/src/testdrive.F90 | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(COMPILE_FLAGS) -c $< -o $@

$(BUILD_DIR)/xnet_types.o: $(SOURCE_DIR)/xnet_types.F90 | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(COMPILE_FLAGS) -c $< -o $@

$(BUILD_DIR)/xnet_timers.o: $(SOURCE_DIR)/xnet_timers.F90 \
		$(BUILD_DIR)/xnet_types.o | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(COMPILE_FLAGS) -c $< -o $@

$(BUILD_DIR)/xnet_fd.o: $(SOURCE_DIR)/xnet_fd.F90 \
		$(BUILD_DIR)/xnet_types.o | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(COMPILE_FLAGS) -c $< -o $@

$(BUILD_DIR)/xnet_parallel_stubs.o: $(SOURCE_DIR)/xnet_parallel_stubs.F90 \
		$(BUILD_DIR)/xnet_types.o | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(COMPILE_FLAGS) -c $< -o $@

$(BUILD_DIR)/xnet_constants.o: $(SOURCE_DIR)/xnet_constants.F90 \
		$(BUILD_DIR)/xnet_types.o | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(COMPILE_FLAGS) -c $< -o $@

$(BUILD_DIR)/xnet_util.o: $(SOURCE_DIR)/xnet_util.F90 \
		$(BUILD_DIR)/xnet_types.o \
		$(BUILD_DIR)/xnet_parallel_stubs.o \
		$(BUILD_DIR)/xnet_constants.o | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(COMPILE_FLAGS) -c $< -o $@

$(BUILD_DIR)/test_types.o: unit/test_types.F90 \
		$(TESTDRIVE_OBJECT) $(BUILD_DIR)/xnet_types.o | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(COMPILE_FLAGS) -c $< -o $@

$(BUILD_DIR)/test_timers.o: unit/test_timers.F90 \
		$(TESTDRIVE_OBJECT) $(BUILD_DIR)/xnet_types.o \
		$(BUILD_DIR)/xnet_timers.o | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(COMPILE_FLAGS) -c $< -o $@

$(BUILD_DIR)/test_fd0h.o: unit/test_fd0h.F90 \
		$(TESTDRIVE_OBJECT) $(BUILD_DIR)/xnet_types.o \
		$(BUILD_DIR)/xnet_fd.o | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(COMPILE_FLAGS) -c $< -o $@

$(BUILD_DIR)/test_safe_exp.o: unit/test_safe_exp.F90 \
		$(TESTDRIVE_OBJECT) $(BUILD_DIR)/xnet_util.o | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(COMPILE_FLAGS) -c $< -o $@

$(BUILD_DIR)/test_input_parsing.o: unit/test_input_parsing.F90 \
		$(TESTDRIVE_OBJECT) $(BUILD_DIR)/xnet_types.o \
		$(BUILD_DIR)/xnet_util.o | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(COMPILE_FLAGS) -c $< -o $@

$(BUILD_DIR)/test_name_ordered.o: unit/test_name_ordered.F90 \
		$(TESTDRIVE_OBJECT) $(BUILD_DIR)/xnet_util.o | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(COMPILE_FLAGS) -c $< -o $@

$(BUILD_DIR)/unit_test_main.o: unit/main.F90 \
		$(TESTDRIVE_OBJECT) $(BUILD_DIR)/test_types.o $(BUILD_DIR)/test_timers.o \
		$(BUILD_DIR)/test_fd0h.o \
		$(BUILD_DIR)/test_safe_exp.o $(BUILD_DIR)/test_input_parsing.o \
		$(BUILD_DIR)/test_name_ordered.o | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(COMPILE_FLAGS) -c $< -o $@

$(BUILD_DIR)/failure_fixture.o: selfcheck/failure_fixture.F90 \
		$(TESTDRIVE_OBJECT) | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(COMPILE_FLAGS) -c $< -o $@

$(BUILD_DIR)/test_constants.o: unit/test_constants.F90 \
		$(TESTDRIVE_OBJECT) $(BUILD_DIR)/xnet_types.o \
		$(BUILD_DIR)/xnet_constants.o | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(COMPILE_FLAGS) -c $< -o $@

clean:
	@$(MAKE) --no-print-directory clean_build BUILD_DIR="$(BUILD_ROOT)"
	@$(MAKE) --no-print-directory clean_build BUILD_DIR="$(UNIT_BUILD_DIR)"
	@$(MAKE) --no-print-directory clean_build BUILD_DIR="$(SELFCHECK_BUILD_DIR)"
	@$(MAKE) --no-print-directory clean_build BUILD_DIR="$(STRICT_KIND_BUILD_DIR)"

clean_build:
	rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/unit_tests $(BUILD_DIR)/failure_fixture \
	  $(BUILD_DIR)/strict_kind_tests
	rm -f $(MOD_DIR)/*.mod $(MOD_DIR)/*.smod
