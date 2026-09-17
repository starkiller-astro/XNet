# XNet production build.
# One top-level GNU Make invocation owns a BUILD_DIR. Do not clean that
# directory while the build is running.
SHELL := /bin/sh
ROOT_DIR := $(CURDIR)
SOURCE_DIR := $(ROOT_DIR)/source

# Read the public defaults before choosing the build directory. Compiler flags,
# machine defaults, library paths, and provider selection are resolved later.
include Makefile.opt

# Build output. Without an explicit BUILD_NAME or BUILD_DIR, use a readable
# name derived from the major public selectors. config.txt remains responsible
# for detecting changes to flags, library paths, and other recorded settings.
BUILD_BASE ?= build
ifeq ($(filter /%,$(BUILD_BASE)),)
  BUILD_BASE := $(abspath $(ROOT_DIR)/$(BUILD_BASE))
else
  BUILD_BASE := $(abspath $(BUILD_BASE))
endif
AUTO_BUILD_NAME := $(PE_ENV)-$(CMODE)
ifeq ($(MPI_MODE),ON)
  AUTO_BUILD_NAME := $(AUTO_BUILD_NAME)-MPI
endif
ifeq ($(OPENMP_MODE),ON)
  AUTO_BUILD_NAME := $(AUTO_BUILD_NAME)-OPENMP
endif
ifeq ($(GPU_MODE),ON)
  AUTO_BUILD_NAME := $(AUTO_BUILD_NAME)-$(GPU_BACKEND)-$(GPU_LAPACK_VER)
  ifeq ($(OPENACC_MODE),ON)
    AUTO_BUILD_NAME := $(AUTO_BUILD_NAME)-OPENACC
  else ifeq ($(OPENMP_OL_MODE),ON)
    AUTO_BUILD_NAME := $(AUTO_BUILD_NAME)-OPENMP-OFFLOAD
  endif
endif
ifneq ($(EOS),STARKILLER)
  AUTO_BUILD_NAME := $(AUTO_BUILD_NAME)-$(EOS)
endif
ifneq ($(MATRIX_SOLVER),dense)
  AUTO_BUILD_NAME := $(AUTO_BUILD_NAME)-$(MATRIX_SOLVER)
endif
BUILD_NAME ?= $(AUTO_BUILD_NAME)
ifneq ($(words $(BUILD_NAME)),1)
  $(error BUILD_NAME must be one readable path component without whitespace)
endif
ifneq ($(findstring /,$(BUILD_NAME)),)
  $(error BUILD_NAME must be one readable path component)
endif
ifneq ($(filter . ..,$(BUILD_NAME)),)
  $(error BUILD_NAME may not be '.' or '..')
endif
BUILD_DIR ?= $(BUILD_BASE)/$(BUILD_NAME)
ifeq ($(filter /%,$(BUILD_DIR)),)
  BUILD_DIR := $(abspath $(ROOT_DIR)/$(BUILD_DIR))
else
  BUILD_DIR := $(abspath $(BUILD_DIR))
endif
ifneq ($(words $(BUILD_DIR)),1)
  $(error BUILD_DIR paths containing whitespace are not supported)
endif
ifeq ($(BUILD_DIR),$(ROOT_DIR))
  $(error BUILD_DIR may not be the repository root)
endif
ifeq ($(BUILD_DIR),$(SOURCE_DIR))
  $(error BUILD_DIR may not be the source directory)
endif
ifeq ($(BUILD_DIR),$(BUILD_BASE))
  $(error BUILD_DIR may not equal BUILD_BASE)
endif

OBJ_DIR := $(BUILD_DIR)/obj
MOD_DIR := $(BUILD_DIR)/mod
PP_DIR := $(BUILD_DIR)/pp
BIN_DIR := $(BUILD_DIR)/bin
CONFIG := $(BUILD_DIR)/config.txt

REQUESTED_GOALS := $(if $(MAKECMDGOALS),$(MAKECMDGOALS),xnet)
CLEAN_GOALS := $(filter clean clean-all,$(REQUESTED_GOALS))
OTHER_GOALS := $(filter-out clean clean-all,$(REQUESTED_GOALS))
ifneq ($(CLEAN_GOALS),)
  ifneq ($(OTHER_GOALS),)
    $(error clean targets cannot be combined with build or print targets)
  endif
  ifneq ($(words $(CLEAN_GOALS)),1)
    $(error select exactly one of clean or clean-all)
  endif
endif

.DEFAULT_GOAL := xnet

# Cleaning is deliberately independent of compiler and configuration validation.
ifneq ($(CLEAN_GOALS),)
.PHONY: clean clean-all
clean:
	@case '$(BUILD_DIR)' in /|'$(ROOT_DIR)'|'$(SOURCE_DIR)'|'$(BUILD_BASE)') echo 'refusing unsafe BUILD_DIR' >&2; exit 2;; esac; \
	 test ! -e '$(BUILD_DIR)' || { test ! -L '$(BUILD_DIR)' && test -f '$(CONFIG)' && test "$$(sed -n '1p' '$(CONFIG)')" = XNET_CONFIG_SCHEMA=1 || { echo 'refusing unmarked BUILD_DIR' >&2; exit 2; }; rm -rf -- '$(BUILD_DIR)'; }
clean-all:
	@test "$(CONFIRM_CLEAN_ALL)" = yes || { echo 'set CONFIRM_CLEAN_ALL=yes' >&2; exit 2; }
	@case '$(BUILD_BASE)' in /|'$(ROOT_DIR)'|'$(SOURCE_DIR)') echo 'refusing unsafe BUILD_BASE' >&2; exit 2;; esac; \
	 test ! -L '$(BUILD_BASE)' || { echo 'refusing symlink BUILD_BASE' >&2; exit 2; }; \
	 for d in '$(BUILD_BASE)'/*; do test -e "$$d" || continue; test -L "$$d" && continue; test -f "$$d/config.txt" && test "$$(sed -n '1p' "$$d/config.txt")" = XNET_CONFIG_SCHEMA=1 && rm -rf -- "$$d" || echo "leaving unrecognized $$d" >&2; done

else

include make/configuration.mk
include make/providers.mk
include make/sources.mk
include make/dependencies.mk
include make/rules.mk

endif
