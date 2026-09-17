# XNet production build.
# One top-level GNU Make invocation owns a BUILD_DIR. Do not clean that
# directory while the build is running.
SHELL := /bin/sh
XNET_DIR := $(abspath .)
ROOT_DIR := $(abspath ..)
BUILD_BASE ?= build
ifeq ($(filter /%,$(BUILD_BASE)),)
  BUILD_BASE := $(abspath $(ROOT_DIR)/$(BUILD_BASE))
else
  BUILD_BASE := $(abspath $(BUILD_BASE))
endif
BUILD_NAME ?= default
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
ifeq ($(BUILD_DIR),$(XNET_DIR))
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

PUBLIC_TARGETS := xnet xnse net_setup all xnet_dense xnet_MA41 xnet_MA48 xnet_PARDISO \
                  frontier_gpu_linalg_probe xinab xnet_gpu clean clean-all \
                  print-XNET_EXE print-XNSE_EXE print-NET_SETUP_EXE print-PROBE_EXE
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
	@case '$(BUILD_DIR)' in /|'$(ROOT_DIR)'|'$(XNET_DIR)'|'$(BUILD_BASE)') echo 'refusing unsafe BUILD_DIR' >&2; exit 2;; esac; \
	 test ! -e '$(BUILD_DIR)' || { test ! -L '$(BUILD_DIR)' && test -f '$(CONFIG)' && test "$$(sed -n '1p' '$(CONFIG)')" = XNET_CONFIG_SCHEMA=1 || { echo 'refusing unmarked BUILD_DIR' >&2; exit 2; }; rm -rf -- '$(BUILD_DIR)'; }
clean-all:
	@test "$(CONFIRM_CLEAN_ALL)" = yes || { echo 'set CONFIRM_CLEAN_ALL=yes' >&2; exit 2; }
	@case '$(BUILD_BASE)' in /|'$(ROOT_DIR)'|'$(XNET_DIR)') echo 'refusing unsafe BUILD_BASE' >&2; exit 2;; esac; \
	 test ! -L '$(BUILD_BASE)' || { echo 'refusing symlink BUILD_BASE' >&2; exit 2; }; \
	 for d in '$(BUILD_BASE)'/*; do test -e "$$d" || continue; test -L "$$d" && continue; test -f "$$d/config.txt" && test "$$(sed -n '1p' "$$d/config.txt")" = XNET_CONFIG_SCHEMA=1 && rm -rf -- "$$d" || echo "leaving unrecognized $$d" >&2; done

else

include $(XNET_DIR)/make/configuration.mk
include $(XNET_DIR)/make/providers.mk
include $(XNET_DIR)/make/sources.mk
include $(XNET_DIR)/make/dependencies.mk
include $(XNET_DIR)/make/rules.mk

endif
