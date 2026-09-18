# Prefer the facility name exported by Lmod; otherwise use the short hostname.
# A command-line or environment MACHINE value takes normal Make precedence.
MACHINE ?= $(if $(LMOD_SYSTEM_NAME),$(LMOD_SYSTEM_NAME),$(shell uname -n | sed 's/\..*//'))

CRAY_PE_HOSTS := frontier perlmutter

ifneq ($(filter $(MACHINE),$(CRAY_PE_HOSTS)),)
  include $(ROOT_DIR)/make/machines/cray-pe.mk
else
  include $(ROOT_DIR)/make/machines/generic.mk
endif
