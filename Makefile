# XNet production GNU Make build.
#
# User-selectable defaults are in Makefile.opt. Compiler flags and library
# settings are in Makefile.internal, with host defaults under make/machines/.
# The make/*.mk files contain build-directory handling, implementation
# selection, source lists, module dependencies, and compile/link rules.
# Production Fortran sources are under source/.
include make/build.mk
