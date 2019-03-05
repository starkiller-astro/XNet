MACHINE=$(shell echo $(HOSTNAME) | sed 's/[0-9]\+$$//')
CRAY_MACHINE = \
	       titan-ext \
	       hopper \
	       darter \
	       edison \
	       cori \
	       chester \
	       beacon \
	       mira
ifeq (,$(filter $(MACHINE),$(CRAY_MACHINE)))
    # If not a Cray machine, set PE_ENV as it would be on a Cray machine, unless already set
    PE_ENV  ?= GNU
#   PE_ENV  ?= PGI
#   PE_ENV  ?= INTEL
#   PE_ENV  ?= CRAY
#   PE_ENV  ?= ABSOFT
else
    F9X = ftn
endif

#CMODE = DEBUG
CMODE ?= OPTIMIZE

ifeq ($(PE_ENV),GNU)
    #GNU
    F9X            ?= gfortran
    FFLAGS_DEBUG    = -O0 -g -fbacktrace -fcheck=all -Wall -ffree-line-length-none
    FFLAGS_OPTIMIZE = -O3 -march=native -ffree-line-length-none
else ifeq ($(PE_ENV),PGI)
    #PGI
    F9X            ?= pgf90
    FFLAGS_DEBUG    = -O0 -g -traceback -Mbounds -Mchkstk -Minform=warn
    FFLAGS_OPTIMIZE = -fast
else ifeq ($(PE_ENV),INTEL)
    #Intel
    F9X            ?= ifort
    FFLAGS_DEBUG    = -O0 -g -debug -traceback -C -CB -debug-parameters all -ftrapuv -fpe0
    FFLAGS_OPTIMIZE = -O2 -xHost
else ifeq ($(PE_ENV),CRAY)
    #Cray
    F9X            ?= crayftn
    FFLAGS_DEBUG    = -O0 -eD -h noomp
    FFLAGS_OPTIMIZE = -O2 -h noomp
endif

FLINKER = $(F9X)
FFLAGS  = $(FFLAGS_$(CMODE))

SUBROUTINES = reaclib_reader.o
MODULES = net_module.o partf_module.o ffn_module.o nnu_module.o file_module.o

EXE = build_net

.DEFAULT_GOAL := $(EXE)

$(EXE): $(MODULES) $(SUBROUTINES)
	$(FLINKER) $(FFLAGS) -o $(EXE) $(SUBROUTINES) $(MODULES) $(FLIBS)

scrub: clean
	rm -rf $(EXE)
clean:
	rm -rf *.o *.mod $(EXE).dSYM

%.o : %.f90
	$(F9X) $(FFLAGS) -c $< -o $@
%.o : %.f
	$(F9X) $(FFLAGS) -c $< -o $@
