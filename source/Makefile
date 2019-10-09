##
## Makefile for XNet and associated programs
##

## Build search path
VPATH = .

## Define XNet's source directory
XNET_DIR = .

## Include options for build configurations defined in Makefile.opt
include $(XNET_DIR)/Makefile.opt

#-------------------------------------------------------------------------------------------------
# The rules below are used by make to build a range of executables
# Users should not need to modify this portion of the Makefile
#

## Common XNet files
BASE = \
       xnet_controls.o \
       xnet_data.o \
       xnet_output.o \
       xnet_abundances.o \
       xnet_conditions.o \
       xnet_constants.o \
       xnet_evolve.o \
       xnet_ffn.o \
       xnet_flux.o \
       xnet_integrate.o \
       xnet_integrate_bdf.o \
       xnet_integrate_be.o \
       xnet_match.o \
       xnet_nnu.o \
       xnet_nse.o \
       xnet_preprocess.o \
       xnet_screening.o \
       xnet_timers.o \
       xnet_types.o \
       xnet_util.o

## Choose Driver
DRIVER = model_input_ascii.o net.o
ifeq ($(MPI_MODE),ON)
  MPI_OBJ = xnet_parallel.o
  MPI_SFX = _mpi
else
  MPI_OBJ = xnet_parallel_stubs.o
  MPI_SFX =
endif
BASE += $(MPI_OBJ)

## Configure Libraries and files to link appropriate EOS
ifeq ($(EOS),BAHCALL)
  EOS_OBJ = xnet_eos_bahcall.o
else ifeq ($(EOS),HELMHOLTZ)
  ifndef HELMHOLTZ_PATH
    $(error HELMHOLTZ_PATH variable must be set for EOS=HELMHOLTZ)
  endif
  EOS_OBJ = xnet_eos_helm.o xnet_fd.o helmholtz.o
  VPATH += $(HELMHOLTZ_PATH)
else ifeq ($(EOS),STARKILLER)
  EOS_OBJ = actual_eos.o eos_type.o xnet_fd.o xnet_eos_starkiller.o
  STARKILLER_HELMHOLTZ_PATH ?= ../tools/starkiller-helmholtz
  VPATH += $(STARKILLER_HELMHOLTZ_PATH)
endif
BASE += $(EOS_OBJ)

.DEFAULT_GOAL := $(EXE)

## Include system configurations and architecture/installation specific MACROS
include $(XNET_DIR)/Makefile.internal
-include $(XNET_DIR)/Makefile.dev

$(EXE): $(BASE) xnet_jacobian_$(SOLVER).o $(SOLVER_OBJ) $(LAPACK_OBJ) $(DRIVER)
	$(LDR) $(LDFLAGS) -o $(EXE) \
	    $(BASE) xnet_jacobian_$(SOLVER).o $(notdir $(SOLVER_OBJ) $(LAPACK_OBJ)) $(DRIVER) \
	    $(SOLVER_LIBDIR) $(SOLVER_LIBS) $(LAPACK_LIBDIR) $(LAPACK_LIBS)

all: xnet_dense xnet_MA48 xnet_PARDISO net_setup

test: test_serial test_parallel test_heat

test_serial: xnet_PARDISO
	@cd ../test ; ./test_xnet.sh 0

test_heat: xnet_PARDISO
	@cd ../test ; ./test_xnet.sh 50

test_simple: xnet_PARDISO
	@cd ../test ; ./test_xnet.sh 51

#
# Rules to build specific configurations of XNet and supporting tools
#

## Recursively call make so that the variable definitions reflect the targeted solver.
## This is the only sensible way to update the object file prerequisites based on the target.
ifndef update-prereq
xnet_dense xnet_MA41 xnet_MA48 xnet_PARDISO:
	$(MAKE) $(MAKECMDGOALS) MATRIX_SOLVER=$(@:xnet_%=%) update-prereq=y --no-print-directory
else
xnet_dense: $(BASE) xnet_jacobian_dense.o $(LAPACK_OBJ) $(DRIVER)
	$(LDR) $(LDFLAGS) -o xnetd$(MPI_SFX) \
	    $(BASE) xnet_jacobian_dense.o $(notdir $(LAPACK_OBJ)) $(DRIVER) \
	    $(LAPACK_LIBDIR) $(LAPACK_LIBS)

xnet_MA41: $(BASE) xnet_jacobian_MA41.o $(MA41_OBJ) $(LAPACK_OBJ) $(DRIVER)
	$(LDR) $(LDFLAGS) -o xnetm$(MPI_SFX) \
	    $(BASE) xnet_jacobian_MA41.o $(notdir $(MA41_OBJ) $(LAPACK_OBJ)) $(DRIVER) \
	    $(MA41_LIBDIR) $(MA41_LIBS) $(LAPACK_LIBDIR) $(LAPACK_LIBS)

xnet_MA48: $(BASE) xnet_jacobian_MA48.o $(MA48_OBJ) $(LAPACK_OBJ) $(DRIVER)
	$(LDR) $(LDFLAGS) -o xnetm$(MPI_SFX) \
	    $(BASE) xnet_jacobian_MA48.o $(notdir $(MA48_OBJ) $(LAPACK_OBJ)) $(DRIVER) \
	    $(MA48_LIBDIR) $(MA48_LIBS) $(LAPACK_LIBDIR) $(LAPACK_LIBS)
xnet_PARDISO: $(BASE) xnet_jacobian_$(PARDISO).o $(PARDISO_OBJ) $(LAPACK_OBJ) $(DRIVER)
	$(LDR) $(LDFLAGS) -o xnetp$(MPI_SFX) \
	    $(BASE) xnet_jacobian_$(PARDISO).o $(notdir $(PARDISO_OBJ) $(LAPACK_OBJ)) $(DRIVER) \
	    $(PARDISO_LIBDIR) $(PARDISO_LIBS) $(LAPACK_LIBDIR) $(LAPACK_LIBS)
endif

net_setup: xnet_conditions.o xnet_controls.o xnet_data.o net_setup.o xnet_constants.o xnet_ffn.o xnet_nnu.o xnet_preprocess.o xnet_util.o $(MPI_OBJ)
	$(LDR) $(LDFLAGS) -o net_setup \
	    xnet_conditions.o xnet_controls.o xnet_data.o net_setup.o xnet_constants.o xnet_ffn.o xnet_nnu.o xnet_preprocess.o xnet_util.o $(MPI_OBJ)

xinab: xnet_data.o init_abund.o 
	$(LDR) $(LDFLAGS) -o xinab \
	    init_abund.o xnet_data.o 

xnse: xnet_conditions.o xnet_controls.o xnet_data.o xnet_constants.o xnet_ffn.o xnet_nse.o xnet_preprocess.o xnet_util.o $(EOS_OBJ) $(MPI_OBJ) nse_slice.o
	$(LDR) $(LDFLAGS) -o xnse$(MPI_SFX) \
	    xnet_conditions.o xnet_controls.o xnet_data.o xnet_constants.o xnet_ffn.o xnet_nse.o xnet_preprocess.o xnet_util.o $(EOS_OBJ) $(MPI_OBJ) nse_slice.o \
	    $(notdir $(LAPACK_OBJ)) $(LAPACK_LIBDIR) $(LAPACK_LIBS)

#
# Rules for compiling individual files.
#
xnet_parallel.o: xnet_parallel.f90
	$(FC_MPI) $(FFLAGS) $(F90FLAGS) -c $< -o $@
xnet_eos_helm.o: xnet_eos_helm.f90
	$(FC) $(FFLAGS) $(F90FLAGS) -I$(HELMHOLTZ_PATH) -c $< -o $@
xnet_nse.o: xnet_nse.f90
	$(FC) $(FFLAGS) $(F90FLAGS) $(LAPACK_INC) -c $< -o $@
xnet_integrate.o: xnet_integrate.f90
	$(FC) $(FFLAGS) $(F90FLAGS) $(LAPACK_INC) -c $< -o $@
xnet_integrate_bdf.o: xnet_integrate_bdf.f90
	$(FC) $(FFLAGS) $(F90FLAGS) $(LAPACK_INC) -c $< -o $@
xnet_jacobian_$(SOLVER).o: xnet_jacobian_$(SOLVER).f90
	$(FC) $(FFLAGS) $(F90FLAGS) $(SOLVER_INC) $(LAPACK_INC) -c $< -o $@

$(LAPACK_OBJ_F90): %.o: %.f90
	$(FC) $(FFLAGS) $(F90FLAGS) $(LAPACK_INC) -c $< -o $(notdir $@)
$(LAPACK_OBJ_F): %.o: %.f
	$(FC) $(FFLAGS) $(F77FLAGS) $(LAPACK_INC) -c $< -o $(notdir $@)

$(SOLVER_OBJ_F90): %.o: %.f90
	$(FC) $(FFLAGS) $(F90FLAGS) $(SOLVER_INC) $(LAPACK_INC) -c $< -o $(notdir $@)
$(SOLVER_OBJ_F): %.o: %.f
	$(FC) $(FFLAGS) $(F77FLAGS) $(SOLVER_INC) $(LAPACK_INC) -c $< -o $(notdir $@)

%.o: %.F90
	$(FC) $(FFLAGS) $(F90FLAGS) -c $< -o $@
%.o: %.f90
	$(FC) $(FFLAGS) $(F90FLAGS) -c $< -o $@
%.o: %.f
	$(FC) $(FFLAGS) $(F77FLAGS) -c $< -o $@
%.o: %.c
	$(CC) $(CCFLAGS) -c $< -o $@
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f core *.o *.oo *.mod *.lst *.cub *.ptx *.i *.T *.diag xref.db
	rm -rf $(INLINE_DB) *.dSYM
tar:
	tar cvf net.tar ReadMe Changes Makefile *.f control th_const
print-%:
	@echo $* = $($*)
#
#  Dependencies for object files
#
#

xnet_controls.o: xnet_constants.o xnet_types.o xnet_util.o $(MPI_OBJ)
xnet_data.o: xnet_constants.o xnet_controls.o xnet_ffn.o xnet_nnu.o xnet_types.o xnet_util.o $(MPI_OBJ)
xnet_ffn.o: xnet_constants.o xnet_controls.o xnet_types.o
xnet_flux.o: xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_match.o xnet_timers.o xnet_types.o
xnet_match.o: xnet_controls.o xnet_data.o xnet_types.o
xnet_abundances.o: xnet_controls.o xnet_data.o xnet_timers.o xnet_types.o
xnet_conditions.o: xnet_controls.o xnet_types.o xnet_util.o
xnet_constants.o: xnet_types.o
xnet_evolve.o: xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_integrate.o xnet_integrate_bdf.o xnet_integrate_be.o xnet_output.o xnet_timers.o xnet_types.o xnet_util.o
xnet_gpu.o: cublasf.o cudaf.o openaccf.o xnet_controls.o
xnet_integrate.o: xnet_abundances.o xnet_conditions.o xnet_constants.o xnet_controls.o xnet_data.o xnet_ffn.o xnet_nnu.o xnet_screening.o xnet_timers.o xnet_types.o xnet_util.o $(EOS_OBJ) $(LAPACK_OBJ)
xnet_integrate_bd.o: xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_integrate.o xnet_timers.o xnet_types.o $(JAC_OBJ)
xnet_integrate_bdf.o: xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_integrate.o xnet_timers.o xnet_types.o $(JAC_OBJ)
xnet_integrate_be.o: xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_integrate.o xnet_timers.o xnet_types.o $(JAC_OBJ)
xnet_nnu.o: xnet_conditions.o xnet_controls.o xnet_types.o xnet_util.o
xnet_nse.o: xnet_constants.o xnet_controls.o xnet_data.o xnet_types.o xnet_util.o $(EOS_OBJ) $(LAPACK_OBJ)
xnet_output.o: xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_flux.o xnet_match.o xnet_timers.o xnet_types.o
xnet_preprocess.o: xnet_constants.o xnet_data.o xnet_types.o xnet_util.o $(MPI_OBJ)
xnet_screening.o: xnet_abundances.o xnet_conditions.o xnet_constants.o xnet_controls.o xnet_data.o xnet_types.o $(EOS_OBJ)
xnet_timers.o: xnet_types.o
xnet_util.o: xnet_constants.o xnet_types.o $(MPI_OBJ)

xnet_fd.o: xnet_types.o
xnet_eos_helm.o: xnet_fd.o helmholtz.o xnet_constants.o xnet_controls.o xnet_data.o xnet_types.o xnet_abundances.o 
xnet_eos_bahcall.o: xnet_constants.o xnet_controls.o xnet_data.o xnet_types.o xnet_abundances.o 
xnet_eos_starkiller.o: actual_eos.o eos_type.o xnet_fd.o xnet_constants.o xnet_controls.o xnet_data.o xnet_types.o xnet_abundances.o 

actual_eos.o: eos_type.o xnet_util.o $(MPI_OBJ)
eos_type.o: xnet_types.o xnet_util.o

xnet_jacobian_dense.o: xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_timers.o xnet_types.o xnet_util.o $(LAPACK_OBJ)
xnet_jacobian_MA41.o: xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_timers.o xnet_types.o xnet_util.o $(LAPACK_OBJ) $(MA41_OBJ) $(MPI_OBJ)
xnet_jacobian_MA48.o: xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_timers.o xnet_types.o xnet_util.o $(LAPACK_OBJ) $(MA48_OBJ) $(MPI_OBJ)
xnet_jacobian_PARDISO.o: xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_timers.o xnet_types.o xnet_util.o $(LAPACK_OBJ) $(PARDISO_OBJ) $(MPI_OBJ)
xnet_jacobian_PARDISO_MKL.o: xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_timers.o xnet_types.o xnet_util.o $(LAPACK_OBJ) $(PARDISO_OBJ) $(MPI_OBJ)
xnet_jacobian_gpu.o: cublasf.o cudaf.o xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_gpu.o xnet_timers.o xnet_types.o xnet_util.o $(LAPACK_OBJ) $(GPU_OBJ)

xnet_parallel_stubs.o: xnet_types.o
xnet_parallel.o: xnet_types.o

model_input_ascii.o: xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_nnu.o xnet_nse.o xnet_types.o xnet_util.o 
net.o: xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_evolve.o xnet_flux.o xnet_match.o xnet_nnu.o xnet_nse.o xnet_preprocess.o xnet_screening.o xnet_timers.o xnet_types.o xnet_util.o $(EOS_OBJ) $(JAC_OBJ) $(MPI_OBJ)
net_setup.o: xnet_controls.o xnet_preprocess.o

nse_slice.o: xnet_controls.o xnet_nse.o xnet_preprocess.o xnet_types.o $(EOS_OBJ)
