##
## Makefile for XNet and associated programs
##

## Build search path
VPATH = .

## Define XNet's source directory
XNET_DIR = .

## Inlcude options for build configurations defined in Makefile.opt
include $(XNET_DIR)/Makefile.opt

#-------------------------------------------------------------------------------------------------
# The rules below are used by make to build a range of executables
# Users should not need to modify this portion of the Makefile
#

## Common XNet files
BASE = control.o data.o conditions.o $(FFN_OBJ) $(NNU_OBJ) common.o net_preprocess.o $(NSE_OBJ) full_net.o solve_be.o solve_bd.o match.o flux.o screening.o
MPI_FILES = data_distribute_mpi.o jacobian_bcast_dense.o jacobian_bcast_PARDISO.o jacobian_bcast_MA28.o net_mpi.o

# Define Driver
DRIVER_SM = net.o
DRIVER_MPI = data_distribute_mpi.o net_mpi.o jacobian_bcast_$(MATRIX_SOLVER).o

## Configure Libraries and files to link appropriate EOS
ifeq ($(EOS),BAHCALL)
  EOS_OBJ = eos_bahcall.o
else ifeq ($(EOS),HELMHOLTZ)
  EOS_OBJ = eos_helm.o helmholtz.o
  VPATH += $(HELMHOLTZ_PATH)
endif
BASE += $(EOS_OBJ)

.DEFAULT_GOAL := $(EXE)

## Inlcude system configurations and architecture/installation specific MACROS
include $(XNET_DIR)/Makefile.internal
-include $(XNET_DIR)/Makefile.dev

$(EXE): $(BASE) jacobian_$(MATRIX_SOLVER).o $(DRIVER) $(EXTRA_OBJ) 
	$(LDRL) $(LDRFLAGS) -o $(EXE) \
	    $(BASE) jacobian_$(MATRIX_SOLVER).o $(notdir $(EXTRA_OBJ)) $(DRIVER) \
	    $(EXTRA_LINK) $(LIBS)

#
# Rules to build specific configurations of XNet and supporting tools
#
xnet_dense: $(BASE) jacobian_dense.o $(DRIVER_SM) 
	$(LDR) $(LDRFLAGS) -o xnetd \
	    $(BASE) jacobian_dense.o $(notdir $(LAPACK_OBJ)) $(DRIVER_SM) \
	    $(LAPACK_LIBDIR) $(LAPACK_LIBS) $(EXTRA_LINK)

xnet_MA28: $(BASE) jacobian_MA28.o MA28.o $(EXTRA_OBJ)
	$(LDR) $(LDRFLAGS) -o xnetm \
	    $(BASE) jacobian_MA28.o $(DRIVER_SM) MA28.o \
	    $(LAPACK_LIBDIR) $(LAPACK_LIBS) $(EXTRA_LINK)

xnet_MA41: $(BASE) jacobian_MA41.o MA41.o $(EXTRA_OBJ)
	$(LDR) $(LDRFLAGS) -o xnetm \
	    $(BASE) jacobian_MA41.o $(DRIVER_SM) MA41.o \
	    $(LAPACK_LIBDIR) $(LAPACK_LIBS) $(EXTRA_LINK)

xnet_MA48: $(BASE) jacobian_MA48.o MA48.o $(EXTRA_OBJ)
	$(LDR) $(LDRFLAGS) -o xnetm \
	    $(BASE) jacobian_MA48.o $(DRIVER_SM) MA48.o \
	    $(LAPACK_LIBDIR) $(LAPACK_LIBS) $(EXTRA_LINK)

xnet_PARDISO: $(BASE) jacobian_PARDISO.o $(DRIVER) $(EXTRA_OBJ)
	$(LDR) $(LDRFLAGS) -o xnetp \
	    $(BASE) jacobian_PARDISO.o $(DRIVER_SM) \
	    $(PARDISO_INC) $(PARDISO_LIBDIR) $(PARDISO_LIBS) $(EXTRA_LINK)

xnet_nse: $(BASE) jacobian_$(MATRIX_SOLVER).o nse.o $(EXTRA_OBJ)
	sed -e 's/^\!\(.*\!NSE\)/ \1/' -e 's/^[^\!]\(.*\!NOTNSE\)/!\1/' net.f90 >| net_nse.f90
	sed -e 's/^\!\(.*\!NSE\)/ \1/' -e 's/^[^\!]\(.*\!NOTNSE\)/!\1/' conditions.f90 >| conditions_nse.f90
	$(LDR) $(LDRFLAGS) -o xnet_nse \
	    $(filter-out conditions.o,$(BASE)) jacobian_$(MATRIX_SOLVER).o nse.o $(notdir $(EXTRA_OBJ)) conditions_nse.f90 net_nse.f90 \
	    $(EXTRA_LINK) $(LIBS)

xnet_MPI: $(BASE) jacobian_$(MATRIX_SOLVER).o $(DRIVER_MPI) $(EXTRA_OBJ)
	$(LDR_MPI) $(LDRFLAGS) -o xnet_mpi \
	    $(BASE) jacobian_$(MATRIX_SOLVER).o $(notdir $(EXTRA_OBJ)) $(DRIVER_MPI) \
	    $(EXTRA_LINK) $(LIBS)

xnet_dense_MPI: $(BASE) jacobian_dense.o $(DRIVER_MPI) $(EXTRA_OBJ)
	$(LDR_MPI) $(LDRFLAGS) -o xnetd_mpi \
	    $(BASE) jacobian_dense.o $(notdir $(LAPACK_OBJ)) $(DRIVER_MPI) \
	    $(LAPACK_LIBDIR) $(LAPACK_LIBS) $(EXTRA_LINK)

xnet_PARDISO_MPI: $(BASE) jacobian_PARDISO.o $(DRIVER_MPI) $(EXTRA_OBJ)
	$(LDR_MPI) $(LDRFLAGS) -o xnetp_mpi \
	    $(BASE) jacobian_PARDISO.o $(DRIVER_MPI) \
	    $(PARDISO_INC) $(PARDISO_LIBDIR) $(PARDISO_LIBS) $(EXTRA_LINK)

net_setup: control.o data.o conditions.o ffn.o net_preprocess.o net_setup.o $(NNU_OBJ)
	$(LDR) $(LDRFLAGS) -o net_setup net_setup.o net_preprocess.o control.o data.o conditions.o ffn.o $(NNU_OBJ)

xinab: data.o init_abund.o 
	$(LDR) $(LDRFLAGS) -o xinab \
	    init_abund.o data.o 

xnse: control.o data.o conditions.o ffn.o net_preprocess.o nse.o $(EOS_OBJ) nse_slice.o
	$(LDR) $(LDRFLAGS) -o xnse \
	    control.o data.o conditions.o ffn.o net_preprocess.o nse.o $(EOS_OBJ) nse_slice.o

#
# Rules for compiling individual files.
#
$(MPI_FILES): %.o: %.f90
	$(FC_MPI) $(FFLAGS) -c $< -o $@
eosnom90.o: eosnom90.f
	$(FC) $(FFLAGS) -c eosnom90.f -o $@
eos_helm.o: eos_helm.f90
	$(FC) $(FFLAGS) -I$(HELMHOLTZ_PATH) -c eos_helm.f90 -o $@
press.o: press.f
	$(FC) $(FFLAGS) -c press.f -o $@
net.o: net.f90
	$(FC) $(FFLAGS) -c net.f90 -o $@
common.o: common.f90
	$(FC) $(FFLAGS) $(LAPACK_INC) -c common.f90 -o $@
jacobian_$(MATRIX_SOLVER).o: jacobian_$(MATRIX_SOLVER).f90
	$(FC) $(FFLAGS) $(LAPACK_INC) $(SOLVER_INC) -c jacobian_$(MATRIX_SOLVER).f90 -o $@
solve_be.o: solve_be.f90
	$(FC) $(FFLAGS) -c solve_be.f90 -o $@

$(LAPACK_OBJ_F90): %.o: %.f90
	$(FC) $(FFLAGS) $(LAPACK_INC) -c $< -o $(notdir $@)
$(LAPACK_OBJ_F): %.o: %.f
	$(FC) $(FFLAGS) $(LAPACK_INC) -c $< -o $(notdir $@)

$(SOLVER_OBJ_F90): %.o: %.f90
	$(FC) $(FFLAGS) $(SOLVER_INC) $(LAPACK_INC) -c $< -o $(notdir $@)
$(SOLVER_OBJ_F): %.o: %.f
	$(FC) $(FFLAGS) $(SOLVER_INC) $(LAPACK_INC) -c $< -o $(notdir $@)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f core *.o *.oo *.mod *.lst *.cub *.ptx *.i *.T *.diag xref.db
	rm -rf $(INLINE_DB)
tar:
	tar cvf net.tar ReadMe Changes Makefile *.f control th_const
#
#  Dependencies for object files
#
net.o:      full_net.o data.o flux.o match.o common.o control.o
flux.o:     full_net.o data.o match.o common.o
match.o:    data.o
full_net.o: data.o common.o
solve_be.o: jacobian_$(MATRIX_SOLVER).o
nu_nucleus.o: conditions.o
jacobian_$(MATRIX_SOLVER).o: $(EXTRA_OBJ) data.o full_net.o common.o
common.o:  ffn.o
