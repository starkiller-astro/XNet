#  Makefile for FullNet and associated programs

# Architecture and installation specific MACROS
#include $(HOME)/soft/Makefile
 include Makefile_ifort
#include /users/Shared/Work/Codes/Makefile_gfortran

FFLAGS 	= $(DEBUG) $(R8) # active compiler options
LDRFLAGS= $(DEBUG) $(R8) # linker flags
LIBS	=  # $(LAPACK) or $(PARDISO) #Libraries need at linking
#
# These are the rules which make uses to build the executable
#
# BASE = data.o common.o ffn.o full_net.o solve_be.o solve_bd.o  match.o net.o flux.o eosnom90.o screen.o 
 BASE =  data.o common.o ffn.o full_net.o screen.o qse_build.o solve_be.o solve_bd.o solve_qse.o match.o net.o flux.o eosnom90.o
 QBASE = data.o common.o ffn.o full_net.o qse_build.o solve_be.o solve_bd.o solve_qse.o match.o net.o flux.o eosnom90.o screen.o
xnet: $(BASE) jacobian_dense.o jacobian_qse.o
	$(LDR) $(LDRFLAGS) -o xnet $(BASE) jacobian_dense.o  jacobian_qse.o  $(LAPACK) $(LIBS)
mqsenet: $(BASE) jacobian_dense.o jacobian_qse.o  
	$(LDR) $(LDRFLAGS) -o mqsenet $(QBASE) jacobian_dense.o jacobian_qse.o $(LAPACK) $(LIBS)
qsenet: $(QBASE)  jacobian_qse.o 
	$(LDR) $(LDRFLAGS) -o qsenet $(QBASE)  jacobian_qse.o $(LAPACK) $(LIBS)
xnet_MA28: $(BASE) jacobian_ma28.o MA28.o
	$(LDR) $(LDRFLAGS) -o xnetm $(BASE) jacobian_ma28.o MA28.o
xnet_PARDISO: $(BASE) jacobian_PARDISO.o
	$(LDR) $(LDRFLAGS) -o xnetp $(BASE) jacobian_PARDISO.o $(PARDISO) $(LIBS)
xnet_nse: data.o ffn.o nse.o full_net.o match.o net_nse.o flux.o eosnom90.o jacobian_dense.o
	$(LDR) $(LDRFLAGS) -o xnetn data.o ffn.o full_net.o match.o net_nse.o flux.o eosnom90.o jacobian_dense.o nse.o $(LAPACK) $(LIBS)
net_setup: data.o ffn.o net_setup.o 
	$(LDR) $(LDRFLAGS) -o net_setup net_setup.o data.o ffn.o
xinab: data.o init_abund.o 
	$(LDR)  $(LDRFLAGS) -o xinab init_abund.o data.o 
.f90.o:
	$(FC) $(FFLAGS) -u -c $<
MA28.o: MA28.f
	$(FC)  $(FFLAGS) -c MA28.f
eosnom90.o: eosnom90.f
	$(FC)  $(FFLAGS) -c eosnom90.f
press.o: press.f
	$(FC)  $(FFLAGS) -c press.f
clean:
	rm -f core *.o *.mod
tar:
	tar cvf net.tar ReadMe Changes Makefile *.f control th_const
#
#  Dependencies for object files
#
net.o:      full_net.o data.o flux.o match.o common.o
flux.o:     full_net.o data.o match.o common.o
match.o:    data.o
full_net.o: data.o common.o
qse_build.o: data.o common.o flux.o match.o  
jacobian_qse.o: data.o full_net.o common.o qse_build.o
jacobian_dense.o: data.o full_net.o common.o
jacobian_ma28.o: data.o full_net.o common.o
jacobian_PARDISO.o: data.o full_net.o common.o
common.o: ffn.o 

