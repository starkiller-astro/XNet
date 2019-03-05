#  Makefile for FullNet and associated programs

# Architecture and installation specific MACROS
include $(HOME)/Codes/Makefile_OSX

FFLAGS 	= $(OPT) -qautodbl=dbl4 -qextname # active compiler options
LDRFLAGS=  $(OPT) -g  # linker flags
LIBS	=  $(LAPACK) # or $(PARDISO) Libraries need at linking
#
# These are the rules which make uses to build the executable
#
xnet: data.o full_net.o match.o net.o flux.o eosnom90.o jacobian_dense.o
	$(LDR) $(LDRFLAGS) -o xnet data.o full_net.o match.o net.o flux.o eosnom90.o jacobian_dense.o $(LIBS)
xnet_MA28: data.o full_net.o match.o net.o flux.o eosnom90.o jacobian_dense.o MA28.o
	$(LDR) $(LDRFLAGS) -o xnet data.o full_net.o match.o net.o flux.o eosnom90.o jacobian_dense.o MA28.o
xnet_PARDISO: data.o full_net.o match.o net.o flux.o eosnom90.o jacobian_PARDISO.o
	$(LDR) $(LDRFLAGS) -o xnet data.o full_net.o match.o net.o flux.o eosnom90.o jacobian_PARDISO.o $(LIBS)
net_setup: net_setup.o data.o
	$(LDR) $(LDRFLAGS) -o net_setup net_setup.o data.o
xinab: data.o init_abund.o 
	$(LDR)  $(LDRFLAGS) -o xinab init_abund.o data.o 
.f.o:
	$(FC) $(FFLAGS) -c $<
eosnom90.o: eosnom90.f
	$(FC)  $(FFLAGS) -c eosnom90.f
press.o: press.f
	$(FC)  $(FFLAGS) -c press.f
clean:
	rm -f core *.o 
tar:
	tar cvf net.tar ReadMe Changes Makefile *.f control th_const
#
#  Dependencies for object files
#
net.o:      full_net.o data.o flux.o match.o
flux.o:     full_net.o data.o match.o
match.o:    data.o
full_net.o: data.o 
jacobian_dense.o: data.o full_net.o
jacobian_ma28.o: data.o full_net.o
jacobian_PARDISO.o: data.o full_net.o
