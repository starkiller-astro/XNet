This discussion covers version 5.0 of the network.

Compared to version 4, there is a division of the build and solve of the
jacobian from the full_net.f file.  For more details, see the repository
for changes.

Required Fortran files
-------------
net.f                The outer wrapper for post-processing use.
                     Contains the main program, which handles I/O and loops
                     over zones, and 3 output routines, one called at the
                     start of the evolution, one called at the end of the
                     evolution, and one called by full_net at the end of
                     each timestep.

data.f 	             The primary network data structures, and the routines
                     which load the data from the Data directory

full_net.f           The main network routines

One of jacobian_*.f  For large networks (ny > ~300), using MA28 or PARDISO
                     (both of which are sparse) solvers can significantly
                     increase the speed of computation and significantly
                     decrease the amount of memory usage by the network.
                     While they are both faster than the old dense LAPACK
                     solver, neither sparse solver is as robust as LAPACK
                     (PARDISO seems to be more robust than MA28).  The sparse
                     solvers  are much more sensitive to problems with the
                     libraries like rates that blow up at low temperatures.
                     To use the sparse solvers,  you must use net_setup5.0
                     as well (this determines the sparsity pattern for the
                     network).

jacobian_dense.f     Contains the dense Jacobian build and solves with LAPACK
                     or other dense solver libraries.  This file contains the
                     same netmatr routine as version 4 of xnet.

jacobian_PARDISO.f   Contains the sparse Jacobian build and solves with
                     PARDISO libraries from Uni Basel.  This is the
                     prefered sparse solver.

                     PARDISO is included in the Intel Math Kernal libraries
                     (MKL) or you can download the PARDISO library direct
                     from its authors at
                     http://www.computational.unibas.ch/cs/scicomp/software/pardiso/

                     If you are using the MKL version of PARDISO, you should
                     not call pardisoinit.  Comment out this call at the
                     end of the read_jacobian_data routine in jacobian_PARDISO.f.

jacobian_MA28.f      Contains the sparse Jacobian build and solves with MA28
                     package from the Harwell Subroutine Library Archive.

                     MA28 is proprietary software from Harwell, but as an
                     obsolete portion of their library, a free license can
                     be obtained from http://hsl.rl.ac.uk/archive/hslarchive.html.
                     If you use the MA28 code, register at the HSL archive.
                     Registration is free, and makes your use of the code legal.
                     Download both the Package and HSL dependencies and put
                     them in a single file called MA28.f.

net_setup.f          Standalone program which preprocesses the REACLIB formated
                     reaction rate libraries.  Makes binary data files read by XNet,
                     so this must be re-run for each architecture.

Optional Fortran files
----------------------
ffn.f               Contains the data structures and routines for using rates
                    tabulated in the format of Fuller, Fowler & Newman 1985.

eosnom90.f          Nomoto's EOS, for use with the reaction screening routines.

match.f	            Contains the data structures and input routines for
                    matching forward and reverse reactions.

flux.f              Contains the data structures and routines for calculating
                    the net flux of each matched reaction pair.

press.f             Numerical Recipes matrix solver.  It is better to use
                    LAPACK or some other library solver which has been
                    optimized for your platform.

nse.f               Contains the routines necessary to solve for Nuclear
                    Statisital equilibrium for a given density, temperature
                    and electron fraction.

nse_slice.f         Wrapper to use NSE alone.

Input files
-----------
control            Sets the values of the various control flags as well as
                   the locations of input files

th_const           An Example thermodynamic trajectory file.  This example
                   sets the conditions to be constant, but time variation
                   is expected.

Data               The input data, as prepared by net_setup.f, is expected
                   to be in a self contained directory.

Common Errors
-------------
The most common errors come from errors in the input of the initial
abundances or thermodynamic trajectory.  By default the network is
configured to output to the screen the time, Temp. and density for each
timestep.  Check that these are correct.

If you see an error message which begins,
"Time beyond thermodynamic range" ...
Then the current time ( and the Stop time) is larger than the last time value
in the thermo trajectory.  Check to see that the time and stop time values in
your thermo input file are being read in correctly.  Also check that the
parameter nhmx (in the thermo_data module in full_net.f) is larger than the
number of lines in your thermo input file.

If you see an error message which contains either,
"OUT OF RANGE IN EQSTEL" ... or
"ERROR IN EXTDEG " ...
Then you are running under conditions where the EOS is invalid.
Chances are, this is an I/0 error, usually too high a temperature.
Note for high temperature T>10GK, NSE is a better representation of the
composition than the network.  If you wish to run at T>10 GK, you'll
need to replace the EOS.

If the network runs for a while then dies, reporting a floating point error,
Then you've probably reached too low a temperature for the network to run.
The temperature at which problems occur depends strongly on the choice of
reactions and the number of isotopes.

If the network gives a segmentation error entering the netmatr routine, the
likely cause is an over-full stack.  In FORTRAN 90 and later, automatic arrays
are located on the stack, in contrast to allocated arrays, which are located
on the heap.  In xnet, the automatic arrays include the Jacobian matrix, which
is ny*ny*8 bytes.  On many systems the stack is limited by default to 8 MB,
so this problem may arise as ny approaches 1000.  One may use the csh limit
command to change the stack size.  Typing "limit" lists the current limits.
Typing "limit stacksize 12M" will set the stack size to 12 MB, sufficient for
ny=1000 or so.  The equivalent bash command is apparently ulimit, but Google
knows more about this than I.  Using one of the sparse solvers can also
significantly decrease the memory usage of xnet for large networks.
