!*******************************************************************************
!  Full_Net 4.10 2/26/07 
!  The subroutines in this file perform nucleosynthesis for a single Lagrangian 
!  mass zone.  The handling of multiple zones, including multi-processing, is 
!  carried out externally.  
!
!  The data.f file contains the nuclear and reaction data structures and the 
!  subroutines which read in the data and allocate the arrays.  
!*******************************************************************************

      module controls
!===============================================================================
!  This module contains the values of the flags and limits which control 
!  the behavior of the network
!===============================================================================
      integer :: iweak     ! If =0, weak reaction are ignored
                           ! If <0, only weak reactions are used
      integer :: iscrn     ! If =0, screening is ignored
      integer :: itso      ! Sets level of time series output
      integer :: idiag     ! Sets level of diagnostic output
      integer :: iconvc    ! Controls type of convergence condition (0=mass)
      integer :: kstmx     ! Max # of timesteps before exit
      integer :: knrmx     ! Max # of Newton-Raphson iterations
      integer :: nzone     ! Number of zones
      real(8) :: tolc      ! The iterative convergence test limit
      real(8) :: tolm      ! Max network mass error
      real(8) :: changemx  ! Relative abundance change used to guess timestep
      real(8) :: ytime     ! Min abundance included in timestep estimate
      real(8) :: ymin      ! Min abundance, y < ymin =0
      real(8) :: tdelmm    ! new timestep <= tdelmm * previous timestep
      common /flag/ iweak,iscrn,itso,idiag,iconvc,kstmx,knrmx,nzone
      common /tol/ changemx,tolm,tolc,ytime,ymin,tdelmm
      end module controls

      module conditions
!===============================================================================
!  This module contains data on the current time and thermodynamic conditions.
!===============================================================================
      real(8) :: t        ! Time at the beginning of the current timestep
      real(8) :: tt       ! Trial time for end of current timestep
      real(8) :: tdel     ! Trial duration of timestep
      real(8) :: t9t,rhot ! Temperature (GK) and Density (g/cc) at trial time
      end module conditions

      module abundances
!===============================================================================
!  This module contains the abundances of the nuclear species, at the previous 
!  time (yo), current time (y), and trial time (yt), as well as the time 
!  derivatives (ydot), at the trial time, and the abundance change due to the 
!  Newton-Raphson iteration.  The size of these arrays is allocated in the 
!  external routine, where the initial values are set.
!===============================================================================
      use nuc_number
      real(8), dimension(:), allocatable :: yo,y,yt,ydot,dy
      end module abundances

      module thermo_data
!===============================================================================
!  This module contains the thermodynamic trajectory which the network follows
!===============================================================================
      integer, parameter   :: nhmx=3000 ! The max number of thermo points
      integer              :: nh        ! The actual number of points
      real(8) :: tstart,tstop,th(nhmx),t9h(nhmx),rhoh(nhmx)
      end module thermo_data

      module constants
!===============================================================================
!  These are fundamental constants in CGS, except energies in MeV, Temp in GK
!===============================================================================
      real(8), parameter :: pi   =3.1415926536, hbar=6.582122E-22
      real(8), parameter :: amu  =1.036427E-18, bok =8.617385E-02
      real(8), parameter :: avn  =6.022137E+23, e2  =1.439830E-13
      real(8), parameter :: epmev=1.602177E-06
      end module constants

      subroutine full_net(izone)
!===============================================================================
!  The abundances are evolved using a Newton-Raphson iteration scheme to solve 
!  the equation yt(i)-y(i)/tdel = ydot(i), where y(i) is the abundance at the 
!  beginning of the iteration, yt(i) is the trial abundance at the end of the 
!  timestep, ydot(i) is the time derivative of the trial abundance  calculated 
!  from reaction rates and tdel is the timestep.  
!===============================================================================
      use controls
      use nuclear_data
      use conditions
      use abundances
      use thermo_data
      real(8) :: reldy(ny)
      real(8) :: toln      !  Network convergence condition
      real(8) :: testc,testc2,testm,testn    ! Convergence tests
      real(8) :: ta,tdela,xtot,xtoto,enm,enb,enold,en0,edot,ye
      real(8) :: ytot,ztot,atot
      integer irdymx(1),idymx(1)
      integer :: izone,k,kstep,kts,ktsmx,knr,kout,idiag0

!  Set reaction controls not read in from control
      idiag0=idiag
      ktsmx=10
      kstep=0
      kout=0

!  Normalize initial abundances and change units if necessary
      yt=y
!     call norm(yt) 

!  Calculate the total energy of the nuclei
      call benuc(yt,enb,enm,ytot,ztot,atot)
      xtot=atot-1.0
      en0=enm
      edot=0.0

!  Start evolution 
      Write(6,"(a,i6,a,i2,2(a,es10.3))") 'Max Step',kstmx,' 
     &  IDiag=',idiag,' Start Time',tstart,' Stop Time',tstop
      t=tstart
      tt=tstart
      call t9rhofind(kstep)
      If(itso>0) call ts_output(kstep,(enm-en0),edot,kts,knr)
!-----------------------------------------------------------------------------  
!  For each step, an initial guess for the timestep is calculated, based on 
!  the abundance variations of the previous step if available.  
!-----------------------------------------------------------------------------  
      Step: Do kstep=1,kstmx
        call timestep(kstep)
        If(idiag>=1) Write(50,*) 'TDel',tt,tdel
        xtoto=xtot

!  Determine if this is an output step
        idiag=idiag0
!       If(mod(kstep,10).eq.0) idiag=2
!       If(kstep==518.or.kstep==718.or.kstep==803) idiag=5
!-----------------------------------------------------------------------------  
!  For each trial timestep, tdel, the Newton-Raphson iteration is attempted.
!  The possible results for a timestep are a converged set of yt or a 
!  failure to converge.  If convergence fails, iteration is retried with 
!  the trial timestep reduced by tdelmm, up to ktsmx times.
!-----------------------------------------------------------------------------  
        TS: Do kts=1,ktsmx

!  Calculate the thermodynamic factors necessary for reaction rates, 
!  including screening, and the reaction rates.
          call cross_sect 

!  The Newton-Raphson iteration occurs for at most knrmx iterations.  
          NR: Do knr=1,knrmx
                   
!  Calculate the changes in abundaces, dy
            call netmatr(kstep)

!  Evolve the abundances and calculate convergence tests
            yt=yt+dy
            Where(yt<ymin) 
              yt=0.0
              reldy=0.0
            ElseWhere
              reldy=abs(dy/yt)
            EndWhere
            If(idiag>=3) Then
              irdymx=maxloc(reldy)
              idymx=maxloc(dy)
              Write(50,"(a2,i5,2i3,2(a5,2es12.4))") 
     &          'dY',kstep,kts,knr,nname(idymx(1)),
     &          dy(idymx(1)),y(idymx(1)),nname(irdymx(1)),
     &          reldy(irdymx(1)),y(irdymx(1))
              If(idiag>=4) Write(50,"(a5,5es12.4)") 
     &          (nname(k),yt(k),dy(k),reldy(k),(aa(k)*dy(k)),
     &          (aa(k)*yt(k)),k=1,ny)
            Endif

!-----------------------------------------------------------------------------  
!  There are 3 included convergence tests: testc, which measures relative
!  changes, testc2 which measures total abundance changes, and testm
!  which tests mass conservation.  
!-----------------------------------------------------------------------------  
            testc=sum(reldy)
            testc2=sum(aa*dy)
            xtot=sum(aa*yt)-1.0
            testm=xtot-xtoto
            If(idiag>=2) Write(50,"(a3,i5,i3,3es14.6)") 
     &        'KNR',kstep,knr,testm,testc,testc2

!-----------------------------------------------------------------------------  
!  testc is the most stringent test, and requires the most iterations.  
!  testm is the most lax, and therefore the fastest, often requiring only one 
!  iteration.  Considering the uncertainties of the reaction rates, it is
!  doubtful that the increased precision of testc is truly increased
!  accuracy. 
!-----------------------------------------------------------------------------  
!  Ordinarily, test for true convergence
            If (iconvc/=0.or.tt>=tstop) Then
              testn=testc
              toln=tolc
  
!  Otherwise, use mass conservation for convergence condition 
            Else
              testn=testm
              toln=tolm           
            Endif
            If(abs(testn)<=toln) Exit TS
          Enddo NR

!  If convergence is not achieved in knrmx iterations, reset abundances 
!  and try again with the timestep reduced.  
          If(idiag>=1) Write(50,*) 'TS Failure',knr,kts,xtot,testn
          tdel=tdel/tdelmm
          tt=t+tdel
          call t9rhofind(kstep)
          yt=y
        Enddo TS

!  If convergence is successful, update time and abundances 
        If(kts<ktsmx) Then
          If(idiag>=1) Write(50,"(a4,i5,i3,3es12.4)") 
     &      'Conv',kstep,knr,xtot,testn,toln
          ta=t
          tdela=tdel
          t=tt
!         call norm(y)
          yo=y
          y=yt
          ye=sum(zz*y)
          If(idiag>=2) Then
            Write(50,"(a)") 'delta Y'
            Write(50,"(a5,4es12.4)") (nname(k),y(k),yo(k),
     &        (yt(k)-y(k)),(tdel*ydot(k)),k=1,ny)
          Endif
          enold=enm
          call benuc(yt,enb,enm,ytot,ztot,atot)
          edot=-(enm-enold)/tdel
          If(idiag>=1) Then 
            Write(50,"(i5,5es14.7)") kstep,t,tdel,t9t,rhot,ye
            Write(50,"(5(a5,es11.4))") (nname(k),y(k),k=1,ny)
          Endif
          If(itso>0) call ts_output(kstep,(enm-en0),edot,kts,knr)
          If(t>=tstop) Then
            Exit STEP
          Endif

!  If reduced timesteps fail to yield convergence, warn and exit
        Else
          Write(6,"(i5,4es12.4,2i3)") kstep,t,tdel,t9t,rhot,knr,kts
          Write(6,*) 'Timestep retrys fail after ',kts,' attempts'
          Exit STEP
        Endif
      Enddo STEP

!  Test that the stop time is reached
      If(t<tstop) Then
        Write(50,"(a,es12.4,a,es12.4,a)") 'Evolution incomplete!!!' 
        Write(6,"(a,es12.4,a,es12.4,a)") 'Evolution stopped at time=',t,
     &  '.  Stop time (',tstop,') not reached!' 
        Write(6,"(a,i6)") 
     &    'Approximately',int((tstop-t)/tdel),'more steps needed' 
      Endif

!  End Post Processing cycle
      call final_output(kstep)
      Return
      End subroutine full_net
 
      subroutine benuc(y,enb,enm,ytot,ztot,atot)
!===============================================================================
!  this routine finds moments of the abundance distribution useful for
!  hydrodynamics, including the total abundance, electron fraction, binding 
!  energy, and mass excess energy and outputs in and mol/g and ergs/g. 
!===============================================================================
      use constants
      use nuclear_data
      real(8), parameter :: mex_p=7.28899, mex_n=8.07144
      real(8)  :: y(ny),enb,enm,ytot,ztot,atot
      ytot=sum(y)
      ztot=sum(zz*y)
      atot=sum(aa*y)
      enb =sum(be*y)
      enm =mex_p*ztot+mex_n*(1.0-ztot)-enb

!  Change units from MeV/nucleon to erg/g
      enb=epmev*avn*enb
      enm=epmev*avn*enm
      Return
      End subroutine benuc

      subroutine t9rhofind(kstep)
!===============================================================================
!  This routine calculates t9 and rho as a function of time, either via 
!  interpolation or via from an analytic expression
!===============================================================================
      use conditions
      use thermo_data
      integer :: n,kstep
      real(8) :: dt,rdt,dt9,drho
!     t9t=t9start
!     rhot=rhostart

!  Calculate T9 and rho by interpolation 
      Do n=1,nh
        If(tt<=th(n)) Exit
      Enddo
      If(n>1.and.n<=nh) Then
        rdt=1.0/(th(n)-th(n-1))
        dt=tt-th(n-1)
        dt9=t9h(n)-t9h(n-1)
        drho=rhoh(n)-rhoh(n-1)
        t9t=dt*rdt*dt9+t9h(n-1)
        rhot=dt*rdt*drho+rhoh(n-1)
      Elseif(n==1) Then
        t9t=t9h(1)
        rhot=rhoh(1)
      Else
        t9t=t9h(nh)
        rhot=rhoh(nh)
        Write(6,*) 'Time beyond thermodynamic range',tt,' >',th(nh)
      Endif

!  Calculate T9 and rho by function
!     chi=1.0
!     thd=chi*446.0/sqrt(rhostart)
!     tint=(tstart-tt)/thd
!     rhot=rhostart*exp(tint)
!     t9t=t9start*exp(tint/3.)

!  Output T9 and rho
!     Write(50,"(a5,i5,3es12.4)") 'T9rho',kstep,tt,t9t,rhot
      Return
      End subroutine t9rhofind

      subroutine timestep(kstep)
!===============================================================================
!  This routine calculates the trial timestep, based on the prior timestep.  
!  For tdel >0, this calculation is based on the relative changes of the 
!               abundances during the previous step in the evolution.  
!           =0, the timestep is calculated using the time derivatives of 
!               the abundances based on reaction rates.  
!           <0, the timestep is held constant, tdel=-tdelo.
!  The timestep is further constrained to be no more than a factor tdelmm 
!  larger than the timestep in the previous step.  There is also the 
!  provision for limiting the timestep in the event that the thermodynamic 
!  conditions are changing too rapidly.
!===============================================================================
      use nuclear_data
      use controls
      use conditions
      use abundances
      use thermo_data
      integer, save :: ints(1),intso(1)   ! nucleus governing timestep
      integer :: i,kstep
      real(8), dimension(ny) :: ydotoy,tdt
      real(8) :: changeth,changest,tdelo,t9old,rhold,dt,dtherm
      real(8) :: tdels,tdelm,tdeln

!  Retain old values of timestep and thermo and calculate remaining time  
      changeth=.1
      changest=.01
      tdelo=tdel
      t9old=t9t
      rhold=rhot
      tdels=tstop-t
      dt=1.0e20

!  If this is not the initial timestep, calculate timestep from changes 
!  in last timestep.
      If(tdelo>0.0) Then
        tdelm=tdelmm*tdelo
        Where(y>ytime)
          ydotoy=abs((y-yo)/y)
        ElseWhere
          ydotoy=0.0
        EndWhere
        ints=maxloc(ydotoy)
        If(ydotoy(ints(1)).ne.0) Then
          tdeln=changemx*tdelo/ydotoy(ints(1))
        Else
          tdeln=tdelm
        Endif
        tdel=min(tdeln,tdels,tdelm)

!  If this is an initial timestep, yo does not exist, so calculate 
!  timestep from derivatives. If derivatives are zero, take a small step.
      Elseif(tdelo==0.0) Then
        tdelm=1.0e-4*tdels
        intso(1)=0
!       If(t9t.le.1.0) Then
!         changest=.0001*changemx
!       Elseif(t9t.le.3.0) Then
!         changest=.001*changemx
!       Else
!         changest=.01*changemx
!       Endif
        call cross_sect
        call yderiv
        Where(y>ytime)
          ydotoy=abs(ydot/y)
        ElseWhere
          ydotoy=0.0
        EndWhere
        ints=maxloc(ydotoy)
        If(ydotoy(ints(1)).ne.0) Then
          tdeln=changest*changemx/ydotoy(ints(1))
        Else
          tdeln=tdelm
        Endif
        tdel=min(tdels,tdeln,tdelm)

!  Keep timestep constant
      Else
        tdel=-tdelo
      Endif

!  Diagnostic Output
      If(idiag>=1) Write(50,"(a4,i5,4es12.4,i5)") 
     &  'tdel',kstep,tdel,tdeln,tdelo,tdels,ints(1)
!     If(idiag>=2) Then
!       Write(50,"(a5,i4,2es12.4)") 
!    &    (nname(k),k,y(k),ydotoy(k),k=1,ny)
!     Endif

!  Retain the index of the species setting the timestep 
      If(ints(1)/=intso(1)) Then
        If(idiag>=1) Write(50,*) 'ITC ',nname(ints(1)),t,tdel
        intso=ints
      Endif

!  Limit timestep if Thermodynamic variation is too large
      Do i=1,10
        tt=t+tdel
        call t9rhofind(kstep)
        If(t9old>0) Then
          dtherm=abs(t9t-t9old)/t9old+abs(rhot-rhold)/rhold
        Endif
        If(dtherm<changeth) Exit
        tdel=.5*tdel
        If(i==10) Write(6,*) 'Error in Thermo variations after ',i,
     &    'reductions',tdel,t9t,rhot
      Enddo

!     If(idiag>=1) Write(50,"(a5,i5,2es12.4)") 
!    &  'T9del',kstep,tdel,dtherm
      Return
      End subroutine timestep                                                                 

      subroutine en(n,y,ye,yps,y2e)  
!===============================================================================
!  This routine calculates moments of the abundance distribution needed 
!  for the EOS.
!===============================================================================
      use nuclear_data
      integer :: n
      real(8) :: ye,yps,y2e
      real(8), dimension(n) :: y
      ye=sum(zz*y)
      y2e=sum(zz**2*y)
      yps=sum(y)+ye 
!     Write(50,"(a3,3es16.8)") 'EN',ye,y2e,yps
      Return                                                                    
      End subroutine en                                                                      

      subroutine netmatr(kstep)
!===============================================================================
!  This routine calculates the Jacobian matrix dYdot/dY, and solves for the 
!  Newton-Raphson iteration, dy.
!===============================================================================
      use controls
      use nuclear_data
      use conditions
      use abundances
      use reac_rate_data
      integer, dimension(ny) :: indx
      integer :: i,j,kstep,i0,la1,le1,la2,le2,la3,le3,j1,l1,l2,l3
      real(8), dimension(ny) :: f,ydot0,um
      real(8), dimension(ny,ny) :: am  !,at ! the Jacobian Matrix
      real(8) :: alph,beta,rdt,d

      alph=0.0                                                                  
      beta=1.0                                                                  
      ydot0=0.0

! Calculate the reaction rates and abundance time derivatives
      call yderiv

!  Build the Jacobian, row by row
      rdt=1.0/tdel/beta 
      Do i0=1,ny 
        um=0.0              
        um(i0)=rdt  
        la1=la(1,i0) 
        le1=le(1,i0)  
        Do j1=la1,le1  
          l1=n11(j1)
          um(l1)=um(l1)-b1(j1)
        Enddo 
        la2=la(2,i0) 
        le2=le(2,i0)  
        Do j1=la2,le2  
          l1=n21(j1)
          l2=n22(j1) 
          um(l1)=um(l1)-b2(j1)*yt(l2) 
          um(l2)=um(l2)-b2(j1)*yt(l1)  
        Enddo       
        la3=la(3,i0) 
        le3=le(3,i0)  
        Do j1=la3,le3 
          l1=n31(j1) 
          l2=n32(j1) 
          l3=n33(j1)
          um(l1)=um(l1)-b3(j1)*yt(l2)*yt(l3)    
          um(l2)=um(l2)-b3(j1)*yt(l1)*yt(l3)     
          um(l3)=um(l3)-b3(j1)*yt(l1)*yt(l2)      
        Enddo                 

!  Tranfer to matrix row
        am(i0,:)=um
!       at(:,i0)=um ! or column if the solver wants the transpose
      Enddo                                                      
!     am=transpose(at)

!  Calculate equation to zero
      f=(y-yt)*rdt+ydot
!     f=y*rdt-yt*rdt+ydot+alph*ydot0/beta  
      If(idiag>=4) Then
        Write(50,"(a2,i5,es14.7)") 'F',kstep,rdt
        Do i=1,ny
          Write(50,"(a5,4es17.9)") nname(i),f(i),ydot(i),yt(i),y(i)
          Write(50,"(5es16.8)") (am(i,j),j=1,ny)
        Enddo
      Endif

!  Test the eigenvalues 
!     If(idiag>=6) Then
!       call eigen_test(kstep,am,rdt)
!     Endif
!-----------------------------------------------------------------------------  
!  The bulk of the computational cost of the network (60-95%) is the solving 
!  of the matrix equation.  Careful selection of the matrix solver is therefore 
!  very important to fast computation.  Generally, hand tuned solvers such as 
!  those supplied by the hardware manufacturer or third-parties like NAG, IMSL,
!  etc. are the fastest.  However for portability, by default we use Numerical 
!  Recipes routines.  
!-----------------------------------------------------------------------------  

!  Use Num Rec LU Decomp.
!     call ludcmp(am,ny,ny,indx,d)
!     call lubksb(am,ny,ny,indx,f)
!     dy=f

!  Use LAPACK solver
!     call sgesv(ny,1,am,ny,indx,f,ny,info) ! Single precision version
      call dgesv(ny,1,am,ny,indx,f,ny,info) ! Double precision version
      dy=f

!  Diagnostic output
!     If(idiag>=4) Then
!       Write(50,"(a5,3es12.4)") (nname(i),dy(i),y(i),yt(i),i=1,ny)
!     Endif
      Return                                                                    
      End subroutine netmatr                                                                       

      subroutine partf(t9) 
!===============================================================================
!  This routine calculates the nuclear partition functions as a function of T9
!===============================================================================
      use part_funct_data
      integer :: i,ii
      real(8) :: t9,rdt9
      Do i=1,24                                                              
        If(t9<=t9i(i)) Exit
      Enddo                                                                     
      ii=i                                                                      
      Select Case (ii)
        Case(2:24)
          rdt9=(t9-t9i(ii-1))/(t9i(ii)-t9i(ii-1))
          gg(1:ny)=exp(rdt9*log(g(ii,1:ny))+(1-rdt9)*log(g(ii-1,1:ny)))
        Case (1)
          gg(1:ny)=g(1,1:ny) 
        Case (25)
          gg(1:ny)=g(24,1:ny)  
      End Select                                                                

!  gg(0) is a placeholder for non-nuclei, gamma-rays, etc. 
      gg(0)=1.0
!     Write(50,"(a5,i3,es14.7)") 'PartF',ii,t9
!     Write(50,"(5(i4,es12.4))") (i, gg(i), i=1,ny) 
      Return                                                                    
      End subroutine partf                                                                      

      subroutine cross_sect
!===============================================================================
!  This routine calculates the cross section for each reaction.
!===============================================================================
      use controls
      use nuclear_data
      use abundances
      use conditions
      use cross_sect_data
      use part_funct_data 
      use ffn_data
      use reac_rate_data
      real(8) :: t09(7)
      real(8) :: ene,ye,yps,y2e,v,pa,ea,dlma,pta,pva,eta,eva,emas,beta
      real(8) :: funct,rat,rhot2,zz1,zz2,zz3,zz12,aa12,h01,h02,fsr,fst
      integer :: j,k,mu,ieos
      emas=4.
      ieos=0
      eta=1.e+30                                                                

!  Call EOS, necessary for screening
      call en(ny,yt,ye,yps,y2e)
      If(iscrn>=1) Then
        V=1./rhot                                                       
        call state(t9t,v,y,pa,ea,dlma,ye,pta,pva,eta,eva,ieos,ny,
     &    emas,yps,beta)
      Endif                                                     
      
!     Write(50,"('CS2:',i2,7es10.3)") (i,rc2(:,i),i=1,nreac(2))
	    
!     Write(50,"('CS1:',i2,7es10.3)") 2,rc1(:,2)
!     Write(50,"('CS1:',i2,7es10.3)") 8,rc1(:,8)
!     Write(50,"('CS2:',i2,7es10.3)") 4,rc2(:,4)
!     Write(50,"('CS2:',i2,7es10.3)") 8,rc2(:,8)
!     Write(50,"('CS2:',i2,7es10.3)") 9,rc2(:,9)

!  Calculate necessary thermodynamic moments
      ene=ye*rhot                                                               
      t09(1)=1.                                                                 
      t09(2)=t9t**(-1) 
      t09(3)=t9t**(-1./3.) 
      t09(4)=t9t**(+1./3.) 
      t09(5)=t9t          
      t09(6)=t9t**(5./3.)
      t09(7)=log(t9t)  
!     Write(50,"(a3,4es12.4)") 'THR',t9t,rhot,ene,ye

!  Calculate partition functions for each nucleus at t9t
      call partf(t9t) 

!  If there are any FFN reactions, calculate their rates
      If(nffn>0) Then
        call ffnrate(t9t,ene)   
      Endif

!  Calculate csects for reactions with one reactant
      funct=1.0
      rat=1.0
      Where(irev1==1)                   ! If it's a reverse reaction
        rpf1=gg(n1i(2,:))*gg(n1i(3,:))/gg(n1i(1,:))
      ElseWhere
        rpf1=1.0
      EndWhere
      If(iweak>0) Then             ! All Rates on
        Where(iwk1/=2.and.iwk1/=3)        ! If it's not an FFN reaction
          csect1=rpf1*exp(matmul(t09,rc1))
        ElseWhere
          iffn=int(rc1(1,:))
          csect1=rat*funct*rf(iffn)
        EndWhere
        Where (iwk1==1) csect1=ene*csect1 ! If it's a non-FFN EC reaction
      ElseIf(iweak<0) Then         ! Only weak rates used
        Where(iwk1==0) csect1=0.0
        Where(iwk1==1) csect1=ene*rpf1*exp(matmul(t09,rc1))
        Where(iwk1==2.or.iwk1==3) 
          iffn=int(rc1(1,:))
          csect1=rat*funct*rf(iffn)
        EndWhere
        Where(iwk1>=4) csect1=rpf1*exp(matmul(t09,rc1))    
      Else                         ! weak interactions are off (iweak==0)
        Where(iwk1==0) 
          csect1=rpf1*exp(matmul(t09,rc1))  
        ElseWhere
          csect1=0.0
        EndWhere
      Endif
      If(idiag>=5) Then
        Write(50,"(a,i5)") 'CSect1',nreac(1)
        Write(50,"(i5,5a5,3i3,es17.9)") 
     &    (k,nname(n1i(1,k)),'-->',(nname(n1i(j,k)),j=2,3),'+++',
     &    iwk1(k),irev1(k),ires1(k),csect1(k), k=1,nreac(1))
      Endif

!  Calculate screening corrections for 2 reactant reactions
      If(iscrn>=1) call screen2

!  Calculate the csect for reactions with 2 reactants
      Where(irev2==1)   ! If it's a reverse reaction mult. by ratio of part. func.
        rpf2=gg(n2i(3,:))*gg(n2i(4,:))/(gg(n2i(1,:))*gg(n2i(2,:)))
      ElseWhere
        rpf2=1.0
      EndWhere
      If(iweak>0) Then        ! All rates are on
        csect2=rhot*rpf2*exp(matmul(t09,rc2)+h2)
        Where(iwk2==1) csect2=ene*csect2  ! If it's a non FFN EC reaction
      ElseIf(iweak<0) Then    ! Only weak rates
        Where(iwk2==0) 
          csect2=0.0
        ElseWhere
          csect2=rhot*rpf2*exp(matmul(t09,rc2)+h2)
        EndWhere
        Where(iwk2==1) csect2=ene*csect2  
      Else                    ! Weak interactions are off  (iweak=0)
        Where(iwk2==0)
          csect2=rhot*rpf2*exp(matmul(t09,rc2)+h2)           
        ElseWhere
          csect2=0.0
        EndWhere
      Endif
      If(idiag>=5) Then
        Write(50,"(a,i5)") 'CSect2',nreac(2)
        Write(50,"(i5,5a5,3i3,es17.9)") 
     &    (k,(nname(n2i(j,k)),j=1,2),'-->',(nname(n2i(j,k)),j=3,4),
     &    iwk2(k),irev2(k),ires2(k),csect2(k),k=1,nreac(2))
      Endif
      If(idiag>=5) Write(50,"(a,i5)") 'CSect3',nreac(3)
      rhot2=rhot**2
      Do mu=1,nreac(3)
!  Compute screening corrections
        If(iscrn>=1) Then
          zz1=zz(n3i(1,mu))                                 
          zz2=zz(n3i(2,mu))                                
          zz3=zz(n3i(3,mu))                               
          zz12=zz1+zz2                              
          aa12=aa(n3i(1,mu))+aa(n3i(2,mu))                   
          If(zz1*zz2.eq.0.0) Then
            h01=0.0
          Else
            call screen(zz1,zz2,aa(n3i(1,mu)),aa(n3i(2,mu)),h01,fsr,fst) 
          Endif
          If(zz12*zz3.eq.0.0) Then                          
            h02=0.0    
          Else
            call screen(zz12,zz3,aa12,aa(n3i(3,mu)),h02,fsr,fst)  
          Endif
        Else
          h01=0.0 ; h02=0.0
        Endif

!  Compute the crossection
        If(iweak>0) Then
          csect3(mu)=t09(1)*rc3(1,mu)+t09(2)*rc3(2,mu)+
     &      t09(3)*rc3(3,mu)+t09(4)*rc3(4,mu)+t09(5)*rc3(5,mu)+
     &      t09(6)*rc3(6,mu)+t09(7)*rc3(7,mu)
          csect3(mu)=rhot2*exp(csect3(mu)+h01+h02)
          If(iwk3(mu).eq.1) csect3(mu)=csect3(mu)*ene        
          If(csect3(mu).lt.1.e-20) csect3(mu)=0.0                   
        ElseIf(iweak<0) Then  ! Only Weak Reactions on
          If(iwk3(mu)==0) csect3(mu)=0. 
        Else                  ! Weak Reactions off (iweak=0)
          If(iwk3(mu)==0) Then
            csect3(mu)=t09(1)*rc3(1,mu)+t09(2)*rc3(2,mu)+
     &        t09(3)*rc3(3,mu)+t09(4)*rc3(4,mu)+t09(5)*rc3(5,mu)+
     &        t09(6)*rc3(6,mu)+t09(7)*rc3(7,mu)
            csect3(mu)=rhot2*exp(csect3(mu)+h01+h02)
            If(csect3(mu).lt.1.e-20) csect3(mu)=0.0 
          Else
            csect3(mu)=0.
          Endif  
        Endif
        If(idiag>=5) Write(50,"(i5,5a5,3i3,es17.9)") 
     &    mu,(nname(n3i(j,mu)),j=1,3),'-->','+++',
     &    iwk3(mu),irev3(mu),ires3(mu),csect3(mu)
      Enddo    
      Return                                  
      End subroutine cross_sect                                   

      subroutine yderiv
!-----------------------------------------------------------------------------  
!  This routine calculates time derivatives for each nuclear species.
!  This calculation is performed by looping over nuclei, and summing the 
!  reaction rates for each reaction which effects that nucleus.
!-----------------------------------------------------------------------------  
      use controls
      use nuclear_data
      use abundances
      use cross_sect_data
      use reac_rate_data
      integer :: i,j,i0,i1,la1,le1,la2,le2,la3,le3
      real(8) :: s1,s11,s2,s22,s3,s33

!  From the cross sections and the counting array, calculate the reaction 
!  rates
      b1=a1*csect1(mu1)
      b2=a2*csect2(mu2)
      b3=a3*csect3(mu3)

!  Calculate Ydot for each nucleus, summing over the reactions which affect it.
      Do i0=1,ny                 
        If(idiag>=5) Write(50,"(a3,a6,i4)") 'NUC',nname(i0),i0

!  Sum over the reactions with 1 reactant
        la1=la(1,i0)                          
        le1=le(1,i0)                         
        s1=0.0                              
        Do i1=la1,le1                    
          s11=b1(i1)*yt(n11(i1))
          s1=s1+s11
          If(idiag>=5) Write(50,"(3a5,'  1  ',4es15.7)") 
     &      (nname(n1i(j,mu1(i1))),j=1,3),b1(i1),yt(n11(i1)),
     &      s11,a1(i1)
        Enddo       
        If(idiag>=5) Write(50,*) '1->',nname(i0),la1,le1,s1

!  Sum over the reactions with 2 reactants
        la2=la(2,i0)  
        le2=le(2,i0) 
        s2=0.0      
        Do i1=la2,le2           
          s22=b2(i1)*yt(n21(i1))*yt(n22(i1))
          s2=s2+s22                                       
          If(idiag>=5) Write(50,"(4a5,4es15.7)")
     &      (nname(n2i(i,mu2(i1))),i=1,4),b2(i1),yt(n21(i1)),
     &      yt(n22(i1)),s22
        Enddo                  
        If(idiag>=5) Write(50,*) '2->',nname(i0),la2,le2,s2

!  Sum over the reactions with 3 reactants
        la3=la(3,i0)                          
        le3=le(3,i0)                         
        s3=0.0                              
        Do i1=la3,le3                  
          s33=b3(i1)*yt(n31(i1))*yt(n32(i1))*yt(n33(i1))
          s3=s3+s33
          If(idiag>=5) Write(50,"(3a5,'  3  ',5es12.4)") 
     &      (nname(n3i(i,mu3(i1))),i=1,3),b3(i1),yt(n31(i1)),
     &      yt(n32(i1)),yt(n33(i1)),s33
        Enddo                                       
        If(idiag>=5) Write(50,*) '3->',nname(i0),la3,le3,s3

!  Sum the 3 components of Ydot
        ydot(i0)=s1+s2+s3                         
        If(idiag>=5) Write(50,"(a4,a5,2es24.16)") 
     &    'YDOT',nname(i0),yt(i0),ydot(i0)
      Enddo                                 
      Return                                  
      End subroutine yderiv                                   

      subroutine ffnrate(t9,ene) 
!-----------------------------------------------------------------------------  
!  This routine calculates the reaction rates for the weak rates of 
!  Fuller, Fowler, & Newman
!-----------------------------------------------------------------------------  
      use ffn_data
      real(8) :: t9,ene,tg(13),egl(11),enl,dt,de,ddt
      integer i1(4)
      integer i,le1,lt1
      data tg/0.01,0.1,0.2,0.4,0.7,1.0,1.5,2.,3.,5.,10.,30.,100./               
      data egl/1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11./                              
      Do i=1,13
        If(t9<=tg(i)) Exit
      Enddo
      lt1=i-1
      If(lt1<1) lt1=1
      enl=log10(ene)                                                           
      le1=int(enl)
      If(le1<1) le1=1
      dt=tg(lt1+1)-tg(lt1)                                                      
      de=egl(le1+1)-egl(le1)                                                    
      ddt=t9-tg(lt1)                                                            
      i1(1)=13*(le1-1)+lt1                                                      
      i1(2)=i1(1)+1                                                             
      i1(3)=13*le1+lt1                                                          
      i1(4)=i1(3)+1                                                             
      Do i=1,nffn                                                              
        r1(i)=ffnsum(i,i1(1))+(ffnsum(i,i1(2))-ffnsum(i,i1(1)))/dt*ddt
        r2(i)=ffnsum(i,i1(3))+(ffnsum(i,i1(4))-ffnsum(i,i1(3)))/dt*ddt
        rf(i)=r1(i)+(r2(i)-r1(i))/de*(enl-egl(le1))    
        If(rf(i).lt.-30.) Then
          rf(i)=0.
        Else
          rf(i)=10.**rf(i)     
        Endif
      Enddo                                                                     
      Return                                                                    
      End subroutine ffnrate

      subroutine screen2
      use controls
      use nuclear_data
      use cross_sect_data
      integer :: mu
      real(8) :: z1,z2,z12,a1,a2,z,zl,dzldt,dzldr
      real(8) :: gsi,gsiz,gsb,f,fsr,fst,gnp,gnz,em2,tau,t9er,gt3,f90
      real(8) :: hw,hi,hs,dhdz
      common /scrn/ gnp,gsi,zl,dzldr,dzldt,t9er                         
!     Write(20,"(a3,3es15.7)") 'EOS',GNP,GSI,ZL
      Do mu=1,nreac(2)                                                      
        z1=zz(n2i(1,mu))                                                  
        z2=zz(n2i(2,mu))                                                  
        a1=aa(n2i(1,mu))                                                  
        a2=aa(n2i(2,mu))                                                  
        z12=z1*z2
        If(z12==0.0.or.iscrn==0) Then
          h2(mu)=0.0                                                        
        Else
  
!  Weak and Intermediate Screening Graboske et al.         
          z=zl*z12                                                           
          gsiz=gsi*z12                                                       
          gsb=gsiz**0.86                                                    
          If (gsiz>1.0) Then
            f=0.38*((1.+gsiz)**1.86-gsb*gsiz-1.)/gsb          
          Else
            f=0.61943*gsiz**0.07                              
          Endif
          hw=z*(1.+z*(log(z)+0.8364))                                      
          hi=f*z**0.86                                                      
          If (hw<=hi) Then
            h2(mu)=hw         
            dhdz=z*(1.+2.*z*(log(z)+1.3364))                
          Else
            h2(mu)=hi                             
            dhdz=0.86*hi                                   
          Endif
          fsr=dhdz*dzldr                                                    
          fst=dhdz*dzldt                                                    

!  Strong screening by Itoh et al.(1990)                             
          gnz=gnp*z12*2.0/(z1**(1./3.)+z2**(1./3.))
          If (gnz>=0.4) Then
            em2=a1*a2*2.0/(a1+a2)                
            tau=3.3722*(em2*z12**2*t9er)**(1./3.) 
            gt3=3.0*gnz/tau
            f90=(.0455*gt3+.348*gt3**3+9.49*gt3**6-.123*gt3**12+
     &          .101*gt3**13)/(1.+100.*gt3**4+.267*gt3**12)
            hs=1.25*gnz-tau*f90 
!           hs=gnz*(1.25-0.285*gt3)    
            If (hs<h2(mu)) Then
                h2(mu)=hs               
                fsr=gnz*(1.25-.57*gt3)/3. 
                fst=-gnz*(1.25-.475*gt3)   
            Endif
          Endif
!         Write(20,"(3a5,i6,4es11.4)") 'SCR2',nname(n2i(1,mu)),
!    &         nname(n2i(2,mu)),dexp(h0(mu)),hw,hi,hs
        Endif
      Enddo                                                             
!     Write(20,*) (nname(n2i(1,mu)),nname(n2i(2,mu)),h2(mu),mu=1,n)
      Return                                                            
      End subroutine screen2                                                              
      
      subroutine norm(yy)
!---------------------------------------------------------------------------- 
!  This routine renormalizes the abundances to guarantee mass conservation.
!---------------------------------------------------------------------------- 
      use nuclear_data
      real(8) :: xtot,rxt
      real(8), dimension(ny) :: yy
      xtot=sum(yy*aa)   
      rxt=1.0/xtot      
      yy=yy*rxt  
      Return                                               
      End subroutine norm                                                

      subroutine ye_norm(yy,ye)             
!-----------------------------------------------------------------------------  
! This routine renormalizes the abundances to guarantee mass and charge
! conservation if Ye is known.
!-----------------------------------------------------------------------------  
      use nuclear_data
      real(8) :: ye,zy,nny,zny,zzy,beta,alph
      real(8) :: yy(ny)
      nny=sum(nn*yy)
      zy=sum(zz*yy)
      zny=sum(nn*zz*yy/aa)
      zzy=sum(zz*zz*yy/aa)
      beta=(ye*nny-zny)/(nny*zzy-zy*zny)
      alph=(1-beta*zy)/nny
      yy=yy*(alph*nn+beta*zz)/aa
      Return
      End subroutine ye_norm

