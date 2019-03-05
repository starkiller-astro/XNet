!*******************************************************************************
!  Common.f 11/24/09
!  This subroutine contains all those pices which are common to the various ways
!  of running xnet.                                                           
!                           
!
!  The data.f file contains the nuclear and reaction data structures and the 
!  subroutines which read in the data and allocate the arrays.  
!*******************************************************************************

      module controls
!===============================================================================
!  This module contains the values of the flags and limits which control 
!  the behavior of the network
!===============================================================================
      integer :: iweak     ! If >0, strong and weak reaction are used
                           ! If =0, weak reaction are ignored
                           ! If <0, only weak reactions are used
      integer :: iweak0    ! Saves input iweak flag
      integer :: iscrn     ! If =0, screening is ignored
      integer :: itso      ! Sets level of time series output
      integer :: idiag     ! Sets level of diagnostic output
      integer :: iconvc    ! Controls type of convergence condition (0=mass)
      integer :: kstmx     ! Max # of timesteps before exit
      integer :: knrmx     ! Max # of Newton-Raphson iterations
      integer :: nzone     ! Number of zones
      integer :: kts       ! Current number of timestep iterations
      integer :: knr       ! Current number of Newton-Raphson iterations
      real(8) :: tolc      ! The iterative convergence test limit
      real(8) :: tolm      ! Max network mass error
      real(8) :: changemx  ! Relative abundance change used to guess timestep
      real(8) :: ytime     ! Min abundance included in timestep estimate
      real(8) :: ymin      ! Min abundance, y < ymin =0
      real(8) :: tdelmm    ! new timestep <= tdelmm * previous timestep
      real(8) :: t9cut=0.01! Temperature cut for turning strong reactions off
      integer :: inout(14) ! List of species to output in condensed form
      integer :: myid, nproc, mythread, nthread ! task & thread ids and counts
      integer :: lun_diag, lun_ev, lun_ts ! logical units
      common /flag/ iweak,iscrn,itso,idiag,iconvc,kstmx,knrmx,nzone
      common /tol/ changemx,tolm,tolc,ytime,ymin,tdelmm

!$OMP THREADPRIVATE(mythread,nthread,lun_diag,lun_ev,lun_ts)

      end module controls

      module conditions
!===============================================================================
!  This module contains data on the current time and thermodynamic conditions.
!===============================================================================
      real(8) :: t        ! Time at the beginning of the current timestep
      real(8) :: tt       ! Trial time for end of current timestep
      real(8) :: tdel     ! Trial duration of timestep
      real(8) :: t9t,rhot ! Temperature (GK) and Density (g/cc) at trial time

! Threading Scope
!$OMP THREADPRIVATE(t,tt,tdel,t9t,rhot)


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

! Threading Scope
!$OMP THREADPRIVATE(yo,y,yt,ydot,dy)

      end module abundances

      module thermo_data
!===============================================================================
!  This module contains the thermodynamic trajectory which the network follows
!===============================================================================
      integer, parameter   :: nhmx=5000 ! The max number of thermo points
      integer              :: nh        ! The actual number of points
      real(8) :: tstart,tstop,th(nhmx),t9h(nhmx),rhoh(nhmx)

! Threading Scope
!$OMP THREADPRIVATE(nh,tstart,tstop,th,t9h,rhoh)

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

      subroutine t9rhofind(kstep,tint,nn)
!===============================================================================
!  This routine calculates t9 and rho as a function of time, either via 
!  interpolation or via from an analytic expression
!===============================================================================
      use conditions
      use thermo_data
      use controls
      
      integer :: n,nn,kstep
      real*8  :: tint
      real(8) :: dt,rdt,dt9,drho
!     t9t=t9start
!     rhot=rhostart

!  Calculate T9 and rho by interpolation 
      Do n=1,nh
        If(tint<=th(n)) Exit
      Enddo
      nn = n
      If(n>1.and.n<=nh) Then
        rdt=1.0/(th(n)-th(n-1))
        dt=tint-th(n-1)
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
!     Write(lun_diag,"(a5,i5,3es12.4)") 'T9rho',kstep,tt,t9t,rhot

! Turn off strong interactions if the temperature is less than t9cut
      If (t9t<t9cut) Then

        If (iweak /= 0) Then
          iweak = -1
        Else
          iweak=iweak0
          Write(6,*)'Warning: Strong reactions ignored for T9< T9cut.'
        Endif
      Endif
      
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
      integer :: i,kstep,nnew,nold
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
      If(idiag>=1) Write(lun_diag,"(a4,i5,4es12.4,i5)") 
     &  'tdel',kstep,tdel,tdeln,tdelo,tdels,ints(1)
!     If(idiag>=2) Then
!       Write(lun_diag,"(a5,i4,2es12.4)") 
!    &    (nname(k),k,y(k),ydotoy(k),k=1,ny)
!     Endif

!  Retain the index of the species setting the timestep 
      If(ints(1)/=intso(1)) Then
        If(idiag>=1) Write(lun_diag,*) 'ITC ',nname(ints(1)),t,tdel
        intso=ints
      Endif

! Make sure to not skip any features in the temperature or density profiles by checking
! for profile monotonicity between t and t+del      
      call t9rhofind(kstep,t,nold)
      tt = t+tdel
      call t9rhofind(kstep,tt,nnew)
      if (nnew-1>nold) then
        do j=nold,nnew-1
          if (t9h(j)>t9old.and.t9h(j)>t9t) then
            tdel = th(j) - t
            exit
          elseif (t9h(j)<t9old.and.t9h(j)<t9t) then
            tdel = th(j) - t
            exit
          elseif (rhoh(j)>rhold.and.rhoh(j)>rhot) then
            tdel = th(j) - t
            exit
          elseif (rhoh(j)<rhold.and.rhoh(j)<rhot) then
            tdel = th(j) - t
            exit    
          endif
        enddo
      endif  
      
!  Limit timestep if Thermodynamic variation is too large
      Do i=1,10
        tt=t+tdel
! reset the iweak flag to its original value
        iweak = iweak0
        call t9rhofind(kstep,tt,n)
        If(t9old>0) Then
          dtherm=abs(t9t-t9old)/t9old+abs(rhot-rhold)/rhold
        Endif
        If(dtherm<changeth) Exit
        tdel=.5*tdel
        If(i==10) Write(6,*) 'Error in Thermo variations after ',i,
     &    'reductions',tdel,t9t,rhot
      Enddo

!     If(idiag>=1) Write(lun_diag,"(a5,i5,2es12.4)") 
!    &  'T9del',kstep,tdel,dtherm
      Return
      End subroutine timestep                                                                 

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
!     Write(lun_diag,"(a5,i3,es14.7)") 'PartF',ii,t9
!     Write(lun_diag,"(5(i4,es12.4))") (i, gg(i), i=1,ny) 
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
      
!     Write(lun_diag,"('CS2:',i2,7es10.3)") (i,rc2(:,i),i=1,nreac(2))
	    
!     Write(lun_diag,"('CS1:',i2,7es10.3)") 2,rc1(:,2)
!     Write(lun_diag,"('CS1:',i2,7es10.3)") 8,rc1(:,8)
!     Write(lun_diag,"('CS2:',i2,7es10.3)") 4,rc2(:,4)
!     Write(lun_diag,"('CS2:',i2,7es10.3)") 8,rc2(:,8)
!     Write(lun_diag,"('CS2:',i2,7es10.3)") 9,rc2(:,9)

!  Calculate necessary thermodynamic moments
      ene=ye*rhot                                                               
      t09(1)=1.                                                                 
      t09(2)=t9t**(-1) 
      t09(3)=t9t**(-1./3.) 
      t09(4)=t9t**(+1./3.) 
      t09(5)=t9t          
      t09(6)=t9t**(5./3.)
      t09(7)=log(t9t)  
!     Write(lun_diag,"(a3,4es12.4)") 'THR',t9t,rhot,ene,ye

!  Calculate partition functions for each nucleus at t9t
      call partf(t9t) 

!  If there are any FFN reactions, calculate their rates
      If(nffn>0) Then
        call ffn_rate(t9t,ene)   
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
        Write(lun_diag,"(a,i5)") 'CSect1',nreac(1)
        Write(lun_diag,"(i5,5a5,3i3,es17.9)") 
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
        Write(lun_diag,"(a,i5)") 'CSect2',nreac(2)
        Write(lun_diag,"(i5,5a5,3i3,es17.9)") 
     &    (k,(nname(n2i(j,k)),j=1,2),'-->',(nname(n2i(j,k)),j=3,4),
     &    iwk2(k),irev2(k),ires2(k),csect2(k),k=1,nreac(2))
      Endif
      If(idiag>=5) Write(lun_diag,"(a,i5)") 'CSect3',nreac(3)
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
        If(idiag>=5) Write(lun_diag,"(i5,5a5,3i3,es17.9)") 
     &    mu,(nname(n3i(j,mu)),j=1,3),'-->','+++',
     &    iwk3(mu),irev3(mu),ires3(mu),csect3(mu)
      Enddo    
      Return                                  
      End subroutine cross_sect                                   

      
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



