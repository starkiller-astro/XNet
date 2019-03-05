!*******************************************************************************
! Common.f90 5/28/10
! This file contains modules and subroutines are common between the various 
! xnet targets, except for those associated with the nuclear data.  The data.f 
! file contains the nuclear and reaction data structures and the subroutines 
! which read in the data and set up the data structures.  
!*******************************************************************************
  
Module controls
!===============================================================================
! This module contains the values of the flags and limits which control 
! the behavior of the network
!===============================================================================
  Integer :: iweak        ! If >0, strong and weak reaction are Used
                          ! If =0, weak reaction are ignored
                          ! If <0, only weak reactions are Used
  Integer :: iweak0       ! Saves input iweak flag
  Integer :: iscrn        ! If =0, screening is ignored
  Integer :: itsout       ! Sets level of time series output
  Integer :: idiag        ! Sets level of diagnostic output
  Integer :: isolv        ! Sets the integration method (1=BE, 2=BD,3=qse)
  Integer :: iadpt        ! Turns adaptive mode on and off(1=on 0=off)         
  Integer :: iconvc       ! Controls type of convergence condition (0=mass)
  Integer :: kstmx        ! Max # of timesteps before exit
  Integer :: kitmx        ! Max # of iterations or substeps within a timestep
  Integer :: nzone        ! Number of zones
  Integer :: kmon(2)      ! Solver-dependent behaviour monitors
  Real(8) :: tolc         ! The iterative convergence test limit
  Real(8) :: tolm         ! Max network mass error
  Real(8) :: changemx     ! Relative abundance change Used to guess timestep
  Real(8) :: yacc         ! Min abundance required to be accuracte, 
                          ! Used in timestep determination.
  Real(8) :: ymin         ! Minimum abundance, y < ymin =0
  Real(8) :: tdel_maxmult ! new timestep <= tdel_maxmult * previous timestep
  Real(8) :: t9cut=0.01   ! Temperature cut for turning strong reactions off
  Integer :: inout(14)    ! List of species to output in condensed form
  Integer :: myid, nproc, mythread, nthread ! task & thread ids and counts
  Integer :: lun_diag, lun_ev, lun_ts ! logical units
! common /flag/ iweak,iscrn,itsout,idiag,iconvc,kstmx,kitmx,nzone
! common /tol/ changemx,tolm,tolc,yacc,ymin,tdel_maxmult
  
! Threading Scope
!$OMP THREADPRIVATE(mythread,nthread,lun_diag,lun_ev,lun_ts,kmon)
  
End Module controls
  
Module conditions
!===============================================================================
! This module contains data on the current time and thermodynamic conditions.
!===============================================================================
  Real(8) :: t         ! Time at the beginning of the current timestep
  Real(8) :: tt        ! Trial time for End of current timestep
  Real(8) :: to        ! Time at the beginning of the previous timestep
  Real(8) :: tdel      ! Trial duration of timestep
  Real(8) :: tdel_next ! Integrator estimate of duration for next timestep
  Real(8) :: t9t,rhot  ! Temperature (GK) and Density (g/cc) at trial time
  
! Threading Scope
!$OMP THREADPRIVATE(t,tt,tdel,tdel_next,t9t,rhot)
  
End Module conditions
  
Module abundances
!===============================================================================
! This module contains the abundances of the nuclear species, at the previous 
! time (yo), current time (y), and trial time (yt), as well as the time 
! derivatives (ydot), at the trial time, and the abundance change due to the 
! Newton-Raphson iteration.  The size of these arrays is Allocated in the 
! external routine, where the initial values are set.
!===============================================================================
  Use nuc_number
  Real(8), Dimension(:), Allocatable :: yo,y,yt,ydot
  
! Threading Scope
!$OMP THREADPRIVATE(yo,y,yt,ydot,dy)
  
End Module abundances
  
Module thermo_data
!===============================================================================
! This module contains the thermodynamic trajectory which the network follows
!===============================================================================
  Integer, Parameter   :: nhmx=5000 ! The max number of thermo points
  Integer              :: nh        ! The actual number of points
  Real(8) :: tstart,tstop,th(nhmx),t9h(nhmx),rhoh(nhmx)
  
! Threading Scope
!$OMP THREADPRIVATE(nh,tstart,tstop,th,t9h,rhoh)
  
End Module thermo_data
  
Module constants
!===============================================================================
! These are fundamental constants in CGS, except energies in MeV, Temp in GK
!===============================================================================
  Real(8), Parameter :: pi   =3.1415926536, hbar=6.582122E-22
  Real(8), Parameter :: amu  =1.036427E-18, bok =8.617385E-02
  Real(8), Parameter :: avn  =6.022137E+23, e2  =1.439830E-13
  Real(8), Parameter :: epmev=1.602177E-06
End Module constants

      module screen_qse
!-----------------------------------------------------------------------------
!  This module contains the corrections to the binding energies due to
!  nuclear screening (AKA Coulomb corrections).
!-----------------------------------------------------------------------------
      integer :: zmax
      integer, dimension(:), allocatable :: intz
      real*8, dimension(:), allocatable :: he
      end module screen_qse
  
   
Subroutine benuc(y,enb,enm,ytot,ztot,atot)
!===============================================================================
! This routine finds moments of the abundance distribution useful for
! hydrodynamics, including the total abundance, electron fraction, binding 
! energy, and mass excess energy and outputs in and mol/g and ergs/g. 
!===============================================================================
  Use constants
  Use nuclear_data
  Real(8), Parameter :: mex_p=7.28899, mex_n=8.07144
  Real(8)  :: y(ny),enb,enm,ytot,ztot,atot
  ytot=sum(y)
  ztot=sum(zz*y)
  atot=sum(aa*y)
  enb =sum(be*y)
  enm =mex_p*ztot+mex_n*(1.0-ztot)-enb
  
! Change units from MeV/nucleon to erg/g
  enb=epmev*avn*enb
  enm=epmev*avn*enm
  Return
End Subroutine benuc
  
Subroutine t9rhofind(kstep,tint,nn)
!===============================================================================
! This routine calculates t9 and rho as a function of time, either via 
! interpolation or from an analytic expression
!===============================================================================
  Use conditions
  Use thermo_data
  Use controls
  
  Integer :: n,nn,kstep
  Real*8  :: tint
  Real(8) :: dt,rdt,dt9,drho
  
! Calculate T9 and rho by interpolation 
  Do n=1,nh
    If(tint<=th(n)) Exit
  EndDo
  nn = n
  If(n>1.and.n<=nh) Then
    rdt=1.0/(th(n)-th(n-1))
    dt=tint-th(n-1)
    dt9=t9h(n)-t9h(n-1)
    drho=rhoh(n)-rhoh(n-1)
    t9t=dt*rdt*dt9+t9h(n-1)
    rhot=dt*rdt*drho+rhoh(n-1)
  ElseIf(n==1) Then
    t9t=t9h(1)
    rhot=rhoh(1)
  Else
    t9t=t9h(nh)
    rhot=rhoh(nh)
    Write(6,*) 'Time beyond thermodynamic range',tt,' >',th(nh)
  EndIf
  
! Calculate T9 and rho by function
! chi=1.0
! thd=chi*446.0/sqrt(rhostart)
! tint=(tstart-tt)/thd
! rhot=rhostart*exp(tint)
! t9t=t9start*exp(tint/3.)
  
! Output T9 and rho
! Write(lun_diag,"(a5,i5,3es12.4)") 'T9rho',kstep,tt,t9t,rhot
  
! Turn off strong interactions if the temperature is less than t9cut
  If (t9t<t9cut) Then
    If (iweak /= 0) Then
      iweak = -1
    Else
      iweak=iweak0
      Write(6,*)'Warning: Strong reactions ignored for T9< T9cut.'
    EndIf
  EndIf
  
  Return
End Subroutine t9rhofind
  
Subroutine timestep(kstep)
!===============================================================================
! This routine calculates the trial timestep.  
! For tdel >0, this calculation is based on the relative changes of the 
!              abundances during the previous step in the evolution and
!              the integrators estimate of the next timestep.  
!          =0, the timestep is calculated using the time derivatives of 
!              the abundances based on reaction rates.  
!          <0, the timestep is held constant, tdel=-tdel_old.
! There is also the provision for limiting the timestep in the event that the 
! thermodynamic conditions are changing too rapidly.
!===============================================================================
  Use nuclear_data
  Use controls
  Use conditions
  Use abundances
  Use thermo_data
  Integer, Save :: ints(1),intso(1) ! nucleus governing timestep
  Integer :: i,j,n,kstep,nnew,nold
  Real(8), Dimension(ny) :: ydotoy
  Real(8) :: changest,changeth,t9old,rhold,dt,dtherm,ttest
  Real(8) :: tdel_old,tdel_stop,tdel_deriv,tdel_fine
  
! Retain old values of timestep and thermo and calculate remaining time  
  changeth=.1
  changest=.01
  tdel_old=tdel
  t9old=t9t
  rhold=rhot
  tdel_stop=tstop-t
  tdel_fine=0.0
  
! If this is not the initial timestep, calculate timestep from changes in last timestep.
  If(tdel_old>0.0) Then
    Where(y>yacc)
      ydotoy=abs((y-yo)/y)
    ElseWhere
      ydotoy=0.0
    EndWhere
    ints=maxloc(ydotoy)
    If(ydotoy(ints(1)).ne.0) Then
      tdel_deriv=changemx*tdel_old/ydotoy(ints(1))
!     tdel_next=tdel_deriv
    Else
      tdel_deriv=tdel_next
    EndIf
    tdel=min(tdel_deriv,tdel_stop,tdel_next)
  
! If this is an initial timestep, yo does not exist, so calculate 
! timestep directly from derivatives.  
  
  ElseIf(tdel_old==0.0) Then
    tdel_stop=1.0e-4*tdel_stop
    intso(1)=0
    Call cross_sect
    Call yderiv
    Where(y>yacc)
      ydotoy=abs(ydot/y)
    ElseWhere
      ydotoy=0.0
    EndWhere
  
! If derivatives are non-zero, Use Y/(dY/dt).
    ints=maxloc(ydotoy)
    If(ydotoy(ints(1)).ne.0) Then
      tdel_deriv=changemx/ydotoy(ints(1))
! For unevolved initial abundances, as often found in test problems, 
! tdel_deriv may produce large jumps from zero abundance in the first 
! timestep. While generally inconsequential, tdel_fine limits these 
! to the accuracy abundance limit.
!    If(t9t.le.1.0) Then
!      changest=.00001
!    ElseIf(t9t.le.3.0) Then
!      changest=.0001
!    Else
!      changest=.001
!    EndIf
      tdel_fine=changest*tdel_deriv
!   tdel_fine=1.0*yacc/maxval(abs(ydot),y<yacc)
! If derivatives are zero, take a small step.
    Else
      tdel_deriv=tdel_stop
      tdel_fine=tdel_stop
    EndIf
    tdel=min(tdel_stop,tdel_deriv,tdel_fine)
  
! Keep timestep constant
  Else
    tdel=-tdel_old
  EndIf
  
! Diagnostic Output
  If(idiag>=1) Write(lun_diag,"(a4,i5,6es12.4,i5)") &
&   'tdel',kstep,tdel,tdel_deriv,tdel_old,tdel_stop,tdel_next,tdel_fine,ints(1)
! If(idiag>=2) Then
!   Write(lun_diag,"(a5,i4,2es12.4)") (nname(k),k,y(k),ydotoy(k),k=1,ny)
! EndIf
  
! Retain the index of the species setting the timestep 
  If(ints(1)/=intso(1)) Then
    If(idiag>=1) Write(lun_diag,*) 'ITC ',nname(ints(1)),t,tdel
    intso=ints
  EndIf
  
! Make sure to not skip any features in the temperature or density profiles by checking
! for profile monotonicity between t and t+del      
  Call t9rhofind(kstep,t,nold)
  ttest = t+tdel
  Call t9rhofind(kstep,ttest,nnew)
  If (nnew-1>nold) then
    do j=nold,nnew-1
      If (t9h(j)>t9old.and.t9h(j)>t9t) then
        tdel = th(j) - t
        exit
      elseIf (t9h(j)<t9old.and.t9h(j)<t9t) then
        tdel = th(j) - t
        exit
      elseIf (rhoh(j)>rhold.and.rhoh(j)>rhot) then
        tdel = th(j) - t
        exit
      elseIf (rhoh(j)<rhold.and.rhoh(j)<rhot) then
        tdel = th(j) - t
        exit    
      EndIf
    EndDo
  EndIf  
  
! Limit timestep if fractional density change is larger than changeth (10% by default)
! or fraction temperature change is larger than 0.1*changeth (1% by default)
  Do i=1,10
    ttest=t+tdel
! Reset the iweak flag to its original value
    iweak = iweak0
    Call t9rhofind(kstep,ttest,n)
    If(t9old>0) Then
      dtherm=10.0*abs(t9t-t9old)/t9old+abs(rhot-rhold)/rhold
    EndIf
    If(dtherm<changeth) Exit
    tdel=.5*tdel
    If(i==10) Write(lun_diag,*) 'Error in Thermo variations after ',i,&
&     'reductions',tdel,t9t,rhot
  EndDo
  tt=t+tdel
  
  If(idiag>=1) Write(lun_diag,"(a5,i5,2es12.4)") 'T9del',kstep,tdel,dtherm
  Return
End Subroutine timestep                                                                 
  
Subroutine partf(t9) 
!===============================================================================
! This routine calculates the nuclear partition functions as a function of T9
!===============================================================================
  Use part_funct_data
  Integer :: i,ii
  Real(8) :: t9,rdt9
  Do i=1,24                                                              
    If(t9<=t9i(i)) Exit
  EndDo                                                                     
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
  
! gg(0) is a placeholder for non-nuclei, gamma-rays, etc. 
  gg(0)=1.0
  
!If(idiag>=1) Then
!  Write(lun_diag,"(a5,i3,es14.7)") 'PartF',ii,t9
!  Write(lun_diag,"(5(i4,es12.4))") (i, gg(i), i=1,ny)
!EndIf
  Return                                                                    
End Subroutine partf                                                                      
  
Subroutine cross_sect
!===============================================================================
! This routine calculates the cross section for each reaction.
!===============================================================================
  Use controls
  Use nuclear_data
  Use abundances
  Use conditions
  USE screen_qse
  Use cross_sect_data
  Use part_funct_data 
  Use ffn_data
  Use reac_rate_data
  Real(8) :: t09(7)
  Real(8) :: ene,ye,yps,y2e,v,pa,ea,dlma,pta,pva,eta,eva,emas,beta
  Real(8) :: funct,rat,rhot2,zz1,zz2,zz3,zz12,aa12,h01,h02,fsr,fst
  Integer :: j,k,mu,ieos,iout
  emas=4.
  ieos=0
  eta=1.e+30                                                                
  
! Call EOS, necessary for screening
 Call en(ny,yt,ye,yps,y2e)
  If(iscrn>=1) Then
    V=1./rhot                                                       
    Call state(t9t,v,y,pa,ea,dlma,ye,pta,pva,eta,eva,ieos,ny,emas,yps,beta)
  EndIf                                                     
  
! Calculate necessary thermodynamic moments
  ene=ye*rhot                                                               
  t09(1)=1.                                                                 
  t09(2)=t9t**(-1) 
  t09(3)=t9t**(-1./3.) 
  t09(4)=t9t**(+1./3.) 
  t09(5)=t9t          
  t09(6)=t9t**(5./3.)
  t09(7)=log(t9t)  
!    Write(lun_diag,"(a3,4es12.4)") 'THR',t9t,rhot,ene,ye
  
! Calculate partition functions for each nucleus at t9t
  Call partf(t9t) 
  
! If there are any FFN reactions, calculate their rates
  If(nffn>0) Then
    Call ffn_rate(t9t,ene)   
  EndIf
  
! Calculate csects for reactions with one reactant
  funct=1.0
  rat=1.0
  Where(irev1==1)                 ! If it's a reverse reaction
    rpf1=gg(n1i(2,:))*gg(n1i(3,:))/gg(n1i(1,:))
  ElseWhere
    rpf1=1.0
  EndWhere
  If(iweak>0) Then           ! All Rates on
    Where(iwk1/=2.and.iwk1/=3)      ! If it's not an FFN reaction
      csect1=rpf1*exp(matmul(t09,rc1))
    ElseWhere
      iffn=int(rc1(1,:))
      csect1=rat*funct*rf(iffn)
    EndWhere
    Where (iwk1==1) csect1=ene*csect1 ! If it's a non-FFN EC reaction
 write(805,*) ene
  ElseIf(iweak<0) Then       ! Only weak rates Used
    Where(iwk1==0) csect1=0.0
    Where(iwk1==1) csect1=ene*rpf1*exp(matmul(t09,rc1))
    Where(iwk1==2.or.iwk1==3) 
      iffn=int(rc1(1,:))
      csect1=rat*funct*rf(iffn)
    EndWhere
    Where(iwk1>=4) csect1=rpf1*exp(matmul(t09,rc1))    
  Else                       ! weak interactions are off (iweak==0)
    Where(iwk1==0) 
      csect1=rpf1*exp(matmul(t09,rc1))  
    ElseWhere
      csect1=0.0
    EndWhere
  EndIf
  If(idiag>=5) Then
    Write(lun_diag,"(a,i5)") 'CSect1',nreac(1)
    Write(lun_diag,"(i5,5a5,3i3,es17.9)") &
&    (k,nname(n1i(1,k)),'-->',(nname(n1i(j,k)),j=2,3),'+++',&
&    iwk1(k),irev1(k),ires1(k),csect1(k), k=1,nreac(1))
  EndIf
 
! Calculate screening corrections for 2 reactant reactions
  If(iscrn>=1) Call screen2
  If(iscrn>=1) Call qse_screen(t9t,rhot,ye,iout)
  
! Calculate the csect for reactions with 2 reactants
  Where(irev2==1) ! If it's a reverse reaction mult. by ratio of part. func.
    rpf2=gg(n2i(3,:))*gg(n2i(4,:))/(gg(n2i(1,:))*gg(n2i(2,:)))
  ElseWhere
    rpf2=1.0
  EndWhere
  If(iweak>0) Then      ! All rates are on
    csect2=rhot*rpf2*exp(matmul(t09,rc2)+h2)
    Where(iwk2==1) csect2=ene*csect2! If it's a non FFN EC reaction
  ElseIf(iweak<0) Then  ! Only weak rates
    Where(iwk2==0) 
      csect2=0.0
    ElseWhere
      csect2=rhot*rpf2*exp(matmul(t09,rc2)+h2)
    EndWhere
    Where(iwk2==1) csect2=ene*csect2  
  Else                  ! Weak interactions are off  (iweak=0)
    Where(iwk2==0)
      csect2=rhot*rpf2*exp(matmul(t09,rc2)+h2)           
    ElseWhere
      csect2=0.0
    EndWhere
  EndIf
  If(idiag>=5) Then
    Write(lun_diag,"(a,i5)") 'CSect2',nreac(2)
    Write(lun_diag,"(i5,5a5,3i3,es17.9)") &
&      (k,(nname(n2i(j,k)),j=1,2),'-->',(nname(n2i(j,k)),j=3,4), &
&      iwk2(k),irev2(k),ires2(k),csect2(k),k=1,nreac(2))
  EndIf
  If(idiag>=5) Write(lun_diag,"(a,i5)") 'CSect3',nreac(3)
  rhot2=rhot**2
  Do mu=1,nreac(3)
!! Compute screening corrections
    If(iscrn>=1) Then
      zz1=zz(n3i(1,mu))                                 
      zz2=zz(n3i(2,mu))                                
      zz3=zz(n3i(3,mu))                               
      zz12=zz1+zz2                              
      aa12=aa(n3i(1,mu))+aa(n3i(2,mu))                   
      If(zz1*zz2.eq.0.0) Then
        h01=0.0
      Else
        Call screen(zz1,zz2,aa(n3i(1,mu)),aa(n3i(2,mu)),h01,fsr,fst) 
      EndIf
      If(zz12*zz3.eq.0.0) Then                          
        h02=0.0    
      Else
        Call screen(zz12,zz3,aa12,aa(n3i(3,mu)),h02,fsr,fst)  
      EndIf
    Else
     h01=0.0 ; h02=0.0
    EndIf
!    write(310,*) zz1,zz2,zz3,aa12,h01,h02  
! Compute the crossection
    If(iweak>0) Then
      csect3(mu)=t09(1)*rc3(1,mu)+t09(2)*rc3(2,mu)+t09(3)*rc3(3,mu)+ &
&       t09(4)*rc3(4,mu)+t09(5)*rc3(5,mu)+t09(6)*rc3(6,mu)+t09(7)*rc3(7,mu)
      csect3(mu)=rhot2*exp(csect3(mu)+h01+h02)
      If(iwk3(mu).eq.1) csect3(mu)=csect3(mu)*ene        
      If(csect3(mu).lt.1.e-20) csect3(mu)=0.0                   
    ElseIf(iweak<0) Then! Only Weak Reactions on
      If(iwk3(mu)==0) csect3(mu)=0. 
    Else                ! Weak Reactions off (iweak=0)
      If(iwk3(mu)==0) Then
        csect3(mu)=t09(1)*rc3(1,mu)+t09(2)*rc3(2,mu)+t09(3)*rc3(3,mu)+ &
&         t09(4)*rc3(4,mu)+t09(5)*rc3(5,mu)+t09(6)*rc3(6,mu)+t09(7)*rc3(7,mu)
        csect3(mu)=rhot2*exp(csect3(mu)+h01+h02)
        If(csect3(mu).lt.1.e-20) csect3(mu)=0.0 
      Else
        csect3(mu)=0.
      EndIf  
    EndIf
    If(idiag>=5) Write(lun_diag,"(i5,5a5,3i3,es17.9)") &
&     mu,(nname(n3i(j,mu)),j=1,3),'-->','+++',iwk3(mu),irev3(mu),ires3(mu),csect3(mu)
  EndDo    
  Return                                  
End Subroutine cross_sect                                   
  
  
Subroutine norm(yy)
!---------------------------------------------------------------------------- 
! This routine renormalizes the abundances to guarantee mass conservation.
!---------------------------------------------------------------------------- 
  Use nuclear_data
  Real(8) :: xtot,rxt
  Real(8), Dimension(ny) :: yy
  xtot=sum(yy*aa)   
  rxt=1.0/xtot      
  yy=yy*rxt  
  Return                                               
End Subroutine norm                                                
  
Subroutine ye_norm(yy,ye)             
!-----------------------------------------------------------------------------  
! This routine renormalizes the abundances to guarantee mass and charge
! conservation if Ye is specified.
!-----------------------------------------------------------------------------  
  Use nuclear_data
  Real(8) :: ye,zy,nny,zny,zzy,beta,alph
  Real(8) :: yy(ny)
  nny=sum(nn*yy)
  zy=sum(zz*yy)
  zny=sum(nn*zz*yy/aa)
  zzy=sum(zz*zz*yy/aa)
  beta=(ye*nny-zny)/(nny*zzy-zy*zny)
  alph=(1-beta*zy)/nny
  yy=yy*(alph*nn+beta*zz)/aa
  Return
End Subroutine ye_norm
  
  
  
