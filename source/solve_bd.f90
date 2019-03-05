!*******************************************************************************
! Bader-Deufelhard solver, part of XNet 7 5/25/10
!
! The routines in this file perform the Bader-Deufelhard time integration for
! thermonuclear reaction network.  Based on routines from Press et al (1992,1996)
! "Numerical Recipes" and inspired by Timmes (1999; ApJS 124 241-263).
!
!*******************************************************************************

Subroutine solve_bd(kstep,its)
!------------------------------------------------------------------------------
! This routine performs semi-implicit extrapolation of integrations returned by
! step_bd.  Success of a timestep is judged against truncation errors from
! successive integrations of increasign order.  Subsequent timesteps are also
! based on these truncation errors.
! Based on stifbs routine of Numerical Recipes.
!------------------------------------------------------------------------------
  Use controls
  Use nuclear_data
  Use conditions
  Use abundances
  Use reac_rate_data

  Integer kstep,its,kts
  Real eps,yscal(ny),dy(ny)
  Integer, Parameter :: imax=8, kmaxx=imax-1, ktsmx=100
  Real, Parameter :: safe1=0.25,safe2=0.7,scalmx=0.1,redmax=1.0e-5,redmin=0.7,tiny=1.0e-30
  Integer :: i,iq,k,km,iminloc(1)
  Integer, Dimension(imax) :: nseq = (/2,6,10,14,22,34,50,70/)
  Integer, Save :: kopt,kmax
  Real, Dimension(kmaxx,kmaxx), Save :: alf
  Real, Dimension(kmaxx) :: err
  Real, Dimension(imax), Save :: a
  Real :: eps1,errmax,fact,red,scale,wrkmin,test
  Real, Dimension(size(y)) :: yerr,yseq
  Real :: t9err,t9seq
  Logical :: reduct

! Set accuracy limits
  eps=tolm
  Do i=1,ny
    yscal(i)=max(yacc,abs(y(i)))
  EndDo

  tt=t+tdel
  If(idiag>=2) Write(lun_diag,"(a,i4,2es10.3)")   'Solve_BD',kstep,t,tt

  If (kstep==1) Then
! Initialize values
    kopt=kmaxx
    eps1=safe1*eps

! Compute work coeffs a(k)
    a(1)=nseq(1)+1
    Do k=1,kmaxx
      a(k+1)=a(k)+nseq(k+1)
    EndDo

! Compute alf
    Do iq=2,kmaxx
      Do k=1,iq-1
        alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.)*(2*k+1)))
      EndDo
    EndDo

    a(1)=ny+a(1)
    Do k=1,kmaxx
      a(k+1)=a(k)+nseq(k+1)
    EndDo

    Do kopt=2,kmaxx-1
      If(a(kopt+1)>a(kopt)*alf(kopt-1,kopt)) exit
    EndDo
    kmax=kopt
  EndIf

  reduct=.false.

TS: Do kts=1,ktsmx
    k=0
    ktot(2)=ktot(2)+k
    Do k=1,kmax
      tt=t+tdel
      If(tt.eq.t) pause 'stepsize zero in solve_bd'
      If(idiag>=2) Write(lun_diag,"(a7,2i4,2es10.3)") 'Step_BD',kts,k,t,tdel
      yt=y
      If(iheat>0) t9t=t9
      Call step_bd(kstep,nseq(k),yseq,t9seq)
      test=(tdel/nseq(k))**2

! Perform Richardson Extrapolation
      If(iheat>0) Then
        Call poly_extr_t9(k,test,yseq,yt,yerr,t9seq,t9t,t9err)
      Else
        Call poly_extr(k,test,yseq,yt,yerr)
      EndIf
      If(idiag>=2) Then
        Write(lun_diag,"(a5,i4,3es10.3)") 'Extrp',nseq(k),t,test,tdel
        Write(lun_diag,"(3es10.3)") (yt(i),yseq(i),yerr(i),i=1,ny)
        If(iheat>0) Write(lun_diag,"(3es10.3)") t9t,t9seq,t9err
      EndIf
      If (k /= 1) Then
        errmax=maxval(abs(yerr(:)/yscal(:)))
        errmax=max(tiny,errmax)/eps
        If(iheat>0) errmax=max(t9err,errmax)
        km=k-1
        err(km)=(errmax/safe1)**(1.0/(2*km+1))
      EndIf
      If (k /= 1 .and. (k >= kopt-1 .or. kstep==1)) Then
        If (errmax < 1.0) exit TS
        If (k == kmax .or. k == kopt+1) Then
          red=safe2/err(km)
          exit
        Else If (k == kopt) Then
          If (alf(kopt-1,kopt) < err(km)) Then
            red=1.0/err(km)
            exit
          EndIf
        Else If (kopt == kmax) Then
          If (alf(km,kmax-1) < err(km)) Then
            red=alf(km,kmax-1)*safe2/err(km)
            exit
          EndIf
        Else If (alf(km,kopt) < err(km)) Then
          red=alf(km,kopt-1)/err(km)
          exit
        EndIf
      EndIf
    EndDo
    ktot(2)=ktot(2)+k
    red=max(min(red,redmin),redmax)
    tdel=tdel*red
    reduct=.true.
  EndDo TS
  iminloc=minloc(a(2:km+1)*max(err(1:km),SCALMX))
  kopt=1+iminloc(1)
  scale=max(err(kopt-1),SCALMX)
  wrkmin=scale*a(kopt)
  tdel_next=tdel/scale
  If (kopt >= k .and. kopt /= kmax .and. .not. reduct) Then
    fact=max(scale/alf(kopt-1,kopt),SCALMX)
    If (a(kopt+1)*fact <= wrkmin) Then
      tdel_next=tdel/fact
      kopt=kopt+1
    EndIf
  EndIf

! Update time, time_step, and abundances for successful timestep
  If(kts<=ktsmx) Then
    to=t
    t=tt
    yo=y
    y=yt
    dy=y-yo
    If(iheat>0) Then
      t9o=t9
      t9=t9t
    EndIf
    its=0
  Else
    its=1
  Endif
  kmon(1)=kts; kmon(2) = k
  ktot(1)=ktot(1)+kts
  Write(lun_diag,"(a,4es10.3)") 'dt',tdel,tdel_next
End Subroutine solve_bd

Subroutine step_bd(kstep,numsteps,yout,t9out)
!-----------------------------------------------------------------------------
! This routine performs one step of the semi-implicit midpoint rule, starting
! from initial abundances, y, to final abundances, yout.
! Based on simpr routine of Numerical Recipes.
!-----------------------------------------------------------------------------
  Use controls
  Use nuclear_data
  Use conditions
  Use abundances
  Use thermo_data

  Integer, Intent(in)                 :: kstep,numsteps
  Real(8), Intent(out), Dimension(ny) :: yout
  Real(8), Intent(out)                :: t9out
  Integer                             :: i,j,step_loop
  Real(8), Dimension(ny)              :: yrhs,del,dy
  Real(8)                             :: d,ttemp,dt
  Real(8)                             :: t9rhs,dt9,delt9,relt9

! Set length of sub-timesteps
  ttemp=t
  dt = tdel/numsteps

! Calculate the thermodynamic factors necessary for reaction rates, including
! screening, and the reaction rates, if thermodynamic conditions are changing.
  If(nh>1.or.iheat>0) Call cross_sect

! Calculate the reaction rates and abundance time derivatives
  Call yderiv

! Build the Jacobian and LU decompose
  Call jacobian_build(1.0,-dt)
  Call jacobian_decomp(kstep)

! First step
! Calculate the right hand side
  yrhs = dt*ydot
  If(iheat>0) t9rhs=0.0

  If(idiag>=2) Write(lun_diag,"(i10,(7es10.3))") 0,ttemp,yt

  If(idiag>=4) Then
    Write(lun_diag,"(a2,i5,es14.7)") 'First',kstep,dt
    Do i=1,ny
      Write(lun_diag,"(a5,4es17.9)") nname(i),yrhs(i),ydot(i),yt(i),y(i)
    EndDo
    If(iheat>0) Write(lun_diag,"(a9,3es17.9)") '   T9',t9rhs,t9dot,t9t,t9
  EndIf

! Do the back substitution
  Call jacobian_bksub(yrhs,dy,t9rhs,dt9)

! Save values for the next iteration
  del = dy
  yt = y+del
  If(iheat>0) Then
    delt9 = dt9
    t9t = t9t+delt9
    Call cross_sect
  EndIf

! Advance time
  ttemp=ttemp+dt
  If(idiag>=2) Write(lun_diag,"(i10,(7es10.3))") 1,ttemp,yt

  Call yderiv

! For intermediate steps
  Do step_loop = 2,numsteps

! Calculate right hand side
    yrhs = dt*ydot-del
    If(iheat>0) t9rhs=dt*t9dot-delt9

! Perform back substitution
    Call jacobian_bksub(yrhs,dy,t9rhs,dt9)

! Save values for the next iteration
    del = del+2.0*dy
    yt = yt+del
    If(iheat>0) Then
      delt9 = delt9+2.0*dt9
      t9t = t9t+delt9
      Call cross_sect
    EndIf

! Advance time
    ttemp=ttemp+dt
    If(idiag>=2) Write(lun_diag,"(i10,(7es10.3))") step_loop,ttemp,yt
    Call yderiv
  EndDo

! Final step
! Calculate right hand side
  yrhs = dt*ydot-del
  If(iheat>0) t9rhs=dt*t9dot-delt9

! Perform back substition
  Call jacobian_bksub(yrhs,dy,t9rhs,dt9)

! Use yout for the final population
  yt = yt + dy
  yout = yt
  If(iheat>0) Then
    t9t = t9t+dt9
    t9out = t9t
  EndIf
  If(idiag>=2) Write(lun_diag,"(i10,(7es10.3))") step_loop,ttemp,yout

End Subroutine step_bd

Subroutine poly_extr(iest,test,yest,yz,dy)
!-----------------------------------------------------------------------------
! This routine performs a Richardson extrapolation, using a polynomial basis,
! of successive estimates of the abundances from finer timesliced integrations.
! Based on the psextr routine from Numerical Recipes.
!-----------------------------------------------------------------------------
  Use nuc_number
  Integer, Intent(IN) :: iest
  Real, Intent(IN) :: test
  Real, Dimension(ny), Intent(IN) :: yest
  Real, Dimension(ny), Intent(OUT) :: yz,dy
  Integer, Parameter :: iest_max=16
  Integer :: j
  Real :: delta,f1,f2
  Real, Dimension(ny) :: d,tmp,q
  Real, Dimension(iest_max), Save :: t
  Real, Dimension(:,:), Allocatable, Save :: qcol
!$OMP THREADPRIVATE(qcol,t)
  If (iest > iest_max) Write(*,*) 'probable misuse, too much extrapolation'
  If (.not.Allocated(qcol)) Allocate(qcol(ny,iest_max))
  t(iest)=test
  dy=yest
  yz=yest
  If (iest == 1) Then
    qcol(:,1)=yest
  Else
    d=yest
    Do j=1,iest-1
      delta=1.0/(t(iest-j)-test)
      f1=test*delta
      f2=t(iest-j)*delta
      q=qcol(:,j)
      qcol(:,j)=dy
      tmp=d-q
      dy=f1*tmp
      d=f2*tmp
      yz=yz+dy
    EndDo
    qcol(:,iest)=dy
  EndIf
End Subroutine poly_extr

Subroutine poly_extr_t9(iest,test,yest,yz,dy,t9est,t9z,dt9)
!-----------------------------------------------------------------------------
! This routine performs a Richardson extrapolation, using a polynomial basis,
! of successive estimates of the abundances from finer timesliced integrations.
! Based on the psextr routine from Numerical Recipes.
!-----------------------------------------------------------------------------
  Use nuc_number
  Integer, Intent(IN) :: iest
  Real, Intent(IN) :: test
  Real, Dimension(ny), Intent(IN) :: yest
  Real, Dimension(ny), Intent(OUT) :: yz,dy
  Real, Intent(IN) :: t9est
  Real, Intent(OUT) :: t9z,dt9
  Integer, Parameter :: iest_max=16
  Integer :: j
  Real :: delta,f1,f2
  Real, Dimension(ny+1) :: d,tmp,q
  Real, Dimension(iest_max), Save :: t
  Real, Dimension(:,:), Allocatable, Save :: qcol
!$OMP THREADPRIVATE(qcol,t)
  If (iest > iest_max) Write(*,*) 'probable misuse, too much extrapolation'
  If (.not.Allocated(qcol)) Allocate(qcol(ny+1,iest_max))
  t(iest)=test
  dy=yest
  yz=yest
  dt9=t9est
  t9z=t9est
  If (iest == 1) Then
    qcol(:,1)=(/yest,t9est/)
  Else
    d=(/yest,t9est/)
    Do j=1,iest-1
      delta=1.0/(t(iest-j)-test)
      f1=test*delta
      f2=t(iest-j)*delta
      q=qcol(:,j)
      qcol(:,j)=(/dy,dt9/)
      tmp=d-q
      dy=f1*tmp(1:ny)
      dt9=f1*tmp(ny+1)
      d=f2*tmp
      yz=yz+dy
      t9z=t9z+dt9
    EndDo
    qcol(:,iest)=(/dy,dt9/)
  EndIf
End Subroutine poly_extr_t9
