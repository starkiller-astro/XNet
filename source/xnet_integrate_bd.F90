!***************************************************************************************************
! xnet_integrate_bd.f90 10/18/17
! Bader-Deufelhard solver
!
! The routines in this file perform the Bader-Deufelhard time integration for thermonuclear reaction
! network. Based on routines from Press et al (1992,1996) "Numerical Recipes" and inspired by
! Timmes (1999; ApJS 124 241-263).
!***************************************************************************************************

Module xnet_integrate_bd
  Implicit None

Contains

  Subroutine solve_bd(kstep,its)
    !-----------------------------------------------------------------------------------------------
    ! This routine performs semi-implicit extrapolation of integrations returned by step_bd. Success
    ! of a timestep is judged against truncation errors from successive integrations of increasing
    ! order. Subsequent timesteps are also based on these truncation errors.
    ! Based on stifbs routine of Numerical Recipes.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny, nname
    Use reaction_data
    Use xnet_abundances, Only: y, ydot, yo, yt
    Use xnet_conditions
    Use xnet_controls, Only: idiag, iheat, kmon, ktot, lun_diag, lun_stdin, lun_stdout, tolm, yacc
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep

    ! Output variables
    Integer, Intent(out) :: its

    ! Local variables
    Real(dp), Parameter :: safe1 = 0.25, safe2 = 0.7, scalmx = 0.1
    Real(dp), Parameter :: redmax = 1.0e-5, redmin = 0.7, errmin = 1.0e-30
    Integer, Parameter :: imax = 8, kmaxx = imax-1, ktsmx = 100
    Integer, Parameter :: nseq(imax) = (/ 2, 6, 10, 14, 22, 34, 50, 70 /)
    Real(dp), Save :: alf(kmaxx,kmaxx), a(imax)
    Integer, Save :: kopt, kmax
    Real(dp) :: yerr(ny), yseq(ny), yscal(ny), err(kmaxx), t9seq, t9err
    Real(dp) :: eps, eps1, errmax, fact, red, scal, wrkmin, test
    Integer :: i, iq, k, km, kts, iminloc
    Logical :: reduct

    ! Set accuracy limits
    eps = tolm
    yscal = max(yacc,abs(y))

    tt = t + tdel
    If ( idiag >= 2 ) Write(lun_diag,"(a,i4,2es10.3)") 'Solve_BD',kstep,t,tt

    If ( kstep == 1 ) Then

      ! Initialize values
      kopt = kmaxx
      eps1 = safe1*eps

      ! Compute work coeffs a(k)
      a(1) = nseq(1) + 1
      Do k = 1, kmaxx
        a(k+1) = a(k) + nseq(k+1)
      EndDo

      ! Compute alf
      Do iq = 2, kmaxx
        Do k = 1, iq-1
          alf(k,iq) = eps1**( (a(k+1) - a(iq+1)) / ((a(iq+1) - a(1) + 1.0) * real((2*k+1),dp)) )
        EndDo
      EndDo

      a(1) = ny + a(1)
      Do k = 1, kmaxx
        a(k+1) = a(k) + nseq(k+1)
      EndDo

      Do kopt = 2, kmaxx-1
        If ( a(kopt+1) > a(kopt)*alf(kopt-1,kopt) ) Exit
      EndDo
      kmax = kopt
    EndIf

    reduct = .false.

  TS: Do kts = 1, ktsmx
      k = 0
      ktot(2) = ktot(2) + k
      Do k = 1, kmax
        tt = t + tdel
        If ( tt <= t ) Then
          Write(lun_stdout,*) 'Stepsize zero in solve_bd'
          Read(lun_stdin,*)
        EndIf
        If ( idiag >=2 ) Write(lun_diag,"(a7,2i4,2es10.3)") 'Step_BD',kts,k,t,tdel

        ! Reset temperature, density, electron fraction, and abundances
        If ( iheat > 0 ) t9t = t9
        yet = ye
        yt = y

        ! Attempt BD integration
        Call step_bd(kstep,nseq(k),yseq,t9seq)
        test = (tdel / nseq(k))**2

        ! Perform Richardson Extrapolation
        If ( iheat > 0 ) Then
          Call poly_extr_t9(k,test,yseq,yt,yerr,t9seq,t9t,t9err)
        Else
          Call poly_extr(k,test,yseq,yt,yerr)
        EndIf

        If ( idiag >= 2 ) Then
          Write(lun_diag,"(a5,i4,3es10.3)") 'Extrp',nseq(k),t,test,tdel
          Write(lun_diag,"(3es10.3)") (yt(i),yseq(i),yerr(i),i=1,ny)
          If ( iheat > 0 ) Write(lun_diag,"(3es10.3)") t9t,t9seq,t9err
        EndIf

        If ( k /= 1 ) Then
          errmax = maxval(abs(yerr / yscal))
          errmax = max(errmin,errmax) / eps
          If ( iheat > 0 ) errmax = max(t9err,errmax)
          km = k - 1
          err(km) = (errmax/safe1)**(1.0/real(2*km+1,dp))
        EndIf

        If ( k /= 1 .and. ( k >= kopt-1 .or. kstep==1 ) ) Then
          If ( errmax < 1.0 ) Exit TS
          If ( k == kmax .or. k == kopt+1 ) Then
            red = safe2 / err(km)
            Exit
          ElseIf ( k == kopt ) Then
            If ( alf(kopt-1,kopt) < err(km) ) Then
              red = 1.0 / err(km)
              Exit
            EndIf
          ElseIf ( kopt == kmax ) Then
            If ( alf(km,kmax-1) < err(km) ) Then
              red = alf(km,kmax-1) * safe2 / err(km)
              Exit
            EndIf
          ElseIf ( alf(km,kopt) < err(km) ) Then
            red = alf(km,kopt-1) / err(km)
            Exit
          EndIf
        EndIf
      EndDo
      ktot(2) = ktot(2) + k
      red = max(min(red,redmin),redmax)
      tdel = tdel*red
      reduct = .true.
    EndDo TS
    iminloc = minloc(a(2:km+1) * max(err(1:km),scalmx),dim=1)
    kopt = 1 + iminloc
    scal = max(err(kopt-1),scalmx)
    wrkmin = scal * a(kopt)
    tdel_next = tdel / scal
    If ( kopt >= k .and. kopt /= kmax .and. .not. reduct ) Then
      fact = max(scal/alf(kopt-1,kopt),scalmx)
      If ( a(kopt+1)*fact <= wrkmin ) Then
        tdel_next = tdel/fact
        kopt = kopt+1
      EndIf
    EndIf

    ! Update time, time_step, and abundances for successful timestep
    If ( kts <= ktsmx ) Then
      its = 0
      nto = nt   ; nt = ntt
      to = t     ; t = tt
      t9o = t9   ; t9 = t9t
      rhoo = rho ; rho = rhot
      yeo = ye   ; ye = yet
      yo = y     ; y = yt
    Else
      its = 1
    EndIf
    kmon(1) = kts; kmon(2) = k
    ktot(1) = ktot(1) + kts
    Write(lun_diag,"(a,4es10.3)") 'dt',tdel,tdel_next

    Return
  End Subroutine solve_bd

  Subroutine step_bd(kstep,numsteps,yout,t9out)
    !-----------------------------------------------------------------------------------------------
    ! This routine performs one step of the semi-implicit midpoint rule, starting from initial
    ! abundances, y, to final abundances, yout. Based on simpr routine of Numerical Recipes.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data
    Use xnet_abundances
    Use xnet_conditions
    Use xnet_controls
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep, numsteps

    ! Output variables
    Real(dp), Intent(out) :: yout(ny), t9out

    ! Local variables
    Real(dp) :: yrhs(ny), del(ny), dy(ny)
    Real(dp) :: t9rhs, dt9, delt9, relt9
    Real(dp) :: d, ttemp, dt
    Integer :: i, j, step_loop

    ! Set length of sub-timesteps
    ttemp = t
    dt = tdel / numsteps

    ! Calculate the thermodynamic factors necessary for reaction rates, including
    ! screening, and the reaction rates, if thermodynamic conditions are changing.
    If ( nh > 1 .or. iheat > 0 ) Call cross_sect

    ! Calculate the reaction rates and abundance time derivatives
    Call yderiv

    ! Build the Jacobian and LU decompose
    Call jacobian_build(1.0,-dt)
    Call jacobian_decomp(kstep)

    ! First step
    ! Calculate the right hand side
    yrhs = dt*ydot
    If ( iheat > 0 ) t9rhs = 0.0

    If ( idiag >= 2 ) Write(lun_diag,"(i10,(7es10.3))") 0,ttemp,yt
    If ( idiag >= 4 ) Then
      Write(lun_diag,"(a2,i5,es14.7)") 'First',kstep,dt
      Write(lun_diag,"(a5,4es23.15)") (nname(i),yrhs(i),ydot(i),yt(i),y(i),i=1,ny)
      If ( iheat > 0 ) Write(lun_diag,"(a9,3es23.15)") '   T9',t9rhs,t9dot,t9t,t9
    EndIf

    ! Do the back substitution
    Call jacobian_bksub(kstep,yrhs,dy,t9rhs,dt9)

    ! Save values for the next iteration
    del = dy
    yt = y + del
    If ( iheat > 0 ) Then
      delt9 = dt9
      t9t = t9t + delt9
      Call cross_sect
    EndIf

    ! Advance time
    ttemp = ttemp + dt
    If ( idiag >= 2 ) Write(lun_diag,"(i10,(7es10.3))") 1,ttemp,yt

    Call yderiv

    ! For intermediate steps
    Do step_loop = 2, numsteps

      ! Calculate right hand side
      yrhs = dt*ydot - del
      If ( iheat > 0 ) t9rhs = dt*t9dot - delt9

      ! Perform back substitution
      Call jacobian_bksub(kstep,yrhs,dy,t9rhs,dt9)

      ! Save values for the next iteration
      del = del + 2.0*dy
      yt = yt + del
      If ( iheat > 0 ) Then
        delt9 = delt9 + 2.0*dt9
        t9t = t9t + delt9
        Call cross_sect
      EndIf

      ! Advance time
      ttemp = ttemp + dt
      If ( idiag >= 2 ) Write(lun_diag,"(i10,(7es10.3))") step_loop,ttemp,yt
      Call yderiv
    EndDo

    ! Final step
    ! Calculate right hand side
    yrhs = dt*ydot - del
    If ( iheat > 0 ) t9rhs = dt*t9dot - delt9

    ! Perform back substition
    Call jacobian_bksub(kstep,yrhs,dy,t9rhs,dt9)

    ! Use yout for the final population
    yt = yt + dy
    yout = yt
    If ( iheat > 0 ) Then
      t9t = t9t + dt9
      t9out = t9t
    EndIf
    If ( idiag >= 2 ) Write(lun_diag,"(i10,(7es10.3))") step_loop,ttemp,yout

    Return
  End Subroutine step_bd

  Subroutine poly_extr(iest,test,yest,yz,dy)
    !-----------------------------------------------------------------------------------------------
    ! This routine performs a Richardson extrapolation, using a polynomial basis, of successive
    ! estimates of the abundances from finer timesliced integrations.
    ! Based on the psextr routine from Numerical Recipes.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use xnet_controls, Only: lun_stdout
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Integer, Intent(in) :: iest
    Real(dp), Intent(in) :: test, yest(ny)

    ! Output variables
    Real(dp), Intent(out) :: yz(ny), dy(ny)
    Integer, Parameter :: iest_max = 16
    Integer :: j
    Real(dp) :: delta, f1, f2
    Real(dp) :: d(ny), tmp(ny), q(ny)
    Real(dp), Save :: t(iest_max)
    Real(dp), Allocatable, Save :: qcol(:,:)
    !$omp threadprivate(qcol,t)

    If (iest > iest_max) Write(lun_stdout,*) 'Probable misuse, too much extrapolation'
    If ( .not. Allocated(qcol) ) Allocate(qcol(ny,iest_max))
    t(iest) = test
    dy = yest
    yz = yest
    If ( iest == 1 ) Then
      qcol(:,1) = yest
    Else
      d = yest
      Do j = 1, iest-1
        delta = 1.0 / (t(iest-j) - test)
        f1 = test * delta
        f2 = t(iest-j) * delta
        q = qcol(:,j)
        qcol(:,j) = dy
        tmp = d - q
        dy = f1 * tmp
        d = f2 * tmp
        yz = yz + dy
      EndDo
      qcol(:,iest) = dy
    EndIf

    Return
  End Subroutine poly_extr

  Subroutine poly_extr_t9(iest,test,yest,yz,dy,t9est,t9z,dt9)
    !-----------------------------------------------------------------------------------------------
    ! This routine performs a Richardson extrapolation, using a polynomial basis, of successive
    ! estimates of the abundances from finer timesliced integrations.
    ! Based on the psextr routine from Numerical Recipes.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use xnet_controls, Only: lun_stdout
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Integer, Intent(in) :: iest
    Real(dp), Intent(in) :: test, yest(ny), t9est

    ! Output variables
    Real(dp), Intent(out) :: yz(ny), dy(ny), t9z, dt9

    ! Local variables
    Integer, Parameter :: iest_max = 16
    Integer :: j
    Real(dp) :: delta, f1, f2
    Real(dp) :: d(ny), tmp(ny), q(ny)
    Real(dp), Save :: t(iest_max)
    Real(dp), Allocatable, Save :: qcol(:,:)
    !$omp threadprivate(qcol,t)

    If ( iest > iest_max ) Write(lun_stdout,*) 'Probable misuse, too much extrapolation'
    If ( .not. Allocated(qcol) ) Allocate(qcol(ny+1,iest_max))
    t(iest) = test
    dy = yest
    yz = yest
    dt9 = t9est
    t9z = t9est
    If ( iest == 1 ) Then
      qcol(:,1) = (/ yest, t9est /)
    Else
      d = (/ yest, t9est /)
      Do j = 1, iest-1
        delta = 1.0 / (t(iest-j) - test)
        f1 = test * delta
        f2 = t(iest-j) * delta
        q = qcol(:,j)
        qcol(:,j) = dy
        tmp = d - q
        dy = f1 * tmp
        dt9 = f1 * tmp(ny+1)
        d = f2 * tmp
        yz = yz + dy
        t9z = t9z + dt9
      EndDo
      qcol(:,iest) = (/ dy, dt9 /)
    EndIf

    Return
  End Subroutine poly_extr_t9

End Module xnet_integrate_bd
