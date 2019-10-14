!***************************************************************************************************
! jacobian_gpu.f90 10/18/17
! The routines in this file assume a dense Jacobian and use a dense linear algebra package.
!
! The bulk of the computational cost of the network (60-95%) is the solving of the matrix equation.
! Careful selection of the matrix solver is therefore very important to fast computation. For
! networks from a few dozen up to a couple hundred species, hand tuned dense solvers such as those
! supplied by the hardware manufacturer (often LAPACK) or third-parties like NAG, IMSL, etc. are
! fastest. However for larger matrices, sparse solvers are faster.
!*******************************************************************************

Module xnet_jacobian
  !-----------------------------------------------------------------------------------------------
  ! The Jacobian matrix for the solver.
  !-----------------------------------------------------------------------------------------------
  Use, Intrinsic :: iso_c_binding
  Use xnet_types, Only: dp
  Implicit None
  Real(dp), Allocatable, Target :: dydotdy(:,:,:) ! dYdot/dY part of jac
  Real(dp), Allocatable, Target :: jac(:,:,:)     ! the Jacobian matrix
  Real(dp), Allocatable, Target :: rhs(:,:)       ! the Jacobian matrix
  Real(dp), Allocatable, Target :: work(:,:)      ! the Jacobian matrix
  Integer, Allocatable, Target :: indx(:,:)       ! Pivots in the LU decomposition
  Integer, Allocatable, Target :: indxinfo(:,:)
  Integer, Allocatable, Target :: info(:)

  Type(C_PTR), Allocatable, Target :: hjac(:), hrhs(:)
  Type(C_PTR), Allocatable, Target :: djac(:), drhs(:), dwork(:), dindx(:), dindxinfo(:)

  ! Array size parameters
  Real(C_DOUBLE), Parameter :: ddummy = 0.0
  Integer(C_INT), Parameter :: idummy = 0
  Integer(C_INTPTR_T), Parameter :: cptr_dummy = 0
  Integer(C_SIZE_T), Parameter :: sizeof_double = c_sizeof(ddummy)
  Integer(C_SIZE_T), Parameter :: sizeof_int = c_sizeof(idummy)
  Integer(C_SIZE_T), Parameter :: sizeof_cptr = c_sizeof(cptr_dummy)
  Integer(C_SIZE_T) :: sizeof_jac, sizeof_rhs, sizeof_indx, sizeof_info, sizeof_batch

  ! Parameters for GPU array dimensions
  Integer :: msize ! Size of linear system to be solved

  Real(dp), Allocatable, Target :: diag0(:)
  Real(dp), Allocatable, Target :: mult1(:)

  Logical, Parameter :: pivot = .false.

  Type(C_PTR), Allocatable, Target :: djacp(:), drhsp(:), dindxp(:)

  !$acc declare &
  !$acc create(dydotdy,jac,rhs,work,indx,indxinfo,info, &
  !$acc        djac,drhs,dwork,dindx,dindxinfo,diag0,mult1, &
  !$acc        djacp,drhsp,dindxp)

Contains

  Subroutine read_jacobian_data(data_dir)
    !-------------------------------------------------------------------------------------------------
    ! Initializes the Jacobian data.
    !-------------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use xnet_controls, Only: iheat, nzevolve, nzbatchmx, tid
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: data_dir

    ! Local variables
    Integer :: istat, izb

    ! Calculate array sizes
    If ( iheat > 0 ) Then
      msize = ny + 1
    Else
      msize = ny
    EndIf

    Allocate (diag0(nzevolve))
    Allocate (mult1(nzevolve))
    diag0 = 0.0
    mult1 = 1.0

    Allocate (dydotdy(msize,msize,nzevolve))
    Allocate (jac(msize,msize,nzevolve))
    Allocate (rhs(msize,nzevolve))
    Allocate (work(msize,nzevolve))
    Allocate (indx(msize,nzevolve))
    Allocate (indxinfo(msize,nzevolve))
    Allocate (info(nzevolve))
    dydotdy = 0.0
    jac = 0.0
    rhs = 0.0
    work = 0.0
    indx = 0
    indxinfo = 0
    info = 0

    Allocate (hjac(nzevolve))
    Allocate (hrhs(nzevolve))
    Do izb = 1, nzevolve
      hjac(izb) = c_loc( jac(1,1,izb) )
      hrhs(izb) = c_loc( rhs(1,izb) )
    EndDo

    Allocate (djac(nzevolve))
    Allocate (djacp(nzevolve))
    Allocate (drhs(nzevolve))
    Allocate (drhsp(nzevolve))
    Allocate (dwork(nzevolve))
    Allocate (dindx(nzevolve))
    Allocate (dindxp(nzevolve))
    Allocate (dindxinfo(nzevolve))
    !$acc host_data use_device(jac,rhs,work,indx,indxinfo)
    Do izb = 1, nzevolve
      djac(izb) = c_loc( jac(1,1,izb) )
      djacp(izb) = c_loc( jac(1,1,izb) )
      drhs(izb) = c_loc( rhs(1,izb) )
      drhsp(izb) = c_loc( rhs(1,izb) )
      dwork(izb) = c_loc( work(1,izb) )
      dindx(izb) = c_loc( indx(1,izb) )
      dindxp(izb) = c_loc( indx(1,izb) )
      dindxinfo(izb) = c_loc( indxinfo(1,izb) )
    EndDo
    !$acc end host_data

    !$acc update async(tid) &
    !$acc device(dydotdy,jac,rhs,work,indx,indxinfo,info,djac, &
    !$acc        djacp,drhs,drhsp,dwork,dindx,dindxp,dindxinfo,diag0,mult1)

    Return
  End Subroutine read_jacobian_data

  Subroutine jacobian_finalize
    !-----------------------------------------------------------------------------------------------
    ! Free the page-locked and device memory used in the dense solver.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Local variables
    Integer :: istat

    Deallocate (diag0,mult1)
    Deallocate (dydotdy,jac,rhs,indx,info)
    Deallocate (hjac,hrhs)
    Deallocate (djac,drhs,dwork,dindx,dindxinfo)
    Deallocate (djacp,drhsp,dindxp)

  End Subroutine jacobian_finalize

  Subroutine jacobian_scale(diag_in,mult_in,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This augments a previously calculation Jacobian matrix by multiplying all elements by mult and
    ! adding diag to the diagonal elements.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny, nname
    Use xnet_controls, Only: idiag, iheat, lun_diag, szbatch, zb_lo, zb_hi, lzactive, tid
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Optional, Target, Intent(in) :: diag_in(zb_lo:zb_hi), mult_in(zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer :: i, j, i0, j1, izb, izone
    Real(dp), Pointer :: diag(:), mult(:)
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    If ( present(diag_in) ) Then
      diag(zb_lo:) => diag_in
    Else
      diag(zb_lo:) => diag0(zb_lo:zb_hi)
    EndIf
    If ( present(mult_in) ) Then
      mult(zb_lo:) => mult_in
    Else
      mult(zb_lo:) => mult1(zb_lo:zb_hi)
    EndIf

    !$acc enter data async(tid) &
    !$acc copyin(mask,diag,mult)

    !$acc parallel loop gang collapse(2) async(tid) &
    !$acc present(mask,mult,diag,jac,dydotdy)
    Do izb = zb_lo, zb_hi
      Do j1 = 1, msize
        If ( mask(izb) ) Then
          !$acc loop vector
          Do i0 = 1, msize
            jac(i0,j1,izb) = mult(izb) * dydotdy(j1,i0,izb)
          EndDo
          jac(j1,j1,izb) = jac(j1,j1,izb) + diag(izb)
        EndIf
      EndDo
    EndDo

    If ( idiag >= 5 ) Then
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          izone = izb + szbatch - zb_lo
          !$acc update wait(tid) &
          !$acc host(diag(izb),mult(izb),jac(:,:,izb))
          Write(lun_diag,"(a9,i5,2es24.16)") 'JAC_SCALE',izone,diag(izb),mult(izb)
          Do i = 1, ny
            Write(lun_diag,"(3a)") 'J(',nname(i),',Y)'
            Write(lun_diag,"(7es24.16)") (jac(i,j,izb),j=1,ny)
          EndDo
          If ( iheat > 0 ) Then
            Write(lun_diag,"(3a)") 'J(Y,T9)'
            Write(lun_diag,"(7es24.16)") (jac(i,ny+1,izb),i=1,ny)
            Write(lun_diag,"(a)") 'J(T9,Y)'
            Write(lun_diag,"(7es24.16)") (jac(ny+1,j,izb),j=1,ny)
            Write(lun_diag,"(a)") 'J(T9,T9)'
            Write(lun_diag,"(es24.16)") jac(ny+1,ny+1,izb)
          EndIf
        EndIf
      EndDo
      Flush(lun_diag)
    EndIf

    !$acc exit data async(tid) &
    !$acc delete(mask,diag,mult)
    
    Return
  End Subroutine jacobian_scale

  Subroutine jacobian_build(diag_in,mult_in,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the reaction Jacobian matrix, dYdot/dY, and augments by multiplying
    ! all elements by mult and adding diag to the diagonal elements.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny, mex, nname
    Use reaction_data, Only: a1, a2, a3, a4, b1, b2, b3, b4, la, le, mu1, mu2, mu3, mu4, n11, n21, &
      & n22, n31, n32, n33, n41, n42, n43, n44, dcsect1dt9, dcsect2dt9, dcsect3dt9, dcsect4dt9, nan, &
      & n10, n20, n30, n40
    Use xnet_abundances, Only: yt
    Use xnet_conditions, Only: cv
    Use xnet_controls, Only: iheat, idiag, ktot, lun_diag, nzbatchmx, szbatch, zb_lo, zb_hi, &
      & lzactive, tid
    Use xnet_timers, Only: xnet_wtime, start_timer, stop_timer, timer_jacob
    Use xnet_types, Only: dp
    Implicit None

    ! Optional variables
    Real(dp), Optional, Target, Intent(in) :: diag_in(zb_lo:zb_hi), mult_in(zb_lo:zb_hi)
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer :: i, j, i0, j1, izb, izone
    Integer :: la1, le1, la2, le2, la3, le3, la4, le4
    Integer :: i11, i21, i22, i31, i32, i33, i41, i42, i43, i44
    Real(dp) :: s1, s2, s3, s4, sdot, r1, r2, r3, r4
    Real(dp) :: y11, y21, y22, y31, y32, y33, y41, y42, y43, y44
    Real(dp) :: dt9dotdy(msize), dr1dt9, dr2dt9, dr3dt9, dr4dt9
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    start_timer = xnet_wtime()
    timer_jacob = timer_jacob - start_timer

    !$acc enter data async(tid) &
    !$acc copyin(mask)

    ! Build the Jacobian

    !$acc parallel loop gang collapse(2) async(tid) &
    !$acc present(mask,dydotdy,yt,b1,b2,b3,b4,la,le,cv,mex,ktot, &
    !$acc         n10,n11,n20,n21,n22,n30,n31,n32,n33,n40,n41,n42,n43,n44, &
    !$acc         dcsect1dt9,dcsect2dt9,dcsect3dt9,dcsect4dt9, &
    !$acc         mu1,mu2,mu3,mu4,a1,a2,a3,a4)
    Do izb = zb_lo, zb_hi
      Do i0 = 1, ny
        If ( mask(izb) ) Then
          !$acc loop vector
          Do j1 = 1, msize
            dydotdy(j1,i0,izb) = 0.0
          EndDo
          !$acc loop vector
          Do j1 = la(1,i0), le(1,i0)
            !$acc atomic
            dydotdy(n11(j1),i0,izb) = dydotdy(n11(j1),i0,izb) + b1(j1,izb)
          EndDo
          !$acc loop vector
          Do j1 = la(2,i0), le(2,i0)
            !$acc atomic
            dydotdy(n21(j1),i0,izb) = dydotdy(n21(j1),i0,izb) + b2(j1,izb) * yt(n22(j1),izb)
            !$acc atomic
            dydotdy(n22(j1),i0,izb) = dydotdy(n22(j1),i0,izb) + b2(j1,izb) * yt(n21(j1),izb)
          EndDo
          !$acc loop vector
          Do j1 = la(3,i0), le(3,i0)
            !$acc atomic
            dydotdy(n31(j1),i0,izb) = dydotdy(n31(j1),i0,izb) + b3(j1,izb) * yt(n32(j1),izb) * yt(n33(j1),izb)
            !$acc atomic
            dydotdy(n32(j1),i0,izb) = dydotdy(n32(j1),i0,izb) + b3(j1,izb) * yt(n33(j1),izb) * yt(n31(j1),izb)
            !$acc atomic
            dydotdy(n33(j1),i0,izb) = dydotdy(n33(j1),i0,izb) + b3(j1,izb) * yt(n31(j1),izb) * yt(n32(j1),izb)
          EndDo
          !$acc loop vector
          Do j1 = la(4,i0), le(4,i0)
            !$acc atomic
            dydotdy(n41(j1),i0,izb) = dydotdy(n41(j1),i0,izb) + b4(j1,izb) * yt(n42(j1),izb) * yt(n43(j1),izb) * yt(n44(j1),izb)
            !$acc atomic
            dydotdy(n42(j1),i0,izb) = dydotdy(n42(j1),i0,izb) + b4(j1,izb) * yt(n43(j1),izb) * yt(n44(j1),izb) * yt(n41(j1),izb)
            !$acc atomic
            dydotdy(n43(j1),i0,izb) = dydotdy(n43(j1),i0,izb) + b4(j1,izb) * yt(n44(j1),izb) * yt(n41(j1),izb) * yt(n42(j1),izb)
            !$acc atomic
            dydotdy(n44(j1),i0,izb) = dydotdy(n44(j1),i0,izb) + b4(j1,izb) * yt(n41(j1),izb) * yt(n42(j1),izb) * yt(n43(j1),izb)
          EndDo

          If ( iheat > 0 ) Then
            s1 = 0.0
            !$acc loop vector &
            !$acc reduction(+:s1)
            Do j1 = la(1,i0), le(1,i0)
              s1 = s1 + a1(j1) * dcsect1dt9(mu1(j1),izb) * yt(n11(j1),izb)
            EndDo
            s2 = 0.0
            !$acc loop vector &
            !$acc reduction(+:s2)
            Do j1 = la(2,i0), le(2,i0)
              s2 = s2 + a2(j1) * dcsect2dt9(mu2(j1),izb) * yt(n21(j1),izb) * yt(n22(j1),izb)
            EndDo
            s3 = 0.0
            !$acc loop vector &
            !$acc reduction(+:s3)
            Do j1 = la(3,i0), le(3,i0)
              s3 = s3 + a3(j1) * dcsect3dt9(mu3(j1),izb) * yt(n31(j1),izb) * yt(n32(j1),izb) * yt(n33(j1),izb)
            EndDo
            s4 = 0.0
            !$acc loop vector &
            !$acc reduction(+:s4)
            Do j1 = la(4,i0), le(4,i0)
              s4 = s4 + a4(j1) * dcsect4dt9(mu4(j1),izb) * yt(n41(j1),izb) * yt(n42(j1),izb) * yt(n43(j1),izb) * yt(n44(j1),izb)
            EndDo
            dydotdy(ny+1,i0,izb) = s1 + s2 + s3 + s4
          EndIf
        EndIf
      EndDo
    EndDo

    If ( iheat > 0 ) Then

      !$acc parallel loop gang collapse(2) async(tid) &
      !$acc present(mask,dydotdy,cv,mex) &
      !$acc private(sdot)
      Do izb = zb_lo, zb_hi
        Do j1 = 1, msize
          If ( mask(izb) ) Then
            sdot = 0.0
            !$acc loop vector &
            !$acc reduction(-:sdot)
            Do i0 = 1, ny
              sdot = sdot - mex(i0)*dydotdy(j1,i0,izb) / cv(izb)
            EndDo
            dydotdy(j1,ny+1,izb) = sdot
          EndIf
        EndDo
      EndDo
    EndIf

    !$acc parallel loop gang async(tid) &
    !$acc present(mask,ktot)
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        ktot(3,izb) = ktot(3,izb) + 1
      EndIf
    EndDo

    ! Apply the externally provided factors
    Call jacobian_scale(diag_in,mult_in,mask_in = mask)

    If ( idiag >= 5 ) Then
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          izone = izb + szbatch - zb_lo
          !$acc update wait(tid) &
          !$acc host(dydotdy(:,:,izb))
          Write(lun_diag,"(a9,i5)") 'JAC_BUILD',izone
          Do i = 1, ny
            Write(lun_diag,"(3a)") 'dYDOT(',nname(i),')/dY'
            Write(lun_diag,"(7es24.16)") (dydotdy(j,i,izb),j=1,ny)
          EndDo
          If ( iheat > 0 ) Then
            Write(lun_diag,"(3a)") 'dYDOT/dT9'
            Write(lun_diag,"(7es24.16)") (dydotdy(ny+1,i,izb),i=1,ny)
            Write(lun_diag,"(a)") 'dT9DOT/dY'
            Write(lun_diag,"(7es24.16)") (dydotdy(j,ny+1,izb),j=1,ny)
            Write(lun_diag,"(a)") 'dT9DOT/dT9'
            Write(lun_diag,"(es24.16)") dydotdy(ny+1,ny+1,izb)
          EndIf
        EndIf
      EndDo
      Flush(lun_diag)
    EndIf

    !$acc exit data async(tid) &
    !$acc delete(mask)

    stop_timer = xnet_wtime()
    timer_jacob = timer_jacob + stop_timer

    Return
  End Subroutine jacobian_build

  Subroutine jacobian_solve(kstep,yrhs,dy,t9rhs,dt9,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine solves the system of equations composed of the Jacobian and RHS vector.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use xnet_controls, Only: idiag, iheat, lun_diag, nzbatch, szbatch, zb_lo, zb_hi, lzactive, &
      & tid
    Use xnet_linalg, Only: LinearSolveBatched_GPU, LinearSolve_CPU
    Use xnet_timers, Only: xnet_wtime, start_timer, stop_timer, timer_solve
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep
    Real(dp), Intent(in) :: yrhs(ny,zb_lo:zb_hi)
    Real(dp), Intent(in) :: t9rhs(zb_lo:zb_hi)

    ! Output variables
    Real(dp), Intent(out) :: dy(ny,zb_lo:zb_hi)
    Real(dp), Intent(out) :: dt9(zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer :: i, izb, izb_p, nzmask, izone, istat
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return
    nzmask = count(mask)

    start_timer = xnet_wtime()
    timer_solve = timer_solve - start_timer

    !$acc enter data async(tid) &
    !$acc copyin(mask)

    !$acc parallel loop gang async(tid) &
    !$acc present(mask,yrhs,t9rhs,rhs)
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        !$acc loop vector
        Do i = 1, ny
          rhs(i,izb) = yrhs(i,izb)
        EndDo
        If ( iheat > 0 ) rhs(ny+1,izb) = t9rhs(izb)!*1.0e9
      Else
        !$acc loop vector
        Do i = 1, msize
          rhs(i,izb) = 0.0
        EndDo
      EndIf
    EndDo

    ! Solve the linear system
#if defined(XNET_GPU)
    !$acc serial async(tid) &
    !$acc present(mask,djacp,djac,drhsp,drhs,dindxp,dindx)
    i = 0
    !$acc loop &
    !$acc private(izb_p) &
    !$acc reduction(+:i)
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        i = i + 1
        izb_p = zb_lo + i - 1
        djacp(izb_p) = djac(izb)
        drhsp(izb_p) = drhs(izb)
        dindxp(izb_p) = dindx(izb)
      EndIf
    EndDo
    !$acc end serial
    call LinearSolveBatched_GPU &
      & ( 'N', msize, 1, djacp(zb_lo), msize, dindxp(zb_lo), dindxinfo(zb_lo), &
      &   drhsp(zb_lo), msize, dwork(zb_lo), info(zb_lo), nzmask )
#else
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        call LinearSolve_CPU &
          & ( 'N', msize, 1, jac(1,1,izb), msize, indx(1,izb), rhs(1,izb), msize, info(izb) )
      EndIf
    EndDo
#endif

    !$acc parallel loop gang async(tid) &
    !$acc present(mask,dy,dt9,rhs)
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        !$acc loop vector
        Do i = 1, ny
          dy(i,izb) = rhs(i,izb)
        EndDo
        If ( iheat > 0 ) dt9(izb) = rhs(ny+1,izb)!*1.0e-9
      EndIf
    EndDo

    If ( idiag >= 6 ) Then
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          izone = izb + szbatch - zb_lo
          !$acc update wait(tid) &
          !$acc host(dy(:,izb),dt9(izb))
          Write(lun_diag,"(a,i5)") 'JAC_SOLVE',izone
          Write(lun_diag,"(14es10.3)") (dy(i,izb),i=1,ny)
          If ( iheat > 0 ) Write(lun_diag,"(es10.3)") dt9(izb)
        EndIf
      EndDo
      Flush(lun_diag)
    EndIf

    !$acc exit data async(tid) &
    !$acc delete(mask)

    stop_timer = xnet_wtime()
    timer_solve = timer_solve + stop_timer

    Return
  End Subroutine jacobian_solve

  Subroutine jacobian_decomp(kstep,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine performs the LU matrix decomposition for the Jacobian.
    !-----------------------------------------------------------------------------------------------
    Use xnet_controls, Only: idiag, lun_diag, nzbatch, szbatch, zb_lo, zb_hi, lzactive, tid
    Use xnet_linalg, Only: LUDecompBatched_GPU, LUDecomp_CPU
    Use xnet_timers, Only: xnet_wtime, start_timer, stop_timer, timer_solve, timer_decmp
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer :: i, j, izb, izone, izb_p, nzmask, istat
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return
    nzmask = count(mask)

    start_timer = xnet_wtime()
    timer_solve = timer_solve - start_timer
    timer_decmp = timer_decmp - start_timer

    ! Calculate the LU decomposition
#if defined(XNET_GPU)
    !$acc serial async(tid) &
    !$acc present(mask,djacp,djac,dindxp,dindx)
    i = 0
    !$acc loop &
    !$acc private(izb_p) &
    !$acc reduction(+:i)
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        i = i + 1
        izb_p = zb_lo + i - 1
        djacp(izb_p) = djac(izb)
        dindxp(izb_p) = dindx(izb)
      EndIf
    EndDo
    !$acc end serial
    call LUDecompBatched_GPU &
      & ( msize, msize, djacp(zb_lo), msize, dindxp(zb_lo), dindxinfo(zb_lo), info(zb_lo), nzmask )
    !$acc update async(tid) &
    !$acc device(info(zb_lo))
#else
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        call LUDecomp_CPU &
          & ( msize, msize, jac(1,1,izb), msize, indx(1,izb), info(izb) )
      EndIf
    EndDo
#endif

    If ( idiag >= 6 ) Then
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          izone = izb + szbatch - zb_lo
          !$acc update wait(tid) &
          !$acc host(jac(:,:,izb))
          Write(lun_diag,"(a3,i5)") 'LUD',izone
          Write(lun_diag,"(14es9.1)") ((jac(i,j,izb),j=1,msize),i=1,msize)
        EndIf
      EndDo
      Flush(lun_diag)
    EndIf

    stop_timer = xnet_wtime()
    timer_solve = timer_solve + stop_timer
    timer_decmp = timer_decmp + stop_timer

    Return
  End Subroutine jacobian_decomp

  Subroutine jacobian_bksub(kstep,yrhs,dy,t9rhs,dt9,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine performs back-substitution for a LU matrix and the RHS vector.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use xnet_controls, Only: idiag, iheat, lun_diag, nzbatch, szbatch, zb_lo, zb_hi, lzactive, &
      & tid
    Use xnet_linalg, Only: LUBksubBatched_GPU, LUBksub_CPU
    Use xnet_timers, Only: xnet_wtime, start_timer, stop_timer, timer_solve, timer_bksub
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep
    Real(dp), Intent(in) :: yrhs(ny,zb_lo:zb_hi)
    Real(dp), Intent(in) :: t9rhs(zb_lo:zb_hi)

    ! Output variables
    Real(dp), Intent(out) :: dy(ny,zb_lo:zb_hi)
    Real(dp), Intent(out) :: dt9(zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer :: i, izb, izone, izb_p, nzmask, istat
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return
    nzmask = count(mask)

    start_timer = xnet_wtime()
    timer_solve = timer_solve - start_timer
    timer_bksub = timer_bksub - start_timer

    !$acc enter data async(tid) &
    !$acc copyin(mask)

    !$acc parallel loop gang async(tid) &
    !$acc present(mask,yrhs,t9rhs,rhs)
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        !$acc loop vector
        Do i = 1, ny
          rhs(i,izb) = yrhs(i,izb)
        EndDo
        If ( iheat > 0 ) rhs(ny+1,izb) = t9rhs(izb)!*1.0e9
      Else
        !$acc loop vector
        Do i = 1, msize
          rhs(i,izb) = 0.0
        EndDo
      EndIf
    EndDo

    ! Solve the LU-decomposed triangular system via back-substitution
#if defined(XNET_GPU)
    !$acc serial async(tid) &
    !$acc present(mask,djacp,djac,drhsp,drhs,dindxp,dindx)
    i = 0
    !$acc loop &
    !$acc private(izb_p) &
    !$acc reduction(+:i)
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        i = i + 1
        izb_p = zb_lo + i - 1
        djacp(izb_p) = djac(izb)
        drhsp(izb_p) = drhs(izb)
        dindxp(izb_p) = dindx(izb)
      EndIf
    EndDo
    !$acc end serial
    call LUBksubBatched_GPU &
      & ( 'N', msize, 1, djacp(zb_lo), msize, dindxp(zb_lo), drhsp(zb_lo), msize, &
      &   dwork(zb_lo), info(zb_lo), nzmask )
#else
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        call LUBksub_CPU &
          & ( 'N', msize, 1, jac(1,1,izb), msize, indx(1,izb), rhs(1,izb), msize, info(izb) )
      EndIf
    EndDo
#endif

    !$acc parallel loop gang async(tid) &
    !$acc present(mask,dy,dt9,rhs)
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        !$acc loop vector
        Do i = 1, ny
          dy(i,izb) = rhs(i,izb)
        EndDo
        If ( iheat > 0 ) dt9(izb) = rhs(ny+1,izb)!*1.0e-9
      EndIf
    EndDo

    If ( idiag >= 6 ) Then
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          izone = izb + szbatch - zb_lo
          !$acc update wait(tid) &
          !$acc host(dy(:,izb),dt9(izb))
          Write(lun_diag,"(a,i5)") 'BKSUB', izone
          Write(lun_diag,"(14es10.3)") (dy(i,izb),i=1,ny)
          If ( iheat > 0 ) Write(lun_diag,"(es10.3)") dt9(izb)
        EndIf
      EndDo
      Flush(lun_diag)
    EndIf

    !$acc exit data async(tid) &
    !$acc delete(mask)

    stop_timer = xnet_wtime()
    timer_solve = timer_solve + stop_timer
    timer_bksub = timer_bksub + stop_timer

    Return
  End Subroutine jacobian_bksub

End Module xnet_jacobian
