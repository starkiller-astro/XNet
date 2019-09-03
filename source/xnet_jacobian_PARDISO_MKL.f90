!***************************************************************************************************
! jacobian_PARDISO.f90 10/18/17
! Sparse Jacobian interface for PARDISO
!
! The routines in this file are used to replace the standard dense Jacobian and associated solver
! with the PARDISO sparse solver package.
!
! The bulk of the computational cost of the network (60-95%) is the solving of the matrix equation.
! Careful selection of the matrix solver is therefore very important to fast computation. For
! networks from a few dozen up to a couple hundred species, hand tuned dense solvers such as those
! supplied by the hardware manufacturer (often LAPACK) or third-parties like NAG, IMSL, etc. are
! fastest. However for larger matrices, sparse solvers are faster. In our experience, PARDISO is
! prefered sparse solver of XNet, because it is generally the fastest and most bulletproof sparse
! package available, with fairly wide distribution as part of the Intel MKL.
!***************************************************************************************************

Module xnet_jacobian
  !-------------------------------------------------------------------------------------------------
  ! Contains data for use in the sparse solver using Compressed-Row Storage (CRS)
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp, i8
  Implicit None

  ! Jacobian data arrays
  Real(dp), Allocatable :: dydotdy(:,:) ! dYdot/dY part of jac
  Real(dp), Allocatable :: tvals(:,:)   ! Jacobian matrix
  Real(dp), Allocatable :: sident(:)    ! Identity matrix
  Integer, Allocatable  :: perm(:)      ! Permutation vector
  Integer, Allocatable  :: ridx(:)      ! Row indices of non-zeroes
  Integer, Allocatable  :: cidx(:)      ! Column indices of non-zeroes
  Integer, Allocatable  :: pb(:)        ! Pointer to first non-zero for each column
  Integer               :: lval         ! Number of non-zeroes from reaction rates only
  Integer               :: nnz          ! Number of non-zeroes (including self-heating terms)
  Integer               :: msize        ! Size of linear system to be solved

  ! nsij(k) is the matrix index in CRS format of the jth reactant in i-reactant reactions for reaction k
  Integer, Allocatable  :: ns11(:), ns21(:), ns22(:), ns31(:), ns32(:), ns33(:)
  Integer, Allocatable  :: ns41(:), ns42(:), ns43(:), ns44(:)
  Integer               :: l1s, l2s, l3s, l4s ! Sizes of ns arrays

  ! PARDISO Attributes
  Integer, Parameter :: mtype = 11 ! PARDISO matrix type code (11 = real and nonsymmetric)
  Integer, Parameter :: nrhs = 1   ! Number of right-hand sides
  Integer, Parameter :: solver = 0 ! Solver method (0 = sparse direct)

  ! Solver controls
  Integer(i8) :: pt(64)            ! PARDISO internal data address pointers
  Integer     :: iparm(64)         ! PARDISO solver parameters
  Real(dp)    :: dparm(64)         ! PARDISO solver parameters for iterative solver
  Integer     :: phase             ! PARDISO execution mode
  Integer     :: msglvl            ! Verbosity of PARDISO diagnostics
  Integer     :: maxfct            ! Maximal number of factors with identical nonzero sparsity structure to keep in memory
  !$omp threadprivate(pt,iparm,dparm,phase)
  Namelist /pardiso_controls/ iparm, dparm

  Interface xnet_pardiso
    Module Procedure mkl_pardiso
  End Interface xnet_pardiso

  Interface xnet_pardisoinit
    Module Procedure mkl_pardisoinit
  End Interface xnet_pardisoinit

Contains

  Subroutine mkl_pardisoinit(error)
    Implicit None
    Integer, Intent(out) :: error
    Interface
      Subroutine pardisoinit(pt,mtype,iparm)
        Use xnet_types, Only: i8
        Implicit None
        Integer(i8), Intent(inout) :: pt(*)
        Integer, Intent(in) :: mtype
        Integer, Intent(inout) :: iparm(*)
      End Subroutine pardisoinit
    End Interface
    Call pardisoinit(pt,mtype,iparm)
    iparm(3) = 1 ! Set # of OpenMP threads for PARDISO to use
    Return
  End Subroutine mkl_pardisoinit

  Subroutine mkl_pardiso(mnum,a,b,x,error)
    Use xnet_types, Only: dp
    Implicit None
    Integer, Intent(in) :: mnum
    Real(dp), Intent(in) :: a(:)
    Real(dp), Intent(inout) :: b(:)
    Real(dp), Intent(out) :: x(:)
    Integer, Intent(out) :: error
    Interface
      Subroutine pardiso(pt,maxfct,mnum,mtype,phase,n,a,ia,ja,perm,nrhs,iparm,msglvl,b,x,error)
        Use xnet_types, Only: dp, i8
        Implicit None
        Integer(i8), Intent(inout) :: pt(*)
        Integer, Intent(in) :: maxfct
        Integer, Intent(in) :: mnum
        Integer, Intent(in) :: mtype
        Integer, Intent(in) :: phase
        Integer, Intent(in) :: n
        Real(dp), Intent(in) :: a(*)
        Integer, Intent(in) :: ia(*)
        Integer, Intent(in) :: ja(*)
        Integer, Intent(inout) :: perm(*)
        Integer, Intent(in) :: nrhs
        Integer, Intent(inout) :: iparm(*)
        Integer, Intent(in) :: msglvl
        Real(dp), Intent(inout) :: b(*)
        Real(dp), Intent(out) :: x(*)
        Integer, Intent(out) :: error
      End Subroutine pardiso
    End Interface
    Call pardiso(pt,maxfct,mnum,mtype,phase,msize,a,pb,cidx,perm,nrhs,iparm,msglvl,b,x,error)
    Return
  End Subroutine mkl_pardiso

  Subroutine read_jacobian_data(data_dir)
    !-----------------------------------------------------------------------------------------------
    ! Reads in data necessary to use sparse solver and initializes the Jacobian data.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use reaction_data, Only: la, le, n11, n21, n22, n31, n32, n33, n41, n42, n43, n44
    Use xnet_controls, Only: idiag, iheat, lun_diag, nzbatchmx, nzevolve, zb_lo, zb_hi
    Use xnet_parallel, Only: parallel_bcast, parallel_IOProcessor
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: data_dir

    ! Local variables
    Integer, Allocatable :: ridxo(:), cidxo(:), pbo(:)
    Integer :: i, i0, i1, la1, le1, la2, le2, la3, le3, la4, le4, j1, l1, l2, l3, l4
    Integer :: rstart, rstarto, rend, rendo
    Integer :: ierr, lun_sparse, lun_solver

    If ( parallel_IOProcessor() ) Then
      Open(newunit=lun_sparse, file=trim(data_dir)//"/sparse_ind", status='old', form='unformatted')
      Read(lun_sparse) lval
    EndIf
    Call parallel_bcast(lval)

    ! Determine number of non-zeroes
    If ( iheat > 0 ) Then
      nnz = lval + 2*ny + 1
      msize = ny+1
    Else
      nnz = lval
      msize = ny
    EndIf

    ! Allocate, read, and broadcast CRS arrays
    Allocate (ridx(nnz),cidx(nnz),sident(nnz),pb(msize+1))
    If ( parallel_IOProcessor() ) Then
      If ( iheat > 0 ) Then
        Allocate (ridxo(lval),cidxo(lval),pbo(ny+1))
        Read(lun_sparse) ridxo, cidxo, pbo

        ! Add indices for self-heating
        pb(1) = pbo(1)
        Do i0 = 1, ny

          ! Shift row pointers to adjust for additional column
          pb(i0+1) = pbo(i0+1) + i0
          rstarto = pbo(i0)
          rendo = pbo(i0+1) - 1
          rstart  = pb(i0)
          rend  = pb(i0+1) - 1

          ! Extra column indices
          cidx(rstart:rend-1) = cidxo(rstarto:rendo)
          cidx(rend) = ny+1
          ridx(rstart:rend-1) = ridxo(rstarto:rendo)
          ridx(rend) = i0
        EndDo
        Deallocate(ridxo,cidxo,pbo)

        ! Extra row indices
        pb(msize+1) = nnz + 1
        rstart = pb(msize)
        rend = pb(msize+1) - 1
        Do i0 = 1, ny
          cidx(rstart+i0-1) = i0
          ridx(rstart+i0-1) = ny+1
        EndDo

        ! dT9dot/dT9 term
        cidx(nnz) = ny + 1
        ridx(nnz) = ny + 1
      Else
        Read(lun_sparse) ridx, cidx, pb
      EndIf
      Read(lun_sparse) l1s, l2s, l3s, l4s
    EndIf
    Call parallel_bcast(ridx)
    Call parallel_bcast(cidx)
    Call parallel_bcast(pb)
    Call parallel_bcast(l1s)
    Call parallel_bcast(l2s)
    Call parallel_bcast(l3s)
    Call parallel_bcast(l4s)

    ! Allocate, read, and broadcast arrays for direct sparse representation Jacobian build
    Allocate (ns11(l1s))
    Allocate (ns21(l2s),ns22(l2s))
    Allocate (ns31(l3s),ns32(l3s),ns33(l3s))
    Allocate (ns41(l4s),ns42(l4s),ns43(l4s),ns44(l4s))

    ! Rebuild arrays mapping reaction rates to locations in CRS Jacobian
    If ( parallel_IOProcessor() ) Then
      If ( iheat > 0 ) Then
        ns11 = 0 ; ns21 = 0 ; ns22 = 0 ; ns31 = 0 ; ns32 = 0 ; ns33 = 0
        ns41 = 0 ; ns42 = 0 ; ns43 = 0 ; ns44 = 0
        Do i0 = 1, ny
          la1 = la(1,i0) ; la2 = la(2,i0) ; la3 = la(3,i0) ; la4 = la(4,i0)
          le1 = le(1,i0) ; le2 = le(2,i0) ; le3 = le(3,i0) ; le4 = le(4,i0)
          Do j1 = la1, le1
            l1 = n11(j1)
            Do i1 = 1, nnz
              If ( cidx(i1) == l1 .and. ridx(i1) == i0 ) ns11(j1) = i1
            EndDo
          EndDo
          Do j1 = la2, le2
            l1 = n21(j1) ; l2 = n22(j1)
            Do i1 = 1, nnz
              If ( cidx(i1) == l1 .and. ridx(i1) == i0 ) ns21(j1) = i1
              If ( cidx(i1) == l2 .and. ridx(i1) == i0 ) ns22(j1) = i1
            EndDo
          EndDo
          Do j1 = la3, le3
            l1 = n31(j1) ; l2 = n32(j1) ; l3 = n33(j1)
            Do i1 = 1, nnz
              If ( cidx(i1) == l1 .and. ridx(i1) == i0 ) ns31(j1) = i1
              If ( cidx(i1) == l2 .and. ridx(i1) == i0 ) ns32(j1) = i1
              If ( cidx(i1) == l3 .and. ridx(i1) == i0 ) ns33(j1) = i1
            EndDo
          EndDo
          Do j1 = la4, le4
            l1 = n41(j1) ; l2 = n42(j1) ; l3 = n43(j1) ; l4 = n44(j1)
            Do i1 = 1, nnz
              If ( cidx(i1) == l1 .and. ridx(i1) == i0 ) ns41(j1) = i1
              If ( cidx(i1) == l2 .and. ridx(i1) == i0 ) ns42(j1) = i1
              If ( cidx(i1) == l3 .and. ridx(i1) == i0 ) ns43(j1) = i1
              If ( cidx(i1) == l4 .and. ridx(i1) == i0 ) ns44(j1) = i1
            EndDo
          EndDo
        EndDo
      Else
        Read(lun_sparse) ns11, ns21, ns22
        Read(lun_sparse) ns31
        Read(lun_sparse) ns32
        Read(lun_sparse) ns33
        Read(lun_sparse) ns41
        Read(lun_sparse) ns42
        Read(lun_sparse) ns43
        Read(lun_sparse) ns44
      EndIf
      Close(lun_sparse)
    EndIf
    Call parallel_bcast(ns11)
    Call parallel_bcast(ns21)
    Call parallel_bcast(ns22)
    Call parallel_bcast(ns31)
    Call parallel_bcast(ns32)
    Call parallel_bcast(ns33)
    Call parallel_bcast(ns41)
    Call parallel_bcast(ns42)
    Call parallel_bcast(ns43)
    Call parallel_bcast(ns44)

    ! Build a CRS format version of the identity matrix
    Where ( ridx == cidx )
      sident = 1.0
    ElseWhere
      sident = 0.0
    EndWhere

    ! Read and broadcast user-defined PARDISO controls
    If ( parallel_IOProcessor() ) Then

      call xnet_pardisoinit(ierr)

      ! Override defaults with user-defined inputs
      Open(newunit=lun_solver, file="sparse_controls.nml", action='read', status='old', iostat=ierr)
      If ( ierr == 0 ) Then
        Read(lun_solver,nml=pardiso_controls)
        Close(lun_solver)
      EndIf
    EndIf
    Call parallel_bcast(iparm)
    Call parallel_bcast(dparm)

    ! Set level of PARDISO verbosity
    If ( idiag >= 5 ) Then
      msglvl = 1
    Else
      msglvl = 0
    EndIf

    ! Set number of factorizations for PARDISO to store with identical sparsity structure
    maxfct = nzbatchmx

    Allocate (dydotdy(nnz,nzevolve),tvals(nnz,nzevolve))

    !$omp parallel default(shared) copyin(pt,iparm,dparm)

    ! Initialize work arrays
    Allocate (perm(msize))
    perm = 0

    !$omp end parallel

    Return
  End Subroutine read_jacobian_data

  Subroutine jacobian_scale(diag,mult,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This augments a previously calculation Jacobian matrix by multiplying all elements by mult and
    ! adding diag to the diagonal elements.
    !-----------------------------------------------------------------------------------------------
    Use xnet_controls, Only: zb_lo, zb_hi, lzactive
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: diag(zb_lo:zb_hi), mult(zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer :: i, i0, izb
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        tvals(:,izb) = mult(izb) * dydotdy(:,izb) + diag(izb) * sident(:)
      EndIf
    EndDo
    
    Return
  End Subroutine jacobian_scale

  Subroutine jacobian_build(diag_in,mult_in,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the reaction Jacobian matrix, dYdot/dY, and augments by multiplying all
    ! elements by mult and adding diag to the diagonal elements.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny, mex
    Use reaction_data, Only: a1, a2, a3, a4, b1, b2, b3, b4, la, le, mu1, mu2, mu3, mu4, n11, n21, &
      & n22, n31, n32, n33, n41, n42, n43, n44, dcsect1dt9, dcsect2dt9, dcsect3dt9, dcsect4dt9, nan
    Use xnet_abundances, Only: yt
    Use xnet_conditions, Only: cv
    Use xnet_controls, Only: iheat, idiag, ktot, lun_diag, szbatch, zb_lo, zb_hi, lzactive
    Use xnet_timers, Only: xnet_wtime, start_timer, stop_timer, timer_jacob
    Use xnet_types, Only: dp
    Implicit None

    ! Optional variables
    Real(dp), Optional, Intent(in) :: diag_in(zb_lo:zb_hi), mult_in(zb_lo:zb_hi)
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer :: i, j, i0, j1, izb, izone
    Integer :: la1, le1, la2, le2, la3, le3, la4, le4
    Integer :: i11, i21, i22, i31, i32, i33, i41, i42, i43, i44
    Integer :: is11, is21, is22, is31, is32, is33, is41, is42, is43, is44
    Real(dp) :: s1, s2, s3, s4, r1, r2, r3, r4
    Real(dp) :: y11, y21, y22, y31, y32, y33, y41, y42, y43, y44
    Real(dp) :: dt9dotdy(msize), dr1dt9, dr2dt9, dr3dt9, dr4dt9
    Real(dp) :: diag(zb_lo:zb_hi), mult(zb_lo:zb_hi)
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    start_timer = xnet_wtime()
    timer_jacob = timer_jacob - start_timer

    ! Build the Jacobian
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) THen
        dydotdy(:,izb) = 0.0
        Do j1 = 1, nan(1)
          r1 = b1(j1,izb)
          is11 = ns11(j1)
          dydotdy(is11,izb) = dydotdy(is11,izb) + r1
        EndDo
        Do j1 = 1, nan(2)
          r2 = b2(j1,izb)
          is21 = ns21(j1)
          is22 = ns22(j1)
          i21 = n21(j1)
          i22 = n22(j1)
          y21 = yt(i21,izb)
          y22 = yt(i22,izb)
          dydotdy(is21,izb) = dydotdy(is21,izb) + r2 * y22
          dydotdy(is22,izb) = dydotdy(is22,izb) + r2 * y21
        EndDo
        Do j1 = 1, nan(3)
          r3 = b3(j1,izb)
          is31 = ns31(j1)
          is32 = ns32(j1)
          is33 = ns33(j1)
          i31 = n31(j1)
          i32 = n32(j1)
          i33 = n33(j1)
          y31 = yt(i31,izb)
          y32 = yt(i32,izb)
          y33 = yt(i33,izb)
          dydotdy(is31,izb) = dydotdy(is31,izb) + r3 * y32 * y33
          dydotdy(is32,izb) = dydotdy(is32,izb) + r3 * y31 * y33
          dydotdy(is33,izb) = dydotdy(is33,izb) + r3 * y31 * y32
        EndDo
        Do j1 = 1, nan(4)
          r4 = b4(j1,izb)
          is41 = ns41(j1)
          is42 = ns42(j1)
          is43 = ns43(j1)
          is44 = ns44(j1)
          i41 = n41(j1)
          i42 = n42(j1)
          i43 = n43(j1)
          i44 = n44(j1)
          y41 = yt(i41,izb)
          y42 = yt(i42,izb)
          y43 = yt(i43,izb)
          y44 = yt(i44,izb)
          dydotdy(is41,izb) = dydotdy(is41,izb) + r4 * y42 * y43 * y44
          dydotdy(is42,izb) = dydotdy(is42,izb) + r4 * y43 * y44 * y41
          dydotdy(is43,izb) = dydotdy(is43,izb) + r4 * y44 * y41 * y42
          dydotdy(is44,izb) = dydotdy(is44,izb) + r4 * y41 * y42 * y43
        EndDo
      EndIf
    EndDo

    If ( iheat > 0 ) Then
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) THen
          Do i0 = 1, ny
            la1 = la(1,i0)
            la2 = la(2,i0)
            la3 = la(3,i0)
            la4 = la(4,i0)
            le1 = le(1,i0)
            le2 = le(2,i0)
            le3 = le(3,i0)
            le4 = le(4,i0)
            s1 = 0.0
            Do j1 = la1, le1
              dr1dt9 = a1(j1)*dcsect1dt9(mu1(j1),izb)
              i11 = n11(j1)
              y11 = yt(i11,izb)
              s1 = s1 + dr1dt9 * y11
            EndDo
            s2 = 0.0
            Do j1 = la2, le2
              dr2dt9 = a2(j1)*dcsect2dt9(mu2(j1),izb)
              i21 = n21(j1)
              i22 = n22(j1)
              y21 = yt(i21,izb)
              y22 = yt(i22,izb)
              s2 = s2 + dr2dt9 * y21 * y22
            EndDo
            s3 = 0.0
            Do j1 = la3, le3
              dr3dt9 = a3(j1)*dcsect3dt9(mu3(j1),izb)
              i31 = n31(j1)
              i32 = n32(j1)
              i33 = n33(j1)
              y31 = yt(i31,izb)
              y32 = yt(i32,izb)
              y33 = yt(i33,izb)
              s3 = s3 + dr3dt9 * y31 * y32 * y33
            EndDo
            s4 = 0.0
            Do j1 = la4, le4
              dr4dt9 = a4(j1)*dcsect4dt9(mu4(j1),izb)
              i41 = n41(j1)
              i42 = n42(j1)
              i43 = n43(j1)
              i44 = n44(j1)
              y41 = yt(i41,izb)
              y42 = yt(i42,izb)
              y43 = yt(i43,izb)
              y44 = yt(i44,izb)
              s4 = s4 + dr4dt9 * y41 * y42 * y43 * y44
            EndDo
            dydotdy(pb(i0+1)-1,izb) = s1 + s2 + s3 + s4
          EndDo

          dt9dotdy = 0.0
          Do i0 = 1, ny
            Do j1 = pb(i0), pb(i0+1)-1
              dt9dotdy(cidx(j1)) = dt9dotdy(cidx(j1)) + mex(i0)*dydotdy(j1,izb)
            EndDo
          EndDo
          dydotdy(pb(ny+1):nnz,izb) = -dt9dotdy(:) / cv(izb)

        EndIf
      EndDo
    EndIf

    ! Apply the externally provided factors
    If ( present(diag_in) ) Then
      diag = diag_in
    Else
      diag = 0.0
    EndIf
    If ( present(mult_in) ) Then
      mult = mult_in
    Else
      mult = 1.0
    EndIf
    Call jacobian_scale(diag,mult,mask_in = mask)

    If ( idiag >= 6 ) Then
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) THen
          izone = izb + szbatch - zb_lo
          Write(lun_diag,"(a9,i5,2es14.7)") 'JAC_BUILD',izone,diag(izb),mult(izb)
          Write(lun_diag,"(14es9.1)") (tvals(i,izb),i=1,nnz)
        EndIf
      EndDo
    EndIf

    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        ktot(3,izb) = ktot(3,izb) + 1
      EndIf
    EndDo

    stop_timer = xnet_wtime()
    timer_jacob = timer_jacob + stop_timer

    Return
  End Subroutine jacobian_build

  Subroutine jacobian_solve(kstep,yrhs,dy,t9rhs,dt9,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine solves the system of equations composed of the Jacobian and RHS vector.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use xnet_controls, Only: idiag, iheat, lun_diag, szbatch, zb_lo, zb_hi, lzactive
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
    Integer :: i, izb, izone
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    Call jacobian_decomp(kstep,mask_in = mask)
    Call jacobian_bksub(kstep,yrhs,dy,t9rhs,dt9,mask_in = mask)

    If ( idiag >= 5 ) Then
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          izone = izb + szbatch - zb_lo
          Write(lun_diag,"(a)") 'JAC_SOLVE', izone
          Write(lun_diag,"(14es10.3)") (dy(i,izb),i=1,ny)
          If ( iheat > 0 ) Write(lun_diag,"(es10.3)") dt9(izb)
        EndIf
      EndDo
    EndIf

    Return
  End Subroutine jacobian_solve

  Subroutine jacobian_decomp(kstep,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine performs the LU matrix decomposition for the Jacobian.
    !-----------------------------------------------------------------------------------------------
    Use xnet_controls, Only: kitmx, kmon, lun_stdout, zb_lo, zb_hi, lzactive
    Use xnet_timers, Only: xnet_wtime, start_timer, stop_timer, timer_solve, timer_decmp
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer :: i, izb, izone, err
    Real(dp) :: rhs(msize), dx(msize)
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    start_timer = xnet_wtime()
    timer_solve = timer_solve - start_timer
    timer_decmp = timer_decmp - start_timer

    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then

        ! Check if this is the first step, or if matrix needs to be reanalyzed
        If ( kstep == 1 .or. kmon(2,izb) > kitmx ) Then
          phase = 12
        Else
          phase = 22
        EndIf
        Call xnet_pardiso(izb,tvals(:,izb),rhs,dx,err)
        If ( err /= 0 ) Write(lun_stdout,*) 'PARDISO error ',err,' phase=',phase
      EndIf
    EndDo

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
    Use xnet_controls, Only: idiag, iheat, lun_diag, lun_stdout, szbatch, zb_lo, zb_hi, lzactive
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
    Real(dp) :: rhs(msize), dx(msize)
    Integer :: i, izb, izone, err
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    start_timer = xnet_wtime()
    timer_solve = timer_solve - start_timer
    timer_bksub = timer_bksub - start_timer

    ! Solve linear system using existing factorization
    phase = 33
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        rhs(1:ny) = yrhs(:,izb)
        If ( iheat > 0 ) rhs(ny+1) = t9rhs(izb)
        Call xnet_pardiso(izb,tvals(:,izb),rhs,dx,err)
        If ( err /= 0 ) Write(lun_stdout,*) 'PARDISO error ',err,' phase=',phase
        dy(:,izb) = dx(1:ny)
        If ( iheat > 0 ) dt9(izb) = dx(ny+1)
      EndIf
    EndDo

    If ( idiag >= 5 ) Then
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          izone = izb + szbatch - zb_lo
          Write(lun_diag,"(a,i5)") 'BKSUB', izone
          Write(lun_diag,"(14es10.3)") (dy(i,izb),i=1,ny)
          If ( iheat > 0 ) Write(lun_diag,"(es10.3)") dt9(izb)
        EndIf
      EndDo
    EndIf

    stop_timer = xnet_wtime()
    timer_solve = timer_solve + stop_timer
    timer_bksub = timer_bksub + stop_timer

    Return
  End Subroutine jacobian_bksub

End Module xnet_jacobian
