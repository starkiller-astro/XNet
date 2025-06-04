!***************************************************************************************************
! jacobian_MA48.f90 10/18/17
! Sparse Jacobian interface for MA48
!
! The routines in this file are used to replace the standard dense Jacobian and associated solver
! with the HSL MA48 sparse solver package.
!
! The bulk of the computational cost of the network (60-95%) is the solving of the matrix equation.
! Careful selection of the matrix solver is therefore very important to fast computation. For
! networks from a few dozen up to a couple hundred species, hand tuned dense solvers such as those
! supplied by the hardware manufacturer (often LAPACK) or third-parties like NAG, IMSL, etc. are
! fastest. However for larger matrices, sparse solvers are faster. MA48 is a solver from HSL, which
! is available under an academic license. We find it to be faster than PARDISO for matrices of size
! 100 < ny < 1000, but slower at larger sizes, but the availablilty of the MA48 source makes it
! valuable for some applications.
!
! Reference:
! HSL(2013). A collection of Fortran codes for large scale scientific computation.
!   http://www.hsl.rl.ac.uk
!***************************************************************************************************

Module xnet_jacobian
  !-------------------------------------------------------------------------------------------------
  ! Contains data for use in the sparse solver.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None

  ! Jacobian data arrays
  Real(dp), Allocatable :: dydotdy(:,:) ! dYdot/dY part of jac
  Real(dp), Allocatable :: tvals(:,:)   ! Jacobian matrix
  Real(dp), Allocatable :: jac(:,:,:)   ! Jacobian matrix
  Real(dp), Allocatable :: jac_cond(:)  ! Jacobian matrix
  Real(dp), Allocatable :: sident(:)    ! Identity matrix
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

  ! MA48 internal working arrays and solver controls
  Real(dp), Allocatable :: vals(:,:), wB(:), wC(:)
  Integer, Allocatable  :: irn(:,:), jcn(:,:)
  Integer, Allocatable  :: keep(:,:), iwA(:), iwB(:), iwC(:)
  Integer, Allocatable :: jobA(:), jobB(:), jobC(:)
  Integer :: lia
  Integer :: icntl(20), info(20)
  Real(dp) :: cntl(10), rinfo(10), maxerr
  !$omp threadprivate(vals,lia,jcn,irn,wB,wC,iwA,iwB,iwC,keep,info,rinfo,cntl,icntl)
  Namelist /ma48_controls/ icntl, cntl, maxerr

  ! Temporary copies for reallocating
  Real(dp), Allocatable ::  vals0(:,:)
  Integer, Allocatable :: irn0(:,:), jcn0(:,:)
  Integer :: lia0
  !$omp threadprivate(vals0,jcn0,irn0,lia0)

  ! Some other solver parameters
  Logical, Parameter :: trans = .false. ! Flag for solving the transposed system
  Integer, Parameter :: kbksubmx = 3, kdecompmx = 5 ! Max loop counts

  Interface
    Subroutine mmwrite(ounit,rep,field,symm,rows,cols,nnz,indx,jndx,ival,rval,cval)
      Integer :: ival(*)
      Double Precision :: rval(*)
      Complex :: cval(*)
      Integer :: indx(*)
      Integer :: jndx(*)
      Integer :: ounit, rows, cols, nnz
      Character(len=*) :: rep, field, symm
    End Subroutine mmwrite
  End Interface

Contains

  Subroutine read_jacobian_data(data_dir)
    !-----------------------------------------------------------------------------------------------
    ! Reads in data necessary to use sparse solver and initializes the Jacobian data.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use xnet_controls, Only: idiag, iheat, lun_diag, nzevolve, zb_lo, zb_hi
    Use xnet_parallel, Only: parallel_bcast, parallel_IOProcessor
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: data_dir

    ! Local variables
    Integer :: i, ierr, lun_sparse, lun_solver

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
      Read(lun_sparse) ridx(1:lval), cidx(1:lval), pb(1:ny+1)
      If ( iheat > 0 ) Then
        ! Add indices for self-heating
        Do i = 1, ny
          cidx(i+lval) = ny + 1     ! Extra column (dYdot/dT9)
          ridx(i+lval) = i
          cidx(i+lval+ny) = i       ! Extra row (dT9dot/dY)
          ridx(i+lval+ny) = ny + 1
        EndDo
        cidx(nnz) = ny + 1          ! dT9dot/dT9 term
        ridx(nnz) = ny + 1
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
    If ( parallel_IOProcessor() ) Then
      Read(lun_sparse) ns11, ns21, ns22
      Read(lun_sparse) ns31
      Read(lun_sparse) ns32
      Read(lun_sparse) ns33
      Read(lun_sparse) ns41
      Read(lun_sparse) ns42
      Read(lun_sparse) ns43
      Read(lun_sparse) ns44
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

    ! Build a compressed row format version of the identity matrix
    Where ( ridx == cidx )
      sident = 1.0
    ElseWhere
      sident = 0.0
    EndWhere

    ! Read and broadcast user-defined MA48 controls
    If ( parallel_IOProcessor() ) Then

      ! Set the value for the maximum allowed error in the call to MA48CD
      maxerr = 1.0d-11

      ! Set default values for MA48 control parameters
      Call MA48ID(cntl,icntl)

      ! Set level of MA48 verbosity (default=2)
      If ( idiag >= 5 ) Then
        icntl(3) = 3
      EndIf

      ! Set block column size (default=32)
      If ( msize < 32 ) Then
        icntl(5) = msize
      ElseIf ( msize < 512 ) Then
        icntl(5) = 32
      ElseIf ( msize < 2048 ) Then
        icntl(5) = 128
      Else
        icntl(5) = 256
      EndIf

      ! Unit numbers for MA48 error/diagnostic output will be set to XNet diagnostic file by default,
      ! but we put in random values here so we can later check for user-defined input
      icntl(1) = -99
      icntl(2) = -99

      ! Override defaults with user-defined inputs
      Open(newunit=lun_solver, file="sparse_controls.nml", action='read', status='old', iostat=ierr)
      If ( ierr == 0 ) Then
        Read(lun_solver,nml=ma48_controls)
        Close(lun_solver)
      EndIf
    EndIf
    Call parallel_bcast(icntl)
    Call parallel_bcast(cntl)
    Call parallel_bcast(maxerr)

    Allocate (jac(msize,msize,nzevolve))
    Allocate (jac_cond(nzevolve))
    Allocate (dydotdy(nnz,nzevolve),tvals(nnz,nzevolve))
    Allocate (jobA(nzevolve))
    Allocate (jobB(nzevolve))
    Allocate (jobC(nzevolve))

    !$omp parallel default(shared) copyin(cntl,icntl)

    ! These are variables to be used by the MA48 solver
    lia = 8*nnz
    Allocate (vals(lia,zb_lo:zb_hi),jcn(lia,zb_lo:zb_hi),irn(lia,zb_lo:zb_hi))
    Allocate (wB(msize),wC(4*msize))
    Allocate (iwA(9*msize),iwB(4*msize),iwC(msize))
    If (icntl(8) == 0) Then
      Allocate (keep(6*msize + 4*msize/icntl(6) + 7 - max(msize/icntl(6),1),zb_lo:zb_hi))
    Else
      Allocate (keep(6*msize + 4*msize/icntl(6) + 7,zb_lo:zb_hi))
    EndIf

    ! Initialize work arrays
    jac = 0.0
    jac_cond = 0.0
    vals = 0.0
    wB = 0.0
    wC = 0.0
    jcn = 0
    irn = 0
    iwA = 0
    iwB = 0
    iwC = 0
    keep = 0
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
    Use matrix_util, Only: matrix_cond
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
        tvals(:,izb) = mult(izb) * dydotdy(:,izb) + diag(izb) * sident
        jac(:,:,izb) = 0.0_dp
        Do i = 1, nnz
          jac(ridx(i),cidx(i),izb) = tvals(i,izb)
        EndDo
        Call matrix_cond(msize,jac(:,:,izb),jac_cond(izb))
      EndIf
    EndDo
    
    Return
  End Subroutine jacobian_scale

  Subroutine jacobian_build(diag_in,mult_in,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the reaction Jacobian matrix, dYdot/dY, and augments by multiplying
    ! all elements by mult and adding diag to the diagonal elements.
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
    Real(dp) :: dt9dotdy(ny), dr1dt9, dr2dt9, dr3dt9, dr4dt9
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
            dydotdy(lval+i0,izb) = s1 + s2 + s3 + s4
          EndDo

          dt9dotdy = 0.0
          Do j1 = 1, lval
            dt9dotdy(cidx(j1)) = dt9dotdy(cidx(j1)) + mex(ridx(j1))*dydotdy(j1,izb)
          EndDo
          dydotdy(lval+ny+1:lval+2*ny,izb) = -dt9dotdy / cv(izb)

          s1 = 0.0
          Do i0 = 1, ny
            s1 = s1 + mex(i0)*dydotdy(lval+i0,izb)
          EndDo
          dydotdy(nnz,izb) = -s1 / cv(izb)
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
    Call jacobian_scale(diag,mult,mask_in = mask_in)

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

  Subroutine jacobian_write_matrix(mask_in)
    Use xnet_controls, Only: zb_lo, zb_hi, lzactive
    Implicit None

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer :: ival(1)
    Complex :: cval(1)
    Integer :: izb
    Integer :: lun_matrix
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        !XDIR XUPDATE XWAIT(tid) &
        !XDIR XHOST(tvals(:,izb))
        Open(newunit=lun_matrix, file="matrix.mtx")
        Call mmwrite(lun_matrix,'coordinate','real','general',msize,msize,nnz, &
          & ridx,cidx,ival,tvals(:,izb),cval)
        Close(lun_matrix)
      EndIf
    EndDo

    Return
  End Subroutine jacobian_write_matrix

  Subroutine jacobian_write_rhs(yrhs,t9rhs,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine solves the system of equations composed of the Jacobian and RHS vector.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use xnet_controls, Only: iheat, zb_lo, zb_hi, lzactive, tid
    Use xnet_types, Only: dp
    Use xnet_util, Only: xnet_terminate
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: yrhs(ny,zb_lo:zb_hi)
    Real(dp), Intent(in) :: t9rhs(zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Real(dp) :: rhs(msize)
    Integer :: ival(1)
    Complex :: cval(1)
    Integer :: izb
    Integer :: lun_rhs
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        !XDIR XUPDATE XWAIT(tid) &
        !XDIR XHOST(yrhs(:,izb),t9rhs(izb))
        rhs(1:ny) = yrhs(:,izb)
        If ( iheat > 0 ) rhs(ny+1) = t9rhs(izb)
        Open(newunit=lun_rhs, file="rhs.mtx")
        Call mmwrite(lun_rhs,'array','real','general',msize,1,msize, &
          & ridx,cidx,ival,rhs,cval)
        Close(lun_rhs)
        Call xnet_terminate('Wrote matrix and RHS')
      EndIf
    EndDo

    Return
  End Subroutine jacobian_write_rhs

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

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Output variables
    Real(dp), Intent(out) :: dy(ny,zb_lo:zb_hi)
    Real(dp), Intent(out) :: dt9(zb_lo:zb_hi)

    ! Local variables
    Integer :: i, izb, izone
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    Call jacobian_decomp(kstep,mask_in = mask_in)
    Call jacobian_bksub(kstep,yrhs,dy,t9rhs,dt9,mask_in = mask_in)

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
    Use xnet_controls, Only: lun_diag, zb_lo, zb_hi, lzactive
    Use xnet_timers, Only: xnet_wtime, start_timer, stop_timer, timer_solve, timer_decmp
    Use xnet_types, Only: dp
    Use xnet_util, Only: xnet_terminate
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer :: i, izb, izone, j, kdecomp, jdecomp
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
        If ( kstep == 1 .or. jobB(izb) == 1 ) Then

          ! Try to restrict pivoting to the diagonal
          jobA(izb) = 3

          ! Perform symbolic analysis
          Do kdecomp = 1, kdecompmx
            jcn(1:nnz,izb) = cidx
            irn(1:nnz,izb) = ridx
            vals(1:nnz,izb) = tvals(:,izb)
            info = 0
            rinfo = 0.0
            Call MA48AD(msize,msize,nnz,jobA(izb),lia,vals(:,izb),irn(:,izb),jcn(:,izb), &
              & keep(:,izb),cntl,icntl,iwA,info,rinfo)
            If ( info(1) == 0 .and. info(4) <= lia .and. info(2) <= 10 ) Then
              jobB(izb) = 1 ! If analysis is successful, proceed to factorization
              Exit
            ElseIf ( kdecomp == kdecompmx ) Then
              Write(lun_diag,"(1x,a,i5,a,i2,a,i2)") 'kstep=',kstep,', Error during MA48AD. kdecomp=',kdecomp,', info(1)=',info(1)
              Call xnet_terminate('Error during MA48AD',info(1))
            Else
              Write(lun_diag,"(1x,a,i5,a,i2,a,i2)") 'kstep=',kstep,', Warning during MA48AD. kdecomp=',kdecomp,', info(1)=',info(1)

              ! Address any error codes
              If ( info(1) >= 4 .and. jobA(izb) == 3 ) Then
                Write(lun_diag,"(1x,a,i3,a,i7)") 'Not possible to choose all pivots from diagonal'
                jobA(izb) = 1
              ElseIf ( info(1) == -3 .or. info(4) > lia ) Then
                Write(lun_diag,"(1x,a,i7,a,i7)") 'Reallocating MA48 arrays: lia=',lia,' < info(4)=',max(info(3),info(4))
                lia0 = lia
                lia = int(1.2*max(info(3),info(4)))

                Allocate (jcn0(lia,zb_lo:zb_hi))
                jcn0(1:lia0,:) = jcn
                Call move_alloc(jcn0,jcn)
                jcn(:,izb) = 0

                Allocate (irn0(lia,zb_lo:zb_hi))
                irn0(1:lia0,:) = irn
                Call move_alloc(irn0,irn)
                irn(:,izb) = 0

                Allocate (vals0(lia,zb_lo:zb_hi))
                vals0(1:lia0,:) = vals
                Call move_alloc(vals0,vals)
                vals(:,izb) = 0.0
              ElseIf ( info(2) > 10 ) Then
                Write(lun_diag,"(1x,a,i3,a,i7)") 'Reallocating MA48 arrays: info(2)=',info(2),', info(4)=',max(info(3),info(4))
                lia0 = lia
                lia = int(2.0*max(info(3),info(4)))

                Allocate (jcn0(lia,zb_lo:zb_hi))
                jcn0(1:lia0,:) = jcn
                Call move_alloc(jcn0,jcn)
                jcn(:,izb) = 0

                Allocate (irn0(lia,zb_lo:zb_hi))
                irn0(1:lia0,:) = irn
                Call move_alloc(irn0,irn)
                irn(:,izb) = 0

                Allocate (vals0(lia,zb_lo:zb_hi))
                vals0(1:lia0,:) = vals
                Call move_alloc(vals0,vals)
                vals(:,izb) = 0.0
              EndIf
            EndIf
          EndDo
        Else
          vals(1:nnz,izb) = tvals(:,izb)
        EndIf

        ! Perform numerical decomposition using previously determined symbolic decomposition
        Do kdecomp = 1, kdecompmx
          info = 0
          rinfo = 0.0
          Call MA48BD(msize,msize,nnz,jobB(izb),lia,vals(:,izb),irn(:,izb),jcn(:,izb), &
            & keep(:,izb),cntl,icntl,wB,iwB,info,rinfo)
          If ( info(1) == 0 .and. info(4) <= lia .and. info(6) == 0 ) Then
            If ( icntl(8) == 0 ) Then
              jobB(izb) = 2 ! Unless using special case of icntl(8)/=0, use the "fast" MA48BD call
            Else
              jobB(izb) = 3 ! For the case of icntl(8)/=0, use the "intermediate" MA48BD all
            EndIf
            If ( kstep /= 0 ) jobC(izb) = 1
            Exit
          ElseIf ( kdecomp == kdecompmx ) Then
            Write(lun_diag,"(1x,a,i5,a,i2,a,i2)") 'kstep=',kstep,', Error during MA48BD. kdecomp=',kdecomp,', info(1)=',info(1)
            Call xnet_terminate('Error during MA48BD',info(1))
          Else
            Write(lun_diag,"(1x,a,i5,a,i2,a,i2)") 'kstep=',kstep,', Warning during MA48BD. kdecomp=',kdecomp,', info(1)=',info(1)

            ! Perform check to see if an error occured in MA48BD.
            ! Most likely a singularity error caused by incorrect symbolic matrix, so run MA48AD again.
            ! Also, make sure that workspaces are large enough.
            If ( info(1) == -3 .or. info(4) > lia ) Then
              Write(lun_diag,"(1x,a,i7,a,i7)") 'Reallocating MA48 arrays: lia=',lia,' < info(4)=',info(4)
              lia0 = lia
              lia = int(1.2*max(info(3),info(4)))

              Allocate (jcn0(lia,zb_lo:zb_hi))
              jcn0(1:lia0,:) = jcn
              Call move_alloc(jcn0,jcn)
              jcn(:,izb) = 0

              Allocate (irn0(lia,zb_lo:zb_hi))
              irn0(1:lia0,:) = irn
              Call move_alloc(irn0,irn)
              irn(:,izb) = 0

              Allocate (vals0(lia,zb_lo:zb_hi))
              vals0(1:lia0,:) = vals
              Call move_alloc(vals0,vals)
              vals(:,izb) = 0.0
            EndIf

            ! Try to restrict pivoting to the diagonal
            jobA(izb) = 3

            ! Perform symbolic analysis
            Do jdecomp = 1, kdecompmx
              jcn(1:nnz,izb) = cidx
              irn(1:nnz,izb) = ridx
              vals(1:nnz,izb) = tvals(:,izb)
              info = 0
              rinfo = 0.0
              Call MA48AD(msize,msize,nnz,jobA(izb),lia,vals(:,izb),irn(:,izb),jcn(:,izb), &
                & keep(:,izb),cntl,icntl,iwA,info,rinfo)
              If ( info(1) == 0 .and. info(4) <= lia .and. info(2) <= 10 ) Then
                jobB(izb) = 1 ! If analysis is successful, proceed to factorization
                Exit
              ElseIf ( jdecomp == kdecompmx ) Then
                Write(lun_diag,"(1x,a,i5,a,i2,a,i2)") 'kstep=',kstep,', Error during MA48AD. jdecomp=',jdecomp,', info(1)=',info(1)
                Call xnet_terminate('Error during MA48AD',info(1))
              Else
                Write(lun_diag,"(1x,a,i5,a,i2,a,i2)") 'kstep=',kstep,', Warning during MA48AD. jdecomp=',jdecomp,', info(1)=',info(1)

                ! Address any error codes
                If ( info(1) >= 4 .and. jobA(izb) == 3 ) Then
                  Write(lun_diag,"(1x,a,i3,a,i7)") 'Not possible to choose all pivots from diagonal'
                  jobA(izb) = 1
                ElseIf ( info(1) == -3 .or. info(4) > lia ) Then
                  Write(lun_diag,"(1x,a,i7,a,i7)") 'Reallocating MA48 arrays: lia=',lia,' < info(4)=',max(info(3),info(4))
                  lia0 = lia
                  lia = int(1.2*max(info(3),info(4)))

                  Allocate (jcn0(lia,zb_lo:zb_hi))
                  jcn0(1:lia0,:) = jcn
                  Call move_alloc(jcn0,jcn)
                  jcn(:,izb) = 0

                  Allocate (irn0(lia,zb_lo:zb_hi))
                  irn0(1:lia0,:) = irn
                  Call move_alloc(irn0,irn)
                  irn(:,izb) = 0

                  Allocate (vals0(lia,zb_lo:zb_hi))
                  vals0(1:lia0,:) = vals
                  Call move_alloc(vals0,vals)
                  vals(:,izb) = 0.0
                ElseIf(info(2) > 10) Then
                  Write(lun_diag,"(1x,a,i3,a,i7)") 'Reallocating MA48 arrays: info(2)=',info(2),', info(4)=',max(info(3),info(4))
                  lia0 = lia
                  lia = int(2.0*max(info(3),info(4)))

                  Allocate (jcn0(lia,zb_lo:zb_hi))
                  jcn0(1:lia0,:) = jcn
                  Call move_alloc(jcn0,jcn)
                  jcn(:,izb) = 0

                  Allocate (irn0(lia,zb_lo:zb_hi))
                  irn0(1:lia0,:) = irn
                  Call move_alloc(irn0,irn)
                  irn(:,izb) = 0

                  Allocate (vals0(lia,zb_lo:zb_hi))
                  vals0(1:lia0,:) = vals
                  Call move_alloc(vals0,vals)
                  vals(:,izb) = 0.0
                EndIf
              EndIf
            EndDo
          EndIf
        EndDo
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
    Use xnet_controls, Only: idiag, iheat, kitmx, kmon, lun_diag, szbatch, zb_lo, zb_hi, lzactive
    Use xnet_timers, Only: xnet_wtime, start_timer, stop_timer, timer_solve, timer_bksub
    Use xnet_types, Only: dp
    Use xnet_util, Only: xnet_terminate
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep
    Real(dp), Intent(in) :: yrhs(ny,zb_lo:zb_hi)
    Real(dp), Intent(in) :: t9rhs(zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Output variables
    Real(dp), Intent(out) :: dy(ny,zb_lo:zb_hi)
    Real(dp), Intent(out) :: dt9(zb_lo:zb_hi)

    ! Local variables
    Real(dp) :: rhs(msize), dx(msize)
    Real(dp) :: relerr(3)
    Integer :: i, izb, izone, kbksub
    Logical :: mask_1zone(zb_lo:zb_hi)
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

    mask_1zone = .false.
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then

        ! Perform back substitution
        If ( kmon(2,izb) > kitmx ) Then
          jobC(izb) = 2 ! Previous NR iteration failed, so estimate error for possible recalculation of data structures
        Else
          jobC(izb) = 1 ! Do not estimate error
        EndIf

        ! Solve linear system using existing factorization
        Do kbksub = 1, kbksubmx

          rhs(1:ny) = yrhs(:,izb)
          If ( iheat > 0 ) rhs(ny+1) = t9rhs(izb)

          relerr = 0.0
          info = 0
          Call MA48CD(msize,msize,trans,jobC(izb),lia,vals(:,izb),irn(:,izb),keep(:,izb), &
            & cntl,icntl,rhs,dx,relerr,wC,iwC,info)
          If ( info(1) == 0 .and. ( jobC(izb) == 1 .or. maxval(relerr) <= maxerr .or. kbksub == kbksubmx ) ) Then
            dy(:,izb) = dx(1:ny)
            If ( iheat > 0 ) dt9(izb) = dx(ny+1)
            Exit
          ElseIf ( kbksub == kbksubmx .and. info(1) /= 0 ) Then
            Write(lun_diag,"(1x,a,i5,a,i2,a,i2)") 'kstep=',kstep,', Error during MA48CD. kbksub=',kbksub,', info(1)=',info(1)
            Call xnet_terminate('Error during MA48CD',info(1))
          Else
            Write(lun_diag,"(1x,a,i5,a,i2,a,i2)") 'kstep=',kstep,', Warning during MA48CD. kbksub=',kbksub,', info(1)=',info(1)

            ! If the relative error becomes sufficiently large, redo analysis and factorization
            If ( jobC(izb) > 1 .and. maxval(relerr) > maxerr ) Then
              Write(lun_diag,"(1x,a,3es12.5,a)") 'Warning: relerr=',(relerr(i),i=1,3),' > maxerr'
              jobB(izb) = 1
              mask_1zone(izb) = .true.
              Call jacobian_decomp(0,mask_in = mask_1zone)
              mask_1zone(izb) = .false.
            EndIf
          EndIf
        EndDo
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
