Module solver_probe
  Use xnet_types, Only: dp
  Implicit None

  Integer, Parameter :: max_calls = 32
  Integer :: analysis_calls = 0
  Integer :: analysis_job_history(max_calls) = 0
  Integer :: factor_calls = 0
  Integer :: factor_job_history(max_calls) = 0
  Integer :: init_abi = 48
  Integer :: init_calls = 0
  Integer :: solve_calls = 0
  Integer :: solve_job_history(max_calls) = 0
  Logical :: coordinate_valid = .True.
  Logical :: options_valid = .True.
  Real(dp) :: factor_matrix(4,4,2) = 0.0_dp

Contains

  Subroutine reset_solver_probe
    Implicit None

    analysis_calls = 0
    analysis_job_history = 0
    factor_calls = 0
    factor_job_history = 0
    solve_calls = 0
    solve_job_history = 0
    coordinate_valid = .True.
    options_valid = .True.
    factor_matrix = 0.0_dp

    Return
  End Subroutine reset_solver_probe

  Logical Function failure_is(expected)
    Implicit None

    Character(*), Intent(in) :: expected
    Character(80) :: requested
    Integer :: length, status

    requested = ''
    Call get_environment_variable('XNET_SPARSE_FAIL',requested,length=length,status=status)
    failure_is = status == 0 .and. requested(1:length) == expected

    Return
  End Function failure_is

  Logical Function recovery_is(expected)
    Implicit None

    Character(*), Intent(in) :: expected
    Character(80) :: requested
    Integer :: length, status

    requested = ''
    Call get_environment_variable('XNET_SPARSE_RECOVERY',requested,length=length,status=status)
    recovery_is = status == 0 .and. requested(1:length) == expected

    Return
  End Function recovery_is

  Subroutine solve_dense(n,matrix,b,x,error)
    Implicit None

    Integer, Intent(in) :: n
    Real(dp), Intent(in) :: matrix(n,n), b(n)
    Real(dp), Intent(out) :: x(n)
    Integer, Intent(out) :: error

    Integer :: info, ipiv(n)
    Real(dp) :: work(n,n)

    work = matrix
    x = b
    Call dgesv(n,1,work,n,ipiv,x,n,info)
    If ( info /= 0 ) Then
      error = -5
      x = 0.0_dp
    Else
      error = 0
    EndIf

    Return
  End Subroutine solve_dense

End Module solver_probe

Subroutine MA48ID(cntl,icntl)
  Use solver_probe, Only: init_calls
  Use xnet_types, Only: dp
  Implicit None

  Real(dp), Intent(out) :: cntl(10)
  Integer, Intent(out) :: icntl(20)

  init_calls = init_calls + 1
  cntl = 0.0_dp
  cntl(1) = 0.5_dp
  cntl(2) = 0.1_dp
  icntl = 0
  icntl(1) = 6
  icntl(2) = 6
  icntl(3) = 2
  icntl(5) = 32
  icntl(6) = 2
  icntl(8) = 0

  Return
End Subroutine MA48ID

Subroutine MA48AD(m,n,ne,job,la,a,irn,jcn,keep,cntl,icntl,iw,info,rinfo)
  Use solver_probe, Only: analysis_calls, analysis_job_history, coordinate_valid, &
    & failure_is, max_calls, options_valid, recovery_is
  Use xnet_types, Only: dp
  Implicit None

  Integer, Intent(in) :: m, n, ne, job, la
  Real(dp), Intent(inout) :: a(la)
  Integer, Intent(inout) :: irn(la), jcn(la)
  Integer, Intent(inout) :: keep(*), iw(*)
  Real(dp), Intent(in) :: cntl(10)
  Integer, Intent(in) :: icntl(20)
  Integer, Intent(out) :: info(20)
  Real(dp), Intent(out) :: rinfo(10)

  analysis_calls = analysis_calls + 1
  If ( analysis_calls <= max_calls ) analysis_job_history(analysis_calls) = job
  coordinate_valid = coordinate_valid .and. all(irn(1:ne) >= 1) .and. &
    & all(irn(1:ne) <= m) .and. all(jcn(1:ne) >= 1) .and. all(jcn(1:ne) <= n)
  options_valid = options_valid .and. m == n .and. ne > 0 .and. la >= ne .and. &
    & (job == 1 .or. job == 3) .and. icntl(5) == n .and. &
    & (abs(cntl(2)-0.25_dp) <= epsilon(1.0_dp) .or. abs(cntl(2)-0.1_dp) <= epsilon(1.0_dp))
  info = 0
  rinfo = 0.0_dp
  info(4) = ne
  If ( recovery_is('ma48_analysis_warning') .and. analysis_calls == 1 ) Then
    info(1) = 4
  ElseIf ( recovery_is('ma48_storage_resize') .and. analysis_calls == 1 ) Then
    info(1) = -3
    info(3) = la + 4
    info(4) = la + 4
  ElseIf ( failure_is('ma48_analysis') ) Then
    info(1) = 4
  ElseIf ( failure_is('ma48_storage') ) Then
    info(1) = -3
    info(3) = la + 4
    info(4) = la + 4
  EndIf

  Return
End Subroutine MA48AD

Subroutine MA48BD(m,n,ne,job,la,a,irn,jcn,keep,cntl,icntl,w,iw,info,rinfo)
  Use solver_probe, Only: factor_calls, factor_job_history, factor_matrix, failure_is, &
    & max_calls, options_valid
  Use xnet_types, Only: dp
  Implicit None

  Integer, Intent(in) :: m, n, ne, job, la
  Real(dp), Intent(inout) :: a(la)
  Integer, Intent(inout) :: irn(la), jcn(la), keep(*), iw(*)
  Real(dp), Intent(in) :: cntl(10)
  Integer, Intent(in) :: icntl(20)
  Real(dp), Intent(inout) :: w(*)
  Integer, Intent(out) :: info(20)
  Real(dp), Intent(out) :: rinfo(10)

  Integer :: entry, zone

  factor_calls = factor_calls + 1
  If ( factor_calls <= max_calls ) factor_job_history(factor_calls) = job
  options_valid = options_valid .and. m == n .and. ne > 0 .and. la >= ne .and. &
    & (job == 1 .or. job == 2 .or. job == 3)
  info = 0
  rinfo = 0.0_dp
  info(4) = ne
  If ( failure_is('ma48_factor') ) Then
    info(1) = -7
    Return
  ElseIf ( failure_is('ma48_singular') ) Then
    info(1) = -5
    Return
  EndIf

  zone = mod(factor_calls-1,2) + 1
  factor_matrix(:,:,zone) = 0.0_dp
  Do entry = 1, ne
    factor_matrix(irn(entry),jcn(entry),zone) = a(entry)
  EndDo

  Return
End Subroutine MA48BD

Subroutine MA48CD(m,n,trans,job,la,a,irn,keep,cntl,icntl,rhs,x,relerr,w,iw,info)
  Use solver_probe, Only: factor_matrix, failure_is, max_calls, solve_calls, &
    & solve_dense, solve_job_history
  Use xnet_types, Only: dp
  Implicit None

  Integer, Intent(in) :: m, n, job, la
  Logical, Intent(in) :: trans
  Real(dp), Intent(in) :: a(la), cntl(10), rhs(n)
  Integer, Intent(in) :: irn(la), keep(*), icntl(20)
  Real(dp), Intent(out) :: x(n), relerr(3)
  Real(dp), Intent(inout) :: w(*)
  Integer, Intent(inout) :: iw(*)
  Integer, Intent(out) :: info(20)

  Integer :: error, zone

  solve_calls = solve_calls + 1
  If ( solve_calls <= max_calls ) solve_job_history(solve_calls) = job
  info = 0
  relerr = 0.0_dp
  If ( failure_is('ma48_solve') ) Then
    info(1) = -8
    x = 0.0_dp
    Return
  EndIf
  zone = mod(solve_calls-1,2) + 1
  Call solve_dense(n,factor_matrix(1:n,1:n,zone),rhs,x,error)
  info(1) = error
  If ( trans ) info(1) = -99

  Return
End Subroutine MA48CD
