Module solver_probe
  Use xnet_types, Only: dp
  Implicit None

  Integer, Parameter :: max_calls = 32
  Integer :: call_count = 0
  Integer :: init_abi = 0
  Integer :: init_calls = 0
  Integer :: handle_history(max_calls) = 0
  Integer :: mnum_history(max_calls) = 0
  Integer :: phase_history(max_calls) = 0
  Logical :: index_base_valid = .True.
  Logical :: options_valid = .True.

Contains

  Subroutine reset_solver_probe
    Implicit None

    call_count = 0
    handle_history = 0
    mnum_history = 0
    phase_history = 0
    index_base_valid = .True.
    options_valid = .True.

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

  Subroutine solve_crs(n,a,ia,ja,b,x,error)
    Implicit None

    Integer, Intent(in) :: n, ia(n+1), ja(*)
    Real(dp), Intent(in) :: a(*), b(n)
    Real(dp), Intent(out) :: x(n)
    Integer, Intent(out) :: error

    Integer :: i, info, ipiv(n), k
    Real(dp) :: matrix(n,n)

    matrix = 0.0_dp
    Do i = 1, n
      Do k = ia(i), ia(i+1)-1
        matrix(i,ja(k)) = a(k)
      EndDo
    EndDo
    x = b
    Call dgesv(n,1,matrix,n,ipiv,x,n,info)
    If ( info /= 0 ) Then
      error = -4
      x = 0.0_dp
    Else
      error = 0
    EndIf

    Return
  End Subroutine solve_crs

End Module solver_probe

Subroutine pardisoinit(pt,mtype,iparm)
  Use solver_probe, Only: init_abi, init_calls
  Use xnet_types, Only: i8
  Implicit None

  Integer(i8), Intent(inout) :: pt(*)
  Integer, Intent(in) :: mtype
  Integer, Intent(inout) :: iparm(*)

  init_calls = init_calls + 1
  init_abi = 64
  pt(1:64) = 0_i8
  iparm(1:64) = 0
  iparm(1) = 1
  iparm(7) = 7
  If ( mtype /= 11 ) iparm(1) = -100

  Return
End Subroutine pardisoinit

Subroutine pardiso(pt,maxfct,mnum,mtype,phase,n,a,ia,ja,perm,nrhs,iparm,msglvl,b,x,error)
  Use solver_probe, Only: call_count, failure_is, handle_history, index_base_valid, max_calls, &
    & mnum_history, options_valid, phase_history, solve_crs
  Use xnet_types, Only: dp, i8
  Implicit None

  Integer(i8), Intent(inout) :: pt(*)
  Integer, Intent(in) :: maxfct, mnum, mtype, phase, n, ia(*), ja(*), nrhs, msglvl
  Integer, Intent(inout) :: iparm(*), perm(*)
  Real(dp), Intent(in) :: a(*)
  Real(dp), Intent(inout) :: b(*)
  Real(dp), Intent(out) :: x(*)
  Integer, Intent(out) :: error

  call_count = call_count + 1
  If ( call_count <= max_calls ) Then
    handle_history(call_count) = pt(1)
    phase_history(call_count) = phase
    mnum_history(call_count) = mnum
  EndIf
  pt(1) = call_count
  index_base_valid = index_base_valid .and. ia(1) == 1 .and. ia(n+1) > ia(1)
  options_valid = options_valid .and. maxfct == 1 .and. mnum == 1 .and. &
    & mtype == 11 .and. nrhs == 1 .and. iparm(3) == 0 .and. &
    & iparm(5) == 0 .and. iparm(6) == 0 .and. iparm(12) == 0 .and. &
    & iparm(31) == 0 .and. iparm(35) == 0 .and. iparm(36) == 0 .and. &
    & (iparm(8) == 0 .or. iparm(8) == 3) .and. &
    & msglvl == 1 .and. all(perm(1:n) == 0)

  If ( .not. index_base_valid ) Then
    error = -201
  ElseIf ( failure_is('pardiso_factor') .and. (phase == 12 .or. phase == 22) ) Then
    error = -4
  ElseIf ( failure_is('pardiso_solve') .and. phase == 33 ) Then
    error = -7
  ElseIf ( phase == 33 ) Then
    Call solve_crs(n,a,ia,ja,b(1:n),x(1:n),error)
  Else
    x(1:n) = 0.0_dp
    error = 0
  EndIf

  Return
End Subroutine pardiso
