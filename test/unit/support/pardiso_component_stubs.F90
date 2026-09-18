Module solver_probe
  Use xnet_types, Only: dp
  Implicit None

  Integer, Parameter :: max_calls = 32
  Integer :: call_count = 0
  Integer :: init_abi = 0
  Integer :: init_calls = 0
  Integer :: mnum_history(max_calls) = 0
  Integer :: phase_history(max_calls) = 0
  Logical :: index_base_valid = .True.
  Logical :: options_valid = .True.
  Real(dp) :: saved_solution(4) = 0.0_dp

Contains

  Subroutine reset_solver_probe
    Implicit None

    call_count = 0
    mnum_history = 0
    phase_history = 0
    index_base_valid = .True.
    options_valid = .True.
    saved_solution = 0.0_dp

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

  Logical Function mutation_is(expected)
    Implicit None

    Character(*), Intent(in) :: expected
    Character(80) :: requested
    Integer :: length, status

    requested = ''
    Call get_environment_variable('XNET_SPARSE_MUTATION',requested,length=length,status=status)
    mutation_is = status == 0 .and. requested(1:length) == expected

    Return
  End Function mutation_is

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

Subroutine pardisoinit(pt,mtype,solver,iparm,dparm,error)
  Use solver_probe, Only: failure_is, init_abi, init_calls
  Use xnet_types, Only: dp, i8
  Implicit None

  Integer(i8), Intent(inout) :: pt(64)
  Integer, Intent(in) :: mtype, solver
  Integer, Intent(inout) :: iparm(64)
  Real(dp), Intent(inout) :: dparm(64)
  Integer, Intent(out) :: error

  init_calls = init_calls + 1
  init_abi = 61
  pt = 0_i8
  iparm = 0
  dparm = 0.0_dp
  iparm(1) = 1
  iparm(7) = 7
  dparm(1) = 0.1_dp
  dparm(8) = 0.8_dp
  If ( mtype /= 11 .or. solver /= 0 ) Then
    error = -100
  ElseIf ( failure_is('pardiso_init') ) Then
    error = -101
  Else
    error = 0
  EndIf

  Return
End Subroutine pardisoinit

Subroutine pardiso(pt,maxfct,mnum,mtype,phase,n,a,ia,ja,perm,nrhs,iparm,msglvl,b,x,error,dparm)
  Use solver_probe, Only: call_count, failure_is, index_base_valid, max_calls, mnum_history, &
    & mutation_is, options_valid, phase_history, saved_solution, solve_crs
  Use xnet_types, Only: dp, i8
  Implicit None

  Integer(i8), Intent(inout) :: pt(64)
  Integer, Intent(in) :: maxfct, mnum, mtype, phase, n, ia(*), ja(*), perm(*), nrhs, msglvl
  Integer, Intent(inout) :: iparm(*)
  Real(dp), Intent(in) :: a(*)
  Real(dp), Intent(inout) :: b(n,nrhs), dparm(64)
  Real(dp), Intent(out) :: x(n,nrhs)
  Integer, Intent(out) :: error

  call_count = call_count + 1
  If ( call_count <= max_calls ) Then
    phase_history(call_count) = phase
    If ( mutation_is('incorrect_phase') ) phase_history(call_count) = phase + 1
    mnum_history(call_count) = mnum
  EndIf
  index_base_valid = index_base_valid .and. ia(1) == 1 .and. ia(n+1) > ia(1)
  options_valid = options_valid .and. maxfct == 2 .and. mnum >= 1 .and. mnum <= 2 .and. &
    & mtype == 11 .and. nrhs == 1 .and. iparm(3) == 1 .and. &
    & iparm(5) == 0 .and. iparm(6) == 0 .and. iparm(12) == 0 .and. &
    & iparm(31) == 0 .and. iparm(35) == 0 .and. msglvl == 1 .and. &
    & all(perm(1:n) == 0) .and. (abs(dparm(8)-8.5_dp) <= epsilon(1.0_dp) .or. &
    & abs(dparm(8)-0.8_dp) <= epsilon(1.0_dp))

  If ( .not. index_base_valid ) Then
    error = -201
  ElseIf ( failure_is('pardiso_factor') .and. (phase == 12 .or. phase == 22) ) Then
    error = -4
  ElseIf ( failure_is('pardiso_solve') .and. phase == 33 ) Then
    error = -7
  ElseIf ( phase == 33 ) Then
    Call solve_crs(n,a,ia,ja,b(:,1),x(:,1),error)
    If ( error == 0 .and. mutation_is('result_offset') ) x(1:n,1) = cshift(x(1:n,1),1)
    If ( error == 0 .and. mutation_is('excessive_residual') ) x(1,1) = x(1,1) + 1.0e-11_dp
    If ( error == 0 .and. mutation_is('cross_zone_copy') ) Then
      If ( mnum == 1 ) Then
        saved_solution(1:n) = x(1:n,1)
      Else
        x(1:n,1) = saved_solution(1:n)
      EndIf
    EndIf
  Else
    x = 0.0_dp
    error = 0
  EndIf

  Return
End Subroutine pardiso
