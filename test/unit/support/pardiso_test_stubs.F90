Subroutine pardisoinit(pt,mtype,solver,iparm,dparm,error)
  Use xnet_types, Only: dp, i8
  Implicit None

  Integer(i8), Intent(inout) :: pt(64)
  Integer, Intent(in) :: mtype, solver
  Integer, Intent(inout) :: iparm(64)
  Real(dp), Intent(inout) :: dparm(64)
  Integer, Intent(out) :: error

  pt = 0_i8
  iparm = 0
  dparm = 0.0_dp
  error = 0

  Return
End Subroutine pardisoinit

Subroutine pardiso(pt,maxfct,mnum,mtype,phase,n,a,ia,ja,perm,nrhs,iparm,msglvl,b,x,error,dparm)
  Use xnet_types, Only: dp, i8
  Implicit None

  Integer(i8), Intent(inout) :: pt(64)
  Integer, Intent(in) :: maxfct, mnum, mtype, phase, n, ia(*), ja(*), perm(*), nrhs, msglvl
  Integer, Intent(inout) :: iparm(*)
  Real(dp), Intent(in) :: a(*)
  Real(dp), Intent(inout) :: b(n,nrhs), dparm(64)
  Real(dp), Intent(out) :: x(n,nrhs)
  Integer, Intent(out) :: error

  x = b
  error = 0

  Return
End Subroutine pardiso
