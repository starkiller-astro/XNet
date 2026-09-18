Module xnet_controls
  Use, Intrinsic :: iso_fortran_env, Only: error_unit
  Implicit None

  Integer :: idiag = 0
  Integer :: iheat = 1
  Integer :: iscrn = 1
  Integer :: lun_diag = error_unit
  Integer :: nzevolve = 0
  Integer :: tid = 1
  Integer :: zb_lo = 1
  Integer :: zb_hi = 1
  Logical, Allocatable, Target :: lzactive(:)
End Module xnet_controls

Module nuclear_data
  Use xnet_types, Only: dp
  Implicit None

  Integer :: ny = 0
  Real(dp), Allocatable :: aa(:)
  Real(dp), Allocatable :: zz(:)
  Real(dp), Allocatable :: zz2(:)
  Real(dp), Allocatable :: zz53(:)
  Real(dp), Allocatable :: zzi(:)
End Module nuclear_data
