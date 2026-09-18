Module xnet_controls
  Use, Intrinsic :: iso_fortran_env, Only: error_unit, output_unit
  Implicit None

  Integer :: idiag = 0
  Integer :: ineutrino = 0
  Integer :: iscrn = 0
  Integer :: itsout = 0
  Integer :: lun_diag = error_unit
  Integer :: lun_stderr = error_unit
  Integer :: lun_stdout = output_unit
  Integer :: nzevolve = 0
  Integer :: szbatch = 1
  Integer :: tid = 1
  Integer :: zb_lo = 1
  Integer :: zb_hi = 1
  Logical, Allocatable, Target :: lzactive(:)
End Module xnet_controls

Module nuclear_data
  Use xnet_types, Only: dp
  Implicit None

  Integer, Parameter :: ng = 24
  Integer :: ny = 0
  Integer :: izmax = 0
  Real(dp), Allocatable :: aa(:)
  Real(dp), Allocatable :: angm(:)
  Real(dp), Allocatable :: be(:)
  Real(dp), Allocatable :: g(:,:)
  Integer, Allocatable :: ia(:)
  Integer, Allocatable :: iz(:)
  Real(dp), Allocatable :: mex(:)
  Real(dp), Allocatable :: mm(:)
  Character(5), Allocatable :: nname(:)
  Real(dp), Allocatable :: nn(:)
  Real(dp), Allocatable :: t9i(:)
  Real(dp), Allocatable :: zz(:)
  Real(dp), Allocatable :: zz2(:)
  Real(dp), Allocatable :: zzi(:)
  Real(dp), Allocatable :: zseq(:)
  Real(dp), Allocatable :: zseq53(:)
  Real(dp), Allocatable :: zseqi(:)
End Module nuclear_data

Module xnet_eos
  Use nuclear_data, Only: zz
  Use xnet_types, Only: dp
  Implicit None

  Integer :: eos_screen_calls = 0

Contains

  Subroutine eos_interface(t9,rho,y,ye,cv,etae,detaedt9,xext,aext,zext)
    Implicit None

    Real(dp), Intent(in) :: aext, rho, t9, xext, y(:), zext
    Real(dp), Intent(out) :: cv, detaedt9, etae, ye

    cv = 1.0_dp
    etae = 0.0_dp
    detaedt9 = 0.0_dp
    ye = sum(zz*y) + xext*zext/aext

    Return
  End Subroutine eos_interface

  Subroutine eos_screen(t9,rho,y,etae,detaedt9,ztilde,zinter,lambda0,gammae,dztildedt9, &
    & xext,aext,zext)
    Implicit None

    Real(dp), Intent(in) :: aext, detaedt9, etae, rho, t9, xext, y(:), zext
    Real(dp), Intent(out) :: dztildedt9, gammae, lambda0, zinter, ztilde

    eos_screen_calls = eos_screen_calls + 1
    ztilde = 1.0_dp
    zinter = 1.0_dp
    lambda0 = 1.0e-2_dp
    gammae = 1.0e-2_dp
    dztildedt9 = 0.0_dp

    Return
  End Subroutine eos_screen

End Module xnet_eos

Module xnet_parallel
  Implicit None

Contains

  Subroutine parallel_abort(str,errorcode)
    Implicit None

    Character(*), Intent(in), Optional :: str
    Integer, Intent(in), Optional :: errorcode

    If ( present(str) ) Write(*,*) trim(str)
    If ( present(errorcode) ) Write(*,*) errorcode
    Stop 1
  End Subroutine parallel_abort

End Module xnet_parallel
