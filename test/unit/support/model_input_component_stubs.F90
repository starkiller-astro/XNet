Module xnet_controls
  Use, Intrinsic :: iso_fortran_env, Only: error_unit
  Implicit None

  Integer, Allocatable :: iaux(:)
  Integer :: idiag = -1
  Integer :: lun_ab = 0
  Integer :: lun_diag = error_unit
  Integer :: lun_stderr = error_unit
  Integer :: lun_th = 0
  Integer :: nzone = 1
  Integer :: nzevolve = 1
  Integer :: szbatch = 1
  Integer :: zb_hi = 1
  Integer :: zb_lo = 1
  Real :: t9nse = 10.0
  Logical, Allocatable, Target :: lzactive(:)

Contains

  Subroutine initialize_controls_fixture
    Implicit None

    nzone = 1
    nzevolve = 1
    szbatch = 1
    zb_lo = 1
    zb_hi = 1
    If ( allocated(iaux) ) Deallocate(iaux)
    If ( allocated(lzactive) ) Deallocate(lzactive)
    Allocate(iaux(nzevolve),lzactive(nzevolve))
    iaux = 1
    lzactive = .True.

    Return
  End Subroutine initialize_controls_fixture
End Module xnet_controls

Module xnet_nnu
  Use xnet_types, Only: dp
  Implicit None

  Integer :: nnuspec = 0
  Real(dp) :: fluxcms(1,1,1) = 0.0_dp
  Real(dp) :: tmevnu(1,1,1) = 0.0_dp
End Module xnet_nnu

Module nuclear_data
  Use xnet_types, Only: dp
  Implicit None

  Integer :: ny = 2
  Real(dp), Allocatable :: aa(:)
  Real(dp), Allocatable :: zz(:)

Contains

  Subroutine initialize_nuclear_fixture
    Implicit None

    ny = 2
    If ( allocated(aa) ) Deallocate(aa)
    If ( allocated(zz) ) Deallocate(zz)
    Allocate(aa(ny),zz(ny))
    aa = (/ 1.0_dp, 4.0_dp /)
    zz = (/ 1.0_dp, 2.0_dp /)

    Return
  End Subroutine initialize_nuclear_fixture

  Subroutine index_from_name(name,index)
    Implicit None

    Character(*), Intent(in) :: name
    Integer, Intent(out) :: index

    Select Case (trim(name))
    Case ('h1')
      index = 1
    Case ('he4')
      index = 2
    Case Default
      index = 0
    End Select

    Return
  End Subroutine index_from_name
End Module nuclear_data

Module xnet_conditions
  Use xnet_types, Only: dp
  Implicit None

  Integer :: nh(1) = 1
  Integer :: nhmx = 1
  Integer :: nstart(1) = 1
  Real(dp) :: rhostart(1) = 1.0_dp
  Real(dp) :: t9start(1) = 1.0_dp
  Real(dp) :: tstart(1) = 0.0_dp
  Real(dp) :: tstop(1) = 0.0_dp
  Real(dp) :: tdelstart(1) = 0.0_dp
  Real(dp) :: th(1,1) = 0.0_dp
  Real(dp) :: t9h(1,1) = 0.0_dp
  Real(dp) :: rhoh(1,1) = 0.0_dp
  Real(dp) :: yeh(1,1) = 0.5_dp
  Real(dp) :: yestart(1) = 0.0_dp
End Module xnet_conditions

Module xnet_abundances
  Use xnet_types, Only: dp
  Implicit None

  Real(dp) :: aext(1) = 1.0_dp
  Real(dp) :: xext(1) = 0.0_dp
  Real(dp) :: ystart(2,1) = 0.0_dp
  Real(dp) :: zext(1) = 0.0_dp

  Interface y_moment
    Module Procedure y_moment_scalar
  End Interface y_moment

Contains

  Subroutine y_moment_scalar(y,ye,ytot,abar,zbar,z2bar,zibar,xext_loc,aext_loc,zext_loc)
    Use nuclear_data, Only: aa, zz
    Implicit None

    Real(dp), Intent(in) :: aext_loc, xext_loc, y(:), zext_loc
    Real(dp), Intent(out) :: abar, ye, ytot, z2bar, zbar, zibar

    ytot = sum(y)
    ye = sum(zz*y) + xext_loc*zext_loc/aext_loc
    abar = 1.0_dp / ytot
    zbar = ye / ytot
    z2bar = sum(zz*zz*y) / ytot
    zibar = 1.0_dp

    Return
  End Subroutine y_moment_scalar
End Module xnet_abundances

Module xnet_nse
  Use xnet_types, Only: dp
  Implicit None

  Real(dp) :: ynse(2) = 0.0_dp

Contains

  Subroutine nse_solve(rho,t9,ye)
    Implicit None

    Real(dp), Intent(in) :: rho, t9, ye

    ynse = 0.0_dp

    Return
  End Subroutine nse_solve
End Module xnet_nse
