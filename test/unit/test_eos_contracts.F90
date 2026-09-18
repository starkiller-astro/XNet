Module eos_component_fixture
  Use nuclear_data, Only: aa, ny, zz, zz2, zz53, zzi
  Use xnet_constants, Only: five3rd, thbim1
  Use xnet_controls, Only: lzactive, nzevolve, zb_hi, zb_lo
  Use xnet_eos, Only: eos_initialize
  Use xnet_types, Only: dp
  Implicit None
  Private

  Integer, Parameter, Public :: nstate = 5
  Integer, Parameter, Public :: nzone = 3
  Character(32), Public :: effectiveness_mutation = ''

  Public :: initialize_component
  Public :: set_interface_fixture
  Public :: set_scientific_fixture

Contains

  Subroutine initialize_component
    Implicit None

    Integer :: length, status

    ny = 4
    nzevolve = nzone
    zb_lo = 1
    zb_hi = nzone
    Allocate (aa(ny),zz(ny),zz2(ny),zz53(ny),zzi(ny),lzactive(nzone))
    aa = (/ 1.0_dp, 4.0_dp, 12.0_dp, 56.0_dp /)
    zz = (/ 1.0_dp, 2.0_dp, 6.0_dp, 26.0_dp /)
    zz2 = zz**2
    zz53 = zz**five3rd
    zzi = zz**thbim1
    lzactive = .True.

    Call get_environment_variable('XNET_EOS_MUTATION',effectiveness_mutation, &
      & length=length,status=status)
    If ( status /= 0 ) effectiveness_mutation = ''
    Call eos_initialize

    Return
  End Subroutine initialize_component

  Subroutine set_interface_fixture(t9,rho,y,xext,aext,zext,mask)
    Implicit None

    Logical, Intent(out) :: mask(nzone)
    Real(dp), Intent(out) :: aext(nzone), rho(nzone), t9(nzone), xext(nzone)
    Real(dp), Intent(out) :: y(4,nzone), zext(nzone)

    t9 = (/ 0.08_dp, 0.9_dp, 3.0_dp /)
    rho = (/ 2.0e4_dp, 4.0e7_dp, 1.0e9_dp /)
    y = 0.0_dp
    y(:,1) = (/ 0.60_dp, 0.36_dp/4.0_dp, 0.04_dp/12.0_dp, 0.0_dp /)
    y(2,2) = 0.25_dp
    y(3,3) = 0.40_dp/12.0_dp
    y(4,3) = 0.60_dp/56.0_dp
    xext = 0.0_dp
    aext = 1.0_dp
    zext = 0.0_dp
    mask = (/ .True., .False., .True. /)
    If ( trim(effectiveness_mutation) == 'mask' ) mask(2) = .True.

    Return
  End Subroutine set_interface_fixture

  Subroutine set_scientific_fixture(t9,rho,y)
    Implicit None

    Real(dp), Intent(out) :: rho(nstate), t9(nstate), y(4,nstate)

    t9 = (/ 0.032_dp, 0.17_dp, 0.8_dp, 4.0_dp, 9.0_dp /)
    rho = (/ 2.5e2_dp, 3.0e6_dp, 2.0e9_dp, 7.0e7_dp, 1.0e10_dp /)
    y = 0.0_dp
    y(:,1) = (/ 0.70_dp, 0.28_dp/4.0_dp, 0.02_dp/12.0_dp, 0.0_dp /)
    y(2,2) = 0.25_dp
    y(3,3) = 1.0_dp/12.0_dp
    y(4,4) = 1.0_dp/56.0_dp
    y(3,5) = 0.50_dp/12.0_dp
    y(4,5) = 0.50_dp/56.0_dp

    Return
  End Subroutine set_scientific_fixture

End Module eos_component_fixture

Module test_eos_contracts
  Use, Intrinsic :: ieee_arithmetic, Only: ieee_is_finite
  Use eos_component_fixture, Only: effectiveness_mutation, nstate, nzone, &
    & set_interface_fixture, set_scientific_fixture
  Use testdrive, Only: check, error_type, new_unittest, unittest_type
  Use xnet_constants, Only: amu
  Use xnet_eos, Only: eos_interface, eos_screen
  Use xnet_types, Only: dp
  Implicit None
  Private

  Real(dp), Parameter :: comparison_absolute_tolerance = 1.0e-12_dp
  Real(dp), Parameter :: comparison_relative_tolerance = 1.0e-11_dp
  Real(dp), Parameter :: sentinel = -9.87654321e99_dp

  Public :: collect_eos_contracts

Contains

  Subroutine collect_eos_contracts(testsuite)
    Implicit None

    Type(unittest_type), Allocatable, Intent(out) :: testsuite(:)

#if defined(TEST_STARKILLER)
    testsuite = [ &
      & new_unittest('EOS scalar/vector contract',test_eos_interface), &
      & new_unittest('screening scalar/vector contract',test_screen_interface), &
      & new_unittest('Helmholtz versus direct Timmes states',test_scientific_states) ]
#else
    testsuite = [ &
      & new_unittest('EOS scalar/vector contract',test_eos_interface), &
      & new_unittest('screening scalar/vector contract',test_screen_interface) ]
#endif

    Return
  End Subroutine collect_eos_contracts

  Subroutine test_eos_interface(error)
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error

    Integer :: zone
    Logical :: mask(nzone)
    Real(dp) :: aext(nzone), cv(nzone), detaedt9(nzone), etae(nzone), rho(nzone)
    Real(dp) :: repeat_cv(nzone), repeat_detaedt9(nzone), repeat_etae(nzone)
    Real(dp) :: repeat_ye(nzone), scalar_cv, scalar_detaedt9, scalar_etae, scalar_ye
    Real(dp) :: t9(nzone), xext(nzone), y(4,nzone), ye(nzone), zext(nzone)

    Call set_interface_fixture(t9,rho,y,xext,aext,zext,mask)
    ye = sentinel
    cv = sentinel
    etae = sentinel
    detaedt9 = sentinel
    Call eos_interface(t9,rho,y,ye,cv,etae,detaedt9,xext,aext,zext,mask)

    Do zone = 1, nzone, 2
      If ( zone == 1 .and. trim(effectiveness_mutation) == 'argument_swap' ) Then
        Call eos_interface(rho(zone),t9(zone),y(:,zone),scalar_ye,scalar_cv, &
          & scalar_etae,scalar_detaedt9,xext(zone),aext(zone),zext(zone))
      Else
        Call eos_interface(t9(zone),rho(zone),y(:,zone),scalar_ye,scalar_cv, &
          & scalar_etae,scalar_detaedt9,xext(zone),aext(zone),zext(zone))
      EndIf
      Call check_close(error,ye(zone),scalar_ye)
      If ( allocated(error) ) Return
      Call check_close(error,cv(zone),scalar_cv)
      If ( allocated(error) ) Return
      Call check_close(error,etae(zone),scalar_etae)
      If ( allocated(error) ) Return
      Call check_close(error,detaedt9(zone),scalar_detaedt9)
      If ( allocated(error) ) Return
      Call check(error,ye(zone) > 0.0_dp .and. ye(zone) <= 1.0_dp)
      If ( allocated(error) ) Return
      Call check(error,cv(zone) > 0.0_dp .and. ieee_is_finite(cv(zone)))
      If ( allocated(error) ) Return
      Call check(error,ieee_is_finite(etae(zone)) .and. ieee_is_finite(detaedt9(zone)))
      If ( allocated(error) ) Return
    EndDo
    Call check_inactive(error,(/ ye(2), cv(2), etae(2), detaedt9(2) /))
    If ( allocated(error) ) Return

    repeat_ye = sentinel
    repeat_cv = sentinel
    repeat_etae = sentinel
    repeat_detaedt9 = sentinel
    Call eos_interface(t9,rho,y,repeat_ye,repeat_cv,repeat_etae,repeat_detaedt9, &
      & xext,aext,zext,mask)
    Call check_repeatable(error,ye,repeat_ye)
    If ( allocated(error) ) Return
    Call check_repeatable(error,cv,repeat_cv)
    If ( allocated(error) ) Return
    Call check_repeatable(error,etae,repeat_etae)
    If ( allocated(error) ) Return
    Call check_repeatable(error,detaedt9,repeat_detaedt9)

    Return
  End Subroutine test_eos_interface

  Subroutine test_screen_interface(error)
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error

    Integer :: zone
    Logical :: all_active(nzone), mask(nzone)
    Real(dp) :: aext(nzone), cv(nzone), detaedt9(nzone), dztildedt9(nzone)
    Real(dp) :: etae(nzone), gammae(nzone), lambda0(nzone), repeat_dztildedt9(nzone)
    Real(dp) :: repeat_gammae(nzone), repeat_lambda0(nzone), repeat_zinter(nzone)
    Real(dp) :: repeat_ztilde(nzone), rho(nzone), scalar_dztildedt9, scalar_gammae
    Real(dp) :: scalar_lambda0, scalar_zinter, scalar_ztilde, t9(nzone), xext(nzone)
    Real(dp) :: y(4,nzone), ye(nzone), zext(nzone), zinter(nzone), ztilde(nzone)

    Call set_interface_fixture(t9,rho,y,xext,aext,zext,mask)
    all_active = .True.
    ye = 0.0_dp
    cv = 0.0_dp
    etae = 0.0_dp
    detaedt9 = 0.0_dp
    Call eos_interface(t9,rho,y,ye,cv,etae,detaedt9,xext,aext,zext,all_active)

    ztilde = sentinel
    zinter = sentinel
    lambda0 = sentinel
    gammae = sentinel
    dztildedt9 = sentinel
    Call eos_screen(t9,rho,y,etae,detaedt9,ztilde,zinter,lambda0,gammae, &
      & dztildedt9,xext,aext,zext,mask)

    Do zone = 1, nzone, 2
      Call eos_screen(t9(zone),rho(zone),y(:,zone),etae(zone),detaedt9(zone), &
        & scalar_ztilde,scalar_zinter,scalar_lambda0,scalar_gammae, &
        & scalar_dztildedt9,xext(zone),aext(zone),zext(zone))
      Call check_close(error,ztilde(zone),scalar_ztilde)
      If ( allocated(error) ) Return
      Call check_close(error,zinter(zone),scalar_zinter)
      If ( allocated(error) ) Return
      Call check_close(error,lambda0(zone),scalar_lambda0)
      If ( allocated(error) ) Return
      Call check_close(error,gammae(zone),scalar_gammae)
      If ( allocated(error) ) Return
      Call check_close(error,dztildedt9(zone),scalar_dztildedt9)
      If ( allocated(error) ) Return
      Call check(error,ztilde(zone) > 0.0_dp .and. zinter(zone) > 0.0_dp)
      If ( allocated(error) ) Return
      Call check(error,lambda0(zone) > 0.0_dp .and. gammae(zone) > 0.0_dp)
      If ( allocated(error) ) Return
      Call check(error,ieee_is_finite(dztildedt9(zone)))
      If ( allocated(error) ) Return
    EndDo
    Call check_inactive(error, &
      & (/ ztilde(2), zinter(2), lambda0(2), gammae(2), dztildedt9(2) /))
    If ( allocated(error) ) Return

    repeat_ztilde = sentinel
    repeat_zinter = sentinel
    repeat_lambda0 = sentinel
    repeat_gammae = sentinel
    repeat_dztildedt9 = sentinel
    Call eos_screen(t9,rho,y,etae,detaedt9,repeat_ztilde,repeat_zinter, &
      & repeat_lambda0,repeat_gammae,repeat_dztildedt9,xext,aext,zext,mask)
    Call check_repeatable(error,ztilde,repeat_ztilde)
    If ( allocated(error) ) Return
    Call check_repeatable(error,zinter,repeat_zinter)
    If ( allocated(error) ) Return
    Call check_repeatable(error,lambda0,repeat_lambda0)
    If ( allocated(error) ) Return
    Call check_repeatable(error,gammae,repeat_gammae)
    If ( allocated(error) ) Return
    Call check_repeatable(error,dztildedt9,repeat_dztildedt9)

    Return
  End Subroutine test_screen_interface

#if defined(TEST_STARKILLER)
  Subroutine test_scientific_states(error)
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error

    Real(dp), Parameter :: cv_reference_cgs(nstate) = (/ &
      & 2.0653230156968951e8_dp, 4.6199576060239784e7_dp, &
      & 1.8932479108972270e7_dp, 9.8736503214137957e7_dp, &
      & 4.3195927677688472e7_dp /)
    Real(dp), Parameter :: detaedt_reference_per_k(nstate) = (/ &
      & -4.9592562974862513e-8_dp, -1.0876372950003569e-7_dp, &
      & -8.4721540555149227e-8_dp, -1.0485292044799550e-9_dp, &
      & -1.2091574746119109e-9_dp /)
    Real(dp), Parameter :: etae_reference(nstate) = (/ &
      & -1.8807503609806382_dp, 1.8319585249344883e1_dp, &
      & 6.7689241964853380e1_dp, 2.8245988707608238_dp, &
      & 1.0296830235239677e1_dp /)
    Real(dp), Parameter :: expected_ye(nstate) = (/ &
      & 0.85_dp, 0.5_dp, 0.5_dp, 26.0_dp/56.0_dp, 27.0_dp/56.0_dp /)
    Integer :: state
    Real(dp) :: cv, cv_reference, detaedt9, detaedt9_reference, etae, etae_expected
    Real(dp) :: rho(nstate), t9(nstate), y(4,nstate), ye

    Call set_scientific_fixture(t9,rho,y)
    Do state = 1, nstate
      Call eos_interface(t9(state),rho(state),y(:,state),ye,cv,etae,detaedt9, &
        & 0.0_dp,1.0_dp,0.0_dp)
      cv_reference = cv_reference_cgs(state)*amu*1.0e9_dp
      detaedt9_reference = detaedt_reference_per_k(state)*1.0e9_dp
      etae_expected = etae_reference(state)
      If ( trim(effectiveness_mutation) == 'unit_conversion' .and. state == 1 ) &
        & cv_reference = cv_reference_cgs(state)*amu
      If ( trim(effectiveness_mutation) == 'reference' .and. state == 1 ) &
        & etae_expected = 1.01_dp*etae_expected

      Call check_close(error,ye,expected_ye(state),1.0e-13_dp,1.0e-13_dp)
      If ( allocated(error) ) Return
      Call check_close(error,cv,cv_reference,1.0e-12_dp,1.0e-5_dp)
      If ( allocated(error) ) Return
      Call check_close(error,etae,etae_expected,5.0e-6_dp,2.0e-5_dp)
      If ( allocated(error) ) Return
      Call check_close(error,detaedt9,detaedt9_reference,1.0e-6_dp,3.0e-4_dp)
      If ( allocated(error) ) Return
    EndDo

    Return
  End Subroutine test_scientific_states
#endif

  Subroutine check_close(error,actual,expected,absolute_tolerance,relative_tolerance)
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error
    Real(dp), Intent(in) :: actual, expected
    Real(dp), Intent(in), Optional :: absolute_tolerance, relative_tolerance

    Real(dp) :: absolute_limit, relative_limit

    absolute_limit = comparison_absolute_tolerance
    relative_limit = comparison_relative_tolerance
    If ( present(absolute_tolerance) ) absolute_limit = absolute_tolerance
    If ( present(relative_tolerance) ) relative_limit = relative_tolerance
    Call check(error,abs(actual-expected) <= absolute_limit + relative_limit*abs(expected))

    Return
  End Subroutine check_close

  Subroutine check_inactive(error,values)
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error
    Real(dp), Intent(in) :: values(:)

    Call check(error,all(abs(values-sentinel) <= 0.0_dp))

    Return
  End Subroutine check_inactive

  Subroutine check_repeatable(error,first,second)
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error
    Real(dp), Intent(in) :: first(:), second(:)

    Call check(error,all(abs(first-second) <= 0.0_dp))

    Return
  End Subroutine check_repeatable

End Module test_eos_contracts

Program eos_contract_test_runner
  Use, Intrinsic :: iso_fortran_env, Only: error_unit
  Use eos_component_fixture, Only: initialize_component
  Use test_eos_contracts, Only: collect_eos_contracts
  Use testdrive, Only: new_testsuite, run_testsuite, testsuite_type
  Implicit None

  Integer :: stat
  Type(testsuite_type), Allocatable :: testsuites(:)

  Call initialize_component
  stat = 0
#if defined(TEST_STARKILLER)
  testsuites = [ new_testsuite('STARKILLER EOS contracts',collect_eos_contracts) ]
#else
  testsuites = [ new_testsuite('Bahcall EOS contracts',collect_eos_contracts) ]
#endif
  Write(error_unit,'("# Testing: ",a)') testsuites(1)%name
  Call run_testsuite(testsuites(1)%collect,error_unit,stat,parallel=.False.)
  If ( stat > 0 ) Then
    Write(error_unit,'(i0,1x,a)') stat,'test(s) failed'
    Stop 1
  EndIf
End Program eos_contract_test_runner
