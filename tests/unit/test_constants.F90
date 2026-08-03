Module test_constants
  Use testdrive, Only: check, error_type, new_unittest, unittest_type
  Use xnet_constants, Only: clt, ln_2, third
  Use xnet_types, Only: dp
  Implicit None
  Private

  Public :: collect_constants

Contains

  Subroutine collect_constants(tests)
    Type(unittest_type), Allocatable, Intent(out) :: tests(:)

    tests = [new_unittest("explicit real(dp) constants", test_dp_constants)]
  End Subroutine collect_constants

  Subroutine test_dp_constants(error)
    Type(error_type), Allocatable, Intent(out) :: error

    Call check(error, Kind(1.0) /= dp, &
      "strict-kind test must preserve the compiler default real kind")
    If ( Allocated(error) ) Return

    Call check(error, third == 1.0_dp / 3.0_dp, &
      "third must be evaluated using real(dp) division")
    If ( Allocated(error) ) Return

    Call check(error, ln_2 == Log(2.0_dp), &
      "ln_2 must select the real(dp) logarithm")
    If ( Allocated(error) ) Return

    Call check(error, clt == 2.99792458e+10_dp, &
      "clt must preserve its real(dp) literal value")
  End Subroutine test_dp_constants

End Module test_constants

Program strict_kind_constants
  Use, Intrinsic :: iso_fortran_env, Only: error_unit
  Use test_constants, Only: collect_constants
  Use testdrive, Only: run_testsuite
  Implicit None

  Integer :: stat

  stat = 0
  Call run_testsuite(collect_constants, error_unit, stat)
  If ( stat > 0 ) Error Stop 1
End Program strict_kind_constants
