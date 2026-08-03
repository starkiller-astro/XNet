Module failure_fixture_suite
  Use testdrive, Only: check, error_type, new_unittest, unittest_type
  Implicit None
  Private

  Public :: collect_failure_fixture

Contains

  Subroutine collect_failure_fixture(tests)
    Type(unittest_type), Allocatable, Intent(out) :: tests(:)

    tests = [new_unittest("deliberate unexpected failure", deliberate_failure)]
  End Subroutine collect_failure_fixture

  Subroutine deliberate_failure(error)
    Type(error_type), Allocatable, Intent(out) :: error

    Call check(error, .False., &
      "deliberate failure used only to verify process-status propagation")
  End Subroutine deliberate_failure

End Module failure_fixture_suite

Program failure_fixture
  Use, Intrinsic :: iso_fortran_env, Only: error_unit
  Use failure_fixture_suite, Only: collect_failure_fixture
  Use testdrive, Only: run_testsuite
  Implicit None

  Integer :: stat

  stat = 0
  Call run_testsuite(collect_failure_fixture, error_unit, stat)
  ! Fall through with zero only when the runner misses the deliberate failure.
  If ( stat /= 0 ) Error Stop 1
End Program failure_fixture
