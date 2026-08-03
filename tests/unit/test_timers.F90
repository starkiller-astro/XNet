Module test_timers
  Use testdrive, Only: check, error_type, new_unittest, unittest_type
  Use xnet_timers, Only: reset_timers, timer_output
  Use xnet_types, Only: dp
  Implicit None
  Private

  Public :: collect_timers

Contains

  Subroutine collect_timers(tests)
    Type(unittest_type), Allocatable, Intent(out) :: tests(:)

    tests = [new_unittest("reset output timer", test_reset_output_timer)]
  End Subroutine collect_timers

  Subroutine test_reset_output_timer(error)
    Type(error_type), Allocatable, Intent(out) :: error

    timer_output = 7.5_dp
    Call reset_timers()

    Call check(error, timer_output == 0.0_dp, "reset_timers must reset timer_output to exact zero")
  End Subroutine test_reset_output_timer

End Module test_timers
