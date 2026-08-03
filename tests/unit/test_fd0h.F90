Module test_fd0h
  Use testdrive, Only: check, error_type, new_unittest, unittest_type
  Use fd, Only: fd0h
  Use xnet_types, Only: dp
  Implicit None
  Private

  Public :: collect_fd0h

Contains

  Subroutine collect_fd0h(tests)
    Type(unittest_type), Allocatable, Intent(out) :: tests(:)

    tests = [ &
      new_unittest("exact order-zero identity", test_exact_order_zero_identity) &
    ]
  End Subroutine collect_fd0h

  Subroutine test_exact_order_zero_identity(error)
    Type(error_type), Allocatable, Intent(out) :: error
    Integer :: i
    Real(dp) :: actual, expected
    Real(dp) :: eta(5)

    eta = [-4.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 4.0_dp]

    Do i = 1, Size(eta)
      actual = fd0h(eta(i))
      expected = Log(1.0_dp + Exp(eta(i)))
      Call check(error, &
        Abs(actual - expected) <= 32.0_dp * Epsilon(1.0_dp) * Abs(expected), &
        "fd0h must agree with the exact order-zero identity")
      If ( Allocated(error) ) Return
    EndDo
  End Subroutine test_exact_order_zero_identity

End Module test_fd0h
