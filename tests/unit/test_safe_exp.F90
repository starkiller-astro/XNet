Module test_safe_exp
  Use, Intrinsic :: ieee_arithmetic, Only: ieee_is_finite
  Use testdrive, Only: check, error_type, new_unittest, unittest_type
  Use xnet_types, Only: dp
  Use xnet_util, Only: exp_max, exp_min, safe_exp
  Implicit None
  Private

  Public :: collect_safe_exp

Contains

  Subroutine collect_safe_exp(tests)
    Type(unittest_type), Allocatable, Intent(out) :: tests(:)

    tests = [ &
      new_unittest("ordinary scalar", test_ordinary_scalar), &
      new_unittest("finite clamps", test_finite_clamps), &
      new_unittest("scalar vector consistency", test_scalar_vector_consistency) &
    ]
  End Subroutine collect_safe_exp

  Subroutine test_ordinary_scalar(error)
    Type(error_type), Allocatable, Intent(out) :: error
    Real(dp) :: actual, bound, expected, x

    x = 1.25_dp
    actual = safe_exp(x)
    expected = Exp(x)
    bound = 8.0_dp * Epsilon(1.0_dp) * Max(1.0_dp, Abs(expected))

    Call check(error, Abs(actual - expected) <= bound, &
      "safe_exp must agree with intrinsic exp in the ordinary range")
  End Subroutine test_ordinary_scalar

  Subroutine test_finite_clamps(error)
    Type(error_type), Allocatable, Intent(out) :: error
    Real(dp) :: actual_high, actual_low, bound, expected

    actual_high = safe_exp(Huge(1.0_dp))
    actual_low = safe_exp(-Huge(1.0_dp))

    Call check(error, ieee_is_finite(actual_high), &
      "positive extreme must remain finite")
    If ( Allocated(error) ) Return
    Call check(error, ieee_is_finite(actual_low), &
      "negative extreme must remain finite")
    If ( Allocated(error) ) Return

    expected = Exp(exp_max)
    bound = 16.0_dp * Epsilon(1.0_dp) * Abs(expected)
    Call check(error, Abs(actual_high - expected) <= bound, &
      "positive extreme must clamp at exp_max")
    If ( Allocated(error) ) Return

    expected = Exp(exp_min)
    bound = 16.0_dp * Epsilon(1.0_dp) * Abs(expected)
    Call check(error, Abs(actual_low - expected) <= bound, &
      "negative extreme must clamp at exp_min")
  End Subroutine test_finite_clamps

  Subroutine test_scalar_vector_consistency(error)
    Type(error_type), Allocatable, Intent(out) :: error
    Integer :: i
    Real(dp) :: bound, scalar_value
    Real(dp) :: inputs(5), vector_values(5)

    inputs = [ &
      -Huge(1.0_dp), -2.0_dp, 0.0_dp, 3.0_dp, Huge(1.0_dp) &
    ]
    vector_values = safe_exp(inputs)

    Do i = 1, Size(inputs)
      scalar_value = safe_exp(inputs(i))
      bound = 16.0_dp * Epsilon(1.0_dp) * &
        Max(Abs(scalar_value), Tiny(1.0_dp))
      Call check(error, Abs(vector_values(i) - scalar_value) <= bound, &
        "scalar and vector safe_exp results must agree")
      If ( Allocated(error) ) Return
    EndDo
  End Subroutine test_scalar_vector_consistency

End Module test_safe_exp
