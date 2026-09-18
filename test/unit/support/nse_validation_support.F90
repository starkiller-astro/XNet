Module nse_validation_support
  Use xnet_types, Only: dp
  Implicit None
  Private

  Integer, Parameter, Public :: nse_fail_identity = 1
  Integer, Parameter, Public :: nse_fail_duplicate = 2
  Integer, Parameter, Public :: nse_fail_finite = 4
  Integer, Parameter, Public :: nse_fail_nonnegative = 8
  Integer, Parameter, Public :: nse_fail_mass = 16
  Integer, Parameter, Public :: nse_fail_charge = 32
  Integer, Parameter, Public :: nse_fail_ye = 64
  Integer, Parameter, Public :: nse_fail_l1 = 128
  Integer, Parameter, Public :: nse_fail_linf = 256
  Integer, Parameter, Public :: nse_fail_dominant = 512

  Type, Public :: nse_validation_tolerances
    Real(dp) :: mass
    Real(dp) :: charge
    Real(dp) :: ye
    Real(dp) :: l1
    Real(dp) :: linf
  End Type nse_validation_tolerances

  Type, Public :: nse_validation_metrics
    Real(dp) :: mass_residual = 0.0_dp
    Real(dp) :: charge_residual = 0.0_dp
    Real(dp) :: ye_error = 0.0_dp
    Real(dp) :: l1 = 0.0_dp
    Real(dp) :: linf = 0.0_dp
    Real(dp) :: dominant_error = 0.0_dp
    Integer :: dominant_index = 0
  End Type nse_validation_metrics

  Public :: apply_nse_validation_gates, evaluate_nse_candidate

Contains

  Subroutine apply_nse_validation_gates(metrics,tolerances,failures)
    Implicit None

    Type(nse_validation_metrics), Intent(in) :: metrics
    Type(nse_validation_tolerances), Intent(in) :: tolerances
    Integer, Intent(inout) :: failures

    If ( abs(metrics%mass_residual) > tolerances%mass ) &
      & failures = ior(failures,nse_fail_mass)
    If ( abs(metrics%charge_residual) > tolerances%charge ) &
      & failures = ior(failures,nse_fail_charge)
    If ( abs(metrics%ye_error) > tolerances%ye ) &
      & failures = ior(failures,nse_fail_ye)
    If ( metrics%l1 > tolerances%l1 ) failures = ior(failures,nse_fail_l1)
    If ( metrics%linf > tolerances%linf ) failures = ior(failures,nse_fail_linf)
    If ( metrics%dominant_error > tolerances%linf ) &
      & failures = ior(failures,nse_fail_dominant)

    Return
  End Subroutine apply_nse_validation_gates

  Subroutine evaluate_nse_candidate(expected_names,candidate_names,expected_a,expected_z, &
    & expected_n,candidate_a,candidate_z,candidate_n,expected,candidate,ye,tolerances, &
    & metrics,failures)
    Implicit None

    Character(5), Intent(in) :: candidate_names(:), expected_names(:)
    Integer, Intent(in) :: expected_a(:), expected_n(:), expected_z(:)
    Real(dp), Intent(in) :: candidate(:), candidate_a(:), candidate_n(:), candidate_z(:)
    Real(dp), Intent(in) :: expected(:), ye
    Type(nse_validation_tolerances), Intent(in) :: tolerances
    Type(nse_validation_metrics), Intent(out) :: metrics
    Integer, Intent(out) :: failures

    Integer :: i, j, species_count
    Real(dp), Allocatable :: difference(:)

    failures = 0
    metrics = nse_validation_metrics()
    species_count = size(expected)
    If ( size(expected_names) /= species_count .OR. &
      & size(candidate_names) /= species_count .OR. &
      & size(candidate) /= species_count .OR. &
      & size(expected_a) /= species_count .OR. &
      & size(expected_z) /= species_count .OR. &
      & size(expected_n) /= species_count .OR. &
      & size(candidate_a) /= species_count .OR. &
      & size(candidate_z) /= species_count .OR. &
      & size(candidate_n) /= species_count ) Then
      failures = ior(failures,nse_fail_identity)
      Return
    EndIf

    If ( any(expected_names /= candidate_names) ) &
      & failures = ior(failures,nse_fail_identity)
    If ( any(real(expected_a,dp) /= candidate_a) .OR. &
      & any(real(expected_z,dp) /= candidate_z) .OR. &
      & any(real(expected_n,dp) /= candidate_n) ) &
      & failures = ior(failures,nse_fail_identity)
    Do i = 1, species_count
      Do j = i + 1, species_count
        If ( expected_names(i) == expected_names(j) .OR. &
          & candidate_names(i) == candidate_names(j) ) Then
          failures = ior(failures,nse_fail_duplicate)
          failures = ior(failures,nse_fail_identity)
        EndIf
      EndDo
    EndDo

    If ( .NOT.all_binary64_finite(expected) .OR. &
      & .NOT.all_binary64_finite(candidate) ) Then
      failures = ior(failures,nse_fail_finite)
      Return
    EndIf
    If ( any(expected < 0.0_dp) .OR. any(candidate < 0.0_dp) ) &
      & failures = ior(failures,nse_fail_nonnegative)

    Allocate(difference(species_count))
    difference = abs(candidate-expected)
    metrics%mass_residual = sum(candidate) - 1.0_dp
    metrics%charge_residual = sum((candidate_z/candidate_a-ye)*candidate)
    metrics%ye_error = sum(candidate_z*candidate/candidate_a) - ye
    metrics%l1 = sum(difference)
    metrics%linf = maxval(difference)
    metrics%dominant_index = maxloc(expected,dim=1)
    metrics%dominant_error = difference(metrics%dominant_index)
    Deallocate(difference)

    Call apply_nse_validation_gates(metrics,tolerances,failures)

    Return
  End Subroutine evaluate_nse_candidate

  Logical Function all_binary64_finite(values) Result(finite)
    Use, Intrinsic :: iso_fortran_env, Only: int64
    Implicit None

    Real(dp), Intent(in) :: values(:)

    Integer(int64), Parameter :: exponent_mask = int(z'7ff0000000000000',int64)
    Integer :: i
    Integer(int64) :: bits

    finite = .True.
    Do i = 1, size(values)
      bits = transfer(values(i),bits)
      If ( iand(bits,exponent_mask) == exponent_mask ) Then
        finite = .False.
        Exit
      EndIf
    EndDo

    Return
  End Function all_binary64_finite

End Module nse_validation_support
