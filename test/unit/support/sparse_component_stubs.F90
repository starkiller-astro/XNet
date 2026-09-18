Module nuclear_data
  Use xnet_types, Only: dp
  Implicit None

  Integer :: ny = 0
  Real(dp), Allocatable :: mex(:)
  Character(5), Allocatable :: nname(:)
End Module nuclear_data

Module reaction_data
  Use xnet_types, Only: dp
  Implicit None

  Integer, Allocatable :: la(:,:), le(:,:), mu1(:), mu2(:), mu3(:), mu4(:), nan(:)
  Integer, Allocatable :: n10(:), n11(:), n20(:), n21(:), n22(:)
  Integer, Allocatable :: n30(:), n31(:), n32(:), n33(:)
  Integer, Allocatable :: n40(:), n41(:), n42(:), n43(:), n44(:)
  Real(dp), Allocatable :: a1(:), a2(:), a3(:), a4(:)
  Real(dp), Allocatable :: b1(:,:), b2(:,:), b3(:,:), b4(:,:)
  Real(dp), Allocatable :: dcsect1dt9(:,:), dcsect2dt9(:,:)
  Real(dp), Allocatable :: dcsect3dt9(:,:), dcsect4dt9(:,:)
End Module reaction_data

Module xnet_abundances
  Use xnet_types, Only: dp
  Implicit None

  Real(dp), Allocatable :: yt(:,:)
End Module xnet_abundances

Module xnet_conditions
  Use xnet_types, Only: dp
  Implicit None

  Real(dp), Allocatable :: cv(:)
End Module xnet_conditions

Module xnet_controls
  Use, Intrinsic :: iso_fortran_env, Only: error_unit, output_unit
  Implicit None

  Integer :: idiag = 0
  Integer :: iheat = 0
  Integer :: kitmx = 5
  Integer :: lun_diag = error_unit
  Integer :: lun_stdout = output_unit
  Integer :: nzbatch = 0
  Integer :: nzbatchmx = 0
  Integer :: nzevolve = 0
  Integer :: szbatch = 1
  Integer :: tid = 1
  Integer :: zb_hi = 1
  Integer :: zb_lo = 1
  Integer, Allocatable :: kmon(:,:), ktot(:,:)
  Logical, Allocatable, Target :: lzactive(:)
End Module xnet_controls

Module xnet_parallel
  Use xnet_types, Only: dp
  Implicit None

  Interface parallel_bcast
    Module Procedure parallel_bcast_i0
    Module Procedure parallel_bcast_i1
    Module Procedure parallel_bcast_r0
    Module Procedure parallel_bcast_r1
  End Interface parallel_bcast

Contains

  Logical Function parallel_IOProcessor()
    Implicit None

    parallel_IOProcessor = .True.

    Return
  End Function parallel_IOProcessor

  Subroutine parallel_bcast_i0(value)
    Implicit None

    Integer, Intent(inout) :: value

    Return
  End Subroutine parallel_bcast_i0

  Subroutine parallel_bcast_i1(value)
    Implicit None

    Integer, Intent(inout) :: value(:)

    Return
  End Subroutine parallel_bcast_i1

  Subroutine parallel_bcast_r0(value)
    Implicit None

    Real(dp), Intent(inout) :: value

    Return
  End Subroutine parallel_bcast_r0

  Subroutine parallel_bcast_r1(value)
    Implicit None

    Real(dp), Intent(inout) :: value(:)

    Return
  End Subroutine parallel_bcast_r1

End Module xnet_parallel

Module xnet_timers
  Use xnet_types, Only: dp
  Implicit None

  Real(dp) :: start_timer = 0.0_dp
  Real(dp) :: stop_timer = 0.0_dp
  Real(dp) :: timer_bksub = 0.0_dp
  Real(dp) :: timer_decmp = 0.0_dp
  Real(dp) :: timer_jacob = 0.0_dp
  Real(dp) :: timer_solve = 0.0_dp

Contains

  Real(dp) Function xnet_wtime()
    Implicit None

    xnet_wtime = 0.0_dp

    Return
  End Function xnet_wtime

End Module xnet_timers

Module xnet_util
  Use, Intrinsic :: iso_fortran_env, Only: error_unit
  Implicit None

Contains

  Subroutine xnet_terminate(c_diagnostic,i_diagnostic)
    Implicit None

    Character(*), Intent(in) :: c_diagnostic
    Integer, Intent(in), Optional :: i_diagnostic

    Write(error_unit,'(a)') trim(c_diagnostic)
    If ( present(i_diagnostic) ) Write(error_unit,'(i0)') i_diagnostic
    Stop 1
  End Subroutine xnet_terminate

End Module xnet_util
