!***************************************************************************************************
! xnet_timers.f90 10/18/17
! This file contains modules and subroutines for internal XNet timers
!***************************************************************************************************

Module xnet_timers
  !-------------------------------------------------------------------------------------------------
  ! This module contains the performance timers for XNet.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None
  Real(dp) :: timer_burner = 0.0  ! Total burner execution time
  Real(dp) :: timer_xnet   = 0.0  ! Total XNet execution time
  Real(dp) :: timer_setup  = 0.0  ! Data loading or preprocessing time
  Real(dp) :: timer_csect  = 0.0  ! Cross section calculation time
  Real(dp) :: timer_deriv  = 0.0  ! Derivative calculation time
  Real(dp) :: timer_jacob  = 0.0  ! Jacobian building time
  Real(dp) :: timer_decmp  = 0.0  ! LU Decomposition time
  Real(dp) :: timer_bksub  = 0.0  ! Backsubstitution time
  Real(dp) :: timer_nraph  = 0.0  ! Newton Raphson iteration timer
  Real(dp) :: timer_tstep  = 0.0  ! Time integration step timer
  Real(dp) :: timer_solve  = 0.0  ! Solution time
  Real(dp) :: timer_scrn   = 0.0  ! Screening and EOS time
  Real(dp) :: timer_eos    = 0.0  ! Screening and EOS time
  Real(dp) :: timer_nse    = 0.0  ! NSE timer
  Real(dp) :: timer_nseinit= 0.0  ! NSE init timer
  Real(dp) :: timer_nsesolv= 0.0  ! NSE solver timer
  Real(dp) :: timer_nsenrap= 0.0  ! NSE Newton Raphson timer
  Real(dp) :: timer_nsels  = 0.0  ! NSE line search timer
  Real(dp) :: timer_nseeval= 0.0  ! NSE function evaluation timer
  Real(dp) :: timer_nsescrn= 0.0  ! NSE screening timer
  Real(dp) :: timer_output = 0.0  ! Output time
  Real(dp) :: start_timer  = 0.0  ! cpu time at the beginning of the timer block
  Real(dp) :: stop_timer   = 0.0  ! cpu time at the end of the timer block
  !$omp threadprivate(timer_burner,timer_xnet,timer_setup,timer_csect,timer_deriv,timer_jacob,timer_decmp, &
  !$omp   timer_bksub,timer_nraph,timer_tstep,timer_solve,timer_scrn,timer_eos,timer_output,start_timer,stop_timer, &
  !$omp   timer_nse,timer_nseinit,timer_nsesolv,timer_nsenrap,timer_nsels,timer_nseeval,timer_nsescrn)

Contains

  Function xnet_wtime()
    !-----------------------------------------------------------------------------------------------
    ! This function returns the wall time in a manner akin to omp_get_wtime().
    !-----------------------------------------------------------------------------------------------
    Use xnet_types, Only: i8
    Implicit None

    ! Function variable
    Real(dp) :: xnet_wtime

    ! Local variables
    Integer(i8) :: clock_read
    Integer(i8) :: clock_rate
    Integer(i8) :: clock_max

    Call system_clock(clock_read,clock_rate,clock_max)
    xnet_wtime = real(clock_read,dp) / real(clock_rate,dp)

    Return
  End Function xnet_wtime

  Subroutine reset_timers
    !-----------------------------------------------------------------------------------------------
    ! This routine resets timers for zone-independent timing.
    !-----------------------------------------------------------------------------------------------
    Implicit None
    timer_xnet  = 0.0
    timer_setup = 0.0
    timer_csect = 0.0
    timer_deriv = 0.0
    timer_jacob = 0.0
    timer_decmp = 0.0
    timer_bksub = 0.0
    timer_nraph = 0.0
    timer_tstep = 0.0
    timer_solve = 0.0
    timer_scrn  = 0.0
    timer_eos   = 0.0
    timer_nse    = 0.0
    timer_nseinit= 0.0
    timer_nsesolv= 0.0
    timer_nsenrap= 0.0
    timer_nsels  = 0.0
    timer_nseeval= 0.0
    timer_nsescrn= 0.0

    Return
  End Subroutine reset_timers

End Module xnet_timers
