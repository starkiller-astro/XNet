!***************************************************************************************************
! xnet_conditions.f90 10/18/17
! This file contains modules and subroutines associated with the thermodynamic conditions in the
! matter undergoing nucleosynthesis.
!***************************************************************************************************

Module xnet_conditions
  !-------------------------------------------------------------------------------------------------
  ! This module contains data on the current time and thermodynamic conditions.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None
  Real(dp), Allocatable :: t(:)                  ! Time at the beginning of the current timestep
  Real(dp), Allocatable :: tt(:)                 ! Trial time for end of current timestep
  Real(dp), Allocatable :: to(:)                 ! Time at the beginning of the previous timestep
  Real(dp), Allocatable :: tdel(:)               ! Trial duration of timestep
  Real(dp), Allocatable :: tdel_next(:)          ! Integrator estimate of duration for next timestep
  Real(dp), Allocatable :: tdel_old(:)           ! Duration of the previous timestep
  Real(dp), Allocatable :: t9t(:),rhot(:),yet(:) ! Temperature, density, and electron fraction at trial time
  Real(dp), Allocatable :: t9(:),rho(:),ye(:)    ! Temperature, density, and electron fraction at current time
  Real(dp), Allocatable :: t9o(:),rhoo(:),yeo(:) ! Temperature, density, and electron fraction at previous time
  Real(dp), Allocatable :: t9dot(:)              ! Time derivative of temperature at trial time

  Real(dp), Allocatable :: cv(:)                 ! Specific heat at constant volume
  Real(dp), Allocatable :: etae(:)               ! Ratio of electron chemical potential to kT
  Real(dp), Allocatable :: detaedt9(:)           ! d(etae)/dT9

  Integer, Allocatable :: nt(:)    ! Point in thermo trajectory for current time
  Integer, Allocatable :: ntt(:)   ! Point in thermo trajectory for trial time
  Integer, Allocatable :: nto(:)   ! Point in thermo trajectory for previous time
  Integer, Allocatable :: ints(:)  ! Nucleus governing timestep
  Integer, Allocatable :: intso(:) ! Nucleus governing last timestep

  ! Arrays for thermodynamic trajectory that the network follows
  Integer, Parameter    :: nhmx = 50000 ! The max number of thermo points
  Integer, Allocatable  :: nh(:)        ! The actual number of points
  Integer, Allocatable  :: nstart(:)
  Real(dp), Allocatable :: tstart(:),tstop(:),tdelstart(:),t9start(:),rhostart(:),yestart(:)
  Real(dp), Allocatable :: th(:,:),t9h(:,:),rhoh(:,:),yeh(:,:)

  Interface t9rhofind
    Module Procedure t9rhofind_scalar
    Module Procedure t9rhofind_vector
  End Interface t9rhofind

Contains

  Subroutine t9rhofind_scalar(kstep,tf,nf,t9f,rhof,ns,ts,t9s,rhos)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates t9 and rho as a function of time, either via interpolation or from an
    ! analytic expression. (scalar interface)
    !-----------------------------------------------------------------------------------------------
    Use, Intrinsic :: iso_fortran_env, Only: lun_stdout=>output_unit
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep, ns
    Real(dp), Intent(in) :: tf, ts(:), t9s(:), rhos(:)

    ! Output variables
    Integer, Intent(out) :: nf
    Real(dp), Intent(out) :: t9f, rhof

    ! Local variables
    Real(dp) :: dt, rdt, dt9, drho
    Integer :: n

    ! For constant conditions (ns = 1), set temperature and density
    If ( ns == 1 ) Then
      t9f = t9s(1)
      rhof = rhos(1)
      nf = 1

    ! Otherwise, calculate T9 and rho by interpolation
    Else
      Do n = 1, ns
        If ( tf <= ts(n) ) Exit
      EndDo
      nf = n
      If ( n > 1 .and. n <= ns ) Then
        rdt = 1.0 / (ts(n)-ts(n-1))
        dt = tf - ts(n-1)
        dt9 = t9s(n) - t9s(n-1)
        drho = rhos(n) - rhos(n-1)
        t9f = dt*rdt*dt9 + t9s(n-1)
        rhof = dt*rdt*drho + rhos(n-1)
      ElseIf ( n == 1 ) Then
        t9f = t9s(1)
        rhof = rhos(1)
      Else
        t9f = t9s(ns)
        rhof = rhos(ns)
        Write(lun_stdout,*) 'Time beyond thermodynamic range',tf,' >',ts(ns)
      EndIf
    EndIf

    Return
  End Subroutine t9rhofind_scalar

  Subroutine t9rhofind_vector(kstep,tf,nf,t9f,rhof,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates t9 and rho as a function of time, either via interpolation or from an
    ! analytic expression. (vector interface)
    !-----------------------------------------------------------------------------------------------
    Use xnet_controls, Only: lun_diag, lun_stdout, zb_lo, zb_hi, lzactive
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep
    Real(dp), Intent(in) :: tf(zb_lo:zb_hi)

    ! Input/Output variables
    Integer, Intent(inout) :: nf(zb_lo:zb_hi)
    Real(dp), Intent(inout) :: t9f(zb_lo:zb_hi), rhof(zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Real(dp) :: dt, rdt, dt9, drho
    Integer :: n, izb
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    End If
    If ( .not. any(mask) ) Return

    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then

        ! For constant conditions (nh = 1), set temperature and density
        If ( nh(izb) == 1 ) Then
          t9f(izb) = t9h(1,izb)
          rhof(izb) = rhoh(1,izb)
          nf(izb) = 1

        ! Otherwise, calculate T9 and rho by interpolation
        Else
          Do n = 1, nh(izb)
            If ( tf(izb) <= th(n,izb) ) Exit
          EndDo
          nf(izb) = n
          If ( n > 1 .and. n <= nh(izb) ) Then
            rdt = 1.0 / (th(n,izb)-th(n-1,izb))
            dt = tf(izb) - th(n-1,izb)
            dt9 = t9h(n,izb) - t9h(n-1,izb)
            drho = rhoh(n,izb) - rhoh(n-1,izb)
            t9f(izb) = dt*rdt*dt9 + t9h(n-1,izb)
            rhof(izb) = dt*rdt*drho + rhoh(n-1,izb)
          ElseIf ( n == 1 ) Then
            t9f(izb) = t9h(1,izb)
            rhof(izb) = rhoh(1,izb)
          Else
            t9f(izb) = t9h(nh(izb),izb)
            rhof(izb) = rhoh(nh(izb),izb)
            Write(lun_stdout,*) 'Time beyond thermodynamic range',tf(izb),' >',th(nh(izb),izb)
          EndIf
        EndIf
      EndIf
    EndDo

    Return
  End Subroutine t9rhofind_vector

End Module xnet_conditions
