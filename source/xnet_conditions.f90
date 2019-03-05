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
  !$omp threadprivate(t,tt,to,tdel,tdel_next,tdel_old,t9t,t9,t9o,t9dot,rhot,rho,rhoo,yet,ye,yeo)

  Real(dp), Allocatable :: cv(:)                 ! Specific heat at constant volume
  Real(dp), Allocatable :: etae(:)               ! Ratio of electron chemical potential to kT
  Real(dp), Allocatable :: detaedt9(:)           ! d(etae)/dT9
  !$omp threadprivate(cv,etae,detaedt9)

  Integer, Allocatable :: nt(:)    ! Point in thermo trajectory for current time
  Integer, Allocatable :: ntt(:)   ! Point in thermo trajectory for trial time
  Integer, Allocatable :: nto(:)   ! Point in thermo trajectory for previous time
  Integer, Allocatable :: ints(:)  ! Nucleus governing timestep
  Integer, Allocatable :: intso(:) ! Nucleus governing last timestep
  !$omp threadprivate(nto,nt,ntt,ints,intso)

  ! Arrays for thermodynamic trajectory that the network follows
  Integer, Parameter    :: nhmx = 50000 ! The max number of thermo points
  Integer, Allocatable  :: nh(:)        ! The actual number of points
  Integer, Allocatable  :: nstart(:)
  Real(dp), Allocatable :: tstart(:),tstop(:),tdelstart(:),t9start(:),rhostart(:),yestart(:)
  Real(dp), Allocatable :: th(:,:),t9h(:,:),rhoh(:,:),yeh(:,:)
  !$omp threadprivate(nh,nstart,tstart,tstop,tdelstart,t9start,rhostart,yestart,th,t9h,rhoh,yeh)

Contains

  Subroutine read_thermo_file( thermo_file, thermo_desc, ierr, mask_in )
    !-----------------------------------------------------------------------------------------------
    ! Read the thermdynamic trajectory
    !-----------------------------------------------------------------------------------------------
    Use, Intrinsic :: iso_fortran_env, Only: iostat_end
!   Use neutrino_data, Only: fluxcms, tmevnu                                                    !NNU
    Use xnet_controls, Only: idiag, lun_diag, lun_th, nzone, szbatch, nzbatchmx, lzactive
    Use xnet_util, Only: replace_tabs, readnext, xnet_terminate
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: thermo_file(nzone)

    ! Output variables
    Character(80), Intent(out) :: thermo_desc(nzbatchmx)
    Integer, Intent(out) :: ierr

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(:)

    ! Local variables
    Integer, Parameter :: max_line_length = 1024
    Integer :: pos, i, n, izb, izone
    Character(max_line_length) :: line
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask => mask_in
    Else
      mask => lzactive
    EndIf
    If ( .not. any(mask(:)) ) Return

    ! Initialize
    tstart = 0.0
    tstop = 0.0
    tdelstart = 0.0
    th = 0.0
    t9h = 0.0
    rhoh = 0.0
    yeh = 0.0
!   fluxcms = 0.0                                                                               !NNU
!   tmevnu = 0.0                                                                                !NNU
    ierr = 0

    !$omp critical(th_read)
    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
        izone = izb + szbatch - 1

        Open(newunit=lun_th, file=trim(thermo_file(izone)), action='read', status='old', iostat=ierr)
        If ( ierr /= 0 ) Then
          Call xnet_terminate('Failed to open input file: '//trim(thermo_file(izone)))
        EndIf

        Read(lun_th,*) thermo_desc(izb)
        Read(lun_th,*) tstart(izb)
        Read(lun_th,*) tstop(izb)
        Read(lun_th,*) tdelstart(izb)
        Do n = 1, nhmx
          line(:) = ' '
          Read(lun_th,"(a)",iostat=ierr) line
          If ( ierr == iostat_end ) Then
            If ( idiag >= 1 ) Write(lun_diag,"(a,i6,a)") 'End of Thermo File Reached after ',n,' records'
            Exit
          ElseIf ( ierr /= 0 ) Then
            Call xnet_terminate('Failed while trying to read input file: '//trim(thermo_file(izone)),ierr)
          Else
            ! Parse the line as space delimited, one value at a time since the format could vary
            Call replace_tabs(line)

            ! Read the required data
            pos = 1
            Call readnext(line,pos,th(n,izb))
            Call readnext(line,pos,t9h(n,izb))
            Call readnext(line,pos,rhoh(n,izb))
            If ( pos == 0 ) Call xnet_terminate('Not enough columns in thermo file: '//trim(thermo_file(izone)))

            ! See if electron fraction is in file, otherwise continue to next line
            Call readnext(line,pos,yeh(n,izb))
            If ( pos == 0 ) Cycle

            ! See if neutrino data is in file, otherwise continue to next line
!           Do i = 1, 4                                                                         !NNU
!             Call readnext(line,pos,fluxcms(n,i,izb))                                          !NNU
!           EndDo                                                                               !NNU
!           If ( pos == 0 ) Cycle                                                               !NNU
!           Do i = 1, 4                                                                         !NNU
!             Call readnext(line,pos,tmevnu(n,i,izb))                                           !NNU
!           EndDo                                                                               !NNU
          EndIf
        EndDo
        nh(izb) = n - 1
        Close(lun_th)

        ! Do not use tdelstart from thermo files
        tdelstart(izb) = min(0.0,tdelstart(izb))

        ! Log thermo description
        If ( idiag >= 0 ) Write(lun_diag,"(a)") thermo_desc(izb)
      EndIf
    EndDo
    !$omp end critical(th_read)

    ! Convert to appropriate units (CGS, except temperature (GK) and neutrino flux)
!   t9h = t9h * 1.0e-9
!   fluxcms = 1.0e-42 * fluxcms                                                                 !NNU

    Return
  End Subroutine read_thermo_file

  Subroutine t9rhofind(kstep,tf,nf,t9f,rhof,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates t9 and rho as a function of time, either via interpolation or from an
    ! analytic expression.
    !-----------------------------------------------------------------------------------------------
    Use xnet_controls, Only: lun_diag, lun_stdout, nzbatchmx, lzactive
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep
    Real(dp), Intent(in) :: tf(:)

    ! Input/Output variables
    Integer, Intent(inout) :: nf(size(tf))
    Real(dp), Intent(inout) :: t9f(size(tf)), rhof(size(tf))

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(:)

    ! Local variables
    Real(dp) :: dt, rdt, dt9, drho
    Integer :: n, izb
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask => mask_in
    Else
      mask => lzactive
    End If
    If ( .not. any(mask(:)) ) Return

    Do izb = 1, nzbatchmx
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
  End Subroutine t9rhofind

End Module xnet_conditions
