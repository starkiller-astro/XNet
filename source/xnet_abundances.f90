!***************************************************************************************************
! xnet_abundances.f90 10/18/17
! This file contains modules and subroutines associated with matter composition.
!***************************************************************************************************

Module xnet_abundances
  !-------------------------------------------------------------------------------------------------
  ! This module contains the abundances of the nuclear species, at the previous time (yo), current
  ! time (y), and trial time (yt), as well as the time derivatives (ydot), at the trial time. The
  ! size of these arrays is allocated in an external routine, where the initial values are set.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None
  Real(dp), Allocatable :: ystart(:,:) ! Abundances at start time
  Real(dp), Allocatable :: yo(:,:)     ! Abundances at previous time
  Real(dp), Allocatable :: y(:,:)      ! Abundances at current time
  Real(dp), Allocatable :: yt(:,:)     ! Abundances at trial time
  Real(dp), Allocatable :: ydot(:,:)   ! Abundance time derivatives at trial time
  !$omp threadprivate(yo,y,yt,ydot,ystart)

Contains

  Subroutine y_moment(y,ye,ytot,abar,zbar,z2bar,zibar)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates moments of the abundance distribution for the EOS.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: aa, zz, zz2, zzi
    Use xnet_controls, Only: idiag, lun_diag
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: y(:)

    ! Output variables
    Real(dp), Intent(out) :: ye, ytot, abar, zbar, z2bar, zibar

    ! Local variables
    Real(dp) :: atot, ztot

    ! Calculate abundance moments
    ytot  = sum(y(:))
    atot  = sum(y(:) * aa(:))
    ztot  = sum(y(:) * zz(:))
    abar  = atot / ytot
    zbar  = ztot / ytot
    z2bar = sum(y(:) * zz2(:)) / ytot
    zibar = sum(y(:) * zzi(:)) / ytot
    ye = ztot / atot
    If ( idiag >= 3 ) Write(lun_diag,"(a4,6es23.15)") 'YMom',ytot,abar,zbar,z2bar,zibar,ye

    Return
  End Subroutine y_moment

End Module xnet_abundances
