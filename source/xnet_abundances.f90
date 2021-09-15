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
  Real(dp), Allocatable :: xext(:)   ! Abundances of auxiliary/external species
  Real(dp), Allocatable :: aext(:)   ! Mass number of auxiliary/external species
  Real(dp), Allocatable :: zext(:)   ! Charge number of auxiliary/external species

Contains

  Subroutine y_moment(y,ye,ytot,abar,zbar,z2bar,zibar,izb)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates moments of the abundance distribution for the EOS.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny, aa, zz, zz2, zzi
    Use xnet_controls, Only: idiag, lun_diag
    Use xnet_constants, Only: thbim1
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: y(ny)
    Integer, Intent(in),Optional :: izb 

    ! Output variables
    Real(dp), Intent(out) :: ye, ytot, abar, zbar, z2bar, zibar

    ! Local variables
    Real(dp) :: atot, ztot, yext, xext_loc,zext_loc, aext_loc
       
    If (present(izb)) Then
        yext = xext(izb)/aext(izb)
        xext_loc=xext(izb)
        aext_loc=aext(izb)
        zext_loc=zext(izb)
    Else
        yext=0.0
        xext_loc=0.0
        aext_loc=1.0
        zext_loc=0.0
    Endif
    ! Calculate abundance moments
    ytot  = sum(y) + yext
    atot  = sum(y(:) * aa(:)) + xext_loc
    ztot  = sum(y * zz) + yext*zext_loc
    abar  = atot / ytot
    zbar  = ztot / ytot
    z2bar = ( sum(y * zz2) + yext * zext_loc * zext_loc ) / ytot 
    zibar = ( sum(y * zzi) + yext * zext_loc**thbim1 ) / ytot
    ye = ztot / atot
    If ( idiag >= 3 ) Write(lun_diag,"(a4,7es23.15)") 'YMom',ytot,abar,zbar,z2bar,zibar,ye, atot

    Return
  End Subroutine y_moment

End Module xnet_abundances
