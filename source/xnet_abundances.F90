!***************************************************************************************************
! xnet_abundances.f90 10/18/17
! This file contains modules and subroutines associated with matter composition.
!***************************************************************************************************

#include "xnet_macros.fh"

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

  Interface y_moment
    Module Procedure y_moment_scalar
    Module Procedure y_moment_vector
  End Interface y_moment

  Private :: y_moment_internal

Contains

  Subroutine y_moment_internal(y,ye,ytot,abar,zbar,z2bar,zibar,xext,yext,zext)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates moments of the abundance distribution for the EOS.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny, aa, zz, zz2, zzi
    Use xnet_constants, Only: thbim1
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: y(ny)
    Real(dp), Intent(in) :: xext, yext, zext

    ! Output variables
    Real(dp), Intent(out) :: ye, ytot, abar, zbar, z2bar, zibar

    ! Local variables
    Real(dp) :: atot, ztot

    ! Calculate abundance moments
    ytot  = sum(y) + yext
    atot  = sum(y * aa) + xext
    ztot  = sum(y * zz) + yext * zext
    abar  = atot / ytot
    zbar  = ztot / ytot
    z2bar = ( sum(y * zz2) + yext * zext * zext ) / ytot 
    zibar = ( sum(y * zzi) + yext * zext**thbim1 ) / ytot
    ye    = ztot / atot

    Return
  End Subroutine y_moment_internal

  Subroutine y_moment_scalar(y,ye,ytot,abar,zbar,z2bar,zibar,xext,aext,zext)
    !-----------------------------------------------------------------------------------------------
    ! This is the interface for the scalar version.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use xnet_controls, Only: idiag, lun_diag
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: y(ny)
    Real(dp), Intent(in) :: xext, aext, zext

    ! Output variables
    Real(dp), Intent(out) :: ye, ytot, abar, zbar, z2bar, zibar

    ! Local variables
    Real(dp) :: yext
       
    ! Calculate abundance moments
    yext = xext / aext
    Call y_moment_internal(y,ye,ytot,abar,zbar,z2bar,zibar,xext,yext,zext)
    If ( idiag >= 3 ) Write(lun_diag,"(a4,6es23.15)") 'YMom',ytot,abar,zbar,z2bar,zibar,ye

    Return
  End Subroutine y_moment_scalar

  Subroutine y_moment_vector(y,ye,ytot,abar,zbar,z2bar,zibar,xext,aext,zext,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This is the interface for the vector version.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny, aa, zz, zz2, zzi
    Use xnet_constants, Only: thbim1
    Use xnet_controls, Only: idiag, lun_diag, zb_lo, zb_hi, lzactive, tid
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: y(ny,zb_lo:zb_hi)
    Real(dp), Intent(in) :: xext(zb_lo:zb_hi), aext(zb_lo:zb_hi), zext(zb_lo:zb_hi)

    ! Output variables
    Real(dp), Intent(out) :: ye(zb_lo:zb_hi), ytot(zb_lo:zb_hi), abar(zb_lo:zb_hi)
    Real(dp), Intent(out) :: zbar(zb_lo:zb_hi), z2bar(zb_lo:zb_hi), zibar(zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer :: izb, k
    Real(dp) :: yext, ntot, atot, ztot, z2tot, zitot
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    !__dir_enter_data &
    !__dir_async(tid) &
    !__dir_create(ye,ytot,abar,zbar,z2bar,zibar) &
    !__dir_copyin(mask,y,xext,aext,zext)
       
    ! Calculate abundance moments
    !__dir_loop_outer(1) &
    !__dir_async(tid) &
    !__dir_present(ye,ytot,abar,zbar,z2bar,zibar) &
    !__dir_present(mask,y,xext,aext,zext) &
    !__dir_private(yext,ntot,atot,ztot,z2tot,zitot)
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        yext = xext(izb) / aext(izb)

        ntot  = 0.0
        atot  = 0.0
        ztot  = 0.0
        z2tot = 0.0
        zitot = 0.0
        !__dir_loop_inner(1) &
        !__dir_reduction(+,ntot,atot,ztot,z2tot,zitot)
        Do k = 1, ny
          ntot  = ntot  + y(k,izb)
          atot  = atot  + y(k,izb) * aa(k)
          ztot  = ztot  + y(k,izb) * zz(k)
          z2tot = z2tot + y(k,izb) * zz2(k)
          zitot = zitot + y(k,izb) * zzi(k)
        EndDo
        ntot  = ntot  + yext
        atot  = atot  + xext(izb)
        ztot  = ztot  + yext * zext(izb)
        z2tot = z2tot + yext * zext(izb) * zext(izb)
        zitot = zitot + yext * zext(izb)**thbim1

        ye(izb)    = ztot  / atot
        ytot(izb)  = ntot
        abar(izb)  = atot  / ntot
        zbar(izb)  = ztot  / ntot
        z2bar(izb) = z2tot / ntot
        zibar(izb) = zitot / ntot
      EndIf
    EndDo
       
    If ( idiag >= 3 ) Then
      !__dir_update &
      !__dir_wait(tid) &
      !__dir_host(ye,ytot,abar,zbar,z2bar,zibar)
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          Write(lun_diag,"(a4,6es23.15)") 'YMom', &
            ytot(izb),abar(izb),zbar(izb),z2bar(izb),zibar(izb),ye(izb)
        EndIf
      EndDo
    EndIf

    !__dir_exit_data &
    !__dir_async(tid) &
    !__dir_copyout(ye,ytot,abar,zbar,z2bar,zibar) &
    !__dir_delete(mask,y,xext,aext,zext)

    Return
  End Subroutine y_moment_vector

End Module xnet_abundances
