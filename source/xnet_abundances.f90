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

Contains

  Subroutine y_moment(y,ye,ytot,abar,zbar,z2bar,zibar)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates moments of the abundance distribution for the EOS.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny, aa, zz, zz2, zzi
    Use xnet_controls, Only: idiag, lun_diag
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: y(ny)

    ! Output variables
    Real(dp), Intent(out) :: ye, ytot, abar, zbar, z2bar, zibar

    ! Local variables
    Real(dp) :: atot, ztot

    ! Calculate abundance moments
    ytot  = sum(y)
    atot  = sum(y * aa)
    ztot  = sum(y * zz)
    abar  = atot / ytot
    zbar  = ztot / ytot
    z2bar = sum(y * zz2) / ytot
    zibar = sum(y * zzi) / ytot
    ye = ztot / atot
    If ( idiag >= 3 ) Write(lun_diag,"(a4,6es23.15)") 'YMom',ytot,abar,zbar,z2bar,zibar,ye

    Return
  End Subroutine y_moment

  Subroutine y_moment2(y,ye,ytot,abar,zbar,z2bar,zibar,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates moments of the abundance distribution for the EOS.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny, aa, zz, zz2, zzi
    Use xnet_controls, Only: idiag, lun_diag, zb_lo, zb_hi, lzactive, tid
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: y(ny,zb_lo:zb_hi)

    ! Output variables
    Real(dp), Intent(out) :: ye(zb_lo:zb_hi), ytot(zb_lo:zb_hi), abar(zb_lo:zb_hi)
    Real(dp), Intent(out) :: zbar(zb_lo:zb_hi), z2bar(zb_lo:zb_hi), zibar(zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer :: izb, k
    Real(dp) :: ntot, atot, ztot, z2tot, zitot
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    !$acc enter data async(tid) &
    !$acc copyin(mask)

    ! Calculate abundance moments
    !$acc kernels loop gang async(tid) &
    !$acc present(mask,y,aa,zz,zz2,zzi,ytot,ye,abar,zbar,z2bar,zibar) &
    !$acc private(ntot,atot,ztot,z2tot,zitot)
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        ntot  = 0.0
        atot  = 0.0
        ztot  = 0.0
        z2tot = 0.0
        zitot = 0.0
        !$acc loop vector &
        !$acc reduction(+:ntot,atot,ztot,z2tot,zitot)
        Do k = 1, ny
          ntot  = ntot  + y(k,izb)
          atot  = atot  + y(k,izb) * aa(k)
          ztot  = ztot  + y(k,izb) * zz(k)
          z2tot = z2tot + y(k,izb) * zz2(k)
          zitot = zitot + y(k,izb) * zzi(k)
        EndDo
        ytot(izb)  = ntot
        ye(izb)    = ztot  / atot
        abar(izb)  = atot  / ntot
        zbar(izb)  = ztot  / ntot
        z2bar(izb) = z2tot / ntot
        zibar(izb) = zitot / ntot
      EndIf
    EndDo
    If ( idiag >= 3 ) Then
      !$acc update wait(tid) &
      !$acc host(ytot,abar,zbar,z2bar,zibar,ye)
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          Write(lun_diag,"(a4,6es23.15)") 'YMom',ytot,abar,zbar,z2bar,zibar,ye
        EndIf
      EndDo
    EndIf

    !$acc exit data async(tid) &
    !$acc delete(mask)

    Return
  End Subroutine y_moment2

End Module xnet_abundances
