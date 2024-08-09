!***************************************************************************************************
! eos_starkiller.f90 10/18/17
! Interface to starkiller
! This file contains routines which calculate EoS quantites needed to calculate screening
! corrections for reaction rates.
!***************************************************************************************************

#include "xnet_macros.fh"

Module xnet_eos
  Use xnet_types, Only: dp
  Implicit None
  Real(dp), Allocatable :: ye(:), ytot(:), abar(:), zbar(:), z2bar(:), zibar(:), sratio(:)

  Interface eos_interface
    Module Procedure eos_interface_scalar
    Module Procedure eos_interface_vector
  End Interface

  Interface eos_screen
    Module Procedure eos_screen_scalar
    Module Procedure eos_screen_vector
  End Interface

  Interface salpeter_ratio
    Module Procedure salpeter_ratio_scalar
    Module Procedure salpeter_ratio_vector
  End Interface

Contains

  Subroutine eos_initialize
    !-----------------------------------------------------------------------------------------------
    ! This routine initializes starkiller
    !-----------------------------------------------------------------------------------------------
    Use xnet_controls, Only: nzevolve

    Use actual_eos_module, Only: actual_eos_init
    Implicit None

    Call actual_eos_init()

    Allocate (ye(nzevolve))
    Allocate (ytot(nzevolve))
    Allocate (abar(nzevolve))
    Allocate (zbar(nzevolve))
    Allocate (z2bar(nzevolve))
    Allocate (zibar(nzevolve))
    Allocate (sratio(nzevolve))
    ye = 0.0
    ytot = 0.0
    abar = 0.0
    zbar = 0.0
    z2bar = 0.0
    zibar = 0.0
    sratio = 0.0

    Return
  End Subroutine eos_initialize

  Subroutine eosx(t9,rho,ye,abar,zbar,cv,etae,detaedt9)
    !-----------------------------------------------------------------------------------------------
    ! This routine interfaces with and calls the underlying EoS.
    !-----------------------------------------------------------------------------------------------
    !__dir_routine_seq
    Use xnet_constants, Only: amu
    Use xnet_controls, Only: iheat, iscrn
    Use xnet_types, Only: dp

    Use actual_eos_module, Only: xnet_actual_eos
    Use eos_type_module, Only: eos_input_rt, eos_t
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9, rho, ye, abar, zbar

    ! Ouput variables
    Real(dp), Intent(out) :: cv, etae, detaedt9

    ! Local variables
    Type(eos_t) :: eos_state

    cv = 0.0
    etae = 0.0
    detaedt9 = 0.0
    If ( iscrn > 0 .or. iheat > 0 ) Then

      ! Load input variables for the eos
      eos_state%rho = rho
      eos_state%T = t9*1e9
      eos_state%y_e = ye
      eos_state%abar = abar
      eos_state%zbar = zbar

      ! Call the eos
      Call xnet_actual_eos(eos_input_rt,eos_state)

      ! Convert units from ergs/g to MeV/nucleon and K to GK
      cv = eos_state%cv * amu * 1e9
      etae = eos_state%eta
      detaedt9 = eos_state%detadt * 1e9
    EndIf

  End Subroutine eosx

  Subroutine eos_interface_scalar(t9,rho,y,ye,cv,etae,detaedt9,xext,aext,zext)
    !-----------------------------------------------------------------------------------------------
    ! This routine updates the equation of state for changes in temperature and density.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use xnet_abundances, Only: y_moment
    Use xnet_controls, Only: idiag, lun_diag
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9, rho, y(ny), xext, aext, zext

    ! Ouput variables
    Real(dp), Intent(out) :: ye, cv, etae, detaedt9

    ! Local variables
    Real(dp) :: ytot, abar, zbar, z2bar, zibar

    ! Calculate Ye
    Call y_moment(y,ye,ytot,abar,zbar,z2bar,zibar,xext,aext,zext)

    ! Call the eos
    Call eosx(t9,rho,ye,abar,zbar,cv,etae,detaedt9)

    If ( idiag >= 3 ) Write(lun_diag,"(a,6es23.15)") 'EOS',t9,rho,ye,cv,etae,detaedt9

    Return
  End Subroutine eos_interface_scalar

  Subroutine eos_interface_vector(t9,rho,y,ye,cv,etae,detaedt9,xext,aext,zext,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine updates the equation of state for changes in temperature and density.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use xnet_abundances, Only: y_moment
    Use xnet_controls, Only: idiag, lun_diag, zb_lo, zb_hi, lzactive
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9(zb_lo:zb_hi), rho(zb_lo:zb_hi), y(ny,zb_lo:zb_hi)
    Real(dp), Intent(in) :: xext(zb_lo:zb_hi), aext(zb_lo:zb_hi), zext(zb_lo:zb_hi)

    ! Ouput variables
    Real(dp), Intent(out) :: ye(zb_lo:zb_hi), cv(zb_lo:zb_hi)
    Real(dp), Intent(out) :: etae(zb_lo:zb_hi), detaedt9(zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer :: izb
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    !__dir_enter_data &
    !__dir_async &
    !__dir_create(ye,cv,etae,detaedt9) &
    !__dir_create(ytot,abar,zbar,z2bar,zibar) &
    !__dir_copyin(mask,t9,rho,y,xext,aext,zext)

    ! Calculate Ye
    Call y_moment(y,ye,ytot(zb_lo:zb_hi), &
      & abar(zb_lo:zb_hi),zbar(zb_lo:zb_hi),z2bar(zb_lo:zb_hi),zibar(zb_lo:zb_hi), &
      & xext(zb_lo:zb_hi),aext(zb_lo:zb_hi),zext(zb_lo:zb_hi),mask_in = mask)

    ! Call the eos
    !__dir_loop_outer(1) &
    !__dir_async &
    !__dir_present(mask,t9,rho,ye,abar,zbar,cv,etae,detaedt9)
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        Call eosx(t9(izb),rho(izb),ye(izb),abar(izb),zbar(izb),cv(izb),etae(izb),detaedt9(izb))
      EndIf
    EndDo

    If ( idiag >= 3 ) Then
      !__dir_update &
      !__dir_wait &
      !__dir_host(t9,rho,ye,cv,etae,detaedt9)
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          Write(lun_diag,"(a,6es23.15)") 'EOS',t9(izb),rho(izb),ye(izb),cv(izb),etae(izb),detaedt9(izb)
        EndIf
      EndDo
    EndIf

    !__dir_exit_data &
    !__dir_async &
    !__dir_copyout(ye,cv,etae,detaedt9) &
    !__dir_copyout(ytot,abar,zbar,z2bar,zibar) &
    !__dir_delete(mask,t9,rho,y,xext,aext,zext)

    !__dir_wait

    Return
  End Subroutine eos_interface_vector

  Subroutine eos_screen_scalar(t9,rho,y,etae,detaedt9,ztilde,zinter,lambda0,gammae,dztildedt9,xext,aext,zext)
    !-----------------------------------------------------------------------------------------------
    ! This routine uses the current composition and prior updates to the Equation of State to
    ! calculate the factors needed for screening.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use xnet_abundances, Only: y_moment
    Use xnet_controls, Only: idiag, iheat, lun_diag
    Use xnet_types, Only: dp
    Use xnet_util, Only: plasma
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9, rho, y(ny), etae, detaedt9, xext, aext, zext

    ! Output variables
    Real(dp), Intent(out) :: ztilde, zinter, lambda0, gammae, dztildedt9

    ! Local variables
    Real(dp) :: ye, ytot, abar, zbar, z2bar, zibar
    Real(dp) :: sratio

    ! Calculate Ye and other needed moments of the abundance distribution
    Call y_moment(y,ye,ytot,abar,zbar,z2bar,zibar,xext,aext,zext)

    ! Calculate ratio f'/f for electrons (Salpeter, Eq. 24; DGC, Eq. 5)
    Call salpeter_ratio(etae,sratio,dztildedt9)
    ztilde = sqrt(z2bar + zbar*sratio) ! DGC, Eq. 4
    dztildedt9 = 0.5*zbar/ztilde * dztildedt9*detaedt9

    ! Calculate plasma quantities
    Call plasma(t9,rho,ytot,ye,zbar,zibar,ztilde,zinter,lambda0,gammae)
    If ( idiag >= 3 ) Write(lun_diag,"(a14,9es23.15)") 'EOS Screen', &
      & t9,rho,ye,z2bar,zbar,sratio,ztilde,ztilde*lambda0,gammae

    Return
  End Subroutine eos_screen_scalar

  Subroutine eos_screen_vector(t9,rho,y,etae,detaedt9,ztilde,zinter,lambda0,gammae,dztildedt9,xext,aext,zext,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine uses the current composition and prior updates to the Equation of State to
    ! calculate the factors needed for screening.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use xnet_abundances, Only: y_moment
    Use xnet_controls, Only: idiag, lun_diag, zb_lo, zb_hi, lzactive
    Use xnet_types, Only: dp
    Use xnet_util, Only: plasma
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9(zb_lo:zb_hi), rho(zb_lo:zb_hi), y(ny,zb_lo:zb_hi)
    Real(dp), Intent(in) :: etae(zb_lo:zb_hi), detaedt9(zb_lo:zb_hi)
    Real(dp), Intent(in) :: xext(zb_lo:zb_hi), aext(zb_lo:zb_hi), zext(zb_lo:zb_hi)

    ! Output variables
    Real(dp), Intent(out) :: ztilde(zb_lo:zb_hi), zinter(zb_lo:zb_hi)
    Real(dp), Intent(out) :: lambda0(zb_lo:zb_hi), gammae(zb_lo:zb_hi)
    Real(dp), Intent(out) :: dztildedt9(zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer :: izb
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    !__dir_enter_data &
    !__dir_async &
    !__dir_create(ztilde,zinter,lambda0,gammae,dztildedt9) &
    !__dir_create(ye,ytot,abar,zbar,z2bar,zibar,sratio) &
    !__dir_copyin(mask,t9,rho,y,etae,detaedt9,xext,aext,zext)

    ! Calculate Ye
    Call y_moment(y,ye(zb_lo:zb_hi),ytot(zb_lo:zb_hi), &
      & abar(zb_lo:zb_hi),zbar(zb_lo:zb_hi),z2bar(zb_lo:zb_hi),zibar(zb_lo:zb_hi), &
      & xext(zb_lo:zb_hi),aext(zb_lo:zb_hi),zext(zb_lo:zb_hi),mask_in = mask)

    ! Calculate ratio f'/f for electrons (Salpeter, Eq. 24; DGC, Eq. 5)
    Call salpeter_ratio(etae,sratio(zb_lo:zb_hi),dztildedt9,mask_in = mask)

    !__dir_loop_outer(1) &
    !__dir_async &
    !__dir_present(mask,t9,rho,y,etae,detaedt9) &
    !__dir_present(ztilde,zinter,lambda0,gammae,dztildedt9) &
    !__dir_present(ye,ytot,abar,zbar,z2bar,zibar,sratio)
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        ztilde(izb) = sqrt(z2bar(izb) + sratio(izb)*zbar(izb)) ! DGC, Eq. 4
        dztildedt9(izb) = 0.5*zbar(izb)/ztilde(izb) * dztildedt9(izb)*detaedt9(izb)

        ! Calculate plasma quantities
        Call plasma(t9(izb),rho(izb),ytot(izb),ye(izb),zbar(izb), &
          & zibar(izb),ztilde(izb),zinter(izb),lambda0(izb),gammae(izb))
      EndIf
    EndDo
    If ( idiag >= 3 ) Then
      !__dir_update &
      !__dir_wait &
      !__dir_host(t9,rho,ye,z2bar,zbar,sratio,ztilde,lambda0,gammae)
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          Write(lun_diag,"(a14,9es23.15)") 'EOS Screen', &
            & t9(izb),rho(izb),ye(izb),z2bar(izb),zbar(izb),sratio(izb), &
            & ztilde(izb),ztilde(izb)*lambda0(izb),gammae(izb)
        EndIf
      EndDo
    EndIf

    !__dir_exit_data &
    !__dir_async &
    !__dir_copyout(ztilde,zinter,lambda0,gammae,dztildedt9) &
    !__dir_copyout(ye,ytot,abar,zbar,z2bar,zibar,sratio) &
    !__dir_delete(mask,t9,rho,y,etae,detaedt9,xext,aext,zext)

    !__dir_wait

    Return
  End Subroutine eos_screen_vector

  Subroutine salpeter_ratio_scalar(eta,ratio,dratiodeta)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the Salpeter (1954) ratio f'/f(eta) needed for electron screening.
    ! eta is the ratio of electron chemical potential to kT.
    ! f'/f is also defined as theta_e in DeWitt+ (1973).
    !
    ! Calculation uses Fermi function relation d/dx f_(k+1) = (k+1) f_k and the rational function
    ! expansions of Fukushima (2015; AMC 259 708) for the F-D integrals of order 1/2, -1/2, and -3/2.
    !-----------------------------------------------------------------------------------------------
    !__dir_routine_seq
    Use fd, Only: fdm1h, fd1h, fdm3h
    Use xnet_controls, Only: iheat
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: eta

    ! Output variables
    Real(dp), Intent(out) :: ratio, dratiodeta

    ! Local variables
    Real(dp) :: fermip, fermim
    Real(dp) :: dfmdeta, dfpdeta

    ! Calculate f_(-1/2) and f_(1/2)
    fermim = fdm1h(eta)
    fermip = fd1h(eta)

    ! Evalutate the Salpeter ratio (extra factor of 1/2 from FD integral definition)
    ratio = 0.5 * fermim/fermip
    If ( iheat > 0 ) Then
      dfmdeta = -0.5 * fdm3h(eta)
      dfpdeta = +0.5 * fermim
      dratiodeta = ratio * (dfmdeta/fermim - dfpdeta/fermip)
    Else
      dratiodeta = 0.0
    EndIf

    Return
  End Subroutine salpeter_ratio_scalar

  Subroutine salpeter_ratio_vector(eta,ratio,dratiodeta,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the Salpeter (1954) ratio f'/f(eta) needed for electron screening.
    ! eta is the ratio of electron chemical potential to kT.
    ! f'/f is also defined as theta_e in DeWitt+ (1973).
    !
    ! Calculation uses Fermi function relation d/dx f_(k+1) = (k+1) f_k and the rational function
    ! expansions of Fukushima (2015; AMC 259 708) for the F-D integrals of order 1/2, -1/2, and -3/2.
    !-----------------------------------------------------------------------------------------------
    Use xnet_controls, Only: zb_lo, zb_hi, lzactive
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: eta(zb_lo:zb_hi)

    ! Output variables
    Real(dp), Intent(out) :: ratio(zb_lo:zb_hi), dratiodeta(zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer :: izb
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    !__dir_enter_data &
    !__dir_async &
    !__dir_create(ratio,dratiodeta) &
    !__dir_copyin(mask,eta)

    !__dir_loop_outer(1) &
    !__dir_async &
    !__dir_present(ratio,dratiodeta) &
    !__dir_present(mask,eta)
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        Call salpeter_ratio_scalar(eta(izb),ratio(izb),dratiodeta(izb))
      EndIf
    EndDo

    !__dir_exit_data &
    !__dir_async &
    !__dir_copyout(ratio,dratiodeta) &
    !__dir_delete(mask,eta)

    !__dir_wait

    Return
  End Subroutine salpeter_ratio_vector

End Module xnet_eos
