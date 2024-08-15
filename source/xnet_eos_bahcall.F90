!***************************************************************************************************
! eos_bahcall.f90 10/18/17
! EoS Replacement based on Bahcall ().
! This file contains routines which calculate EoS quantites needed to calculate screening
! corrections for reaction rates.
!***************************************************************************************************

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
  End Interface eos_screen

  Interface salpeter_ratio
    Module Procedure salpeter_ratio_scalar
    Module Procedure salpeter_ratio_vector
  End Interface

Contains

  Subroutine eos_initialize
    !-----------------------------------------------------------------------------------------------
    ! This routine initializes the EoS.
    !-----------------------------------------------------------------------------------------------
    Use xnet_controls, Only: nzevolve
    Implicit None

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

  Subroutine eosx(t9,rho,ye,ytot,abar,zbar,z53bar,cv,efermkt,defermktdt9)
    !-----------------------------------------------------------------------------------------------
    ! This routine interfaces with and calls the underlying EoS.
    !-----------------------------------------------------------------------------------------------
    Use xnet_constants, Only: asig, avn, bok, e2, pi, pi2, third, emass, clt, ele_en, hbar
    Use xnet_controls, Only: iheat, iscrn
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9, rho, ye, ytot, abar, zbar, z53bar

    ! Ouput variables
    Real(dp), Intent(out) :: cv, efermkt, defermktdt9

    ! Local variables
    Real(dp), Parameter :: a1 = 0.898004, b1 = 0.96786, c1 = 0.220703, d1 = -0.86097
    Real(dp), Parameter :: a2 = -0.86602540378, b2 = 0.29561, c2 = 1.9885
    Real(dp) :: bkt
    Real(dp) :: etae_mb, ae, gam, gam14, gam32, gamc2, rel_ef
    Real(dp) :: eion, erad, ecoul
    Real(dp) :: deiondt9, deraddt9, decouldt9

    cv = 0.0
    efermkt = 0.0
    defermktdt9 = 0.0
    If ( iscrn > 0 .or. iheat > 0 ) Then

      bkt = bok*t9
      etae_mb = log(rho*avn*ye * (2.0*pi*hbar**2 / (emass*bkt))**1.5)
      rel_ef = hbar * (3.0*pi2*rho*avn*ye)**third / (emass*clt)
      efermkt = ele_en * (sqrt(1.0 + rel_ef**2) - 1.0) / bkt
      defermktdt9 = -efermkt/t9

      ! Calculate energy derivatives (assume etot = eion + erad + ecoul) (ignore degenerate electron energy)
      eion = 1.5*bkt*ytot
      deiondt9 = 1.5*ytot*bok

      erad = asig*bkt*bkt*bkt*bkt / (rho*avn)
      deraddt9 = 4.0*erad / t9

      ae = (3.0 / (4.0*pi*avn*rho*ye))**third ! electron-sphere radius
      gam = z53bar*e2 / (ae*bkt) ! ionic Coulomb coupling parameter
      If ( gam >= 1.0 ) Then
        gam14 = gam**0.25
        ecoul = (a1*gam + b1*gam14 + c1/gam14 + d1) * bkt * ytot
        decouldt9 = ecoul/t9 - bok * ytot * (a1*gam + 0.25 * (b1*gam14 - c1/gam14))
      Else
        gam32 = gam*sqrt(gam)
        gamc2 = gam**c2
        ecoul = (a2*gam32 + b2*gamc2) * bkt * ytot
        decouldt9 = ecoul/t9 - bok * ytot * (1.5*a2*gam32 + b2*c2*gamc2)
      EndIf
      cv = deiondt9 + deraddt9 + decouldt9
    EndIf

  End Subroutine eosx

  Subroutine eos_interface_scalar(rho,t9,y,ye,cv,efermkt,defermktdt9,xext,aext,zext)
    !-----------------------------------------------------------------------------------------------
    ! This routine updates the Equation of State for changes in temperature and density.
    !
    ! Bahcall introduced a fit for the plasma parameters needed to calculate the screening factors
    ! for nuclear reaction based on the easier to calculate efermkt (the ratio of the electron
    ! Fermi Energy to the themal energy) instead of the more general expresions based on etae (the
    ! ratio of the electron chemical potential to the thermal energy), which generally require a
    ! complex Equation of State to calculate.   For basic cases, XNet uses Bahcall's approach
    ! to simplify the
    !
    ! Note: the factors efermkt and defermktdt9 (the derivative of efermkt wrt to temperature) in the
    ! subroutine interface are stored in XNet's conditions module as etae and detaedt9.  For more
    ! general EoS, etae and detaedt9 are their proper selves.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny, zz53
    Use xnet_abundances, Only: y_moment
    Use xnet_controls, Only: idiag, lun_diag
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9, rho, y(ny), xext, aext, zext

    ! Output variables
    Real(dp), Intent(out) :: ye, cv, efermkt, defermktdt9

    ! Local variables
    Real(dp) :: ytot, abar, zbar, z2bar, zibar, z53bar

    ! Calculate Ye
    Call y_moment(y,ye,ytot,abar,zbar,z2bar,zibar,xext,aext,zext)
    z53bar = sum(zz53 * y) / ytot

    ! Call the eos
    Call eosx(t9,rho,ye,ytot,abar,zbar,z53bar,cv,efermkt,defermktdt9)

    If ( idiag >= 3 ) Write(lun_diag,"(a,6es23.15)") 'EOS',t9,rho,ye,cv,efermkt,defermktdt9

    Return
  End Subroutine eos_interface_scalar

  Subroutine eos_interface_vector(t9,rho,y,ye,cv,efermkt,defermktdt9,xext,aext,zext,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine updates the Equation of State for changes in temperature and density.
    !
    ! Bahcall introduced a fit for the plasma parameters needed to calculate the screening factors
    ! for nuclear reaction based on the easier to calculate efermkt (the ratio of the electron
    ! Fermi Energy to the themal energy) instead of the more general expresions based on etae (the
    ! ratio of the electron chemical potential to the thermal energy), which generally require a
    ! complex Equation of State to calculate.   For basic cases, XNet uses Bahcall's approach
    ! to simplify the
    !
    ! Note: the factors efermkt and defermktdt9 (the derivative of efermkt wrt to temperature) in the
    ! subroutine interface are stored in XNet's conditions module as etae and detaedt9.  For more
    ! general EoS, etae and detaedt9 are their proper selves.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny, zz53
    Use xnet_abundances, Only: y_moment
    Use xnet_controls, Only: idiag, lun_diag, zb_lo, zb_hi, lzactive
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9(zb_lo:zb_hi), rho(zb_lo:zb_hi), y(ny,zb_lo:zb_hi)
    Real(dp), Intent(in) :: xext(zb_lo:zb_hi), aext(zb_lo:zb_hi), zext(zb_lo:zb_hi)

    ! Ouput variables
    Real(dp), Intent(out) :: ye(zb_lo:zb_hi), cv(zb_lo:zb_hi)
    Real(dp), Intent(out) :: efermkt(zb_lo:zb_hi), defermktdt9(zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer :: izb
    Real(dp) :: z53bar
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    ! Calculate Ye
    Call y_moment(y,ye,ytot(zb_lo:zb_hi), &
      & abar(zb_lo:zb_hi),zbar(zb_lo:zb_hi),z2bar(zb_lo:zb_hi),zibar(zb_lo:zb_hi), &
      & xext(zb_lo:zb_hi),aext(zb_lo:zb_hi),zext(zb_lo:zb_hi),mask_in = mask_in)

    ! Call the eos
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        z53bar = sum(zz53 * y(:,izb)) / ytot(izb)
        Call eosx(t9(izb),rho(izb),ye(izb),ytot(izb),abar(izb),zbar(izb), &
          & z53bar,cv(izb),efermkt(izb),defermktdt9(izb))
      EndIf
    EndDo

    If ( idiag >= 3 ) Then
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          Write(lun_diag,"(a,6es23.15)") 'EOS',t9(izb),rho(izb),ye(izb),cv(izb),efermkt(izb),defermktdt9(izb)
        EndIf
      EndDo
    EndIf

    Return
  End Subroutine eos_interface_vector

  Subroutine eos_screen_scalar(t9,rho,y,efermkt,defermktdt9,ztilde,zinter,lambda0,gammae,dztildedt9,xext,aext,zext)
    !-----------------------------------------------------------------------------------------------
    ! This routine Bahcall's approach to calculate the factors needed for screening from the input
    ! temperature, density and composition.
    !
    ! Note: the efermkt and defermktdt9 in the subroutine interface are stored in XNet's conditions
    ! module as etae and detaedt9.  For more general EoS, etae and detaedt9 are their proper selves.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use xnet_abundances, Only: y_moment
    Use xnet_controls, Only: idiag, lun_diag
    Use xnet_types, Only: dp
    Use xnet_util, Only: plasma
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9, rho, y(ny), efermkt, defermktdt9, xext, aext, zext

    ! Output variables
    Real(dp), Intent(out) :: ztilde, zinter, lambda0, gammae, dztildedt9

    ! Local variables
    Real(dp) :: ye, ytot, bkt, abar, zbar, z2bar, zibar
    Real(dp) :: sratio, ae, dsratiodefermkt

    ! Calculate Ye and other needed moments of the abundance distribution
    Call y_moment(y,ye,ytot,abar,zbar,z2bar,zibar,xext,aext,zext)

    ! Calculate ratio f'/f for electrons (Salpeter, Eq. 24)
    Call salpeter_ratio(efermkt,sratio,dztildedt9)
    ztilde = sqrt(z2bar + zbar*sratio)
    dztildedt9 = 0.5*zbar/ztilde * dztildedt9*defermktdt9

    ! Calculate plasma quantities
    Call plasma(t9,rho,ytot,ye,zbar,zibar,ztilde,zinter,lambda0,gammae)
    If ( idiag >= 3 ) Write(lun_diag,"(a14,9es23.15)") 'EOS Screen', &
      & t9,rho,ye,z2bar,zbar,sratio,ztilde,ztilde*lambda0,gammae

    Return
  End Subroutine eos_screen_scalar

  Subroutine eos_screen_vector(t9,rho,y,efermkt,defermktdt9,ztilde,zinter,lambda0,gammae,dztildedt9,xext,aext,zext,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine Bahcall's approach to calculate the factors needed for screening from the input
    ! temperature, density and composition.
    !
    ! Note: the efermkt and defermktdt9 in the subroutine interface are stored in XNet's conditions
    ! module as etae and detaedt9.  For more general EoS, etae and detaedt9 are their proper selves.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use xnet_abundances, Only: y_moment
    Use xnet_controls, Only: idiag, lun_diag, zb_lo, zb_hi, lzactive
    Use xnet_types, Only: dp
    Use xnet_util, Only: plasma
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9(zb_lo:zb_hi), rho(zb_lo:zb_hi), y(ny,zb_lo:zb_hi)
    Real(dp), Intent(in) :: efermkt(zb_lo:zb_hi), defermktdt9(zb_lo:zb_hi)
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

    ! Calculate Ye
    Call y_moment(y,ye(zb_lo:zb_hi),ytot(zb_lo:zb_hi), &
      & abar(zb_lo:zb_hi),zbar(zb_lo:zb_hi),z2bar(zb_lo:zb_hi),zibar(zb_lo:zb_hi), &
      & xext(zb_lo:zb_hi),aext(zb_lo:zb_hi),zext(zb_lo:zb_hi),mask_in = mask_in)

    ! Calculate ratio f'/f for electrons (Salpeter, Eq. 24; DGC, Eq. 5)
    Call salpeter_ratio(efermkt,sratio(zb_lo:zb_hi),dztildedt9,mask_in = mask_in)

    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        ztilde(izb) = sqrt(z2bar(izb) + sratio(izb)*zbar(izb)) ! DGC, Eq. 4
        dztildedt9(izb) = 0.5*zbar(izb)/ztilde(izb) * dztildedt9(izb)*defermktdt9(izb)

        ! Calculate plasma quantities
        Call plasma(t9(izb),rho(izb),ytot(izb),ye(izb),zbar(izb), &
          & zibar(izb),ztilde(izb),zinter(izb),lambda0(izb),gammae(izb))
      EndIf
    EndDo
    If ( idiag >= 3 ) Then
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          Write(lun_diag,"(a14,9es23.15)") 'EOS Screen', &
            & t9(izb),rho(izb),ye(izb),z2bar(izb),zbar(izb),sratio(izb), &
            & ztilde(izb),ztilde(izb)*lambda0(izb),gammae(izb)
        EndIf
      EndDo
    EndIf

    Return
  End Subroutine eos_screen_vector

  Subroutine salpeter_ratio_scalar(efmkt,ratio,dratiodefmkt)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the Salpeter (1954) ratio f'/f(eta) needed for electron screening, using
    ! a fit to Figure 24 in that paper. efmkt is the ratio of electron chemical potential to kT.
    !-----------------------------------------------------------------------------------------------
    Use xnet_constants, Only: ln_10
    Use xnet_controls, Only: iheat
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: efmkt

    ! Output variables
    Real(dp), Intent(out) :: ratio, dratiodefmkt

    ! Local variables
    Real(dp), Parameter :: a0 = 1.578347, c0 = 0.75793, c1 = -0.54621
    Real(dp), Parameter ::  c2 = -0.30964,  c3 = 0.12535,  c4 = 0.1203,  c5 = -0.012857,  c6 = -0.014768
    Real(dp), Parameter :: dc2 = -0.61928, dc3 = 0.37605, dc4 = 0.4812, dc5 = -0.064285, dc6 = -0.088608
    Real(dp) :: lefmkt

    lefmkt = log10(efmkt)
    If ( lefmkt <= -a0 ) Then ! Bahcall uses limit of -2, but this gives ratio slightly above 1.
      ratio = 1.0
    ElseIf ( lefmkt >= 1.5 ) Then
      ratio = 0.0
    Else
      ratio = c0 + (c1*lefmkt) + (c2*lefmkt**2) + (c3*lefmkt**3) + (c4*lefmkt**4) + (c5*lefmkt**5) + (c6*lefmkt**6)
    Endif

    If ( iheat > 0 .and. lefmkt > a0 .and. lefmkt < 1.5 ) Then
      dratiodefmkt = c1 + dc2*lefmkt + dc3*lefmkt**2 + dc4*lefmkt**3 + dc5*lefmkt**4 + dc6*lefmkt**5
      dratiodefmkt = dratiodefmkt / (ln_10*efmkt)
    Else
      dratiodefmkt = 0.0
    EndIf

    Return
  End Subroutine salpeter_ratio_scalar

  Subroutine salpeter_ratio_vector(efmkt,ratio,dratiodefmkt,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the Salpeter (1954) ratio f'/f(eta) needed for electron screening, using
    ! a fit to Figure 24 in that paper. efmkt is the ratio of electron chemical potential to kT.
    !-----------------------------------------------------------------------------------------------
    Use xnet_controls, Only: zb_lo, zb_hi, lzactive
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: efmkt(zb_lo:zb_hi)

    ! Output variables
    Real(dp), Intent(out) :: ratio(zb_lo:zb_hi), dratiodefmkt(zb_lo:zb_hi)

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

    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        Call salpeter_ratio_scalar(efmkt(izb),ratio(izb),dratiodefmkt(izb))
      EndIf
    EndDo

    Return
  End Subroutine salpeter_ratio_vector

End Module xnet_eos
