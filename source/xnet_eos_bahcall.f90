!***************************************************************************************************
! eos_bahcall.f90 10/18/17
! EoS Replacement based on Bahcall ().
! This file contains routines which calculate EoS quantites needed to calculate screening
! corrections for reaction rates.
!***************************************************************************************************

Module xnet_eos
  Implicit None

Contains

  Subroutine eos_initialize
    !-----------------------------------------------------------------------------------------------
    ! This routine initializes the EoS.
    !-----------------------------------------------------------------------------------------------
    Implicit None
    Return
  End Subroutine eos_initialize

  Subroutine eos_interface(rho,t9,y,ye,cv,efermkt,defermktdt9)
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
    Use nuclear_data, Only: zz53
    Use xnet_constants, Only: asig, avn, bok, e2, pi, pi2, third, emass, clt, ele_en, hbar
    Use xnet_controls, Only: idiag, iheat, iscrn, lun_diag
    Use xnet_types, Only: dp
    Use xnet_abundances, Only: y_moment
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9, rho, y(:)

    ! Output variables
    Real(dp), Intent(out) :: ye, cv, efermkt, defermktdt9

    ! Local variables
    Real(dp), Parameter :: a1 = 0.898004, b1 = 0.96786, c1 = 0.220703, d1 = -0.86097
    Real(dp), Parameter :: a2 = -0.86602540378, b2 = 0.29561, c2 = 1.9885
    Real(dp) :: ytot, bkt, abar, zbar, z2bar, zibar, z53bar
    Real(dp) :: etae_mb, ae, gam, gam14, gam32, gamc2, rel_ef
    Real(dp) :: eion, erad, ecoul
    Real(dp) :: deiondt9, deraddt9, decouldt9

    ! Calculate Ye and other needed moments of the abundance distribution
    Call y_moment(y,ye,ytot,abar,zbar,z2bar,zibar)

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
      z53bar = sum(zz53 * y) / ytot
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
    Else
      efermkt = 0.0
      defermktdt9 = 0.0
      cv = 0.0
    EndIf

    If ( idiag >= 3 ) Write(lun_diag,"(a,6es23.15)") 'EOS',t9,rho,ye,cv,efermkt,defermktdt9

    Return
  End Subroutine eos_interface

  Subroutine eos_screen(t9,rho,y,efermkt,defermktdt9,ztilde,zinter,lambda0,gammae,dztildedt9)
    !-----------------------------------------------------------------------------------------------
    ! This routine Bahcall's approach to calculate the factors needed for screening from the input
    ! temperature, density and composition.
    !
    ! Note: the efermkt and defermktdt9 in the subroutine interface are stored in XNet's conditions
    ! module as etae and detaedt9.  For more general EoS, etae and detaedt9 are their proper selves.
    !-----------------------------------------------------------------------------------------------
    Use xnet_constants, Only: avn, bok, clt, e2, ele_en, emass, hbar, pi, pi2, third, two3rd, &
      & thbim2, twm2bi
    Use xnet_controls, Only: idiag, iheat, lun_diag
    Use xnet_types, Only: dp
    Use xnet_abundances, Only: y_moment
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9, rho, y(:), efermkt, defermktdt9

    ! Output variables
    Real(dp), Intent(out) :: ztilde, zinter, lambda0, gammae, dztildedt9

    ! Local variables
    Real(dp) :: ye, ytot, bkt, abar, zbar, z2bar, zibar
    Real(dp) :: sratio, ae, dsratiodefermkt

    ! Calculate Ye and other needed moments of the abundance distribution
    Call y_moment(y,ye,ytot,abar,zbar,z2bar,zibar)

    ! Calculate ratio f'/f for electrons (Salpeter, Eq. 24)
    Call salpeter_ratio(efermkt,sratio,dsratiodefermkt)
    ztilde = sqrt(z2bar + zbar*sratio)
    If ( iheat > 0 ) Then
      dztildedt9 = 0.5*zbar/ztilde * dsratiodefermkt*defermktdt9
    Else
      dztildedt9 = 0.0
    EndIf

    ! Calculate plasma quantities
    bkt = bok*t9
    lambda0 = sqrt(4.0*pi*rho*avn*ytot) * (e2/bkt)**1.5 ! DGC, Eq. 3
    ae = (3.0 / (4.0*pi*avn*rho*ye))**third ! electron-sphere radius
    gammae = e2 / (ae*bkt) ! electron Coulomb coupling parameter
    zinter = zibar / (ztilde**thbim2 * zbar**twm2bi)
    If ( idiag >= 2 ) Write(lun_diag,"(a14,9es23.15)") 'EOS Screen', &
      & t9,rho,ye,z2bar,zbar,sratio,ztilde,ztilde*lambda0,gammae

    Return
  End Subroutine eos_screen

  Subroutine salpeter_ratio(efmkt,ratio,dratiodefmkt)
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
  End Subroutine salpeter_ratio

End Module xnet_eos
