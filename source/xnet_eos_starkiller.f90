!***************************************************************************************************
! eos_starkiller.f90 10/18/17
! Interface to starkiller
! This file contains routines which calculate EoS quantites needed to calculate screening
! corrections for reaction rates.
!***************************************************************************************************

Module xnet_eos
  Use eos_type_module, Only: eos_t
  Implicit None
  Type(eos_t) :: eos_state
  !$omp threadprivate(eos_state)

Contains

  Subroutine eos_initialize
    !-----------------------------------------------------------------------------------------------
    ! This routine initializes starkiller
    !-----------------------------------------------------------------------------------------------
    Use actual_eos_module, Only: actual_eos_init
    Implicit None
    Call actual_eos_init()

    Return
  End Subroutine eos_initialize

  Subroutine eos_interface(t9,rho,y,ye,cv,etae,detaedt9)
    !-----------------------------------------------------------------------------------------------
    ! This routine updates the equation of state for changes in temperature and density.
    !-----------------------------------------------------------------------------------------------
    Use xnet_constants, Only: avn, epmev
    Use xnet_controls, Only: idiag, iheat, iscrn, lun_diag
    Use xnet_types, Only: dp
    Use xnet_abundances, Only: y_moment

    Use actual_eos_module, Only: actual_eos
    Use eos_type_module, Only: eos_input_rt
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9, rho, y(:)

    ! Ouput variables
    Real(dp), Intent(out) :: ye, cv, etae, detaedt9

    ! Local variables
    Real(dp) :: ytot, abar, zbar, z2bar, zibar

    ! Calculate Ye
    Call y_moment(y,ye,ytot,abar,zbar,z2bar,zibar)

    If ( iscrn > 0 .or. iheat > 0 ) Then

      ! Load input variables for the eos
      eos_state%rho = rho
      eos_state%T = t9*1e9
      eos_state%y_e = ye
      eos_state%abar = abar
      eos_state%zbar = ye*abar

      ! Call the eos
      Call actual_eos(eos_input_rt,eos_state)

      ! Convert units from ergs/g to MeV/nucleon and K to GK
      etae = eos_state%eta
      detaedt9 = eos_state%detadt * 1e9
      cv = eos_state%cv * 1e9/epmev/avn
    Else
      etae = 0.0
      detaedt9 = 0.0
      cv = 0.0
    EndIf

    If ( idiag >= 3 ) Write(lun_diag,"(a,6es23.15)") 'EOS',t9,rho,ye,cv,etae,detaedt9

    Return
  End Subroutine eos_interface

  Subroutine eos_screen(t9,rho,y,etae,detaedt9,ztilde,zinter,lambda0,gammae,dztildedt9)
    !-----------------------------------------------------------------------------------------------
    ! This routine uses the current composition and prior updates to the Equation of State to
    ! calculate the factors needed for screening.
    !-----------------------------------------------------------------------------------------------
    Use xnet_constants, Only: avn, bok, clt, e2, ele_en, emass, hbar, pi, pi2, third, two3rd, &
      & thbim2, twm2bi
    Use xnet_controls, Only: idiag, iheat, lun_diag
    Use xnet_types, Only: dp
    Use xnet_abundances, Only: y_moment
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9, rho, y(:), etae, detaedt9

    ! Output variables
    Real(dp), Intent(out) :: ztilde, zinter, lambda0, gammae, dztildedt9

    ! Local variables
    Real(dp) :: ye, ytot, bkt, abar, zbar, z2bar, zibar
    Real(dp) :: sratio, ae, dsratiodeta

    ! Calculate Ye and other needed moments of the abundance distribution
    Call y_moment(y,ye,ytot,abar,zbar,z2bar,zibar)

    ! Calculate ratio f'/f for electrons (Salpeter, Eq. 24)
    Call salpeter_ratio(etae,sratio,dsratiodeta)
    ztilde = sqrt(z2bar + zbar*sratio)
    If ( iheat > 0 ) Then
      dztildedt9 = 0.5*zbar/ztilde * dsratiodeta*detaedt9
    Else
      dztildedt9 = 0.0
    EndIf

    ! Calculate plasma quantities
    bkt = bok*t9
    lambda0 = sqrt(4.0*pi*rho*avn*ytot) * (e2/bkt)**1.5 ! DGC, Eq. 3
    ae = (3.0 / (4.0*pi*avn*rho*ye))**third ! electron-sphere radius
    gammae = e2 / (ae*bkt) ! electron Coulomb coupling parameter
    zinter = zibar / (ztilde**thbim2 * zbar**twm2bi)
    If ( idiag >= 3 ) Write(lun_diag,"(a14,9es23.15)") 'EOS Screen', &
      & t9,rho,ye,z2bar,zbar,sratio,ztilde,ztilde*lambda0,gammae

    Return
  End Subroutine eos_screen

  Subroutine salpeter_ratio(eta,ratio,dratiodeta)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the Salpeter (1954) ratio f'/f(eta) needed for electron screening.
    ! eta is the ratio of electron chemical potential to kT.
    !
    ! Calculation uses Fermi function relation d/dx f_(k+1) = (k+1) f_k and the rational function
    ! expansions of Fukushima (2015; AMC 259 708) for the F-D integrals of order 1/2, -1/2, and -3/2.
    !-----------------------------------------------------------------------------------------------
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

    ! Evalutate the Salpeter ratio
    ratio = 0.5 * fermim/fermip
    If ( iheat > 0 ) Then
      dfmdeta = -0.5 * fdm3h(eta)
      dfpdeta = +0.5 * fermim
      dratiodeta = ratio * (dfmdeta/fermim - dfpdeta/fermip)
    Else
      dratiodeta = 0.0
    EndIf

    Return
  End Subroutine salpeter_ratio

End Module xnet_eos
