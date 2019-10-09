!***************************************************************************************************
! eos_starkiller.f90 10/18/17
! Interface to starkiller
! This file contains routines which calculate EoS quantites needed to calculate screening
! corrections for reaction rates.
!***************************************************************************************************

Module xnet_eos
  Use xnet_types, Only: dp
  Implicit None
  Real(dp), Allocatable :: ye(:), ytot(:), abar(:), zbar(:), z2bar(:), zibar(:), sratio(:)

Contains

  Subroutine eos_initialize
    !-----------------------------------------------------------------------------------------------
    ! This routine initializes starkiller
    !-----------------------------------------------------------------------------------------------
    Use xnet_controls, Only: nzevolve, tid

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

    !$acc enter data async(tid) &
    !$acc copyin(ye,ytot,abar,zbar,z2bar,zibar,sratio)

    Return
  End Subroutine eos_initialize

  Subroutine eos_interface(t9,rho,y,ye,cv,etae,detaedt9)
    !-----------------------------------------------------------------------------------------------
    ! This routine updates the equation of state for changes in temperature and density.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use xnet_abundances, Only: y_moment
    Use xnet_constants, Only: avn, epmev, amu
    Use xnet_controls, Only: idiag, iheat, iscrn, lun_diag
    Use xnet_types, Only: dp

    Use actual_eos_module, Only: xnet_actual_eos
    Use eos_type_module, Only: eos_input_rt, eos_t
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9, rho, y(ny)

    ! Ouput variables
    Real(dp), Intent(out) :: ye, cv, etae, detaedt9

    ! Local variables
    Real(dp) :: ytot, abar, zbar, z2bar, zibar
    Type(eos_t) :: eos_state

    ! Calculate Ye
    Call y_moment(y,ye,ytot,abar,zbar,z2bar,zibar)

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
      etae = eos_state%eta
      detaedt9 = eos_state%detadt * 1e9
      cv = eos_state%cv * amu * 1e9
    Else
      etae = 0.0
      detaedt9 = 0.0
      cv = 0.0
    EndIf

    If ( idiag >= 3 ) Write(lun_diag,"(a,6es24.16)") 'EOS',t9,rho,ye,cv,etae,detaedt9

    Return
  End Subroutine eos_interface

  Subroutine eos_interface2(t9,rho,y,yeout,cv,etae,detaedt9,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine updates the equation of state for changes in temperature and density.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use xnet_abundances, Only: y_moment2
    Use xnet_constants, Only: avn, epmev, amu
    Use xnet_controls, Only: idiag, iheat, iscrn, lun_diag, zb_lo, zb_hi, lzactive, tid
    Use xnet_types, Only: dp

    Use actual_eos_module, Only: xnet_actual_eos
    Use eos_type_module, Only: eos_input_rt, eos_t
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9(zb_lo:zb_hi), rho(zb_lo:zb_hi), y(ny,zb_lo:zb_hi)

    ! Ouput variables
    Real(dp), Intent(out) :: yeout(zb_lo:zb_hi), cv(zb_lo:zb_hi)
    Real(dp), Intent(out) :: etae(zb_lo:zb_hi), detaedt9(zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Type(eos_t) :: eos_state
    Integer :: izb
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    !$acc enter data async(tid) &
    !$acc copyin(mask)

    ! Calculate Ye
    Call y_moment2(y,ye(zb_lo:zb_hi),ytot(zb_lo:zb_hi),abar(zb_lo:zb_hi),zbar(zb_lo:zb_hi), &
      & z2bar(zb_lo:zb_hi),zibar(zb_lo:zb_hi),mask_in = mask)

    If ( iscrn > 0 .or. iheat > 0 ) Then

      !$acc parallel loop gang async(tid) &
      !$acc present(mask,t9,rho,y,ye,abar,zbar,cv,etae,detaedt9,yeout) &
      !$acc private(eos_state)
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then

          ! Load input variables for the eos
          eos_state%rho = rho(izb)
          eos_state%T = t9(izb)*1e9
          eos_state%y_e = ye(izb)
          eos_state%abar = abar(izb)
          eos_state%zbar = zbar(izb)

          ! Call the eos
          Call xnet_actual_eos(eos_input_rt,eos_state)

          ! Convert units from ergs/g to MeV/nucleon and K to GK
          etae(izb) = eos_state%eta
          detaedt9(izb) = eos_state%detadt * 1e9
          cv(izb) = eos_state%cv * amu * 1e9
          yeout(izb) = ye(izb)
        EndIf
      EndDo
    Else

      !$acc parallel loop gang async(tid) &
      !$acc present(mask,cv,etae,detaedt9,ye,yeout)
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          etae(izb) = 0.0
          detaedt9(izb) = 0.0
          cv(izb) = 0.0
          yeout(izb) = ye(izb)
        EndIf
      EndDo
    EndIf
    If ( idiag >= 3 ) Then
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          !$acc update wait(tid) &
          !$acc host(t9(izb),rho(izb),ye(izb),cv(izb),etae(izb),detaedt9(izb))
          Write(lun_diag,"(a,6es24.16)") 'EOS',t9(izb),rho(izb),ye(izb),cv(izb),etae(izb),detaedt9(izb)
        EndIf
      EndDo
      Flush(lun_diag)
    EndIf

    !$acc exit data async(tid) &
    !$acc delete(mask)

    Return
  End Subroutine eos_interface2

  Subroutine eos_screen(t9,rho,y,etae,detaedt9,ztilde,zinter,lambda0,gammae,dztildedt9)
    !-----------------------------------------------------------------------------------------------
    ! This routine uses the current composition and prior updates to the Equation of State to
    ! calculate the factors needed for screening.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use xnet_abundances, Only: y_moment
    Use xnet_constants, Only: avn, bok, clt, e2, ele_en, emass, hbar, pi, pi2, third, two3rd, &
      & thbim2, twm2bi
    Use xnet_controls, Only: idiag, iheat, lun_diag
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9, rho, y(ny), etae, detaedt9

    ! Output variables
    Real(dp), Intent(out) :: ztilde, zinter, lambda0, gammae, dztildedt9

    ! Local variables
    Real(dp) :: ye, ytot, abar, zbar, z2bar, zibar
    Real(dp) :: sratio, ae, bkt, nb, ni, ne

    ! Calculate Ye and other needed moments of the abundance distribution
    Call y_moment(y,ye,ytot,abar,zbar,z2bar,zibar)

    ! Calculate ratio f'/f for electrons (Salpeter, Eq. 24; DGC, Eq. 5)
    Call salpeter_ratio(etae,sratio,dztildedt9)
    ztilde = sqrt(z2bar + zbar*sratio) ! DGC, Eq. 4
    If ( iheat > 0 ) Then
      dztildedt9 = 0.5*zbar/ztilde * dztildedt9*detaedt9
    EndIf

    ! Calculate plasma quantities
    bkt = bok*t9
    nb = avn*rho
    ni = nb*ytot
    ne = nb*ye
    lambda0 = sqrt(4.0*pi*ni) * (e2/bkt)**1.5 ! DGC, Eq. 3
    ae = (3.0 / (4.0*pi*ne))**third ! electron-sphere radius
    gammae = e2 / (ae*bkt) ! electron Coulomb coupling parameter
    zinter = zibar / (ztilde**thbim2 * zbar**twm2bi) ! GDC, Table 4
    If ( idiag >= 3 ) Write(lun_diag,"(a14,9es24.16)") 'EOS Screen', &
      & t9,rho,ye,z2bar,zbar,sratio,ztilde,ztilde*lambda0,gammae

    Return
  End Subroutine eos_screen

  Subroutine eos_screen2(t9,rho,y,etae,detaedt9,ztilde,zinter,lambda0,gammae,dztildedt9,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine uses the current composition and prior updates to the Equation of State to
    ! calculate the factors needed for screening.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use xnet_abundances, Only: y_moment2
    Use xnet_constants, Only: avn, bok, clt, e2, ele_en, emass, hbar, pi, pi2, third, two3rd, &
      & thbim2, twm2bi
    Use xnet_controls, Only: idiag, iheat, lun_diag, zb_lo, zb_hi, lzactive, tid
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9(zb_lo:zb_hi), rho(zb_lo:zb_hi), y(ny,zb_lo:zb_hi)
    Real(dp), Intent(in) :: etae(zb_lo:zb_hi), detaedt9(zb_lo:zb_hi)

    ! Output variables
    Real(dp), Intent(out) :: ztilde(zb_lo:zb_hi), zinter(zb_lo:zb_hi)
    Real(dp), Intent(out) :: lambda0(zb_lo:zb_hi), gammae(zb_lo:zb_hi)
    Real(dp), Intent(out) :: dztildedt9(zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Real(dp) :: ae, bkt, nb, ni, ne
    Integer :: izb
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    !$acc enter data async(tid) &
    !$acc copyin(mask)

    ! Calculate ratio f'/f for electrons (Salpeter, Eq. 24; DGC, Eq. 5)
    Call salpeter_ratio2(etae,sratio(zb_lo:zb_hi),dztildedt9,mask_in = mask)

    !$acc parallel loop gang async(tid) &
    !$acc present(mask,t9,rho,y,etae,detaedt9,ztilde,zinter,lambda0,gammae,dztildedt9, &
    !$acc         ye,ytot,abar,zbar,z2bar,zibar,sratio) &
    !$acc private(ae,bkt,nb,ni,ne)
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        ztilde(izb) = sqrt(z2bar(izb) + sratio(izb)*zbar(izb)) ! DGC, Eq. 4
        If ( iheat > 0 ) Then
          dztildedt9(izb) = 0.5*zbar(izb)/ztilde(izb) * dztildedt9(izb)*detaedt9(izb)
        EndIf

        ! Calculate plasma quantities
        bkt = bok*t9(izb)
        nb = avn*rho(izb)
        ni = nb*ytot(izb)
        ne = nb*ye(izb)
        lambda0(izb) = sqrt(4.0*pi*ni) * (e2/bkt)**1.5 ! DGC, Eq. 3
        ae = (3.0 / (4.0*pi*ne))**third ! electron-sphere radius
        gammae(izb) = e2 / (ae*bkt) ! electron Coulomb coupling parameter
        zinter(izb) = zibar(izb) / (ztilde(izb)**thbim2 * zbar(izb)**twm2bi) ! GDC, Table 4
      EndIf
    EndDo
    If ( idiag >= 3 ) Then
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          !$acc update wait(tid) &
          !$acc host(t9(izb),rho(izb),ye(izb),zbar(izb),z2bar(izb),ztilde(izb), &
          !$acc      lambda0(izb),gammae(izb),sratio(izb))
          Write(lun_diag,"(a14,9es24.16)") 'EOS Screen', &
            & t9(izb),rho(izb),ye(izb),z2bar(izb),zbar(izb),sratio(izb), &
            & ztilde(izb),ztilde(izb)*lambda0(izb),gammae(izb)
        EndIf
      EndDo
      Flush(lun_diag)
    EndIf

    !$acc exit data async(tid) &
    !$acc delete(mask)

    Return
  End Subroutine eos_screen2

  Subroutine salpeter_ratio(eta,ratio,dratiodeta)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the Salpeter (1954) ratio f'/f(eta) needed for electron screening.
    ! eta is the ratio of electron chemical potential to kT.
    ! f'/f is also defined as theta_e in DeWitt+ (1973).
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
  End Subroutine salpeter_ratio

  Subroutine salpeter_ratio2(eta,ratio,dratiodeta,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the Salpeter (1954) ratio f'/f(eta) needed for electron screening.
    ! eta is the ratio of electron chemical potential to kT.
    ! f'/f is also defined as theta_e in DeWitt+ (1973).
    !
    ! Calculation uses Fermi function relation d/dx f_(k+1) = (k+1) f_k and the rational function
    ! expansions of Fukushima (2015; AMC 259 708) for the F-D integrals of order 1/2, -1/2, and -3/2.
    !-----------------------------------------------------------------------------------------------
    Use fd, Only: fdm1h, fd1h, fdm3h
    Use xnet_controls, Only: iheat, zb_lo, zb_hi, lzactive, tid
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: eta(zb_lo:zb_hi)

    ! Output variables
    Real(dp), Intent(out) :: ratio(zb_lo:zb_hi), dratiodeta(zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Real(dp) :: fermip, fermim
    Real(dp) :: dfmdeta, dfpdeta
    Integer :: izb
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    !$acc enter data async(tid) &
    !$acc copyin(mask)

    !$acc parallel loop gang async(tid) &
    !$acc present(mask,eta,ratio,dratiodeta) &
    !$acc private(fermim,fermip,dfmdeta,dfpdeta)
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then

        ! Calculate f_(-1/2) and f_(1/2)
        fermim = fdm1h(eta(izb))
        fermip = fd1h(eta(izb))

        ! Evalutate the Salpeter ratio (extra factor of 1/2 from FD integral definition)
        ratio(izb) = 0.5 * fermim/fermip
        If ( iheat > 0 ) Then
          dfmdeta = -0.5 * fdm3h(eta(izb))
          dfpdeta = +0.5 * fermim
          dratiodeta(izb) = ratio(izb) * (dfmdeta/fermim - dfpdeta/fermip)
        Else
          dratiodeta(izb) = 0.0
        EndIf
      EndIf
    EndDo

    !$acc exit data async(tid) &
    !$acc delete(mask)

    Return
  End Subroutine salpeter_ratio2

End Module xnet_eos
