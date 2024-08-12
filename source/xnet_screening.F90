!***************************************************************************************************
! xnet_screening.f90 10/18/17
! This file contains the routines needed to calculate screening corrections for reaction rates.
!***************************************************************************************************

#include "xnet_macros.fh"

Module xnet_screening
  !-------------------------------------------------------------------------------------------------
  ! This module contains data and routines used to calculate the screening corrections that appear
  ! in the reaction rate exponents.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None
  Real(dp), Allocatable :: h1(:,:), h2(:,:), h3(:,:), h4(:,:)                 ! Screening factors
  Real(dp), Allocatable :: dh1dt9(:,:), dh2dt9(:,:), dh3dt9(:,:), dh4dt9(:,:) ! d(screening factors)/dT9

  ! Proton (charge) numbers for individual reactants in 2-, 3-, and 4-reactant reactions
  Real(dp), Allocatable :: z21(:), z22(:), z31(:), z32(:), z33(:)
  Real(dp), Allocatable :: z41(:), z42(:), z43(:), z44(:)
  Integer, Allocatable  :: iz21(:), iz22(:), iz31(:), iz32(:), iz33(:)
  Integer, Allocatable  :: iz41(:), iz42(:), iz43(:), iz44(:)

  ! Composite proton (charge) numbers for reactants in 2-, 3-, and 4-reactant reactions
  Real(dp), Allocatable :: z2c(:), z3c(:), z4c(:)
  Integer, Allocatable  :: iz2c(:), iz3c(:), iz4c(:)

  ! Screening factors from Table 4 of Graboske+ (1973)
  Real(dp), Allocatable :: zeta2w(:), zeta2i(:), zeta3w(:), zeta3i(:), zeta4w(:), zeta4i(:) ! Reaction charge parameter

  Real(dp), Allocatable :: ztilde(:), zinter(:), lambda0(:), gammae(:), dztildedt9(:)

  Real(dp), Parameter :: lam_1 = 0.1_dp
  Real(dp), Parameter :: lam_2 = 0.125_dp
  Real(dp), Parameter :: lam_3 = 2.0_dp
  Real(dp), Parameter :: lam_4 = 2.15_dp
  Real(dp), Parameter :: lam_5 = 4.85_dp
  Real(dp), Parameter :: lam_6 = 5.0_dp
  Real(dp), Parameter :: rdlam12 = 1.0_dp / (lam_1 - lam_2)
  Real(dp), Parameter :: rdlam34 = 1.0_dp / (lam_3 - lam_4)
  Real(dp), Parameter :: rdlam56 = 1.0_dp / (lam_6 - lam_5)

Contains

  Subroutine screening_init
    !-----------------------------------------------------------------------------------------------
    ! This routine allocates and initializes screening arrays
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: zz
    Use reaction_data, Only: nreac, n2i, n3i, n4i
    Use xnet_constants, Only: bip1, thbim1, five3rd
    Use xnet_controls, Only: iheat, iscrn, nzevolve, tid
    Use xnet_types, Only: dp
    Implicit None

    ! Local variables
    Integer :: i, nr1, nr2, nr3, nr4

    nr1 = nreac(1)
    nr2 = nreac(2)
    nr3 = nreac(3)
    nr4 = nreac(4)

    If ( iscrn > 0 ) Then

      ! 2-reactant screening terms
      Allocate (z21(nr2),z22(nr2))
      Allocate (iz21(nr2),iz22(nr2))
      Allocate (iz2c(nr2),z2c(nr2),zeta2w(nr2),zeta2i(nr2))
      z21 = zz(n2i(1,:))
      z22 = zz(n2i(2,:))
      iz21 = nint(z21)
      iz22 = nint(z22)
      iz2c = iz21 + iz22
      z2c = real(iz2c,dp)
      zeta2w = z21*z22
      zeta2i = z2c**bip1 - z21**bip1 - z22**bip1

      ! 3-reactant screening terms
      Allocate (z31(nr3),z32(nr3),z33(nr3))
      Allocate (iz31(nr3),iz32(nr3),iz33(nr3))
      Allocate (iz3c(nr3),z3c(nr3),zeta3w(nr3),zeta3i(nr3))
      z31 = zz(n3i(1,:))
      z32 = zz(n3i(2,:))
      z33 = zz(n3i(3,:))
      iz31 = nint(z31)
      iz32 = nint(z32)
      iz33 = nint(z33)
      iz3c = iz31 + iz32 + iz33
      z3c = real(iz3c,dp)
      zeta3w = z31*z32 + z31*z33 + z32*z33
      zeta3i = z3c**bip1 - z31**bip1 - z32**bip1 - z33**bip1

      ! 4-reactant screening terms
      Allocate (z41(nr4),z42(nr4),z43(nr4),z44(nr4))
      Allocate (iz41(nr4),iz42(nr4),iz43(nr4),iz44(nr4))
      Allocate (iz4c(nr4),z4c(nr4),zeta4w(nr4),zeta4i(nr4))
      z41 = zz(n4i(1,:))
      z42 = zz(n4i(2,:))
      z43 = zz(n4i(3,:))
      z44 = zz(n4i(4,:))
      iz41 = nint(z41)
      iz42 = nint(z42)
      iz43 = nint(z43)
      iz44 = nint(z44)
      iz4c = iz41 + iz42 + iz43 + iz44
      z4c = real(iz4c,dp)
      zeta4w = z41*(z42 + z43 + z44) + z42*(z43 + z44) + z43 * z44
      zeta4i = z4c**bip1 - z41**bip1 - z42**bip1 - z43**bip1 - z44**bip1

      Allocate (ztilde(nzevolve))
      Allocate (zinter(nzevolve))
      Allocate (lambda0(nzevolve))
      Allocate (gammae(nzevolve))
      Allocate (dztildedt9(nzevolve))
      ztilde = 0.0
      zinter = 0.0
      lambda0 = 0.0
      gammae = 0.0
      dztildedt9 = 0.0
      !__dir_enter_data &
      !__dir_async(tid) &
      !__dir_copyin(iz21,iz22,iz31,iz32,iz33,iz41,iz42,iz43,iz44) &
      !__dir_copyin(iz2c,iz3c,iz4c,zeta2w,zeta3w,zeta4w,zeta2i,zeta3i,zeta4i) &
      !__dir_copyin(ztilde,zinter,lambda0,gammae,dztildedt9)
    EndIf

    Allocate (h1(nr1,nzevolve))
    Allocate (h2(nr2,nzevolve))
    Allocate (h3(nr3,nzevolve))
    Allocate (h4(nr4,nzevolve))
    Allocate (dh1dt9(nr1,nzevolve))
    Allocate (dh2dt9(nr2,nzevolve))
    Allocate (dh3dt9(nr3,nzevolve))
    Allocate (dh4dt9(nr4,nzevolve))
    h1 = 0.0
    h2 = 0.0
    h3 = 0.0
    h4 = 0.0
    dh1dt9 = 0.0
    dh2dt9 = 0.0
    dh3dt9 = 0.0
    dh4dt9 = 0.0

    !__dir_enter_data &
    !__dir_async(tid) &
    !__dir_copyin(h1,h2,h3,h4,dh1dt9,dh2dt9,dh3dt9,dh4dt9)

    Return
  End Subroutine screening_init

  Subroutine screening(mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the screening factors necessary for XNet. An equation of state,
    ! typically Helmholtz (Timmes & Swesty 1999) is used to determine the electron distribution and
    ! chemical potential.
    !
    ! References:
    ! Weak Screening:
    !    Salpeter (1954) Aust J Phys 7 353.
    ! Intermediate Screening:
    !    DeWitt, Graboske & Cooper (1973) ApJ 181 439.
    !    Graboske, DeWitt, Grossman & Cooper (1973) ApJ 181 457.
    ! Strong Screening:
    !    DeWitt & Slattery (2003) Contrib Plasma Phys 43 279.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: izmax, nname, zseq, zseq53, zseqi
    Use reaction_data, Only: n2i, n3i, n4i, nreac
    Use xnet_abundances, Only: yt, xext, aext, zext
    Use xnet_constants, Only: bi, bip1, cds, kbi, thbim2
    Use xnet_conditions, Only: rhot, t9t, etae, detaedt9
    Use xnet_controls, Only: idiag, iheat, iscrn, lun_diag, szbatch, zb_lo, zb_hi, lzactive, tid
    Use xnet_eos, Only: eos_screen
    Use xnet_timers, Only: xnet_wtime, start_timer, stop_timer, timer_prescrn, timer_scrn
    Use xnet_types, Only: dp
    Implicit None

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Real(dp) :: fhs(0:izmax+2), dfhsdt9(0:izmax+2)
    Real(dp) :: gammaz, gammaz5, lgammaz, lambda
    Real(dp) :: hw0, hi0, dlnhw0dt9, dlnhi0dt9
    Real(dp) :: h, hw, hi, hs
    Real(dp) :: dhdt9, dhwdt9, dhidt9, dhsdt9
    Integer :: j, mu, izb, izone
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    start_timer = xnet_wtime()
    timer_prescrn = timer_prescrn - start_timer

    !__dir_enter_data &
    !__dir_async(tid) &
    !__dir_copyin(mask)

    ! Call EOS to get plasma quantities
    Call eos_screen(t9t(zb_lo:zb_hi),rhot(zb_lo:zb_hi),yt(:,zb_lo:zb_hi),etae(zb_lo:zb_hi), &
      & detaedt9(zb_lo:zb_hi),ztilde(zb_lo:zb_hi),zinter(zb_lo:zb_hi),lambda0(zb_lo:zb_hi), &
      & gammae(zb_lo:zb_hi),dztildedt9(zb_lo:zb_hi),xext(zb_lo:zb_hi),aext(zb_lo:zb_hi), &
      & zext(zb_lo:zb_hi),mask_in = mask)

    stop_timer = xnet_wtime()
    timer_prescrn = timer_prescrn + stop_timer

    start_timer = xnet_wtime()
    timer_scrn = timer_scrn - start_timer

    !__dir_loop_outer(1) &
    !__dir_async(tid) &
    !__dir_present(nreac,zseq,zseq53,zseqi) &
    !__dir_present(iz21,iz22,iz31,iz32,iz33,iz41,iz42,iz43,iz44) &
    !__dir_present(iz2c,iz3c,iz4c,zeta2w,zeta3w,zeta4w,zeta2i,zeta3i,zeta4i) &
    !__dir_present(h1,h2,h3,h4,dh1dt9,dh2dt9,dh3dt9,dh4dt9) &
    !__dir_present(ztilde,zinter,lambda0,gammae,dztildedt9) &
    !__dir_present(mask,t9t) &
    !__dir_private(fhs,dfhsdt9,hw0,hi0,dlnhw0dt9,dlnhi0dt9)
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then

        ! Calculate screening energies as a function of Z, for prescriptions that follow this approach
        fhs(0) = 0.0
        dfhsdt9(0) = 0.0
        !__dir_loop_inner(1) &
        !__dir_private(gammaz,gammaz5,lgammaz)
        Do j = 1, izmax+2
          gammaz = gammae(izb) * zseq53(j)
          gammaz5 = gammaz**cds(5)
          lgammaz = log(gammaz)
          fhs(j) = cds(1)*gammaz + cds(2)/cds(5)*gammaz5 + cds(3)*lgammaz + cds(4)
          dfhsdt9(j) = -(cds(1)*gammaz + cds(2)*gammaz5 + cds(3))/t9t(izb)
        EndDo

        ! No screening term for 1-reactant reactions
        !__dir_loop_inner(1)
        Do mu = 1, nreac(1)
          h1(mu,izb) = 0.0
          dh1dt9(mu,izb) = 0.0
        EndDo

        hw0 = ztilde(izb) * lambda0(izb)
        hi0 = kbi * zinter(izb) * lambda0(izb)**bi
        dlnhw0dt9 = dztildedt9(izb)/ztilde(izb) - 1.5/t9t(izb)
        dlnhi0dt9 = - thbim2*dztildedt9(izb)/ztilde(izb) - 1.5*bi/t9t(izb)

        ! 2-reactant screening
        !__dir_loop_inner(1) &
        !__dir_private(lambda,h,hw,hi,hs,dhdt9,dhwdt9,dhidt9,dhsdt9)
        Do mu = 1, nreac(2)
          If ( zeta2w(mu) > 0.0 ) Then

            lambda = hw0 * zeta2w(mu)
            hw = lambda
            hi = hi0 * zeta2i(mu)
            hs = fhs(iz21(mu)) + fhs(iz22(mu)) - fhs(iz2c(mu))

            dhwdt9 = hw * dlnhw0dt9
            dhidt9 = hi * dlnhi0dt9
            dhsdt9 = dfhsdt9(iz21(mu)) + dfhsdt9(iz22(mu)) - dfhsdt9(iz2c(mu))

            ! Select Screening factor for 2 reactant reactions
            Call screen_blend(lambda,hw,hi,hs,dhwdt9,dhidt9,dhsdt9,h,dhdt9)
            h2(mu,izb) = h
            dh2dt9(mu,izb) = dhdt9
          Else
            h2(mu,izb) = 0.0
            dh2dt9(mu,izb) = 0.0
          EndIf
        EndDo

        ! 3-reactant screening
        !__dir_loop_inner(1) &
        !__dir_private(lambda,h,hw,hi,hs,dhdt9,dhwdt9,dhidt9,dhsdt9)
        Do mu = 1, nreac(3)
          If ( zeta3w(mu) > 0.0 ) Then

            lambda = hw0 * zeta3w(mu)
            hw = lambda
            hi = hi0 * zeta3i(mu)
            hs = fhs(iz31(mu)) + fhs(iz33(mu)) - fhs(iz3c(mu))

            dhwdt9 = hw * dlnhw0dt9
            dhidt9 = hi * dlnhi0dt9
            dhsdt9 = dfhsdt9(iz31(mu)) + dfhsdt9(iz32(mu)) + dfhsdt9(iz33(mu)) - dfhsdt9(iz3c(mu))

            ! Select Screening factor for 3 reactant reactions
            Call screen_blend(lambda,hw,hi,hs,dhwdt9,dhidt9,dhsdt9,h,dhdt9)
            h3(mu,izb) = h
            dh3dt9(mu,izb) = dhdt9
          Else
            h3(mu,izb) = 0.0
            dh3dt9(mu,izb) = 0.0
          EndIf
        EndDo

        ! 4-reactant screening
        !__dir_loop_inner(1) &
        !__dir_private(lambda,h,hw,hi,hs,dhdt9,dhwdt9,dhidt9,dhsdt9)
        Do mu = 1, nreac(4)
          If ( zeta4w(mu) > 0.0 ) Then

            lambda = hw0 * zeta4w(mu)
            hw = lambda
            hi = hi0 * zeta4i(mu)
            hs = fhs(iz41(mu)) + fhs(iz42(mu)) + fhs(iz43(mu)) + fhs(iz44(mu)) - fhs(iz4c(mu))

            dhwdt9 = hw * dlnhw0dt9
            dhidt9 = hi * dlnhi0dt9
            dhsdt9 = dfhsdt9(iz41(mu)) + dfhsdt9(iz42(mu)) + dfhsdt9(iz43(mu)) + dfhsdt9(iz44(mu)) - dfhsdt9(iz4c(mu))

            ! Select Screening factor for 4 reactant reactions
            Call screen_blend(lambda,hw,hi,hs,dhwdt9,dhidt9,dhsdt9,h,dhdt9)
            h4(mu,izb) = h
            dh4dt9(mu,izb) = dhdt9
          Else
            h4(mu,izb) = 0.0
            dh4dt9(mu,izb) = 0.0
          EndIf
        EndDo
      EndIf
    EndDo

    If ( idiag >= 5 ) Then
      !__dir_update &
      !__dir_wait(tid) &
      !__dir_host(h1,h2,h3,h4,dh1dt9,dh2dt9,dh3dt9,dh4dt9) &
      !__dir_host(ztilde,zinter,lambda0,gammae,dztildedt9)
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          izone = izb + szbatch - zb_lo
          Write(lun_diag,"(a,i5)") 'SCREEN',izone
          hw0 = ztilde(izb) * lambda0(izb)
          hi0 = kbi * zinter(izb) * lambda0(izb)**bi
          dlnhw0dt9 = dztildedt9(izb)/ztilde(izb) - 1.5/t9t(izb)
          dlnhi0dt9 = - thbim2*dztildedt9(izb)/ztilde(izb) - 1.5*bi/t9t(izb)
          fhs(0) = 0.0
          dfhsdt9(0) = 0.0
          Do j = 1, izmax+2
            gammaz = gammae(izb) * zseq53(j)
            gammaz5 = gammaz**cds(5)
            lgammaz = log(gammaz)
            fhs(j) = cds(1)*gammaz + cds(2)/cds(5)*gammaz5 + cds(3)*lgammaz + cds(4)
            dfhsdt9(j) = -(cds(1)*gammaz + cds(2)*gammaz5 + cds(3))/t9t(izb)
          EndDo
          Do mu = 1, nreac(2)
            lambda = hw0 * zeta2w(mu)
            hw = lambda
            hi = hi0 * zeta2i(mu)
            hs = fhs(iz21(mu)) + fhs(iz22(mu)) - fhs(iz2c(mu))
            Write(lun_diag,"(3a5,i6,5es23.15)") &
              & 'H2',(nname(n2i(j,mu)),j=1,2),mu,lambda,h2(mu,izb),hw,hi,hs
          EndDo
          Do mu = 1, nreac(3)
            lambda = hw0 * zeta3w(mu)
            hw = lambda
            hi = hi0 * zeta3i(mu)
            hs = fhs(iz31(mu)) + fhs(iz32(mu)) + fhs(iz33(mu)) - fhs(iz3c(mu))
            Write(lun_diag,"(4a5,i6,5es23.15)") &
              & 'H3',(nname(n3i(j,mu)),j=1,3),mu,lambda,h3(mu,izb),hw,hi,hs
          EndDo
          Do mu = 1, nreac(4)
            lambda = hw0 * zeta4w(mu)
            hw = lambda
            hi = hi0 * zeta4i(mu)
            hs = fhs(iz41(mu)) + fhs(iz42(mu)) + fhs(iz43(mu)) + fhs(iz44(mu)) - fhs(iz4c(mu))
            Write(lun_diag,"(5a5,i6,5es23.15)") &
              & 'H4',(nname(n4i(j,mu)),j=1,4),mu,lambda,h4(mu,izb),hw,hi,hs
          EndDo
          If ( iheat > 0 ) Then
            Do mu = 1, nreac(2)
              lambda = hw0 * zeta2w(mu)
              hw = lambda
              hi = hi0 * zeta2i(mu)
              hs = fhs(iz21(mu)) + fhs(iz22(mu)) - fhs(iz2c(mu))
              dhwdt9 = hw * dlnhw0dt9
              dhidt9 = hi * dlnhi0dt9
              dhsdt9 = dfhsdt9(iz21(mu)) + dfhsdt9(iz22(mu)) - dfhsdt9(iz2c(mu))
              Write(lun_diag,"(a7,2a5,i6,4es23.15)") &
                & 'dH2/dT9',(nname(n2i(j,mu)),j=1,2),mu,dh2dt9(mu,izb),dhwdt9,dhidt9,dhsdt9
            EndDo
            Do mu = 1, nreac(3)
              lambda = hw0 * zeta3w(mu)
              hw = lambda
              hi = hi0 * zeta3i(mu)
              hs = fhs(iz31(mu)) + fhs(iz32(mu)) + fhs(iz33(mu)) - fhs(iz3c(mu))
              dhwdt9 = hw * dlnhw0dt9
              dhidt9 = hi * dlnhi0dt9
              dhsdt9 = dfhsdt9(iz31(mu)) + dfhsdt9(iz32(mu)) + dfhsdt9(iz33(mu)) - dfhsdt9(iz3c(mu))
              Write(lun_diag,"(a7,3a5,i6,4es23.15)") &
                & 'dH3/dT9',(nname(n3i(j,mu)),j=1,3),mu,dh3dt9(mu,izb),dhwdt9,dhidt9,dhsdt9
            EndDo
            Do mu = 1, nreac(4)
              lambda = hw0 * zeta4w(mu)
              hw = lambda
              hi = hi0 * zeta4i(mu)
              hs = fhs(iz41(mu)) + fhs(iz42(mu)) + fhs(iz43(mu)) + fhs(iz44(mu)) - fhs(iz4c(mu))
              dhwdt9 = hw * dlnhw0dt9
              dhidt9 = hi * dlnhi0dt9
              dhsdt9 = dfhsdt9(iz41(mu)) + dfhsdt9(iz42(mu)) + dfhsdt9(iz43(mu)) + dfhsdt9(iz44(mu)) - dfhsdt9(iz4c(mu))
              Write(lun_diag,"(a7,4a5,i6,4es23.15)") &
                & 'dH4/dT9',(nname(n4i(j,mu)),j=1,4),mu,dh4dt9(mu,izb),dhwdt9,dhidt9,dhsdt9
            EndDo
          EndIf
        EndIf
      EndDo
    EndIf

    !__dir_exit_data &
    !__dir_async(tid) &
    !__dir_delete(mask)

    stop_timer = xnet_wtime()
    timer_scrn = timer_scrn + stop_timer

    Return
  End Subroutine screening

  Subroutine screen_blend(lambda,hw,hi,hs,dhwdt9,dhidt9,dhsdt9,h,dhdt9)
    !-----------------------------------------------------------------------------------------------
    ! This function linearly blends screening prescriptions
    !-----------------------------------------------------------------------------------------------
    !__dir_routine_seq
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: lambda, hw, hi, hs, dhwdt9, dhidt9, dhsdt9

    ! Output variables
    Real(dp), Intent(out) :: h, dhdt9

    ! Local variables
    Real(dp) :: theta

    If ( lambda <= lam_1 ) Then
      h = hw
      dhdt9 = dhwdt9
    ElseIf ( lambda <= lam_2 ) Then
      theta = (lambda - lam_2) * rdlam12
      h = theta*hw + (1.0 - theta)*hi
      dhdt9 = theta*dhwdt9 + (1.0 - theta)*dhidt9
    ElseIf ( lambda <= lam_3 ) Then
      h = hi
      dhdt9 = dhidt9
    ElseIf ( lambda > lam_6 ) Then
      h = hs
      dhdt9 = dhsdt9
    ElseIf ( hi >= hs ) Then
      If ( lambda <= lam_4 ) Then
        theta = (lambda - lam_4) * rdlam34
        h = theta*hi + (1.0 - theta)*hs
        dhdt9 = theta*dhidt9 + (1.0 - theta)*dhsdt9
      Else
        h = hs
        dhdt9 = dhsdt9
      EndIf
    Else
      If ( lambda <= lam_5 ) Then
        h = hi
        dhdt9 = dhidt9
      Else
        theta = (lambda - lam_5) * rdlam56
        h = theta*hs + (1.0 - theta)*hi
        dhdt9 = theta*dhsdt9 + (1.0 - theta)*dhidt9
      EndIf
    EndIf

    Return
  End Subroutine screen_blend

End Module xnet_screening
