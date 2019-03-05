!***************************************************************************************************
! xnet_screening.f90 10/18/17
! This file contains the routines needed to calculate screening corrections for reaction rates.
!***************************************************************************************************

Module xnet_screening
  !-------------------------------------------------------------------------------------------------
  ! This module contains data and routines used to calculate the screening corrections that appear
  ! in the reaction rate exponents.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None
  Real(dp), Allocatable :: h1(:,:), h2(:,:), h3(:,:), h4(:,:)                 ! Screening factors
  Real(dp), Allocatable :: dh1dt9(:,:), dh2dt9(:,:), dh3dt9(:,:), dh4dt9(:,:) ! d(screening factors)/dT9
  !$omp threadprivate(h1,h2,h3,h4,dh1dt9,dh2dt9,dh3dt9,dh4dt9)

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

Contains

  Subroutine screening_init
    !-----------------------------------------------------------------------------------------------
    ! This routine allocates and initializes screening arrays
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: zz
    Use reaction_data, Only: nreac, n2i, n3i, n4i
    Use xnet_constants, Only: bip1, thbim1, five3rd
    Use xnet_controls, Only: iheat, iscrn, nzbatchmx
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
      z21(:) = zz(n2i(1,:))
      z22(:) = zz(n2i(2,:))
      iz21(:) = nint(z21(:))
      iz22(:) = nint(z22(:))
      iz2c(:) = iz21(:) + iz22(:)
      z2c(:) = real(iz2c(:),dp)
      zeta2w(:) = z21(:)*z22(:)
      zeta2i(:) = z2c(:)**bip1 - z21(:)**bip1 - z22(:)**bip1

      ! 3-reactant screening terms
      Allocate (z31(nr3),z32(nr3),z33(nr3))
      Allocate (iz31(nr3),iz32(nr3),iz33(nr3))
      Allocate (iz3c(nr3),z3c(nr3),zeta3w(nr3),zeta3i(nr3))
      z31(:) = zz(n3i(1,:))
      z32(:) = zz(n3i(2,:))
      z33(:) = zz(n3i(3,:))
      iz31(:) = nint(z31(:))
      iz32(:) = nint(z32(:))
      iz33(:) = nint(z33(:))
      iz3c(:) = iz31(:) + iz32(:) + iz33(:)
      z3c(:) = real(iz3c(:),dp)
      zeta3w(:) = z31(:)*z32(:) + z31(:)*z33(:) + z32(:)*z33(:)
      zeta3i(:) = z3c(:)**bip1 - z31(:)**bip1 - z32(:)**bip1 - z33(:)**bip1

      ! 4-reactant screening terms
      Allocate (z41(nr4),z42(nr4),z43(nr4),z44(nr4))
      Allocate (iz41(nr4),iz42(nr4),iz43(nr4),iz44(nr4))
      Allocate (iz4c(nr4),z4c(nr4),zeta4w(nr4),zeta4i(nr4))
      z41(:) = zz(n4i(1,:))
      z42(:) = zz(n4i(2,:))
      z43(:) = zz(n4i(3,:))
      z44(:) = zz(n4i(4,:))
      iz41(:) = nint(z41(:))
      iz42(:) = nint(z42(:))
      iz43(:) = nint(z43(:))
      iz44(:) = nint(z44(:))
      iz4c(:) = iz41(:) + iz42(:) + iz43(:) + iz44(:)
      z4c(:) = real(iz4c(:),dp)
      zeta4w(:) = z41(:)*(z42(:) + z43(:) + z44(:)) + z42(:)*(z43(:) + z44(:)) + z43(:) * z44(:)
      zeta4i(:) = z4c(:)**bip1 - z41(:)**bip1 - z42(:)**bip1 - z43(:)**bip1 - z44(:)**bip1
    EndIf

    !$omp parallel default(shared)
    Allocate (h1(nr1,nzbatchmx))
    Allocate (h2(nr2,nzbatchmx))
    Allocate (h3(nr3,nzbatchmx))
    Allocate (h4(nr4,nzbatchmx))
    h1(:,:) = 0.0
    h2(:,:) = 0.0
    h3(:,:) = 0.0
    h4(:,:) = 0.0
    If ( iheat > 0 ) Then
      Allocate (dh1dt9(nr1,nzbatchmx))
      Allocate (dh2dt9(nr2,nzbatchmx))
      Allocate (dh3dt9(nr3,nzbatchmx))
      Allocate (dh4dt9(nr4,nzbatchmx))
      dh1dt9(:,:) = 0.0
      dh2dt9(:,:) = 0.0
      dh3dt9(:,:) = 0.0
      dh4dt9(:,:) = 0.0
    EndIf
    !$omp end parallel

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
    Use xnet_abundances, Only: yt
    Use xnet_constants, Only: bi, bip1, cds, kbi, thbim2
    Use xnet_conditions, Only: rhot, t9t, etae, detaedt9
    Use xnet_controls, Only: idiag, iheat, iscrn, lun_diag, nzbatchmx, szbatch, lzactive
    Use xnet_eos, Only: eos_screen
    Use xnet_types, Only: dp
    Implicit None

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(:)

    ! Local variables
    Integer :: j, mu, izb, izone
    Real(dp), Dimension(nreac(2)) :: h2w, h2i, h2s, lambda12
    Real(dp), Dimension(nreac(2)) :: dh2wdt9, dh2idt9, dh2sdt9
    Real(dp), Dimension(nreac(3)) :: h3w, h3i, h3s, lambda123
    Real(dp), Dimension(nreac(3)) :: dh3wdt9, dh3idt9, dh3sdt9
    Real(dp), Dimension(nreac(4)) :: h4w, h4i, h4s, lambda1234
    Real(dp), Dimension(nreac(4)) :: dh4wdt9, dh4idt9, dh4sdt9
    Real(dp), Dimension(0:izmax+2) :: gammaz, fhs, fhi, dfhsdt9
    Real(dp) :: ztilde, zinter, lambda0, gammae, dztildedt9
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask => mask_in
    Else
      mask => lzactive
    EndIf
    If ( .not. any(mask(:)) ) Return

    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then

        ! No screening term for 1-reactant reactions
        h1(:,izb) = 0.0
        If ( iheat > 0 ) dh1dt9(:,izb) = 0.0

        ! Call EOS to get plasma quantities
        call eos_screen(t9t(izb),rhot(izb),yt(:,izb),etae(izb),detaedt9(izb), &
          & ztilde,zinter,lambda0,gammae,dztildedt9)

        ! Calculate screening energies as a function of Z, for prescriptions that follow this approach
        gammaz(0) = 0.0
        gammaz(1:izmax+2) = gammae * zseq53(1:izmax+2)
        fhi(0) = 0.0
        fhi(1:izmax+2) = kbi * zinter * lambda0**bi * zseqi(1:izmax+2)
        fhs(0) = 0.0
        fhs(1:izmax+2) = + cds(1) * gammaz(1:izmax+2) &
          &              + cds(2) * gammaz(1:izmax+2)**cds(5) / cds(5) &
          &              + cds(3) * log(gammaz(1:izmax+2)) &
          &              + cds(4)
        dfhsdt9(0) = 0.0
        dfhsdt9(1:izmax+2) = + cds(1) * gammaz(1:izmax+2) &
          &                  + cds(2) * gammaz(1:izmax+2)**cds(5) &
          &                  + cds(3) * log(gammaz(1:izmax+2))
        dfhsdt9(1:izmax+2) = -dfhsdt9(1:izmax+2) / t9t(izb)

        ! Weak and intermediate screening factors, Table 4 of Graboske et al. (1973)
        lambda12 = zeta2w * ztilde * lambda0
        h2w = lambda12
        h2i = kbi * zinter * lambda0**bi * zeta2i
        !h2i = fhi(iz2c) - fhi(iz21) - fhi(iz22)

        ! Strong screening from Dewitt & Slattery (2003) using linear mixing.
        h2s = fhs(iz21) + fhs(iz22) - fhs(iz2c)

        ! Select Screening factor for 2 reactant reactions
        Where ( iz2c == 0 )
          h2(:,izb) = 0.0
        ElseWhere ( lambda12 < 0.1 )
          h2(:,izb) = h2w
        ElseWhere ( lambda12 < 2.0 )
          h2(:,izb) = h2i
        ElseWhere ( lambda12 > 5.0 )
          h2(:,izb) = h2s
        ElseWhere ( h2i < h2s )
          h2(:,izb) = h2i
        ElseWhere
          h2(:,izb) = h2s
        EndWhere
        If ( iheat > 0 ) Then
          dh2wdt9 = +h2w * (dztildedt9/ztilde - 1.5/t9t(izb))
          dh2idt9 = -h2i * (thbim2*dztildedt9/ztilde + bi*1.5/t9t(izb))
          dh2sdt9 = dfhsdt9(iz21) + dfhsdt9(iz22) - dfhsdt9(iz2c)
          Where ( iz2c == 0 )
            dh2dt9(:,izb) = 0.0
          ElseWhere ( lambda12 < 0.1 )
            dh2dt9(:,izb) = dh2wdt9
          ElseWhere ( lambda12 < 2.0 )
            dh2dt9(:,izb) = dh2idt9
          ElseWhere ( lambda12 > 5.0 )
            dh2dt9(:,izb) = dh2sdt9
          ElseWhere ( h2i < h2s )
            dh2dt9(:,izb) = dh2idt9
          ElseWhere
            dh2dt9(:,izb) = dh2sdt9
          EndWhere
        EndIf

        ! Weak and intermediate screening factors, Table 4 of Graboske+ (1973)
        lambda123 = zeta3w * ztilde * lambda0
        h3w = lambda123
        !h3i = kbi * zinter * lambda0**bi * zeta3i
        h3i = fhi(iz3c) - fhi(iz31) - fhi(iz32) - fhi(iz33)

        ! Strong screening from Dewitt & Slattery (2003) using linear mixing.
        h3s = fhs(iz31) + fhs(iz32) + fhs(iz33) - fhs(iz3c)

        ! Select screening factor for 3 reactant reactions
        Where ( iz3c == 0 )
          h3(:,izb) = 0.0
        ElseWhere ( lambda123 < 0.1 )
          h3(:,izb) = h3w
        ElseWhere ( lambda123 < 2.0 )
          h3(:,izb) = h3i
        ElseWhere ( lambda123 > 5.0 )
          h3(:,izb) = h3s
        ElseWhere ( h3i < h3s )
          h3(:,izb) = h3i
        ElseWhere
          h3(:,izb) = h3s
        EndWhere

        If ( iheat > 0 ) Then
          dh3wdt9 = +h3w * (dztildedt9/ztilde - 1.5/t9t(izb))
          dh3idt9 = -h3i * (thbim2*dztildedt9/ztilde + bi*1.5/t9t(izb))
          dh3sdt9 = dfhsdt9(iz31) + dfhsdt9(iz32) + dfhsdt9(iz33) - dfhsdt9(iz3c)
          Where ( iz3c == 0 )
            dh3dt9(:,izb) = 0.0
          ElseWhere ( lambda123 < 0.1 )
            dh3dt9(:,izb) = dh3wdt9
          ElseWhere ( lambda123 < 2.0 )
            dh3dt9(:,izb) = dh3idt9
          ElseWhere ( lambda123 > 5.0 )
            dh3dt9(:,izb) = dh3sdt9
          ElseWhere ( h3i < h3s )
            dh3dt9(:,izb) = dh3idt9
          ElseWhere
            dh3dt9(:,izb) = dh3sdt9
          EndWhere
        EndIf

        ! Weak and intermediate screening factors, Table 4 of Graboske+ (1973)
        lambda1234 = zeta4w * ztilde * lambda0
        h4w = lambda1234
        !h4i = kbi * zinter * lambda0**bi * zeta4i
        h4i = fhi(iz4c) - fhi(iz41) - fhi(iz42) - fhi(iz43) - fhi(iz44)

        ! Strong screening from Dewitt & Slattery (2003) using linear mixing.
        h4s = fhs(iz41) + fhs(iz42) + fhs(iz43) + fhs(iz44) - fhs(iz4c)

        ! Select screening factor for 4 reactant reactions
        Where ( iz4c == 0 )
          h4(:,izb) = 0.0
        ElseWhere ( lambda1234 < 0.1 )
          h4(:,izb) = h4w
        ElseWhere ( lambda1234 < 2.0 )
          h4(:,izb) = h4i
        ElseWhere ( lambda1234 > 5.0 )
          h4(:,izb) = h4s
        ElseWhere ( h4i < h4s )
          h4(:,izb) = h4i
        ElseWhere
          h4(:,izb) = h4s
        EndWhere

        If ( iheat > 0 ) Then
          dh4wdt9 = +h4w * (dztildedt9/ztilde - 1.5/t9t(izb))
          dh4idt9 = -h4i * (thbim2*dztildedt9/ztilde + bi*1.5/t9t(izb))
          dh4sdt9 = dfhsdt9(iz41) + dfhsdt9(iz42) + dfhsdt9(iz43) + dfhsdt9(iz44) - dfhsdt9(iz4c)
          Where ( iz4c == 0 )
            dh4dt9(:,izb) = 0.0
          ElseWhere ( lambda1234 < 0.1 )
            dh4dt9(:,izb) = dh4wdt9
          ElseWhere ( lambda1234 < 2.0 )
            dh4dt9(:,izb) = dh4idt9
          ElseWhere ( lambda1234 > 5.0 )
            dh4dt9(:,izb) = dh4sdt9
          ElseWhere ( h4i < h4s )
            dh4dt9(:,izb) = dh4idt9
          ElseWhere
            dh4dt9(:,izb) = dh4sdt9
          EndWhere
        EndIf

        If ( idiag >= 5 ) Then
          izone = izb + szbatch - 1
          Write(lun_diag,"(a,i5)") 'SCREEN',izone
          Write(lun_diag,"(3a5,i6,5es23.15)") &
            & ('H2',(nname(n2i(j,mu)),j=1,2),mu,lambda12(mu),  h2(mu,izb),h2w(mu),h2i(mu),h2s(mu),mu=1,nreac(2))
          Write(lun_diag,"(4a5,i6,5es23.15)") &
            & ('H3',(nname(n3i(j,mu)),j=1,3),mu,lambda123(mu), h3(mu,izb),h3w(mu),h3i(mu),h3s(mu),mu=1,nreac(3))
          Write(lun_diag,"(5a5,i6,5es23.15)") &
            & ('H4',(nname(n4i(j,mu)),j=1,4),mu,lambda1234(mu),h4(mu,izb),h4w(mu),h4i(mu),h4s(mu),mu=1,nreac(4))
          If ( iheat > 0 ) Then
            Write(lun_diag,"(a7,2a5,i6,4es23.15)") &
              & ('dH2/dT9',(nname(n2i(j,mu)),j=1,2),mu,dh2dt9(mu,izb),dh2wdt9(mu),dh2idt9(mu),dh2sdt9(mu),mu=1,nreac(2))
            Write(lun_diag,"(a7,3a5,i6,4es23.15)") &
              & ('dH3/dT9',(nname(n3i(j,mu)),j=1,3),mu,dh3dt9(mu,izb),dh3wdt9(mu),dh3idt9(mu),dh3sdt9(mu),mu=1,nreac(3))
            Write(lun_diag,"(a7,4a5,i6,4es23.15)") &
              & ('dH4/dT9',(nname(n4i(j,mu)),j=1,4),mu,dh4dt9(mu,izb),dh4wdt9(mu),dh4idt9(mu),dh4sdt9(mu),mu=1,nreac(4))
          EndIf
        EndIf
      EndIf
    EndDo

    Return
  End Subroutine screening

End Module xnet_screening
