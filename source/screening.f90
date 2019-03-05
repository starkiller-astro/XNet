!***********************************************************************
! Screening for XNet 6.0 7/6/2010
! This file contains the routines needed to calculate screening
! corrections for reaction rates.
!***********************************************************************

Subroutine screening()
!-----------------------------------------------------------------------------
! This routine calculates the screening factors necessary for Xnet.
! The HELMHOLTZ Equation of State (Timmes & Swesty 1999) is Used to
! determine the electron distribution and chemical potential.
!
! References:
! Weak Screening:
!    Salpeter (1954) Aust J Phys 7 353.
! Intermediate Screening:
!    DeWitt, Graboske & Cooper (1973) ApJ 181 439.
!    Graboske, DeWitt, Grossman & Cooper (1973) ApJ 181 457.
! Strong Screening:
!    DeWitt & Slattery (1999) Contrib Plasma Phys 39 97.
!-----------------------------------------------------------------------------
  Use controls
  Use constants
  Use abundances
  Use conditions
  Use nuclear_data
  Use cross_sect_data

  Integer :: j,mu
  Real(8), Dimension(nreac(2)) :: h2w,h2i,h2s,lambda12!,gamma12
  Real(8), Dimension(nreac(2)) :: dh2wdt9,dh2idt9,dh2sdt9
  Real(8), Dimension(nreac(3)) :: h3w,h3i,h3s,lambda123!,gamma123
  Real(8), Dimension(nreac(3)) :: dh3wdt9,dh3idt9,dh3sdt9
  Real(8) :: ztilde,zinter,lambda0,gammae,dztildedt9
  Real(8) :: gamma,z
  Real(8) :: onethird=1./3.,twothird=2./3.,fivethird=5./3.,bkt
  Real(8) :: fhs(0:izmax+2),fhi(0:izmax+2),dfhsdt9(0:izmax+2)
  Real(8) :: cds(5)=(/-.899172,0.602249,-.274823,-1.401915,0.3230064/)

! Call EOS to get plasma quantities
  call eos_interface(t9t,rhot,yt,yet,ztilde,zinter,lambda0,gammae,dztildedt9)

!-----------------------------------------------------------------------------
! Calculate screening energies as a function of Z, for prescriptions that
! follow this approach
!-----------------------------------------------------------------------------
  fhi(0)=0.0
  fhs(0)=0.0
  dfhsdt9(0)=0.0
  Do j=1,izmax+2
    z=Real(j)
    fhi(j)=0.38*zinter*lambda0**.86*(z)**1.86
    gamma=(z)**fivethird*gammae
    fhs(j)=cds(1)*gamma+cds(2)*gamma**cds(5)/cds(5)+cds(3)*log(gamma)+cds(4)
    dfhsdt9(j)=-(cds(1)*gamma+cds(2)*gamma**cds(5)+cds(3)*log(gamma))/t9t
  Enddo

! Loop over 1 reactanct reactions to build screening factors.
  h1=0.0
  If(iheat>0) dh1dt9=0.0

  If(iscrn>0) Then

! Weak and intermediate screening factors, Table 4 of Graboske et al.
    lambda12=z12*ztilde*lambda0
    h2w=lambda12
    h2i=0.38*zinter*lambda0**.86*z12i

! Strong screening from Dewitt & Slattery using linear mixing.
    h2s=fhs(iz21)+fhs(iz22)-fhs(iz12)

! Strong Screening from Itoh, non-radial component
!   gamma12=2.0*z21*z22*gammae/(z21**onethird+z22**onethird)
!   h2s=1.25*gamma12

! Select Screening factor for 2 reactant reactions
!   h2=min(h2w,h2i,h2s)
    Where(iz21==0.or.iz22==0)
      h2=0.0
    ElseWhere(lambda12<0.1)
      h2=h2w
    ElseWhere(lambda12<2.0)
      h2=h2i
    ElseWhere(lambda12>5.0)
      h2=h2s
    ElseWhere(h2i<h2s)
      h2=h2i
    ElseWhere
      h2=h2s
    EndWhere
    If(idiag>3) Then
      Do mu=1,nreac(2)
        Write(lun_diag,'(3a5,i6,6es12.5)') 'H2',nname(n2i(1:2,mu)),mu,lambda12(mu),h2(mu),h2w(mu),h2i(mu),h2s(mu)
      EndDo
    EndIf

! Weak and intermediate screening factors, Table 4 of Graboske et al.
    lambda123=z123*ztilde*lambda0
    h3w=lambda123
    h3i=fhi(iz123)-fhi(iz31)-fhi(iz32)-fhi(iz33)

! Strong screening from Dewitt & Slattery using linear mixing.
    h3s=fhs(iz31)+fhs(iz32)+fhs(iz33)-fhs(iz123)

! Strong Screening from Itoh, non-radial component
!   gamma123=2.0*gammae* &
! & ((z31*z32/(z31**onethird+z32**onethird))+(z31*z32*z33/((z31*z32)**onethird+z33**onethird)))
!   h3s=1.25*gamma123

! Select Screening factor for 3 reactant reactions
!   h3=min(h3w,h3i,h3s)
    Where(iz31==0.or.iz32==0.or.iz33==0)
      h3=0.0
    ElseWhere(lambda123<0.1)
      h3=h3w
    ElseWhere(lambda123<2.0)
      h3=h3i
    ElseWhere(lambda123>5.0)
      h3=h3s
    ElseWhere(h3i<h3s)
      h3=h3i
    ElseWhere
      h3=h3s
    EndWhere
    If(idiag>3) Then
      Do mu=1,nreac(3)
        Write(lun_diag,'(4a5,i6,7es12.5)') 'H3',nname(n3i(1:3,mu)),mu,lambda123(mu),h3(mu),h3w(mu),h3i(mu),h3s(mu)
      EndDo
    EndIf

    If(iheat>0) Then
      dh2wdt9=h2w*(dztildedt9/ztilde-1.5/t9t)
      dh2idt9=-h2i*(0.58*dztildedt9/ztilde+0.86*1.5/t9t)
      dh2sdt9=dfhsdt9(iz21)+dfhsdt9(iz22)-dfhsdt9(iz12)
      Where(iz21==0.or.iz22==0)
        dh2dt9=0.0
      ElseWhere(lambda12<0.1)
        dh2dt9=dh2wdt9
      ElseWhere(lambda12<2.0)
        dh2dt9=dh2idt9
      ElseWhere(lambda12>5.0)
        dh2dt9=dh2sdt9
      ElseWhere(h2i<h2s)
        dh2dt9=dh2idt9
      ElseWhere
        dh2dt9=dh2sdt9
      EndWhere
      If(idiag>3) Then
        Do mu=1,nreac(2)
          Write(lun_diag,'(a7,2a5,i6,5es12.5)') 'dH2/dT9',nname(n2i(1:2,mu)),mu,dh2dt9(mu),dh2wdt9(mu),dh2idt9(mu),dh2sdt9(mu)
        EndDo
      EndIf

      dh3wdt9=h3w*(dztildedt9/ztilde-1.5/t9t)
      dh3idt9=-h3i*(0.58*dztildedt9/ztilde+0.86*1.5/t9t)
      dh3sdt9=dfhsdt9(iz31)+dfhsdt9(iz32)+dfhsdt9(iz33)-dfhsdt9(iz123)
      Where(iz31==0.or.iz32==0.or.iz33==0)
        dh3dt9=0.0
      ElseWhere(lambda123<0.1)
        dh3dt9=dh3wdt9
      ElseWhere(lambda123<2.0)
        dh3dt9=dh3idt9
      ElseWhere(lambda123>5.0)
        dh3dt9=dh3sdt9
      ElseWhere(h3i<h3s)
        dh3dt9=dh3idt9
      ElseWhere
        dh3dt9=dh3sdt9
      EndWhere
      If(idiag>3) Then
        Do mu=1,nreac(3)
          Write(lun_diag,'(a7,3a5,i6,6es12.5)') 'dH3/dT9',nname(n3i(1:3,mu)),mu,dh3dt9(mu),dh3wdt9(mu),dh3idt9(mu),dh3sdt9(mu)
        EndDo
      EndIf
    EndIf
  Else
    h2=0.0
    h3=0.0
    If(iheat>0) Then
      dh2dt9=0.0
      dh3dt9=0.0
    EndIf
  Endif

  Return
End Subroutine screening
