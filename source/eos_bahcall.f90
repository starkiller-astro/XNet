!***********************************************************************
! EoS Replacement based on Bahcall () for XNet 6.0 7/6/2010
! This file contains routines which calculate EoS quantites needed to
! calculate screening corrections for reaction rates.
!***********************************************************************

Subroutine eos_initialize
!-----------------------------------------------------------------------
! This routine initializes the EoS
!-----------------------------------------------------------------------
End Subroutine eos_initialize


Subroutine eos_cv(rho,t9,y,cv)
  Use nuclear_data
  Use constants
  Use controls
  Real(8), Intent(in) :: t9,rho,y(ny)
  Real(8), Intent(out) :: cv
  Real(8) :: ye,ytot,bkt,abar,zbar,z2bar,zibar,z53bar,z53(ny)
  Real(8) :: ae,gam,gam14,gam32,gamc2
  Real(8) :: eion,erad,ecoul
  Real(8) :: deiondt9,deraddt9,decouldt9
  Real(8), Parameter :: asig=8.56345d31
  Real(8), Parameter :: onethird=1./3.,twothird=2./3.,fivethird=5./3.
  Real(8), Parameter :: a1=0.898004,b1=0.96786,c1=0.220703,d1=-0.86097
  Real(8), Parameter :: a2=-0.86602540378,b2=0.29561,c2=1.9885

! Calculate Ye and other needed moments of the abundance distribution
  call y_moment(y,ye,ytot,abar,zbar,z2bar,zibar)  

! Calculate energy derivatives (assume etot = eion + erad + ecoul) (ignore degenerate electron energy)
  bkt=bok*t9

  eion=1.5*bkt*ytot
  deiondt9=1.5*ytot*bok

  erad=asig*bkt*bkt*bkt*bkt/(rho*avn)
  deraddt9=4.0*erad/t9

  ae=(3./(4.*pi*avn*rho*ye))**onethird ! electron-sphere radius
  z53(1:ny)=zz(1:ny)**fivethird
  z53bar=sum(z53*y)/ytot
  gam=z53bar*e2/(ae*bkt) ! ionic Coulomb coupling parameter 
  If(gam>=1.0) Then
    gam14=gam**0.25
    ecoul=(a1*gam+b1*gam14+c1/gam14+d1)*bkt*ytot
    decouldt9=ecoul/t9-bok*ytot*(a1*gam+0.25*(b1*gam14-c1/gam14))
  Else
    gam32=gam*sqrt(gam)
    gamc2=gam**c2
    ecoul=(a2*gam32+b2*gamc2)*bkt*ytot
    decouldt9=ecoul/t9-bok*ytot*(1.5*a2*gam32+b2*c2*gamc2)
  EndIf

  cv=deiondt9+deraddt9+decouldt9

  If(idiag>0) Write(lun_diag,'(a,5es12.5)') 'CV',t9,rho,ye,cv

  Return
End Subroutine eos_cv

Subroutine eos_interface(t9,rho,y,ye,ztilde,zinter,lambda0,gammae,dztildedt9)
!-----------------------------------------------------------------------
! This routine Bahcall's approach to calculate the factors needed for 
! screening from the input temperature, density and composition.
!-----------------------------------------------------------------------
  Use constants
  Use controls
  Real(8), Intent(in) :: t9, rho, y(*)
  Real(8), Intent(out) :: ztilde, zinter, lambda0, gammae, ye, dztildedt9
  Real(8) :: ytot,bkt,abar,zbar,z2bar,zibar
  Real(8) :: etae_mb,sratio,efermkt,rel_ef,emass,efc,ae
  Real(8) :: defermktdt9,dsratiodefermkt
  Real(8) :: onethird=1./3.,twothird=2./3.

! Calculate Ye and other needed moments of the abundance distribution
  call y_moment(y,ye,ytot,abar,zbar,z2bar,zibar)  

! Calculate electon distribution
  bkt=bok*t9
  emass=ele_en/clt**2
  etae_mb=log(rho*avn*ye*(2.0*pi*hbar**2/(emass*bkt))**1.5)
  rel_ef=hbar*(3.0*pi**2*rho*avn*ye)**onethird/(emass*clt)
  efermkt=ele_en*(sqrt(1+rel_ef**2)-1)/bkt
  efc=.5*hbar**2*(3.0*pi**2*rho*avn*ye)**twothird/(emass*bkt)
 !Write(lun_diag,'(a4,6es12.5)') 'MUb',bkt,etae_mb,efermkt,efc,rel_ef

! Calculate ratio f'/f for electrons (Salpeter, Eq. 24)
  call salpeter_ratio(efermkt,sratio,dsratiodefermkt)
  ztilde=sqrt(z2bar+zbar*sratio)
  If(iheat>0) Then
    defermktdt9=-efermkt/t9
    dztildedt9=0.5*zbar/ztilde * dsratiodefermkt*defermktdt9
  Else
    dztildedt9=0.0
  EndIf
  
! Calculate plasma quantities
  lambda0=sqrt(4*pi*rho*avn*ytot)*(e2/bkt)**1.5 ! DGC, Eq. 3
  ae=(3./(4.*pi*avn*rho*ye))**onethird ! electron-sphere radius
  gammae=e2/(ae*bkt) ! electron Coulomb coupling parameter 
  zinter=zibar/(ztilde**.58*zbar**.28)
  If(idiag>0) Write(lun_diag,'(a14,9es12.5)') 'Bahcall SCRN', t9, rho, ye, z2bar, zbar, sratio, ztilde,ztilde*lambda0, gammae
  
  Return
End Subroutine eos_interface

Subroutine salpeter_ratio(efmkt,ratio,dratiodefmkt)
!-----------------------------------------------------------------------------
! This routine calculates the Salpeter (1954) ratio f'/f(eta) needed for
! electron screening, using a fit to Figure 24 in that paper.  
! efmkt is the ratio of electron chemical potential to kT.
!-----------------------------------------------------------------------------
  Use controls
  Real(8) :: efmkt,lefmkt,ratio,dratiodefmkt
  lefmkt=log10(efmkt)
  If(lefmkt<=-1.578347) Then ! Bahcall uses limit of -2, but this gives ratio slightly above 1.
    ratio=1.0
  Elseif(lefmkt>=1.5) Then
    ratio=0.0
  Else
    ratio=0.75793-(0.54621*lefmkt)-(0.30964*lefmkt**2)+(0.12535*lefmkt**3) &
&    +(0.1203*lefmkt**4)-(0.012857*lefmkt**5)-(0.014768*lefmkt**6)
  Endif

  If(iheat>0.and.lefmkt>1.578347.and.lefmkt<1.5) Then
    dratiodefmkt=(-0.54621-0.61928*lefmkt+0.37605*lefmkt**2+0.4812*lefmkt**3 &
&    -0.064285*lefmkt**4-0.088608*lefmkt**5)/(log(10.0)*efmkt)
  Else
    dratiodefmkt=0.0
  EndIf

  Return
End Subroutine salpeter_ratio

Subroutine y_moment(y,ye,ytot,abar,zbar,z2bar,zibar)  
!------------------------------------------------------------------------------  
! This routine calculates moments of the abundance distribution for the EOS.
!------------------------------------------------------------------------------  
  Use nuclear_data
  Use controls
  Real(8), Intent(in)  :: y(ny)
  Real(8), Intent(out) :: ye,ytot,abar,zbar,z2bar,zibar
  Real(8)              :: atot,ztot

! Calculate abundance moments
  ytot =sum(y(:))
  atot =sum(aa*y)
  ztot =sum(zz*y)
  abar =atot/ytot
  zbar =ztot/ytot
  z2bar=sum(zz*zz*y)/ytot
  zibar=sum(zz**1.58*y)/ytot
  ye=ztot
  If(idiag>0) Write(lun_diag,'(a4,6es12.5)') 'YMom',ytot,abar,zbar,z2bar,zibar,ye
  
  Return 
End Subroutine y_moment
