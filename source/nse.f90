!***********************************************************************
! NSE for XNet
! These routines calculate the Nuclear Statistical Equilibrium distribution 
! for matter of a given temperature, density and Ye. 
! 
! Y(A,Z) can be expressed as a product of the free proton and neutron 
! abundances, Yp and Yn, and Cnse which contains all of the hydrodynamic and 
! nuclear parameter dependencies.  
!***********************************************************************


Module nse_abundance
!-----------------------------------------------------------------------------
! This module contains the abundances of all species.
!-----------------------------------------------------------------------------
  Use nuc_number
  Real(8), Dimension(:), allocatable :: ynse
!$OMP THREADPRIVATE(ynse)
End Module nse_abundance
  
Module cnse_data
!-----------------------------------------------------------------------------
! This module contains the NSE coefficients, which contain the nuclear 
! physics of each species.
!-----------------------------------------------------------------------------
  Real(8), Dimension(:), allocatable :: cnse
!$OMP THREADPRIVATE(cnse)
End Module cnse_data

Module screen_nse
!-----------------------------------------------------------------------------
! This module contains the corrections to the binding energies due to 
! nuclear screening (AKA Coulomb corrections).
!-----------------------------------------------------------------------------
  Integer :: zmax
  Integer, Dimension(:), allocatable :: intz
  Real(8), Dimension(:), allocatable :: he
!$OMP THREADPRIVATE(zmax,intz,he)
End Module screen_nse

Subroutine nse_initialize
!-----------------------------------------------------------------------------
! This routine allocates and initializes arrays related to NSE abundances.
! Thses tasks are only performed the first time the routine is called.
!-----------------------------------------------------------------------------
  Use nuclear_data
  Use nse_abundance                                       
  Use cnse_data                                       
  Use screen_nse         
! Logical, save :: first_call=.true.
  
! If this is the first call to this routine  
! If(first_call) Then
    Allocate(ynse(ny),cnse(ny),intz(ny))
    zmax=nint(maxval(zz))
    Allocate(he(zmax))
    Where(zz>1.0)
      intz=nint(zz)
    ElseWhere
      intz=1
    EndWhere

!   first_call=.false.                             
! Endif
End Subroutine nse_initialize

Subroutine nse_solv(kmax,iout,t9,rho,ye,yp,yn,kit,ubind,ytst,ztst,atst)
!-----------------------------------------------------------------------------
! Solution occcurs via a two Dimensional Newton-Raphson method of the 
! equations for charge and nucleon number conservation yields new values 
! for Yp and Yn.  The routine uses previous values of Yp and Yn to 
! accelerate convergence.  
!
! The routine also reports successful convergence.  Failure to converge, 
! requires a more finer step in T, rho or Ye.
!-----------------------------------------------------------------------------
  Use controls
  Use nuclear_data
  Use nse_abundance
  Use screen_nse
  Use cnse_data
  Integer :: kmax,iout,kit,k,i,imax(1)
  Real(8), Intent(in) :: t9,rho,ye
  Real(8), Intent(inout) :: yp,yn
  Real(8), Intent(out) :: ubind,ytst,ztst,atst
  Real(8) :: yp0,yn0,ypl,ynl,ypr,ynr,tol,delyp,delyn,ypt,ynt,testyp,testyn
  Real(8) :: f,g,dfdp,dfdn,dgdp,dgdn,det,ddet,testk,bkt,bktr,ronth
  Real(8) :: abar,zbar,benuc
  tol=1.0e-6
  If(iout==2) tol=1.0e-6

! Save initial Yp and Yn
  yp0=yp
  yn0=yn

! We compute logrithmic abundances to avoid overflow when calculating abund
!    ypl=dlog(yp)
!    ynl=dlog(yn)

! Calculate the screening corrections
  If(iscrn>0) Then
    Call nse_screen(t9,rho,ye)
  Else
    he=0.0
  Endif
  
! Calculate the thermodynamic dependendancies for each iterative T and rho
  Call cnsecalc(t9,rho)
  kit=0


! Interively solve the system of equations, f, the expression for charge 
! conservation, and g, the expression for nucleon number conservation.

  Do k=1,kmax
    ypl=log(yp)
    ynl=log(yn)
    If(maxval(cnse+zz*ypl+nn*ynl)>7.0e+2) Then
      imax=maxloc(abs(cnse+zz*ypl+nn*ynl))
      If(idiag>0) Write(lun_diag,'(a,4es15.7,2i4)') 'Exponent too large ',&
&       cnse(imax(1))+zz(imax(1))*ypl+nn(imax(1))*ynl,&
&       cnse(imax(1)),zz(imax(1))*ypl,nn(imax(1))*ynl,iz(imax(1)),in(imax(1))
      kit=kmax+1
      Exit
    EndIf
    Where(cnse+zz*ypl+nn*ynl<-7.0e+2)
      ynse=exp(-7.0e+2)
    ElseWhere
      ynse=exp(cnse+zz*ypl+nn*ynl)
    EndWhere
    If(maxval(ynse)>1.0e+100) Then
      If(idiag>0) Write(lun_diag,'(a,es10.3)') 'Max(y)=',maxval(ynse)
      kit=kmax+1
      Exit
    Endif
    ypr=1.0/yp
    ynr=1.0/yn
    f=sum(zz*ynse)-ye
    g=sum(aa*ynse)-1.0

! df/dp, df/dn, dg/dn, and dg/dp, the derivatives of f and g with respect to yp and yn.
    dfdp=sum(zz**2*ynse)*ypr
    dfdn=sum(zz*nn*ynse)*ynr
    dgdp=sum(aa*zz*ynse)*ypr
    dgdn=sum(aa*nn*ynse)*ynr
    det=dfdp*dgdn-dfdn*dgdp

! Solve the matrix equation for the changes in yp and yn
    If(det/=0.0)Then
      ddet=1.0/det
      delyp=(dgdn*f-dfdn*g)*ddet
      delyn=(dfdp*g-dgdp*f)*ddet
    Else
      If(idiag>0) Write(lun_diag,'(a)') 'Zero Determinant'
      kit=kmax+2
      Exit
    Endif

! Update the abundances, if positive
    ypt=yp-delyp
    ynt=yn-delyn
    If(ypt>0.0.and.ynt>0.0) Then
      yp=ypt
      yn=ynt
    
! Otherwise exit
    Else
      If(idiag>0) Write(lun_diag,'(a,2es10.3)') 'Non-positive Yp, Yn=',ypt,ynt
      kit=kmax+3
      Exit
    Endif
    If(idiag>=3) Then
      Write(lun_diag,'(a,i3,a,es16.8,a,es16.8)') 'Iteration ',k,' Yp=',yp,' Yn=',yn
      Write(lun_diag,'(a,es10.3,a,es10.3)') 'Mass Cons.= ',f,' Charge Cons.=',g
      Write(lun_diag,'(a,4es10.3)') 'Deriv',dfdp,dfdn,dgdp,dgdn
      Write(lun_diag,'(a,es16.8,a,es16.8,a,es10.3)') 'delYp=',delyp,' delYn=',delyn,' Det=',det
    Endif

! The iteritive solution continues until either kmax iterations are run, or the change in solution is less
! than some minimum tolerance, tol.
    testyp=delyp/yp
    testyn=delyn/yn
    testk=sqrt(testyp**2+testyn**2+f**2+g**2)
    kit=k
    If(testk<tol) Exit
  Enddo

! If NSE converges
  If(kit<kmax) Then
    ypl=log(yp)
    ynl=log(yn)
    Where(cnse+zz*ypl+nn*ynl<-7.0e+2)
      ynse=exp(-7.0e+2)
    ElseWhere
      ynse=exp(cnse+zz*ypl+nn*ynl)
    EndWhere
    If(idiag>=2) Then
      Write(lun_diag,'(a,i3,a,es16.8,a,es16.8)') 'Converged in ',kit,' steps, yp=',yp,', yn=',yn
      If(idiag>=5) Then
        Write(lun_diag,'(3(a4,es10.3))') 'T9',t9,'Rho',rho,'Ye',ye
        Write(lun_diag,'(5(a6,es10.3))') (nname(i),(aa(i)*ynse(i)),i=1,ny) 
      Endif
    Endif

! atst, ztst, ytst are the total mass and electron fraction, 
! abar, zbar are the average baryon number and proton number,
! benuc, ubind are the binding energy, in Mev/nuc and erg/g.
    atst=sum(aa*ynse)
    ztst=sum(zz*ynse)
    ytst=sum(ynse)
    abar=atst/ytst
    zbar=ztst/ytst
    benuc=sum(be*ynse)
    ubind=-9.616221376e17*benuc
!   Where(ynse<1e-99) ynse=0.0
  Else

! NSE fails to converge in Kmax steps
    If(idiag>=3) Write(lun_diag,'(a,i3,a,es10.3,a,es10.3)') &
&  'NSE Does not converge in',kit,' steps at T9=',t9,' and density=',rho 
    yp=yp0
    yn=yn0
  Endif

! Having completed the loop for T,rho (and set iflag=1) or discovered that 
! the solution will not converge (and set iflag=0), theSubroutine returns
! for the next set of T, rho and Ye.  
  Return
End Subroutine nse_solv

Subroutine cnsecalc(t9,rho)      
  Use constants
  Use controls
  Use nuclear_data 
  Use part_funct_data
  Use cnse_data
  Use screen_nse
  Real(8), Intent(in)  :: t9,rho
  Integer :: i,ii
  Real(8) :: rdt9,bkt,bktr,ronth
!--------------------------------------------------------------------------
! The first step in calculating Cnse is to interpolate the partition 
! function, g, for the temperature.           
!-------------------------------------------------------------------------- 
  Do i=1,24
    If(t9<=t9i(i)) Exit
  Enddo
  ii=i
  Select Case (ii)
    Case(2:24)
      rdt9=(t9-t9i(ii-1))/(t9i(ii)-t9i(ii-1))
      gg(1:ny)=(rdt9*log(g(ii,:))+(1-rdt9)*log(g(ii-1,:)))
    Case (1)
      gg(1:ny)=log(g(1,:))
    Case (25)
      gg(1:ny)=log(g(24,:))
  End Select
 If(idiag>=5) Then
    Write(lun_diag,"(a5,i3,es14.7)") 'PartF',ii,t9
    Write(lun_diag,"(4(i4,es16.8))") (i, gg(i), i=1,ny)
  Endif
  bkt=bok*t9 
  bktr=1.0/bkt
  ronth=log(avn*rho)+1.5*log((2.0*pi*hbar**2)/(bkt*amu))
  cnse(1)=0.0
  cnse(2)=0.0
  cnse(3:ny)=log(angm(3:ny)*aa(3:ny)*sqrt(aa(3:ny)))+gg(3:ny) &
&    -.69314720*aa(3:ny)+(aa(3:ny)-1.0)*ronth+be(3:ny)*bktr+he(intz(3:ny))
  If(idiag>=4) Then
    Write(lun_diag,'(a5,2es16.8)') 'Cnse',bktr,ronth
    Write(lun_diag,'(a5,4es16.8)') (nname(i),gg(i),be(i),he(intz(i)),cnse(i),i=1,ny)
  Endif
  Return
End Subroutine cnsecalc  
!
Subroutine nse_screen(t9,rho,ye)
  Use controls
  Use nse_abundance
  Use screen_nse
  Use nuclear_data
  Real(8), Intent(in)  :: t9,rho,ye
  Integer :: j,mu
  Integer, PARAMETER :: iz1=1
  Real(8), PARAMETER :: z1=1.0
  Real(8), Dimension(zmax-1) :: z2,h0,hw,hi,hs,lambda12!,hi2,gamma12
  Integer, Dimension(zmax-1) :: iz2
  Real(8) :: ztot,ztilde,zinter,lambda0,gammae,dztildedt9
  Real(8) :: gamma,z
  Real(8) :: onethird=1./3.,twothird=2./3.,fivethird=5./3.,bkt
  Real(8) :: fhs(0:zmax+1),fhi(0:zmax+1)
  Real(8) :: cds(5)=(/-.899172,0.602249,-.274823,-1.401915,0.3230064/)

! Call EOS to get plasma quantities
  call eos_interface(t9,rho,ynse,ztot,ztilde,zinter,lambda0,gammae,dztildedt9)

!-----------------------------------------------------------------------------
! Calculate screening energies as a function of Z, for prescriptions that 
! follow this approach 
!-----------------------------------------------------------------------------
  fhi(0)=0.0
  fhs(0)=0.0
  Do j=1,zmax+1
    z=Real(j)
    fhi(j)=0.38*zinter*lambda0**.86*(z)**1.86
    gamma=(z)**fivethird*gammae                            
    fhs(j)=cds(1)*gamma+cds(2)*gamma**cds(5)/cds(5)+cds(3)*log(gamma)+cds(4)
  Enddo

! Loop over proton capture reactions to build screening factors.
  h0=0.0

  If(iscrn>0) Then

    Do j=1,zmax-1
      iz2(j)=j
      z2(j)=real(j)
    EndDo

! Weak and intermediate screening factors, Table 4 of Graboske et al.         
    lambda12=z1*z2*ztilde*lambda0    
    hw=lambda12
!   hw=lambda12*(1.+lambda12*(log(lambda12)+0.8364))
    hi=0.38*zinter*lambda0**.86*((z1+z2)**1.86-z1**1.86-z2**1.86)
!   hi2=fhi(iz1+iz2)-fhi(iz1)-fhi(iz2)                                                     

! Strong screening from Dewitt & Slattery using linear mixing.
    hs=fhs(iz1)+fhs(iz2)-fhs(iz1+iz2)

! Strong Screening from Itoh, non-radial component
!   gamma12=2.0*z1*z2*gammae/(z1**onethird+z2**onethird)
!   hs=1.25*gamma12

! Select Screening factor
!   h0=min(hw,hi,hs)
    Where(iz2==0)
      h0=0.0
    ElseWhere(lambda12<0.1)
      h0=hw
    ElseWhere(lambda12<2.0)
      h0=hi
    ElseWhere(lambda12>5.0)
      h0=hs
    ElseWhere
      h0=min(hi,hs)
    EndWhere
    If(idiag>3) Then
      Do mu=2,zmax
        Write(lun_diag,'(a5,3i6,6es12.5)') 'HNSE',iz1,iz2(mu-1),mu,lambda12(mu-1),h0(mu-1),hw(mu-1),hi(mu-1),hs(mu-1)
      EndDo
    EndIf

! Add succeeding screening factors
    he(1)=0.0
    Do mu=2,zmax
      he(mu)=sum(h0(1:(mu-1)))
    EndDo
  Else
    he=0.0
  EndIf

  Return
End Subroutine nse_screen


Subroutine nse_descend(rho,ye,t9fnl,t9out)
!-----------------------------------------------------------------------------
! This routine progress to the desired conditions from an initially high 
! temperature, at which the initial NSE abundances can be accurately predicted.
!-----------------------------------------------------------------------------
  Use controls
  Use nse_abundance
  Integer, parameter :: kmax=10
  Real(8), Intent(in) :: rho,ye,t9fnl,t9out(*)
  Integer :: kit,iout,iw
  Real(8) :: t9start,t9,t9t,delt9,yp,yn,ubind,ztst,atst,ytst,zbar,abar,y_save(ny)
  
! Start with large T9 because for large T9, free nucleons dominate
  t9start=100.0
  yp=ye
  ynse(2)=yp
  yn=1.0-ye
  ynse(1)=yn
  t9=t9start
  delt9=sign(10.d0,(t9fnl-t9start))
  iw=1
  iout=2

! Descend in T9 
  Do  
!   If(yp<1.0e-30.or.yn<1.0e-30) Then
!     Write(6,'(a,2es11.3)') 'Yp,Yn too small',yp,yn
!     t9final=t9
!   Endif

! Solve the NSE
    If(idiag>=2) Write(lun_diag,'(a,4es13.6)') 'NSE trying',t9,delt9,yp,yn
    Call nse_solv(kmax,iout,t9,rho,ye,yp,yn,kit,ubind,ytst,ztst,atst)
    Select Case (kit) 

! If Convergence successful, continue with same step
      Case (3:kmax-1)
        zbar=ztst/ytst
        abar=atst/ytst
        If(idiag>=2) Write(lun_diag,'(a,es12.5,5es13.6,i3)') 'NSE success',t9,yn,yp,zbar,abar,-ubind,kit

! If convergence easy, double step
      Case (1:2)
        delt9=2.0*delt9
        zbar=ztst/ytst
        abar=atst/ytst
        If(idiag>=2) Write(lun_diag,'(a,es12.5,5es13.6,i3)') 'NSE easy',t9,yn,yp,zbar,abar,-ubind,kit

! If convergence fails, retry with delta T9 halved
      Case Default 
        t9=t9-delt9
        delt9=.5*delt9
        ynse = y_save ! Restore abundances so EOS calculation functions
        If(idiag>=2) Write(lun_diag,'(a,3es12.5,i3)') 'NSE Failed, revert to', t9,yn,yp,kit
    End Select
    delt9=sign(min(abs(t9-t9out(iw)),10.0,abs(delt9)),delt9)
    t9t=t9+delt9

! If Final T9, step back to T9Fnl
    If(t9<=t9fnl) Then
      t9=t9fnl
      iout=2
      Exit

! If trial t9 below Output T9, step back and solve for output T9
!   Elseif(t9t<=1.001*t9out(iw)) Then
!     t9=t9out(iw)
!     iw=iw+1
!     iout=2

! Else continue descent
    Else 
      t9=t9+delt9
      iout=1
    Endif
    y_save= ynse
  Enddo

! Solve for T9fnl
  Call nse_solv(kmax,iout,t9,rho,ye,yp,yn,kit,ubind,ytst,ztst,atst)
  zbar=ztst/ytst
  abar=atst/ytst
  Write(lun_diag,'(es12.5,5es13.6,i3)') t9,yn,yp,zbar,abar,-ubind,kit
  Return
End Subroutine nse_descend
  
