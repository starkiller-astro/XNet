!***********************************************************************
!  NSE for XNet
!  These routines calculate the Nuclear Statistical Equilibrium distribution 
!  for matter of a given temperature, density and Ye. 
! 
!  Y(A,Z) can be expressed as a product of the free proton and neutron 
!  abundances, Yp and Yn, and Cnse which contains all of the hydrodynamic and 
!  nuclear parameter dependencies.  
!***********************************************************************


      module nse_abundance
!-----------------------------------------------------------------------------
!  This module contains the abundances of all species.
!-----------------------------------------------------------------------------
      use nuc_number
      real*8, dimension(:), allocatable :: ynse
      end module nse_abundance
      
      module cnse_data
!-----------------------------------------------------------------------------
!  This module contains the NSE coefficients, which contain the nuclear 
!  physics of each species.
!-----------------------------------------------------------------------------
      real*8, dimension(:), allocatable :: cnse
      end module cnse_data

      module screen_nse
!-----------------------------------------------------------------------------
!  This module contains the corrections to the binding energies due to 
!  nuclear screening (AKA Coulomb corrections).
!-----------------------------------------------------------------------------
      integer :: zmax
      integer, dimension(:), allocatable :: intz
      real*8, dimension(:), allocatable :: he
      end module screen_nse

      subroutine nse_solv (kmax,iout,t9,rho,ye,yp,yn,kit,ubind,
     &                   ytst,ztst,atst,iscrn)
!-----------------------------------------------------------------------------
!  Solution occcurs via a two dimensional Newton-Raphson method of the 
!  equations for charge and nucleon number conservation yields new values 
!  for Yp and Yn.  The subroutine uses previous values of Yp and Yn to 
!  accelerate convergence.  
!
!  The subroutine also reports successful convergence.  Failure to converge, 
!  requires a more finer step in T, rho or Ye.
!-----------------------------------------------------------------------------
      use nuclear_data
      use nse_abundance
      use screen_nse
      use cnse_data
      integer :: kmax,iout,kit,iscrn
      real(8) :: t9,rho,ye,yp,yn,ubind,ztst,atst
      tol=1.0e-6
      If(iout==2) tol=1.0e-6

!  Save initial Yp and Yn
      yp0=yp
      yn0=yn

!  We compute logrithmic abundances to avoid overflow when calculating abund
!     ypl=dlog(yp)
!     ynl=dlog(yn)

!  Calculate the screening corrections
      If(iscrn>0) Then
        If(T9>1.0) Then
          call nse_screen(T9,rho,ye,iout)
        Else
          call nse_screen2(T9,rho,ye,iout)
        Endif
      Else
        he=0.0
      Endif
  
!  Calculate the thermodynamic dependendancies for each iterative T and rho
      call cnsecalc(t9,rho,iout)
      kit=0


!  Interively solve the system of equations, f, the expression for charge 
!  conservation, and g, the expression for nucleon number conservation.

      Do k=1,kmax
        ypl=log(yp)
        ynl=log(yn)
        ynse=exp(cnse+zz*ypl+nn*ynl)
        If(maxval(ynse)>1.0e+100) Then
          Write(15,'(a,es10.3)') 'Max(y)=',maxval(ynse)
          kit=kmax+1
          Exit
        Endif
        ypr=1.0/yp
        ynr=1.0/yn
        f=sum(zz*ynse)-ye
        g=sum(aa*ynse)-1.0

!  df/dp, df/dn, dg/dn, and dg/dp, the derivatives of f and g with respect to yp and yn.
        dfdp=sum(zz**2*ynse)*ypr
        dfdn=sum(zz*nn*ynse)*ynr
        dgdp=sum(aa*zz*ynse)*ypr
        dgdn=sum(aa*nn*ynse)*ynr
        det=dfdp*dgdn-dfdn*dgdp

!  Solve the matrix equation for the changes in yp and yn
        If(det/=0.0)Then
          ddet=1.0/det
          delyp=(dgdn*f-dfdn*g)*ddet
          delyn=(dfdp*g-dgdp*f)*ddet
        Else
          Write(15,'(a)') 'Zero Determinant'
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
          Write(15,'(a,2es10.3)') 'Yp, Yn=',ypt,ynt
          kit=kmax+3
          Exit
        Endif
        If(iout>0) Then
          Write(15,'(a,i3,a,es16.8,a,es16.8)')
     &         'Iteration ',k,' Yp=',yp,' Yn=',yn
          Write(15,'(a,es10.3,a,es10.3)')
     &         'Mass Cons.= ',f,' Charge Cons.=',g
          Write(15,'(a,4es10.3)') 'Deriv',dfdp,dfdn,dgdp,dgdn
          Write(15,'(a,es16.8,a,es16.8,a,es10.3)')
     &          'delYp=',delyp,' delYn=',delyn,' Det=',det
        Endif

! The iteritive solution continues until either kmax iterations are run, or the change in solution is less
!  than some minimum tolerance, tol.
        testyp=delyp/yp
        testyn=delyn/yn
        testk=sqrt(testyp**2+testyn**2+f**2+g**2)
        kit=k
        If(testk<tol) Exit
      Enddo

!  If NSE converges
      If(kit<kmax) Then
        ypl=log(yp)
        ynl=log(yn)
        ynse=exp(cnse+ypl*zz+ynl*nn)
        If(iout>0) Then
          Write(15,'(a,i3,a,es16.8,a,es16.8)') 
     &      'Converged in ',kit,' steps, yp=',yp,', yn=',yn

!  atst, ztst, ytst are the total mass and electron fraction, 
!  abar, zbar are the average baryon number and proton number,
!  benuc, ubind are the binding energy, in Mev/nuc and erg/g.
          atst=sum(aa*ynse)
          ztst=sum(zz*ynse)
          ytst=sum(ynse)
          abar=atst/ytst
          zbar=ztst/ytst
          benuc=sum(be*ynse)
          ubind=-9.616221376e17*benuc
          If(iout>1) Then
            Write(14,'(3(a4,es10.3))') 'T9',t9,'Rho',rho,'Ye',ye
!           Where(ynse<1e-99) ynse=0.0
            Write(14,'(5(a6,es10.3))')
     &        (nname(i),(aa(i)*ynse(i)),i=1,ny) 
          Endif
        Endif
      Else

!  NSE fails to converge in Kmax steps
        Write(15,'(a,i3,a,es10.3,a,es10.3)')
     &    'Does not converge in',kit,' steps at T9=',t9,
     &    ' and density=',rho 
        yp=yp0
        yn=yn0
      Endif

!  Having completed the loop for T,rho (and set iflag=1) or discovered that 
!  the solution will not converge (and set iflag=0), the subroutine returns
!  for the next set of T, rho and Ye.  
      Return
      End

      subroutine cnsecalc (t9,rho,iout)              
      use constants
      use nuclear_data 
      use part_funct_data
      use cnse_data
      use screen_nse
      integer :: i,ii
      real(8)  :: t9,rdt9,rho
!--------------------------------------------------------------------------
!  The first step in calculating Cnse is to interpolate the partition 
!  function, g, for the temperature.                           
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
!     If(iout>3) Then
!         Write(15,"(a5,i3,es14.7)") 'PartF',ii,t9
!         Write(15,"(4(i4,es16.8))") (i, gg(i), i=1,ny)
!     Endif
      bkt=bok*t9 
      bktr=1.0/bkt
      ronth=log(avn*rho)+1.5*log((2.0*pi*hbar**2)/(bkt*amu))
      cnse(1)=0.0
      cnse(2)=0.0
      cnse(3:ny)=log(angm(3:ny)*aa(3:ny)*sqrt(aa(3:ny)))+gg(3:ny)
     &          -.69314720*aa(3:ny)+(aa(3:ny)-1.0)*ronth+be(3:ny)*bktr
     &                +he(intz(3:ny))
!      If(iout>3) Then
          Write(15,'(a5,2es16.8)') 'Cnse',bktr,ronth
          Write(15,'(a5,4es16.8)') 
     &         (nname(i),gg(i),be(i),he(intz(i)),cnse(i),i=1,ny)
!      Endif
      RETURN
      END  
!                    
      subroutine nse_screen(T9,rho,ye,iout)
      use screen_nse
c
      third=1.0/3.0
      gnp=2.27472e-4*(rho*ye)**third/T9
      If(iout>3) Write(15,'(a5,5es13.6)') 'SCR',t9,rho,ye,gnp
      he(1)=0.0
      Do 100 j=2,zmax
        z1=1.0
        z2=float(j-1)
        a1=1.0
        a2=2.*z2
        If(z2.eq.1.0) a2=1.0                         
        zz=z1*z2
        z13=z1**third
        z23=z2**third

!  STRONG SCREENING BY ITOH ET AL.(1990)
        gnz=gnp*zz*2.0/(z13+z23)                                           
        em2=a1*a2*2.0/(a1+a2)                                              
        tau=3.37220*(em2*zz*zz/t9)**third
        gt=gnz/tau                                                        
        gt3=3.0*gt
        h0=1.25*gnz

!  Add succeeding screening factors
        he(j)=he(j-1)+h0
         If(iout>3) Write(15,'(a5,5es13.6)') 'H',he(j),h0,gnz,tau
  100 Continue                        
      RETURN 
      END 

      subroutine nse_screen2(t9,rho,ye,iout)
      use nse_abundance
      use screen_nse
      common /scrn/ gnp,gsi,zl,dzldr,dzldt,t9er                         
      he(1)=0.0
      v=1.0/rho
      write(301,*) t9, y
      call state(t9,v,y,pa,ea,dlma,fey,pta,pva,eta,eva,ieos,ny,
     &          emas,ysp,beta)
      write(301,*) t9, y
      Do j=2,zmax
        z1=1.0
        z2=float(j-1)
        a1=1.0
        If(z2>1.0) Then
            a2=2.*z2
        Else
            a2=1.0                         
        Endif

!  Weak and intermediate screening factors, from  Graboske et al. (1973)        
        zz=z1*z2    
        z=zl*zz                                                           
        gsiz=gsi*zz                                                       
        gsb=gsiz**0.86                                                    
        If (gsiz>1.) Then 
            f=0.38*((1.+gsiz)**1.86-gsb*gsiz-1.)/gsb          
        Else
            f=0.61943*gsiz**0.07                              
        Endif
        hw=z*(1.+z*(log(z)+0.8364))                                      
        hi=f*z**0.86                                                      
        If (hw>hi) Then
            h0=hi                                                         
        Else                                                              
            h0=hw                                                         
        Endif

!  Strong screening from Itoh et al.(1990)                             
        z13=z1**(1./3.)                                                   
        z23=z2**(1./3.)                                                   
        gnz=gnp*zz*2./(z13+z23)                                           
        If (gnz>=0.4) Then
            em2=a1*a2*2./(a1+a2)                                              
            tau=3.3722*(em2*zz**2*t9er)**(1./3.)                              
            gt=gnz/tau                                                        
            gt3=3.*gt                                                         
            numf=.0455*gt3+.348*gt3**3.+9.49*gt3**6.-.123*gt3**12.+
     #        .101*gt3**13.   
            denf=1.+100.*gt3**4.+.267*gt3**12.
            f90=numf/denf
            scr=1.25*gnz-tau*f90 
            If (scr<h0) Then
                h0=scr   
            Endif
        Endif

!  Add succeeding screening factors
        he(j)=he(j-1)+h0
         If(iout>3) Write(15,'(a,4es13.6)') 'H2',he(j),h0,gnz,tau
      Enddo                                                             
      Return                                                            
      End                                                               


      subroutine nse_descend(Rho,Ye,t9fnl,t9out,is)
!-----------------------------------------------------------------------------
!  This routine progress to the desired conditions from an initially high 
!  temperature, at which the initial NSE abundances can be accurately predicted.
!-----------------------------------------------------------------------------
      integer, parameter :: kmax=10
      real(8) :: rho,ye,t9fnl
      real(8), dimension(*) :: t9out      
      integer :: kit,iout,iw,is
      real(8) :: t9start,t9,delt9,yp,yn,ubind,ztst,atst,zbar,abar
      
!  Start with large T9 because for large T9, free nucleons dominate
      T9start=100.0
      yp=ye
      yn=1.0-ye
      T9=T9start
      DelT9=sign(10.d0,(t9fnl-t9start))
      iw=2
      iout=2

!  Descend in T9 
      Do  
        Write(6,'(4es13.6)') t9,delt9,yp,yn
!       If(yp<1.0e-30.or.yn<1.0e-30) Then
!          Write(6,'(a,2es11.3)') 'Yp,Yn too small',yp,yn
!          t9final=t9
!       Endif

!  Solve the NSE
        If(t9==5.0) iout=4
        call nse_solv(kmax,iout,t9,rho,ye,yp,yn,kit,ubind,ytst,ztst,
     &                atst,is)
        Select Case (kit) 

!  If Convergence successful, continue with same step
          Case (4:kmax-1)
              zbar=ztst/ytst
              abar=atst/ytst
              Write(12,'(es12.5,5es13.6,i3)') 
     &              t9,yn,yp,zbar,abar,-ubind,kit

!  If convergence easy, double step
          Case (1:3)
              delT9=2.0*delT9
              zbar=ztst/ytst
              abar=atst/ytst
              Write(12,'(es12.5,5es13.6,i3)') 
     &              t9,yn,yp,zbar,abar,-ubind,kit

!  If convergence fails, retry with delta T9 halved
          Case Default 
              T9=T9-delT9
              delT9=.5*delT9
        End Select
        delT9=sign(min(abs(t9-t9out(iw)),abs(delT9)),delT9)
        t9t=t9+delt9

!  If Final T9, step back to T9Fnl
        If(t9<=t9fnl) Then
            t9=t9fnl
            iout=2
            Exit

!  If trial t9 below Output T9, step back and solve for output T9
        Elseif(t9t<=1.001*t9out(iw)) Then
            t9=t9out(iw)
            iw=iw+1
            iout=2

!  Else continue descent
        Else 
            T9=T9+delT9
            iout=1
        Endif
      Enddo

!  Solve for T9fnl
      call nse_solv(kmax,iout,t9,rho,ye,yp,yn,kit,ubind,ytst,ztst,
     &                atst,is)
      zbar=ztst/ytst
      abar=atst/ytst
      Write(12,'(es12.5,5es13.6,i3)') t9,yn,yp,zbar,abar,-ubind,kit
      Return
      end subroutine nse_descend
      
