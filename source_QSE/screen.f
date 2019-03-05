!*******************************************************************************
!  screen.f               
!
!  
!
!*******************************************************************************
!      module screen_qse
!-----------------------------------------------------------------------------
!  This module contains the corrections to the binding energies due to
!  nuclear screening (AKA Coulomb corrections).
!-----------------------------------------------------------------------------
!      integer :: zmax
!      integer, dimension(:), allocatable :: intz
!      real*8, dimension(:), allocatable :: he
!      end module screen_qse


      subroutine en(n,y,ye,yps,y2e)  
!===============================================================================
!  This routine calculates moments of the abundance distribution needed 
!  for the EOS.
!===============================================================================
      use nuclear_data
      integer :: n
      real(8) :: ye,yps,y2e
      real(8), dimension(n) :: y
      ye=sum(zz*y)
      y2e=sum(zz**2*y)
      yps=sum(y)+ye 
!     Write(lun_diag,"(a3,3es16.8)") 'EN',ye,y2e,yps
      Return                                                                    
      End subroutine en                                                                      

      subroutine screen2
!===============================================================================
! This routine calculates screenign factors for reactions with two reactants.
!===============================================================================

      use controls
      use nuclear_data
      use cross_sect_data
      integer :: mu
      real(8) :: z1,z2,z12,a1,a2,z,zl,dzldt,dzldr
      real(8) :: gsi,gsiz,gsb,f,fsr,fst,gnp,gnz,em2,tau,t9er,gt3,f90
      real(8) :: hw,hi,hs,dhdz
      common /scrn/ gnp,gsi,zl,dzldr,dzldt,t9er                         
!     Write(20,"(a3,3es15.7)") 'EOS',GNP,GSI,ZL
      Do mu=1,nreac(2)                                                      
        z1=zz(n2i(1,mu))                                                  
        z2=zz(n2i(2,mu))                                                  
        a1=aa(n2i(1,mu))                                                  
        a2=aa(n2i(2,mu))                                                  
        z12=z1*z2
        If(z12==0.0.or.iscrn==0) Then
          h2(mu)=0.0                                                        
        Else
  
!  Weak and Intermediate Screening Graboske et al.         
          z=zl*z12                                                           
          gsiz=gsi*z12                                                       
          gsb=gsiz**0.86                                                    
          If (gsiz>1.0) Then
            f=0.38*((1.+gsiz)**1.86-gsb*gsiz-1.)/gsb          
          Else
            f=0.61943*gsiz**0.07                              
          Endif
          hw=z*(1.+z*(log(z)+0.8364))                                      
          hi=f*z**0.86                                                      
          If (hw<=hi) Then
            h2(mu)=hw         
            dhdz=z*(1.+2.*z*(log(z)+1.3364))                
          Else
            h2(mu)=hi                             
            dhdz=0.86*hi                                   
          Endif
          fsr=dhdz*dzldr                                                    
          fst=dhdz*dzldt                                                    

!  Strong screening by Itoh et al.(1990)                             
          gnz=gnp*z12*2.0/(z1**(1./3.)+z2**(1./3.))
!          If (gnz>=0.4) Then
            em2=a1*a2*2.0/(a1+a2)                
            tau=3.3722*(em2*z12**2*t9er)**(1./3.) 
            gt3=3.0*gnz/tau
             f90 =0.0
!            f90=(.0455*gt3+.348*gt3**3+9.49*gt3**6-.123*gt3**12+
!     &          .101*gt3**13)/(1.+100.*gt3**4+.267*gt3**12)
            hs=1.25*gnz-tau*f90 
!           hs=gnz*(1.25-0.285*gt3)    
            If (hs<h2(mu)) Then
                h2(mu)=hs               
                fsr=gnz*(1.25-.57*gt3)/3. 
                fst=-gnz*(1.25-.475*gt3)   
            Endif
!          Endif
!         Write(20,"(3a5,i6,4es11.4)") 'SCR2',nname(n2i(1,mu)),
!    &         nname(n2i(2,mu)),dexp(h0(mu)),hw,hi,hs
        Endif
      Enddo                                                             
!     Write(20,*) (nname(n2i(1,mu)),nname(n2i(2,mu)),h2(mu),mu=1,n)
      Return                                                            
      End subroutine screen2   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                           
      subroutine qse_screen(T9,rho,ye,iout)
      use screen_qse
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine qse_screen2(t9,rho,ye,iout)
      use abundances
      use screen_qse
      common /scrn/ gnp,gsi,zl,dzldr,dzldt,t9er
      he(1)=0.0
      v=1.0/rhot
      write(301,*) t9, y
      call state(t9t,v,y,pa,ea,dlma,fey,pta,pva,eta,eva,ieos,ny,
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

      subroutine qse_screen3
!===============================================================================
! This routine calculates screenign factors for reactions with two reactants.
!===============================================================================

      use controls
      use nuclear_data
      use cross_sect_data
      use screen_qse
      use abundances
      integer :: mu
      real(8) :: z1,z2,z12,a1,a2,z,zl,dzldt,dzldr
      real(8) :: gsi,gsiz,gsb,f,fsr,fst,gnp,gnz,em2,tau,t9er,gt3,f90
      real(8) :: hw,hi,hs,dhdz
      common /scrn/ gnp,gsi,zl,dzldr,dzldt,t9er                         
       Write(*,"(a3,3es15.7)") 'EOS',GNP,GSI,ZL
      he(1)=0.0
!      v=1.0/rhot
!        write(*,*) "stuck here"
!       Call state(t9t,v,y,pa,ea,dlma,ye,pta,pva,eta,eva,
!     & ieos,ny,emas,yps,beta)
!      call state(t9,v,y,pa,ea,dlma,fey,pta,pva,eta,eva,ieos,ny,
 !    &          emas,ysp,beta)
      Do j=2,zmax
         write(*,*) j
        z1=1.0                                                
        z2=float(j-1)                                                 
        a1=1.0                                                  
        If(z2>1.0) Then
            a2=2.*z2
        Else
            a2=1.0
        Endif
                                              
        z12=z1*z2
        If(z12==0.0.or.iscrn==0) Then
          h0=0.0                                                        
        Else
  
!  Weak and Intermediate Screening Graboske et al.         
          z=zl*z12                                                           
          gsiz=gsi*z12                                                       
          gsb=gsiz**0.86                                                    
          If (gsiz>1.0) Then
            f=0.38*((1.+gsiz)**1.86-gsb*gsiz-1.)/gsb          
          Else
            f=0.61943*gsiz**0.07                              
          Endif
          hw=z*(1.+z*(log(z)+0.8364))                                      
          hi=f*z**0.86                                                      
          If (hw<=hi) Then
            h0=hw         
            dhdz=z*(1.+2.*z*(log(z)+1.3364))                
          Else
            h0=hi                             
            dhdz=0.86*hi                                   
          Endif
          fsr=dhdz*dzldr                                                    
          fst=dhdz*dzldt                                                    

!  Strong screening by Itoh et al.(1990)                             
          gnz=gnp*z12*2.0/(z1**(1./3.)+z2**(1./3.))
          If (gnz>=0.4) Then
            em2=a1*a2*2.0/(a1+a2)                
            tau=3.3722*(em2*z12**2*t9er)**(1./3.) 
            gt3=3.0*gnz/tau
            f90=(.0455*gt3+.348*gt3**3+9.49*gt3**6-.123*gt3**12+
     &          .101*gt3**13)/(1.+100.*gt3**4+.267*gt3**12)
            hs=1.25*gnz-tau*f90 
!           hs=gnz*(1.25-0.285*gt3)    
            If (hs<h0) Then
                h0=hs               
                fsr=gnz*(1.25-.57*gt3)/3. 
                fst=-gnz*(1.25-.475*gt3)   
            Endif
          Endif
!         Write(20,"(3a5,i6,4es11.4)") 'SCR2',nname(z1),
!    &         nname(z2),dexp(h0),hw,hi,hs
        Endif
!  Add succeeding screening factors
        he(j)=he(j-1)+h0
         If(iout>3) Write(15,'(a,4es13.6)') 'H2',he(j),h0,gnz,tau

       Enddo                                                             
      Return                                                            
      End subroutine qse_screen3                                                              
