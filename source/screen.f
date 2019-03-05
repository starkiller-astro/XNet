!*******************************************************************************
!  screen.f               
!
!  
!
!*******************************************************************************

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
!     Write(lun_diag,"(a3,3es15.7)") 'EOS',GNP,GSI,ZL
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
          hs=0.0
          If (gnz>=0.4) Then
            em2=a1*a2*2.0/(a1+a2)                
            tau=3.3722*(em2*z12**2*t9er)**(1./3.) 
            gt3=3.0*gnz/tau
            f90=(.0455*gt3+.348*gt3**3+9.49*gt3**6-.123*gt3**12+
     &          .101*gt3**13)/(1.+100.*gt3**4+.267*gt3**12)
            hs=1.25*gnz-tau*f90
!           hs=gnz*(1.25-0.285*gt3)    
            If (hs<h2(mu)) Then
                h2(mu)=hs               
                fsr=gnz*(1.25-.57*gt3)/3. 
                fst=-gnz*(1.25-.475*gt3)   
            Endif
          Endif
          If(idiag>3) Write(lun_diag,"(3a5,i6,6es12.5)") 'SCR2',
     &      nname(n2i(1,mu)),nname(n2i(2,mu)),mu,z,gnz,h2(mu),hw,hi,hs
        Endif
      Enddo
!     Write(lun_diag,"(2a6,1es12.5)") (nname(n2i(1,mu)),nname(n2i(2,mu)),h2(mu),mu=1,nreac(2))
      Return                                                            
      End subroutine screen2                                                              

      
