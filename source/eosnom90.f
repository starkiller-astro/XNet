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
      Return                                                                    
      End subroutine en                                                                      

      SUBROUTINE STATE(T,V,X,P,E,DLM,FYE,PT,PV,ET,EV,IEOS,M,EMAS,       
     *  YPS,BETA)                                                       
c------------------------------------------------------------------     
      implicit double precision (a-h,o-z)                               
C.....INCLUDES IDEAL ION GAS,RADIATION PRESSURE,AND NONINTERACTING FERMI
C.....-DIRAC GAS OF ELECTRONS AND POSITRONS                             
      DIMENSION X(m)
C                                                                       
C.....USES CHI OVER A                                                   
C.....CALLED BY MAIN                                                    
C.....EQUATION OF STATE AND LUMINOSITY                                  
C                                                                       
      COMMON/CONST/ PI,GRAV,ARAD,RGAS,CLN10                             
      DATA PI/3.14     /,GRAV/6.673231E-8/,ARAD/7.5641E-15/,            
     *  RGAS/8.314E7/,CLN10/2.302585093/                                
      FACT=0.0 
      CALL EN(m,X,FYE,YPS,F2Y)
      FACT=YPS-FYE                                                      
c     write(6,*) 'EN',M,FYE,YPS,F2Y
      TT=9.+dLOG10(T)                                                   
      RHO=dLOG10(1./V)                                                  
      P=0.0                                                             
      CALL EOS(RHO,TT,FACT,FYE,F2Y,P,S,E,ER,ET,PR,PT)                   
      P=10.**P                                                          
      TT=10.**TT                                                        
      EV=-ER/V                                                          
      PV=-PR*P/V                                                        
      ET=ET/TT                                                          
      PT=PT*P/TT                                                        
C      LUMINOSITY                                                       
      BETA=BET(YPS,EMAS)                                                
      DLM=DELDM(FYE,BETA,YPS,EMAS)                                      
    1 CONTINUE                                                          
      IEOS=1                                                            
      RETURN                                                            
       END                                                              
      SUBROUTINE EOS(RHO,T,EMNI,EMEI,Z2AJ,P,S,UT,DUDR,DUTR,DPDRJ,DPDTJ) 
C C   SUBROUTINE EOS(RHO,T,EMNI,EMEI,P,S,UT,CV,DSDR,DRDT,DPDRJ,         
C C  1               EPNU,J,IDEF)                                       
      Use controls
      implicit double precision(A-H,O-Z)
      COMMON /SCRN/ GNP,GSI,ZL,DZLDR,DZLDT,T9ER
      COMMON /EG/ EPG,EPGT,PB,ROB,UTB,SB,TBD                            
      COMMON /PHY/ CP,PSI,UMB                                           
      DEX(DEXARG)=DEXP(2.302585093*DEXARG)                              
      DATA SNCONS /-38.30841/, CONST /2.302585093/
      CONSTA=1./CONST
      RHOJ=RHO                                                          
      TJ=T                                                              
      EMNIL=DLOG10(EMNI)                                                
      EMEIL=DLOG10(EMEI)                                                
C     Z2AJ=EMEI**2/EMNI                                                 
  110 XI=RHOJ+TJ+7.919827                                               
      XIN  =XI+EMNIL                                                    
      SN=(1.5*XIN-2.5*RHOJ)*CONST+SNCONS                                
      PNE=DEX(XIN)                                                      
C EFFECTS OF COULOMB INTERACTIONS: SLATTERY ET AL.(1982)                
      GML=5.35692+RHOJ/3.-TJ+2.*EMEIL-5./3.*EMNIL                       
      G=10**GML                                                         
      GB=G                                                              
      If (GB.ge.0.1) Then
          CALL COULOM(RHOJ,TJ,G,GB,EMNI,EMEI,DPCL,DSCL,DPDGM,DSDG)          
          PNE=PNE*DPCL                                                      
          SN=SN+DSCL                                                        
      Else
          dpcl=0.0
          dscl=0.0
      Endif
   10 PRE=DEX(TJ*4.-14.59833)                                           
      SR=DEX(TJ*3.-RHOJ-21.91618)                                       
      CALL EQSTEL(RHOJ,TJ,EMEIL,ELAM,SPSI,SE,DPEDR,DSPSDR               
     1,DSEDR,DPEDT,DSPSDT,DSEDT,PRR,SRR,PTR,STR,PTT,STT,                
     2PSI,PSRA,DPSDT,DRDAPS,RPRA,RPAR,PSIP,0)                           
      SPSI=DEX(SPSI)                                                    
      XIE=XI+ELAM+EMEIL                                                 
      PEE=DEX(XIE)                                                      
   15 PTE=PEE+PRE+PNE                                                   
      P= DLOG10(PTE)                                                    
      DPDRJ=         (PNE+PEE*DPEDR)/PTE                                
      DPDTJ  =(PNE+PEE*DPEDT+PRE*4.)/PTE                                
      IF (GB.GE.0.1) GO TO 19                                           
      DSDT=1.5*CONST*EMNI                                               
      DSDR=-CONST*EMNI                                                  
      GO TO 20                                                          
   19 DPDGM=DPDGM*PNE/PTE                                               
      DPDRJ=DPDRJ+DPDGM/3.                                              
      DPDTJ=DPDTJ-DPDGM                                                 
      DSDT=(1.5-DSDG)*CONST*EMNI                                        
      DSDR=(-1.0+DSDG/3.)*CONST*EMNI                                    
   20 DSDT=DSDT+SR*3.*CONST+EMEI*DSEDT                                  
      DSDR=DSDR-SR*CONST+EMEI*DSEDR                                     
      DRDP=1./DPDRJ                                                     
      DRDT=-DPDTJ/DPDRJ                                                 
      UMB=PRE/PTE                                                       
      S=SN*EMNI+SE*EMEI+SR                                              
      CV=DSDT/CONST                                                     
      CP=CV+DSDR*DRDT/CONST                                             
C*****SCREENING                                                         
      T9ER=10.**(9.0-T)                                                 
      T6=(T-6.0)*CONST                                                  
      DD=RHO*CONST-T6*3.                                                
      GNP=DEXP((DD+EMEIL*CONST)/3.)*0.22747                             
      ZZ=EMEI/EMNI                                                      
      ZZ3=ZZ**(1./3.)                                                   
      ES=CONST/PSRA*EMEI                                                
      Z2ES=Z2AJ+ES                                                      
      GSI=EMNI/Z2ES                                                     
      ZLAM=-1.671845+0.5*(DD+DLOG(Z2ES))                                
      ZL=DEXP(ZLAM)
      If(idiag>=0) Write(lun_diag,'(a14,9es12.5)') 'Nomoto EOS',
     & 10.0**(t-9.0), 10.0**rho, emei, z2aj/EMNI, zz, const/psra,
     & z2es/sqrt(emni), zl, gnp
      DZLDR=0.5+0.5*ES/Z2ES*RPRA
      DZLDT=-1.5+0.5*ES/Z2ES*RPAR                                       
C-----ZONE FLASH                                                        
      RHOI=DEX(-RHOJ)                                                   
      PER=PEE*RHOI                                                      
      TEE=DEX(TJ+7.919827+EMEIL)                                        
      PNN=DEX(RHOJ+TJ+7.919827+EMNIL)                                   
      DU1=(DPCL-1.D0)*3.+1.                                             
      UN=1.5*PNN*RHOI                                                   
      IF (GB.GT.0.1) UN=UN*DU1                                          
      UE=TEE*SPSI-PER                                                   
      UR=3.0*PRE*RHOI                                                   
      UT=UN+UE+UR                                                       
      DUEDR=TEE*DSEDR*CONSTA+PER                                        
      DUEDT=TEE*DSEDT*CONSTA                                            
      DUDR=DUEDR-UR                                                     
      DUTR=UN+DUEDT+4.*UR                                               
C     PME=10.**PB                                                       
C     DTIME=1.D0                                                        
C     EG1=(UT-UTB)/DTIME                                                
C     EG2=PME*(RHOI-DEX(-ROB))/DTIME                                    
      DUDP=DUDR*DRDP*CONST                                              
      DUDT=(DUTR+DUDR*DRDT)*CONST                                       
C     EPG=-(EG1+EG2)                                                    
C     EGPR=-PME*RHOI/DTIME*CONST                                        
C     EPGP=-DUDP/DTIME-EGPR*DRDP-0.5*EG2*CONST                          
C     EPGT=-DUDT/DTIME-EGPR*DRDT                                        
C     BKHT=10.**((T+TBD)*0.5)*8.31433E7                                 
C     EPG=-(S-SB)*BKHT                                                  
C     EPGT=-DSDT*BKHT+0.5*EPG*CONST-DSDR*DRDT*BKHT                      
      RETURN                                                            
      END                                                               
      SUBROUTINE SCREEN(Z1,Z2,A1,A2,H0,FSR,FST)                         
C-----WEAK AND INTERMEDIATE SCREENING FACTOR BY GRABOSKE ET AL.         
      Use controls
      implicit double precision (a-h,o-z)
      COMMON /SCRN/ GNP,GSI,ZL,DZLDR,DZLDT,T9ER                         
      ZZ=Z1*Z2  
c     write(51,*) z1,z2,zl          
      IF(ZZ.GT.0.0) GOTO 5                                   
      H0=0.0                                                            
      GOTO 30                                                           
    5 Z=ZL*ZZ                                                           
      GSIZ=GSI*ZZ                                                       
      GSB=GSIZ**0.86                                                    
      IF (GSIZ.GT.1.) F=0.38*((1.+GSIZ)**1.86-GSB*GSIZ-1.)/GSB          
      IF (GSIZ.LE.1.) F=0.61943*GSIZ**0.07                              
      HW=Z*(1.+Z*(dLOG(Z)+0.8364))                                      
      HI=F*Z**0.86                                                      
      IF (HW.GT.HI) GO TO 10
c     write(51,*) 'h0=hw',z,f,hw                                           
      H0=HW                                                             
      DHDZ=Z*(1.+2.*Z*(dLOG(Z)+1.3364))                                 
      GO TO 20                                                          
   10 H0=HI 
c     write(51,*) 'ho=hi',z,f,hi                                                  
      DHDZ=0.86*HI                                                      
   20 FSR=DHDZ*DZLDR                                                    
      FST=DHDZ*DZLDT                                                    
C-----STRONG SCREENING BY ITOH ET AL.(1990)                             
      Z13=Z1**(1./3.)                                                   
      Z23=Z2**(1./3.)
      SCR=0.0
      GNZ=GNP*ZZ*2./(Z13+Z23)                                           
      IF (GNZ.LT.0.4) GO TO 30                                          
      EM2=A1*A2*2./(A1+A2)                                              
      TAU=3.3722*(EM2*ZZ**2*T9ER)**(1./3.)                              
      GT=GNZ/TAU                                                        
      GT3=3.*GT                           
      NUMF=.0455*GT3+.348*GT3**3.+9.49*GT3**6.-.123*GT3**12.+
     #  .101*GT3**13.   
      DENF=1.+100.*GT3**4.+.267*GT3**12.
      F90=NUMF/DENF
      SCR=1.25*GNZ !-TAU*F90
      SCRO=GNZ*(1.25-0.855*GT)
      tst=dexp(scr)/dexp(scro)
      IF (SCR.GE.H0) GO TO 30 
      H0=SCR                   
      FSR=GNZ*(1.25-1.71*GT)/3.
      FST=-GNZ*(1.25-1.425*GT)                                          
30    If(idiag>3) Write(lun_diag,"(a5,3i6,6es12.5)") 'SCRN',
     &      int(z1),int(z2),0,z,gnz,h0,hw,hi,SCR
      RETURN
      END                                                               
      SUBROUTINE SCREENV(ZP,ZA,N2I,H0,inam,N,ny)                                
C-----WEAK AND INTERMEDIATE SCREENING FACTOR BY GRABOSKE ET AL.         
      implicit double precision (a-h,o-z)                               
      DIMENSION ZP(ny),ZA(ny),H0(n),N2I(4,n)                              
      character*5 inam(0:ny)
      COMMON /SCRN/ GNP,GSI,ZL,DZLDR,DZLDT,T9ER                         
      common /flag/ incwk,iscrn,iwrt,idbg,iconvc,kstmx,knrmx,nzone
c     write(lun_diag,*) 'EOS',GNP,GSI,ZL
      DO 30 MU=1,N                                                      
        Z1=ZP(N2I(1,MU))                                                  
        Z2=ZP(N2I(2,MU))                                                  
        A1=ZA(N2I(1,MU))                                                  
        A2=ZA(N2I(2,MU))                                                  
        ZZ=Z1*Z2    
        IF(ZZ.GT.0.0.and.iscrn.eq.1) GOTO 5              
      H0(MU)=0.0                                                        
      GOTO 30                                                           
    5 Z=ZL*ZZ                                                           
      GSIZ=GSI*ZZ                                                       
      GSB=GSIZ**0.86                                                    
      IF (GSIZ.GT.1.) F=0.38*((1.+GSIZ)**1.86-GSB*GSIZ-1.)/GSB          
      IF (GSIZ.LE.1.) F=0.61943*GSIZ**0.07                              
      HW=Z*(1.+Z*(dLOG(Z)+0.8364))                                      
      HI=F*Z**0.86                                                      
      IF (HW.GT.HI) GO TO 10                                            
      H0(MU)=HW                                                         
      DHDZ=Z*(1.+2.*Z*(dLOG(Z)+1.3364))                                 
      GO TO 20                                                          
   10 H0(MU)=HI                                                         
      DHDZ=0.86*HI                                                      
   20 FSR=DHDZ*DZLDR                                                    
      FST=DHDZ*DZLDT                                                    
C-----STRONG SCREENING BY ITOH ET AL.(1990)                             
      Z13=Z1**(1./3.)                                                   
      Z23=Z2**(1./3.)                                                   
      GNZ=GNP*ZZ*2./(Z13+Z23)                                           
      IF (GNZ.LT.0.4) GO TO 30                                          
      EM2=A1*A2*2./(A1+A2)                                              
      TAU=3.3722*(EM2*ZZ**2*T9ER)**(1./3.)                              
      GT=GNZ/TAU                                                        
      GT3=3.*GT                                                         
      NUMF=.0455*GT3+.348*GT3**3.+9.49*GT3**6.-.123*GT3**12.+
     #  .101*GT3**13.   
      DENF=1.+100.*GT3**4.+.267*GT3**12.
      F90=NUMF/DENF
      SCR=1.25*GNZ-TAU*F90 
c     SCR=GNZ*(1.25-0.855*GT)                                           
      IF (SCR.GE.H0(MU)) GO TO 30                                       
      H0(MU)=SCR                                                        
      FSR=GNZ*(1.25-1.71*GT)/3.                                         
      FST=-GNZ*(1.25-1.425*GT)                                          
c  30 write(18,991) 'srnv',mu,z1,a1,z2,a2,dexp(h0(mu))
c  30 write(18,992) 'srnm',mu,inam(N2I(1,MU)),inam(N2I(2,MU)),
c    &              dexp(h0(mu)),hw,hi,scr
      CONTINUE                                                          
   30 CONTINUE                                                          
c     DO 40 MU=1,N                                                      
c       Z1=ZP(N2I(1,MU))                                                  
c       Z2=ZP(N2I(2,MU))
c       tstscr=dexp(h0(mu))  
c       write(lun_diag,*) inam(n2i(1,mu)),inam(n2i(2,mu)),h0(mu)
c  40 Continue                                              
      RETURN                                                            
  991 Format(a5,1x,i5,1x,4(f3.0,1x),4(1x,1pd11.4))
  992 Format(a5,1x,i5,1x,2(a5,1x),4(1x,1pd11.4))
  993 Format(a5,1x,i5,1x,2(i5,1x),4(1x,1pd11.4))
      END                                                               
      SUBROUTINE COULOM(RHOJ,TJ,G,GB,EMNI,EMEI,DPCL,DSCL,DPDGM,DSDGM)   
      implicit double precision(A-H,O-Z)                                
      DATA A1,A2,A3,A4,A5/-0.897744,0.95043,0.18956,-0.81487,2.58204/   
      DATA  B1,B2,B3,B4,B5/-0.895929,-1612.5,3225.0,5.25E3,0.0/         
      DATA CONST/2.302585/ 
      IF (G.LT.1.0) GO TO 10                                            
C EFFECTS OF COULOMB INTERACTIONS: LIQUID: SLATTERY ET AL.(1982)        
      GML=DLOG10(G)                                                     
       G4=G**0.25                                                       
       DPCL=1.+(A1*G+A2*G4+A3/G4+A4)/3.                                 
       DSCL=-3.*A2*G4+5.*A3/G4+A4*(1.-GML*CONST)+A5                     
C      DSCL=DSCL+0.01105671D0                                           
       DPDGM=(4.*A1*G+A2*G4-A3/G4)/(12.*DPCL)                           
       DSDGM=(-3.*A2*G4-5.*A3/G4-A4*4.)*0.25                            
      RETURN                                                            
C EFFECTS OF COULOMB INTERACTIONS: VAN HORN (1970)                      
   10 G12=DSQRT(G)                                                      
      G32=G12*G                                                         
      GF1=1./(1.+0.142*G)                                               
      GF2=1./(1.+0.575*G)                                               
      SF1=1./(1.+G12)                                                   
      SF2=1./(1.+1.308*G32)                                             
      S3=DSQRT(3.0D+0)                                                  
      DPCL=1.-0.113*G32*(DSQRT(GF1)+1.54*GF2**1.5)                      
      SF=1.+G32/(2.*S3)*(0.015+0.585*SF1+0.4*SF2)                       
      DSCL=-DLOG(SF)                                                    
      DPDGM=-0.0565/DPCL*G32*((3.+0.284*G)*GF1**1.5+4.62*GF2**2.5)      
      DSDGM=-G32/(SF*4.)*S3*(0.015+0.195*(3.+2.*G12)*SF1**2+0.4*SF2**2) 
      RETURN                                                            
       END                                                              
      SUBROUTINE EQSTEL (R,T,EM,EL,U,S,PR,UR,SR,PT,UT,ST,PRR,SRR,PTR,   
     1STR,PTT,STT,PSI,PSRA,PSAR,RAPS,RPRA,RPAR,PSIP,MODE)               
      implicit double precision(A-H,O-Z)                                
      DIMENSION X(4),Y(4,17),W(17),XP(4),YP(4,17),EPM(9)                
      COMMON / EQSDAT / AL(10),PS(18),FT(18,10),YT(9,18,10),AX(18,10),  
     1AG(18,10),AP1(18,10),AP2(18,10),KM,JM                             
      COMMON /IPA/ IXA,IXF,IXR,IST                                      
      DATA C/2.302585/,CR/ .4342945/,J,K/3,2/ 
      DEX(ARG)=DEXP(2.302585*ARG)                                       
C-----ELECTRON PAIRS ARE INCLUDED                                       
      M=MODE                                                            
      F=R+EM-T*1.5+8.043244                                             
      A=DEX(T-9.773064)                                                 
  110 IF (AL(K).GT.A) GO TO 100                                         
      IF (K.GE.KM) GO TO 210                                            
      K=K+1                                                             
      GO TO 110                                                         
  210 IEQ=IEQ+1                                                         
      IF (IEQ.LE.100) WRITE (6,900) R,T                                 
      IF (IEQ.EQ.100) WRITE (6,901)                                     
  901 FORMAT (53H0THE ABOVE MESSAGE WILL BE TERMINATED AFTER 100 LINES, 
     1//)                                                               
  900 FORMAT (23H OUT OF RANGE IN EQSTEL,10X,10HLOG(RHO) =,F10.6,10X,   
     18HLOG(T) =,F10.6)                                                 
  100 K=K-1                                                             
      IF (K.LE.2) GO TO 120                                             
      IF (AL(K).GT.A) GO TO 100                                         
      K4=K+2                                                            
      IF (K4.GT.KM) K4=KM                                               
      K0=K4-4                                                           
      GO TO 130                                                         
  120 K0=0                                                              
  130 IF (FT(J,K).GT.F) GO TO 200                                       
      IF(J.GE.JM) GO TO 240                                             
      J=J+1                                                             
      GO TO 130                                                         
  240 CALL EXTDEG          (F,A,EL,U,S,PR,UR,SR,PT,UT,ST,PRR,SRR,PTR,   
     1STR,PTT,STT,PSI,PSRA,PSAR,RAPS,RPRA,RPAR,M)                       
      PSIP=-(PSI+2./A)                                                  
      RETURN                                                            
  200 IF (J.LE.3) GO TO 220                                             
      J=J-1                                                             
      IF (FT(J,K).GT.F) GO TO 200                                       
      FLIN=(FT(J,K+1)-FT(J,K))*(A-AL(K))/(AL(K+1)-AL(K))+FT(J,K)        
      IF (FLIN.GT.F) J=J-1                                              
      J4=J+2                                                            
      IF (J4.GT.JM) J4=JM                                               
      J0=J4-4                                                           
      GO TO 230                                                         
  220 J0=0                                                              
  230 IXB=-1                                                            
      LYI=1                                                             
      LYE=9                                                             
      MI=1                                                              
      IF (M.EQ.0) GO TO 1230                                            
      LYI=2                                                             
      LYE=3                                                             
      MI=3                                                              
 1230 DO 300 LF=1,4                                                     
      LFJ=LF+J0                                                         
      DO 301 LA=1,4                                                     
      LAK=LA+K0                                                         
      X(LA)=AX(LFJ,LAK)                                                 
      Y(LA,10)=FT(LFJ,LAK)                                              
      DO 301 LY=LYI,LYE                                                 
      Y(LA,LY)=YT(LY,LFJ,LAK)                                           
  301  CONTINUE                                                         
      P=PS(LFJ)                                                         
      IXF=0                                                             
      IXA=0                                                             
      IF (P.GE.1.) GO TO 13                                             
      IXA=1                                                             
      DO 11 LA=1,4                                                      
      LAK=LA+K0                                                         
      Y(LA,11)=Y(LA,7)                                                  
      ALLAK=AL(LAK)                                                     
      CALL PMOD(P,ALLAK,  EPM) 
      DO 10 JP=1,7                                                      
      IF (JP/3*3.EQ.JP) GO TO 10                                        
      Y(LA,JP)= Y(LA,JP)-EPM(JP)                                        
   10 CONTINUE                                                          
      Y(LA,3)=DLOG10( Y(LA,3)/EPM(3))                                   
      Y(LA,6)=DLOG10( -Y(LA,6)*CR+EPM(6))-EPM(2)                        
   11 CONTINUE                                                          
   13 IX=0                                                              
      IF (P.LT.1.999999) GO TO 350                                      
      IX=1                                                              
      AP=A*P                                                            
      APU=AP+1.                                                         
      APUH=AP*0.5+1.                                                    
      V=              DLOG10(APUH)                                      
      GO TO 309                                                         
  350 IF (IXB.EQ.IX) GO TO 309                                          
      V=              DLOG10(1.+A)                                      
  309 IY=0                                                              
      IF (A.LT.0.2) IY=1                                                
      IF((IY.EQ.0).OR.(IX.EQ.0)) GO TO 310                              
      DO 311 LA=1,4                                                     
      LAK=LA+K0                                                         
      Y(LA,3)=Y(LA,3)*AP1(LFJ,LAK)                                      
      IF (M.NE.0) GO TO 311                                             
      Y(LA,4)=Y(LA,4)/AG(LFJ,LAK)                                       
      Y(LA,5)=Y(LA,5)/AP1(LFJ,LAK)                                      
      Y(LA,6)=Y(LA,6)*AP2(LFJ,LAK)                                      
      Y(LA,9)=Y(LA,9)*AP1(LFJ,LAK)                                      
  311 CONTINUE                                                          
      AP1L=APUH/APU                                                     
      IF (M.NE.0) GO TO 330                                             
      AP2L=AP1L/(1.+1./APU**2)                                          
      XI2=AP*(2.+AP)                                                    
      IF (XI2.LT.0.1 ) GO TO 320                                        
      XI=DSQRT(XI2)                                                     
      PHI= DLOG(APU+XI)*3.+XI*(XI2+XI2-3.)*APU                          
      AGL=2.666667*XI*XI2*XI2/(APU*PHI)-1.                              
      GO TO 330                                                         
  320 AGL=.6666667       +((0.134653*XI2+.1927438)*XI2-.2380952)*XI2    
  330 CALL INT4 (X,Y,V,W,MI)                                            
      W(3)=W(3)/AP1L                                                    
      IF (M.NE.0) GO TO 340                                             
      W(4)=W(4)*AGL                                                     
      W(5)=W(5)*AP1L                                                    
      W(6)=W(6)/AP2L                                                    
      W(9)=W(9)/AP1L                                                    
      W(12)=W(12)-C*W(3)/APU                                            
      GO TO 340                                                         
  310 CALL INT4 (X,Y,V,W,MI)                                            
      IF (IXA.EQ.0) GO TO 340                                           
      CALL PMOD(P,A,EPM)                                                
      DO 20 JP=1,7                                                      
      IF (JP/3*3.EQ.JP) GO TO 20                                        
      W(JP)=W(JP)+EPM(JP)                                               
   20 CONTINUE                                                          
      W(3)=10**W(3)*EPM(3)                                              
      W(6)=-C*(10**(W(6)+EPM(2))-EPM(6))                                
  340 XP(LF)=W(10)                                                      
      DO 341 LY=LYI,LYE                                                 
  341 YP(LF,LY)=W(LY)                                                   
      IF (M.NE.0) GO TO 300                                             
      IF (IX.EQ.0) DXDA=A/(A+1.)                                        
      IF (IX.EQ.1) DXDA=AP*0.5/APUH                                     
      YP(LF,10)=W(11)*DXDA                                              
      YP(LF,11)=W(12)*DXDA                                              
      YP(LF,12)=W(13)*DXDA                                              
  300 IXB=IX                                                            
      MI=MI+1                                                           
      PS2=PS(J0+2)                                                      
      IXA=0                                                             
      IXF=0                                                             
      IF (PS2.GE.3.9) GO TO 400                                         
      IXF=1                                                             
      IXR=0                                                             
      PSAJ=PS(J0+4)+1./A                                                
      IF (PSAJ.GT.1.5) IXR=1                                            
      IST=0                                                             
      IF ((PS2.GT.-4.01).OR.(A.GT.0.20001)) IST=1                       
      DO 32 LF=1,4                                                      
      LFJ=LF+J0                                                         
      P=PS(LFJ)                                                         
      YP(LF,13)=YP(LF,4)                                                
      YP(LF,14)=YP(LF,6)                                                
      YP(LF,15)=YP(LF,7)                                                
      YP(LF,16)=YP(LF,9)                                                
      YP(LF,17)=P                                                       
      CALL PMOD(P,A,EPM)                                                
      DO 31 JP=1,8                                                      
      IF (JP/3*3.EQ.JP) GO TO 31                                        
      IF ((JP.GE.7).AND.(IXR.EQ.0)) GO TO 30                            
      YP(LF,JP)=YP(LF,JP)-EPM(JP)                                       
      GO TO 31                                                          
   30 YP(LF,JP)=YP(LF,JP)/EPM(JP)                                       
   31 CONTINUE                                                          
      YP(LF,3)=DLOG10(YP(LF,3)+EPM(3))-EPM(2)                           
      YP(LF,6)=DLOG10(-YP(LF,6)*CR+EPM(6))-EPM(2)                       
      IF (IST.EQ.1) YP(LF,9)=DLOG10(YP(LF,9))-EPM(9)                    
   32 CONTINUE                                                          
  400 CALL INT4 (XP,YP,F,W,MI)                                          
      IF (IXF.EQ.0) GO TO 420                                           
      CALL PMOD(W(17),A,EPM)                                            
      DO 41 JP=1,8                                                      
      IF (JP/3*3.EQ.JP) GO TO 41                                        
      IF ((JP.GE.7).AND.(IXR.EQ.0)) GO TO 40                            
      W(JP)=W(JP)+EPM(JP)                                               
      GO TO 41                                                          
   40 W(JP)=W(JP)*EPM(JP)                                               
   41 CONTINUE                                                          
      W(3)=10**(W(3)+EPM(2))-EPM(3)                                     
      W(6)=-C*(10**(W(6)+EPM(2))-EPM(6))                                
      IF (IST.EQ.1) W(9)=10**(W(9)+EPM(9))                              
  420 IF (M.NE.0) GO TO 1420                                            
      EL=W(1)                                                           
 1420 U=W(2)                                                            
      S=W(3)                                                            
      IF (M.NE.0) RETURN                                                
      PL=W(4)+1.                                                        
      IF (DABS(PL).LT.1.E-7) PL=1.E-7                                   
      PM=W(7)+1.                                                        
      UR=W(5)                                                           
      UT=W(8)                                                           
      SR=W(6)                                                           
      ST=W(9)                                                           
      PSR=DEX(W(1))*PL*C                                                
      PSRA=PSR                                                          
      EXW2=DEX(W(2))                                                    
      PSI=EXW2-W(3)                                                     
      PSIP=-(PSI+2./A)                                                  
      PSA=C*EXW2*W(8)-W(9)                                              
      AT=PSA/PSR                                                        
      RAPS=-AT                                                          
      PLI=CR/PL                                                         
      RPRA=-PLI*W(13)-W(4)                                              
      RPAR=-PLI*W(15)-W(7)                                              
      PR=PL                                                             
      PSAR=PSA                                                          
      AT0=-W(12)-1.5                                                    
      ELTT=W(10)+AT0*W(15)                                              
      PT=PM                                                             
      STT=W(11)+AT0*W(16)                                               
      PRR=C*(C*PL*PL+W(13))                                             
      PTT=C*(C*PM*PM+ELTT)                                              
      PTR=C*(C*PL*PM+W(15))                                             
      STR=W(16)                                                         
      SRR=W(14)                                                         
      RETURN                                                            
      END                                                               
      SUBROUTINE INT4 (X,Y,P,Q,MODE)                                    
      implicit double precision (A-H,O-Z)                               
      COMMON /IPA/ IXA,IXF,IXR,IST                                      
      DIMENSION X(4),Y(4,17),Q(17)  
      D2(Z)= (A*Z*3.+B+B)*Z+C                                           
      X2=X(2)                                                           
      X32=X(3)-X2                                                       
      X12=X(1)-X2                                                       
      X42=X(4)-X2                                                       
      X31=X32-X12                                                       
      X41=X42-X12                                                       
      X43=X42-X32                                                       
      R=P-X2                                                            
      M=MODE                                                            
      IF (M.GE.3) GO TO 3                                               
      KI=1                                                              
      KE=10                                                             
      IF (IXA.EQ.1) KE=11                                               
      IF (M.EQ.2) KE=12                                                 
      IF (IXF.EQ.1) KE=17                                               
   10 DO 20 K=KI,KE                                                     
      Y12=(Y(1,K)-Y(2,K))/X12                                           
      Y32= (Y(3,K)-Y(2,K))/X32                                          
      E32=(Y32-Y12)/X31                                                 
      E42=((Y(4,K)-Y(2,K))/X42-Y12)/X41                                 
      A=(E42-E32)/X43                                                   
      B=E32-(X32+X12)*A                                                 
      C=Y32-(A*X32+B)*X32                                               
      Q(K)=((A*R+B)*R+C)*R+Y(2,K)                                       
      IF (M.GE.3) GO TO 20                                              
      IF (M.EQ.2) GO TO 11                                              
      IF ((IXA.EQ.0).AND.(K.EQ.7)) Q(11)=D2(R)                          
      IF ((IXA.EQ.1).AND.(K.EQ.11))Q(11)=D2(R)                          
      IF (K.EQ.9) Q(12)=D2(R)                                           
      IF (K.EQ.10)Q(13)=D2(R)                                           
      GO TO 20                                                          
   11 IF (K.LT.4) GO TO 20                                              
      IF (IXF.EQ.1) GO TO 30                                            
      IF (K.NE.4) GO TO 12                                              
      Q(13)=D2(R)                                                       
      GO TO 20                                                          
   12 IF (K.NE.6) GO TO 13                                              
      Q(14)=D2(R)                                                       
      GO TO 20                                                          
   13 IF (K.EQ.7) Q(15)=D2(R)                                           
      IF (K.EQ.9) Q(16)=D2(R)                                           
      GO TO 20                                                          
   30 IF ((K.GE.13).AND.(K.LE.16)) Q(K)=D2(R)                           
   20 CONTINUE                                                          
      IF (M.NE.3) RETURN                                                
      IF (KE.NE.3) RETURN                                               
      KI=10                                                             
      KE=10                                                             
      GO TO 10                                                          
    3 KI=2                                                              
      KE=3                                                              
      GO TO 10                                                          
      END                                                               
      SUBROUTINE PMOD (PS,A,EPM)                                        
      implicit double precision(A-H,O-Z)                                
      COMMON /IPA/ IXA,IXF,IXR,IST                                      
      DIMENSION EPM(9)  
      IF (A.LT.1.E-6) GO TO 30                                          
      AI=1./A                                                           
      PSA=2.*(PS+AI)                                                    
      IF (PSA.GT.170.) GO TO 30                                         
      R =DSQRT(1.+4.*DEXP(-PSA))                                        
      GO TO 40                                                          
   30 R=1.                                                              
      AI=1.                                                             
   40 W =R *(4.+3.*A)*5./8.+(R -1.)*AI                                  
      IF ((IXF.EQ.1).AND.(A.GT.0.8)) W=2.5*R                            
      S0=   PS +DLOG((R +1.)*0.5)                                       
      S =W-S0                                                           
      RR=1.-1./(R*R)                                                    
      FA=(W/R-1.)/R-RR*AI                                               
      WA=3.*(1.-0.5/(R*R))-3.75/W                                       
      EPM(1)=DLOG10(R)                                                  
      EPM(2)=DLOG10(W)                                                  
      EPM(4)=-RR                                                        
      EPM(5)=-RR*(1.+A*W)/W*AI                                          
      EPM(6)=FA                                                         
      EPM(7)=(W-R-FA)/R                                                 
      IF ( IXA .EQ. 1 )  GO TO 10                                       
      EPM(3)=S0                                                         
      IF (IXR.EQ.0) EPM(7)=RR                                           
      IF (IXR.EQ.0) EPM(8)=RR*R/W                                       
      IF (IXR.EQ.1) EPM(8)=WA                                           
      IF (IST.EQ.1) EPM(9)=DLOG10((W*WA+FA)*2.302585)                   
      GO TO 20                                                          
   10 EPM(3)=S                                                          
   20 RETURN                                                            
      END                                                               
      SUBROUTINE EXTDEG (LF1,A1,LG,LU,S1,PR1,LULR,SR1,PT1,LULT,SLT,PRR, 
     1SRR,PTR,STR,PTT,STT,PS,PSRA,PSAR,RAPS,RPRA,RPAR,MODE)             
      implicit double precision(A-H,L,O-Z)                              
      DATA PAI2/9.869604/,C/2.302585/,C2/5.301898/ 
      LF=LF1                                                            
      A=A1                                                              
      LX2A=(LF+0.1761)/1.5                                              
      N=0                                                               
      A2=A+A                                                            
   30 X2AR=DEXP(1.151293*LX2A)                                          
      X2A=X2AR*X2AR                                                     
      X2=X2A*A2                                                         
      X2AR3=X2AR*X2A                                                    
      X22=X2+X2                                                         
      X2A2 =X2A*X2A                                                     
      X2A28=X2A2*8.                                                     
      X22U=X22+1.                                                       
      X22UP=PAI2*X22U                                                   
      W=1.+X22UP/X2A28                                                  
      F=X2AR3*W/1.5                                                     
      LFR= DLOG10(F)                                                    
      LFLX=(3.+PAI2*(X22-1.)/X2A28)/W                                   
      IF (DABS(LFR-LF).LT.1.E-5) GO TO 10                               
      IF (N.GT.100) GO TO 20                                            
      N=N+1                                                             
      LX2A=LX2A+2.*(LF-LFR)/LFLX                                        
      GO TO 30                                                          
   20 IEND=1                                                            
      WRITE (6,21) LF,A                                                 
   21 FORMAT (22H0ERROR IN EXTDEG, LF =,F10.6,10X,3HA =,F10.6)          
      RETURN                                                            
   10 X2U=X2+1.                                                         
      UPX2=DSQRT(X2U)                                                   
      X2A2A=X2A+X2A                                                     
      S=PAI2*UPX2*(1.-PAI2*(X2*8.+11.)/(X2A2*60.))/X2A2A                
      IF(X2.GT.1.E-3) GO TO 11                                          
      PSI=X2A*(1.-X2*0.25)                                              
      GO TO 12                                                          
   11 PSI=(UPX2-1.)/A                                                   
   12 U=S+PSI                                                           
      LU= DLOG10(U)                                                     
      S1=S                                                              
      PS=PSI                                                            
      IF (MODE.NE.0) RETURN                                             
      LFLA=-X2AR3*(1.-X22UP/(24.*X2A2))/F                               
      X=DSQRT(X2)                                                       
      X2AR5=X2AR3*X2A                                                   
      IF (X2.LT.0.06 ) GO TO 13                                         
      SFA=(3.* DLOG(UPX2+X)+X*(X22-3.)*UPX2)/A2**2.5                    
      GO TO 14                                                          
   13 SFA=1.6*X2AR5*(1.-X2/2.8*(1.-X2*7./12.))                          
   14 W=X2AR*UPX2*PAI2                                                  
      G=(SFA+W)/6.                                                      
      LG= DLOG10(G)-LF                                                  
      LXLA=-X22UP/(X2A2*12.)                                            
      W=W/SFA                                                           
      LGLX=8.*X2AR5*(1.+X22UP/X2A28-W)/(SFA*UPX2)                       
      LGLA25=W+W                                                        
      LGLR=LGLX/LFLX-1.                                                 
      PT=LGLA25+LXLA*LGLX                                               
      SLX=-PAI2*(X2+2.-PAI2*((24.*X2+87.)*X2+66.)/(X2A2*60.))/(X2A2A    
     1*UPX2)                                                            
      SLA=PAI2*UPX2*(1.-PAI2*(X2*8.+11.)/(20.*X2A2))/X2A2A              
      ULX=SLX+X2A2A/UPX2                                                
      ULA=SLA-PSI                                                       
      LULX=ULX/U                                                        
      LULA=ULA/U                                                        
      SLR=C*SLX/LFLX                                                    
      LULR=LULX/LFLX                                                    
      SLT=C*(SLA+LXLA*SLX)                                              
      LULT=LULA+LXLA*LULX                                               
      PSX=X2A2A/UPX2                                                    
      PSRA=C*PSX/LFLX                                                   
      PSAR=C*(LXLA*PSX-PSI)                                             
      RAPS=-PSAR/PSRA                                                   
      LRPLA=1.-PAI2*X2U/(X2A2*3.)                                       
      LRPLX=-(X2+2.)*LRPLA/X2U                                          
      RPRA=LRPLX/LFLX                                                   
      RPAR=LRPLA+LXLA*LRPLX                                             
      STT=C2*S                                                          
      STR=C*SLR                                                         
      SRR=-C*((X2+2.)*SLR-C*S*X2/(1.5*X2U))/(X2U*3.)                    
      SR1=SLR                                                           
      W=SFA*SFA*X2U                                                     
      LAX=PAI2*X2AR*(SFA*((X2+X2+3.)*X2+2.)/UPX2-X2AR5*(X2+2.)*8.)      
     1/(1.5*W)                                                          
      LAA=PT+PT                                                         
      LRT=LAX/LFLX                                                      
      LTT=LAA+LXLA*LAX                                                  
      PR=LGLR+1.                                                        
      PRR=C2*PR*PR                                                      
      PTT=C2*(PT*PT+LTT)                                                
      PTR=C2*(PT*PR+LRT)                                                
      PR1=PR                                                            
      PT1=PT                                                            
      RETURN                                                            
      END                                                               
      function bet(yps,emas)
      implicit double precision (a-h,o-z)
c ..... solves the equation for beta (1-beta=.....)
c .................................................................
      bet=0.
   10 bet1=bet-eff(bet,yps,emas)/efs(bet,yps,emas)
      if(dabs((bet1-bet)/bet1).le.1.e-12) goto 20
      bet=bet1
      goto 10
   20 bet=bet1
      return
      end
      function deldm(ye,bet,yps,emas)
      implicit double precision (a-h,o-z)
c ...  luminosity l/m eddington limit
c ...............................................................
      deldm=1.88e2/ye*emas**2*(bet/yps)**4
      return
      end
      function eff(bet,yps,emas)
      implicit double precision (a-h,o-z)
      eff=0.00298*(bet/yps)**4*emas**2+bet-1.
      return
      end
      function efs(bet,yps,emas)
      implicit double precision (a-h,o-z)
      efs=0.00298*4.*bet**3/yps**4*emas**2+1.
      return   
      end
C-----BLOCK DATA FOR ELECTRON PAIRS                                     
      BLOCK DATA                                                        
      implicit double precision (a-h,o-z)
      COMMON /EQSDAT/ EQ1(28),EQ2(90), EQ3(90),                         
     1 EQ4 (108), EQ5 (108), EQ6 (108), EQ7 (108), EQ8 (108), EQ9 (108),
     1 EQ10(108), EQ11(108), EQ12(108), EQ13(108), EQ14(108), EQ15(108),
     1 EQ16(108), EQ17(108), EQ18(108), EQ19( 90), EQ20( 90), EQ21( 90),
     1 EQ22( 90), EQ23( 90), EQ24( 90), EQ25( 90), EQ26( 90), KMJM(2)   
      DATA EQ1  /0., 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1., 1.5, 2.,       
     1-9.,-8.,-6.,-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.,6.,8.,11.,15.,20.,30. /
      DATA EQ2  /                                                       
     1 -3.961120, -3.526860, -2.658600, -1.792420, -1.362840, -0.940860,
     2 -0.536850, -0.168710,  0.145000,  0.398360,  0.599550,  0.761220,
     3  1.006220,  1.186960,  1.390430,  1.590420,  1.776790,  2.040180,
     4 -3.921400, -3.487140, -2.618860, -1.752580, -1.322790, -0.900270,
     5 -0.494930, -0.123890,  0.195010,  0.456050,  0.666870,  0.839510,
     6  1.108430,  1.314030,  1.554510,  1.801900,  2.043730,  2.406370,
     7 -3.883600, -3.449330, -2.581040, -1.714670, -1.284700, -0.861690,
     8 -0.455170, -0.081600,  0.241830,  0.509390,  0.728230,  0.909710,
     9  1.197050,  1.420560,  1.685670,  1.961450,  2.232650,  2.639440,
     O -3.813210, -3.378950, -2.510640, -1.644100, -1.213820, -0.790000,
     1 -0.381500, -0.003700,  0.327150,  0.605210,  0.836540,  1.031380,
     2  1.345120,  1.592660,  1.888460,  2.196600,  2.498360,  2.946270,
     3 -3.748980, -3.314710, -2.446380, -1.579720, -1.149180, -0.724710,
     4 -0.314630,  0.066530,  0.403200,  0.689260,  0.929860,  1.134250,
     5  1.465960,  1.728850,  2.042980,  2.369080,  2.686580,  3.154030/
      DATA EQ3  /                                                       
     1 -3.635650, -3.201380, -2.333020, -1.466160, -1.035250, -0.609800,
     2 -0.197360,  0.188770,  0.533880,  0.831300,  1.084610,  1.301710,
     3  1.656230,  1.937560,  2.272570,  2.617970,  2.951550,  3.437900,
     4 -3.515890, -3.081610, -2.213240, -1.346210, -0.914980, -0.488700,
     5 -0.074250,  0.316060,  0.668130,  0.974690,  1.237910,  1.464620,
     6  1.835620,  2.129630,  2.478310,  2.835750,  3.178950,  3.676140,
     7 -3.414070, -2.979790, -2.111400, -1.244260, -0.812800, -0.385950,
     8  0.029880,  0.423060,  0.779850,  1.092480,  1.362170,  1.595000,
     9  1.976210,  2.277810,  2.634450,  2.998690,  3.347180,  3.850210,
     O -3.247590, -2.813320, -1.944900, -1.077620, -0.645870, -0.218270,
     1  0.199320,  0.596150,  0.958880,  1.279060,  1.556640,  1.796800,
     2  2.189980,  2.500250,  2.865750,  3.237410,  3.591560,  4.100720,
     3 -3.114580, -2.680300, -1.811880, -0.944500, -0.512580, -0.084530,
     4  0.334140,  0.733170,  1.099480,  1.424140,  1.706340,  1.950750,
     5  2.350700,  2.665780,  3.036110,  3.411730,  3.768870,  4.281210/
      DATA EQ4  /                                                       
     1 0.00001, 0.39795,  11.50006    , 0.00001, 0.00001, -2.302580    ,
     2-0.00002,-0.00002,  3.453870    , 0.00002, 0.39796,  10.50015    ,
     3 0.00005, 0.00005, -2.302530    ,-0.00008,-0.00008,  3.453800    ,
     4 0.00019, 0.39813,  8.501100    , 0.00043, 0.00043, -2.302100    ,
     5-0.00064,-0.00064,  3.453150    , 0.00139, 0.39933,  6.508040    ,
     6 0.00319, 0.00319, -2.298900    ,-0.00479,-0.00479,  3.448360    ,
     7 0.00373, 0.40167,  5.521610    , 0.00855, 0.00855, -2.292680    ,
     8-0.01283,-0.01283,  3.439020    , 0.00979, 0.40773,  4.557010    ,
     9 0.02223, 0.02223, -2.276540    ,-0.03335,-0.03335,  3.414810    ,
     O 0.02432, 0.42226,  3.643990    , 0.05412, 0.05412, -2.237490    ,
     1-0.08118,-0.08118,  3.356240    , 0.05437, 0.45231,  2.833440    ,
     2 0.11605, 0.11605, -2.155380    ,-0.17408,-0.17408,  3.233080    ,
     3 0.10406, 0.50200,  2.176890    , 0.20725, 0.20725, -2.016380    ,
     4-0.31087,-0.31087,  3.024570    , 0.16883, 0.56677,  1.687870    ,
     5 0.30724, 0.30724, -1.831260    ,-0.46086,-0.46086,  2.746890    ,
     6 0.23945, 0.63739,  1.339010    , 0.39498, 0.39498, -1.628630    ,
     7-0.59247,-0.59247,  2.442950    , 0.30887, 0.70681,  1.091120    ,
     8 0.46280, 0.46280, -1.433880    ,-0.69420,-0.69420,  2.150820    /
      DATA EQ5  /                                                       
     1 0.43322, 0.83116, 0.7790200    , 0.54788, 0.54788, -1.112400    ,
     2-0.82183,-0.82183,  1.668600    , 0.53649, 0.93443, 0.5988300    ,
     3 0.59215, 0.59215,-0.8851500    ,-0.88823,-0.88823,  1.327730    ,
     4 0.66055, 1.05849, 0.4418100    , 0.62453, 0.62453,-0.6660100    ,
     5-0.93679,-0.93679, 0.9990200    , 0.78749, 1.18543, 0.3263200    ,
     6 0.64322, 0.64322,-0.4963200    ,-0.96483,-0.96483, 0.7444900    ,
     7 0.90839, 1.30633, 0.2456200    , 0.65325, 0.65325,-0.3751300    ,
     8-0.97988,-0.97988, 0.5627000    , 1.08155, 1.47949, 0.1641600    ,
     9 0.66063, 0.66063,-0.2514400    ,-0.99094,-0.99094, 0.3771600    ,
     O 0.00001, 0.41317,  11.58926    , 0.00001, 0.00001, -2.302570    ,
     1-0.00002, 0.03272,  3.854510    , 0.00002, 0.41319,  10.58936    ,
     2 0.00005, 0.00005, -2.302520    ,-0.00008, 0.03266,  3.854450    ,
     3 0.00018, 0.41335,  8.590310    , 0.00041, 0.00042, -2.301990    ,
     4-0.00067, 0.03206,  3.853680    , 0.00133, 0.41452,  6.597300    ,
     5 0.00305, 0.00310, -2.298130    ,-0.00499, 0.02772,  3.848150    ,
     6 0.00357, 0.41679,  5.610950    , 0.00817, 0.00831, -2.290610    ,
     7-0.01335, 0.01929,  3.837330    , 0.00935, 0.42268,  4.646580    ,
     8 0.02120, 0.02157, -2.271150    ,-0.03468,-0.00220,  3.809380    /
      DATA EQ6  /                                                       
     1 0.02324, 0.43681,  3.734110    , 0.05145, 0.05240, -2.224260    ,
     2-0.08434,-0.05230,  3.741930    , 0.05194, 0.46607,  2.924660    ,
     3 0.10953, 0.11177, -2.126680    ,-0.18051,-0.14946,  3.601070    ,
     4 0.09934, 0.51452,  2.269850    , 0.19313, 0.19776, -1.964440    ,
     5-0.32129,-0.29210,  3.365350    , 0.16101, 0.57782,  1.782870    ,
     6 0.28115, 0.28929, -1.754110    ,-0.47418,-0.44768,  3.056730    ,
     7 0.22802, 0.64697,  1.435810    , 0.35366, 0.36598, -1.531210    ,
     8-0.60663,-0.58324,  2.725400    , 0.29366, 0.71509,  1.189150    ,
     9 0.40485, 0.42153, -1.324010    ,-0.70757,-0.68725,  2.412780    ,
     O 0.41045, 0.83743, 0.8775700    , 0.45768, 0.48242,-0.9965500    ,
     1-0.83179,-0.81664,  1.907080    , 0.50657, 0.93932, 0.6960500    ,
     2 0.47416, 0.50551,-0.7755900    ,-0.89508,-0.88362,  1.554530    ,
     3 0.62084, 1.06203, 0.5353900    , 0.47325, 0.51197,-0.5704300    ,
     4-0.94068,-0.93281,  1.214140    , 0.73633, 1.18792, 0.4143000    ,
     5 0.45891, 0.50399,-0.4173400    ,-0.96673,-0.96157, 0.9462800    ,
     6 0.84500, 1.30807, 0.3268900    , 0.43931, 0.48888,-0.3117000    ,
     7-0.98065,-0.97733, 0.7491100    , 0.99876, 1.48049, 0.2342500    ,
     8 0.40919, 0.46178,-0.2079600    ,-0.99094,-0.98928, 0.5381200    /
      DATA EQ7  /                                                       
     1 0.09396, 0.75814,  14.61590    ,-0.35121,-0.96418, -14.57590    ,
     2 4.09866, 2.43478,  30.74150    , 0.01537, 0.49499,  11.10817    ,
     3-0.06822,-0.28662, -4.285940    , 0.79659, 0.97928,  9.960740    ,
     4 0.00046, 0.42805,  8.679160    ,-0.00094,-0.00594, -2.339530    ,
     5 0.01495, 0.07821,  4.312840    , 0.00128, 0.42783,  6.678150    ,
     6 0.00290, 0.00290, -2.298210    ,-0.00484, 0.05267,  4.193280    ,
     7 0.00342, 0.43002,  5.691700    , 0.00782, 0.00805, -2.289060    ,
     8-0.01368, 0.04358,  4.178950    , 0.00897, 0.43574,  4.727380    ,
     9 0.02030, 0.02093, -2.266910    ,-0.03563, 0.02124,  4.146510    ,
     O 0.02228, 0.44948,  3.815040    , 0.04912, 0.05074, -2.213980    ,
     1-0.08657,-0.03066,  4.069120    , 0.04981, 0.47796,  3.005830    ,
     2 0.10399, 0.10778, -2.104950    ,-0.18489,-0.13112,  3.908690    ,
     3 0.09525, 0.52520,  2.351220    , 0.18163, 0.18934, -1.927040    ,
     4-0.32792,-0.27806,  3.643730    , 0.15435, 0.58704,  1.864060    ,
     5 0.26113, 0.27432, -1.702400    ,-0.48180,-0.43746,  3.302670    ,
     6 0.21854, 0.65476,  1.516160    , 0.32406, 0.34342, -1.471510    ,
     7-0.61362,-0.57545,  2.943250    , 0.28136, 0.72163,  1.267920    ,
     8 0.36625, 0.39160, -1.263190    ,-0.71299,-0.68072,  2.609510    /
      DATA EQ8  /                                                       
     1 0.39310, 0.84208, 0.9515400    , 0.40571, 0.44087,-0.9454200    ,
     2-0.83392,-0.81115,  2.077260    , 0.48510, 0.94271, 0.7642400    ,
     3 0.41490, 0.45674,-0.7378600    ,-0.89512,-0.87879,  1.708000    ,
     4 0.59475, 1.06427, 0.5950100    , 0.41025, 0.45783,-0.5487900    ,
     5-0.93940,-0.92894,  1.348150    , 0.70620, 1.18932, 0.4641400    ,
     6 0.39731, 0.44790,-0.4085900    ,-0.96509,-0.95878,  1.058920    ,
     7 0.81200, 1.30893, 0.3672900    , 0.38301, 0.43381,-0.3111600    ,
     8-0.97916,-0.97543, 0.8409600    , 0.96390, 1.48089, 0.2621800    ,
     9 0.36468, 0.41192,-0.2131400    ,-0.98993,-0.98830, 0.6020000    ,
     O 2.03780, 2.92758,  851.4021    ,-0.99991,-1.00582, -1960.300    ,
     1 6.80374,-0.12065, -264.0530    , 1.60365, 2.48900,  313.2954    ,
     2-0.99937,-1.01558, -721.0559    , 6.80008,-0.02069, -25.81689    ,
     3 0.74246, 1.58127,  42.94679    ,-0.96715,-1.09404, -96.47340    ,
     4 6.58107, 0.98494,  77.80141    , 0.09548, 0.67325,  8.597500    ,
     5-0.35035,-0.72538, -9.734900    , 2.39338, 2.23562,  25.37520    ,
     6 0.01874, 0.49369,  6.098310    ,-0.06241,-0.17424, -3.504520    ,
     7 0.45770, 0.69160,  8.635050    , 0.01055, 0.46342,  4.904170    ,
     8 0.00839,-0.00870, -2.437260    , 0.03307, 0.15849,  5.316780    /
      DATA EQ9  /                                                       
     1 0.02100, 0.47121,  3.959030    , 0.04373, 0.04349, -2.226000    ,
     2-0.07890, 0.01645,  4.700500    , 0.04631, 0.49749,  3.144000    ,
     3 0.09488, 0.10007, -2.080270    ,-0.18793,-0.10321,  4.412010    ,
     4 0.08856, 0.54235,  2.486230    , 0.16412, 0.17513, -1.880980    ,
     5-0.33380,-0.25838,  4.072200    , 0.14369, 0.60150,  1.994840    ,
     6 0.23268, 0.25092, -1.643430    ,-0.48732,-0.42268,  3.666950    ,
     7 0.20377, 0.66663,  1.641210    , 0.28503, 0.31065, -1.410550    ,
     8-0.61681,-0.56328,  3.256510    , 0.26279, 0.73127,  1.386120    ,
     9 0.31902, 0.35105, -1.208680    ,-0.71338,-0.66988,  2.885440    ,
     O 0.36851, 0.84848,  1.054810    , 0.35026, 0.39098,-0.9122300    ,
     1-0.83041,-0.80195,  2.304910    , 0.45647, 0.94710, 0.8533200    ,
     2 0.35861, 0.40354,-0.7227900    ,-0.89026,-0.87118,  1.903130    ,
     3 0.56269, 1.06693, 0.6664600    , 0.35833, 0.40490,-0.5495000    ,
     4-0.93467,-0.92346,  1.506430    , 0.67241, 1.19085, 0.5185200    ,
     5 0.35348, 0.39852,-0.4177200    ,-0.96142,-0.95526,  1.180870    ,
     6 0.77826, 1.30978, 0.4074800    , 0.34813, 0.38959,-0.3229100    ,
     7-0.97663,-0.97329, 0.9319300    , 0.93269, 1.48124, 0.2865700    ,
     8 0.34186, 0.37624,-0.2236600    ,-0.98865,-0.98736, 0.6576600    /
      DATA EQ10 /                                                       
     1 2.75988, 3.55482,  3591.100    ,-0.99999,-1.00092, -8268.801    ,
     2 5.24227, 1.25493,  10359.60    , 2.32562, 3.11986,  1321.190    ,
     3-0.99997,-1.00250, -3042.100    , 5.24217, 1.26996,  3846.080    ,
     4 1.45755, 2.24473,  178.9840    ,-0.99876,-1.01772, -411.7839    ,
     5 5.23585, 1.41871,  566.6631    , 0.60468, 1.33808,  24.85699    ,
     6-0.93724,-1.08120, -54.80800    , 4.91480, 2.29894,  110.6460    ,
     7 0.24636, 0.88412,  10.33480    ,-0.67073,-0.96638, -18.37801    ,
     8 3.52599, 2.79883,  48.60970    , 0.06262, 0.58319,  5.762500    ,
     9-0.20920,-0.40729, -5.695200    , 1.14125, 1.51452,  16.47990    ,
     O 0.02816, 0.50479,  4.186320    , 0.00096,-0.04019, -2.755120    ,
     1 0.12141, 0.35044,  7.187140    , 0.04480, 0.51518,  3.272910    ,
     2 0.08151, 0.08063, -2.152930    ,-0.15664,-0.02404,  5.206270    ,
     3 0.08356, 0.55585,  2.595930    , 0.15036, 0.16155, -1.873020    ,
     4-0.32893,-0.23195,  4.487030    , 0.13562, 0.61232,  2.095610    ,
     5 0.21341, 0.23313, -1.619450    ,-0.48532,-0.40912,  3.952860    ,
     6 0.19286, 0.67520,  1.733730    , 0.26081, 0.28795, -1.387470    ,
     7-0.61350,-0.55340,  3.480340    , 0.24950, 0.73800,  1.470260    ,
     8 0.29186, 0.32484, -1.192070    ,-0.70853,-0.66137,  3.073140    /
      DATA EQ11 /                                                       
     1 0.35199, 0.85267,  1.123250    , 0.32228, 0.36204,-0.9094500    ,
     2-0.82438,-0.79512,  2.450800    , 0.43833, 0.94982, 0.9088200    ,
     3 0.33305, 0.37503,-0.7284000    ,-0.88470,-0.86590,  2.022240    ,
     4 0.54381, 1.06847, 0.7077200    , 0.33750, 0.37884,-0.5598800    ,
     5-0.93048,-0.92000,  1.596650    , 0.65400, 1.19166, 0.5475800    ,
     6 0.33790, 0.37597,-0.4286200    ,-0.95870,-0.95323,  1.245660    ,
     7 0.76116, 1.31021, 0.4274900    , 0.33701, 0.37058,-0.3322400    ,
     8-0.97499,-0.97216, 0.9770700    , 0.91834, 1.48140, 0.2976300    ,
     9 0.33553, 0.36195,-0.2298800    ,-0.98794,-0.98691, 0.6828800    ,
     O 3.33390, 4.03917,  10946.10    ,-0.99999,-1.00018, -25204.30    ,
     1 4.07405, 2.25026,  56701.60    , 2.89963, 3.60477,  4027.100    ,
     2-0.99999,-1.00049, -9272.602    , 4.07405, 2.25305,  20876.90    ,
     3 2.03129, 2.73505,  545.3101    ,-0.99990,-1.00359, -1255.570    ,
     4 4.07369, 2.28101,  2849.150    , 1.16548, 1.85896,  74.20000    ,
     5-0.99510,-1.02271, -170.3580    , 4.05429, 2.47211,  407.4460    ,
     6 0.74115, 1.41414,  27.75800    ,-0.96536,-1.04034, -62.60300    ,
     7 3.93433, 2.76545,  162.3960    , 0.35754, 0.97944,  11.03890    ,
     8-0.79471,-0.96528, -22.27600    , 3.24697, 2.99888,  65.54620    /
      DATA EQ12 /                                                       
     1 0.11662, 0.66103,  5.451400    ,-0.34137,-0.50725, -7.335200    ,
     2 1.43542, 1.93517,  23.63071    , 0.05761, 0.56285,  3.629470    ,
     3-0.00921,-0.05262, -3.047810    , 0.15920, 0.48574,  9.455220    ,
     4 0.07890, 0.57874,  2.785970    , 0.11802, 0.12123, -2.029010    ,
     5-0.26520,-0.10302,  5.800570    , 0.12469, 0.62794,  2.244560    ,
     6 0.18672, 0.20468, -1.640310    ,-0.46541,-0.36843,  4.533860    ,
     7 0.17800, 0.68682,  1.861810    , 0.23233, 0.25791, -1.387660    ,
     8-0.59999,-0.53395,  3.829840    , 0.23187, 0.74678,  1.581840    ,
     9 0.26252, 0.29281, -1.194750    ,-0.69577,-0.64801,  3.329250    ,
     O 0.33131, 0.85782,  1.208120    , 0.29539, 0.32978,-0.9227700    ,
     1-0.81311,-0.78600,  2.628910    , 0.41672, 0.95299, 0.9741200    ,
     2 0.31051, 0.34509,-0.7464400    ,-0.87581,-0.85933,  2.160280    ,
     3 0.52264, 1.07016, 0.7534300    , 0.32078, 0.35301,-0.5780600    ,
     4-0.92464,-0.91597,  1.695870    , 0.63452, 1.19251, 0.5780000    ,
     5 0.32644, 0.35466,-0.4436300    ,-0.95530,-0.95102,  1.313150    ,
     6 0.74396, 1.31063, 0.4474300    , 0.32940, 0.35326,-0.3434300    ,
     7-0.97310,-0.97099,  1.021970    , 0.90476, 1.48155, 0.3080400    ,
     8 0.33155, 0.34944,-0.2364300    ,-0.98721,-0.98648, 0.7065700    /
      DATA EQ13 /                                                       
     1 3.61829, 4.27758,  18950.00    ,-0.99999,-1.00006, -43634.00    ,
     2 3.56376, 2.64102,  115228.0    , 3.18401, 3.84325,  6971.699    ,
     3-0.99999,-1.00018, -16053.00    , 3.56375, 2.64212,  42402.90    ,
     4 2.31564, 2.97435,  943.9900    ,-0.99997,-1.00138, -2173.590    ,
     5 3.56367, 2.65323,  5755.969    , 1.44891, 2.10369,  128.2630    ,
     6-0.99862,-1.00913, -295.1160    , 3.55892, 2.73076,  795.5969    ,
     7 1.01954, 1.66652,  47.63100    ,-0.99009,-1.01870, -109.0800    ,
     8 3.52881, 2.86095,  303.4351    , 0.60608, 1.23252,  18.14000    ,
     9-0.93307,-1.00707, -40.23199    , 3.32782, 3.05609,  119.3000    ,
     O 0.26334, 0.84562,  7.659700    ,-0.66442,-0.79734, -14.28420    ,
     1 2.38292, 2.74914,  46.21851    , 0.09668, 0.63353,  4.216100    ,
     2-0.19131,-0.26944, -4.994600    , 0.73618, 1.21168,  16.90700    ,
     3 0.08078, 0.60386,  2.999320    , 0.06394, 0.05503, -2.441620    ,
     4-0.11959, 0.11554,  7.875830    , 0.11700, 0.64141,  2.375680    ,
     5 0.16176, 0.17552, -1.732220    ,-0.42537,-0.30857,  5.240120    ,
     6 0.16666, 0.69590,  1.963900    , 0.21250, 0.23444, -1.417810    ,
     7-0.58049,-0.51189,  4.162080    , 0.21862, 0.75331,  1.666260    ,
     8 0.24442, 0.27024, -1.214340    ,-0.68121,-0.63587,  3.536550    /
      DATA EQ14 /                                                       
     1 0.31658, 0.86142,  1.268170    , 0.28075, 0.30900,-0.9416500    ,
     2-0.80271,-0.77914,  2.754410    , 0.40202, 0.95511,  1.018140    ,
     3 0.29921, 0.32675,-0.7644300    ,-0.86844,-0.85471,  2.252340    ,
     4 0.50897, 1.07124, 0.7826700    , 0.31311, 0.33792,-0.5928200    ,
     5-0.92024,-0.91331,  1.758940    , 0.62254, 1.19302, 0.5965500    ,
     6 0.32158, 0.34266,-0.4543100    ,-0.95294,-0.94964,  1.354190    ,
     7 0.73380, 1.31088, 0.4591400    , 0.32637, 0.34378,-0.3507600    ,
     8-0.97188,-0.97030,  1.048270    , 0.89711, 1.48164, 0.3138700    ,
     9 0.33008, 0.34278,-0.2403700    ,-0.98676,-0.98624, 0.7198400    ,
     O 3.75974, 4.39785,  24996.00    ,-1.00000,-1.00003, -57555.00    ,
     1 3.34636, 2.79204,  160689.0    , 3.32546, 3.96354,  9195.898    ,
     2-0.99999,-1.00010, -21174.40    , 3.34636, 2.79270,  59125.20    ,
     3 2.45708, 3.09485,  1245.110    ,-0.99998,-1.00078, -2866.940    ,
     4 3.34632, 2.79930,  8016.801    , 1.59009, 2.22564,  169.1000    ,
     5-0.99927,-1.00522, -389.2200    , 3.34395, 2.84562,  1099.550    ,
     6 1.15963, 1.79077,  62.69200    ,-0.99470,-1.01087, -143.9520    ,
     7 3.32879, 2.92499,  414.2939    , 0.73973, 1.35916,  23.66100    ,
     8-0.96322,-1.00587, -53.42200    , 3.22440, 3.05741,  160.1920    /
      DATA EQ15 /                                                       
     1 0.36626, 0.95836,  9.568000    ,-0.78976,-0.87983, -19.53200    ,
     2 2.64973, 2.92854,  62.65640    , 0.13746, 0.69241,  4.778600    ,
     3-0.33875,-0.41825, -6.832700    , 1.16229, 1.65412,  23.26630    ,
     4 0.08614, 0.62462,  3.181550    , 0.00990,-0.00588, -2.892700    ,
     5 0.03025, 0.30203,  9.738980    , 0.11270, 0.65103,  2.470700    ,
     6 0.14194, 0.15228, -1.838540    ,-0.38404,-0.25761,  5.815300    ,
     7 0.15946, 0.70187,  2.031980    , 0.20013, 0.21886, -1.452770    ,
     8-0.56296,-0.49493,  4.401040    , 0.21024, 0.75744,  1.720290    ,
     9 0.23434, 0.25642, -1.234290    ,-0.66965,-0.62747,  3.672750    ,
     O 0.30755, 0.86361,  1.304860    , 0.27335, 0.29704,-0.9565200    ,
     1-0.79539,-0.77484,  2.830650    , 0.39328, 0.95637,  1.044210    ,
     2 0.29381, 0.31650,-0.7769800    ,-0.86356,-0.85193,  2.306440    ,
     3 0.50111, 1.07186, 0.7994200    , 0.30965, 0.32972,-0.6022000    ,
     4-0.91750,-0.91177,  1.794880    , 0.61585, 1.19331, 0.6068900    ,
     5 0.31950, 0.33627,-0.4607400    ,-0.95153,-0.94886,  1.376990    ,
     6 0.72826, 1.31102, 0.4655100    , 0.32513, 0.33880,-0.3549800    ,
     7-0.97117,-0.96991,  1.062580    , 0.89305, 1.48168, 0.3169700    ,
     8 0.32951, 0.33934,-0.2425400    ,-0.98652,-0.98611, 0.7268900    /
      DATA EQ16 /                                                       
     1 3.90098, 4.52093,  33185.00    ,-1.00000,-1.00001, -76412.00    ,
     2 3.16831, 2.90551,  222011.0    , 3.46670, 4.08664,  12208.70    ,
     3-0.99999,-1.00005, -28111.60    , 3.16831, 2.90587,  81682.69    ,
     4 2.59829, 3.21808,  1652.900    ,-0.99999,-1.00039, -3806.000    ,
     5 3.16829, 2.90946,  11067.60    , 1.73109, 2.34976,  224.4000    ,
     6-0.99961,-1.00259, -516.5901    , 3.16709, 2.93476,  1510.620    ,
     7 1.29988, 1.91634,  83.08900    ,-0.99716,-1.00525, -191.0420    ,
     8 3.15940, 2.97837,  564.5090    , 0.87604, 1.48662,  31.17900    ,
     9-0.98001,-1.00150, -71.05800    , 3.10537, 3.05371,  215.1590    ,
     O 0.48258, 1.07879,  12.26670    ,-0.87580,-0.92564, -26.42270    ,
     1 2.77722, 2.99133,  83.76520    , 0.19748, 0.76963,  5.643500    ,
     2-0.50520,-0.56701, -9.476700    , 1.61191, 2.04367,  31.75650    ,
     3 0.09757, 0.65435,  3.455300    ,-0.07362,-0.09304, -3.637100    ,
     4 0.26172, 0.54536,  12.41760    , 0.10841, 0.66374,  2.598200    ,
     5 0.11332, 0.11956, -2.021130    ,-0.31615,-0.18775,  6.601610    ,
     6 0.15091, 0.70925,  2.117010    , 0.18528, 0.19948, -1.511570    ,
     7-0.53622,-0.47297,  4.701430    , 0.20026, 0.76238,  1.785380    ,
     8 0.22347, 0.24035, -1.265380    ,-0.65346,-0.61718,  3.834810    /
      DATA EQ17 /                                                       
     1 0.29710, 0.86613,  1.347360    , 0.26606, 0.28380,-0.9765200    ,
     2-0.78602,-0.76981,  2.917680    , 0.38338, 0.95778,  1.073680    ,
     3 0.28874, 0.30543,-0.7927200    ,-0.85759,-0.84875,  2.367300    ,
     4 0.49242, 1.07253, 0.8178900    , 0.30657, 0.32103,-0.6133000    ,
     5-0.91428,-0.91005,  1.834340    , 0.60861, 1.19362, 0.6180400    ,
     6 0.31773, 0.32961,-0.4680000    ,-0.94994,-0.94801,  1.401550    ,
     7 0.72237, 1.31116, 0.4722800    , 0.32411, 0.33366,-0.3596100    ,
     8-0.97040,-0.96950,  1.077750    , 0.88881, 1.48173, 0.3202000    ,
     9 0.32905, 0.33582,-0.2448500    ,-0.98626,-0.98597, 0.7342300    ,
     O 3.97171, 4.58437,  38404.00    ,-1.00000,-1.00001, -88428.00    ,
     1 3.09887, 2.94643,  260544.0    , 3.53743, 4.15008,  14128.50    ,
     2-0.99999,-1.00003, -32532.20    , 3.09887, 2.94669,  95857.69    ,
     3 2.66900, 3.28156,  1912.800    ,-0.99999,-1.00025, -4404.398    ,
     4 3.09886, 2.94918,  12985.00    , 1.80169, 2.41352,  259.6101    ,
     5-0.99971,-1.00164, -597.7000    , 3.09800, 2.96673,  1769.150    ,
     6 1.37016, 1.98055,  96.07401    ,-0.99793,-1.00316, -220.9870    ,
     7 3.09250, 2.99684,  659.0339    , 0.94486, 1.55145,  35.97099    ,
     8-0.98532,-0.99925, -82.21001    , 3.05358, 3.04794,  249.6140    /
      DATA EQ18 /                                                       
     1 0.54424, 1.14138,  14.01120    ,-0.90606,-0.93934, -30.70900    ,
     2 2.80894, 2.99344,  96.62511    , 0.23494, 0.81494,  6.230400    ,
     3-0.58637,-0.63370, -11.16490    , 1.82290, 2.19435,  36.86819    ,
     4 0.10659, 0.67362,  3.642100    ,-0.12852,-0.14705, -4.161900    ,
     5 0.41410, 0.68303,  14.11630    , 0.10642, 0.67165,  2.679000    ,
     6 0.09434, 0.09851, -2.154500    ,-0.26767,-0.14542,  7.084520    ,
     7 0.14602, 0.71362,  2.167870    , 0.17663, 0.18797, -1.553780    ,
     8-0.51788,-0.46019,  4.874180    , 0.19454, 0.76522,  1.823150    ,
     9 0.21772, 0.23129, -1.286610    ,-0.64298,-0.61137,  3.924930    ,
     O 0.29123, 0.86754,  1.371230    , 0.26251, 0.27664,-0.9889700    ,
     1-0.78034,-0.76700,  2.965570    , 0.37792, 0.95855,  1.089900    ,
     2 0.28639, 0.29954,-0.8020000    ,-0.85410,-0.84699,  2.400320    ,
     3 0.48773, 1.07290, 0.8278600    , 0.30521, 0.31649,-0.6195700    ,
     4-0.91247,-0.90912,  1.855490    , 0.60477, 1.19379, 0.6239600    ,
     5 0.31698, 0.32616,-0.4720000    ,-0.94907,-0.94755,  1.414620    ,
     6 0.71928, 1.31124, 0.4758200    , 0.32370, 0.33103,-0.3621000    ,
     7-0.96998,-0.96929,  1.085690    , 0.88662, 1.48175, 0.3218600    ,
     8 0.32888, 0.33403,-0.2460700    ,-0.98612,-0.98590, 0.7380200    /
      DATA EQ19 /                                                       
     1  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     2  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     3  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     4  0.021189,  0.021189,  0.021189,  0.021189,  0.021189,  0.021189,
     5  0.021189,  0.021189,  0.021189,  0.021189,  0.031408,  0.041392,
     6  0.060698,  0.079181,  0.105510,  0.138303,  0.176091,  0.243038,
     7  0.041392,  0.041392,  0.041392,  0.041392,  0.041392,  0.041392,
     8  0.041392,  0.041392,  0.041392,  0.041392,  0.060698,  0.079181,
     9  0.113943,  0.146128,  0.190331,  0.243038,  0.301030,  0.397940,
     O  0.079181,  0.079181,  0.079181,  0.079181,  0.079181,  0.079181,
     1  0.079181,  0.079181,  0.079181,  0.079181,  0.113943,  0.146128,
     2  0.204120,  0.255272,  0.322219,  0.397940,  0.477121,  0.602060,
     3  0.113943,  0.113943,  0.113943,  0.113943,  0.113943,  0.113943,
     4  0.113943,  0.113943,  0.113943,  0.113943,  0.161368,  0.204120,
     5  0.278753,  0.342423,  0.423246,  0.511883,  0.602060,  0.740362/
      DATA EQ20 /                                                       
     1  0.176091,  0.176091,  0.176091,  0.176091,  0.176091,  0.176091,
     2  0.176091,  0.176091,  0.176091,  0.176091,  0.243038,  0.301030,
     3  0.397940,  0.477121,  0.574031,  0.676693,  0.778151,  0.929419,
     4  0.243038,  0.243038,  0.243038,  0.243038,  0.243038,  0.243038,
     5  0.243038,  0.243038,  0.243038,  0.243038,  0.327359,  0.397940,
     6  0.511883,  0.602060,  0.709694,  0.821186,  0.929419,  1.088136,
     7  0.301030,  0.301030,  0.301030,  0.301030,  0.301030,  0.301030,
     8  0.301030,  0.301030,  0.301030,  0.301030,  0.397940,  0.477121,
     9  0.602060,  0.698970,  0.812913,  0.929419,  1.041392,  1.204120,
     O  0.397940,  0.397940,  0.397940,  0.397940,  0.397940,  0.397940,
     1  0.397940,  0.397940,  0.397940,  0.397940,  0.511883,  0.602060,
     2  0.740362,  0.845098,  0.966141,  1.088136,  1.204120,  1.371067,
     3  0.477121,  0.477121,  0.477121,  0.477121,  0.477121,  0.477121,
     4  0.477121,  0.477121,  0.477121,  0.477121,  0.602060,  0.698970,
     5  0.845098,  0.954242,  1.079181,  1.204120,  1.322219,  1.491361/
      DATA EQ21 /                                                       
     1  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     2  0.000000,  0.000000,  0.000000,  0.666667,  0.666667,  0.666667,
     3  0.666667,  0.666667,  0.666667,  0.666667,  0.666667,  0.666667,
     4  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     5  0.000000,  0.000000,  0.000000,  0.623886,  0.605621,  0.589062,
     6  0.560297,  0.536239,  0.506850,  0.476762,  0.449094,  0.413365,
     7  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     8  0.000000,  0.000000,  0.000000,  0.589063,  0.560304,  0.536237,
     9  0.498485,  0.470478,  0.440242,  0.413365,  0.391959,  0.368634,
     O  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     1  0.000000,  0.000000,  0.000000,  0.536239,  0.498485,  0.470479,
     2  0.432372,  0.408235,  0.385738,  0.368635,  0.356867,  0.345899,
     3  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     4  0.000000,  0.000000,  0.000000,  0.498485,  0.459103,  0.432372,
     5  0.399351,  0.380458,  0.364320,  0.353071,  0.345899,  0.339687/
      DATA EQ22 /                                                       
     1  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     2  0.000000,  0.000000,  0.000000,  0.449094,  0.413365,  0.391959,
     3  0.368634,  0.356868,  0.347780,  0.342017,  0.338617,  0.335876,
     4  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     5  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     6  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     7  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     8  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     9  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     O  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     1  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     2  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     3  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     4  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     5  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000/
      DATA EQ23 /                                                       
     1  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     2  0.000000,  0.000000,  0.000000,  1.000000,  1.000000,  1.000000,
     3  1.000000,  1.000000,  1.000000,  1.000000,  1.000000,  1.000000,
     4  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     5  0.000000,  0.000000,  0.000000,  0.954545,  0.934783,  0.916666,
     6  0.884616,  0.857143,  0.822581,  0.785714,  0.750000,  0.700000,
     7  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     8  0.000000,  0.000000,  0.000000,  0.916666,  0.884616,  0.857143,
     9  0.812500,  0.777777,  0.738095,  0.700000,  0.666667,  0.625000,
     O  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     1  0.000000,  0.000000,  0.000000,  0.857143,  0.812500,  0.777778,
     2  0.727273,  0.692308,  0.656250,  0.625000,  0.600000,  0.571428,
     3  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     4  0.000000,  0.000000,  0.000000,  0.812500,  0.763158,  0.727273,
     5  0.678571,  0.647059,  0.616279,  0.590909,  0.571429,  0.550000/
      DATA EQ24 /                                                       
     1  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     2  0.000000,  0.000000,  0.000000,  0.750000,  0.700000,  0.666667,
     3  0.625000,  0.600000,  0.576923,  0.558824,  0.545455,  0.531250,
     4  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     5  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     6  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     7  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     8  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     9  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     O  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     1  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     2  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     3  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     4  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     5  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000/
      DATA EQ25 /                                                       
     1  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     2  0.000000,  0.000000,  0.000000,  0.500000,  0.500000,  0.500000,
     3  0.500000,  0.500000,  0.500000,  0.500000,  0.500000,  0.500000,
     4  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     5  0.000000,  0.000000,  0.000000,  0.522624,  0.532293,  0.540983,
     6  0.555762,  0.567567,  0.580823,  0.592308,  0.600000,  0.603448,
     7  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     8  0.000000,  0.000000,  0.000000,  0.540983,  0.555762,  0.567567,
     9  0.584269,  0.594340,  0.601664,  0.603448,  0.600000,  0.588235,
     O  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     1  0.000000,  0.000000,  0.000000,  0.567567,  0.584269,  0.594340,
     2  0.602740,  0.603093,  0.597865,  0.588235,  0.576923,  0.560000,
     3  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     4  0.000000,  0.000000,  0.000000,  0.584269,  0.597614,  0.602740,
     5  0.601810,  0.595541,  0.584659,  0.572000,  0.560000,  0.544555/
      DATA EQ26 /                                                       
     1  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     2  0.000000,  0.000000,  0.000000,  0.600000,  0.603448,  0.600000,
     3  0.588235,  0.576923,  0.563584,  0.551195,  0.540984,  0.529183,
     4  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     5  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     6  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     7  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     8  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     9  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     O  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     1  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     2  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     3  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     4  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     5  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000/
      DATA KMJM / 10, 18/                                               
      END                                                               
