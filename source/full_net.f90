!*******************************************************************************
! Full_Net, part of XNet 7 5/25/10 
! The subroutines in this file perform nucleosynthesis for a single Lagrangian 
! mass zone.  The handling of multiple zones, including multi-processing, is 
! carried out externally.  
!*******************************************************************************
  
Subroutine full_net
!===============================================================================
! The abundance evolution is performed over a series of timesteps, with the 
! duration of the timestep determined by the integration scheme and the changing 
! thermodynamic conditions. Integration is performed by a choice of methods 
! controlled by the isolv flag. 
!===============================================================================
  Use controls
  Use nuclear_data
  Use abundances
  Use conditions
  Use thermo_data
  Use timers
  Real(8) :: xtot,enm,enb,enold,en0,edot
  Real(8) :: ytot,ztot,atot
  Integer :: kstep,idiag0,its
  Integer :: i,j,k 
  
! Initiate timer
  start_timer = xnet_wtime()
  timer_total = timer_total - start_timer  

! Set reaction controls not read in from control
! iweak0=iweak
  idiag0=idiag
  kstep=0
  
! Normalize initial abundances and change units if necessary
  yt=y
! Call norm(yt) 

  If(iheat>0) Then
    t9t=t9
    t9o=t9
  EndIf
  
! Calculate the total energy of the nuclei
  Call benuc(yt,enb,enm,ytot,ztot,atot)
  xtot=atot-1.0
  en0=enm
  edot=0.0
  
! Start evolution 
  Write(6,"(a,i6,a,i2,2(a,es10.3))") & 
&    'Max Step',kstmx,' IDiag=',idiag,' Start Time',tstart,' Stop Time',tstop
  t=tstart
  tt=tstart
  Call t9rhofind(kstep,tt,j)

! Output initial abundances and conditions
  kmon=0
  ktot=0
  If(itsout>0) Call ts_output(kstep,(enm-en0),edot)

  Do kstep=1,kstmx
    If(idiag>=3) Write(lun_diag,*) 'KStep',kstep
!-----------------------------------------------------------------------------  
! For each step, an initial guess for the timestep is calculated, based on 
! successof the prior integration step and the abundance variations of the 
! previous step, If available.  
!-----------------------------------------------------------------------------  
    Call timestep(kstep)
    If(idiag>=1) Write(lun_diag,*) 'TDel',tt,tdel
  
! Determine If this is an output step
    idiag=idiag0
!   If(mod(kstep,10).eq.0) idiag=2
!   If(kstep==518.or.kstep==718.or.kstep==803) idiag=5
  
! Take Integration step
    Select Case (isolv)
      Case(2)
        Call solve_bd(kstep,its)
      Case Default
        Call solve_be(kstep,its)
    End Select
    
! If convergence is successful, update energy and other abundance dependent values 
    If(its==0) Then
      If(idiag>=2) Then
        Write(lun_diag,"(a,es23.16)") 'delta Y',tdel
        Write(lun_diag,"(a5,4es12.4)") &
&         (nname(k),y(k),yo(k),(y(k)-yo(k)),(tdel*ydot(k)),k=1,ny)
        If(iheat>0) Write(lun_diag,"(a,4es12.4)") 'delta T9',t9,t9o,t9-t9o,tdel*t9dot
      EndIf
      enold=enm
      Call benuc(yt,enb,enm,ytot,ztot,atot)
      edot=-(enm-enold)/tdel
      If(idiag>=1) Then 
        Write(lun_diag,"(i5,5es14.7)") kstep,t,tdel,t9t,rhot,yet
!       Write(lun_diag,"(5(a5,es11.4))") (nname(k),y(k),k=1,ny)
      EndIf
      If(itsout>0) Call ts_output(kstep,(enm-en0),edot)
  
! If Stop time is reached, exit
      If(t>=tstop) Exit 
      
! If reduced timesteps fail to successfully integrate, warn and exit
    Else
      Write(lun_diag,"(a,es12.4,a,es12.4,a)") 'Inter!!!' 
      Exit
    EndIf
  EndDo 
  
! Test that the stop time is reached
  If(t<tstop.or.its/=0) Then
    Write(6,"(a,es12.4,a,es12.4,a)") 'Evolution stopped at time=',t, &
&     ', stop time (',tstop,') not reached!' 
    Write(6,"(a,i12,a)") 'Approximately',int((tstop-t)/tdel),'more steps needed' 
  EndIf
  
! Stop timer
  stop_timer = xnet_wtime()
  timer_total = timer_total + stop_timer

! End Post Processing cycle
  iweak=iweak0
  Call final_output(kstep)
  Return
End Subroutine full_net
   
  
Subroutine yderiv
!-----------------------------------------------------------------------------  
! This routine calculates time derivatives for each nuclear species.
! This calculation is performed by looping over nuclei, and summing the 
! reaction rates for each reaction which effects that nucleus.
!-----------------------------------------------------------------------------  
  Use controls
  Use nuclear_data
  Use abundances
  Use conditions
  Use cross_sect_data
  Use reac_rate_data
  Use timers
  Integer :: i,j,i0,i1,la1,le1,la2,le2,la3,le3
  Real(8) :: s1,s11,s2,s22,s3,s33
  
! Initiate timer
  start_timer = xnet_wtime()
  timer_deriv = timer_deriv - start_timer  

! From the cross sections and the counting array, calculate the reaction 
! rates
  b1=a1*csect1(mu1)
  b2=a2*csect2(mu2)
  b3=a3*csect3(mu3)
  
! Calculate Ydot for each nucleus, summing over the reactions which affect it.
  Do i0=1,ny                 
    If(idiag>=5) Write(lun_diag,"(a3,a6,i4)") 'NUC',nname(i0),i0
  
! Sum over the reactions with 1 reactant
    la1=la(1,i0)                          
    le1=le(1,i0)                         
    s1=0.0                              
    Do i1=la1,le1                    
      s11=b1(i1)*yt(n11(i1))
      s1=s1+s11
      If(idiag>=5) Write(lun_diag,"(3a5,'  1  ',4es15.7)") &
&       (nname(n1i(j,mu1(i1))),j=1,3),b1(i1),yt(n11(i1)),s11,a1(i1)
    EndDo       
    If(idiag>=5) Write(lun_diag,*) '1->',nname(i0),la1,le1,s1
  
! Sum over the reactions with 2 reactants
    la2=la(2,i0)  
    le2=le(2,i0) 
    s2=0.0      
    Do i1=la2,le2           
      s22=b2(i1)*yt(n21(i1))*yt(n22(i1))
      s2=s2+s22                                       
      If(idiag>=5) Write(lun_diag,"(4a5,4es15.7)") &
&       (nname(n2i(i,mu2(i1))),i=1,4),b2(i1),yt(n21(i1)),yt(n22(i1)),s22
    EndDo                  
    If(idiag>=5) Write(lun_diag,*) '2->',nname(i0),la2,le2,s2
  
! Sum over the reactions with 3 reactants
    la3=la(3,i0)                          
    le3=le(3,i0)                         
    s3=0.0                              
    Do i1=la3,le3                  
      s33=b3(i1)*yt(n31(i1))*yt(n32(i1))*yt(n33(i1))
      s3=s3+s33
      If(idiag>=5) Write(lun_diag,"(3a5,'  3  ',5es12.4)") &
&       (nname(n3i(i,mu3(i1))),i=1,3),b3(i1),yt(n31(i1)),yt(n32(i1)),yt(n33(i1)),s33
    EndDo                                       
    If(idiag>=5) Write(lun_diag,*) '3->',nname(i0),la3,le3,s3
  
! Sum the 3 components of Ydot
    ydot(i0)=s1+s2+s3                         
    If(idiag>=3) Write(lun_diag,"(a4,a5,2es24.16,3es11.3)") 'YDOT',nname(i0),yt(i0),ydot(i0),s1,s2,s3
  EndDo                                 

  If(iheat>0) Then
    t9dot=-sum(mex*ydot)/cv
    If(idiag>=3) Write(lun_diag,"(5x,a5,2es24.16)") '   T9',t9t,t9dot
  EndIf

! Stop timer
  stop_timer = xnet_wtime()
  timer_deriv = timer_deriv + stop_timer

  Return                                  
End Subroutine yderiv                                   
  
