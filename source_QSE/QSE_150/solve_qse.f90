!*******************************************************************************
! Backward Euler solver, part of XNet 7 5/25/10 
!
! These routines in this file perform the Backward Euler time integration for the 
! thermonuclear reaction network.
!*******************************************************************************
  
Subroutine solve_qse(kstep,its)
!===============================================================================
! This routine handles the potential failure of an individual Backward Euler step. 
! Repeated Calls to the BE step integrator are made with successivley smaller 
! timesteps until success is achieved or ktsmx trials fail.
! The abundances are evolved using a Newton-Raphson iteration scheme to solve 
! the equation yt(i)-y(i)/tdel = ydot(i), where y(i) is the abundance at the 
! beginning of the iteration, yt(i) is the trial abundance at the End of the 
! timestep, ydot(i) is the time derivative of the trial abundance  calculated 
! from reaction rates and tdel is the timestep.  
!===============================================================================
  Use controls
  Use conditions
  Use abundances
  Use qse_abundances
  Use qse_data
  Use screen_qse 
  Use nuclear_data
  Integer, Intent(in) :: kstep  
  Integer, Intent(out) :: its 
  Integer, Parameter :: ktsmx=10
  Integer kts,inr,kit,j,i,dummy1
  Integer :: nsi,zsi,nfe,zfe,ncr,zcr,iout,ieos
!  Real(8) :: ene,ye,yps,y2e,v,pa,ea,dlma,pta,pva,eta,eva,emas,beta
  logical, save :: firstCall = .TRUE.
   if(firstCall)then
    write(*,*) "qse called" 
         firstcaller=1
      firstCall=.false.
   end if
 
   if(firstcaller>0) then
!    write(*,*), "firstcall top", firstcaller
         call group_number()
     if(.not.allocated(dlCqse)) Allocate(dlCqse(ny)) ! here is where we set the size of the QSE_coeff arra
     !   write(*,*) he(int(zz(43))), "hesi"
!Set up step 0 of QSE
!----------------------------------------------------------------------
! Build QSE coeff. from 4 focal
      call qse_coeff !
!repalce initial abund. with group members.
        call qse_make (yt)


!Solve for QSE this is step 0.
        call update(yt) ! y is getting screwed up
!        ye=sum(zz*yt)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Build yg. This is the set with the qse group abundances. Note: it is
! the same size as yr and has the same single nuc except for the focal
! abundances.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           yq =0.0

           do i=1,gp1i
              yq(1) =yq(1) + nn(i)*yt(i)
              yq(2) =yq(2)+  zz(i)*yt(i)
            enddo
            do ig=1,gp2i
               i=gp2number(ig)
               nsi=nn(i)-nn(f3)
               zsi=zz(i)-zz(f3)
               yq(1)=yq(1)+nsi*yt(i)
               yq(2)=yq(2)+zsi*yt(i)
               yq(3)=yq(3)+yt(i)
             enddo
             do ig=1,gp3i
               i=gp3number(ig)
               nfe=nn(i)-nn(f4)
               zfe=zz(i)-zz(f4)
               yq(1)=yq(1)+nfe*yt(i)
               yq(2)=yq(2)+zfe*yt(i)
               yq(4)=yq(4)+yt(i)
             enddo
             do ig=gn+1,gp0i
                i=gp0number(ig)
                yq(ig)=yt(i)
             enddo

                yqt=yq
            firstcaller=0
 endif ! firstCall
!-----------------------------------------------------------------------------  
! For each trial timestep, tdel, the Newton-Raphson iteration is attempted.
! The possible results for a timestep are a converged set of yt or a 
! failure to converge.  If convergence fails, iteration is retried with 
! the trial timestep reduced by tdel_maxmult, up to ktsmx times.
!-----------------------------------------------------------------------------  
 Do kts=1,ktsmx

  
! Attempt Backward Euler integration over desired timestep
!        write(300,*) kstep
        call qse_coeff 
        call qse_make(yt)
       Call step_qse(kstep,inr,kit)
    
! If integration fails, reset abundances, reduce timestep and retry.  
    If(inr/=0) Then
      tdel=tdel/tdel_maxmult
      tt=t+tdel

      If(idiag>=1) Write(lun_diag,"(a,i3,3es12.4)") 'BE TS Reduce',kts,tt,tdel,tdel_maxmult
! Reset abundances
         yqt=yq
          yt=y
  
! Set iweak to the original value before possibly changing it again in t9rhofind
      iweak=iweak0
      Call t9rhofind(kstep,tt,j)
      Call qse_coeff
      Call update2(yqt,yt)
  
    Else
      Exit
    EndIf
  EndDo
! If TS loop exits early, step is successful.  Prepare for next step.
  If(kts<=ktsmx) Then
!     write(*,*) "tag success"
    its=0
    If(idiag>=1) Write(lun_diag,"(a,i5,2i3)") 'TS Success',kstep,kts,kit
    to=t
    t=tt
    tdel_next=tdel*tdel_maxmult
!   Call norm(y)
    yo=y
    yqo=yq
    yq=yqt
    call update2(yq,yt)
    yo=y 
    y=yt
     call group_ab_test(kstep)
  If (iadpt==1)then
  call switcher(kstep)
  endif
  Else
    Write(6,"(i5,4es12.4,2i3)") kstep,t,tdel,t9t,rhot,kit,kts
    Write(6,*) 'Timestep retrys fail after ',kts,' attempts',tdel
    its=1
  EndIf
  kmon(1)=kts ; kmon(2) = kit
End Subroutine solve_qse
  
Subroutine step_qse(kstep,inr,kit)
!===============================================================================
! This routine attempts to integrate a single Backward Euler step for the 
! timestep tdel.  If successful, inr=0
!===============================================================================
  Use controls
  Use nuclear_data
  Use screen_qse 
  Use conditions
  Use abundances
  Use reac_rate_data
  Use qse_data
  Use qse_abundances
  Integer, Intent(in) :: kstep
  Integer, Intent(out) :: inr, kit
  Integer irdymx(1),idymx(1)
  Integer :: i,j
  Real(8) :: toln    ! Network convergence condition
  Real(8) :: testc,testc2,testm,testn  ! Convergence tests
  Real(8) :: xtot,xtot_init,rdt
  Integer :: k,idiag0
  Real(8), Dimension(gp0i) :: rhs,dy,reldy
  Integer :: nsi,zsi,nfe,zfe,ncr,zcr,iout,ieos
  Real(8) :: ene,ye,yps,y2e,v,pa,ea,dlma,pta,pva,eta,eva,emas,beta
  
! Calculate the thermodynamic factors necessary for reaction rates, 
! including screening, and the reaction rates.
           If(gn>2)then
  Call cross_sect 
           Endif
 Call en(ny,y,ye,yps,y2e)
  
  xtot_init=sum(aa*y)-1.0
  
! The Newton-Raphson iteration occurs for at most kitmx iterations.  
  Do kit=1,kitmx
               
! Calculate the reaction rates and abundance time derivatives
!write(*,*) "before qyder" 
   Call qyderiv
! write(*,*) "after qyder" 
! Build the Jacobian
    rdt=1.0/tdel
 !write(*,*) "before jacobian build"
    Call qjacobian_build(rdt,-1.0)
!write(*,*) "after jacobian build"   
! Calculate equation to zero
    rhs=((yq-yqt)*rdt+ydot)

    If(idiag>=1) Then
      Write(lun_diag,"(a3,i5,es14.7)") 'RHS',kstep,rdt
      Do i=1,gp0i
        Write(lun_diag,"(a5,4es17.9)") nname(i),rhs(i),ydot(i),yt(i),y(i)
      EndDo
    EndIf
!   write(*,*) "before jacobian_solve"
! Solve the jacobian and calculate the changes in abundances, dy        
    Call qjacobian_solve(kstep,rhs,dy)
!   Call qjacobian_decomp(kstep)
!   Call qjacobian_bksub(rhs,dy)
!   write(*,*) "after jacobian_solve"
! Evolve the abundances and calculate convergence tests
!    yt=yt+dy
!    Where(yt<ymin) 
!      yt=0.0
!      reldy=0.0
!    ElseWhere
!      reldy=abs(dy/yt)
!    EndWhere
!    If(idiag>=3) Then
!      irdymx=maxloc(reldy)
!      idymx=maxloc(dy)
            testc=0.0
            yqt=yqt+dy
            call update2(yqt,yt)
            Where(yqt<ymin)
              yqt=0.0
              reldy=0.0
            ElseWhere
              reldy=abs(dy/yqt)
            EndWhere
            If(idiag>=3) Then
              irdymx=maxloc(reldy)
              idymx=maxloc(dy)
      Write(lun_diag,"(a2,i5,i3,2(a5,2es12.4))") &
&       'dY',kstep,kit,nname(idymx(1)),dy(idymx(1)),y(idymx(1)), &
&       nname(irdymx(1)),reldy(irdymx(1)),y(irdymx(1))
      If(idiag>=4) Write(lun_diag,"(a5,5es12.4)") &
&       (nname(k),yt(k),dy(k),reldy(k),(aa(k)*dy(k)),(aa(k)*yt(k)),k=1,ny)
    EndIf
  
!-----------------------------------------------------------------------------  
! There are 3 included convergence tests: testc, which measures relative
! changes, testc2 which measures total abundance changes, and testm
! which tests mass conservation.  
!-----------------------------------------------------------------------------  
    testc=sum(reldy)
    testc2=sum(aaq*dy)
    xtot=sum(aaq*yqt)-1.0
    testm=xtot-xtot_init
    If(idiag>=2) Write(lun_diag,"(a,2i5,3es14.6)") 'NR',kstep,kit,testm,testc,testc2
  
!-----------------------------------------------------------------------------  
! testc is the most stringent test, and requires the most iterations.  
! testm is the most lax, and therefore the fastest, often requiring only one 
! iteration.  Considering the uncertainties of the reaction rates, it is
! doubtful that the increased precision of testc is truly increased
! accuracy. 
!-----------------------------------------------------------------------------  
! Ordinarily, test for true convergence
    If (iconvc/=0) Then
      testn=testc
      toln=tolc
    
! Otherwise, Use mass conservation for convergence condition 
    Else
      testn=testm
      toln=tolm           
    EndIf
  
! If converged, exit NR loop
    If(abs(testn)<=toln) Exit
        
  EndDo
  
! If NR loop exits early, step is successful.
  If(kit<=kitmx) Then
    If(idiag>=1) Write(lun_diag,"(a,2i5,3es12.4)") 'Conv',kstep,kit,xtot,testn,toln
    inr=0

! Otherwise, convergence failed
  Else    
    If(idiag>=1) Write(lun_diag,"(a,2i5,2es10.2)") 'BE Failure',kit,kitmx,xtot,testn
    inr=1
  EndIf
  
End Subroutine step_qse
