!*******************************************************************************
! Backward Euler solver, part of XNet 7 5/25/10 
!
! These routines in this file perform the Backward Euler time integration for the 
! thermonuclear reaction network.
!*******************************************************************************
  
Subroutine solve_be(kstep,its)
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
  Integer, Intent(in) :: kstep  
  Integer, Intent(out) :: its 
  Integer, Parameter :: ktsmx=10
  Integer kts,inr,kit,j
!-----------------------------------------------------------------------------  
! For each trial timestep, tdel, the Newton-Raphson iteration is attempted.
! The possible results for a timestep are a converged set of yt or a 
! failure to converge.  If convergence fails, iteration is retried with 
! the trial timestep reduced by tdel_maxmult, up to ktsmx times.
!-----------------------------------------------------------------------------  
  Do kts=1,ktsmx
! Attempt Backward Euler integration over desired timestep
    Call step_be(kstep,inr,kit)
    
! If integration fails, reset abundances, reduce timestep and retry.  
    If(inr/=0) Then
      tdel=tdel/tdel_maxmult
      tt=t+tdel
      If(idiag>=1) Write(lun_diag,"(a,i3,3es12.4)") 'BE TS Reduce',kts,tt,tdel,tdel_maxmult
  
! Set iweak to the original value before possibly changing it again in t9rhofind
      iweak=iweak0
      Call t9rhofind(kstep,tt,j)
  
! Reset abundances
      yt=y
    Else
      Exit
    EndIf
  EndDo
! If TS loop exits early, step is successful.  Prepare for next step.
  If(kts<=ktsmx) Then
    its=0
    If(idiag>=1) Write(lun_diag,"(a,i5,2i3)") 'TS Success',kstep,kts,kit
    to=t
    t=tt
    tdel_next=tdel*tdel_maxmult
!   Call norm(y)
    yo=y
    y=yt
!    call group_number
!    call group_ab_test(kstep) 
  If (iadpt==1)then
!  write(*,*) "si28=",y(77)
  call switcher(kstep) 
  EndIF
  Else
    Write(6,"(i5,4es12.4,2i3)") kstep,t,tdel,t9t,rhot,kit,kts
    Write(6,*) 'Timestep retrys fail after ',kts,' attempts'
    its=1
  EndIf
  kmon(1)=kts ; kmon(2) = kit
End Subroutine solve_be
  
Subroutine step_be(kstep,inr,kit)
!===============================================================================
! This routine attempts to integrate a single Backward Euler step for the 
! timestep tdel.  If successful, inr=0
!===============================================================================
  Use controls
  Use nuclear_data
  Use conditions
  Use abundances
  Use reac_rate_data
  Integer, Intent(in) :: kstep
  Integer, Intent(out) :: inr, kit
  Integer irdymx(1),idymx(1)
  Integer :: i,j
  Real(8) :: toln    ! Network convergence condition
  Real(8) :: testc,testc2,testm,testn  ! Convergence tests
  Real(8) :: xtot,xtot_init,rdt
  Integer :: k,idiag0
  Real(8), Dimension(ny) :: rhs,dy,reldy
  
! Calculate the thermodynamic factors necessary for reaction rates, 
! including screening, and the reaction rates.
  Call cross_sect 
  
  xtot_init=sum(aa*y)-1.0
  
! The Newton-Raphson iteration occurs for at most kitmx iterations.  
  Do kit=1,kitmx
               
! Calculate the reaction rates and abundance time derivatives
    Call yderiv
  
! Build the Jacobian
    rdt=1.0/tdel
    Call jacobian_build(rdt,-1.0)
  
! Calculate equation to zero
    rhs=(y-yt)*rdt+ydot
    If(idiag>=1) Then
      Write(lun_diag,"(a3,i5,es14.7)") 'RHS',kstep,rdt
      Do i=1,ny
        Write(lun_diag,"(a5,4es17.9)") nname(i),rhs(i),ydot(i),yt(i),y(i)
      EndDo
    EndIf
   
! Solve the jacobian and calculate the changes in abundances, dy        
    Call jacobian_solve(kstep,rhs,dy)
!   Call jacobian_decomp(kstep)
!   Call jacobian_bksub(rhs,dy)
  
! Evolve the abundances and calculate convergence tests
    yt=yt+dy
    Where(yt<ymin) 
      yt=0.0
      reldy=0.0
    ElseWhere
      reldy=abs(dy/yt)
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
    testc2=sum(aa*dy)
    xtot=sum(aa*yt)-1.0
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
  
End Subroutine step_be
