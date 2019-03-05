!*******************************************************************************
!  Full_Net 4.10.1 11/24/09 
!  The subroutines in this file perform nucleosynthesis for a single Lagrangian 
!  mass zone.  The handling of multiple zones, including multi-processing, is 
!  carried out externally.  
!
!  The data.f file contains the nuclear and reaction data structures and the 
!  subroutines which read in the data and allocate the arrays.  
!*******************************************************************************

      subroutine full_net(izone)
!===============================================================================
!  The abundances are evolved using a Newton-Raphson iteration scheme to solve 
!  the equation yt(i)-y(i)/tdel = ydot(i), where y(i) is the abundance at the 
!  beginning of the iteration, yt(i) is the trial abundance at the end of the 
!  timestep, ydot(i) is the time derivative of the trial abundance  calculated 
!  from reaction rates and tdel is the timestep.  
!===============================================================================
      use controls
      use nuclear_data
      use conditions
      use abundances
      use thermo_data
      real(8) :: reldy(ny)
      real(8) :: toln      !  Network convergence condition
      real(8) :: testc,testc2,testm,testn    ! Convergence tests
      real(8) :: ta,tdela,xtot,xtoto,enm,enb,enold,en0,edot,ye
      real(8) :: ytot,ztot,atot
      integer irdymx(1),idymx(1)
      integer :: izone,k,kstep,ktsmx,kout,idiag0

!  Set reaction controls not read in from control
      iweak0=iweak
      idiag0=idiag
      ktsmx=10
      kstep=0
      kout=0

!  Normalize initial abundances and change units if necessary
      yt=y
!     call norm(yt) 

!  Calculate the total energy of the nuclei
      call benuc(yt,enb,enm,ytot,ztot,atot)
      xtot=atot-1.0
      en0=enm
      edot=0.0

!  Start evolution 
      Write(6,"(a,i6,a,i2,2(a,es10.3))") 'Max Step',kstmx,' 
     &  IDiag=',idiag,' Start Time',tstart,' Stop Time',tstop
      t=tstart
      tt=tstart
      call t9rhofind(kstep,tt,j)
      
      If(itso>0) call ts_output(kstep,(enm-en0),edot)
!-----------------------------------------------------------------------------  
!  For each step, an initial guess for the timestep is calculated, based on 
!  the abundance variations of the previous step if available.  
!-----------------------------------------------------------------------------  
      Step: Do kstep=1,kstmx
        call timestep(kstep)
        If(idiag>=1) Write(lun_diag,*) 'TDel',tt,tdel
        xtoto=xtot

!  Determine if this is an output step
        idiag=idiag0
!       If(mod(kstep,10).eq.0) idiag=2
!       If(kstep==518.or.kstep==718.or.kstep==803) idiag=5
!-----------------------------------------------------------------------------  
!  For each trial timestep, tdel, the Newton-Raphson iteration is attempted.
!  The possible results for a timestep are a converged set of yt or a 
!  failure to converge.  If convergence fails, iteration is retried with 
!  the trial timestep reduced by tdelmm, up to ktsmx times.
!-----------------------------------------------------------------------------  
        TS: Do kts=1,ktsmx

!  Calculate the thermodynamic factors necessary for reaction rates, 
!  including screening, and the reaction rates.
          call cross_sect 

!  The Newton-Raphson iteration occurs for at most knrmx iterations.  
          NR: Do knr=1,knrmx
                   
! Calculate the reaction rates and abundance time derivatives
            call yderiv

!  Build the Jacobian and calculate the changes in abundances, dy
            call netmatr(kstep)

!  Evolve the abundances and calculate convergence tests
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
              Write(lun_diag,"(a2,i5,2i3,2(a5,2es12.4))") 
     &          'dY',kstep,kts,knr,nname(idymx(1)),
     &          dy(idymx(1)),y(idymx(1)),nname(irdymx(1)),
     &          reldy(irdymx(1)),y(irdymx(1))
              If(idiag>=4) Write(lun_diag,"(a5,5es12.4)") 
     &          (nname(k),yt(k),dy(k),reldy(k),(aa(k)*dy(k)),
     &          (aa(k)*yt(k)),k=1,ny)
            Endif

!-----------------------------------------------------------------------------  
!  There are 3 included convergence tests: testc, which measures relative
!  changes, testc2 which measures total abundance changes, and testm
!  which tests mass conservation.  
!-----------------------------------------------------------------------------  
            testc=sum(reldy)
            testc2=sum(aa*dy)
            xtot=sum(aa*yt)-1.0
            testm=xtot-xtoto
            If(idiag>=2) Write(lun_diag,"(a3,i5,i3,3es14.6)") 
     &        'KNR',kstep,knr,testm,testc,testc2

!-----------------------------------------------------------------------------  
!  testc is the most stringent test, and requires the most iterations.  
!  testm is the most lax, and therefore the fastest, often requiring only one 
!  iteration.  Considering the uncertainties of the reaction rates, it is
!  doubtful that the increased precision of testc is truly increased
!  accuracy. 
!-----------------------------------------------------------------------------  
!  Ordinarily, test for true convergence
            If (iconvc/=0.or.tt>=tstop) Then
              testn=testc
              toln=tolc
  
!  Otherwise, use mass conservation for convergence condition 
            Else
              testn=testm
              toln=tolm           
            Endif
            If(abs(testn)<=toln) Exit TS
          Enddo NR

!  If convergence is not achieved in knrmx iterations, reset abundances 
!  and try again with the timestep reduced.  
          If(idiag>=1) Write(lun_diag,*) 'TS Failure',knr,kts,xtot,testn
          tdel=tdel/tdelmm
          tt=t+tdel

! Set iweak to the original value before possibly changing it again in t9rhofind
          iweak=iweak0
          call t9rhofind(kstep,tt,j)
          yt=y
        Enddo TS

!  If convergence is successful, update time and abundances 
        If(kts<ktsmx) Then
          If(idiag>=1) Write(lun_diag,"(a4,i5,i3,3es12.4)") 
     &      'Conv',kstep,knr,xtot,testn,toln
          ta=t
          tdela=tdel
          t=tt
!         call norm(y)
          yo=y
          y=yt
          ye=sum(zz*y)
          If(idiag>=2) Then
            Write(lun_diag,"(a)") 'delta Y'
            Write(lun_diag,"(a5,4es12.4)") (nname(k),y(k),yo(k),
     &        (yt(k)-y(k)),(tdel*ydot(k)),k=1,ny)
          Endif
          enold=enm
          call benuc(yt,enb,enm,ytot,ztot,atot)
          edot=-(enm-enold)/tdel
          If(idiag>=1) Then 
            Write(lun_diag,"(i5,5es14.7)") kstep,t,tdel,t9t,rhot,ye
!           Write(lun_diag,"(5(a5,es11.4))") (nname(k),y(k),k=1,ny)
          Endif
          If(itso>0) call ts_output(kstep,(enm-en0),edot)
          If(t>=tstop) Then
            Exit STEP
          Endif

!  If reduced timesteps fail to yield convergence, warn and exit
        Else
          Write(6,"(i5,4es12.4,2i3)") kstep,t,tdel,t9t,rhot,knr,kts
          Write(6,*) 'Timestep retrys fail after ',kts,' attempts'
          Exit STEP
        Endif
      Enddo STEP

!  Test that the stop time is reached
      If(t<tstop) Then
        Write(lun_diag,"(a,es12.4,a,es12.4,a)") 
     &       'Evolution incomplete!!!' 
        Write(6,"(a,es12.4,a,es12.4,a)") 'Evolution stopped at time=',t,
     &  '.  Stop time (',tstop,') not reached!' 
        Write(6,"(a,i6)") 
     &    'Approximately',int((tstop-t)/tdel),'more steps needed' 
      Endif

!  End Post Processing cycle
      iweak=iweak0
      call final_output(kstep)
      Return
      End subroutine full_net
 

      subroutine yderiv
!-----------------------------------------------------------------------------  
!  This routine calculates time derivatives for each nuclear species.
!  This calculation is performed by looping over nuclei, and summing the 
!  reaction rates for each reaction which effects that nucleus.
!-----------------------------------------------------------------------------  
      use controls
      use nuclear_data
      use abundances
      use cross_sect_data
      use reac_rate_data
      integer :: i,j,i0,i1,la1,le1,la2,le2,la3,le3
      real(8) :: s1,s11,s2,s22,s3,s33

!  From the cross sections and the counting array, calculate the reaction 
!  rates
      b1=a1*csect1(mu1)
      b2=a2*csect2(mu2)
      b3=a3*csect3(mu3)

!  Calculate Ydot for each nucleus, summing over the reactions which affect it.
      Do i0=1,ny                 
        If(idiag>=5) Write(lun_diag,"(a3,a6,i4)") 'NUC',nname(i0),i0

!  Sum over the reactions with 1 reactant
        la1=la(1,i0)                          
        le1=le(1,i0)                         
        s1=0.0                              
        Do i1=la1,le1                    
          s11=b1(i1)*yt(n11(i1))
          s1=s1+s11
          If(idiag>=5) Write(lun_diag,"(3a5,'  1  ',4es15.7)") 
     &      (nname(n1i(j,mu1(i1))),j=1,3),b1(i1),yt(n11(i1)),
     &      s11,a1(i1)
        Enddo       
        If(idiag>=5) Write(lun_diag,*) '1->',nname(i0),la1,le1,s1

!  Sum over the reactions with 2 reactants
        la2=la(2,i0)  
        le2=le(2,i0) 
        s2=0.0      
        Do i1=la2,le2           
          s22=b2(i1)*yt(n21(i1))*yt(n22(i1))
          s2=s2+s22                                       
          If(idiag>=5) Write(lun_diag,"(4a5,4es15.7)")
     &      (nname(n2i(i,mu2(i1))),i=1,4),b2(i1),yt(n21(i1)),
     &      yt(n22(i1)),s22
        Enddo                  
        If(idiag>=5) Write(lun_diag,*) '2->',nname(i0),la2,le2,s2

!  Sum over the reactions with 3 reactants
        la3=la(3,i0)                          
        le3=le(3,i0)                         
        s3=0.0                              
        Do i1=la3,le3                  
          s33=b3(i1)*yt(n31(i1))*yt(n32(i1))*yt(n33(i1))
          s3=s3+s33
          If(idiag>=5) Write(lun_diag,"(3a5,'  3  ',5es12.4)") 
     &      (nname(n3i(i,mu3(i1))),i=1,3),b3(i1),yt(n31(i1)),
     &      yt(n32(i1)),yt(n33(i1)),s33
        Enddo                                       
        If(idiag>=5) Write(lun_diag,*) '3->',nname(i0),la3,le3,s3

!  Sum the 3 components of Ydot
        ydot(i0)=s1+s2+s3                         
        If(idiag>=5) Write(lun_diag,"(a4,a5,2es24.16)") 
     &    'YDOT',nname(i0),yt(i0),ydot(i0)
      Enddo                                 
      Return                                  
      End subroutine yderiv                                   

