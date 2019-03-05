!*******************************************************************************
! This is an example driver for running FullNet on serial machines, 
! including the potential Use of OpenMP threads. 
!
! This driver reads in the controlling flags, the nuclear and 
! reaction data. It also loops through reading the initial abundances 
! and the thermodynamic trajectory before Calling full_net for each zone.  
!
! This file also contains the ts_output routine, which controls the 
! information output at the End of each timestep, and the final_output 
! routine, which controls what is output an the End of the run. 
!
! If you wish to Use full_net in concert with a hydrodynmics code, 
! you will need to supply these services from within the hydro. 
!*******************************************************************************
   
  Program net
!===============================================================================
!===============================================================================
  Use controls
  Use conditions
  Use abundances
  Use thermo_data
  Use screen_qse
  Use nuclear_data
!$ Use omp_lib
  Integer :: i,j,k,n,nz,izone! Loop indices
  Integer :: szone,kstep,kstart,ii(14),ierr,index
  Integer, Parameter   :: nzmx= 100 ! Max number of zones
  Integer :: nhL(nzmx)
  Real(8), Dimension(:,:), Allocatable :: yinL
  Real(8) :: tdelstart(nzmx)
  Real(8) :: tbegin,tbeginL(nzmx)
  Real(8) :: tstartL(nzmx)
  Real(8) :: tEnd,tEndL(nzmx)
  Real(8) :: tcut, tcutL(nzmx)
  Real(8) :: thL(nhmx,nzmx),t9hL(nhmx,nzmx)
  Real(8) :: rhohL(nhmx,nzmx),t9oldL(nzmx)
  Real(8) :: t9cutL(nzmx)
  Character (LEN=5)  :: output_nuc(14)
  Character (LEN=80) :: descript(3),data_descL(nzmx),data_desc
  Character (LEN=80) :: abund_desc(nzmx),thermo_desc(nzmx)
  Character (LEN=80) :: data_dir,thermo_file(nzmx),inab_file(nzmx)
  Character (LEN=80) :: ev_file_base,ev_file,bin_file_base,bin_file
  Character (LEN=80) :: diag_file,diag_file_base='net_diag'
!$OMP THREADPRIVATE(ev_file,bin_file,diag_file,tEnd,tcut)
  
!-----------------------------------------------------------------------
! Read control file, which contains the Parameters which control the 
! action of the network.  
!-----------------------------------------------------------------------
  Open(5,FILE='control')       
  write(*,*) "1"                                  
  Read(5,"(72x)")
  Read(5,"(a80)") (descript(i),i=1,3) ! text description of the problem.
  Read(5,"(72x)")
  write(*,*) "1"                                  
  Read(5,*) szone! number of the zone with which to begin
  Read(5,*) nzone! total # of zones
  Read(5,*) isolv ! Choice of integrations scheme
  Read(5,*) iadpt !  Adaptive or not! 
  Read(5,*) kstmx! max # of timesteps for each zone
  Read(5,*) kitmx! max # of iterations before retry
  Read(5,*) idiag! sets diagnostic output level
  Read(5,*) itsout ! sets per timestep output level
  Read(5,*) iweak! controls the treatment of weak reactions 
  Read(5,*) iscrn! controls the treatment of nuclear screening
  Read(5,*) iconvc! determines which convergence condition is Used
  Read(5,*) changemx! allowed abundance change Used to set the timestep.
  Read(5,*) yacc! abundances > yacc Used for timestep calculation
  Read(5,*) tolm! mass conservation convergence criterion
  Read(5,*) tolc ! convergence limit on the iterative abundance change
  Read(5,*) ymin ! abundance < ymin is set to 0.0
  Read(5,*) tdel_maxmult ! max factor by which the timestep is changed
!-----------------------------------------------------------------------
! control also indicates the base of the filenames to which ASCII and 
! binary output are written at the End of each timestep, the relative 
! directory from which the nuclear data should be loaded, as well as 
! the names of the files containing the initial abundances and
! thermodynamic trajectories.
!-----------------------------------------------------------------------
  Read(5,"(72x)")
  Read(5,"(14a5)") output_nuc            
  Read(5,"(72x)")
  Read(5,"(a80)") ev_file_base            
  Read(5,"(72x)")
  Read(5,"(a80)") bin_file_base            
  Read(5,"(72x)")
  Read(5,"(a80)",End=10) data_dir
  Read(5,"(72x)")
  Do nz=1,nzone
    Read(5,"(a80)",End=10) inab_file(nz)
    Read(5,"(a80)",End=10) thermo_file(nz)
  EndDo
     10 nz=nz-1
  If((nz)/=nzone) Write(6,*) nz,' datafiles for ',nzone,' zones!'
! Write(6,*) inab_file(nz),thermo_file(nz)
  
! IdentIfy threads
!$OMP PARALLEL DEFAULT(SHARED)
  mythread=1
!$ mythread=OMP_get_thread_num()      
  nthread=1
!$ nthread=OMP_get_num_threads()    
  
!-----------------------------------------------------------------------
! Output files:
! The diagnostic output is sorted by thread, lun_diag = 50+mythread
!
! General output is per zone, however logical units are spaced per thread
! lun_ev = 50+nthread+myid ; lun_ts ; unit 50+2*nthread+myid     
!-----------------------------------------------------------------------
! Open diagnositic output file, per thread If OMP
  If(idiag>=0) Then
    diag_file=trim(diag_file_base)
    Call name_ordered(diag_file,mythread,nthread)
    Open(lun_diag,file=diag_file)
    Write(lun_diag,"(a,i4,a,i4)") 'Thread ',mythread,' of ',nthread
  EndIf
!$OMP End PARALLEL
  
! Retain the value of iweak
  iweak0=iweak
  
! Read nuclear data (proton and atomic number, mass, partition functions)
  Call read_nuclear_data(data_dir,data_desc)
  
! Read reaction rate data
  Call read_reaction_data(data_dir)
  If(idiag>=0) Write(lun_diag,"(a)") (descript(i),i=1,3),data_desc
  
! Read jacobian matrix data 
  Call read_jacobian_data(data_dir)
   
! Read data on matching forward and reverse reactions 
  Call read_match_data(data_dir)
  
! Initialize flux tracking       
  Call flux_init
  
! Convert output_nuc names into indices
  Do i=1,14
    Call index_from_name(output_nuc(i),index)
    If(index<1.or.index>ny) Then
      Write(6,*) 'Output Nuc:',i,output_nuc(i),' not found'
      inout(i)=ny
    Else
      inout(i)=index
    EndIf
  EndDo
  
! Set sizes of abundance arrays
!$OMP PARALLEL DEFAULT(SHARED)
  Allocate (y(ny),yo(ny),yt(ny),ydot(ny))
       zmax=maxval(zz)
      allocate(he(zmax))

!$OMP End PARALLEL
  
!-----------------------------------------------------------------------
! The initial abundances and thermo-data are read into a series of local
! arrays that are functions of the zone. All of the variables associated 
! with the local arrays End in "L". The local arrays are then loaded by 
! a loop over zones into global varaibles and dIfferent threads according 
! to zone in the OMP parallel region.  The local arrays are Allocated to 
! have the zone size set to the maximum number of zones, nzmx. 
!-----------------------------------------------------------------------
  Allocate (yinL(ny,nzone))
  
! For each Zone
  Do izone=szone,nzone
  
! Read Initial Abundances
    Open(8,file=inab_file(izone))
    kstart=1
    Read(8,"(a)") abund_desc(izone)
    Read(8,"(4(5x,es14.7,1x))") (yinL(i,izone),i=1,ny)
    Close(8)
!   Write(lun_ts) inab_file(izone),abund_desc
!   Write(lun_diag,"(a)") abund_desc
!   Call sum_test(yinL)
! Read the thermdynamic trajectory 
    tstartL=0.0
    tcutL=0.0
    t9oldL=0.0
    Open(8,file=thermo_file(izone))
    Read(8,"(a)") thermo_desc(izone)
    Read(8,*) tbeginL(izone)
    Read(8,*) tEndL(izone)
    Read(8,*) tdelstart(izone)
    
    Do n=1,nhmx
      Read(8,*,IOSTAT=ierr) thL(n,izone),t9hL(n,izone),rhohL(n,izone)
      If(ierr==-1) Then
!       Write(6,*) 'End of File Reached after',i,' records'
        Exit
      EndIf
    EndDo
    nh=n-1
    nhL(izone)=nh
    Close(8)
  EndDo
  
! For each zone
!$OMP PARALLEL DEFAULT(SHARED) !COPYIN(nh,tstart,tstop,th,t9h,rhoh)
!$OMP DO
  Do izone=szone,nzone

! Load in zone trajectory's timestep, temperature and density
    nh=nhL(izone)
    Do i=1,nh
      th(i)=thL(i,izone)
      t9h(i)=T9hL(i,izone)
      rhoh(i)=rhohL(i,izone)
    EndDo
    tstart=max(tbeginL(izone),tstartL(izone))
    tcut=tcutL(izone)
    tbegin=tbeginL(izone)
    tEnd=tEndL(izone)
    data_desc=data_descL(izone)
  
! Open the evolution file
    lun_ev=50+nthread+mythread
    If(itsout>=2) Then
      ev_file=ev_file_base
      Call name_ordered(ev_file,izone,nzone)
      If(idiag>=0) Write(lun_diag,"(a,i5,7es10.3)") trim(ev_file),&
&       nh,th(nh),t9h(nh),rhoh(nh),tstart,tstop,tbegin,tEnd
      Open(lun_ev,file=ev_file)
  
! Log abundance description
      If(idiag>=0) Write(lun_diag,"(a)") abund_desc(izone)
! Log thermo description
      If(idiag>=0) Write(lun_diag,"(a)") thermo_desc(izone)
  
      ii=(/5,21,29,36,43,53,63,73,89,102,115,117,129,144/)
      Write(lun_ev,"(a4,a15,4a10,15a9,a4)") &
&       'k ',' Time ',' T(GK) ',' Density ',' dE/dt ',' Timestep ',nname(inout), ' It '
    EndIf
  
! Open the time series file
    lun_ts=50+2*nthread+mythread 
    If(itsout>=1) Then
      bin_file=bin_file_base
      Call name_ordered(bin_file,izone,nzone)
      Open(lun_ts,file=bin_file,form='unformatted')
! Write Control Parameters to ts file
      Write(lun_ts) (descript(i),i=1,3),data_desc
      Write(lun_ts) kstmx,kitmx,iweak,iscrn,iconvc,changemx,tolm,tolc,yacc,ymin,tdel_maxmult
! Write abundance description to ts file
      Write(lun_ts) inab_file(izone),abund_desc(izone)
! Write thermo description to ts file
      Write(lun_ts) thermo_file(izone),thermo_desc(izone)
    EndIf
  
! Load initial abundances, time and timestep
    tdel=0.0
    Do i=1,ny
      y(i)=yinL(i,izone)
    EndDo
    If(tcut>0.0) Then
      tstop=min(tEnd,tcut)
    Else
      tstop=tEnd
    EndIf
       
! Evolve abundance from tstart to tstop, using the original iweak
    iweak=iweak0
    Call full_net(izone)
  
! For evolution at temperatures below T9=t9min, we Use only weak reactions
    If(tstop<tEnd) Then
      iweak=-1   
      tstart=t   
      tstop=tEnd 
      Call full_net(izone)
    EndIf
!  yo=yinL
    Close(lun_ev)
    Close(lun_ts)
  EndDo
  
!$OMP End DO
!$OMP End PARALLEL
  
  Stop
  End
  
Subroutine ts_output(kstep,enuc,edot) 
!===============================================================================
! The per timestep output routine.  If the flag itsout is > 0, 
! full_net Calls this routine to handle stepwise output.  
!===============================================================================
  Use nuclear_data
  Use controls
  Use conditions
  Use abundances
  Use match_data
  Use flux_data
    use qse_data
    use screen_qse
    use part_funct_data
    real, parameter :: pi=3.14159,hbr=6.58217D-22,amu=1.036435E-18
    real, parameter :: bok=8.6174D-02,avn=6.02205D+23

     double precision xCqse(ny), xdlCqse(ny), tmpa,Yqse(ny),Yqse2(ny)
     double precision Fesum,t9oldy,Ynse(ny)
     real(8) :: xg(7),tmp,bkt,dummy,ratiosi,ratiofe, sumsi,sumfe,rsumsi,rsumfe
  
  Real(8) :: enuc,edot,yout(0:ny)
  Integer :: i,j,k,l,kstep,kout,aa1(ny), dummycount
!$OMP THREADPRIVATE(yout)
  
!$OMP PARALLEL DEFAULT(SHARED)
  yout(0)=0.0
  yout(1:ny)= y
!$OMP End PARALLEL
  
! Calculate reaction fluxes
  If(kstep>0) Then
     Call flux
  Else
     flx_int=0.0
  EndIf
  
! An abundance snapshot is written to the binary file on unit 24.
  Write(lun_ts) kstep,t,t9t,rhot,tdel,y,flx
!  Write(208, *) kstep,t,t9t,rhot,tdel
  
! For itsout>=2, output important mass fractions to the ASCII file on unit 22
  If(itsout>=2) Then
    xg=0.0
!  Do i=1,ny
!    If(zz(i)<=1) Then
!      ig=1
!     ElseIf(zz(i)<=2) Then
!       ig=2
!     ElseIf(zz(i)<=8) Then
!       ig=3
!     ElseIf(zz(i)<=10) Then
!       ig=4
!     ElseIf(zz(i)<=12) Then
!       ig=5
!     ElseIf(zz(i)<=14) Then
!       ig=6
!     Else
!       ig=7
!     EndIf
!     xg(ig)=xg(ig)+aa(i)*y(i)
!   EndDo
!   Write(lun_ev,"(i4,1es15.8,2es10.3,2es10.2,8es9.2,i4)") &
!&     kstep,t,t9t,rhot,edot,enuc,tdel,(xg(i),i=1,7),kout
!   ii=(/6,21,33,45,61,77,95,115,139,162,184,206,231,257/)
    Write(lun_ev,"(i4,1es15.8,2es10.3,2es10.2,14es9.2,2i2)") &
&     kstep,t,t9t,rhot,edot,tdel,(aa(inout)*yout(inout)),kmon
!  Write(lun_ev,"(i4,1es15.8,2es10.3,2es10.2,14es9.2,i4)") kstep,t,t9t,rhot,edot,tdel,(aa*y),kout
    Write(901,"(1es15.8,14es9.2)") &
&     t,(aa(inout)*yout(inout))
  
! For itsout>=3, output time and thermo evolution to the screen
    If(itsout>=3) Write(6,"(i5,4es12.4,2i3)") kstep,t,tdel,t9t,rhot,kmon
  
! For itsout>=4, output abundances for each timestep to the diagnostic file
    If(itsout>=4) Write(lun_diag,"(4(a5,es14.7,1x))") (nname(i),aa(i)*y(i),i=1,ny)
  
! For itsout>=5, output fluxes for each timestep to the diagnostic file
    If(itsout>=5) Write(lun_diag,'(i5,8a5,i5,es11.3)') &
&     (k,nname(nflx(1:3,k)),' <-> ',nname(nflx(4:7,k)),iwflx(k),flx(k),k=1,mflx)
  EndIf
!--------------------------------------------------------------------------
! For itsout>=3 print out diaignotic information for QSE
        If(itsout>=3)then !!!!!!!!!!!QSE itsout if

! temp brackets
      do i=1,ny
       aa1(i)=int(aa(i))
       enddo


       OPEN(406,file='T6.dat')
       OPEN(405,file='T5.dat')
       OPEN(404,file='T4.dat')
       OPEN(435,file='T35.dat')
       OPEN(403,file='T3.dat')
       OPEN(407,file='T45.dat')
       OPEN(408,file='25.dat')
       OPEN(409,file='37.dat')
       OPEN(410,file='T10.dat')
       If(t9t.ge.5.95.and.t9t.le.6.05)then
         Write(406,'(i5,1x,1es13.6,1x,2es13.6)') kstep,t,t9t,rhot
         Write(406,"((a5,i4,es14.7,1x))") (nname(i),aa1(i),y(i),i=1,ny)
       Endif
       If(t9t.ge.9.95.and.t9t.le.10.05)then
         Write(410,'(i5,1x,1es13.6,1x,2es13.6)') kstep,t,t9t,rhot
         Write(410,"((a5,i4,es14.7,1x))") (nname(i),aa1(i),y(i),i=1,ny)
       Endif

       If(t9t.ge.4.95.and.t9t.le.5.05)then
         Write(405,'(i5,1x,1es13.6,1x,2es13.6)') kstep,t,t9t,rhot
         Write(405,"((a5,i4,es14.7,1x))") (nname(i),aa1(i),y(i),i=1,ny)
       Endif

       If(t9t.ge.3.95.and.t9t.le.4.05)then
         Write(404,'(i5,1x,1es13.6,1x,2es13.6)') kstep,t,t9t,rhot
         Write(404,"((a5,i4,es14.7,1x))") (nname(i),aa1(i),y(i),i=1,ny)
       Endif
       If(t9t.ge.3.47.and.t9t.le.3.525)then
         Write(435,'(i5,1x,1es13.6,1x,2es13.6)') kstep,t,t9t,rhot
         Write(435,"((a5,i4,es14.7,1x))") &
     & (nname(i),aa1(i),y(i),i=1,ny)
       Endif
       If(t9t.ge.2.99.and.t9t.le.3.015)then
         Write(403,'(i5,1x,1es13.6,1x,2es13.6)') kstep,t,t9t,rhot
         Write(403,"((a5,i4,es14.7,1x))") (nname(i),aa1(i),y(i),i=1,ny)
       Endif

       If(t9t.ge.4.39.and.t9t.le.4.515)then
         Write(407,'(i5,1x,1es13.6,1x,2es13.6)') kstep,t,t9t,rhot
         Write(407,"((a5,i4,es14.7,1x))") (nname(i),aa1(i),y(i),i=1,ny)
       Endif
       If(t9t.ge.2.20.and.t9t.le.2.55)then
          Write(408,'(i5,1x,1es13.6,1x,2es13.6)') kstep,t,t9t,rhot
          Write(408,"((a5,i4,es14.7,1x))") (nname(i),aa1(i),y(i),i=1,ny)
       Endif
       If(t9t.ge.3.65.and.t9t.le.3.75)then
         Write(409,'(i5,1x,1es13.6,1x,2es13.6)') kstep,t,t9t,rhot
         Write(409,"((a5,i4,es14.7,1x))") (nname(i),aa1(i),y(i),i=1,ny)
       Endif
!Makes inital abundance profs.
!       if(kstep>100)then
!        Write(417,'(i5,1x,1es13.6,1x,2es13.6)') kstep,t,t9t,rhot
!        Write(417,"(4(a5,es14.7,1x))") (nname(i),y(i),i=1,ny)
!       endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!-----------------------------------------------------------
!! The following calculates QSE cofficents around si28.
! It is used to determine which kstep to yank an initial
! abundance profile for qse.exe from.
!-----------------------------------------------------------

         bkt=bok*t9t
         tmp=dlog(((avn*rhot)**(2./3.)*2.*pi*hbr**2)/(bkt*amu))
         xCqse(1)=1.
         xdlCqse(1)=0.
         xCqse(2)=1.
         xdlCqse(2)=0.
         DO l=3,ny
            tmpa=dlog(angm(l)*gg(l)*(aa(l)**1.5)*(0.5**aa(l)))
            xdlCqse(l)=tmpa+(1.5*(aa(l)-1))*tmp+(be(l)/bkt)+he(int(zz(l)))
            xCqse(l)=dexp(xdlCqse(l))

         Enddo
!220    CONTINUE

!---------------------------------------------------------------------
! Now that the QSE coefficent has been calcuated, the Qse abundances can
! be calucated near Si28.Si29 is nucli # 77 so Cqse(77) is the qse
! cofficinet for si28.
!----------------------------------------------------------------------
       write(929,*) t, kstep,rhot,t9t
!       write(937,*) t, kstep
       write(923,*) "step",t, kstep,t9t
       DO i=1,ny
        Yqse(i)=dexp(xdlCqse(i)-xdlCqse(43))*y(43) &
     &  *y(2)**(zz(i)-14)*y(1)**(nn(i)-14)
        Yqse2(i)=dexp(xdlCqse(i)-xdlCqse(104))*y(104) &
     &  *y(2)**(zz(i)-24)*y(1)**(nn(i)-26)
        Ynse(i)=dexp(xdlCqse(i)+dlog(y(1))*nn(i)+dlog(y(2))*zz(i))
       Enddo
!----------------------------------------------------------------------
! this routine tries to establish the groups by comparing YQSE/T.
! Where YQSE/y is within 10% of 1 we consider group membership.
! For the 299 network si spans from 60-130 fe from 147-299.
!----------------------------------------------------------------------
       write(987,*) kstep
       write(988,*) kstep
       sumsi=0
       sumfe=0
       dummycount=0
       Do i=1,ny
         ratiosi=Yqse(i)/y(i)
         ratiofe=Yqse2(i)/y(i)
      
         !If(ratiosi>=.80.and.ratiosi<=1.20)then
         ! write(987,*) i, nname(i),Yqse(i),y(i),ratiosi
         !Endif

         !If(ratiofe>=.80.and.ratiofe<=1.20)then
         write(987,*) "fe",i,nname(i),y(i),Yqse2(i),ratiofe
      !   Endif
       Enddo
       Do i =43,89 
          If(Yqse(i)>0.and.y(i)>0)then
          sumsi=sumsi+ y(i)*aa(i)
           else
           dummycount=dummycount+1
         Endif
       Enddo
         rsumsi=sumsi/70
       Do i=89,ny  
          If(Yqse2(i)>0.and.y(i)>0)then
          sumfe=sumfe+ y(i)*aa(i)
          else 
          dummycount=dummycount+1
          Endif
       enddo
          rsumfe=sumfe/152
        write(987,*)"sumsi",rsumsi, "sumfe", rsumfe, dummycount
   
    If(t>9.8e-11.and.t<1.2e-10)then
          do i =1,ny
          write(990,'(i4,3es14.7)') aa1(i),Yqse(i)/y(i),Yqse2(i)/y(i),Ynse(i)/y(i)
          enddo
          write(989,*) t, sumsi,sumfe
        Endif
  
    If(t>9.8e-10.and.t<1.2e-9)then
          do i =1,ny
          write(991,'(i4,3es14.7)') aa1(i),Yqse(i)/y(i),Yqse2(i)/y(i),Ynse(i)/y(i)
          enddo
          write(989,*) t, sumsi,sumfe
        Endif
     If(t>9.8e-9.and.t<1.2e-8)then
          do i =1,ny
          write(992,'(i4,3es14.7)') aa1(i),Yqse(i)/y(i),Yqse2(i)/y(i),Ynse(i)/y(i)
          enddo
          write(989,*) t, sumsi,sumfe
        Endif

          If(t>9.6e-8.and.t<1.4e-7)then
          do i =1,ny
          write(993,'(i4,3es14.7)') aa1(i),Yqse(i)/y(i),Yqse2(i)/y(i),Ynse(i)/y(i)
          enddo
          write(989,*) t, sumsi,sumfe
        Endif

  If(t>9.4e-7.and.t<1.4e-6)then
          do i =1,ny
          write(994,'(i4,3es14.7)') aa1(i),Yqse(i)/y(i),Yqse2(i)/y(i),Ynse(i)/y(i)
          enddo
          write(989,*) t, sumsi,sumfe
        Endif
       If(t>9.8e-6.and.t<1.2e-5)then
          do i =1,ny
          write(995,'(i4,3es14.7)') aa1(i),Yqse(i)/y(i),Yqse2(i)/y(i),Ynse(i)/y(i)
          enddo
          write(989,*) t, sumsi,sumfe
        Endif

       If(t>9.8e-5.and.t<1.3e-4)then
          do i =1,ny
          write(996,'(i4,3es14.7)') aa1(i),Yqse(i)/y(i),Yqse2(i)/y(i),Ynse(i)/y(i)
          enddo
          write(989,*) t, sumsi,sumfe
        Endif
  
       If(t>9.8e-4.and.t<1.3e-3)then
          do i =1,ny
          write(997,'(i4,3es14.7)') aa1(i),Yqse(i)/y(i),Yqse2(i)/y(i),Ynse(i)/y(i)
          enddo
          write(989,*) t, sumsi,sumfe
        Endif
       If(t>9.8e-3.and.t<1.3e-2)then
          do i =1,ny
          write(998,'(i4,3es14.7)') aa1(i),Yqse(i)/y(i),Yqse2(i)/y(i),Ynse(i)/y(i)
          enddo
          write(989,*) t, sumsi,sumfe
        Endif
       If(t>9.9e-2.and.t<1.1e-1)then
          do i =1,ny
          write(999,'(i4,3es14.7)') aa1(i),Yqse(i)/y(i),Yqse2(i)/y(i),Ynse(i)/y(i)
          enddo
          write(989,*) t, sumsi,sumfe
        Endif
            Endif !!!!!QSE itsout if
	Return
End Subroutine ts_output
  
Subroutine final_output(kstep)
!===============================================================================
! full_net Calls this routine to handle output at the End of its
! execution.  
!===============================================================================
  Use controls
  Use conditions
  Use nuclear_data
  Use thermo_data
  Use abundances
  Use match_data
  Use flux_data
  Integer, Dimension(1) :: Iflx_mx
  Integer :: i,k,kstep
  Real(8), Dimension(0:ny) :: dyf
  
! Write final abundances to diagnostic output (in ASCII)
  Write(lun_diag,"(a3,2i6,2es14.7)") 'End',kstep,kstmx,tt,tstop
  Write(lun_diag,"(4(a5,es14.7,1x))") (nname(i),aa(i)*y(i),i=1,ny)
  
! Write integrated fluxes to diagnotic output
! Iflx_mx=maxloc(abs(flx_int))
! Write(lun_diag,'(2x,a8,i5,es11.3,8a5)') 'Flux Max',Iflx_mx(1),flx_int(Iflx_mx(1)),&
!&   nname(nflx(1:3,Iflx_mx(1))),' <-> ',nname(nflx(4:7,Iflx_mx(1)))
! Write(lun_diag,'(i5,8a5,i5,es11.3)') (k,nname(nflx(1:3,k)),' <-> ',&
!&   nname(nflx(4:7,k)),iwflx(k),flx_int(k),k=1,mflx)
  
! Write flux output formatted for Matlab
  Write(lun_diag,'(i5,4f6.1,es13.5,a5)') &
&   (k,zz(Ifl_orig(k)),nn(Ifl_orig(k)),zz(Ifl_term(k)),nn(Ifl_term(k)),flx_int(k),descx(k),k=1,mflx)
  
! Test how well sums of fluxes match abbundances changes
! dyf=0.0
! Do k=1,mflx 
!   dyf(nflx(1:3,k))=dyf(nflx(1:3,k))+flx_int(k)
!   dyf(nflx(4:7,k))=dyf(nflx(4:7,k))-flx_int(k)
! EndDo
! dy=y-yo+dyf(1:ny)
! Write(lun_diag,'(a5,es10.3)') 'DYF',dyf(0)
! Write(lun_diag,'(a5,4es10.3)') (nname(k),dyf(k),dy(k),y(k),yo(k), k=1,ny)
  Return
End Subroutine final_output
  
Subroutine sum_test(y)
  Use nuclear_data
  Real(8), Dimension(ny) :: y
  Real(8), Dimension(ny,2) :: stest
  Real(8), Dimension(ny,2) :: xtot
  stest(:,1)=y
  stest(:,2)=10.0*y
! xtot=sum((aa*stest),dim=2)
! xtot(:,2)=aa*stest(:,2)
  Write(6,*) xtot
  Return
End Subroutine sum_test
  
