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
  Use nuclear_data
  Use match_data
  Use flux_data
  Use timers
! Use nse_abundance                                                     !NSE
! Use neutrino_data                                                     !NNU
!$ Use omp_lib
  Integer :: i,j,k,n,nz,izone! Loop indices
  Integer :: kstep,ii(14),ierr,index,lun_control
  Integer :: nstart,kstart
  Real(8) :: tdelstart,t9start,rhostart
  Integer, Parameter :: nzmx=1000
  Real(8), Dimension(:), Allocatable :: yin
  Real(8) :: dt,rdt,dye,yestart,ytot,abar,zbar,z2bar,zibar              !NSE
  Real(8), Dimension(:), Allocatable, Save :: dyf,flx_diff
  Character (LEN=5)  :: output_nuc(14)
  Character (LEN=80) :: descript(3),data_desc
  Character (LEN=80) :: data_dir,thermo_file(nzmx),inab_file(nzmx)
  Character (LEN=80) :: ev_file_base,bin_file_base,diag_file_base='net_diag'
  Character (LEN=80) :: ev_file,bin_file,diag_file
  Character (LEN=80) :: abund_desc,thermo_desc
  
! Identify threads
!$OMP PARALLEL DEFAULT(SHARED)
  mythread=1
!$ mythread=OMP_get_thread_num()      
!$OMP SINGLE
  nthread=1
!$ nthread=OMP_get_num_threads()    
!$OMP END SINGLE
!$OMP END PARALLEL

! Initiate setup timer
  start_timer = xnet_wtime()
  timer_setup = timer_setup - start_timer  

!-----------------------------------------------------------------------
! The control file contains the parameters which determine the actions 
! of XNet.  
!-----------------------------------------------------------------------
  lun_control=5
  Open(lun_control,FILE='control')                                         

! Read Problem Description
  call find_controls_block(lun_control,'Problem Description',ierr)  
  Read(lun_control,"(a80)") (descript(i),i=1,3) ! text description of the problem.
  
! Read Job Controls
  call find_controls_block(lun_control,'Job Controls',ierr)  
  Read(lun_control,*) szone        ! number of the zone with which to begin
  Read(lun_control,*) nzone        ! total # of zones
  Read(lun_control,*) iweak0       ! controls the treatment of weak reactions 
  Read(lun_control,*) iscrn        ! controls the treatment of nuclear screening
  Read(lun_control,*) iprocess     ! controls the runtime pre-processing of the network data 

! Read Integration Controls
  call find_controls_block(lun_control,'Integration Controls',ierr)  
  Read(lun_control,*) isolv        ! Choice of integrations scheme
  Read(lun_control,*) kstmx        ! max # of timesteps for each zone
  Read(lun_control,*) kitmx        ! max # of iterations before retry
  Read(lun_control,*) ijac         ! rebuild jacobian every ijac iterations after the first
  Read(lun_control,*) iconvc       ! determines which convergence condition is Used
  Read(lun_control,*) changemx     ! allowed abundance change Used to set the timestep.
  Read(lun_control,*) yacc         ! abundances > yacc Used for timestep calculation
  Read(lun_control,*) tolm         ! mass conservation convergence criterion
  Read(lun_control,*) tolc         ! convergence limit on the iterative abundance change
  Read(lun_control,*) ymin         ! abundance < ymin is set to 0.0
  Read(lun_control,*) tdel_maxmult ! max factor by which the timestep is changed

! Read Self-heating Controls
  call find_controls_block(lun_control,'Self-heating Controls',ierr)  
  If(ierr/=0) Then
    Read(lun_control,*) iheat      ! controls the coupling of the network to temperature
    Read(lun_control,*) changemxt  ! allowed temperature change used to set the timestep.
    Read(lun_control,*) tolt9      ! convergence limit on the iterative temperature change
  Else
    Write(6,*) 'Using Default Self-heating behavior'
    iheat = 0
    changemxt = 1.0e-2
    tolt9 = 1.0e-4
  EndIf

! Read NSE Initial Abundance Controls 
! temperature at which NSE initial conditions are used, t9nse =8 is default.
! call find_controls_block(lun_control,'NSE Initial Conditions',ierr)   !NSE
! If(ierr/=0) Then                                                      !NSE
!     Read(lun_control,*) t9nse                                         !NSE
! Else                                                                  !NSE
!   Write(6,*) 'Using Default NSE behavior'                             !NSE 
!   t9nse = 8.0                                                         !NSE
! Endif                                                                 !NSE

! Read Neutrino Controls
! call find_controls_block(lun_control,'Neutrinos',ierr)                !NNU
! If(ierr/=0) Then                                                      !NNU
!   Read(lun_control,*) ineutrino   ! controls neutrino reactions       !NNU
! Else                                                                  !NNU
!   Write(6,*) 'Using Default Neutrino behavior'                        !NNU 
!   ineutrino=0                                                         !NNU
! Endif                                                                 !NNU

!-------------------------------------------------------------------------------
! XNet output controls include the base of the filenames to which ASCII and 
! binary output are written at the end of each timestep, and a subset of nuclei
! to be included in the per timestep ASCII output file.
!-------------------------------------------------------------------------------
! Read Output Controls
  call find_controls_block(lun_control,'Output Controls',ierr)  
  Read(lun_control,*) idiag        ! sets diagnostic output level
  Read(lun_control,*) itsout       ! sets per timestep output level
  Read(lun_control,"(72x)")
  Read(lun_control,"(a80)") ev_file_base            
  Read(lun_control,"(72x)")
  Read(lun_control,"(a80)") bin_file_base            
  Read(lun_control,"(72x)")
  Read(lun_control,"(14a5)") output_nuc         

!-------------------------------------------------------------------------------
! XNet input controls include the relative directory from which the nuclear data 
! should be loaded, as well as the names of the files containing the initial 
! abundances and thermodynamic trajectories.
!-------------------------------------------------------------------------------
! Read Input Controls
  call find_controls_block(lun_control,'Input Controls',ierr)  
  Read(lun_control,"(72x)")
  Read(lun_control,"(a80)") data_dir
  Read(lun_control,"(72x)")
  Do nz=1,nzone
    Read(lun_control,"(a80)",IOSTAT=ierr) inab_file(nz)
    Read(lun_control,"(a80)",IOSTAT=ierr) thermo_file(nz)
    If(ierr < 0) Then
      Exit
    ElseIf(ierr > 0) Then
      Write(6,*) 'Problem reading Input Filenames'
    Endif
  EndDo
  If(nz<nzone) Write(6,*) nz,' datafiles for ',nzone,' zones!'
! Write(6,*) inab_file(nz),thermo_file(nz)
  Close(lun_control)
  
!-----------------------------------------------------------------------
! Output files:
! The diagnostic output is sorted by thread, lun_diag = 50+mythread
!
! General output is per zone, however logical units are spaced per thread
! lun_ev = 50+nthread+myid ; lun_ts ; unit 50+2*nthread+myid     
!-----------------------------------------------------------------------
! Open diagnositic output file, per thread if OMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(diag_file)
  If(idiag>=0) Then
    lun_diag=50+mythread
    diag_file=trim(diag_file_base)
    Call name_ordered(diag_file,mythread,nthread)
    Open(lun_diag,file=diag_file)
    FLUSH(lun_diag)
!$  Write(lun_diag,"(a,i4,a,i4)") 'Thread ',mythread,' of ',nthread
  EndIf

! Set iweak to original value
  iweak=iweak0
!$OMP End PARALLEL
  
! Retain the value of iweak
! iweak0=iweak
  
! In requested, pre-process the nuclear and reaction data.
  If(iprocess>0) Call net_preprocess(6,data_dir,data_dir)

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

! Initialize EoS for screening
  If(iscrn>0.or.iheat>0) call eos_initialize
  
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
  
! Stop setup timer
  stop_timer = xnet_wtime()
  timer_setup = timer_setup + stop_timer
  
! Set sizes of abundance arrays
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP   PRIVATE(yin,dyf,flx_diff,kstart,nstart,tdelstart,t9start, &
!$OMP     rhostart,yestart,dt,rdt,dye,abund_desc,thermo_desc, &
!$OMP     ev_file,bin_file,izone,n,ierr,ytot,abar,zbar,z2bar,zibar) &
!$OMP   COPYIN(timer_setup)

! Start setup timer
  start_timer = xnet_wtime()
  timer_setup = timer_setup - start_timer

  Allocate (y(ny),yo(ny),yt(ny),ydot(ny),yin(ny),dyf(0:ny),flx_diff(ny))
  
!-----------------------------------------------------------------------
! The initial abundances and thermo-data are read into a series of local
! arrays that are functions of the zone. All of the variables associated 
! with the local arrays End in "L". The local arrays are then loaded by 
! a loop over zones into global varaibles and dIfferent threads according 
! to zone in the OMP parallel region.  The local arrays are Allocated to 
! have the zone size set to the maximum number of zones, nzmx. 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Loop over all zones, reading each zone and calling full_net to perform 
! each integration
!-----------------------------------------------------------------------

! call nse_initialize                                               !NSE                                
  
! Stop setup timer
  stop_timer = xnet_wtime()
  timer_setup = timer_setup + stop_timer

!$OMP DO
  Do izone=szone,nzone

! Start setup timer
    start_timer = xnet_wtime()
    timer_setup = timer_setup - start_timer

! Load the zone ID quadruplet
    zone_id(1)=izone
    zone_id(2)=1
    zone_id(3)=1
    zone_id(4)=1
  
! Read the thermdynamic trajectory
    lun_th=50+4*nzmx+izone
    Open(lun_th,file=thermo_file(izone))
    Read(lun_th,"(a)") thermo_desc
    Read(lun_th,*) tstart
    Read(lun_th,*) tstop
    Read(lun_th,*) tdelstart
  
!   yeh = 0.0                                                       !NSE
    Do n=1,nhmx
      Read(lun_th,*,IOSTAT=ierr) th(n),t9h(n),rhoh(n) !NOTNSE !NOTNNU
!     Read(lun_th,*,IOSTAT=ierr) th(n),t9h(n),rhoh(n),yeh(n) !NSE !NOTNNU
!     Read(lun_th,*,IOSTAT=ierr) th(n),t9h(n),rhoh(n),yeh(n),fluxcms(n,:),tmevnu(n,:) !NNU

      If(ierr<0) Then
        If(idiag>2) Write(lun_diag,"(a,i6,a)") 'End of Thermo File Reached after',n,' records'
        Exit
      EndIf
    EndDo
    nh=n-1
    Close(lun_th)
    
! Convert to appropriate units (CGS, except temperature (GK) and neutrino flux)
!   t9h=t9h/1e9
!   fluxcms=1.d-42*fluxcms                                              !NNU

! Determine thermodynamic conditions at tstart.
    Call t9rhofind(0,tstart,nstart)
    t9start=t9t ; rhostart =rhot    
    If(idiag>0) Write(lun_diag,"(a,i6,a,f6.3,a,es10.3)") &
 &  'Start',nstart,' T9=',t9start,' Rho=',rhostart

! Load initial Abundances.    
! For High temperatures, use NSE initial abundance.
!   If(t9start>t9nse .and. yeh(1)>0.0d0 .and. yeh(1)<=1.0d0) Then       !NSE

! Interpolate initial Ye
!     If(nstart>1.and.nstart<=nh) Then                                  !NSE
!       rdt=1.0/(th(nstart)-th(nstart-1))                               !NSE
!       dt=tstart-th(nstart-1)                                          !NSE
!       dye=yeh(nstart)-yeh(nstart-1)                                   !NSE
!       yestart=dt*rdt*dye+yeh(nstart-1)                                !NSE
!     ElseIf(nstart==1) Then                                            !NSE
!       yestart=yeh(1)                                                  !NSE
!     Else                                                              !NSE
!       yestart=yeh(nh)                                                 !NSE
!     EndIf                                                             !NSE

! Calculate NSE abundances
!     Write(lun_diag,fmt='(a,es10.4,a,es10.4,a,f5.4)') &                !NSE 
!&      'NSE abundances for T9=',t9start,', rho=',rhostart,', Ye=',yestart !NSE
!     call nse_descend(rhostart,yestart,t9start,t9start)                !NSE
!     yin=ynse                                                          !NSE

! Read Initial Abundances
!   Else                                                                !NSE
      lun_ab=50+3*nzmx+izone
      Open(lun_ab,file=inab_file(izone))
      Read(lun_ab,"(a)") abund_desc
      Read(lun_ab,"(4(5x,es14.7,1x))") (yin(i),i=1,ny)
      Close(lun_ab)
      If(idiag>=0) Then
! Log abundance file and description
        Write(lun_diag,"(a)") inab_file(izone)
        Write(lun_diag,"(a)") abund_desc
      Endif
      call y_moment(yin,yestart,ytot,abar,zbar,z2bar,zibar)
!     If( t9start>t9nse ) Then                                          !NSE
! Calculate NSE abundances
!       Write(lun_diag,fmt='(a,es10.4,a,es10.4,a,f5.4)') &              !NSE
!       & 'NSE abundances for T9=',t9start,', rho=',rhostart,', Ye=',yestart !NSE
!       call nse_descend(rhostart,yestart,t9start,t9start)              !NSE
!       yin=ynse                                                        !NSE
!     EndIf                                                             !NSE
!   Endif                                                               !NSE
!   Call sum_test(yin)
  
! Load initial abundances, time and timestep
    kstart=1
    tdel= 0.0
    y   = yin
    yet = yestart
    t9  = t9t
 
! Log initial abundance.
    If(idiag>=0) Then
      Write(lun_diag,"(5(a6,1es10.3))") (nname(i),yin(i),i=1,ny)
    Endif
                   
 ! Log thermo description
    If(idiag>=0) Write(lun_diag,"(a)") thermo_desc
                   
! Open the evolution file
    lun_ev=50+nzmx+izone
    If(itsout>=2) Then
      ev_file=ev_file_base
      Call name_ordered(ev_file,izone,nzone)
      If(idiag>=0) Write(lun_diag,"(a,i5,7es10.3)") trim(ev_file),&
&       nh,th(nh),t9h(nh),rhoh(nh),tstart,tstop
      Open(lun_ev,file=ev_file)
      ii=(/5,21,29,36,43,53,63,73,89,102,115,117,129,144/)
      Write(lun_ev,"(a4,a15,4a10,15a9,a4)") &
&       'k ',' Time ',' T(GK) ',' Density ',' dE/dt ',' Timestep ',nname(inout), ' It '
    EndIf
  
! Open the time series file
    lun_ts=50+2*nzmx+izone 
    If(itsout>=1) Then
      bin_file=bin_file_base
      Call name_ordered(bin_file,izone,nzone)
      Open(lun_ts,file=bin_file,form='unformatted')
! Write Control Parameters to ts file
      Write(lun_ts) (descript(i),i=1,3),data_desc
      Write(lun_ts) kstmx,kitmx,iweak,iscrn,iconvc,changemx,tolm,tolc,yacc,ymin,tdel_maxmult
! Write abundance description to ts file
      Write(lun_ts) inab_file(izone),abund_desc
! Write thermo description to ts file
      Write(lun_ts) thermo_file(izone),thermo_desc
! Write species identifiers to ts file
      Write(lun_ts) ny,zz,aa
! Write flux identifiers to ts file
      Write(lun_ts) mflx,ifl_orig,ifl_term
    EndIf

! Stop setup timer
    stop_timer = xnet_wtime()
    timer_setup = timer_setup + stop_timer
  
! Evolve abundance from tstart to tstop, using the original iweak
    iweak=iweak0
    Call full_net

! Test how well sums of fluxes match abundances changes
    If(idiag>=3) Then 
      dyf=0.0
      Do k=1,mflx 
       Do j=1,3 
        dyf(nflx(j,k))=dyf(nflx(j,k))+flx_int(k)
      Enddo
      Do j=4,7
        dyf(nflx(j,k))=dyf(nflx(j,k))-flx_int(k)
      Enddo
      EndDo
      flx_diff=y-yin+dyf(1:ny)
      Write(lun_diag,'(a,es11.3)') 'Compare Integrated flux to abundance change',dyf(0)
      Write(lun_diag,'(a)') 'Species Flux Sum + Y Final - Y Initial = Flux Diff'
      Write(lun_diag,'(a5,4es11.3)') (nname(k),dyf(k),y(k),yin(k),flx_diff(k), k=1,ny)
    Endif
   
    Close(lun_ev)
    Close(lun_ts)
  EndDo
  
!$OMP End DO
!$OMP End PARALLEL
  
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
  Use timers
  Real(8) :: xg(7),enuc,edot
  Real(8), Dimension(ny) :: yout
  Integer :: i,j,k,ig,kstep,kout
  
! Initiate output timer
  start_timer = xnet_wtime()
  timer_output = timer_output - start_timer  
  
  yout(:) = y(:)
  
! Calculate reaction fluxes
  If(kstep>0) Then
    Call flux
  Else
    flx_int=0.0
  EndIf
  
! An abundance snapshot is written to the binary file on unit 24.
  Write(lun_ts) kstep,t,t9t,rhot,tdel,edot,y,flx
  
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
  
! For itsout>=3, output time and thermo evolution to the screen
    If(itsout>=3) Write(6,"(i5,4es12.4,2i3)") kstep,t,tdel,t9t,rhot,kmon
  
! For itsout>=4, output abundances for each timestep to the diagnostic file
    If(itsout>=4) Write(lun_diag,"(4(a5,es14.7,1x))") (nname(i),aa(i)*y(i),i=1,ny)
  
! For itsout>=5, output fluxes for each timestep to the diagnostic file
    If(itsout>=5) Write(lun_diag,'(i5,8a5,i5,es11.3)') &
&     (k,nname(nflx(1:3,k)),' <-> ',nname(nflx(4:7,k)),iwflx(k),flx(k),k=1,mflx)
  EndIf
  
! Stop output timer
  stop_timer = xnet_wtime()
  timer_output = timer_output + stop_timer

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
  Use timers
  Integer, Dimension(1) :: iflx_mx
  Integer :: i,k,kstep
  
! Initiate output timer
  start_timer = xnet_wtime()
  timer_output = timer_output - start_timer  

! Write final abundances to diagnostic output (in ASCII)
  Write(lun_diag,"(a3,2i6,2es14.7)") 'End',kstep,kstmx,tt,tstop
  Write(lun_diag,"(4(a5,es14.7,1x))") (nname(i),aa(i)*y(i),i=1,ny)
  Write(lun_diag,"(a3,3i6,4es10.3)") 'CNT',isolv,iconvc,kitmx,yacc,tolm,tolc,changemx
  
! Write integrated fluxes to diagnotic output
  Iflx_mx=maxloc(abs(flx_int))
  Write(lun_diag,'(2x,a8,i5,es11.3,8a5)') 'Flux Max',Iflx_mx(1),flx_int(Iflx_mx(1)),&
 &   nname(nflx(1:3,Iflx_mx(1))),' <-> ',nname(nflx(4:7,Iflx_mx(1)))
  Write(lun_diag,'(i5,8a5,i5,es11.3)') (k,nname(nflx(1:3,k)),' <-> ',&
 &   nname(nflx(4:7,k)),iwflx(k),flx_int(k),k=1,mflx)
  
! Write flux output formatted for Matlab
! Write(lun_diag,'(i5,4f6.1,es13.5,a5)') &
!& (k,zz(ifl_orig(k)),nn(ifl_orig(k)),zz(ifl_term(k)),nn(ifl_term(k)),flx_int(k),descx(k),k=1,mflx)
  
! Stop output timer
  stop_timer = xnet_wtime()
  timer_output = timer_output + stop_timer

! Write timers
  Write(lun_diag,"(a8,8a10)") 'Timers: ','Total','Solver','Jacobian','Deriv', 'CrossSect', 'Screening', 'Setup','Output'
  Write(lun_diag,"(8x,8es10.3)") timer_total,timer_solve,timer_jacob,timer_deriv,timer_csect,timer_scrn, timer_setup,timer_output
  Write(lun_diag,"(a8,2a10)") ' Counts: ','Timestep','Iteration'
  Write(lun_diag,"(8x,2i10)") ktot(1),ktot(2)
! Uncomment the following line to restart the timers  
! timer_total =0.0 ; timer_solve =0.0 ; timer_jacob =0.0 ; timer_deriv =0.0 ; timer_csect =0.0 ; timer_scrn =0.0 ;  timer_setup =0.0 

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
  
