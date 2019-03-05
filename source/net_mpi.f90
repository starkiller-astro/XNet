!***********************************************************************
! net MPI (for XNet v6)
!
! This program is a driver to run XNet on multiple processors using MPI, 
! by allocating successive zones to each processor
!***********************************************************************

Program net_mpi
!=======================================================================
! This is the shell routine for running the network under MPI.  
!=======================================================================
  Use controls
  Use file_data
  Use conditions
  Use abundances
  Use thermo_data
  Use nuclear_data
  Use match_data
  Use flux_data
  Use timers
  Use mpi
! Use nse_abundance                                                     !NSE
! Use neutrino_data                                                     !NNU
!$ Use omp_lib
  Integer :: i,j,k,n,izone
  Integer :: kstep,ierr,index
  Integer :: nstart,kstart

! Input Data descriptions
  Character (LEN=80) :: data_dir,data_desc
  Character (LEN=80) :: abund_desc,thermo_desc

! Thermodynamic input data
  Real(8) :: tdelstart,t9start,rhostart
  Real(8) :: dt,rdt,dye,yestart,ytot,abar,zbar,z2bar,zibar              !NSE
  Real(8), Dimension(:), Allocatable :: dyf,flx_diff

! Input abundance data
  Real(8), Dimension(:), Allocatable :: yin

! Filenames
  Character (LEN=80) :: ev_file,bin_file,diag_file
  Character (LEN=80) :: diag_file_base='mp_diag'
  Character (LEN=80) :: inab_file,thermo_file

! Combined Output data
! Real(8), Dimension(:,:), Allocatable :: yfinal,flxfinal
  
! Identify number of MPI nodes and ID number for this PE
  call mpi_init(ierr)
  call mpi_comm_rank( MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size( MPI_COMM_WORLD, nproc, ierr)

! Identify OpenMP threads
!$OMP PARALLEL DEFAULT(SHARED)
  mythread=1
!$ mythread=OMP_get_thread_num()      
!$OMP SINGLE
  nthread=1
!$ nthread=OMP_get_num_threads()    
!$OMP End SINGLE
!$OMP End PARALLEL

! Initiate setup timer
  start_timer = xnet_wtime()
  timer_setup = timer_setup - start_timer

! Open Diagnositic file
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(diag_file)
  lun_diag=50+mythread
  diag_file=trim(diag_file_base)
  call name_ordered(diag_file,myid,nproc)
  If(nthread>1) call name_ordered(diag_file,mythread,nthread)
  Open(lun_diag,file=diag_file)
  FLUSH(lun_diag)
  Write(lun_diag,"(a5,4i5)") 'MyId',myid,nproc,mythread,nthread
!$OMP End PARALLEL

! Read and distribute control file data
  call control_bcast(data_dir)

! Read in and distribute nuclear and reaction data
  call netdata_bcast(data_dir,data_desc)
  If(idiag>=0) Write(lun_diag,"(a)") (descript(i),i=1,3),data_desc

! Read in and distribute data for reaction matching and fluxes
  call match_bcast(data_dir)

! BCast jacobian data
  call jacobian_bcast(data_dir)
  
! Initialize EoS for screening
  If(iscrn>0.or.iheat>0) call eos_initialize

! Convert output_nuc names into indices
  Do i=1,14
    Call index_from_name(output_nuc(i),index)
    If(index<1.or.index>ny) Then
      Write(lun_diag,*) 'Output Nuc:',i,output_nuc(i),' not found'
      inout(i)=ny
    Else
      inout(i)=index
    EndIf
  EndDo

! Stop setup timer
  stop_timer = xnet_wtime()
  timer_setup = timer_setup + stop_timer

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP   PRIVATE(yin,dyf,flx_diff,kstart,nstart,tdelstart,t9start, &
!$OMP     rhostart,yestart,dt,rdt,dye,abund_desc,thermo_desc, &
!$OMP     ev_file,bin_file,izone,n,ierr,inab_file,thermo_file, &
!$OMP     ytot,abar,zbar,z2bar,zibar) &
!$OMP   COPYIN(timer_setup)

! Set sizes of abundance arrays
  Allocate (y(ny),yo(ny),yt(ny),ydot(ny),yin(ny),dyf(0:ny),flx_diff(ny))

! call nse_initialize                                               !NSE                                

! Loop over zones assiging zones to MPI tasks in order
!$OMP DO
  
  Do izone=myid+1,nzone,nproc

! Load the zone ID quadruplet
    zone_id(1)=izone
    zone_id(2)=1
    zone_id(3)=1
    zone_id(4)=1

! Build thermo and input abundance filenames
    thermo_file=trim(thermo_file_base)
    Call name_ordered(thermo_file,izone,nzone)
    inab_file=trim(inab_file_base)
    Call name_ordered(inab_file,izone,nzone)
  
! Read the thermdynamic trajectory
!   yeh(1) = 0.0                                                    !NSE
    lun_th=50+4*nzone+izone
    Open(lun_th,file=thermo_file)
    Read(lun_th,"(a)") thermo_desc
    Read(lun_th,*) tstart
    Read(lun_th,*) tstop
    Read(lun_th,*) tdelstart
    Do n=1,nhmx
      Read(lun_th,*,IOSTAT=ierr) th(n),t9h(n),rhoh(n) !NOTNSE !NOTNNU
!     Read(lun_th,*,IOSTAT=ierr) th(n),t9h(n),rhoh(n),yeh(n) !NSE !NOTNNU
!     Read(lun_th,*,IOSTAT=ierr) th(n),t9h(n),rhoh(n),yeh(n),fluxcms(n,:),tmevnu(n,:) !NNU
      If(ierr<0) Then
!       Write(*,*) 'End of File Reached after',i,' records'
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
!   If(t9start>t9nse .AND. yeh(1)>0.0d0 .AND. yeh(1)<=1.0d0) Then       !NSE

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
!     Write(lun_diag,fmt='(a,es10.4,a,es10.4,a,f5.4)') &                    !NSE 
!&      'NSE abundances for T9=',t9start,', rho=',rhostart,', Ye=',yestart  !NSE
!     call nse_descend(rhostart,yestart,t9start,t9start)                    !NSE
!     yin=ynse                                                              !NSE

! Read Initial Abundances
!   Else                                                                    !NSE
      lun_ab=50+3*nzone+izone
      Open(lun_ab,file=inab_file)
      Read(lun_ab,"(a)") abund_desc
      Read(lun_ab,"(4(5x,es14.7,1x))") (yin(i),i=1,ny)
      Close(lun_ab)
! Log abundance file and description
      If(idiag>=0) Write(lun_diag,"(a)") inab_file
      If(idiag>=0) Write(lun_diag,"(a)") abund_desc
      call y_moment(yin,yestart,ytot,abar,zbar,z2bar,zibar)
!     If( t9start>t9nse ) Then                                          !NSE
!         ! Calculate NSE abundances
!         Write(lun_diag,fmt='(a,es10.4,a,es10.4,a,f5.4)') &            !NSE 
!         & 'NSE abundances for T9=',t9start,', rho=',rhostart,', Ye=',yestart !NSE
!         call nse_descend(rhostart,yestart,t9start,t9start)            !NSE
!         yin=ynse                                                      !NSE
!     EndIf                                                             !NSE
!   Endif                                                               !NSE

! Log thermo description
    If(idiag>=0) Write(lun_diag,"(a)") thermo_desc
  
! Load initial abundances, time and timestep
    kstart=1
    tdel=0.0
    y(:)   = yin(:)
    yet    = yestart
    t9     = t9t
  
! Open the evolution file
    lun_ev=50+nzone+izone
    If(itsout>=2) Then
      ev_file=ev_file_base
      Call name_ordered(ev_file,izone,nzone)
      If(idiag>=0) Write(lun_diag,"(a,i5,7es10.3)") trim(ev_file),&
&       nh,th(nh),t9h(nh),rhoh(nh),tstart,tstop
      Open(lun_ev,file=ev_file)

! Write evolution file header
      Write(lun_ev,"(a4,a15,4a10,15a9,a4)") &
&       'k ',' Time ',' T(GK) ',' Density ',' dE/dt ',' Timestep ',nname(inout), ' It '
    EndIf
  
! Open the time series file
    lun_ts=50+2*nzone+izone
    If(itsout>=1) Then
      bin_file=bin_file_base
      Call name_ordered(bin_file,izone,nzone)
      Open(lun_ts,file=bin_file,form='unformatted')
! Write Control Parameters to ts file
      Write(lun_ts) (descript(i),i=1,3),data_desc
      Write(lun_ts) kstmx,kitmx,iweak,iscrn,iconvc,changemx,tolm,tolc,yacc,ymin,tdel_maxmult
! Write abundance description to ts file
      Write(lun_ts) inab_file,abund_desc
! Write thermo description to ts file
      Write(lun_ts) thermo_file,thermo_desc
! Write species identifiers to ts file
      Write(lun_ts) ny,zz,aa
! Write flux identifiers to ts file
      Write(lun_ts) mflx,ifl_orig,ifl_term
    EndIf
  
! Call network routine
    iweak=iweak0
    call full_net

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

! Save abundances and reactions fluxes
!   Write(lun_diag,"(4(a5,es14.7,1x))") (nname(i),y(i),i=1,ny)
!   yfinal(:,iz)=y
!   flxfinal(:,iz)=flx_int

! Close zone output files
    If(itsout>=2) Close(lun_ev)
    If(itsout>=1) Close(lun_ts)

  Enddo

!$OMP End DO
!$OMP End PARALLEL

! Write combined output
! Write(30) itrial,yfinal(:,1:nzone),flxfinal(:,1:nzone)

! Perform Network Calculations
! call net_execute

!$OMP PARALLEL DEFAULT(SHARED)
  Close(lun_diag)
!$OMP End PARALLEL

! Wait for all nodes to finish
  call mpi_finalize(ierr)
  Stop 'driver'
  End

Subroutine ts_output(kstep,enuc,edot) 
!=======================================================================
! If the flag itso is > 0, full_net calls this routine to handle stepwise 
! output.  
!=======================================================================
  Use controls
  Use nuclear_data
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
  iflx_mx=maxloc(abs(flx_int))
  Write(lun_diag,'(2x,a8,i5,es11.3,8a5)') 'Flux Max',iflx_mx(1),flx_int(iflx_mx(1)),&
 &   nname(nflx(1:3,iflx_mx(1))),' <-> ',nname(nflx(4:7,iflx_mx(1)))
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
  xtot(:,2)=aa*stest(:,2)
  Write(6,*) xtot
  Return
End Subroutine sum_test
