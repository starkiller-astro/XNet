!*******************************************************************************
! This is an example driver for running FullNet on serial machines. 
!
! This driver reads in the controlling flags, the nuclear and 
! reaction data. It also loops through reading the initial abundances 
! and the thermodynamic trajectory before calling full_net for each zone.  
!
! This file also contains the ts_output routine, which controls the 
! information output at the end of each timestep, and the final_output 
! routine, which controls what is output an the end of the run. 
!
! If you wish to use full_net in concert with a hydrodynmics code, 
! you'll need to supply these services from within the hydro. 
!*******************************************************************************

      Program net
!===============================================================================
!===============================================================================
      use controls
      use conditions
      use abundances
      use thermo_data
      use nuclear_data
      integer :: i,j,k,n,nz,izone ! Loop indices
      integer :: szone,kstep,kstart,iweak0
      integer, parameter   :: nzmx= 100 ! Max number of zones
      real(8), dimension(:), allocatable :: yin
      real(8) :: tdelstart,tbegin,tend,tcut,t9min
      character (LEN=80) :: descript(3),data_desc,abund_desc,thermo_desc
      character (LEN=80) :: data_dir,thermo_file(nzmx),inab_file(nzmx)
      character (LEN=80) :: ev_file_base,ev_file,bin_file_base,bin_file
      character (LEN=5)  :: output_nuc(14)

! Open diagnositic output file
      Open(50,file='net_diag')

!-----------------------------------------------------------------------
! Read control file, which contains the parameters which control the 
! action of the network.  
!-----------------------------------------------------------------------
      Open(5,FILE='control')                                         
      Read(5,"(72x)")
      Read(5,"(a80)") (descript(i),i=1,3) ! text description of the problem.
      Read(5,"(72x)")
      Read(5,*) szone    ! number of the zone with which to begin
      Read(5,*) nzone    ! total # of zones
      Read(5,*) kstmx    ! max # of timesteps for each zone
      Read(5,*) knrmx    ! max # of Newton-Raphson iterations before retry
      Read(5,*) idiag    ! sets diagnostic output level
      Read(5,*) itso     ! sets per timestep output level
      Read(5,*) iweak    ! controls the treatment of weak reactions 
      Read(5,*) iscrn    ! controls the treatment of nuclear screening
      Read(5,*) iconvc  ! determines which convergence condition is used
      Read(5,*) changemx ! allowed abundance change used to set the timestep.
      Read(5,*) ytime    ! abundances > ytime used for timestep calculation
      Read(5,*) tolm     ! mass conservation convergence criterion
      Read(5,*) tolc    ! convergence limit on the iterative abundance change
      Read(5,*) ymin     ! abundance < ymin is set to 0.0
      Read(5,*) tdelmm   ! max factor by which the timestep is changed
!-----------------------------------------------------------------------
! control also indicates the base of the filenames to which ASCII and 
! binary output are written at the end of each timestep, the relative 
! directory from which the nuclear data should be loaded, as well as 
! the names of the files containing the initial abundances and \
! thermodynamic trajectories.
!-----------------------------------------------------------------------
      Read(5,"(72x)")
      Read(5,"(14a5)") output_nuc            
      Read(5,"(72x)")
      Read(5,"(a80)") ev_file_base            
      Read(5,"(72x)")
      Read(5,"(a80)") bin_file_base            
      Read(5,"(72x)")
      Read(5,"(a80)",END=10) data_dir
      Read(5,"(72x)")
      Do nz=1,nzone
        Read(5,"(a80)",END=10) inab_file(nz)
        Read(5,"(a80)",END=10) thermo_file(nz)
      Enddo
   10 nz=nz-1
      If((nz)/=nzone) Write(6,*) nz,' datafiles for ',nzone,' zones!'
      Write(6,*) inab_file(nz+10)

! Retain the value of iweak
      iweak0=iweak

! Read nuclear data (proton and atomic number, mass, partition functions)
      call read_nuclear_data(data_dir,data_desc)

! Read reaction rate data
      call read_reaction_data(data_dir)
      Write(50,"(a)") (descript(i),i=1,3),data_desc

! Read jacobian matrix data 
      call read_jacobian_data(data_dir)
 
! Read data on matching forward and reverse reactions 
      call read_match_data(data_dir)

! Initialize flux tracking (uncomment 1 line below, also need match)      
      call flux_init

! Convert output_nuc names into indices
      Do i=1,14
        call index_from_name(output_nuc(i),index)
        If(index<1.or.index>ny) Then
          Write(6,*) 'Output Nuc:',i,name,' not found'
          inout(i)=0
        Else
          inout(i)=index
        Endif
      Enddo

! Set minimum temperature for evolution of strong/EM reactions
      t9min=0.01

! Set sizes of abundance arrays
      Allocate (yo(ny),y(ny),yt(ny),ydot(ny),dy(ny),yin(ny))

! For each Zone
      Do izone=szone,nzone

! Open the output files
!       ev_file='ev'
        ev_file=ev_file_base
        call name_ordered(ev_file,izone,nzone)
        Open(52,file=ev_file)

        Write(52,"(a4,a15,4a10,14(2x,a5,2x),a4)") 
     &    'k',' Time   ',' T(GK) ',' Density ',' dE/dt ',
     &    ' Timestep',(nname(inout)), ' It'
!    &    ' Timestep ',
!    &    nname((/6,21,33,45,61,77,95,115,139,162,184,206,231,257/)), 
!    &    ' It '        
        bin_file='tso'
        bin_file=bin_file_base
        call name_ordered(bin_file,izone,nzone)
        Open(54,file=bin_file,form='unformatted')

! Write Control parameters to bin_file
        Write(54) (descript(i),i=1,3),data_desc
        Write(54) kstmx,knrmx,iweak,iscrn,iconvc,changemx,tolm,
     &    tolc,ytime,ymin,tdelmm

! Read Initial Abundances
        Open(8,file=inab_file(izone))
        kstart=1
        Read(8,"(a)") abund_desc
        Read(8,"(4(5x,es14.7,1x))") (yin(i),i=1,ny)
!       Read(8,"(5es15.7)") (y(i),i=1,ny)
        Close(8)
        Write(54) inab_file(izone),abund_desc
        Write(50,"(a)") abund_desc
!       call sum_test(yin)

!  Read the thermdynamic trajectory 
        tstart=0.0
        tcut=0.0
        t9old=0.0
        Open(8,file=thermo_file(izone))
        Read(8,"(a)") thermo_desc

!       Read(8,"(3(12x,1pd20.12))") tbegin,tend,tdelstart
        Read(8,*)tbegin
        Read(8,*)tend
        Read(8,*)tdelstart
        
        Do n=1,nhmx
!         Read(8,"(3es15.8)",end=20,err=20) th(n),t9,rhoh(n)
!         t9h(n)=t9/1.0e9
          Read(8,*,end=20,err=20) th(n),t9h(n),rhoh(n)
!  Find the temperature cutoffs
          If(t9h(n)<t9min) Then
            If(t9old<t9min) Then
!             tstart=th(n)
            Else 
              tcut=th(n)
            Endif
          Endif
          t9old=t9h(n)
          Write(50,'(i5,5es11.3)') n,th(n),t9h(n),rhoh(n),tcut,tstart
        Enddo
   20   nh=n-1
        Close(8)

! Log thermo description
        Write(54) thermo_file(izone),thermo_desc
        Write(50,"(a)") thermo_desc

! Load initial abundances, time and timestep
        tdel=0.0
        y=yin
        tstart=max(tbegin,tstart)
        If(tcut>0.0) Then
          tstop=min(tend,tcut)
        Else
          tstop=tend
        Endif

! Evolve abundance from tstart to tstop, using the original iweak
        iweak=iweak0
        call full_net(izone)

! For evolution at temperatures below T9=t9min, we use only weak reactions
        If(tstop<tend) Then
          iweak=-1   
          tstart=t   
          tstop=tend 
          call full_net(izone)
        Endif
!       yo=yin
        Close(52)
        Close(54)
      Enddo
      Stop
      End

      subroutine ts_output(kstep,enuc,edot) 
!===============================================================================
!  The per timestep output routine.  If the flag itso is > 0, 
!  full_net calls this routine to handle stepwise output.  
!===============================================================================
      use nuclear_data
      use controls
      use conditions
      use abundances
      use match_data
      use flux_data
      real(8) :: xg(7),enuc,edot,yout(0:ny)
      integer :: i,j,k,ig,kstep,kout

      yout(0)=0.0
      yout(1:ny)= y

!  Calculate reaction fluxes
      If(kstep>0) Then
         call flux
      Else
         flx_int=0.0
      Endif

! An abundance snapshot is written to the binary file on unit 24.
      Write(54) kstep,t,t9t,rhot,tdel,y,flx

! For itso>=2, output important mass fractions to the ASCII file on unit 22
      If(itso>=2) Then
        kout=10*kts+knr
        xg=0.0
!       Do i=1,ny
!         If(zz(i)<=1) Then
!           ig=1
!          Elseif(zz(i)<=2) Then
!            ig=2
!          Elseif(zz(i)<=8) Then
!            ig=3
!          Elseif(zz(i)<=10) Then
!            ig=4
!          Elseif(zz(i)<=12) Then
!            ig=5
!          Elseif(zz(i)<=14) Then
!            ig=6
!          Else
!            ig=7
!          Endif
!          xg(ig)=xg(ig)+aa(i)*y(i)
!        Enddo
!       Write(52,"(i4,1es15.8,2es10.3,2es10.2,8es9.2,i4)") 
!    &     kstep,t,t9t,rhot,edot,enuc,tdel,(xg(i),i=1,7),kout
!       ii=(/6,21,33,45,61,77,95,115,139,162,184,206,231,257/)
        Write(52,"(i4,1es15.8,2es10.3,2es10.2,14es9.2,i4)") 
     &    kstep,t,t9t,rhot,edot,tdel,(aa(inout)*yout(inout)),kout
!       Write(52,"(i4,1es15.8,2es10.3,2es10.2,14es9.2,i4)") 
!    &    kstep,t,t9t,rhot,edot,tdel,(aa*y),kout

! For itso>=3, output time and thermo evolution to the screen
        If(itso>=3) Write(6,"(i5,4es12.4,2i3)") 
     &    kstep,t,tdel,t9t,rhot,knr,kts

! For itso>=4, output abundances for each timestep to the diagnostic file
        If(itso>=4) Write(50,"(4(a5,es14.7,1x))") 
     &    (nname(i),aa(i)*y(i),i=1,ny)

! For itso>=5, output fluxes for each timestep to the diagnostic file
        If(itso>=5) Write(50,'(i5,8a5,i5,es11.3)')(k,nname(nflx(1:3,k)),
     &    ' <-> ',nname(nflx(4:7,k)),iwflx(k),flx(k),k=1,mflx)
      Endif
      Return
      End

      subroutine final_output(kstep)
!===============================================================================
! full_net calls this routine to handle output at the end of its
! execution.  
!===============================================================================
      use controls
      use conditions
      use nuclear_data
      use thermo_data
      use abundances
      use match_data
      use flux_data
      integer, dimension(1) :: iflx_mx
      integer :: i,kstep
      real(8), dimension(0:ny) :: dyf

! Write final abundances to diagnostic output (in ASCII)
      Write(50,"(a3,2i6,2es14.7)") 'End',kstep,kstmx,tt,tstop
      Write(50,"(4(a5,es14.7,1x))") (nname(i),aa(i)*y(i),i=1,ny)

! Write integrated fluxes to diagnotic output
!     iflx_mx=maxloc(abs(flx_int))
!     Write(50,'(2x,a8,i5,es11.3,8a5)') 'Flux Max',iflx_mx(1),
!    &    flx_int(iflx_mx(1)),nname(nflx(1:3,iflx_mx(1))),' <-> ',
!    &    nname(nflx(4:7,iflx_mx(1)))
!     Write(50,'(i5,8a5,i5,es11.3)') (k,nname(nflx(1:3,k)),' <-> ',
!    &    nname(nflx(4:7,k)),iwflx(k),flx_int(k),k=1,mflx)
! Write flux output formatted for Matlab
      Write(50,'(i5,4f6.1,es13.5,a5)') 
     &  (k,zz(ifl_orig(k)),nn(ifl_orig(k)),zz(ifl_term(k)),
     &  nn(ifl_term(k)),flx_int(k),descx(k),k=1,mflx)
        
! Test how well sums of fluxes match abbundances changes
!     dyf=0.0
!     Do k=1,mflx 
!        dyf(nflx(1:3,k))=dyf(nflx(1:3,k))+flx_int(k)
!        dyf(nflx(4:7,k))=dyf(nflx(4:7,k))-flx_int(k)
!     Enddo
!     dy=y-yo+dyf(1:ny)
!     Write(50,'(a5,es10.3)') 'DYF',dyf(0)
!     Write(50,'(a5,4es10.3)') 
!    &    (nname(k),dyf(k),dy(k),y(k),yo(k), k=1,ny)
      Return
      End

      subroutine sum_test(y)
      use nuclear_data
      real(8), dimension(ny) :: y
      real(8), dimension(ny,2) :: stest
      real(8), dimension(ny,2) :: xtot
      stest(:,1)=y
      stest(:,2)=10.0*y
!     xtot=sum((aa*stest),dim=2)
!     xtot(:,2)=aa*stest(:,2)
      Write(6,*) xtot
      Return
      End
