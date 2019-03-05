Program nse_slice
  Use controls
  Character (LEN=80) :: descript(3),data_desc,data_dir
  Character (LEN=80) :: diag_file_base='nse_diag'
  Character (LEN=80) :: diag_file
  Real(8) :: rho,ye,t9fnl
  Real(8) :: t9w(27)=(/100.,75.,50.,40.,35.,30.,28.,26.,24.,22.,20., &
&   18.,16.,15.,14.,13.,12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2./)
  Integer :: lun_control,ierr

! For a given density and Ye, this program solves for NSE as a funtion of T9
! down to T9fnl.

  Write(6,'(a)') 'Rho?'
  Read(5,*) rho
  Write(6,'(a)') 'Ye?'
  Read(5,*) ye
  Write(6,'(a)') 'T9 stop?'
  Read(5,*) t9fnl
  Write(6,'(a,f7.4,a,es11.3,a,f7.3)') 'NSE for Ye=',ye,'and Rho=',rho,', stopping at T9=',t9fnl

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
  Read(lun_control,*) iconvc       ! determines which convergence condition is Used
  Read(lun_control,*) changemx     ! allowed abundance change Used to set the timestep.
  Read(lun_control,*) yacc         ! abundances > yacc Used for timestep calculation
  Read(lun_control,*) tolm         ! mass conservation convergence criterion
  Read(lun_control,*) tolc         ! convergence limit on the iterative abundance change
  Read(lun_control,*) ymin         ! abundance < ymin is set to 0.0
  Read(lun_control,*) tdel_maxmult ! max factor by which the timestep is changed

! Read Output Controls
  call find_controls_block(lun_control,'Output Controls',ierr)
  Read(lun_control,*) idiag        ! sets diagnostic output level
  Read(lun_control,*) itsout       ! sets per timestep output level

! Read Input Controls
  call find_controls_block(lun_control,'Input Controls',ierr)
  Read(lun_control,"(72x)")
  Read(lun_control,"(a80)") data_dir
  Read(lun_control,"(72x)")

  Close(lun_control)

! Open diagnositic output file, per thread if OMP
  If(idiag>=0) Then
    lun_diag=50
    diag_file=trim(diag_file_base)
    Open(lun_diag,file=diag_file)
  EndIf

! In requested, pre-process the nuclear and reaction data.
  If(iprocess>0) Call net_preprocess(6,data_dir,data_dir)

! Initialize EoS for screening
  If(iscrn>0) call eos_initialize

! Read nuclear dataset
  Call read_nuclear_data(data_dir,data_desc)

! Allocate and initialize NSE arrays
  Call nse_initialize

! Calculate NSE abundances
  Call nse_descend(rho,ye,t9fnl,t9fnl)

  Close(lun_diag)

End Program nse_slice
