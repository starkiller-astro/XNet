!***************************************************************************************************
! xnet_controls.f90 10/18/17
! This file contains modules and subroutines to control the execution of XNet.
!***************************************************************************************************

Module xnet_controls
  !-------------------------------------------------------------------------------------------------
  ! This module contains the values of the flags and limits which control the behavior of the
  ! network. Generally read from the control file.
  !-------------------------------------------------------------------------------------------------
  Use, Intrinsic :: iso_fortran_env, Only: lun_stderr=>error_unit, lun_stdin=>input_unit, lun_stdout=>output_unit
  Use xnet_types, Only: dp
  Implicit None

  ! Problem Description
  Character(80) :: descript(3)

  ! Job Controls
  Integer              :: szone      ! starting zone
  Integer              :: nzone      ! Number of zones
  Integer              :: zone_id(4) ! Current zone index quadruplet (ix,iy,iz,block)
  Integer, Allocatable :: iweak(:)   ! If >0, strong and weak reactions are used
                                     ! If =0, weak reactions are ignored
                                     ! If <0, only weak reactions are used
  Integer              :: iweak0     ! Saves input iweak flag
  Integer              :: iscrn      ! If =0, screening is ignored
  Integer              :: iprocess   ! If =0, assume network has been pre-processed.
                                     ! If >0, then run network pre-processing (slows calculation)
  !$omp threadprivate(zone_id,iweak)

  ! Zone Batching Controls
  Integer  :: nzbatchmx               ! Maximum number of zones in a batch
  Integer  :: nzbatch                 ! Active number of zones in a batch
  Integer  :: szbatch                 ! Starting zone for a batch
  Logical, Allocatable, Target :: lzactive(:) ! Mask for active zones
  !$omp threadprivate(nzbatch,szbatch,lzactive)

  ! Integration Controls
  Integer  :: isolv                 ! Sets the integration method (1=BE, 2=BD)
  Integer  :: kstmx                 ! Max # of timesteps before exit
  Integer  :: kitmx                 ! Max # of iterations or substeps within a timestep
  Integer  :: ijac                  ! Rebuild the Jacobian every ijac iterations after the first
  Integer  :: iconvc                ! Controls type of convergence condition (0=mass)
  Real(dp) :: changemx              ! Relative abundance change used to guess timestep
  Real(dp) :: yacc                  ! Min abundance required to be accuracte, used in timestep determination.
  Real(dp) :: tolc                  ! The iterative convergence test limit
  Real(dp) :: tolm                  ! Max network mass error
  Real(dp) :: ymin                  ! Minimum abundance, y < ymin = 0
  Real(dp) :: tdel_maxmult          ! new timestep <= tdel_maxmult * previous timestep
  Integer, Allocatable :: kmon(:,:) ! Solver-dependent behaviour monitors
  Integer, Allocatable :: ktot(:,:) ! Cumulative solver-dependent behaviour monitors
  Real(dp) :: t9min = 0.01          ! Temperature minimum for turning strong reactions off
  !$omp threadprivate(kmon,ktot)

  ! Self-heating Controls
  Integer  :: iheat        ! If >0, then implicitly couple network to temperature (self-heating)
  Real(dp) :: changemxt    ! Relative temperature change used to guess timestep
  Real(dp) :: tolt9        ! The iterative convergence test limit

  ! NSE Initial Conditions Controls
  Real(dp) :: t9nse = 8.0  ! Temperature maximum for switching to NSE

  ! Neutrinos Controls
  Integer  :: ineutrino    ! If >0, the include neutrino capture reactions

  ! Output controls
  Integer                   :: nnucout         ! Number of species to be output in condensed form
  Character(4)              :: nnucout_string  ! For output formatting
  Integer, Allocatable      :: inucout(:)      ! List of species to output in condensed form
  Character(5), Allocatable :: output_nuc(:)   ! Names of nuclei to be output in condensed form
  Integer       :: idiag                       ! Sets level of diagnostic output
  Integer       :: itsout                      ! Sets level of time series output
  Character(80) :: ev_file_base, bin_file_base ! Output filename bases
  Integer       :: lun_diag                    ! Logical units for per-thread diagnostic output file
  Integer, Allocatable :: lun_ts(:), lun_ev(:) ! Logical units for per-zone output files
  !$omp threadprivate(lun_diag,lun_ev,lun_ts)

  ! Input controls
  Character(80) :: inab_file_base, thermo_file_base ! Input filename bases
  Character(80), Allocatable :: inab_file(:)        ! Initial abundance files for each zone
  Character(80), Allocatable :: thermo_file(:)      ! Thermo trajectory files for each zone
  Integer       :: lun_th, lun_ab                   ! Logical units for input files
  !$omp threadprivate(lun_th,lun_ab)

  ! Job indentifiers
  Integer :: myid, nproc, mythread, nthread ! task & thread ids and counts
  !$omp threadprivate(mythread)

Contains

  Subroutine find_controls_block(lun_control,block_string,ifcb)
    !-----------------------------------------------------------------------------------------------
    ! This routine scans the control file for the input control block identifier string and
    ! positions the file to read that subset of the controls.
    !-----------------------------------------------------------------------------------------------
    Use, Intrinsic :: iso_fortran_env, Only: iostat_end
    Implicit None

    ! Input variables
    Integer, Intent(in)      :: lun_control
    Character(*), Intent(in) :: block_string

    ! Output variables
    Integer, Intent(out)     :: ifcb

    ! Local variables
    Integer, Parameter       :: line_limit = 1000
    Character(2)             :: test_marker, block_marker="##"
    Character(80)            :: test_read
    Integer                  :: ierr, itest
    Integer                  :: ii ! Loop Indicies

    ! Initialize return status
    ifcb = 1

    ! Rewind control file to the beginning
    Rewind(lun_control,iostat=ierr)

    ! Cycle through file
    Do ii = 1, line_limit
      Read(lun_control,"(a2,a80)",iostat=ierr) test_marker,test_read

      ! Test status of file
      If ( ierr == iostat_end ) Then ! EOF reached
        ifcb = 0
        Exit
      ElseIf ( ierr /= 0 ) Then ! Format Error
        Write(lun_stdout,*) 'Problem Finding Control Block ',trim(block_string)
      EndIf

      ! For lines begining with the block marker
      If ( test_marker == block_marker ) Then

        ! Search for the desired control block identifier string
        itest = index(test_read,block_string)

        ! If the string is found, exit loop, leaving file at desired location
        If ( itest > 0 ) Exit
      EndIf
    EndDo
    If ( ii >= line_limit .or. ifcb == 0 ) Then
      Write(lun_stdout,*) block_string,' not found!'
    EndIf

    Return
  End Subroutine find_controls_block

  Subroutine read_controls(data_dir)
    !-----------------------------------------------------------------------------------------------
    ! The control file contains the parameters which determine the actions of XNet.
    !-----------------------------------------------------------------------------------------------
    Use, Intrinsic :: iso_fortran_env, Only: iostat_end
    Use xnet_parallel, Only: parallel_bcast, parallel_IOProcessor
    Use xnet_util, Only: name_ordered, xnet_terminate
    Implicit None

    ! Output variables
    Character(80), Intent(out) :: data_dir ! Nuclear data directory

    ! Local variables
    Integer :: lun_control, i, ierr, izone, nzone_read

    ! Read Problem Description
    If ( parallel_IOProcessor() ) Then
      Open(newunit=lun_control, file='control', status='old')
      Call find_controls_block(lun_control,'Problem Description',ierr)
      Read(lun_control,"(a80)") (descript(i), i=1,3) ! text description of the problem.
    EndIf
    Call parallel_bcast(descript)

    ! Read Job Controls
    If ( parallel_IOProcessor() ) Then
      Call find_controls_block(lun_control,'Job Controls',ierr)
      Read(lun_control,*) szone        ! number of the zone with which to begin
      Read(lun_control,*) nzone        ! total # of zones
      Read(lun_control,*) iweak0       ! controls the treatment of weak reactions
      Read(lun_control,*) iscrn        ! controls the treatment of nuclear screening
      Read(lun_control,*) iprocess     ! controls the runtime pre-processing of the network data
    EndIf
    Call parallel_bcast(szone)
    Call parallel_bcast(nzone)
    Call parallel_bcast(iweak0)
    Call parallel_bcast(iscrn)

    ! Read Zone Batching Controls
    If ( parallel_IOProcessor() ) Then
      Call find_controls_block(lun_control,'Zone Batching Controls',ierr)
      If ( ierr /= 0 ) Then
        Read(lun_control,*) nzbatchmx  ! maximum number of zones in a batch
      Else
        Write(lun_stdout,*) 'Using Default Zone Batching behavior'
        nzbatchmx = 1
      EndIf
    EndIf
    Call parallel_bcast(nzbatchmx)
    !$omp parallel default(shared)
    Allocate (lzactive(nzbatchmx))
    Allocate (iweak(nzbatchmx),lun_ev(nzbatchmx),lun_ts(nzbatchmx))
    Allocate (kmon(2,nzbatchmx),ktot(5,nzbatchmx))
    !$omp end parallel

    ! Read Integration Controls
    If ( parallel_IOProcessor() ) Then
      Call find_controls_block(lun_control,'Integration Controls',ierr)
      Read(lun_control,*) isolv        ! Choice of integrations scheme
      Read(lun_control,*) kstmx        ! max # of timesteps for each zone
      Read(lun_control,*) kitmx        ! max # of iterations before retry
      Read(lun_control,*) ijac         ! rebuild jacobian every ijac iterations after the first
      Read(lun_control,*) iconvc       ! determines which convergence condition is Used
      Read(lun_control,*) changemx     ! allowed abundance change used to set the timestep.
      Read(lun_control,*) yacc         ! abundances > yacc Used for timestep calculation
      Read(lun_control,*) tolm         ! mass conservation convergence criterion
      Read(lun_control,*) tolc         ! convergence limit on the iterative abundance change
      Read(lun_control,*) ymin         ! abundance < ymin is set to 0.0
      Read(lun_control,*) tdel_maxmult ! max factor by which the timestep is changed
    EndIf
    Call parallel_bcast(isolv)
    Call parallel_bcast(kstmx)
    Call parallel_bcast(kitmx)
    Call parallel_bcast(ijac)
    Call parallel_bcast(iconvc)
    Call parallel_bcast(changemx)
    Call parallel_bcast(yacc)
    Call parallel_bcast(tolm)
    Call parallel_bcast(tolc)
    Call parallel_bcast(ymin)
    Call parallel_bcast(tdel_maxmult)

    ! Read Self-heating Controls (off by default)
    If ( parallel_IOProcessor() ) Then
      Call find_controls_block(lun_control,'Self-heating Controls',ierr)
      If ( ierr /= 0 ) Then
        Read(lun_control,*) iheat      ! controls the coupling of the network to temperature
        Read(lun_control,*) changemxt  ! allowed temperature change used to set the timestep.
        Read(lun_control,*) tolt9      ! convergence limit on the iterative temperature change
      Else
        Write(lun_stdout,*) 'Using Default Self-heating behavior'
        iheat = 0
        changemxt = 1.0e-2
        tolt9 = 1.0e-4
      EndIf
    EndIf
    Call parallel_bcast(iheat)
    Call parallel_bcast(changemxt)
    Call parallel_bcast(tolt9)

    ! If using higher-order solver, XNet does not need to further limit the timestep size
    If ( isolv == 3 ) Then
      changemx = 1.0e+10
      changemxt = 1.0e+10
    EndIf

    ! Read NSE Initial Conditions Controls (t9nse = 8.0 by default)
    If ( parallel_IOProcessor() ) Then
      Call find_controls_block(lun_control,'NSE Initial Conditions',ierr)
      If ( ierr /= 0 ) Then
        Read(lun_control,*) t9nse
      Else
        Write(lun_stdout,*) 'Using Default NSE behavior'
        t9nse = 8.0
      EndIf
    EndIf
    Call parallel_bcast(t9nse)

    ! Read Neutrino Controls
    If ( parallel_IOProcessor() ) Then
      Call find_controls_block(lun_control,'Neutrinos',ierr)
      If ( ierr /= 0 ) Then
        Read(lun_control,*) ineutrino
      Else
        Write(lun_stdout,*) 'Using Default Neutrinos behavior'
        ineutrino = 0
      EndIf
    EndIf
    Call parallel_bcast(ineutrino)

    !-----------------------------------------------------------------------------------------------
    ! XNet output controls include the base of the filenames to which ASCII and binary output are
    ! written at the end of each timestep, and a subset of nuclei to be included in the per timestep
    ! ASCII output file.
    !-----------------------------------------------------------------------------------------------
    ! Read Output Controls
    If ( parallel_IOProcessor() ) Then
      Call find_controls_block(lun_control,'Output Controls',ierr)
      Read(lun_control,*) idiag        ! sets diagnostic output level
      Read(lun_control,*) itsout       ! sets per timestep output level
      Read(lun_control,*)
      Read(lun_control,"(a80)") ev_file_base
      Read(lun_control,*)
      Read(lun_control,"(a80)") bin_file_base
      Read(lun_control,"(50x,i4)") nnucout
    EndIf
    Call parallel_bcast(idiag)
    Call parallel_bcast(itsout)
    Call parallel_bcast(ev_file_base)
    Call parallel_bcast(bin_file_base)
    Call parallel_bcast(nnucout)
    Allocate(output_nuc(nnucout),inucout(nnucout))
    Write(nnucout_string,"(i4)") nnucout
    nnucout_string = adjustl(nnucout_string)
    If ( parallel_IOProcessor() ) Then
      Read(lun_control,"("//nnucout_string//"a5)") output_nuc
      Write(lun_stdout,*) 'nnucout',nnucout
    EndIf
    Call parallel_bcast(output_nuc)

    !-----------------------------------------------------------------------------------------------
    ! XNet input controls include the relative directory from which the nuclear data should be
    ! loaded, as well as the names of the files containing the initial abundances and thermodynamic
    ! trajectories.
    !-----------------------------------------------------------------------------------------------
    ! Read Input Controls
    Allocate (inab_file(nzone),thermo_file(nzone))
    If ( parallel_IOProcessor() ) Then
      Call find_controls_block(lun_control,'Input Controls',ierr)
      Read(lun_control,*)
      Read(lun_control,"(a80)") data_dir
      Read(lun_control,*)
      nzone_read = nzone
      Do izone = 1, nzone
        Read(lun_control,"(a80)",iostat=ierr) inab_file(izone)
        If ( ierr == iostat_end ) Then
          Exit
        ElseIf ( ierr /= 0 ) Then
          Write(lun_stdout,*) 'Problem reading Input Filenames'
        EndIf
        Read(lun_control,"(a80)",iostat=ierr) thermo_file(izone)
        If ( ierr == iostat_end ) Then
          Exit
        ElseIf ( ierr /= 0 ) Then
          Write(lun_stdout,*) 'Problem reading Input Filenames'
        EndIf
      EndDo
      Close(lun_control)
      nzone_read = izone - 1
      If ( nzone_read == 1 .and. nzone_read < nzone ) Then
        inab_file_base = inab_file(1)
        thermo_file_base = thermo_file(1)
        Do izone = 1, nzone
          thermo_file(izone) = trim(thermo_file_base)
          Call name_ordered(thermo_file(izone),izone,nzone)
          inab_file(izone) = trim(inab_file_base)
          Call name_ordered(inab_file(izone),izone,nzone)
        EndDo
      ElseIf ( nzone_read /= nzone ) Then
        Write(lun_stdout,*) 'Number of datafiles does not match number of zones!'
!       Call xnet_terminate('Number of datafiles does not match number of zones!')
      Else
        inab_file_base = ''
        thermo_file_base = ''
      EndIf
    EndIf
    Call parallel_bcast(data_dir)
    Call parallel_bcast(inab_file)
    Call parallel_bcast(thermo_file)

    Return
  End Subroutine read_controls

End Module xnet_controls
