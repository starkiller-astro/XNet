!*******************************************************************************
! Control.f90 6/14/11
! This file contains modules and subroutines to manage the parameters that 
! control the execution of XNet.  
!*******************************************************************************
  
Module controls
!===============================================================================
! This module contains the values of the flags and limits which control 
! the behavior of the network.  Generally read from the control file.
!===============================================================================

! Job Controls
  Integer :: szone        ! starting zone
  Integer :: nzone        ! Number of zones
  Integer :: zone_id(4)   ! Current zone index quadruplet (ix,iy,iz,block)
  Integer :: iweak        ! If >0, strong and weak reaction are Used
                          ! If =0, weak reaction are ignored
                          ! If <0, only weak reactions are Used
  Integer :: iweak0       ! Saves input iweak flag
  Integer :: iscrn        ! If =0, screening is ignored
  Integer :: iprocess     ! If =0, assume network has been pre-processed.
                          ! If >0, then run network pre-processing (slows calculation)
! Read Integration Controls
  Integer :: isolv        ! Sets the integration method (1=BE, 2=BD)
  Integer :: kstmx        ! Max # of timesteps before exit
  Integer :: kitmx        ! Max # of iterations or substeps within a timestep
  Integer :: ijac         ! Rebuild the Jacobian every ijac iterations after the first
  Integer :: iconvc       ! Controls type of convergence condition (0=mass)
  Real(8) :: changemx     ! Relative abundance change Used to guess timestep
  Real(8) :: yacc         ! Min abundance required to be accuracte, 
                          ! Used in timestep determination.
  Real(8) :: tolc         ! The iterative convergence test limit
  Real(8) :: tolm         ! Max network mass error
  Real(8) :: ymin         ! Minimum abundance, y < ymin =0
  Real(8) :: tdel_maxmult ! new timestep <= tdel_maxmult * previous timestep
  Integer :: kmon(2)      ! Solver-dependent behaviour monitors
  Integer :: ktot(2)      ! Cumulative solver-dependent behaviour monitors
  Real(8) :: t9min=0.01   ! Temperature minimum for turning strong reactions off
  Real(8) :: t9nse=8.0    ! Temperature maximum for switching to NSE

! Self-heating Controls
  Integer :: iheat        ! If >0, then implicitly couple network to temperature (self-heating)
  Real(8) :: changemxt    ! Relative temperature change used to guess timestep
  Real(8) :: tolt9        ! The iterative convergence test limit

! Output controls
  Integer :: idiag        ! Sets level of diagnostic output
  Integer :: itsout       ! Sets level of time series output
  Integer :: lun_diag, lun_ev, lun_ts, lun_th, lun_ab ! logical units
  Integer :: inout(14)    ! List of species to output in condensed form

! Job indentifiers
  Integer :: myid, nproc, mythread, nthread ! task & thread ids and counts

! Threading Scope
!$OMP THREADPRIVATE(mythread,lun_diag,lun_ev,lun_ts,lun_th,lun_ab,kmon,ktot,zone_id)
!$OMP THREADPRIVATE(iweak)
  
End Module controls

Module timers
!===============================================================================
! This module contains the performance counters for XNet.
!===============================================================================
Real(8)  :: timer_total  = 0.0  ! Total XNet execution time
Real(8)  :: timer_setup  = 0.0  ! Data loading or preprocessing time
Real(8)  :: timer_csect  = 0.0  ! Cross section calculation time 
Real(8)  :: timer_deriv  = 0.0  ! Derivative calculation time
Real(8)  :: timer_jacob  = 0.0  ! Jacobian building time
Real(8)  :: timer_solve  = 0.0  ! Solution time
Real(8)  :: timer_scrn   = 0.0  ! Screening and EOS time 
Real(8)  :: timer_output = 0.0  ! Output time 
Real(8)  :: start_timer  = 0.0  ! cpu time at the beginning of the timer block
Real(8)  :: stop_timer   = 0.0  ! cpu time at the end of the timer block

! Threading Scope
!$OMP THREADPRIVATE(timer_total,timer_setup,timer_csect,timer_deriv,timer_jacob,&
!$OMP   timer_solve,timer_scrn,timer_output,start_timer,stop_timer)

Contains

Function xnet_wtime()
!===============================================================================
! This function returns the wall time in a manner akin to OMP_GET_WTIME()
!===============================================================================
Integer(8) :: clock_read
Integer(8) :: clock_rate
Integer(8) :: clock_max
Real(8)    :: xnet_wtime

!$OMP CRITICAL(wtime)
Call SYSTEM_CLOCK( clock_read, clock_rate, clock_max )
xnet_wtime = DBLE( clock_read ) / DBLE( clock_rate )
!$OMP END CRITICAL(wtime)

Return
End Function xnet_wtime

End Module timers
  
Subroutine find_controls_block(lun_control,block_string,ifcb)
!===============================================================================
! This routine scans the control file for the input control block identifier 
! string and positions the file to read that subset of the controls.
!===============================================================================
  Integer, Intent(in)          :: lun_control
  Character(LEN=*), Intent(in) :: block_string
  Integer, Intent(out)         :: ifcb
  Character(LEN=2)             :: test_marker, block_marker="##"
  Character(LEN=80)            :: test_read
  Integer                      :: ierr,itest,line_limit=1000
  Integer                      :: ii ! Loop Indicies
  
! Rewind control file to the beginning
  Rewind (lun_control,IOSTAT=ierr)  
  ifcb=1

! Cycle through file
  Do ii=1,line_limit
    Read(lun_control,'(a2,a80)',IOSTAT=ierr) test_marker,test_read

! Test status of file.
    If(ierr < 0) Then ! EOF reached
      ifcb=0
      Exit
    ElseIf(ierr > 0) Then ! Format Error
      Write(6,*) 'Problem Finding Control Block ',trim(block_String)
    Endif

! For lines begining with the block marker
    If(test_marker==block_marker) Then
    
! Search for the desired control block identifier string    
      itest=index(test_read,block_string)
      
! If the string is found, exit loop, leaving file at desired location       
      If(itest>0) Exit
    Endif
  Enddo
  If(ii>=line_limit.or.ifcb==0) Then
    Write(6,*) block_string,' not found!'
  EndIf
End subroutine find_controls_block
