!*******************************************************************************
! Net Setup 6.0 6/13/2011
!
! This program processes the nuclear data directory to speed XNet execution.
! It re-orders and reformats the reaction data.
!
! The pre-processing can be performed at the begining of an XNet calculation, 
! using the routine net_preprocess, or via this standalone wrapper.
!*******************************************************************************

program net_setup
!===============================================================================
!  This program takes the human readable, ASCII data files sunet, netsu, 
!  netweak and netwinv and prepares binary versions for the reaction network. 
!  This involves translating characters to indicies and organizing the data  
!  for more efficient computation.  Several additional tests are performed.
!===============================================================================
  Character(LEN=80)  :: data_dir  ! The location of the data directory
  Character(LEN=80)  :: data_desc ! Description of data directory

! Process local directory
  data_dir='.'

! Solicit Data description
  Write(6,*) "Provide a one line description of this data set"
  Read(5,"(a80)") data_desc

  Call net_preprocess(6,data_dir,data_desc)

  Stop
End program net_setup

