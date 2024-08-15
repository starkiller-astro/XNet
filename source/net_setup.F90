!***************************************************************************************************
! net_setup.f90 10/18/17
! This program processes the nuclear data directory to speed XNet execution. It re-orders and
! reformats the reaction data.
!
! The pre-processing can be performed at the begining of an XNet calculation,
! using the routine net_preprocess or via this standalone wrapper.
!***************************************************************************************************

Program net_setup
  !-------------------------------------------------------------------------------------------------
  ! This program takes the human readable, ASCII data files sunet, netsu, netweak and netwinv and
  ! prepares binary versions for the reaction network. This involves translating characters to
  ! indicies and organizing the data for more efficient computation.
  !-------------------------------------------------------------------------------------------------
  Use xnet_controls, Only: lun_stdin, lun_stdout
  Use xnet_preprocess, Only: net_preprocess
  Implicit None

  ! Local variables
  Character(80) :: data_dir  ! The location of the data directory
  Character(80) :: data_desc ! Description of data directory

  ! Process local directory
  data_dir = '.'

  ! Solicit Data description
  Write(lun_stdout,*) "Provide a one line description of this data set"
  Read(lun_stdin,"(a80)") data_desc

  Call net_preprocess(lun_stdout,data_dir,data_desc)

End Program net_setup

