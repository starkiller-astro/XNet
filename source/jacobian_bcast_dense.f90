!*******************************************************************************
! Jacobian Broadcast for Dense matrix solvers, part of XNet 6, 1/31/12
!
! Needed for MPI execution with dense matrix solvers.
! This routine broadcasts the jacobian data between MPI tasks.
!
!*******************************************************************************

Subroutine jacobian_bcast(data_dir)
!===============================================================================
! This routine distributes Jacobian data for dense solver.
!===============================================================================
  Use controls
  Use file_data
  Use jacobian_data
  Use mpi
  Character (LEN=80), INTENT(IN) :: data_dir  

! For Dense Solvers, no additional data must be read or broadcast, only 
! allocations performed 
  Call read_jacobian_data(data_dir)
  
End Subroutine jacobian_bcast

