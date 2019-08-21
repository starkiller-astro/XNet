!***************************************************************************************************
! xnet_match.f90 10/18/17
! This file containes the data structures and input routines for matching reactions which involve
! the same nuclei (forward and reverse reactions as well as reactions with multiple components).
!***************************************************************************************************

Module xnet_match
  !-------------------------------------------------------------------------------------------------
  ! This module contains the data necessary to match up reactions.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None
  Integer                   :: mflx                               ! Number of unique reaction pairs
  Integer, Allocatable      :: nflx(:,:)                          ! The nuclei in each unique reaction
  Integer, Allocatable      :: iwflx(:)                           ! Reaction pair weak flag
  Integer, Allocatable      :: ifl1(:), ifl2(:), ifl3(:), ifl4(:) ! Maps reaction rate arrays to reaction pair arrays
  Real(dp), Allocatable     :: qflx(:)                            ! Reaction pair Q values
  Character(4), Allocatable :: descx(:)                           ! Descriptor for the reaction pair

Contains

  Subroutine read_match_data(data_dir)
    !-----------------------------------------------------------------------------------------------
    ! This routine reads in the reaction matching data and allocates the necessary arrays.
    !-----------------------------------------------------------------------------------------------
    Use reaction_data, Only: nreac
    Use xnet_controls, Only: idiag, lun_diag, lun_stdout
    Use xnet_parallel, Only: parallel_bcast, parallel_IOProcessor
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: data_dir

    ! Local variables
    Integer :: i, lun_match, nr(4)

    ! Open and read the matching data arrays
    If ( parallel_IOProcessor() ) Then
      Open(newunit=lun_match, file=trim(data_dir)//"/match_data", form="unformatted", status="old")
      Read(lun_match) mflx, nr

      ! Make sure match_data agrees with reaction_data
      Do i = 1, 4
        If ( nr(i) /= nreac(i) ) Write(lun_stdout,*) 'NR mismatch',i,nr(i),nreac(i)
      EndDo
    EndIf
    Call parallel_bcast(mflx)

    Allocate (ifl1(nreac(1)),ifl2(nreac(2)),ifl3(nreac(3)),ifl4(nreac(4)))
    Allocate (nflx(8,mflx),qflx(mflx),iwflx(mflx),descx(mflx))
    If ( parallel_IOProcessor() ) Then
      Read(lun_match) ifl1, ifl2, ifl3, ifl4
      Read(lun_match) nflx,qflx,iwflx,descx
      Close(lun_match)
    EndIf
    Call parallel_bcast(ifl1)
    Call parallel_bcast(ifl2)
    Call parallel_bcast(ifl3)
    Call parallel_bcast(ifl4)
    Call parallel_bcast(nflx)
    Call parallel_bcast(qflx)
    Call parallel_bcast(iwflx)
    Call parallel_bcast(descx)

    !$omp parallel default(shared)
    If ( idiag >= 0 ) Write(lun_diag,*) 'Match',mflx
    !$omp end parallel

    Return
  End Subroutine read_match_data

End Module xnet_match
