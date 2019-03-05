!-----------------------------------------------------------------------------
! This file containes the data structures and input routines for matching
! reactions which involve the same nuclei (forward and reverse reactions as 
! well as reactions with multiple components)
!-----------------------------------------------------------------------------
  
Module match_data
!-----------------------------------------------------------------------------
! This module contains the data necessary to match up reactions.
!-----------------------------------------------------------------------------
  Integer                               :: mflx
  Integer, Dimension(:,:), allocatable  :: nflx ! nflx(7,mflx)
  Integer, Dimension(:),   allocatable  :: iwflx,Ifl1,Ifl2,Ifl3
  Real(8),  Dimension(:),   allocatable  :: qflx
  Character(4), Dimension(:), allocatable :: descx
End Module match_data
  
Subroutine read_match_data(data_dir)
!-----------------------------------------------------------------------------
! This routine reads in the reaction matching data and Allocates the 
! necessary arrays.
!-----------------------------------------------------------------------------
  Use controls
  Use cross_sect_data
  Use match_data
  Integer  :: i,nr(3)
  Character (LEN=*) :: data_dir
  Open(10,FILE=trim(data_dir)//"/match_data",FORM="unformatted")
  Read(10) mflx,nr(1:3)
!$OMP PARALLEL DEFAULT(SHARED)
  Write(lun_diag,*) 'Match',mflx
!$OMP END PARALLEL
  Do i=1,3
    If(nr(i)/=nreac(i)) Write(6,*) 'NR mismatch',i,nr(i),nreac(i)
  Enddo
  Allocate (Ifl1(nr(1)),Ifl2(nr(2)),Ifl3(nr(3)))
  Read(10) Ifl1,Ifl2,Ifl3
  Allocate (nflx(7,mflx),qflx(mflx),iwflx(mflx),descx(mflx))
  Read(10) nflx,qflx,iwflx,descx
  Return
End Subroutine read_match_data
  
