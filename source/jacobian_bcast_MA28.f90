!*******************************************************************************
! Jacobian Broadcast for PARDISO matrix solver, part of XNet 6, 1/31/12
!
! Needed for MPI execution with PARDISO matrix solver.
! This routine broadcasts the jacobian data between MPI tasks.
!
!*******************************************************************************

Subroutine jacobian_bcast(data_dir)
!===============================================================================
! This routine distributes sparse Jacobian data for PARDISO solver
!===============================================================================
  Use controls
  Use file_data
  Use jacobian_data
  Use mpi
  Character (LEN=80) , INTENT(IN):: data_dir  
  Integer :: sparse_int(4),i,ierr

! MA28 Data structures
  Common /MA28ED/ LP,MP,LBLOCK,GROW
  Common /MA28FD/EPS,RMIN,RESID,IRNCP,ICNCP,MINIRN,MINICN,IRANK,ABORT1,ABORT2
  Real(8) EPS,MRESID,RESID,RMIN
  Integer LP,MP,ICNCP,IRANK,IRNCP,MINICN,MINIRN
  Logical LBLOCK,GROW,ABORT1,ABORT2
  
! For PARDISO solver, parameters for the compressed row storage must be
! read and broadcast, along with allocations being performed 
! On PE0 ...  
  If(myid==0) Then
    Call read_jacobian_data(data_dir)

! Pack sparse array size passing arrays
    sparse_int(1)  = lval ; sparse_int(2)  = l1s
    sparse_int(3)  = l2s  ; sparse_int(4)  = l3s

  EndIf

! Broadcast Jacobian array sizes
  call mpi_bcast(sparse_int,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! Unpack Jacobian array sizes and allocate on PE /=0
  If(myid/=0) Then

! Unpack sparse array size passing arrays
    lval = sparse_int(1) ; l1s = sparse_int(2)
    l2s  = sparse_int(3) ; l3s = sparse_int(4)
    If(iheat>0) Then
      nnz=lval+2*ny+1
      msize=ny+1
    Else
      nnz=lval
      msize=ny
    EndIf

    licn = 8*nnz
    lirn = 7*nnz

! Allocate sparse data arrays   
    Allocate(ridx(nnz),cidx(nnz),sident(nnz))
    Allocate(ns11(l1s),ns21(l2s),ns22(l2s))
    Allocate(ns31(l3s),ns32(l3s),ns33(l3s))

  EndIf ! PE/=0      

! Broadcast Jacobian data arrays
  call mpi_bcast(ridx, nnz,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(cidx, nnz,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ns11, l1s,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ns21, l2s,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ns22, l2s,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ns31, l3s,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ns32, l3s,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ns33, l3s,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  
! Complete MA28 Initilization on PE /=0
  If(myid/=0) Then

! Allocate variables used by the MA28 solver
    Allocate(tvals(nnz),vals(licn),icn(licn),irn(lirn))
    Allocate(ma28w(msize))
    Allocate(iwA(8*msize),iwB(5*msize),ikeep(5*msize))

! Build a compressed row format version of the identity matrix
    sident=0.0
    Do i=1,nnz
      If (ridx(i)==cidx(i)) sident(i)=1.0
    Enddo

! Do not abort MA28AD on numerically singular matrices
    ABORT1 = .false.
    ABORT2 = .false. 
  
    LP = lun_diag
    MP = lun_diag
  
  Endif ! PE/=0      

End Subroutine jacobian_bcast

