!*******************************************************************************
! Jacobian Broadcast for MA41 matrix solver, part of XNet 6, 1/31/12
!
! Needed for MPI execution with MA41 matrix solver.
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
  Character (LEN=80), INTENT(IN) :: data_dir  
  Integer :: sparse_int(4),i,ierr

! For MA41 solver, parameters for the compressed row storage must be
! read and broadcast, along with allocations being performed 
! On PE0 ...  
  If(myid==0) Then
    Call read_jacobian_data(data_dir)

! Pack sparse array size passing arrays
    sparse_int(1)  = lval
    sparse_int(2)  = l1s
    sparse_int(3)  = l2s
    sparse_int(4)  = l3s

  EndIf

! Broadcast Jacobian array sizes
  Call mpi_bcast(sparse_int,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! Unpack Jacobian array sizes and allocate on PE /=0
  If(myid/=0) Then

! Unpack sparse array size passing arrays
    lval = sparse_int(1)
    l1s = sparse_int(2)
    l2s = sparse_int(3)
    l3s = sparse_int(4)

    If(iheat>0) Then
      nnz=lval+2*ny+1
      msize=ny+1
    Else
      nnz=lval
      msize=ny
    EndIf

! Allocate sparse data arrays   
    Allocate(ridx(nnz),cidx(nnz),sident(nnz))
    Allocate(ns11(l1s),ns21(l2s),ns22(l2s))
    Allocate(ns31(l3s),ns32(l3s),ns33(l3s))

  EndIf ! PE/=0      

! Broadcast Jacobian data arrays
  call mpi_bcast(ridx,nnz,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(cidx,nnz,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ns11,l1s,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ns21,l2s,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ns22,l2s,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ns31,l3s,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ns32,l3s,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ns33,l3s,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! Complete MA41 Initilization on PE /=0
  If(myid/=0) Then

! Build a compressed row format version of the identity matrix
    sident=0.0
    Do i=1,nnz
      If (ridx(i)==cidx(i)) sident(i)=1.0
    EndDo
  
!$OMP PARALLEL DEFAULT(SHARED)
! Set default values for MA41 control parameters
    Call MA41ID(cntl,icntl,keep)

! Set error/diagnostic output to XNet diagnostic file (default=6)
    icntl(1) = lun_diag
    If(idiag>=1) icntl(2) = lun_diag

! Set statistics output (default=0)
    If(idiag>=2) icntl(3) = lun_diag

! Set level of MA41 error/diagnostic output (default=2)
    If(idiag>=4) icntl(4) = 3

! Set the number of processors (default=1)
!   icntl(5) = 2

! Use a maximum transversal algorithm to get a zero-free diagonal (default=0)
!   icntl(6) = 1

! Use custom pivot order (must set the variable iwork manually) (default=0)
!   icntl(7) = 1

! Scaling strategy (default=0 (no scaling))
!   icntl(8) = 6

! Solve A*x = b (default=1), or A^T*x = b
!   icntl(9) = 0

! Maximum number of steps of iterative refinement (default=0)
!   icntl(10) = 10

! Return the infinity norm, solution, scaled residual, backward error estimates, forward error estimates to rinfo(4:9) (default=0)
!   icntl(11) = 1

! Set the pivoting threshold (near 0.0 emphasizes sparsity; near 1.0 emphasizes stability) (default=0.01)
!   cntl(1) = 0.0

! This holds the jacobian matrix
    Allocate(tvals(nnz))
  
! These are variables to be used by the MA41 solver
    Allocate(jcn(nnz),irn(nnz))
    If(icntl(7)==1) Then
      maxiw = nnz + 12*msize + 1
    Else
      maxiw = 2*nnz + 12*msize + 1
    EndIf
    maxw = 2*maxiw
    Allocate(iwork(maxiw),work(maxiw))

! These are variables to be used by MC29 if scaling is deemed necessary
    If(icntl(8)/=0) Then
      Allocate(mc29r(msize),mc29c(msize))
    Else
      Allocate(mc29r(1),mc29c(1))
    EndIf
    mc29r = 0.0
    mc29c = 0.0
!$OMP END PARALLEL

! Set the value for the maximum allowed error in the call to MA41AD w/ job = 3,5,6 and icntl(11)>0
    maxerr = 1.0d-11

  Endif ! PE/=0      

End Subroutine jacobian_bcast
