!*******************************************************************************
! Jacobian Sparce for PARDISO, part of XNet 7, 6/2/10
!
! The routines in this file are Used to replace the standard dense Jacobian 
! and associated solver with a sparse Jacobian build and the PARDISO sparse 
! solver package. 
!
! The bulk of the computational cost of the network (60-95%) is the solving 
! of the matrix equation.  Careful selection of the matrix solver is therefore 
! very important to fast computation.  For networks from a few dozen up to a 
! couple hundred species, hand tuned dense solvers such as those supplied by the 
! hardware manufacturer (often LAPACK) or third-parties like NAG, IMSL, etc. are 
! the fastest. However for larger matrices, sparse solvers are faster.  In our 
! experience, PARDISO is prefered sparse solver of XNet, because it is 
! generally the fastest and most bulletproof sparse package available, 
! with fairly wide distribution as part of the Intel Math Kernal Libraries. 
!*******************************************************************************

Module jacobian_data
!===============================================================================
! Contains data for use in the sparse solver.
!===============================================================================
  Use nuclear_data
  Real(8),dimension(:), Allocatable :: tvals,sident
  Integer,dimension(:), Allocatable :: cols,pb,perm
  Integer,dimension(:), Allocatable :: ridx,cidx
  Integer,dimension(:), Allocatable :: ns11,ns21,ns22
  Integer,dimension(:), Allocatable :: ns31,ns32,ns33
  Integer :: PARDISO_mem_stat,lval,pt(64),iparm(64)
  Integer :: l1s,l2s,l3s
  Integer :: mtype,maxfct,mnum,phase,nrhs,msglvl
End Module jacobian_data
  
Subroutine read_jacobian_data(data_dir)
!===============================================================================
! Reads in data necessary to use sparse solver and initializes the Jacobian data.
!===============================================================================
  Use controls
  Use nuclear_data
  Use reac_rate_data
  Use jacobian_data
       
  character (LEN=*),  INTENT(in)  :: data_dir
  Integer :: i
  
  Open(600,file=trim(data_dir)//"/sparse_ind",status='old',form='unformatted')
  
  Read(600) lval
  
  Allocate(ridx(lval),cidx(lval))
  Allocate(tvals(lval),sident(lval))
  Allocate(pb(ny+1),perm(ny))
  
! Read compressed row format storage arrays 
  Read(600) ridx,cidx,pb
  
! PARDISO Attributes
  mtype = 11
  maxfct= 1
  mnum  = 1
  nrhs  = 1
  msglvl= 0
  
! Build a compressed row format version of the identity matrix
  Do i=1,lval
    if (ridx(i)==cidx(i))sident(i)=1.0
  Enddo
  
! Read the sizes of the ns arrays
  Read(600)l1s,l2s,l3s
  
! Allocate arrays for direct sparse representation Jacobian build
  Allocate(ns11(l1s))
  Allocate(ns21(l2s),ns22(l2s))
  Allocate(ns31(l3s),ns32(l3s),ns33(l3s))
  
  ns11 = 0
  ns21 = 0
  ns22 = 0
  ns31 = 0
  ns32 = 0
  ns33 = 0
  
  Read(600) ns11,ns21,ns22
  Read(600) ns31
  Read(600) ns32
  Read(600) ns33
  close(600)  

! Initialize PARDISO 
  PARDISO_mem_stat = 0

! For Basel PARDISO libraries, call pardisoinit to initialize

! call pardisoinit(pt,mtype,iparm)

! For MKL, pardisoinit is broken, instead set by hand.
  
  pt=0
  iparm=0
  iparm(1) = 0 ! Due default IPARM values
  iparm(3) = 1 ! numbers of processors

! Write(lun_diag,"(a)") " i  pt   iparm" 
! Write(lun_diag,"(i2,2i8)") (i,pt(i),iparm(i),i=1,64)
      
End Subroutine read_jacobian_data
  
Subroutine jacobian_build(diag,mult)
!===============================================================================
! This routine calculates the Jacobian matrix dYdot/dY, and and augments 
! by multiplying all elements by mult and adding diag to the diagonal elements.
!===============================================================================
  Use controls
  Use nuclear_data
  Use conditions
  Use abundances
  Use reac_rate_data
  Use jacobian_data
  
  Real(8), Intent(in) :: diag, mult
  Integer :: err,knrcut
  Integer :: i,j,kstep,i0,la1,le1,la2,le2,la3,le3,j1,l1,l2,l3
  Integer :: ls1,ls2,ls3
  
  knrcut = int(kitmx*0.5)
  
! The quick solution to taking advantage of sparseness is to create a values array that has the maximum
! number of non-zero elements as well as a 2-by-#non-zero elements matrix and the other vectors
! required by the sparse solver.  The second matrix will contain the ordered pairs that map to a 
! particular place in the Jacobian.  The Jacobian is built by looping over species with a nested loop over 
! reactions which effect that row.  The sparse tvals array is built directly by using sparse mapping arrays
! determined in Read_sparse_data.  This saves significant time, but more importantly saves vast amounts of space.

  tvals = 0.0
! Build the Jacobian, species by species
  
  Do i0=1,ny
    la1=la(1,i0)
    le1=le(1,i0)
    Do j1=la1,le1
      ls1 = ns11(j1) ! ns11(j1) gives the index affected reaction j1 by in the compressed row storage scheme
      tvals(ls1)=tvals(ls1)+b1(j1)
    Enddo
    la2=la(2,i0) 
    le2=le(2,i0)  
    Do j1=la2,le2
      ls1=ns21(j1) ! ns21(j1) gives the first index affected reaction j1 by in the compressed row storage scheme
      ls2=ns22(j1) ! ns22(j1) gives the second index affected reaction j1 by in the compressed row storage scheme
      l1=n21(j1)   ! n21(k) gives the index of first reactant in reaction mu2(k)
      l2=n22(j1)   ! n22(k) gives the index of second reactant in reaction mu2(k)
      tvals(ls1)=tvals(ls1)+b2(j1)*yt(l2)
      tvals(ls2)=tvals(ls2)+b2(j1)*yt(l1)
    Enddo
    la3=la(3,i0)
    le3=le(3,i0)
    Do j1=la3,le3
      ls1=ns31(j1) ! ns31(j1) gives the first index affected reaction j1 by in the compressed row storage scheme
      ls2=ns32(j1) ! ns32(j1) gives the second index affected reaction j1 by in the compressed row storage scheme
      ls3=ns33(j1) ! ns33(j1) gives the third index affected reaction j1 by in the compressed row storage scheme
      l1=n31(j1)   ! n31(j1) gives the index of first reactant in reaction mu2(k)
      l2=n32(j1)   ! n32(j1) gives the index of second reactant in reaction mu2(k)
      l3=n33(j1)   ! n33(j1) gives the index of third reactant in reaction mu2(k)
      tvals(ls1)=tvals(ls1)+b3(j1)*yt(l2)*yt(l3)
      tvals(ls2)=tvals(ls2)+b3(j1)*yt(l1)*yt(l3)
      tvals(ls3)=tvals(ls3)+b3(j1)*yt(l1)*yt(l2)
    Enddo
  Enddo                
  

! Augment matrix with externally provided factors  
  tvals = mult * tvals
  tvals = tvals + sident * diag 
  
  If(idiag>=5) Then
    Write(lun_diag,"(a9,2es14.7)") 'JAC_build',diag,mult
    Write(lun_diag,"(14es9.1)") tvals
  EndIf

  Return   
End Subroutine jacobian_build
  
Subroutine jacobian_solve(kstep,rhs,dy) 
!===============================================================================
! This routine solves the system of abundance equations composed of the jacobian
! matrix and rhs vector.
!===============================================================================
  Use controls
  Use jacobian_data
  Integer, Intent(in)  :: kstep
  Real(8), Intent(IN)  :: rhs(ny)
  Real(8), Intent(out) ::  dy(ny)
  Integer :: i,err

! Use PARDISO sparse solver from MKL

! Check if this is the first step, or if convergence is taking to long
  If(kstep == 1.or.kmon(1)>1) Then

! If PARDISO memory allocated, free memory
! Comment this block for non-MKL PARDISO
    If(PARDISO_mem_stat==1) Then
      iparm(1) = 0
      iparm(3) = 1
      phase=-1
      Call pardiso(pt,maxfct,mnum,mtype,phase,ny,tvals,pb,cidx,perm,nrhs,iparm,msglvl,rhs,dy,err)
      If (err/=0) Write(6,*)'PARDISO error ',err
    Endif

! Analyze matrix
    iparm(1) = 0
    iparm(3) = 1
    phase = 11
    Call pardiso(pt,maxfct,mnum,mtype,phase,ny,tvals,pb,cidx,perm,nrhs,iparm,msglvl,rhs,dy,err)
    PARDISO_mem_stat=1
    If (err/=0) Write(6,*)'PARDISO error ',err
!   Write(lun_daig,"(a)") " i  pt   iparm" 
!   Write(lun_daig,"(i2,2i8)") (i,pt(i),iparm(i),i=1,64)
  Endif
   
! Solve matrix     
  iparm(1) = 0 
  iparm(3) = 1
  phase = 23
  Call pardiso(pt,maxfct,mnum,mtype,phase,ny,tvals,pb,cidx,perm,nrhs,iparm,msglvl,rhs,dy,err)
  If (err/=0) Write(6,*)'PARDISO error ',err
    
! Diagnostic output
  If(idiag>=4) Then
    Write(lun_diag,"(a)") 'JAC_SOLV'
    Write(lun_diag,"(14es10.3)") dy
  EndIf

  Return
End Subroutine jacobian_solve                                                                       
  
Subroutine jacobian_decomp(kstep) 
!===============================================================================
! This routine performs a matrix decomposition for the jacobian
!===============================================================================
  Use controls
  Use jacobian_data
  Integer, Intent(in)  :: kstep
  Integer :: err 
! Use PARDISO sparse solver from MKL

! Check if this is the first step, or if convergence is taking to long
  If(kstep == 1.or.kmon(1)>1) Then

! If PARDISO memory allocated, free memory
! Comment this block for non-MKL PARDISO
    If(PARDISO_mem_stat==1) Then
      iparm(1) = 0
      iparm(3) = 1
      phase=-1
      Call pardiso(pt,maxfct,mnum,mtype,phase,ny,tvals,pb,cidx,perm,nrhs,iparm,msglvl,rhs,dy,err)
      If (err/=0) Write(6,*)'PARDISO error ',err
    Endif

! Analyze matrix
    iparm(1) = 0
    iparm(3) = 1
    phase = 11
    Call pardiso(pt,maxfct,mnum,mtype,phase,ny,tvals,pb,cidx,perm,nrhs,iparm,msglvl,rhs,dy,err)
    PARDISO_mem_stat=1
    If (err/=0) Write(6,*)'PARDISO error ',err
!   Write(lun_daig,"(a)") " i  pt   iparm" 
!   Write(lun_daig,"(i2,2i8)") (i,pt(i),iparm(i),i=1,64)
  Endif
   
! Factor the matrix     
  iparm(1) = 0 
  iparm(3) = 1
  phase = 22
  Call pardiso(pt,maxfct,mnum,mtype,phase,ny,tvals,pb,cidx,perm,nrhs,iparm,msglvl,rhs,dy,err)
  If (err/=0) Write(6,*)'PARDISO error ',err
    
  Return
End Subroutine jacobian_decomp                                                                       
  
Subroutine jacobian_bksub(rhs,dy) 
!===============================================================================
! This routine performs back-substitution for a previously factored matrix and 
! the vector rhs.   
!===============================================================================
  Use controls
  Use jacobian_data
  Real(8), Intent(in) :: rhs(ny)
  Real(8), Intent(out) ::  dy(ny)
  Integer :: err 

! Solve matrix, using existing factorization     
  iparm(1) = 0 
  iparm(3) = 1
  phase = 33
  Call pardiso(pt,maxfct,mnum,mtype,phase,ny,tvals,pb,cidx,perm,nrhs,iparm,msglvl,rhs,dy,err)
  If (err/=0) Write(6,*)'PARDISO error ',err
    
! Diagnostic output
  If(idiag>=4) Then
    Write(lun_diag,"(a)") 'BKSUB'
    Write(lun_diag,"(14es10.3)") dy
  EndIf

  Return
End Subroutine jacobian_bksub

