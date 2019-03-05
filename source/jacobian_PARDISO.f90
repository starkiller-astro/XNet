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
  Integer, parameter :: jacobian_type = 2
  Real(8),dimension(:), Allocatable :: tvals,sident
  Integer,dimension(:), Allocatable :: cols,perm
  Integer,dimension(:), Allocatable :: ridx,cidx,pb
  Integer,dimension(:), Allocatable :: ridxo,cidxo,pbo
  Integer,dimension(:), Allocatable :: ns11,ns21,ns22
  Integer,dimension(:), Allocatable :: ns31,ns32,ns33
  Integer :: PARDISO_mem_stat,lval,nnz,msize,iparm(64)
  Integer*8 :: pt(64)
  Integer :: l1s,l2s,l3s
  Integer :: mtype,maxfct,mnum,phase,nrhs,msglvl

! Threading Scope
!$OMP THREADPRIVATE(tvals,perm,PARDISO_mem_stat,iparm,pt,phase)
End Module jacobian_data
  
Subroutine read_jacobian_data(data_dir)
!===============================================================================
! Reads in data necessary to use sparse solver and initializes the Jacobian data.
!===============================================================================
  Use controls
  Use reac_rate_data
  Use jacobian_data
       
  character (LEN=*),  INTENT(in)  :: data_dir
  Integer :: i,i0,i1,la1,le1,la2,le2,la3,le3,j1,l1,l2,l3
  Integer :: rstart,rstarto,rend,rendo
  
  Open(600,file=trim(data_dir)//"/sparse_ind",status='old',form='unformatted')
  
  Read(600) lval
  If(iheat>0) Then
    nnz=lval+2*ny+1
    msize=ny+1
  Else
    nnz=lval
    msize=ny
  EndIf
  
  Allocate(ridx(nnz),cidx(nnz),sident(nnz),pb(msize+1))
!$OMP PARALLEL
  Allocate(tvals(nnz),perm(msize))
!$OMP End PARALLEL
  
! Read compressed row format storage arrays 
  If(iheat>0) Then
    Allocate(ridxo(lval),cidxo(lval),pbo(ny+1))
    Read(600) ridxo,cidxo,pbo
! Add jacobian indices for temperature coupling
    pb(1)=pbo(1)
    Do i0=1,ny
      ! Shift row pointers to adjust for additional column
      pb(i0+1)=pbo(i0+1)+i0
      rstarto=pbo(i0)
      rendo=pbo(i0+1)-1
      rstart=pb(i0)
      rend=pb(i0+1)-1

      ! Extra column
      cidx(rstart:rend-1)=cidxo(rstarto:rendo)
      ridx(rstart:rend-1)=ridxo(rstarto:rendo)
      cidx(rend)=ny+1
      ridx(rend)=i0
    EndDo
    Deallocate(ridxo,cidxo,pbo)

    ! Extra row
    pb(msize+1)=nnz+1
    rstart=pb(msize)
    rend=pb(msize+1)-1
    Do i0=1,ny
      cidx(rstart+i0-1)=i0
      ridx(rstart+i0-1)=ny+1
    EndDo

    ! dT9dot/dT9 term
    cidx(nnz) = ny+1
    ridx(nnz) = ny+1
  Else
    Read(600) ridx,cidx,pb
  EndIf
  
! PARDISO Attributes
  mtype = 11
  maxfct= 1
  mnum  = 1
  nrhs  = 1
  msglvl= 0
  
! Build a compressed row format version of the identity matrix
  sident=0.0
  Do i=1,nnz
    If(ridx(i)==cidx(i)) sident(i)=1.0
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
  
  If(iheat<=0) Then
    Read(600) ns11,ns21,ns22
    Read(600) ns31
    Read(600) ns32
    Read(600) ns33
  Else
! Rebuild arrays mapping reaction rates to locations in CRS Jacobian
    Do i0=1,ny
      la1=la(1,i0)
      le1=le(1,i0)
      Do j1=la1,le1
        l1=n11(j1)
        Do i1=1,nnz
          If(cidx(i1)==l1.and.ridx(i1)==i0) ns11(j1)=i1
        EndDo
      EndDo
      la2=la(2,i0)
      le2=le(2,i0)
      Do j1=la2,le2
        l1=n21(j1)
        l2=n22(j1)
        Do i1=1,nnz
          If(cidx(i1)==l1.and.ridx(i1)==i0) ns21(j1)=i1
          If(cidx(i1)==l2.and.ridx(i1)==i0) ns22(j1)=i1
        EndDo
      EndDo
      la3=la(3,i0)
      le3=le(3,i0)
      Do j1=la3,le3
        l1=n31(j1)
        l2=n32(j1)
        l3=n33(j1)
        Do i1=1,nnz
          If(cidx(i1)==l1.and.ridx(i1)==i0) ns31(j1)=i1
          If(cidx(i1)==l2.and.ridx(i1)==i0) ns32(j1)=i1
          If(cidx(i1)==l3.and.ridx(i1)==i0) ns33(j1)=i1
        EndDo
      EndDo
    EndDo
  EndIf
  Close(600)  

! Initialize PARDISO 
!$OMP PARALLEL DEFAULT(SHARED)
  PARDISO_mem_stat = 0

! For Basel PARDISO libraries, call pardisoinit to initialize

! call pardisoinit(pt,mtype,iparm)

! For MKL, pardisoinit is broken, instead set by hand.
  
  pt=0
  iparm=0
  iparm(1) = 0 ! Due default IPARM values
  iparm(3) = 1

! Write(lun_diag,"(a)") " i  pt   iparm" 
! Write(lun_diag,"(i2,2i16)") (i,pt(i),iparm(i),i=1,64)
!$OMP END PARALLEL
      
End Subroutine read_jacobian_data

Subroutine jacobian_build(diag,mult)
!===============================================================================
! This routine calculates the Jacobian matrix dYdot/dY, and and augments 
! by multiplying all elements by mult and adding diag to the diagonal elements.
!===============================================================================
  Use controls
  Use cross_sect_data
  Use nuclear_data
  Use conditions
  Use abundances
  Use reac_rate_data
  Use jacobian_data
  Use timers
  
  Real(8), Intent(in) :: diag, mult
  Integer :: err,knrcut
  Integer :: i,j,kstep,i0,la1,le1,la2,le2,la3,le3,j1,l1,l2,l3
  Integer :: ls1,ls2,ls3
  Real(8) :: dydotdt9(ny),dt9dotdy(msize)
  
! Initiate timer
  start_timer = xnet_wtime()
  timer_jacob = timer_jacob - start_timer  
  
  knrcut = int(kitmx*0.5)
  
! The quick solution to taking advantage of sparseness is to create a values array that has the maximum
! number of non-zero elements as well as a 2-by-#non-zero elements matrix and the other vectors
! required by the sparse solver.  The second matrix will contain the ordered pairs that map to a 
! particular place in the Jacobian.  The Jacobian is built by looping over species with a nested loop over 
! reactions which effect that row.  The sparse tvals array is built directly by using sparse mapping arrays
! determined in Read_sparse_data.  This saves significant time, but more importantly saves vast amounts of space.

! Build the Jacobian, species by species
  tvals = 0.0
  
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

  If(iheat>0) Then

    dr1dt9=a1*dcsect1dt9(mu1)*yt(n11)
    dr2dt9=a2*dcsect2dt9(mu2)*yt(n21)*yt(n22)
    dr3dt9=a3*dcsect3dt9(mu3)*yt(n31)*yt(n32)*yt(n33)

    dydotdt9=0.0
    Do i0=1,ny
      la1=la(1,i0)
      le1=le(1,i0)
      Do j1=la1,le1
        dydotdt9(i0)=dydotdt9(i0)+dr1dt9(j1)
      EndDo
      la2=la(2,i0)
      le2=le(2,i0)
      Do j1=la2,le2
        dydotdt9(i0)=dydotdt9(i0)+dr2dt9(j1)
      EndDo
      la3=la(3,i0)
      le3=le(3,i0)
      Do j1=la3,le3
        dydotdt9(i0)=dydotdt9(i0)+dr3dt9(j1)
      EndDo
    EndDo
    tvals(pb(2:ny+1)-1)=dydotdt9

    dt9dotdy=0.0
    Do i0=1,ny
      Do j1=pb(i0),pb(i0+1)-1
        dt9dotdy(cidx(j1))=dt9dotdy(cidx(j1))+mex(i0)*tvals(j1)/cv
      EndDo
    EndDo
    tvals(pb(ny+1):nnz)=-dt9dotdy(1:ny+1)

  EndIf

! Augment matrix with externally provided factors  
  tvals = mult * tvals
  tvals = tvals + sident * diag 
  
  If(idiag>=5) Then
    Write(lun_diag,"(a9,2es14.7)") 'JAC_build',diag,mult
    Write(lun_diag,"(14es9.1)") tvals
  EndIf

! Stop timer
  stop_timer = xnet_wtime()
  timer_jacob = timer_jacob + stop_timer

  Return   
End Subroutine jacobian_build
  
Subroutine jacobian_solve(kstep,yrhs,dy,t9rhs,dt9) 
!===============================================================================
! This routine solves the system of abundance equations composed of the jacobian
! matrix and rhs vector.
!===============================================================================
  Use controls
  Use jacobian_data
  Use timers
  Integer, Intent(in)  :: kstep
  Real(8), Intent(in)  :: yrhs(ny)
  Real(8), Intent(out) :: dy(ny)
  Real(8), Intent(in)  :: t9rhs
  Real(8), Intent(out) :: dt9
  Real(8) :: rhs(msize),dx(msize)
  Integer :: i,err
  
! Initiate timer
  start_timer = xnet_wtime()
  timer_solve = timer_solve - start_timer  

  rhs(1:ny) = yrhs
  If(iheat>0) rhs(ny+1) = t9rhs

! Use PARDISO sparse solver from MKL

! Check if this is the first step, or if convergence is taking to long
  If(kstep == 1.or.kmon(1)>1) Then

! If PARDISO memory allocated, free memory
! Comment this block for non-MKL PARDISO
    If(PARDISO_mem_stat==1) Then
      iparm(1) = 0
      iparm(3) = 1
      phase=-1
      Call pardiso(pt,maxfct,mnum,mtype,phase,msize,tvals,pb,cidx,perm,nrhs,iparm,msglvl,rhs,dx,err)
      If (err/=0) Write(6,*)'PARDISO error ',err
    Endif

! Analyze matrix
    iparm(1) = 0
    iparm(3) = 1
    phase = 11
    Call pardiso(pt,maxfct,mnum,mtype,phase,msize,tvals,pb,cidx,perm,nrhs,iparm,msglvl,rhs,dx,err)
    PARDISO_mem_stat=1
    If (err/=0) Write(6,*)'PARDISO error ',err
!   Write(lun_diag,"(a)") " i  pt   iparm" 
!   Write(lun_diag,"(i2,2i16)") (i,pt(i),iparm(i),i=1,64)
  Endif
   
! Solve matrix     
  iparm(1) = 0 
  iparm(3) = 1
  phase = 23
  Call pardiso(pt,maxfct,mnum,mtype,phase,msize,tvals,pb,cidx,perm,nrhs,iparm,msglvl,rhs,dx,err)
  If (err/=0) Write(6,*)'PARDISO error ',err
  dy = dx(1:ny)
  If(iheat>0) dt9 = dx(ny+1)
    
! Diagnostic output
  If(idiag>=4) Then
    Write(lun_diag,"(a)") 'JAC_SOLV'
    Write(lun_diag,"(14es10.3)") dy
    If(iheat>0) Write(lun_diag,"(es10.3)") dt9
  EndIf

! Stop timer
  stop_timer = xnet_wtime()
  timer_solve = timer_solve + stop_timer

  Return
End Subroutine jacobian_solve                                                                       
  
Subroutine jacobian_decomp(kstep) 
!===============================================================================
! This routine performs a matrix decomposition for the jacobian
!===============================================================================
  Use controls
  Use jacobian_data
  Use timers
  Integer, Intent(in)  :: kstep
  Integer :: err 
  real(8) ::ddum
  
! Initiate timer
  start_timer = xnet_wtime()
  timer_solve = timer_solve - start_timer  

! Use PARDISO sparse solver from MKL

! Check if this is the first step, or if convergence is taking to long
  If(kstep == 1.or.kmon(1)>1) Then

! If PARDISO memory allocated, free memory
! Comment this block for non-MKL PARDISO
    If(PARDISO_mem_stat==1) Then
      iparm(1) = 0
      iparm(3) = 1
      phase=-1
      Call pardiso(pt,maxfct,mnum,mtype,phase,msize,ddum,pb,cidx,perm,nrhs,iparm,msglvl,ddum,ddum,err)
      If (err/=0) Write(6,*)'PARDISO error ',err
    Endif

! Analyze matrix
    iparm(1) = 0
    iparm(3) = 1
    phase = 11
    Call pardiso(pt,maxfct,mnum,mtype,phase,msize,tvals,pb,cidx,perm,nrhs,iparm,msglvl,ddum,ddum,err)
    PARDISO_mem_stat=1
    If (err/=0) Write(6,*)'PARDISO error ',err
!   Write(lun_daig,"(a)") " i  pt   iparm" 
!   Write(lun_daig,"(i2,2i8)") (i,pt(i),iparm(i),i=1,64)
  Endif
   
! Factor the matrix     
  iparm(1) = 0 
  iparm(3) = 1
  phase = 22
  Call pardiso(pt,maxfct,mnum,mtype,phase,msize,tvals,pb,cidx,perm,nrhs,iparm,msglvl,ddum,ddum,err)
  If (err/=0) Write(6,*)'PARDISO error ',err
    
! Stop timer
  stop_timer = xnet_wtime()
  timer_solve = timer_solve + stop_timer

  Return
End Subroutine jacobian_decomp                                                                       
  
Subroutine jacobian_bksub(yrhs,dy,t9rhs,dt9) 
!===============================================================================
! This routine performs back-substitution for a previously factored matrix and 
! the vector rhs.   
!===============================================================================
  Use controls
  Use jacobian_data
  Use timers
  Real(8), Intent(in)  :: yrhs(ny)
  Real(8), Intent(out) :: dy(ny)
  Real(8), Intent(in)  :: t9rhs
  Real(8), Intent(out) :: dt9
  Real(8) :: rhs(msize),dx(msize)
  Integer :: err 
  
! Initiate timer
  start_timer = xnet_wtime()
  timer_solve = timer_solve - start_timer  

  rhs(1:ny) = yrhs
  If(iheat>0) rhs(ny+1) = t9rhs

! Solve matrix, using existing factorization     
  iparm(1) = 0 
  iparm(3) = 1
  phase = 33
  Call pardiso(pt,maxfct,mnum,mtype,phase,msize,tvals,pb,cidx,perm,nrhs,iparm,msglvl,rhs,dx,err)
  If (err/=0) Write(6,*)'PARDISO error ',err
  dy = dx(1:ny)
  If(iheat>0) dt9 = dx(ny+1)
    
! Diagnostic output
  If(idiag>=4) Then
    Write(lun_diag,"(a)") 'BKSUB'
    Write(lun_diag,"(14es10.3)") dy
    If(iheat>0) Write(lun_diag,"(es10.3)") dt9
  EndIf

! Stop timer
  stop_timer = xnet_wtime()
  timer_solve = timer_solve + stop_timer

  Return
End Subroutine jacobian_bksub
