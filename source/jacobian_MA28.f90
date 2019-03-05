!*******************************************************************************
! Jacobian Sparce for MA28, part of XNet 7, 6/2/10
!
! The routines in this file are used to replace the standard dense Jacobian 
! and associated solver with the Harwell MA28 sparse solver package.  
!
! The bulk of the computational cost of the network (60-95%) is the solving 
! of the matrix equation.  Careful selection of the matrix solver is therefore 
! very important to fast computation.  For networks from a few dozen up to a 
! couple hundred species, hand tuned dense solvers such as those supplied by the 
! hardware manufacturer (often LAPACK) or third-parties like NAG, IMSL, etc. are 
! the fastest. However for larger matrices, sparse solvers are faster.  MA28 is 
! an obsolescent solver from the Harwell Subroutine Library, which is therefore 
! publically available.  We find it to be less bulletproof and slower at large 
! network sizes than PARDISO, but the availablilty of the MA28 source makes it 
! valuable for some applications. 
!*******************************************************************************

Module jacobian_data
!===============================================================================
! Contains data for use in the sparse solver.
!===============================================================================
  Use nuclear_data
  Real(8), Dimension(:), Allocatable :: vals,ma28w,tvals,sident
  Integer, Dimension(:), Allocatable :: ridx,cidx
  Integer, Dimension(:), Allocatable :: irn,icn
  Integer, Dimension(:), Allocatable :: ns11,ns21,ns22
  Integer, Dimension(:), Allocatable :: ns31,ns32,ns33
  Integer, Dimension(:), Allocatable :: ikeep,iwA,iwB
  Integer :: lval,nnz,msize,lirn,licn,l1s,l2s,l3s
  Real(8)  :: pivoting = 0.1
End Module jacobian_data
  
Subroutine read_jacobian_data(data_dir)
!===============================================================================
! Reads in data necessary to use sparse solver and initializes the Jacobian data.
!===============================================================================
  Use controls
  Use reac_rate_data
  Use jacobian_data
       
  Character (LEN=*),  Intent(in)  :: data_dir
  Integer :: i, pb(ny+1)
  
  Common /MA28ED/ LP,MP,LBLOCK,GROW
  Common /MA28FD/EPS,RMIN,RESID,IRNCP,ICNCP,MINIRN,MINICN,IRANK,ABORT1,ABORT2
  Real(8) EPS,MRESID,RESID,RMIN
  Integer LP,MP,ICNCP,IRANK,IRNCP,MINICN,MINIRN
  Logical LBLOCK,GROW,ABORT1,ABORT2
  
! Do not abort MA28AD on numeriCally singular matrices
  ABORT1 = .false.
  ABORT2 = .false. 
  
  LP = lun_diag
  MP = lun_diag
  
  Open(600,file=trim(data_dir)//"/sparse_ind",status='old',form='unformatted')
  
  Read(600) lval

  If(iheat>0) Then
    nnz=lval+2*ny+1
    msize=ny+1
  Else
    nnz=lval
    msize=ny
  EndIf
  
  licn = 8*nnz
  lirn = 7*nnz

  Allocate(ridx(nnz),cidx(nnz))
  Allocate(tvals(nnz),sident(nnz))
  
! These are variables to be used by the MA28 solver
  Allocate(vals(licn),icn(licn),irn(lirn))
  Allocate(ma28w(msize))
  Allocate(iwA(8*msize),iwB(5*msize),ikeep(5*msize))
  
  Read(600) ridx(1:lval),cidx(1:lval),pb
  If(iheat>0) Then
! Add jacobian indices for temperature coupling
    Do i=1,ny
      cidx(i+lval) = ny+1    ! Extra column
      ridx(i+lval) = i
     
      cidx(i+lval+ny) = i    ! Extra row
      ridx(i+lval+ny) = ny+1
    EndDo
    cidx(lval+2*ny+1) = ny+1 ! dT9dot/dT9 term
    ridx(lval+2*ny+1) = ny+1
  EndIf
  
  ! Build a compressed row format version of the identity matrix
  sident=0.0
  Do i=1,nnz
    If (ridx(i)==cidx(i)) sident(i)=1.0
  EndDo
  
  Read(600) l1s,l2s,l3s
  
! Build  arrays for direct sparse representation Jacobian build
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
  Close(600)  
  
End Subroutine read_jacobian_data

Subroutine jacobian_build(diag,mult)
!===============================================================================
! This routine calculates the Jacobian matrix dYdot/dY, and and augments 
! by multiplying all elements by mult and adding diag to the diagonal elements.
!===============================================================================
  Use controls
  Use conditions
  Use abundances
  Use reac_rate_data
  Use jacobian_data
  Use timers
  
  Real(8), Intent(in) :: diag, mult
  Integer :: err
  Integer :: i,j,kstep,i0,la1,le1,la2,le2,la3,le3,j1,l1,l2,l3
  Integer :: ls1,ls2,ls3
  Real(8) :: dydotdt9(ny),dt9dotdy(ny)
  
! Initiate timer
  start_timer = xnet_wtime()
  timer_jacob = timer_jacob - start_timer  

! The quick solution to taking advantage of sparseness is to create a values array that has the maximum
! number of non-zero elements as well as a 2-by-#non-zero elements matrix and the other vectors
! required by the sparse solver.  The second matrix will contain the ordered pairs that map to a 
! particular place in the Jacobian.  Build the Jacobian as it is built now, Then pull the values 
! from it using the ordered pairs to fill the values array. 

! Build the Jacobian, species by species
  tvals = 0.0

  Do i0=1,ny
    la1=la(1,i0)
    le1=le(1,i0)
    Do j1=la1,le1
      ls1 = ns11(j1) ! ns11(j1) gives the index effected reaction j1 by in the compressed row storage scheme
      tvals(ls1)=tvals(ls1)+b1(j1)
    EndDo
    la2=la(2,i0) 
    le2=le(2,i0)  
    Do j1=la2,le2
      ls1=ns21(j1) ! ns21(j1) gives the first index effected reaction j1 by in the compressed row storage scheme
      ls2=ns22(j1) ! ns22(j1) gives the second index effected reaction j1 by in the compressed row storage scheme
      l1=n21(j1)   ! n21(k) gives the index of first reactant in reaction mu2(k)
      l2=n22(j1)   ! n22(k) gives the index of second reactant in reaction mu2(k)
      tvals(ls1)=tvals(ls1)+b2(j1)*yt(l2)
      tvals(ls2)=tvals(ls2)+b2(j1)*yt(l1)
    EndDo
    la3=la(3,i0)
    le3=le(3,i0)
    Do j1=la3,le3
      ls1=ns31(j1) ! ns31(j1) gives the first index effected reaction j1 by in the compressed row storage scheme
      ls2=ns32(j1) ! ns32(j1) gives the second index effected reaction j1 by in the compressed row storage scheme
      ls3=ns33(j1) ! ns33(j1) gives the third index effected reaction j1 by in the compressed row storage scheme
      l1=n31(j1)   ! n21(k) gives the index of first reactant in reaction mu2(k)
      l2=n32(j1)   ! n22(k) gives the index of second reactant in reaction mu2(k)
      l3=n33(j1)   ! n22(k) gives the index of third reactant in reaction mu2(k)
      tvals(ls1)=tvals(ls1)+b3(j1)*yt(l2)*yt(l3)
      tvals(ls2)=tvals(ls2)+b3(j1)*yt(l1)*yt(l3)
      tvals(ls3)=tvals(ls3)+b3(j1)*yt(l1)*yt(l2)
    EndDo
  EndDo

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
    tvals(lval+1:lval+ny)=dydotdt9

    dt9dotdy=0.0
    Do j1=1,lval
      dt9dotdy(cidx(j1))=dt9dotdy(cidx(j1))+mex(ridx(j1))*tvals(j1)/cv
    EndDo
    tvals(lval+ny+1:lval+2*ny)=-dt9dotdy
    tvals(nnz)=-sum(mex*dydotdt9)/cv

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
  Integer :: i,info
  
  Call jacobian_decomp(kstep)
  Call jacobian_bksub(yrhs,dy,t9rhs,dt9)
  
! Diagnostic output
  If(idiag>=4) Then
    Write(lun_diag,"(a)") 'JAC_SOLV'
    Write(lun_diag,"(14es10.3)") dy
    If(iheat>0) Write(lun_diag,"(es10.3)") dt9
  EndIf

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
  Integer :: i,j,err,kdecomp 
  
  Common /MA28FD/EPS,RMIN,RESID,IRNCP,ICNCP,MINIRN,MINICN,IRANK,ABORT1,ABORT2
  Common /MA28ED/ LP,MP,LBLOCK,GROW
  Common /MA28HD/TOL,THEMAX,BIG,DXMAX,ERRMAX,DRES,CGCE,NDROP,MAXIT,NOITER,NSRCH,ISTART,LBIG
  Real(8) EPS,MRESID,RESID,RMIN
  Integer ICNCP,IRANK,IRNCP,MINICN,MINIRN,LP,MP
  Logical ABORT1,ABORT2,LBLOCK,GROW
  Real(8) TOL,THEMAX,BIG,DXMAX,ERRMAX,DRES,CGCE 
  Integer NDROP,MAXIT,NOITER,NSRCH,ISTART
  Logical LBIG 

! Initiate timer
  start_timer = xnet_wtime()
  timer_solve = timer_solve - start_timer  

! Reduce the pivot ratio  
  EPS = 1.0e-6
  
! Use Harwell MA28 sparse solver (this one is portable)
  If (kstep == 1) Then
    icn = 0
    irn = 0
    vals = 0.0
    iwA = 0
    ma28w = 0.0
    ikeep = 0
    iwB = 0
    err = 0
    Do i=1,nnz
      icn(i) = cidx(i)
      irn(i) = ridx(i)
    EndDo
    vals(1:nnz) = tvals

! Perform symbolic decomposition
    Call MA28AD(msize,nnz,vals,licn,irn,lirn,icn,pivoting,ikeep,iwA,ma28w,err)
     
!   fs = 0
  EndIf

! Build sparse value array    
  vals(1:nnz) = tvals


! Perform numerical decomposition using previously determined symbolic decomposition
  Call MA28BD(msize,nnz,vals,licn,ridx,cidx,icn,ikeep,iwB,ma28w,err)
  If (err > 0) Then
    Write(lun_diag,*)'Warning during MA28BD.  Code: ',err,' Name: ',nname(err)
  ElseIf (err /= 0) Then
    Write(lun_diag,*)'Error during MA28BD.  Code: ',err
  EndIf
      
! Perform check to see if an error occured in MA28BD.  Most likely a singularity error
! caused by incorrect symbolic matrix, so run MA28AD again.  Also, make sure that 
! workspaces are large enough.
  If (ICNCP > msize/10) err = -66
  If (IRNCP > msize/10) err = -66
    
  Do kdecomp = 1,10
    If (err == 0) exit
    If (kdecomp>1) Write(lun_diag,*)'2. Error during MA28AD.  Code: ',err
    If (err < -2) Then
      Write(lun_diag,*)'Trying to reAllocate MA28 arrays'
      If (MINICN > licn) Then
         Write(lun_diag,*)'MINICN = ',MINICN,' > licn = ',licn
         Deallocate(icn,vals)
         licn = int(1.2*MINICN)
         Allocate(icn(licn),vals(licn))
      EndIf
      If (MINIRN > lirn) Then
        Write(lun_diag,*)'MINIRN = ',MINIRN,' > lirn = ',lirn
        Deallocate(irn)
        lirn = int(1.2*MINIRN)
        Allocate(irn(lirn))
      EndIf   
      If (IRNCP > msize/10) Then
        Write(lun_diag,*)'IRNCP = ',IRNCP,' > msize/10 = ',msize/10
        Deallocate(irn)
        lirn = lirn + int(lirn*0.5)
        Allocate(irn(lirn))
      EndIf  
      If (ICNCP > msize/10) Then
        Write(lun_diag,*)'ICNCP = ',ICNCP,' > msize/10 = ',msize/10
        Deallocate(icn,vals)
        licn = licn + int(licn*0.5)
        Allocate(icn(licn),vals(licn))
      EndIf
    EndIf
      
    ! Rebuild vals array
    icn = 0
    irn = 0
    vals = 0.0
    iwA = 0
    ma28w = 0.0
    ikeep = 0
    iwB = 0
    err = 0
    Do i=1,nnz
      icn(i) = cidx(i)
      irn(i) = ridx(i)
    EndDo
    vals(1:nnz)=tvals
    ikeep = 0
    iwA = 0
    ma28w = 0.0
    Call MA28AD(msize,nnz,vals,licn,irn,lirn,icn,pivoting,ikeep,iwA,ma28w,err)
    If(idiag>=1) Write(lun_diag,"(a,i6)") 'Recomputing symbolic decomposition at step ',kstep
     
  EndDo
    
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
  Real(8) :: rhs(msize)
  
! Initiate timer
  start_timer = xnet_wtime()
  timer_solve = timer_solve - start_timer  

  rhs(1:ny) = yrhs
  If(iheat>0) rhs(ny+1)=t9rhs

! Perform back substitution 
  Call MA28CD(msize,vals,licn,icn,ikeep,rhs,ma28w,1)
  dy = rhs(1:ny)
  If(iheat>0) dt9=rhs(ny+1)
    
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
  
