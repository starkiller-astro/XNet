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
  Integer, Dimension(:), Allocatable :: ridx,cidx,pb
  Integer, Dimension(:), Allocatable :: irn,icn
  Integer, Dimension(:), Allocatable :: ns11,ns21,ns22
  Integer, Dimension(:), Allocatable :: ns31,ns32,ns33
  Integer, Dimension(:), Allocatable :: ikeep,iwA,iwB
  Integer :: lval,lirn,licn
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
  Integer :: la1,la2,la3,le1,le2,le3,l1,l2,l3,l1s,l2s,l3s
  Integer :: i,i0,i1,j1
  
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
  
  licn = 8*lval
  lirn = 7*lval

  Allocate(ridx(lval),cidx(lval))
  Allocate(tvals(lval),sident(lval))
  Allocate(pb(ny+1))
  
! These are variables to be used by the MA28 solver
  Allocate(vals(licn),icn(licn),irn(lirn))
  Allocate(ma28w(ny))
  Allocate(iwA(8*ny),iwB(5*ny),ikeep(5*ny))
  
  Read(600) ridx,cidx,pb
  
  ! Build a compressed row format version of the identity matrix
  Do i=1,lval
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
  
  Real(8), Intent(in) :: diag, mult
  Integer :: err
  Integer :: i,j,kstep,i0,la1,le1,la2,le2,la3,le3,j1,l1,l2,l3
  Integer :: ls1,ls2,ls3

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
  Real(8) :: d 
  Integer :: i,info
  
  Call jacobian_decomp
  Call jacobian_bksub(rhs,dy)
  
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
    Do i=1,lval
      icn(i) = cidx(i)
      irn(i) = ridx(i)
    EndDo
    vals(1:lval) = tvals

! Perform symbolic decomposition
    Call MA28AD(ny,lval,vals,licn,irn,lirn,icn,pivoting,ikeep,iwA,ma28w,err)
     
!   fs = 0
  EndIf

! Build sparse value array    
  vals(1:lval) = tvals


! Perform numerical decomposition using previously determined symbolic decomposition
  Call MA28BD(ny,lval,vals,licn,ridx,cidx,icn,ikeep,iwB,ma28w,err)
  If (err > 0) Then
    Write(LP,*)'Error during MA28BD.  Code: ',err,' Name: ',nname(err)
  ElseIf (err /= 0) Then
    Write(LP,*)'Error during MA28BD.  Code: ',err
  EndIf
      
! Perform check to see if an error occured in MA28BD.  Most likely a singularity error
! caused by incorrect symbolic matrix, so run MA28AD again.  Also, make sure that 
! workspaces are large enough.
  If (ICNCP > ny/10) err = -66
  If (IRNCP > ny/10) err = -66
    
  Do kdecomp = 1,10
    If (err == 0) exit
    If (kdecomp>1) Write(LP,*)'2. Error during MA28AD.  Code: ',err
    If (err < -2) Then
      Write(LP,*)'Trying to reAllocate MA28 arrays'
      If (MINICN > licn) Then
         Write(LP,*)'MINICN = ',MINICN,' > licn = ',licn
         Deallocate(icn,vals)
         licn = int(1.2*MINICN)
         Allocate(icn(licn),vals(licn))
      EndIf
      If (MINIRN > lirn) Then
        Write(LP,*)'MINIRN = ',MINIRN,' > lirn = ',lirn
        Deallocate(irn)
        lirn = int(1.2*MINIRN)
        Allocate(irn(lirn))
      EndIf   
      If (IRNCP > ny/10) Then
        Write(LP,*)'IRNCP = ',IRNCP,' > ny/10 = ',ny/10
        Deallocate(irn)
        lirn = lirn + int(lirn*0.5)
        Allocate(irn(lirn))
      EndIf  
      If (ICNCP > ny/10) Then
        Write(LP,*)'ICNCP = ',ICNCP,' > ny/10 = ',ny/10
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
    Do i=1,lval
      icn(i) = cidx(i)
      irn(i) = ridx(i)
    EndDo
    vals(1:lval)=tvals
    ikeep = 0
    iwA = 0
    ma28w = 0.0
    Call MA28AD(ny,lval,vals,licn,irn,lirn,icn,pivoting,ikeep,iwA,ma28w,err)
     
  EndDo
  Return
    
End Subroutine jacobian_decomp                                                                       
  
Subroutine jacobian_bksub(rhs,dy) 
!===============================================================================
! This routine performs back-substitution for a previously factored matrix and 
! the vector rhs.   
!===============================================================================
  Use controls
  Use jacobian_data
  Real(8), Intent(IN) :: rhs(ny)
  Real(8), Intent(out) ::  dy(ny)

! Perform back substitution 
  Call MA28CD(ny,vals,licn,icn,ikeep,rhs,ma28w,1)
  dy = rhs
    
! Diagnostic output
  If(idiag>=4) Then
    Write(lun_diag,"(a)") 'BKSUB'
    Write(lun_diag,"(14es10.3)") dy
  EndIf

  Return
End Subroutine jacobian_bksub
  
