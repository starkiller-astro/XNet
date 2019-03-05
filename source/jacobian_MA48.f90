!*******************************************************************************
! Jacobian Sparce for MA48, part of XNet 7, 6/2/10
!
! The routines in this file are used to replace the standard dense Jacobian 
! and associated solver with the Harwell MA48 sparse solver package.  
!
! The bulk of the computational cost of the network (60-95%) is the solving 
! of the matrix equation.  Careful selection of the matrix solver is therefore 
! very important to fast computation.  For networks from a few dozen up to a 
! couple hundred species, hand tuned dense solvers such as those supplied by the 
! hardware manufacturer (often LAPACK) or third-parties like NAG, IMSL, etc. are 
! the fastest. However for larger matrices, sparse solvers are faster.  MA48 is 
! a solver from the Harwell Subroutine Library, which is available under an
! academic license.  We find it to be faster than PARDISO for matrices of size
! 100<ny<1000, but slower at larger sizes, but the availablilty of the MA48 source
! makes it valuable for some applications. 
!*******************************************************************************

Module jacobian_data
!===============================================================================
! Contains data for use in the sparse solver.
!===============================================================================
  Use nuclear_data
  Real(8), Dimension(:), Allocatable :: vals,tvals,sident,wB,wC
  Integer, Dimension(:), Allocatable :: ridx,cidx
  Integer, Dimension(:), Allocatable :: irn,jcn
  Integer, Dimension(:), Allocatable :: ns11,ns21,ns22
  Integer, Dimension(:), Allocatable :: ns31,ns32,ns33
  Integer, Dimension(:), Allocatable :: keep,iwA,iwB,iwC
  Integer :: lval,lia,nnz,msize
  Integer :: l1s,l2s,l3s
  Integer :: jobA,jobB,jobC
  Integer :: icntl(20),info(20)
  Real(8) :: cntl(10),rinfo(10),maxerr
  Logical :: trans = .false.
!$OMP THREADPRIVATE(tvals,vals,lia,jcn,irn,wB,wC,iwA,iwB,iwC, &
!$OMP   keep,info,rinfo,cntl,icntl,jobA,jobB,jobC)

  Real(8), Dimension(:), Allocatable :: mc29r,mc29c,mc29w
  Integer :: ifail
  Logical :: scale_jac
!$OMP THREADPRIVATE(mc29r,mc29c,mc29w,ifail,scale_jac)

End Module jacobian_data
  
Subroutine read_jacobian_data(data_dir)
!===============================================================================
! Reads in data necessary to use sparse solver and initializes the Jacobian data.
!===============================================================================
  Use controls
  Use reac_rate_data
  Use jacobian_data
       
  Character (LEN=*),  Intent(in)  :: data_dir
  Integer :: i,pb(ny+1)

  If(iheat>0) Then
    msize=ny+1
  Else
    msize=ny
  EndIf
  
!$OMP PARALLEL DEFAULT(SHARED)
! Set default values for MA48 control parameters
  Call MA48ID(cntl,icntl)

! Set error/diagnostic output to XNet diagnostic file (default=6)
  icntl(1) = lun_diag
  icntl(2) = lun_diag

! Set level of MA48 error/diagnostic output (default=2)
  If(idiag>=4) Then
    icntl(3) = 3
  EndIf

! Limit pivot search to maximum icntl(4) columns or use special search technique for 0 (default=3)
! icntl(4) = 0

! Set block size (default=32)
  If (msize<32) Then
    icntl(5) = msize
  ElseIf (msize<512) Then
    icntl(5) = 32
  ElseIf (msize<2048) Then
    icntl(5) = 128
  Else
    icntl(5) = 256
  EndIf

! Set the minimum size of a block of the block triangular form other than the final block (defualt=1)
! icntl(6) = msize

! Set option to move each column for which iwA(j)=0 to the end of the pivot sequence within its block (default=0)
! icntl(8) = 1

! Set option to automatically reset to job=1 if MA48BD fails for job=2
! icntl(11) = 1

! Set the pivoting threshold (near 0.0 emphasizes sparsity; near 1.0 emphasizes stability) (default=0.1)
! cntl(2) = 1.0
!$OMP END PARALLEL
  
  Open(600,file=trim(data_dir)//"/sparse_ind",status='old',form='unformatted')
  
  Read(600) lval
  If(iheat>0) Then
    nnz=lval+2*ny+1
  Else
    nnz=lval
  EndIf

  Allocate(ridx(nnz),cidx(nnz),sident(nnz))
!$OMP PARALLEL DEFAULT(SHARED)
  Allocate(tvals(nnz))
!$OMP END PARALLEL
  
!$OMP PARALLEL DEFAULT(SHARED)
  lia = 8*nnz
! These are variables to be used by the MA48 solver
  Allocate(vals(lia),jcn(lia),irn(lia))
  Allocate(wB(msize),wC(4*msize))
  Allocate(iwA(9*msize),iwB(4*msize),iwC(msize))
  If (icntl(8) == 0) Then
    Allocate(keep(6*msize+4*msize/icntl(6)+7-max(msize/icntl(6),1)))
  Else
    Allocate(keep(6*msize+4*msize/icntl(6)+7))
  EndIf

! These are variables to be used by MC29 if scaling is deemed necessary
  Allocate(mc29r(msize),mc29c(msize),mc29w(5*msize))
  mc29r = 0.0
  mc29c = 0.0
  mc29w = 0.0
  scale_jac = .false.
!$OMP END PARALLEL

! Set the value for the maximum allowed error in the call to MA48CD
  maxerr = 1.0d-11
  
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
  Use cross_sect_data
  Use jacobian_data
  Use timers
  
  Real(8), Intent(in) :: diag, mult
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
      l1=n31(j1)   ! n31(k) gives the index of first reactant in reaction mu3(k)
      l2=n32(j1)   ! n32(k) gives the index of second reactant in reaction mu3(k)
      l3=n33(j1)   ! n33(k) gives the index of third reactant in reaction mu3(k)
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
  Integer :: i,j,kdecomp,jdecomp
  
! Initiate timer
  start_timer = xnet_wtime()
  timer_solve = timer_solve - start_timer  
  
! Use Harwell MA48 sparse solver (this one is portable)
  If(kstep == 1 .or. jobB == 1) Then

    If(scale_jac) Then
      If(kstep == 1) Then
        ifail = 0
        Call MC29AD(msize,msize,nnz,tvals,ridx,cidx,mc29r,mc29c,mc29w,lun_diag,ifail)
        If(ifail<0) Write(lun_diag,"(1x,a,i2)") 'Error during MC29AD.  Code: ',ifail
        mc29r(1:msize) = exp(mc29r(1:msize))
        mc29c(1:msize) = exp(mc29c(1:msize))
      EndIf
    EndIf

    jobA = 3 ! Try to restrict pivoting to the diagonal

    Do kdecomp=1,5

      jcn(1:nnz) = cidx
      irn(1:nnz) = ridx
      vals(1:nnz) = tvals
      If(scale_jac) vals(1:nnz) = vals(1:nnz)*mc29r(irn(1:nnz))*mc29c(jcn(1:nnz))

! Perform symbolic analysis
      info = 0
      rinfo = 0.0
      Call MA48AD(msize,msize,nnz,jobA,lia,vals,irn,jcn,keep,cntl,icntl,iwA,info,rinfo)

      If(info(1) == 0 .and. info(4) <= lia .and. info(2) <= 10) Then
        jobB = 1 ! If analysis is successful, proceed to factorization
        Exit
      ElseIf(kdecomp == 5) Then
        Write(lun_diag,"(1x,a,i5,a,i2,a,i2)") 'kstep=',kstep,', Error during MA48AD. kdecomp=',kdecomp,', info(1)=',info(1)
        Stop
      Else
        Write(lun_diag,"(1x,a,i5,a,i2,a,i2)") 'kstep=',kstep,', Warning during MA48AD. kdecomp=',kdecomp,', info(1)=',info(1)

        If(info(1) == -3 .or. info(4) > lia) Then
          Write(lun_diag,"(1x,a,i7,a,i7)") 'Reallocating MA48 arrays: lia=',lia,' < info(4)=',max(info(3),info(4))
          Deallocate(irn,jcn,vals)
          lia = max(info(3),info(4))
          Allocate(irn(lia),jcn(lia),vals(lia))
        ElseIf(info(2) > 10) Then
          Write(lun_diag,"(1x,a,i3,a,i7)") 'Reallocating MA48 arrays: info(2)=',info(2),', info(4)=',max(info(3),info(4))
          Deallocate(irn,jcn,vals)
          lia = 2*max(info(3),info(4))
          Allocate(irn(lia),jcn(lia),vals(lia))
        EndIf  

        If(info(1) == 4 .and. jobA == 3) Then
          Write(lun_diag,"(1x,a,i3,a,i7)") 'Not possible to choose all pivots from diagonal'
          jobA = 1
        EndIf

      EndIf
    EndDo
  Else
! Build sparse value array    
    vals(1:nnz) = tvals
    If(scale_jac) vals(1:nnz) = vals(1:nnz)*mc29r(irn(1:nnz))*mc29c(jcn(1:nnz))
  EndIf

! Perform numerical decomposition using previously determined symbolic decomposition
  Do kdecomp=1,5
    info = 0
    rinfo = 0.0
    Call MA48BD(msize,msize,nnz,jobB,lia,vals,irn,jcn,keep,cntl,icntl,wB,iwB,info,rinfo)

    If(info(1) == 0 .and. info(4) <= lia .and. info(6) == 0) Then
      If( icntl(8)==0 ) Then
        jobB = 2 ! Unless using special case of icntl(8)/=0, use the "fast" MA48BD call
      Else
        jobB = 3 ! For the case of icntl(8)/=0, use the "intermediate" MA48BD all
      EndIf
      If(kstep/=0) jobC = 1
      Exit
    ElseIf(kdecomp == 5) Then
      Write(lun_diag,"(1x,a,i5,a,i2,a,i2)") 'kstep=',kstep,', Error during MA48BD. kdecomp=',kdecomp,', info(1)=',info(1)
      Stop
    Else
      Write(lun_diag,"(1x,a,i5,a,i2,a,i2)") 'kstep=',kstep,', Warning during MA48BD. kdecomp=',kdecomp,', info(1)=',info(1)
      
! Perform check to see if an error occured in MA48BD.  Most likely a singularity error
! caused by incorrect symbolic matrix, so run MA48AD again.  Also, make sure that 
! workspaces are large enough.
      If(info(1) == -3 .or. info(4) > lia) Then
        Write(lun_diag,"(1x,a,i7,a,i7)") 'Reallocating MA48 arrays: lia=',lia,' < info(4)=',info(4)
        Deallocate(irn,jcn,vals)
        lia = max(info(3),info(4))
        Allocate(irn(lia),jcn(lia),vals(lia))
      EndIf
      
      jobA = 3 ! Try to restrict pivoting to the diagonal

      Do jdecomp=1,5

! Rebuild vals array
        jcn(1:nnz) = cidx
        irn(1:nnz) = ridx
        vals(1:nnz) = tvals
        If(scale_jac) vals(1:nnz) = vals(1:nnz)*mc29r(irn(1:nnz))*mc29c(jcn(1:nnz))

! Perform symbolic analysis
        info = 0
        rinfo = 0.0
        Call MA48AD(msize,msize,nnz,jobA,lia,vals,irn,jcn,keep,cntl,icntl,iwA,info,rinfo)

        If(info(1) == 0 .and. info(4) <= lia .and. info(2) <= 10) Then
          jobB = 1 ! If analysis is successful, proceed to factorization
          Exit
        ElseIf(jdecomp == 5) Then
          Write(lun_diag,"(1x,a,i5,a,i2,a,i2)") 'kstep=',kstep,', Error during MA48AD. jdecomp=',jdecomp,', info(1)=',info(1)
          Stop
        Else
          Write(lun_diag,"(1x,a,i5,a,i2,a,i2)") 'kstep=',kstep,', Warning during MA48AD. jdecomp=',jdecomp,', info(1)=',info(1)

          If(info(1) == -3 .or. info(4) > lia) Then
            Write(lun_diag,"(1x,a,i7,a,i7)") 'Reallocating MA48 arrays: lia=',lia,' < info(4)=',max(info(3),info(4))
            Deallocate(irn,jcn,vals)
            lia = max(info(3),info(4))
            Allocate(irn(lia),jcn(lia),vals(lia))
          ElseIf(info(2) > 10) Then
            Write(lun_diag,"(1x,a,i3,a,i7)") 'Reallocating MA48 arrays: info(2)=',info(2),', info(4)=',max(info(3),info(4))
            Deallocate(irn,jcn,vals)
            lia = 2*max(info(3),info(4))
            Allocate(irn(lia),jcn(lia),vals(lia))
          EndIf  

          If(info(1) >= 4 .and. jobA == 3) Then
            Write(lun_diag,"(1x,a,i3,a,i7)") 'Not possible to choose all pivots from diagonal'
            jobA = 1
          EndIf

        EndIf
      EndDo
    EndIf
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
  Use abundances
  Real(8), Intent(in)  :: yrhs(ny)
  Real(8), Intent(out) :: dy(ny)
  Real(8), Intent(in)  :: t9rhs
  Real(8), Intent(out) :: dt9
  Real(8) :: rhs(msize),dx(msize)
  Real(8) :: relerr(3)
  Integer :: kbksub
  
! Initiate timer
  start_timer = xnet_wtime()
  timer_solve = timer_solve - start_timer  

! Perform back substitution 
  If(kmon(2) > kitmx) Then
    jobC = 2 ! Previous NR iteration failed, so estimate error for possible recalculation of data structures
  Else
    jobC = 1 ! Do not estimate error
  EndIf

  Do kbksub=1,5

    rhs(1:ny)=yrhs
    If(iheat>0) rhs(ny+1)=t9rhs

    If(scale_jac) rhs(1:msize) = rhs(1:msize)*mc29r(1:msize)

    relerr = 0.0
    info = 0
    Call MA48CD(msize,msize,trans,jobC,lia,vals,irn,keep,cntl,icntl,rhs,dx,relerr,wC,iwC,info)

    If(info(1) == 0 .and. (jobC == 1 .or. maxval(relerr) <= maxerr)) Then
      If(scale_jac) dx(1:msize) = dx(1:msize)*mc29c(1:msize)
      dy=dx(1:ny)
      If(iheat>0) dt9=dx(ny+1)
      Exit
    ElseIf(kbksub == 5) Then
      Write(lun_diag,"(1x,a,i5,a,i2,a,i2)") 'kstep=',kstep,', Error during MA48CD. kbksub=',kbksub,', info(1)=',info(1)
      Stop
    Else
      Write(lun_diag,"(1x,a,i5,a,i2,a,i2)") 'kstep=',kstep,', Warning during MA48CD. kbksub=',kbksub,', info(1)=',info(1)

    ! If the relative error becomes sufficiently large, regenerate data structures
      If( jobC>1 .and. maxval(relerr)>maxerr ) Then
        Write(lun_diag,"(1x,a,3es12.5,a)") 'Warning: relerr=',(relerr(i),i=1,3),' > maxerr'
        jobB = 1
        Call jacobian_decomp(0)
      EndIf

    EndIf
  EndDo
    
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
