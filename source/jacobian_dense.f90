!*******************************************************************************
! Jacobian Dense, part of XNet 6, 5/25/10 
!
! The routines in this file assume a dense Jacobian and use a variety of linear 
! algebra packages.  Choice of matrix solver is determined by uncommented lines.
!
! The bulk of the computational cost of the network (60-95%) is the solving 
! of the matrix equation.  Careful selection of the matrix solver is therefore 
! very important to fast computation.  For networks from a few dozen up to a 
! couple hundred species, hand tuned dense solvers such as those supplied by the 
! hardware manufacturer (often LAPACK) or third-parties like NAG, IMSL, etc. are 
! the fastest. However for larger matrices, sparse solvers are faster.  
!*******************************************************************************
  
Module jacobian_data
!===============================================================================
! The jacobian matrix for the solver.      
!===============================================================================
  Use nuclear_data
  Integer, parameter :: jacobian_type = 1
  Real(8), Dimension(:,:), Allocatable :: jac !,jact ! the Jacobian Matrix
  Integer, Dimension(:), Allocatable :: indx
  Integer :: msize

! Threading Scope
!$OMP THREADPRIVATE(jac,indx)
End Module jacobian_data
  
Subroutine read_jacobian_data(data_dir)
!===============================================================================
! Initializes the Jacobian data. For the dense solver, the only jacobian data is 
! the matrix dimension, ny.
!===============================================================================
  Use controls
  Use jacobian_data
  Character (LEN=*),  Intent(in)  :: data_dir

  If(iheat>0) Then
    msize=ny+1
  Else
    msize=ny
  EndIf

!$OMP PARALLEL DEFAULT(SHARED)
  Allocate(jac(msize,msize),indx(msize))
!$OMP END PARALLEL
  
End Subroutine read_jacobian_data
  
Subroutine jacobian_build(diag,mult)
!===============================================================================
! This routine calculates the reaction Jacobian matrix, dYdot/dY, and augments 
! by multiplying all elements by mult and adding diag to the diagonal elements.
!===============================================================================
  Use controls
  Use jacobian_data
  Use abundances
  Use conditions
  Use reac_rate_data
  Use cross_sect_data
  Use timers
  Real(8), Intent(IN) :: diag, mult
  Integer :: i,j,kstep,i0,la1,le1,la2,le2,la3,le3,j1!,l1,l2,l3
  Real(8) :: dydotdt9(ny),dt9dotdy(msize)
!  Real(8), Dimension(ny) :: um
  
! Initiate timer
  start_timer = xnet_wtime()
  timer_jacob = timer_jacob - start_timer  

! Build the Jacobian
! Building the jacobian directly is faster with sufficient memory
! Use um to build row-by-row
  jac=0.0
  Do i0=1,ny 
!    um=0.0              

! Include reaction terms
    la1=la(1,i0) 
    le1=le(1,i0)  
    Do j1=la1,le1  
!      l1=n11(j1)
!      um(l1)=um(l1)+b1(j1)
      jac(i0,n11(j1))=jac(i0,n11(j1))+mult*b1(j1)
    EndDo 
    la2=la(2,i0) 
    le2=le(2,i0)  
    Do j1=la2,le2  
!      l1=n21(j1)
!      l2=n22(j1) 
!      um(l1)=um(l1)+b2(j1)*yt(l2) 
!      um(l2)=um(l2)+b2(j1)*yt(l1)  
      jac(i0,n21(j1))=jac(i0,n21(j1))+mult*b2(j1)*yt(n22(j1))
      jac(i0,n22(j1))=jac(i0,n22(j1))+mult*b2(j1)*yt(n21(j1))
    EndDo       
    la3=la(3,i0) 
    le3=le(3,i0)  
    Do j1=la3,le3 
!      l1=n31(j1) 
!      l2=n32(j1) 
!      l3=n33(j1)
!      um(l1)=um(l1)+b3(j1)*yt(l2)*yt(l3)    
!      um(l2)=um(l2)+b3(j1)*yt(l1)*yt(l3)     
!      um(l3)=um(l3)+b3(j1)*yt(l1)*yt(l2)      
      jac(i0,n31(j1))=jac(i0,n31(j1))+mult*b3(j1)*yt(n32(j1))*yt(n33(j1))
      jac(i0,n32(j1))=jac(i0,n32(j1))+mult*b3(j1)*yt(n31(j1))*yt(n33(j1))
      jac(i0,n33(j1))=jac(i0,n33(j1))+mult*b3(j1)*yt(n31(j1))*yt(n32(j1))
    EndDo                 
!   jac(i0,i0)=jac(i0,i0)+diag

! Augment matrix with externally provided factors  
!    um=um*mult
!    um(i0)=um(i0)+diag  

! Tranfer to matrix row
!    jac(i0,:)=um
!   jact(:,i0)=um ! or column if the solver wants the transpose
  
  EndDo                                                      
! jac=transpose(jact)

  If(iheat>0) Then

    dr1dt9=mult*a1*dcsect1dt9(mu1)*yt(n11)
    dr2dt9=mult*a2*dcsect2dt9(mu2)*yt(n21)*yt(n22)
    dr3dt9=mult*a3*dcsect3dt9(mu3)*yt(n31)*yt(n32)*yt(n33)

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
!     jac(ny+1,i0)=-sum(mex*jac(1:ny,i0))/cv
    EndDo
    jac(1:ny,ny+1)=dydotdt9
!   jac(ny+1,ny+1)=-sum(mex*dydotdt9)/cv

    Call dgemv('T',ny,msize,1.0/cv,jac,msize,mex,1,0.0,dt9dotdy,1)
    jac(ny+1,1:ny+1)=-dt9dotdy(1:ny+1)

  EndIf

  Do i0=1,msize
    jac(i0,i0)=jac(i0,i0)+diag
  EndDo

  If(idiag>=5) Then
    Write(lun_diag,"(a9,2es14.7)") 'JAC_build',diag,mult
    Do i=1,msize
      Write(lun_diag,"(14es9.1)") (jac(i,j),j=1,msize)
    EndDo
  EndIf
  
! Test the eigenvalues 
! If(idiag>=6) Then
!   Call eigen_test(kstep,jac,rdt)
! EndIf
  
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
  Real(8) :: rhs(msize)
  Real(8) :: d
  Integer :: info
  
! Initiate timer
  start_timer = xnet_wtime()
  timer_solve = timer_solve - start_timer  

  rhs(1:ny)=yrhs
  If(iheat>0) rhs(ny+1)=t9rhs

! To Use Num Rec LU Decomp, uncomment 3 lines below
! Call ludcmp(jac,msize,msize,indx,d)
! Call lubksb(jac,msize,msize,indx,rhs)
! dy=rhs(1:ny)
  
! To Use LAPACK solver,  uncomment one of the first two lines below and the third line 
! Call sgesv(msize,1,jac,msize,indx,rhs,msize,info) ! Single precision version
  Call dgesv(msize,1,jac,msize,indx,rhs,msize,info) ! Double precision version
  dy=rhs(1:ny)
  If(iheat>0) dt9=rhs(ny+1)
  
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
  Real(8) :: d 
  Integer :: i,j,info
  
! Initiate timer
  start_timer = xnet_wtime()
  timer_solve = timer_solve - start_timer  

! To Use Num Rec LU Decomp, uncomment line below
! Call ludcmp(jac,msize,msize,indx,d)
  
! To Use LAPACK ,  uncomment one of the first two lines below 
! Call sgetrf(msize,msize,jac,msize,indx,info)! Single precision version
  Call dgetrf(msize,msize,jac,msize,indx,info)
  
! Diagnostic output
  If(idiag>=4) Then
    Write(lun_diag,"(a3,i4)") 'LUD',info
    Do i=1,msize
      Write(lun_diag,"(14es9.1)") (jac(i,j),j=1,msize)
    EndDo
  EndIf
  
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
  Real(8), Intent(IN) :: yrhs(ny)
  Real(8), Intent(out) ::  dy(ny)
  Real(8), Intent(in) :: t9rhs
  Real(8), Intent(out) :: dt9
  Real(8) :: rhs(msize)
  Real(8) :: d 
  Integer :: i,info
  
! Initiate timer
  start_timer = xnet_wtime()
  timer_solve = timer_solve - start_timer  

  rhs(1:ny)=yrhs
  If(iheat>0) rhs(ny+1)=t9rhs

! To Use Num Rec LU Decomp, uncomment 2 lines below
! Call lubksb(jac,msize,msize,indx,rhs)
! dy=rhs(1:ny)
  
! To use LAPACK solver, uncomment one of the first two lines below and the third line 
! Call sgetrs('No transpose',msize,1,jac,msize,indx,rhs,msize,info) ! Single precision version
  Call dgetrs('No transpose',msize,1,jac,msize,indx,rhs,msize,info) ! Double precision version
  dy=rhs(1:ny)
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
  
