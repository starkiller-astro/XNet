!*******************************************************************************
! Jacobian Dense, part of XNet 7, 5/25/10 
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
  Real(8), Dimension(:,:), Allocatable :: jac !,jact ! the Jacobian Matrix
  Integer, Dimension(:), Allocatable :: indx
End Module jacobian_data
  
Subroutine read_jacobian_data(data_dir)
!===============================================================================
! Initializes the Jacobian data. For the dense solver, the only jacobian data is 
! the matrix dimension, ny.
!===============================================================================
  Use jacobian_data
  Character (LEN=*),  Intent(in)  :: data_dir

  Allocate (jac(ny,ny),indx(ny))
  
End Subroutine read_jacobian_data
  
Subroutine jacobian_build(diag,mult)
!===============================================================================
! This routine calculates the reaction Jacobian matrix, dYdot/dY, and augments 
! by multiplying all elements by mult and adding diag to the diagonal elements.
!===============================================================================
  Use controls
  Use jacobian_data
  Use abundances
  Use reac_rate_data
  Real(8), Intent(IN) :: diag, mult
  Integer :: i,j,kstep,i0,la1,le1,la2,le2,la3,le3,j1,l1,l2,l3
  Real(8), Dimension(ny) :: um
  
! Build the Jacobian, row by row
  Do i0=1,ny 
    um=0.0              

! Include reaction terms
    la1=la(1,i0) 
    le1=le(1,i0)  
    Do j1=la1,le1  
      l1=n11(j1)
      um(l1)=um(l1)+b1(j1)
    EndDo 
    la2=la(2,i0) 
    le2=le(2,i0)  
    Do j1=la2,le2  
      l1=n21(j1)
      l2=n22(j1) 
      um(l1)=um(l1)+b2(j1)*yt(l2) 
      um(l2)=um(l2)+b2(j1)*yt(l1)  
    EndDo       
    la3=la(3,i0) 
    le3=le(3,i0)  
    Do j1=la3,le3 
      l1=n31(j1) 
      l2=n32(j1) 
      l3=n33(j1)
      um(l1)=um(l1)+b3(j1)*yt(l2)*yt(l3)    
      um(l2)=um(l2)+b3(j1)*yt(l1)*yt(l3)     
      um(l3)=um(l3)+b3(j1)*yt(l1)*yt(l2)    
    EndDo                 

! Augment matrix with externally provided factors  
    um=um*mult
    um(i0)=um(i0)+diag  

! Tranfer to matrix row
    jac(i0,:)=um
!   jact(:,i0)=um ! or column if the solver wants the transpose
  
  EndDo                                                      
! jac=transpose(jact)
  
  If(idiag>=5) Then
    Write(lun_diag,"(a9,2es14.7)") 'JAC_build',diag,mult
    Do i=1,ny
      Write(lun_diag,"(14es9.1)") (jac(i,j),j=1,ny)
    EndDo
  EndIf
  
! Test the eigenvalues 
! If(idiag>=6) Then
!   Call eigen_test(kstep,jac,rdt)
! EndIf

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
  Real(8), Intent(in)  :: rhs(ny)
  Real(8), Intent(out) ::  dy(ny)
  Real(8) :: d 
  Integer :: info
  
! To Use Num Rec LU Decomp, uncomment 3 lines below
! Call ludcmp(jac,ny,ny,indx,d)
! Call lubksb(jac,ny,ny,indx,rhs)
! dy=rhs
  
! To Use LAPACK solver,  uncomment one of the first two lines below and the third line 
! Call sgesv(ny,1,jac,ny,indx,rhs,ny,info) ! Single precision version
  Call dgesv(ny,1,jac,ny,indx,rhs,ny,info) ! Double precision version
  dy=rhs
  
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
  Real(8) :: d 
  Integer :: i,j,info
  
! To Use Num Rec LU Decomp, uncomment line below
! Call ludcmp(jac,ny,ny,indx,d)
  
! To Use LAPACK ,  uncomment one of the first two lines below 
! Call sgetrf(ny,ny,jac,ny,indx,info)! Single precision version
  Call dgetrf(ny,ny,jac,ny,indx,info)
  
! Diagnostic output
  If(idiag>=4) Then
    Write(lun_diag,"(a3,i4)") 'LUD',info
    Do i=1,ny
      Write(lun_diag,"(14es9.1)") (jac(i,j),j=1,ny)
    EndDo
  EndIf
  
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
  Real(8) :: d 
  Integer :: i,info
  
! To Use Num Rec LU Decomp, uncomment 2 lines below
! Call lubksb(jac,ny,ny,indx,rhs)
! dy=rhs
  
! To use LAPACK solver, uncomment one of the first two lines below and the third line 
! Call sgetrs('No transpose',ny,1,jac,ny,indx,rhs,ny,info) ! Single precision version
  Call dgetrs('No transpose',ny,1,jac,ny,indx,rhs,ny,info) ! Double precision version
  dy=rhs
  
! Diagnostic output
  If(idiag>=4) Then
    Write(lun_diag,"(a)") 'BKSUB'
    Write(lun_diag,"(14es10.3)") dy
  EndIf

  Return                                                                    
End Subroutine jacobian_bksub                                                                       
  
