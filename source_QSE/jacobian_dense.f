!*******************************************************************************
!  Jacobian Dense, part of Full_Net 4.10 4/20/07 
!
!  The routines in this file assume a dense Jacobian and use a variety of 
!  solvers.  Choice of solver is determined by uncommented lines.
!  The routines in  
!*******************************************************************************

      subroutine read_jacobian_data(data_dir)
!===============================================================================
!  Placeholder for sparse solvers
!  For the dense solver, the only jacobian data is the matrix dimension, ny
!===============================================================================
      end subroutine read_jacobian_data

      subroutine netmatr(kstep)
!===============================================================================
!  This routine calculates the Jacobian matrix dYdot/dY, and solves for the 
!  Newton-Raphson iteration, dy.
!===============================================================================
      use controls
      use nuclear_data
      use conditions
      use abundances
      use reac_rate_data
      integer, dimension(ny) :: indx
      integer :: i,j,kstep,i0,la1,le1,la2,le2,la3,le3,j1,l1,l2,l3
      real(8), dimension(ny) :: f,ydot0,um
      real(8), dimension(ny,ny) :: am  !,at ! the Jacobian Matrix
      real(8) :: alph,beta,rdt,d

      alph=0.0                                                                  
      beta=1.0                                                                  
      ydot0=0.0

!  Build the Jacobian, row by row
      rdt=1.0/tdel/beta 
      Do i0=1,ny 
        um=0.0              
        um(i0)=rdt  
        la1=la(1,i0) 
        le1=le(1,i0)  
        Do j1=la1,le1  
          l1=n11(j1)
          um(l1)=um(l1)-b1(j1)
        Enddo 
        la2=la(2,i0) 
        le2=le(2,i0)  
        Do j1=la2,le2  
          l1=n21(j1)
          l2=n22(j1) 
          um(l1)=um(l1)-b2(j1)*yt(l2) 
          um(l2)=um(l2)-b2(j1)*yt(l1)  
        Enddo       
        la3=la(3,i0) 
        le3=le(3,i0)  
        Do j1=la3,le3 
          l1=n31(j1) 
          l2=n32(j1) 
          l3=n33(j1)
          um(l1)=um(l1)-b3(j1)*yt(l2)*yt(l3)    
          um(l2)=um(l2)-b3(j1)*yt(l1)*yt(l3)     
          um(l3)=um(l3)-b3(j1)*yt(l1)*yt(l2)      
        Enddo                 

!  Tranfer to matrix row
        am(i0,:)=um
!       at(:,i0)=um ! or column if the solver wants the transpose
      Enddo                                                      
!     am=transpose(at)

!  Calculate equation to zero
      f=(y-yt)*rdt+ydot
!     f=y*rdt-yt*rdt+ydot+alph*ydot0/beta  
      If(idiag>=4) Then
        Write(lun_diag,"(a2,i5,es14.7)") 'F',kstep,rdt
        Do i=1,ny
          Write(lun_diag,"(a5,4es17.9)") 
     &          nname(i),f(i),ydot(i),yt(i),y(i)
          Write(lun_diag,"(5es16.8)") (am(i,j),j=1,ny)
        Enddo
      Endif

!  Test the eigenvalues 
!     If(idiag>=6) Then
!       call eigen_test(kstep,am,rdt)
!     Endif
!-----------------------------------------------------------------------------  
!  The bulk of the computational cost of the network (60-95%) is the solving 
!  of the matrix equation.  Careful selection of the matrix solver is therefore 
!  very important to fast computation.  Generally, hand tuned solvers such as 
!  those supplied by the hardware manufacturer or third-parties like NAG, IMSL,
!  etc. are the fastest.  However for portability, by default we use Numerical 
!  Recipes routines.  
!-----------------------------------------------------------------------------  

!  To use Num Rec LU Decomp, uncomment 3 lines below
!     call ludcmp(am,ny,ny,indx,d)
!     call lubksb(am,ny,ny,indx,f)
!     dy=f

!  To use LAPACK solver,  uncomment one of the first two lines below and the third line 
!     call sgesv(ny,1,am,ny,indx,f,ny,info) ! Single precision version
      call dgesv(ny,1,am,ny,indx,f,ny,info) ! Double precision version
      dy=f

!  Diagnostic output
!     If(idiag>=4) Then
!       Write(lun_diag,"(a5,3es12.4)") (nname(i),dy(i),y(i),yt(i),i=1,ny)
!     Endif
      Return                                                                    
      End subroutine netmatr                                                                       

