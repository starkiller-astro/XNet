!*******************************************************************************
!  Jacobian Sparce for PARDISO, part of Full_Net 4.10 4/20/07 
!
!  The routines in this file are used to replace the standard dense Jacobian 
!  and associated solver with a sparse Jacobian build andthe PARDISO sparse 
!  solver package.  
!  The modifications include a additional module with information about 
!  the sparcity pattern, a routine to read this data, and a modified netmatr. 
!*******************************************************************************

      module jacobian_data
!===============================================================================
!  Contains data for use in the sparse solver
!===============================================================================
      real*8, dimension(:), allocatable :: tvals,sident
      integer,dimension(:), allocatable :: cols,pb,perm
      integer,dimension(:), allocatable :: ridx,cidx
      integer,dimension(:), allocatable :: ns11,ns21,ns22
      integer,dimension(:), allocatable :: ns31,ns32,ns33
      integer :: lval,pt(64),iparm(64)
      integer :: l1s,l2s,l3s
      integer :: mtype,maxfct,mnum,phase,nrhs,msglvl
      end module jacobian_data
      
      subroutine read_jacobian_data(data_dir)
!===============================================================================
! Reads in data necessary to use sparse solver
!===============================================================================
      use nuclear_data
      use reac_rate_data
      use jacobian_data
           
      character (LEN=*),  INTENT(in)  :: data_dir
      integer :: i
      
      
      Open(600,file=trim(data_dir)//"/sparse_ind",
     &         status='old',form='unformatted')
      
      read(600)lval
      

      allocate(ridx(lval),cidx(lval))
      allocate(tvals(lval),sident(lval))
      allocate(pb(ny+1),perm(ny))
      
      ! Read compressed row format storage arrays 
      read(600)ridx,cidx,pb
        
      
      
      ! PARDISO Attributes
      mtype = 11
      maxfct= 1
      mnum  = 1
      nrhs  = 1
      msglvl= 0
      
      ! Build a compressed row format version of the identity matrix
      do i=1,lval
        if (ridx(i)==cidx(i))sident(i)=1.0
      enddo
      
      ! Read the sizes of the ns arrays
      read(600)l1s,l2s,l3s
      
      
      !Allocate arrays for direct sparse representation Jacobian build
      allocate(ns11(l1s))
      allocate(ns21(l2s),ns22(l2s))
      allocate(ns31(l3s),ns32(l3s),ns33(l3s))
      
      
      ns11 = 0
      ns21 = 0
      ns22 = 0
      ns31 = 0
      ns32 = 0
      ns33 = 0
      
      read(600)ns11,ns21,ns22
      read(600)ns31
      read(600)ns32
      read(600)ns33
      close(600)      
      
! If using the Basel PARDISO libraries, pardisoinit must be called before the first pardiso call
! If using MKL, do not call pardisoinit!

      iparm(1) = 0
      call pardisoinit(pt,mtype,iparm)
          
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
      use jacobian_data
      
      
      integer :: err,knrcut
      integer :: i,j,kstep,i0,la1,le1,la2,le2,la3,le3,j1,l1,l2,l3
      integer :: ls1,ls2,ls3
      real(8) :: f(ny)
      real(8) :: alph,beta,rdt,d,det_sm,det_fm,difference
      
      knrcut = int(knrmx*0.5)
      alph=0.0                                                                  
      beta=1.0                                                                  
      
! Calculate the reaction rates and abundance time derivatives
      call yderiv

!  The quick solution to taking advantage of sparseness is to create a values array that has the maximum
!  number of non-zero elements as well as a 2-by-#non-zero elements matrix and the other vectors
!  required by the sparse solver.  The second matrix will contain the ordered pairs that map to a 
!  particular place in the Jacobian.  The Jacobian is built by looping over rows with a nested loop over 
!  reactions which effect that row.  The sparse tvals array is built directly by using sparse mapping arrays
!  determined in read_sparse_data.  This saves significant time, but more importantly saves vast amounts of space.

      tvals = 0.0
!  Build the Jacobian, row by row
      rdt=1.0/tdel/beta                        
      
      ! Set kronecker delta portion of the Jacobian in the sparse matrix representation
      tvals = tvals + sident * rdt
      
     
      do i0=1,ny
        la1=la(1,i0)
        le1=le(1,i0)
        do j1=la1,le1
          ls1 = ns11(j1) !ns11(j1) gives the index effected reaction j1 by in the compressed row storage scheme
          tvals(ls1)=tvals(ls1)-b1(j1)
        enddo
        la2=la(2,i0) 
        le2=le(2,i0)  
        do j1=la2,le2
          ls1=ns21(j1) !ns21(j1) gives the first index effected reaction j1 by in the compressed row storage scheme
          ls2=ns22(j1) !ns22(j1) gives the second index effected reaction j1 by in the compressed row storage scheme
          l1=n21(j1) !n21(k) gives the index of first reactant in reaction mu2(k)
          l2=n22(j1) !n22(k) gives the index of second reactant in reaction mu2(k)
          tvals(ls1)=tvals(ls1)-b2(j1)*yt(l2)
          tvals(ls2)=tvals(ls2)-b2(j1)*yt(l1)
        enddo
        la3=la(3,i0)
        le3=le(3,i0)
        do j1=la3,le3
          ls1=ns31(j1) !ns31(j1) gives the first index effected reaction j1 by in the compressed row storage scheme
          ls2=ns32(j1) !ns32(j1) gives the second index effected reaction j1 by in the compressed row storage scheme
          ls3=ns33(j1) !ns33(j1) gives the third index effected reaction j1 by in the compressed row storage scheme
          l1=n31(j1) !n21(k) gives the index of first reactant in reaction mu2(k)
          l2=n32(j1) !n22(k) gives the index of second reactant in reaction mu2(k)
          l3=n33(j1) !n22(k) gives the index of third reactant in reaction mu2(k)
          tvals(ls1)=tvals(ls1)-b3(j1)*yt(l2)*yt(l3)
          tvals(ls2)=tvals(ls2)-b3(j1)*yt(l1)*yt(l3)
          tvals(ls3)=tvals(ls3)-b3(j1)*yt(l1)*yt(l2)
        enddo
      enddo                                            
      

!  Calculate equation to zero
      f=(y-yt)*rdt+ydot      
 

! Use PARDISO sparse solver from MKL
      ! Check if this is the first step, or if convergence is taking to long
      if (kstep == 1.or.kts>1.or.knr>knrcut) then
          iparm(1) = 0
          phase = 11
          call pardiso(pt,maxfct,mnum,mtype,phase,
     &                ny,tvals,pb,cidx,perm,nrhs, 
     &                iparm,msglvl,f,dy,err)
          if (err/=0) print *,'PARDISO error ',err
      endif
        
        phase = 23
        iparm(1) = 0 
        call pardiso(pt,maxfct,mnum,mtype,phase,
     &                ny,tvals,pb,cidx,perm,nrhs, 
     &                iparm,msglvl,f,dy,err)
        if (err/=0) print *,'PARDISO error ',err
        
      Return
      End subroutine netmatr

