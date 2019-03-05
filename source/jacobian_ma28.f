!*******************************************************************************
!  Jacobian Sparce for MA28, part of Full_Net 4.10 4/20/07 
!
!  The routines in this file are used to replace the standard dense Jacobian 
!  and associated solver with the Harwell MA28 sparse solver package.  
!  The modifications include a additional module with information about 
!  the sparcity pattern, a routine to read this data, and a modified netmatr. 
!*******************************************************************************

      module jacobian_data
!===============================================================================
!  Contains data for use in the sparse solver
!===============================================================================
      real*8, dimension(:), allocatable :: vals,ma28w,tvals,sident
      integer,dimension(:), allocatable :: ridx,cidx,pb
      integer,dimension(:), allocatable :: irn,icn
      integer,dimension(:), allocatable :: ns11,ns21,ns22
      integer,dimension(:), allocatable :: ns31,ns32,ns33
      integer,dimension(:), allocatable :: ikeep,iwA,iwB
      integer :: lval,lirn,licn
      real*8  :: pivoting = 0.1
      end module jacobian_data
      
      subroutine read_jacobian_data(data_dir)
!----------------------------------------------------------------------------
! Reads in data necessary to use sparse solver
!----------------------------------------------------------------------------
      use nuclear_data
      use reac_rate_data
      use jacobian_data
           
      character (LEN=*),  INTENT(in)  :: data_dir
      integer :: la1,la2,la3,le1,le2,le3,l1,l2,l3
      integer :: i0,i1,j1
      
      COMMON /MA28ED/ LP,MP,LBLOCK,GROW
      COMMON /MA28FD/EPS,RMIN,RESID,IRNCP,ICNCP,MINIRN,MINICN,IRANK,
     &             ABORT1,ABORT2
      DOUBLE PRECISION EPS,MRESID,RESID,RMIN
      INTEGER LP,MP,ICNCP,IRANK,IRNCP,MINICN,MINIRN
      LOGICAL LBLOCK,GROW,ABORT1,ABORT2
      
! Do not abort MA28AD on numerically singular matrices
      ABORT1 = .false.
      ABORT2 = .false. 
      
      LP = 99
      MP = 99
      open(99,file='MA28.err')
      
      Open(600,file=trim(data_dir)//"/sparse_ind",
     &         status='old',form='unformatted')
      
      read(600)lval
      
      licn = 8*lval
      lirn = 7*lval

      allocate(ridx(lval),cidx(lval))
      allocate(tvals(lval),sident(lval))
      allocate(pb(ny+1))
      
! These are variables to be used by the MA28 solver
      allocate(vals(licn),icn(licn),irn(lirn))
      allocate(ma28w(ny))
      allocate(iwA(8*ny),iwB(5*ny),ikeep(5*ny))
      
      read(600)ridx,cidx,pb
      
      ! Build a compressed row format version of the identity matrix
      do i=1,lval
        if (ridx(i)==cidx(i))sident(i)=1.0
      enddo
      
      read(600)l1s,l2s,l3s
      
      
      !Build  arrays for direct sparse representation Jacobian build
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
      
      end subroutine read_jacobian_data

      subroutine netmatr(kstep)
!===============================================================================
!  This routine calculates the Jacobian matrix dYdot/dY, and solves for the 
!  Newton-Raphson iteration, dy using the MA28 sparse matrix package. 
!===============================================================================
      use controls
      use nuclear_data
      use conditions
      use abundances
      use reac_rate_data
      use jacobian_data
      
      
      integer :: err
      integer :: i,j,kstep,i0,la1,le1,la2,le2,la3,le3,j1,l1,l2,l3
      integer :: ls1,ls2,ls3
      real(8) :: f(ny)
      real(8) :: alph,beta,rdt,d,det_sm,det_fm,difference
      character(6) :: matdescra
      COMMON /MA28FD/EPS,RMIN,RESID,IRNCP,ICNCP,MINIRN,MINICN,IRANK,
     &             ABORT1,ABORT2
      COMMON /MA28ED/ LP,MP,LBLOCK,GROW
      COMMON /MA28HD/TOL,THEMAX,BIG,DXMAX,ERRMAX,DRES,CGCE,NDROP, 
     &             MAXIT,NOITER,NSRCH,ISTART,LBIG
      DOUBLE PRECISION EPS,MRESID,RESID,RMIN
      INTEGER ICNCP,IRANK,IRNCP,MINICN,MINIRN,LP,MP
      LOGICAL ABORT1,ABORT2,LBLOCK,GROW
      DOUBLE PRECISION TOL,THEMAX,BIG,DXMAX,ERRMAX,DRES,CGCE 
      INTEGER NDROP,MAXIT,NOITER,NSRCH,ISTART
      LOGICAL LBIG 


! Reduce the pivot ratio      
      EPS = 1.0e-6
      
      alph=0.0                                                                  
      beta=1.0                                                                  
      
! Calculate the reaction rates and abundance time derivatives
      call yderiv

!  The quick solution to taking advantage of sparseness is to create a values array that has the maximum
!  number of non-zero elements as well as a 2-by-#non-zero elements matrix and the other vectors
!  required by the sparse solver.  The second matrix will contain the ordered pairs that map to a 
!  particular place in the Jacobian.  Build the Jacobian as it is built now, then pull the values 
!  from it using the ordered pairs to fill the values array.  Now use the sparse solver to solve the 
!  linear set of equations.  Doing by just checking for zero values in the Jacobian is hard to do because
!  one must continually allocate and deallocate arrays, and save multiple arrays for transfering purposes.
!  I would guess that this removes all the advantage of a sparse solver in all but the most extreme cases. 

!  Build the Jacobian, row by row
      rdt=1.0/tdel/beta                        
      
      tvals = 0.0
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
      

! Use Harwell MA28 sparse solver (this one is portable)
        if (kstep == 1) then
          icn = 0
          irn = 0
          vals = 0.0
          iwA = 0
          ma28w = 0.0
          ikeep = 0
          iwB = 0
          err = 0
          do i=1,lval
            icn(i) = cidx(i)
            irn(i) = ridx(i)
          enddo
          vals(1:lval) = tvals

          call MA28AD(ny,lval,vals,licn,irn,lirn,icn,
     &      pivoting,ikeep,iwA,ma28w,err)
     
          fs = 0
        endif

! Build sparse value array        
        vals(1:lval) = tvals


! Do the numerical decomposition using the previously determined symbolic decomposition
          call MA28BD(ny,lval,vals,licn,ridx,cidx,
     &       icn,ikeep,iwB,ma28w,err)
          if (err > 0) then
            write(LP,*)'Error during MA28BD.  Code: ',err,
     &      ' Name: ',nname(err)
          elseif (err /= 0) then
            write(LP,*)'Error during MA28BD.  Code: ',err
          endif
          
! Perform check to see if an error occured in MA28BD.  Most likely a singularity error
! caused by incorrect symbolic matrix, so run MA28AD again.  Also, make sure that 
! workspaces are large enough.
        if (ICNCP > ny/10) err = -66
        if (IRNCP > ny/10) err = -66
        
        do iter = 1,10
          if (err == 0) exit
          if (iter>1) write(LP,*)'2. Error during MA28AD.  Code: ',err
          if (err < -2) then
            write(LP,*)'Trying to reallocate MA28 arrays'
            if (MINICN > licn) then
               write(LP,*)'MINICN = ',MINICN,' > licn = ',licn
               deallocate(icn,vals)
               licn = int(1.2*MINICN)
               allocate(icn(licn),vals(licn))
             endif
             if (MINIRN > lirn) then
              write(LP,*)'MINIRN = ',MINIRN,' > lirn = ',lirn
              deallocate(irn)
              lirn = int(1.2*MINIRN)
              allocate(irn(lirn))
            endif   
            if (IRNCP > ny/10) then
              write(LP,*)'IRNCP = ',IRNCP,' > ny/10 = ',ny/10
              deallocate(irn)
              lirn = lirn + int(lirn*0.5)
              allocate(irn(lirn))
            endif      
            if (ICNCP > ny/10) then
              write(LP,*)'ICNCP = ',ICNCP,' > ny/10 = ',ny/10
              deallocate(icn,vals)
              licn = licn + int(licn*0.5)
              allocate(icn(licn),vals(licn))
            endif
        endif
          
        ! Rebuild vals array
        icn = 0
        irn = 0
        vals = 0.0
        iwA = 0
        ma28w = 0.0
        ikeep = 0
        iwB = 0
        err = 0
        do i=1,lval
          icn(i) = cidx(i)
          irn(i) = ridx(i)
        enddo
        vals(1:lval)=tvals
        ikeep = 0
        iwA = 0
        ma28w = 0.0
        call MA28AD(ny,lval,vals,licn,irn,lirn,icn,
     &     pivoting,ikeep,iwA,ma28w,err)
     
      enddo
        
! Perform back substitution 
      call MA28CD(ny,vals,licn,icn,ikeep,f,ma28w,1)
      dy = f
        
      Return
      End subroutine netmatr
      
