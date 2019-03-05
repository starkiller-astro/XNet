      SUBROUTINE ludcmp (a,n,np,indx,d)

c   This subroutine from numerical recipes, vol. 2, pp. 38-39.
c
c   Given a matrix a(1:n,1:n), with physical dimension np by np, this 
c   routine replaces it by
c   the LU decomposition of a rowwise permutation of itself.  a and n 
c   are input.  a is output,
c   arranged as in equation (2.3.14) above; indx(1:n) is an output vector 
c   that records the
c   row permutation effected by the partial pivoting; d is output as +/- 1 
c   depending on whether
c   the number of row interchanges was even or odd, respectively.  This 
c   routine is used in
c   combination with lubksb to solve linear equations or invert a 
c   matrix.
      
      INTEGER(4) n,np,indx(n),NMAX
      REAL(8) d,a(np,np),TINY
      PARAMETER (NMAX=1000,TINY=1.0e-20) ! Largest expected n, & a small no.
      INTEGER(4) i,imax,j,k
      REAL(8) aamax,dum,sum,vv(NMAX)  !vv stores implicit scaling of each row.

      d=1.                         !No row interchanges yet.

      do 12 i=1,n                  !Loop over rows to get the implicit 
                                   !scaling information.
         aamax=0.

         do 11 j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
 11      continue

         if (aamax.eq.0.) pause 'singular matrix in ludcmp' 
                                   !No nonzero largest element.

         vv(i)=1./aamax            !Save the scaling.

 12   continue

      do 19 j=1,n          !This is the loop over columns of Crout's method.

         do 14 i=1,j-1             !This is equation (2.3.12) except  i = j.
            sum=a(i,j)

            do 13 k=1,i-1
               sum=sum-a(i,k)*a(k,j)
 13         continue

            a(i,j)=sum

 14      continue

         aamax=0.            !Initialize for search for largest pivot element.

         do 16 i=j,n         !This is i = j of eq. (2.3.12) and i = j + 1...N
            sum=a(i,j)       !of equation (2.3.13).

            do 15 k=1,j-1
               sum=sum-a(i,k)*a(k,j)
 15         continue

            a(i,j)=sum
            dum=vv(i)*abs(sum)          !Figure of merit for the pivot.
            if (dum.ge.aamax) then      !Is it better than the best so far?
               imax=i
               aamax=dum
            endif

 16      continue

         if (j.ne.imax)then            !Do we need to interchange rows?

            do 17 k=1,n                !Yes, do so...
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
 17         continue

            d=-d                        !...and change the parity of d.
            vv(imax)=vv(j)              !Also interchange the scale factor.
         endif

         indx(j)=imax
         if(a(j,j).eq.0.)a(j,j)=TINY

c   If the pivot element is zero the matrix is singular (at least to the 
c   precision of the algorithm).
c   For some applications on singular matrices, it is desirable to 
c   substitute TINY for zero.

         if(j.ne.n) then      !Now, finally, divide by the pivot element.
            dum=1.0d0/a(j,j)

            do 18 i=j+1,n
               a(i,j)=a(i,j)*dum
 18         continue

         endif

 19   continue                  !Go back for the next column in the reduction.

      return
      END

      SUBROUTINE lubksb(a,n,np,indx,b)

c   This subroutine from numerical recipes, vol. 2, pp. 39.
c
c   Solves the set of n linear equations A * X = B.  Here a is input, not 
c   as the matrix A but
c   rather as its LU decomposition, determined by the routine ludcmp.  
c   indx is input as the 
c   permutation vector returned by ludcmp.  
c   b(1:n) is input as the right-hand side vector B,
c   and returns with the solution vector X.  a,n,np, and indx are not 
c   modified by this routine
c   and can be left in place for successive calls with different right-
c   hand sides b.  This routine
c   takes into account the possibility that b will begin with many zero 
c   elements, so it is efficient
c   for use in matrix inversion.

      INTEGER(4) n,np,indx(n)
      REAL(8) a(np,np),b(n)
      INTEGER(4) i,ii,j,ll
      REAL(8) sum

      ii=0                        ! When ii is set to a posititve value, 
                                  ! it will become the index
      do 12 i=1,n                 ! of the first nonvanishing 
                                  ! element of b.  We now do
         ll=indx(i)               ! the forward substitution, 
                                  ! equation (2.3.6).  The only new
         sum=b(ll)                ! wrinkle is to unscramble the 
                                  ! permutation as we go.
         b(ll)=b(i)

         if (ii.ne.0)then

            do 11 j=ii,i-1
               sum=sum-a(i,j)*b(j)
 11         continue

         else if (sum.ne.0.) then
            ii=i                  !A nonzero element was encountered, 
                                  ! so from now on we will
         endif                    !have to do the sums in the loop above.

         b(i)=sum

 12   continue

      do 14 i=n,1,-1            !Now we do the backsubstitution, eq. (2.3.7).

         sum=b(i)

         do 13 j=i+1,n
            sum=sum-a(i,j)*b(j)
 13      continue

         b(i)=sum/a(i,i)        !Store a component of the solution vector X.

 14   continue

      return                    !All done!
      END

