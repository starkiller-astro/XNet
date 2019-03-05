!***************************************************************************************************
! xnet_preprocess.f90 10/18/17
! Reaction network data preprocessing
!
! These routines take the human readable, ASCII data files sunet, netsu, netweak, netneutr, and
! netwinv and prepares binary versions for the reaction network. This involves translating
! characters to indicies and organizing the data for more efficient computation. Several additional
! tests are performed.
!***************************************************************************************************

Module xnet_preprocess
  !-------------------------------------------------------------------------------------------------
  ! This module contains the reaction rate and flux data, including the indices which type the
  ! reaction and the cross-references, as well as the routines to preprocess the ASCII data files.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None

  ! Reaction data
  Integer, Parameter   :: nchap = 11 ! Number of REACLIB chapters
  Integer, Parameter   :: nnucr(nchap)      = (/ 2, 3, 4, 3, 4, 5, 6, 4, 5, 6, 5 /) ! Number of nuclei for each REACLIB rate chapter
  Integer, Parameter   :: nreactant(nchap)  = (/ 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 1 /) ! Number of reactants for each REACLIB rate chapter
  Integer, Parameter   :: n1m = 5, n2m = 6, n3m = 6, n4m = 6 ! Maximum number of nuclei for i-reactant chapters
  Integer, Parameter   :: nmax = max(n1m,n2m,n3m,n4m)        ! Maximum number of nuclei for any REACLIB rate
  Integer, Parameter   :: nrm1 = 100000             ! Maximum number of 1-reactant reactions
  Integer, Parameter   :: nrm2 = 100000             ! Maximum number of 2-reactant reactions
  Integer, Parameter   :: nrm3 = 100                ! Maximum number of 3-reactant reactions
  Integer, Parameter   :: nrm4 = 100                ! Maximum number of 3-reactant reactions
  Integer, Parameter   :: nrmax = nrm1 + nrm2 + nrm3 + nrm4
  Real(dp)             ::      q1(nrm1),      q2(nrm2),      q3(nrm3),      q4(nrm4)
  Integer              ::    iwk1(nrm1),    iwk2(nrm2),    iwk3(nrm3),    iwk4(nrm4)
  Integer              ::   ires1(nrm1),   ires2(nrm2),   ires3(nrm3),   ires4(nrm4)
  Integer              ::   irev1(nrm1),   irev2(nrm2),   irev3(nrm3),   irev4(nrm4)
  Character(4)         ::   desc1(nrm1),   desc2(nrm2),   desc3(nrm3),   desc4(nrm4)
  Integer              ::  n1(n1m,nrm1),  n2(n2m,nrm2),  n3(n3m,nrm3),  n4(n4m,nrm4)
  Real(dp)             :: an1(n1m,nrm1), an2(n2m,nrm2), an3(n3m,nrm3), an4(n4m,nrm4)
  Integer              :: mu1(n1m*nrm1), mu2(n2m*nrm2), mu3(n3m*nrm3), mu4(n4m*nrm4)
  Real(dp)             ::  a1(n1m*nrm1),  a2(n2m*nrm2),  a3(n3m*nrm3),  a4(n4m*nrm4)
  Real(dp)             ::   rc1(7,nrm1),   rc2(7,nrm2),   rc3(7,nrm3),   rc4(7,nrm4)
  Integer              :: la(nchap), le(nchap)
  Integer, Allocatable :: l1a(:), l2a(:), l3a(:), l4a(:)
  Integer, Allocatable :: l1e(:), l2e(:), l3e(:), l4e(:)
  Integer              :: nreac(4), nan(4)
  Integer              :: lun_s3, lun_s4, lun_ndiag, lun_desc

  ! Flux data
  Integer, Parameter :: mxflx = 60000 ! Maximum number of unique reaction pairs
  Integer            :: nflx(8,mxflx) ! The nuclei in each unique reaction pair
  Integer            :: iwflx(mxflx)  ! Reaction pair weak flag
  Real(dp)           :: qflx(mxflx)   ! Reaction pair Q values
  Character(4)       :: descx(mxflx)  ! Descriptor for the reaction pair

Contains

  Subroutine flux_search(iflx,iw,i8,q,desc,mflx)
    !-----------------------------------------------------------------------------------------------
    ! This routine searches for a given reaction among the matched reaction pairs (fluxes),
    ! returning the flux index
    !-----------------------------------------------------------------------------------------------
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Integer, Intent(in)      :: iw    ! Weak reaction flag
    Integer, Intent(in)      :: i8(8) ! Nuclear indicies to be matched
    Real(dp), Intent(in)     :: q     ! Reaction Q value
    Character(4), Intent(in) :: desc  ! Reaction source descriptor

    ! Input/Output variables
    Integer, Intent(inout) :: mflx  ! Number of reaction pairs

    ! Output variables
    Integer, Intent(out) :: iflx  ! Search result

    ! Local variables
    Integer :: m

    ! Check reaction against existing fluxes
    iflx = 0
    !Write(lun_ndiag,"(a6,9i5)") 'Find',iflx,i8
    If ( iw == 0 .and. mflx /= 0 ) Then
      Do m = 1, mflx
        If ( all( i8 == nflx(:,m) ) ) Then
          iflx = m
          Exit
        EndIf
      EndDo
    EndIf

    ! If flux does not match existing flux, create new flux
    If ( iflx == 0 ) Then
      mflx         = mflx + 1
      nflx(:,mflx) = i8
      qflx(mflx)   = abs(q)
      descx(mflx)  = desc
      iwflx(mflx)  = iw
      iflx         = mflx
      !Write(lun_ndiag,"(a6,9i5)") 'NFLX+',mflx,nflx(:,mflx)
    Else
      !Write(lun_ndiag,"(a6,9i5)") 'Found',iflx,nflx(:,iflx)
    EndIf

    Return
  End Subroutine flux_search

  Subroutine match_react(data_dir)
    !-----------------------------------------------------------------------------------------------
    ! In this data set, the forward and reverse reactions are separate. For the purposes of this
    ! routine, forward reactions are defined as those with positive Q values. Additionally, some
    ! reactions have multiple components. This routine matches these multiple components of forward
    ! and reverse rates. While this data is unnecessary for the basic functioning of the network, it
    ! is useful for some analysis. For example, to calculate the net flux through a reaction
    ! channel, one must sum the components.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: nname, nn, zz
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: data_dir

    ! Local variables
    Integer      :: mflx ! the number of reaction sets
    Integer      :: iflx ! the reaction set index
    Integer      :: i,j,k,ii ! Loop variables
    Integer      :: itest(8),nprod,nreact,nnuc,iw
    Integer      :: ifl1(nrm1),ifl2(nrm2),ifl3(nrm3),ifl4(nrm4)
    Integer      :: revtest,comtest
    Integer      :: ifl_orig,ifl_term ! indicies of originating and termination species
    Real(dp)     :: q,rvar
    Character(5) :: blank5,arrow
    Character(4) :: desc
    Integer      :: lun_matchd, lun_matchr

    mflx = 0

    Do j = 1, nchap
      nprod  = nnucr(j) - nreactant(j)
      nreact = nreactant(j)
      nnuc   = nnucr(j)

      ! Match the one reactant reactions
      If ( nreact == 1 ) Then
        Do k = la(j), le(j)
          q     = q1(k)
          iw    = iwk1(k)
          desc  = desc1(k)
          itest = 0
          If ( q > 0.0 ) Then
            itest(1:nreact)  = n1(1:nreact,k)
            itest(5:nprod+4) = n1(nreact+1:nnuc,k)
            Call flux_search(iflx,iw,itest,q,desc,mflx)
            ifl1(k) = iflx
          Else
            itest(1:nprod)    = n1(nreact+1:nnuc,k)
            itest(5:nreact+4) = n1(1:nreact,k)
            Call flux_search(iflx,iw,itest,q,desc,mflx)
            ifl1(k) = -iflx
          EndIf
        EndDo

        ! Match the two reactant reactions
      ElseIf ( nreact == 2 ) Then
        Do k = la(j), le(j)
          q     = q2(k)
          iw    = iwk2(k)
          desc  = desc2(k)
          itest = 0
          If ( q > 0.0 ) Then
            itest(1:nreact)  = n2(1:nreact,k)
            itest(5:nprod+4) = n2(nreact+1:nnuc,k)
            Call flux_search(iflx,iw,itest,q,desc,mflx)
            ifl2(k) = iflx
          Else
            itest(1:nprod)    = n2(nreact+1:nnuc,k)
            itest(5:nreact+4) = n2(1:nreact,k)
            Call flux_search(iflx,iw,itest,q,desc,mflx)
            ifl2(k) = -iflx
          EndIf
        EndDo

        ! Match three reactant reactions
      ElseIf ( nreact == 3 ) Then
        Do k = la(j), le(j)
          q     = q3(k)
          iw    = iwk3(k)
          desc  = desc3(k)
          itest = 0
          If ( q > 0.0 ) Then
            itest(1:nreact)  = n3(1:nreact,k)
            itest(5:nprod+4) = n3(nreact+1:nnuc,k)
            Call flux_search(iflx,iw,itest,q,desc,mflx)
            ifl3(k) = iflx
          Else
            itest(1:nprod)=n3(nreact+1:nnuc,k)
            itest(5:nreact+4)=n3(1:nreact,k)
            Call flux_search(iflx,iw,itest,q,desc,mflx)
            ifl3(k) = -iflx
          EndIf
        EndDo

        ! Match four reactant reactions
      ElseIf ( nreact == 4 ) Then
        Do k = la(j), le(j)
          q     = q4(k)
          iw    = iwk4(k)
          desc  = desc4(k)
          itest = 0
          If ( q > 0.0 ) Then
            itest(1:nreact)  = n4(1:nreact,k)
            itest(5:nprod+4) = n4(nreact+1:nnuc,k)
            Call flux_search(iflx,iw,itest,q,desc,mflx)
            ifl4(k) = iflx
          Else
            itest(1:nprod)=n4(nreact+1:nnuc,k)
            itest(5:nreact+4)=n4(1:nreact,k)
            Call flux_search(iflx,iw,itest,q,desc,mflx)
            ifl4(k) = -iflx
          EndIf
        EndDo
      EndIf
    EndDo

    ! Output reaction matching
    Open(newunit=lun_matchd,file=trim(data_dir)//'/match_data',FORM='unformatted')
    Write(lun_matchd) mflx,(nreac(k),k=1,4)
    Write(lun_matchd) ifl1(1:nreac(1)),ifl2(1:nreac(2)),ifl3(1:nreac(3)),ifl4(1:nreac(4))
    Write(lun_matchd) nflx(:,1:mflx),qflx(1:mflx),iwflx(1:mflx),descx(1:mflx)
    Close(lun_matchd)

    ! Output ASCII match info
    arrow = ' --> '
    rvar = 0.0
    Open(newunit=lun_matchr, file=trim(data_dir)//'/match_read')
    Write(lun_matchr,"(10a5,1es10.3)") ((nname(nflx(j,i)),j=1,4),arrow,(nname(nflx(k,i)),k=5,8),descx(i),rvar,i=1,mflx)

    ! Test for matching components, reverse reactions and output
    blank5 = '     '
    Do ii = 1, mflx
      revtest = 0
      comtest = 0
      Write(lun_ndiag,"(a2)") '--'
      Do k = 1, nreac(1)
        If ( abs(ifl1(k)) == ii ) Then
          Write(lun_ndiag,"(9a5,2i2,i6,1es12.4)") &
            & nname(n1(1,k)),blank5,blank5,arrow,(nname(n1(j,k)),j=2,4),blank5,desc1(k),ires1(k),irev1(k),ifl1(k),q1(k)
          comtest = comtest + 1
          revtest = revtest + irev1(k)
        EndIf
      EndDo
      Do k = 1, nreac(2)
        If ( abs(ifl2(k)) == ii ) Then
          Write(lun_ndiag,"(9a5,2i2,i6,1es12.4)") &
            & (nname(n2(i,k)),i=1,2),blank5,arrow,(nname(n2(j,k)),j=3,6),desc2(k),ires2(k),irev2(k),ifl2(k),q2(k)
          comtest = comtest + 1
          revtest = revtest + irev2(k)
        EndIf
      EndDo
      Do k = 1, nreac(3)
        If ( abs(ifl3(k)) == ii ) Then
          Write(lun_ndiag,"(9a5,2i2,i6,1es12.4)") &
            & (nname(n3(i,k)),i=1,3),arrow,(nname(n3(j,k)),j=4,5),blank5,blank5,desc3(k),ires3(k),irev3(k),ifl3(k),q3(k)
          comtest = comtest + 1
          revtest = revtest + irev3(k)
        EndIf
      EndDo
      Do k = 1, nreac(4)
        If ( abs(ifl4(k)) == ii ) Then
          Write(lun_ndiag,"(10a5,2i2,i6,1es12.4)") &
            & (nname(n4(i,k)),i=1,4),arrow,(nname(n4(j,k)),j=5,6),blank5,blank5,desc4(k),ires4(k),irev4(k),ifl4(k),q4(k)
          comtest = comtest + 1
          revtest = revtest + irev4(k)
        EndIf
      EndDo

      ! Write test results
      If ( mod(comtest,2) /= 0 .and. iwflx(ii) == 0 ) Write(lun_ndiag,"(a,i4)") 'Unbalanced reaction',comtest
      If ( comtest > 1 .and. revtest == 0 ) Write(lun_ndiag,"(a,i4)") 'No reverse defined',revtest

      ! Write flux output formatted for Matlab
      ifl_orig = nflx(count(nflx(1:4,ii) /= 0 ),  ii)
      ifl_term = nflx(count(nflx(5:8,ii) /= 0 )+4,ii)
      Write(lun_matchr,'(i5,4f6.1,es13.5,a5)') &
        & ii,zz(ifl_orig),nn(ifl_orig),zz(ifl_term),nn(ifl_term),1.0,descx(ii)
    EndDo
    Close(lun_matchr)

    Return
  End Subroutine match_react

  Subroutine sparse_check(data_dir)
    !-----------------------------------------------------------------------------------------------
    ! This routine maps the network Jacobian and calculates parameters for the sparseness. The
    ! network Jacobian is reasonably described as a doubly bordered band diagonal, hough even this
    ! structure contains considerable zero elements.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny, nname
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: data_dir

    ! Local variables
    Real(dp) :: am(ny,ny) ! Dense Jacobian
    Real(dp) :: um(ny)    ! Jacobian Row
    Integer :: col(ny*ny), row(ny*ny) ! Column and row positions of non-zero elements
    Integer :: llum(ny) ! List of non-zero row elements
    Integer :: pb(ny+1) ! First column index for each row in row vector
    Integer :: pe(ny)   ! Last column index for each row in row vector
    Integer :: ns11(nan(1)), ns21(nan(2)), ns22(nan(2)) ! Map of sparse matrix location
    Integer :: ns31(nan(3)), ns32(nan(3)), ns33(nan(3)) ! for each reaction
    Integer :: ns41(nan(4)), ns42(nan(4)), ns43(nan(4)), ns44(nan(4))
    Integer :: vl, i0, i1
    Integer :: i, j, j1, l1, l2, l3, l4
    Integer :: nlum, leftm, imin, imax, iwhere, jmin, jmax, jwhere, ijd
    Integer :: lun_shape, lun_sparse

    ! The file matr_shape lists the coordinates of the non-zero elements in a form suitable for plotting.
    ! Output the matrix dimension and trial right hand side
    Open(newunit=lun_shape, file=trim(data_dir)//'/matr_shape')
    Write(lun_shape,"(a3,i4)") 'NY=',ny
    um = 1.0
    Write(lun_shape,"(8es10.3)") um

    ! Build the test matrix
    am = 0.0
    Do i0 = 1, ny
      um = 0.0
      um(i0) = 10.0
      Do j1 = l1a(i0), l1e(i0)
        l1 = n1(1,mu1(j1))
        um(l1) = um(l1) + 1.0
      EndDo
      Do j1 = l2a(i0), l2e(i0)
        l1 = n2(1,mu2(j1))
        l2 = n2(2,mu2(j1))
        um(l1) = um(l1) + 1.0
        um(l2) = um(l2) + 1.0
      EndDo
      Do j1 = l3a(i0), l3e(i0)
        l1 = n3(1,mu3(j1))
        l2 = n3(2,mu3(j1))
        l3 = n3(3,mu3(j1))
        um(l1) = um(l1) + 1.0
        um(l2) = um(l2) + 1.0
        um(l3) = um(l3) + 1.0
      EndDo
      Do j1 = l4a(i0), l4e(i0)
        l1 = n4(1,mu4(j1))
        l2 = n4(2,mu4(j1))
        l3 = n4(3,mu4(j1))
        l4 = n4(4,mu4(j1))
        um(l1) = um(l1) + 1.0
        um(l2) = um(l2) + 1.0
        um(l3) = um(l3) + 1.0
        um(l4) = um(l4) + 1.0
      EndDo
      Write (lun_ndiag,*) i0
      am(i0,:) = um

      ! Check for non-zero elements
      nlum = 0
      Do i1 = 1, ny
        If ( abs(um(i1)) > tiny(0.0) ) Then
          nlum = nlum + 1
          llum(nlum) = i1
          Write(lun_shape,"(i4,2x,i4,es11.3)") i0,i1,um(i1)
        EndIf
      EndDo
      Write(lun_ndiag,"(18i4)") (llum(i),i=1,nlum)
    EndDo
    Close(lun_shape)

    ! Find heaviest incident nucleus, which sets border width
    leftm = 0
    Do i = 1, ny
      If ( nname(i) == '    n' .or. nname(i) == '    p' .or. nname(i) == '    d' .or. nname(i) == '    t' .or. &
        &  nname(i) == '  he3' .or. nname(i) == '  he4' .or. nname(i) == '  c12' .or. nname(i) == '  n14' .or. &
        &  nname(i) == '  o16' .or. nname(i) == ' ne20' ) leftm = i
    EndDo

    ! Find band widths
    imax = 0
    jmax = 0
    imin = 0
    jmin = 0
    ijd  = leftm + 1

    ! From the diagonal, for each row beyond the border
    Do i = ijd, ny

      ! Search rightmost (lowermost) elements of the band diagonal
      Do j = 1, ny-i
        If ( abs(am(i,i+j)) > tiny(0.0) .and. j > jmax ) Then
          jmax = j
          jwhere = i
        EndIf
        If ( abs(am(i+j,i)) > tiny(0.0) .and. j > imax ) Then
          imax = j
          iwhere = i
        EndIf
      EndDo

      ! Search for the leftmost (uppermost) element of the band diagonal
      Do j = 1, i-1
        If ( i-j > leftm .and. abs(am(i,i-j)) > tiny(0.0) .and. j > jmin ) Then
          jmin = j
          If (i-j <= leftm ) jmin = i - leftm - 1
        EndIf
        If ( i-j > leftm .and. abs(am(i-j,i)) > tiny(0.0) .and. j > imin ) Then
          imin = j
          If ( i-j <= leftm ) imin = i - leftm - 1
        EndIf
      EndDo
    EndDo

    ! Output parameters describing doubly-bordered band diagonal form.
    Write(lun_s4) leftm,leftm,ijd,imin,imax
    Write(lun_desc,"('Matrix Sparseness parameters')")
    Write(lun_desc,"('Border Widths   ',2i5)") leftm,leftm
    Write(lun_desc,"('Diagonal Widths ',2i5,'(j<i,j>i)')") jmin,jmax

    !-------------------------------------------------------------------------------------------------
    ! Sparse solver methods like PARDISO and MA48 use Compressed Row storage (CRS) to efficiently
    ! store the matrix to be solved.
    !-------------------------------------------------------------------------------------------------
    ! Build the Jacobian in CRS
    col(:) = 0
    row(:) = 0
    pb(:) = 0
    pe(:) = 0
    vl = 0
    Do i0 = 1, ny
      pb(i0) = vl + 1
      Do i1 = 1, ny
        If ( abs(am(i0,i1)) > tiny(0.0) ) Then
          vl = vl + 1
          col(vl) = i1
          row(vl) = i0
        EndIf
      EndDo
      pe(i0) = vl
    EndDo
    pb(ny+1) = vl + 1

    ! Build arrays to map reaction rates into locations in the CRS Jacobian
    Do i0 = 1, ny
      Do j1 = l1a(i0), l1e(i0)
        l1 = n1(1,mu1(j1))
        Do i1 = 1, vl
          If ( col(i1) == l1 .and. row(i1) == i0 ) ns11(j1) = i1
        EndDo
      EndDo
      Do j1 = l2a(i0), l2e(i0)
        l1 = n2(1,mu2(j1))
        l2 = n2(2,mu2(j1))
        Do i1 = 1, vl
          If ( col(i1) == l1 .and. row(i1) == i0 ) ns21(j1) = i1
          If ( col(i1) == l2 .and. row(i1) == i0 ) ns22(j1) = i1
        EndDo
      EndDo
      Do j1 = l3a(i0), l3e(i0)
        l1 = n3(1,mu3(j1))
        l2 = n3(2,mu3(j1))
        l3 = n3(3,mu3(j1))
        Do i1 = 1, vl
          If ( col(i1) == l1 .and. row(i1) == i0 ) ns31(j1) = i1
          If ( col(i1) == l2 .and. row(i1) == i0 ) ns32(j1) = i1
          If ( col(i1) == l3 .and. row(i1) == i0 ) ns33(j1) = i1
        EndDo
      EndDo
      Do j1 = l4a(i0), l4e(i0)
        l1 = n4(1,mu4(j1))
        l2 = n4(2,mu4(j1))
        l3 = n4(3,mu4(j1))
        l4 = n4(4,mu4(j1))
        Do i1 = 1, vl
          If ( col(i1) == l1 .and. row(i1) == i0 ) ns41(j1) = i1
          If ( col(i1) == l2 .and. row(i1) == i0 ) ns42(j1) = i1
          If ( col(i1) == l3 .and. row(i1) == i0 ) ns43(j1) = i1
          If ( col(i1) == l4 .and. row(i1) == i0 ) ns44(j1) = i1
        EndDo
      EndDo
    EndDo

    ! Output the CRS jacobian format data
    Open(newunit=lun_sparse, file=trim(data_dir)//'/sparse_ind', form='unformatted')
    Write(lun_sparse) vl
    Write(lun_sparse) row(1:vl),col(1:vl),pb
    Write(lun_sparse) (nan(i),i=1,4)
    Write(lun_sparse) ns11,ns21,ns22
    Write(lun_sparse) ns31
    Write(lun_sparse) ns32
    Write(lun_sparse) ns33
    Write(lun_sparse) ns41
    Write(lun_sparse) ns42
    Write(lun_sparse) ns43
    Write(lun_sparse) ns44
    Close(lun_sparse)

    Return
  End Subroutine sparse_check

  Subroutine net_preprocess(lun_out,data_dir,data_desc)
    !-----------------------------------------------------------------------------------------------
    ! This is the outer shell of net_setup. This routine also reformats the reaction data.
    !-----------------------------------------------------------------------------------------------
    Use, Intrinsic :: iso_fortran_env, Only: iostat_end
    Use nuclear_data, Only: ny, nname, aa, zz, nn, be, mex, iz, in, mex_n, mex_p, &
      & index_from_name, format_nuclear_data
    Use xnet_util, Only: ifactorial, xnet_terminate
    Implicit None

    ! Input variables
    Integer, Intent(in) :: lun_out    ! The logical unit for progress output
    Character(*), Intent(in) :: data_dir  ! The location of the data directory
    Character(*), Intent(in) :: data_desc ! Description of data directory

    ! Local variables
    Integer, Parameter :: ndesc = 169 ! Number of unique REACLIB rate labels
    Integer      :: n, i, j, k, l, jj ! Loop indicies
    Integer      :: nk1, nk2, nk3, nk4
    Integer      :: ni(nmax), ki(nmax), ierr
    Integer      :: nffn, nnnu
    Integer      :: krt, kl, lo
    Integer      :: iec, ier, iev    ! Reaction flags
    Integer      :: idesc, nnew_desc ! Reaction descriptor indices
    Integer      :: lun_blank, lun_su
    Real(dp)     :: p(7), q, qtemp, r
    Character(5) :: nucname(nmax), blank5 = '     '
    Character(4) :: desc, desc_new(100)
    Character(1) :: nr, vw, nrr = 'r', nrn = 'n', nrw = 'w', vvw = 'v', blank1 = ' '
    Character(4) :: desc_known(ndesc) = &
      & (/'  ec','bet+','betp','bet-','betm','btyk',' bec','bkmo', &
      &   'mo92',' ffn',' lmp','fkth','cf88','wies','laur','bb92', &
      &   'baka','rolf','wawo','mafo','ms01',' wag','bb93','ca85', &
      &   'ornl',' wfh','rath','rpsm','hg94','ro95','wien','iucf', &
      &   'smwa','smmo','wi2p','aexp','pexp','bex+','bex-','bqa+', &
      &   'bhi+','bec ','bklq','nuba','bede','bqa+','gamm','ja01', &
      &   'nosm','nacr','tang','ta04','ch04','ku02','il04','go00', &
      &   'most','il03','de04','da03','vanc','il01','ha00','ku92', &
      &   'ro04','sc83','he98','po00','il00','he95','sc98','vo00', &
      &   'il99','beau','da04','ca00','da77','bl03','ch01','fynb', &
      &   'De99','bu96','ba00','Ha96','sh03','tr04','sc05','nfis', &
      &   'GW95','ch05','re98','nb03','mo03','wc07','an06','ka02', &
      &   'gl07','GG90','lg06','im05','GG92','ct07','ga00','ww02', &
      &   'wh87','dh03','ha01','cb01','HF99','hg95','pg06','jm06', &
      &   'sb05','SM93','SM86','thra','fy05','ha04','pt05','bk91', &
      &   'bk92','kg91','SM91','iwam',' nen','nebp','lznu','lzan', &
      &   'ths8','wc12','nk06','mb11','jz10','cd08','nw00','il10', &
      &   'li10','chw0','co10','kd02','st08','pkrF','br10','hi09', &
      &   'ld11','hd08','ua08','ol11','ls09','ia08','dc11','mb07', &
      &   'wc17','mo97','ks03','nac2','mp17','cb09','li12','ma10', &
      &   'mm11','si13','sa12','ks12','hg12','rk12','ac12','gl12', &
      &   'mv09'/)

    nnew_desc = 1
    desc_new(1) = '    '

    Open(newunit=lun_su, file=trim(data_dir)//'/netsu')
    Open(newunit=lun_s3, file=trim(data_dir)//'/nets3', form='unformatted')
    Open(newunit=lun_ndiag, file=trim(data_dir)//'/net_diag')

    ! Read and Reformat the nuclear data
    Call format_nuclear_data(lun_out,data_dir)
    Write(lun_out,"('Listing Nuclei')")
    Write(lun_out,"(10a6)") (nname(i),i=1,ny)
    Write(lun_out,"(1x,'Number of species=',i4)") ny

    ! Read and reformat the reaction data
    Write(lun_out,"('Reading and reformating reactions data')")
    Read(lun_su,"(2i5)") nffn, nnnu
    la(:) = 1
    le(:) = 0
    nk1 = 1
    nk2 = 1
    nk3 = 1
    nk4 = 1
    n = 1
    krt = 1
    Do jj = 1, nrmax

      ! Read in each entry from the reaction data file
      Read(lun_su,"(i2,3x,6a5,8x,a4,a1,a1,3x,1pe12.5)",iostat=ierr) &
        & k, (nucname(j),j=1,nmax), desc, nr, vw, q
      If ( ierr == iostat_end ) Then
        Exit
      ElseIf ( ierr /= 0 ) Then
        Write(lun_out,*) 'Problem reading reaction ',jj
        Call xnet_terminate('Problem reading reaction',jj)
      EndIf
      Read(lun_su,"(4e13.6)",iostat=ierr) (p(j),j=1,7)
      If ( ierr == iostat_end ) Then
        Exit
      ElseIf ( ierr /= 0 ) Then
        Write(lun_out,*) 'Problem reading reaction ',jj
        Call xnet_terminate('Problem reading reaction',jj)
      EndIf

      ! If the entry is a reaction
      ni(:) = 0
      If ( k == 0 ) Then

        ! Convert nuclear names to indicies
        Do i = 1, nmax
          Call index_from_name(nucname(i),ni(i))
        EndDo

        If ( any( ni > ny ) ) Then
          Write(lun_out,*) 'Dropping reaction ',(nucname(i),i=1,nmax),'Some nuclei not included in netsu.'
          Cycle
        EndIf

        ! Identify the reaction source against the list of known identifiers
        Do i = 1, ndesc
          If ( desc == desc_known(i) ) Exit
        EndDo
        If ( i > ndesc ) Then

          ! The reaction source is unknown, treat it as a normal strong rate.
          idesc = i

          ! Check if this is the first time this new reaction source has showed up in this particular netsu.
          Do i = 1, nnew_desc
            If ( desc == desc_new(i) ) Exit
          EndDo

          ! If this is the first time, remember the name so the user can be informed of unknown codes at the end.
          If ( i > nnew_desc ) Then
            nnew_desc = nnew_desc + 1
            desc_new(i) = desc
            Write(lun_out,"(3a,i6,a)") 'Reaction Desc. ',desc,' not found for entry ',jj,', treating as strong rate.'
          EndIf
        Else
          idesc = i
        EndIf

        !---------------------------------------------------------------------------------------------
        ! The reaction rates are calculated via a number of templates, depending on the type of
        ! reaction and its source. In this section, flags are set for correct treatment of these
        ! different types of reactions.
        !---------------------------------------------------------------------------------------------
        ! For weak interactions;
        ! (iec=1,4), Reaclib style weak rates
        ! (iec=2,3), Fuller, Fowler, Neuman rates
        ! (iec=7,8), Neutrino capture rates
        Select Case (idesc)
        Case (1)              ! electron capture rates from FFN
          iec = 1
        Case (10,11)          ! FFN type tabulated rate
          If ( iz(ni(1)) > iz(ni(2)) .and. iz(ni(1)) /= 1 ) Then
            iec = 3
          Else
            iec = 2
          EndIf
        Case (2:9,38:45,92:94,130,153:154,157) ! beta decay and non-tabulated ec
          iec = 4
        Case (125,127)        ! electron neutrino capture
          iec = 7
        Case (126,128)        ! electron antineutrino capture
          iec = 8
        Case Default
          iec = 0             ! not a weak rate
        End Select

        ! Check resonant nature of rate, resonant rate (ier=2), nonresonant (ier=1),
        ! or none of the above (ier=0).  JINA REALCIB also includes a weak flag.
        If ( nr == blank1 ) Then
          ier = 0
        ElseIf ( nr == nrw ) Then ! weak reaction
          ier = 0
          If ( iec == 0 ) Then
            Write(lun_out,"(5a,i2)") 'Reaction with description ',desc,' has resonant flag ',nr,' but iec=',iec
          EndIf
        ElseIf ( nr == nrn ) Then
          ier = 1
        ElseIf ( nr == nrr ) Then
          ier = 2
        Else
          ier = 2
        EndIf

        ! Check for inverse rate (iev=1). Inverse rates must be multiplied with
        ! Temperature dependant ratio of partition functions. iev=0 is a forward rate.
        If ( vw == vvw ) Then
          iev = 1
        Else
          iev = 0
        EndIf

        ! Divide nuclei name list into reactants (ki=-1), products (ki=+1) and blanks (ki=0)
        Do i = 1, nmax
          If ( i > nreactant(krt) ) Then
            If ( nucname(i) == blank5 ) Then
              ki(i) = 0
            Else
              ki(i) = 1
            EndIf
          Else
            ki(i) = -1
          EndIf
        EndDo

        ! Don't trust the Q-values in the netsu file
        qtemp = 0.0
        Do i = 1, nmax
          If ( i > nreactant(krt) ) Then
            If ( nucname(i) == blank5 ) Then
              qtemp = qtemp
            Else
              qtemp = qtemp + be(ni(i)) - mex_n*nn(ni(i)) - mex_p*zz(ni(i))
            EndIf
          Else
            qtemp = qtemp - be(ni(i)) + mex_n*nn(ni(i)) + mex_p*zz(ni(i))
          EndIf
        EndDo
        If ( abs(qtemp-q) > abs(0.1*q) ) Then
          Write(lun_ndiag,"(a,7a6,a,es9.2,a,es9.2)") &
            & 'Inconsistent q-value for ',desc,(nucname(i),i=1,nmax),'  netsu: ',q,' netwinv: ',qtemp
        EndIf
        q = qtemp

        ! Count identical reactants
        kl = nreactant(krt)
        j = 1
        Do i = 2, kl
          If ( nucname(i-1) /= nucname(i) ) Then
            j = i
          Else
            ki(j) = ki(j) + ki(i)
            ki(i) = 0
          EndIf
        EndDo

        ! Count identical products
        kl = nreactant(krt) + 2
        j = kl - 1
        Do i = kl, nmax
          If ( nucname(i-1) /= nucname(i) ) Then
            j = i
          Else
            ki(j) = ki(j) + ki(i)
            ki(i) = 0
          EndIf
        EndDo

        ! Calculate double counting corrections and output reaction data to nets3
        If ( nreactant(krt) == 1 ) Then
          r = 1.0
          Do i = 1, n1m
            n1(i,n) = ni(i)
            an1(i,n) = real(ki(i),dp) * r
          EndDo
          q1(n)    = q
          iwk1(n)  = iec
          ires1(n) = ier
          irev1(n) = iev
          desc1(n) = desc
          rc1(:,n) = p
          nk1 = nk1 + 1
        Else
          lo = 0
          Do i = 1, nreactant(krt)
            If ( ki(i) < lo ) lo = ki(i)
          EndDo
          r = 1.0 / real(ifactorial(-lo),dp)
          If ( nreactant(krt) == 2 ) Then
            Do i = 1, n2m
              n2(i,n) = ni(i)
              an2(i,n) = real(ki(i),dp) * r
            EndDo
            q2(n)    = q
            iwk2(n)  = iec
            ires2(n) = ier
            irev2(n) = iev
            desc2(n) = desc
            rc2(:,n) = p
            nk2 = nk2 + 1
          ElseIf ( nreactant(krt) == 3 ) Then
            Do i = 1, n3m
              n3(i,n) = ni(i)
              an3(i,n) = real(ki(i),dp) * r
            EndDo
            q3(n)    = q
            iwk3(n)  = iec
            ires3(n) = ier
            irev3(n) = iev
            desc3(n) = desc
            rc3(:,n) = p
            nk3 = nk3 + 1
          ElseIf ( nreactant(krt) == 4 ) Then
            Do i = 1, n4m
              n4(i,n) = ni(i)
              an4(i,n) = real(ki(i),dp) * r
            EndDo
            q4(n)    = q
            iwk4(n)  = iec
            ires4(n) = ier
            irev4(n) = iev
            desc4(n) = desc
            rc4(:,n) = p
            nk4 = nk4 + 1
          EndIf
        EndIf
        n = n + 1

      ! If data entry is not a reaction, it is a type or sub-type marker
      ! For type marker, reset type and sub-type counters
      ElseIf ( k == 1 ) Then
        krt = k
        la(k) = 1
        n = 1
      ElseIf ( k == 4 .or. k == 8 .or. k == 10 ) Then
        krt = k
        la(k) = 1
        le(k-1) = n - 1
        n = 1
      ElseIf ( k == 11 ) Then
        krt = k
        la(k) = le(3) + 1
        le(k-1) = n - 1
        n = le(3) + 1

      ! For sub-type marker, restart sub-type counters
      Else
        krt = k
        la(k) = n
        le(k-1) = n - 1
      EndIf
      !Write(lun_ndiag,"(1x,i3,6a5,6i4,6i2,3i2,2e12.4,7e12.4)") n,(nname(ni(i)),i=1,nmax),ni(:),ki(:),iec,ier,iev,r,q,p(:)
    EndDo
    le(krt) = n - 1
    nreac(1) = nk1 - 1
    nreac(2) = nk2 - 1
    nreac(3) = nk3 - 1
    nreac(4) = nk4 - 1
    Close(lun_su)

    Do n = 1, nreac(1)
      Write(lun_s3) n,(n1(i,n),i=1,n1m),iwk1(n),ires1(n),irev1(n),(rc1(j,n),j=1,7),q1(n)
    EndDo
    Do n = 1, nreac(2)
      Write(lun_s3) n,(n2(i,n),i=1,n2m),iwk2(n),ires2(n),irev2(n),(rc2(j,n),j=1,7),q2(n)
    EndDo
    Do n = 1, nreac(3)
      Write(lun_s3) n,(n3(i,n),i=1,n3m),iwk3(n),ires3(n),irev3(n),(rc3(j,n),j=1,7),q3(n)
    EndDo
    Do n = 1, nreac(4)
      Write(lun_s3) n,(n4(i,n),i=1,n4m),iwk4(n),ires4(n),irev4(n),(rc4(j,n),j=1,7),q4(n)
    EndDo

    !-------------------------------------------------------------------------------------------------
    ! In order to build the Jacobian efficiently (row or column wise rather than scattered),
    ! information about which reactions affect each nucleus is necessary.
    !-------------------------------------------------------------------------------------------------
    Write(lun_out,"('Registering reactions to species')")

    Allocate (l1a(ny),l1e(ny),l2a(ny),l2e(ny),l3a(ny),l3e(ny),l4a(ny),l4e(ny))

    ! For each nucleus, register each single reactant reaction which involves it.
    n = 0
    Do i = 1, ny
      l1a(i) = n + 1
      Do j = 1, nchap
        If ( nreactant(j) == 1 ) Then
          Do k = la(j), le(j)
            Do l = 1, nnucr(j)
              If ( n1(l,k) == i .and. abs(an1(l,k)) > tiny(0.0) ) Then
                n = n + 1
                a1(n) = an1(l,k)
                mu1(n) = k
                Write(lun_s3) a1(n),mu1(n)
              EndIf
            EndDo
          EndDo
        EndIf
      EndDo
      l1e(i) = n
    EndDo

    ! For each nucleus, register each two reactant reaction which involves it.
    n = 0
    Do i = 1, ny
      l2a(i) = n + 1
      Do j = 1, nchap
        If ( nreactant(j) == 2 ) Then
          Do k = la(j), le(j)
            Do l = 1, nnucr(j)
              If ( n2(l,k) == i .and. abs(an2(l,k)) > tiny(0.0) ) Then
                n = n + 1
                a2(n) = an2(l,k)
                mu2(n) = k
                Write(lun_s3) a2(n),mu2(n)
              EndIf
            EndDo
          EndDo
        EndIf
      EndDo
      l2e(i) = n
    EndDo

    ! For each nucleus, register each three reactant reaction which involves it.
    n = 0
    Do i = 1, ny
      l3a(i) = n + 1
      Do j = 1, nchap
        If ( nreactant(j) == 3 ) Then
          Do k = la(j), le(j)
            Do l = 1, nnucr(j)
              If ( n3(l,k) == i .and. abs(an3(l,k)) > tiny(0.0) ) Then
                n = n + 1
                a3(n) = an3(l,k)
                mu3(n) = k
                Write(lun_s3) a3(n),mu3(n)
              EndIf
            EndDo
          EndDo
        EndIf
      EndDo
      l3e(i) = n
    EndDo

    ! For each nucleus, register each four reactant reaction which involves it.
    n = 0
    Do i = 1, ny
      l4a(i) = n + 1
      Do j = 1, nchap
        If ( nreactant(j) == 4 ) Then
          Do k = la(j), le(j)
            Do l = 1, nnucr(j)
              If ( n4(l,k) == i .and. abs(an4(l,k)) > tiny(0.0) ) Then
                n = n + 1
                a4(n) = an4(l,k)
                mu4(n) = k
                Write(lun_s3) a4(n),mu4(n)
              EndIf
            EndDo
          EndDo
        EndIf
      EndDo
      l4e(i) = n
    EndDo
    nan(1) = l1e(ny)
    nan(2) = l2e(ny)
    nan(3) = l3e(ny)
    nan(4) = l4e(ny)
    Close(lun_s3)

    ! Output the nucleus to reaction registration to nets4
    Open(newunit=lun_s4, file=trim(data_dir)//'/nets4', form='unformatted')
    Write(lun_s4) ny
    Write(lun_s4) (nname(i),i=1,ny)
    Write(lun_s4) nffn,nnnu
    Write(lun_s4) (nreac(k),k=1,4)
    Do i = 1, ny
      Write(lun_s4) i,l1a(i),l1e(i),l2a(i),l2e(i),l3a(i),l3e(i),l4a(i),l4e(i)
    EndDo
    Write(lun_out,"(10x,'number of reactions of different types')")
    Write(lun_out,"(3(10x,i8))") (j,la(j),le(j),j=1,nchap)
    Write(lun_out,"(10x,'necessary dimensions')")
    Write(lun_out,"(4(10x,i8))") (nreac(k),k=1,4)
    Write(lun_out,"(4(10x,i8))") (nan(k),k=1,4)
    Write(lun_out,"(3(10x,i8))") nffn, nnnu

    ! Prepare network description file
    Open(newunit=lun_desc, file=trim(data_dir)//'/net_desc')
    Write(lun_desc,"(a80)") data_desc
    Write(lun_desc,"('Number of Nuclear Species=',i5)") ny
    Write(lun_desc,"('Reaction Count for the 11 different types')")
    Write(lun_desc,"(3i10)") (j,la(j),le(j),j=1,nchap)
    Write(lun_desc,"('Necessary dimensions')")
    Write(lun_desc,"(a17,4i8)") 'nreac(1,2,3,4)=',(nreac(k),k=1,4)
    Write(lun_desc,"(a17,4i8)") 'nan(1,2,3,4)=',(nan(k),k=1,4)
    Write(lun_desc,"('Reaction Count for non-REACLIB rates')")
    Write(lun_desc,"(a7,i8,a7,i8)") 'nffn=',nffn,' nnu=',nnnu

    ! Match reactions
    Write(lun_out,*) "Matching forward and reverse reactions"
    Call match_react(data_dir)

    ! Calculate sparseness
    Write(lun_out,*) "Determining network sparsity pattern"
    Call sparse_check(data_dir)
    Close(lun_s4)

    ! Create inab template
    Open(newunit=lun_blank, file=trim(data_dir)//'/ab_blank')
    Write(lun_blank,"(a)") 'Abundance Description'
    Write(lun_blank,"(4(a5,es14.7,1x))") (nname(i),0.0,i=1,ny)
    Close(lun_blank)

    ! Tell user to add new descriptions to net_setup
    If ( nnew_desc > 1 ) Then
      Write(lun_out,*) "These descriptions need to be added to net_setup:"
      Write(lun_out,"(8a4)") (desc_new(i),i=2,nnew_desc)
    EndIf

    ! Deallocate Nuclear Data arrays
    Deallocate (nname,aa,zz,nn,be,mex,iz,in)

    ! Close other files
    Close(lun_desc)
    Close(lun_ndiag)

    Return
  End Subroutine net_preprocess

End Module xnet_preprocess
