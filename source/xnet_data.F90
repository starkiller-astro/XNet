!***************************************************************************************************
! xnet_data.f90 10/18/17
! This file contains the nuclear and reaction data structures and the subroutines which read in the
! data and allocate the arrays.
!***************************************************************************************************

#include "xnet_macros.fh"

Module nuclear_data
  !-------------------------------------------------------------------------------------------------
  ! This module contains the essential data for each included species. Their array sizes are set in
  ! the routine read_nuclear_data.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None
  Integer                   :: ny                  ! The number of nuclear species evolved by the network
  Character(5), Allocatable :: nname(:)            ! Nuclei names (e.g. he4) (nname(0) is placeholder for non-nuclei)
  Real(dp), Allocatable     :: aa(:), zz(:), nn(:) ! Mass (aa), proton (zz), and neutron (nn) numbers
  Real(dp), Allocatable     :: be(:), mex(:)       ! Binding energy and mass excess [MeV c^{-2}]
  Real(dp), Allocatable     :: mm(:)               ! Mass of nuclei [g]
  Integer, Allocatable      :: ia(:), iz(:), in(:) ! Integer copies of mass, proton, and neutron numbers
  Integer                   :: inmin, inmax        ! Min and max neutron numbers
  Integer                   :: izmin, izmax        ! Min and max proton numbers

  ! Commonly used powers of Z (for, e.g., screening)
  Real(dp), Allocatable :: zz2(:), zz53(:), zzi(:)      ! zz^2, zz^{5/3}, and zz^{3b-1}
  Real(dp), Allocatable :: zseq(:), zseq53(:), zseqi(:) ! Sequence of numbers spanning the range of Z

  ! Get the neutron and proton mass excess (and masses) from netwinv for consistency
  Integer  :: ineut, iprot        ! Indices for neutron and proton in network
  Real(dp) :: mex_n ! = 8.0713179 ! Neutron mass excess [MeV c^{-2}]
  Real(dp) :: mex_p ! = 7.2889848 ! Proton mass excess [MeV c^{-2}]

  !-------------------------------------------------------------------------------------------------
  ! This module contains the nuclear partition function data. g is the interpolation data (ng,ny),
  ! gg is the current parition function, and angm is the J. gg(0) and angm(0) are placeholders for
  ! non-nuclei. The array sizes are set in read_nuclear_data.
  !-------------------------------------------------------------------------------------------------
  Integer, Parameter    :: ng = 24      ! Number of grid points for partition function data
  Integer, Allocatable  :: it9i(:)      ! Raw partition function temperature grid from file
  Real(dp), Allocatable :: t9i(:)       ! Temperature grid for partition function data
  Real(dp), Allocatable :: g(:,:)       ! Partition function data
  Real(dp), Allocatable :: gg(:,:)      ! Interpolated partition function
  Real(dp), Allocatable :: angm(:)      ! Angular momentum
  Real(dp), Allocatable :: dlngdt9(:,:) ! d(ln(partition functions))/dT9

  Interface benuc
    Module Procedure benuc_scalar
    Module Procedure benuc_vector
  End Interface benuc

Contains

  Subroutine index_from_name(nuc_name,inuc)
    !-----------------------------------------------------------------------------------------------
    ! This routine takes a nuclear name and finds the corresponding index for the current data set.
    ! inuc = 0      indicates a blank name.
    ! inuc = ny + 1 indicates that the nuclear name is not found in the current set.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: nuc_name

    ! Output variables
    Integer, Intent(out) :: inuc

    ! Local variables
    Integer :: i, name_len

    name_len = len_trim(nuc_name)
    If ( name_len > 0 ) then
      Do i = 1, ny
        If ( trim(adjustl(nuc_name)) == trim(adjustl(nname(i))) ) Exit
      EndDo
      inuc = i
    Else
      inuc = 0
    EndIf

    Return
  End Subroutine index_from_name

  Subroutine benuc_scalar(y,enb,enm)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the binding energy and mass excess energy [ergs/g] of the abundance
    ! distribution.
    !-----------------------------------------------------------------------------------------------
    Use xnet_constants, Only: epmev, avn
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: y(ny)

    ! Output variables
    Real(dp), Intent(out) :: enb ! Binding energy [ergs g^{-1}]
    Real(dp), Intent(out) :: enm ! Mass excess [ergs g^{-1}]

    ! Local variables
    Real(dp) :: ztot ! Total proton number
    Integer :: k

    ztot = 0.0
    enb  = 0.0
    Do k = 1, ny
      ztot = ztot + y(k) * zz(k)
      enb  = enb  + y(k) * be(k)
    EndDo
    enm  = mex_p*ztot + mex_n*(1.0-ztot) - enb

    ! Change units from MeV/nucleon to erg/g
    enb = epmev * avn * enb
    enm = epmev * avn * enm

    Return
  End Subroutine benuc_scalar

  Subroutine benuc_vector(y,enb,enm,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the binding energy and mass excess energy [ergs/g] of the abundance
    ! distribution.
    !-----------------------------------------------------------------------------------------------
    Use xnet_constants, Only: epmev, avn
    Use xnet_controls, Only: zb_lo, zb_hi, lzactive, tid
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: y(ny,zb_lo:zb_hi)

    ! Input/Output variables
    Real(dp), Intent(inout) :: enb(zb_lo:zb_hi) ! Binding energy [ergs g^{-1}]
    Real(dp), Intent(inout) :: enm(zb_lo:zb_hi) ! Mass excess [ergs g^{-1}]

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Real(dp) :: ztot, btot
    Integer :: k, izb
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    !XDIR XENTER_DATA XASYNC(tid) &
    !XDIR XCOPYIN(mask,y,enb,enm)

    !XDIR XLOOP_OUTER(1) XASYNC(tid) &
    !XDIR XPRESENT(mask,y,enb,enm,zz,be) &
    !XDIR XPRIVATE(ztot,btot)
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then

        ztot = 0.0
        btot = 0.0
        !XDIR XLOOP_INNER(1) &
        !XDIR XREDUCTION(+,ztot,btot)
        Do k = 1, ny
          ztot = ztot + y(k,izb) * zz(k)
          btot = btot + y(k,izb) * be(k)
        EndDo
        enb(izb) = btot
        enm(izb) = mex_p*ztot + mex_n*(1.0-ztot) - btot

        ! Change units from MeV/nucleon to erg/g
        enb(izb) = epmev * avn * enb(izb)
        enm(izb) = epmev * avn * enm(izb)
      EndIf
    EndDo

    !XDIR XEXIT_DATA XASYNC(tid) &
    !XDIR XCOPYOUT(enb,enm) &
    !XDIR XDELETE(mask,y)

    Return
  End Subroutine benuc_vector

  Subroutine partf(t9,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the nuclear partition functions as a function of temperature.
    !-----------------------------------------------------------------------------------------------
    Use xnet_controls, Only: idiag, iheat, lun_diag, zb_lo, zb_hi, lzactive, tid
    Use xnet_types, Only: dp
    Use xnet_util, Only: safe_exp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9(zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer :: i, ii, k, izb
    Real(dp) :: rdt9
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    !XDIR XENTER_DATA XASYNC(tid) &
    !XDIR XCOPYIN(mask,t9)

    !XDIR XLOOP_OUTER(1) XASYNC(tid) &
    !XDIR XPRESENT(mask,t9,t9i,gg,g,dlngdt9) &
    !XDIR XPRIVATE(ii,rdt9)
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then

        !XDIR XLOOP_SERIAL(1)
        Do i = 1, ng
          If ( t9(izb) <= t9i(i) ) Exit
        EndDo
        ii = i

        ! Linear interpolation in log-space
        gg(0,izb) = 1.0 ! placeholder for non-nuclei, gamma-rays, etc.
        Select Case (ii)
        Case (1)
          !XDIR XLOOP_INNER(1)
          Do k = 1, ny
            gg(k,izb) = g(1,k)
          EndDo
        Case (ng+1)
          !XDIR XLOOP_INNER(1)
          Do k = 1, ny
            gg(k,izb) = g(ng,k)
          EndDo
        Case Default
          rdt9 = (t9(izb)-t9i(ii-1)) / (t9i(ii)-t9i(ii-1))
          !XDIR XLOOP_INNER(1)
          Do k = 1, ny
            gg(k,izb) = safe_exp( rdt9*log(g(ii,k)) + (1.0-rdt9)*log(g(ii-1,k)) )
          EndDo
        End Select

        If ( iheat > 0 ) Then
          dlngdt9(0,izb) = 0.0
          Select Case (ii)
          Case (1)
            !XDIR XLOOP_INNER(1)
            Do k = 1, ny
              dlngdt9(k,izb) = log(g(2,k)/g(1,k)) / (t9i(2)-t9i(1))
            EndDo
          Case (ng+1)
            !XDIR XLOOP_INNER(1)
            Do k = 1, ny
              dlngdt9(k,izb) = log(g(ng,k)/g(ng-1,k)) / (t9i(ng)-t9i(ng-1))
            EndDo
          Case Default
            !XDIR XLOOP_INNER(1)
            Do k = 1, ny
              dlngdt9(k,izb) = log(g(ii,k)/g(ii-1,k)) / (t9i(ii)-t9i(ii-1))
            EndDo
          End Select
        EndIf
      EndIf
    EndDo

    !If ( idiag >= 1 ) Then
    !  Do izb = zb_lo, zb_hi
    !    If ( mask(izb) ) Then
    !      Write(lun_diag,"(a5,i3,es14.7)") 'PartF',ii,t9(izb)
    !      Write(lun_diag,"(5(i4,es12.4))") (i, gg(i,izb), i=1,ny)
    !    EndIf
    !  EndDo
    !EndIf

    !XDIR XEXIT_DATA XASYNC(tid) &
    !XDIR XDELETE(mask,t9)

    Return
  End Subroutine partf

  Subroutine read_sunet(data_dir)
    Use, Intrinsic :: iso_fortran_env, Only: iostat_end
    Use xnet_util, Only: xnet_terminate
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: data_dir

    ! Local variables
    Character(256) :: filename
    Integer :: lun_sunet, inuc, ierr

    filename = trim(data_dir)//'/sunet'
    Open(newunit=lun_sunet, file=trim(filename), status='old', action='read', iostat=ierr)
    If ( ierr /= 0 ) Call xnet_terminate('Failed to open sunet file',ierr)

    If ( .not. allocated(nname) ) Then
      inuc = 1
      Do
        Read (lun_sunet,*,iostat=ierr)
        If ( ierr == iostat_end ) Then
          Exit
        ElseIf ( ierr /= 0 ) Then
          Call xnet_terminate('Error reading sunet file',ierr)
        EndIf
        inuc = inuc + 1
      end do
      ny = inuc - 1
      Allocate (nname(0:ny))
      Rewind(lun_sunet)
    EndIf
    nname(0) = ' === '
    Do inuc = 1, ny
      Read(lun_sunet,"(a5)",iostat=ierr) nname(inuc)
      If ( ierr /= 0 ) Call xnet_terminate('Error reading sunet file',ierr)
    EndDo
    Close(lun_sunet)

    Return
  End Subroutine read_sunet

  Subroutine read_netwinv(data_dir)
    Use, Intrinsic :: iso_fortran_env, Only: iostat_end
    Use xnet_types, Only: dp
    Use xnet_util, Only: xnet_terminate
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: data_dir

    ! Local variables
    Character(256) :: filename
    Character(5) :: nam
    Real(dp) :: spin ! Ground state spin
    Integer :: inuc, j, ierr, lun_winv

    filename = trim(data_dir)//'/netwinv'
    Open(newunit=lun_winv, file=trim(filename), status='old', action='read', iostat=ierr)
    If ( ierr /= 0 ) Call xnet_terminate('Failed to open netwinv file',ierr)

    Read(lun_winv,"(i5)",iostat=ierr) ny
    If ( ierr /= 0 ) Call xnet_terminate('Error reading netwinv file',ierr)

    ! Read in the partition function iteration grid, and fix endpoints
    If ( .not. allocated(it9i) ) Allocate (it9i(ng))
    If ( .not. allocated(t9i) )  Allocate (t9i(ng))
    Read(lun_winv,"(24i3)",iostat=ierr) (it9i(j), j=1,ng)
    If ( ierr /= 0 ) Call xnet_terminate('Error reading netwinv file',ierr)

    ! Convert to the proper units
    t9i = real(it9i,dp)
    t9i(1:ng-1) = 0.01*t9i(1:ng-1)
    t9i(ng) = 0.1*t9i(ng)

    ! Read the nuclear names
    If ( .not. allocated(nname) ) Allocate (nname(0:ny))
    nname(0) = ' === '
    Do inuc = 1, ny
      Read(lun_winv,"(a5)",iostat=ierr) nname(inuc)
      If ( ierr /= 0 ) Call xnet_terminate('Error reading netwinv file',ierr)
    EndDo

    ! Set size of nuclear parameter arrays and read in nuclear parameters
    ! and partition function interpolation table.
    If ( .not. allocated(aa) )   Allocate (aa(ny))
    If ( .not. allocated(iz) )   Allocate (iz(ny))
    If ( .not. allocated(in) )   Allocate (in(ny))
    If ( .not. allocated(mex) )  Allocate (mex(ny))
    If ( .not. allocated(g) )    Allocate (g(ng,ny))
    If ( .not. allocated(angm) ) Allocate (angm(0:ny))
    angm(0) = 0.0
    Do inuc = 1, ny
      Read(lun_winv,*,iostat=ierr) nam, aa(inuc), iz(inuc), in(inuc), spin, mex(inuc)
      If ( ierr /= 0 ) Call xnet_terminate('Error reading netwinv file',ierr)
      Read(lun_winv,*,iostat=ierr) (g(j,inuc), j=1,ng)
      If ( ierr /= 0 ) Call xnet_terminate('Error reading netwinv file',ierr)
      angm(inuc) = 2.0*spin + 1.0

      ! Check that data entry name matches header
      If ( adjustl(nam) /= adjustl(nname(inuc)) ) Call xnet_terminate('netwinv data entry name does not match header for inuc=',inuc)
    EndDo
    Close(lun_winv)

    Return
  End Subroutine read_netwinv

  Subroutine read_nuclear_data(data_dir,data_desc)
    !-----------------------------------------------------------------------------------------------
    ! This routine reads, from the file netwinv, the nuclei included along with the nuclear data
    ! which will be needed for later calculations. This data includes the atomic number, the number
    ! of protons and neutrons, and the binding energy (calculated from the tabulated mass excess).
    ! Also the tabulations of the partition functions, g, are read in for later interpolation. Once
    ! the set of nuclear data is read in, it is assigned to the proper nuclei.
    !-----------------------------------------------------------------------------------------------
    Use xnet_constants, Only: avn, bip1, m_e, m_n, m_p, m_u, five3rd, thbim1
    Use xnet_controls, Only: iheat, nzevolve, tid
    Use xnet_parallel, Only: parallel_bcast, parallel_IOProcessor
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: data_dir ! Network data directory

    ! Output variables
    Character(80), Intent(out) :: data_desc ! Brief network description

    ! Local variables
    Real(dp) :: spin ! Ground state spin
    Integer :: lun_desc
    Integer :: i, j, ierr

    ! Read in the data description
    If ( parallel_IOProcessor() ) Then
      Open(newunit=lun_desc, file=trim(data_dir)//"/net_desc", status='old', action='read', iostat=ierr)
      If ( ierr == 0 ) Then
        Read(lun_desc,"(a80)") data_desc
        Close(lun_desc)
      Else
        data_desc = ""
      EndIf
    EndIf
    Call parallel_bcast(data_desc)

    ! Read the size of the network and partition function data
    If ( parallel_IOProcessor() ) Call read_sunet(data_dir)
    Call parallel_bcast(ny)

    ! Set size of nuclear data arrays and read in nuclear data and partition function interpolation table.
    ! nname(0), gg(0) and angm(0) are placeholders for non-nuclei.
    If ( .not. allocated(nname) ) Allocate (nname(0:ny))
    If ( .not. allocated(aa) )    Allocate (aa(ny))
    If ( .not. allocated(zz) )    Allocate (zz(ny))
    If ( .not. allocated(nn) )    Allocate (nn(ny))
    If ( .not. allocated(be) )    Allocate (be(ny))
    If ( .not. allocated(mex) )   Allocate (mex(ny))
    If ( .not. allocated(mm) )    Allocate (mm(ny))
    If ( .not. allocated(ia) )    Allocate (ia(ny))
    If ( .not. allocated(iz) )    Allocate (iz(ny))
    If ( .not. allocated(in) )    Allocate (in(ny))
    If ( .not. allocated(g) )     Allocate (g(ng,ny))
    If ( .not. allocated(angm) )  Allocate (angm(0:ny))
    If ( .not. allocated(it9i) )  Allocate (it9i(ng))
    If ( .not. allocated(t9i) )   Allocate (t9i(ng))
    If ( parallel_IOProcessor() ) Call read_netwinv(data_dir)
    Call parallel_bcast(it9i)
    Call parallel_bcast(t9i)
    Call parallel_bcast(nname)
    Call parallel_bcast(aa)
    Call parallel_bcast(iz)
    Call parallel_bcast(in)
    Call parallel_bcast(mex)
    Call parallel_bcast(angm)
    Call parallel_bcast(g)

    ia = nint(aa)
    zz = real(iz,dp)
    izmin = minval(iz)
    izmax = maxval(iz)
    nn = real(in,dp)
    inmin = minval(in)
    inmax = maxval(in)

    ! Some commonly used factors of Z
    Allocate (zz2(ny),zz53(ny),zzi(ny))
    zz2 = zz*zz
    zz53 = zz**five3rd
    zzi = zz**thbim1

    Allocate (zseq(0:izmax+2),zseq53(0:izmax+2),zseqi(0:izmax+2))
    zseq = (/ (real(i,dp), i=0,izmax+2) /)
    zseq53 = zseq**five3rd
    zseqi = zseq**bip1

    ! Get neutron and proton indices
    ineut = 0
    iprot = 0
    Do i = 1, ny
      If ( iz(i) == 0 .and. in(i) == 1 ) ineut = i
      If ( iz(i) == 1 .and. in(i) == 0 ) iprot = i
    EndDo

    ! For consistency, use neutron and proton mass excess from netwinv to
    ! calculate binding energies. Otherwise, use CODATA recommended 2014 values.
    If ( ineut > 0 ) Then
      mex_n = mex(ineut)
    Else
      mex_n = m_n - m_u
    EndIf
    If ( iprot > 0 ) Then
      mex_p = mex(iprot)
    Else
      mex_p = m_p + m_e - m_u
    EndIf
    be = mex_n*nn + mex_p*zz - mex

    ! Uncomment the commented end of the line below to use the actual mass instead of A*m_u
    mm = aa / avn! + mex(:)*epmev/(clt*clt)
    !mm(:) = zz(:)*(m_p+m_e) + nn(:)*m_n - be(:)*epmev/(clt*clt)

    Allocate (gg(0:ny,nzevolve),dlngdt9(0:ny,nzevolve))

    !XDIR XENTER_DATA XASYNC(tid) &
    !XDIR XCOPYIN(aa,zz,nn,be,mex,mm,ia,iz,in) &
    !XDIR XCOPYIN(zz2,zz53,zzi,zseq,zseq53,zseqi) &
    !XDIR XCOPYIN(it9i,t9i,g,angm) &
    !XDIR XCREATE(gg,dlngdt9)

    Return
  End Subroutine read_nuclear_data

  Subroutine format_nuclear_data(lun_out,data_dir)
    !-----------------------------------------------------------------------------------------------
    ! This routine reads from the file netwinv, the nuclei included along with the nuclear data
    ! which will be needed for later calculations. This data includes the atomic number, the number
    ! of protons and neutrons, and the binding energy (calculated from the tabulated mass excess).
    ! Also, the tabulations of the  partition functions, g, are read in for later interpolation. A
    ! binary data file is then prepared for faster input.
    !-----------------------------------------------------------------------------------------------
    Use, Intrinsic :: iso_fortran_env, Only: iostat_end
    Use xnet_constants, Only: m_e, m_n, m_p, m_u
    Use xnet_types, Only: dp
    Use xnet_util, Only: xnet_terminate
    Implicit None

    ! Input variables
    Integer, Intent(in) :: lun_out
    Character(*), Intent(in) :: data_dir

    ! Local variables
    Character(5), Allocatable :: nname_test(:)
    Integer :: ny_test
    Real(dp) :: spin ! Ground state spin
    Integer :: i, n, l, m, ntest, ierr
    Integer :: lun_data, lun_sunet, lun_winv

    ! Read in sunet
    Call read_sunet(data_dir)
    Write(lun_out,*) ny
    Allocate (nname_test(0:ny))
    nname_test = nname
    ny_test = ny

    ! Read nuclear data from netwinv and check for consistency with sunet
    Call read_netwinv(data_dir)
    If ( ny /= ny_test ) Then
      Write(lun_out,*) 'netwinv /= sunet',i,ny,ny_test
      Call xnet_terminate('netwinv /= sunet')
    EndIf
    Do n = 1, ny
      If ( adjustl(nname(n)) /= adjustl(nname_test(n)) ) Then
        Write(lun_out,*) 'netwinv /= sunet',n,nname(n),nname_test(n)
        Call xnet_terminate('netwinv /= sunet')
      EndIf
    EndDo
    If ( .not. allocated(zz) )   Allocate (zz(ny))
    If ( .not. allocated(nn) )   Allocate (nn(ny))
    If ( .not. allocated(be) )   Allocate (be(ny))
    zz = real(iz,dp)
    nn = real(in,dp)

    ! Get neutron and proton indices
    ineut = 0
    iprot = 0
    Do i = 1, ny
      If ( iz(i) == 0 .and. in(i) == 1 ) ineut = i
      If ( iz(i) == 1 .and. in(i) == 0 ) iprot = i
    EndDo

    ! For consistency, use neutron and proton mass excess from netwinv to
    ! calculate binding energies. Otherwise, use CODATA recommended 2014 values.
    If ( ineut > 0 ) Then
      mex_n = mex(ineut)
    Else
      mex_n = m_n - m_u
    EndIf
    If ( iprot > 0 ) Then
      mex_p = mex(iprot)
    Else
      mex_p = m_p + m_e - m_u
    EndIf
    be = mex_n*nn + mex_p*zz - mex

    ! Write binary data file
    Open(newunit=lun_data, file=trim(data_dir)//'/nuc_data', form='unformatted', action='write')
    Write(lun_data) ny
    Write(lun_data) t9i
    Write(lun_data) (nname(i),aa(i),zz(i),nn(i),be(i),(g(n,i),n=1,ng),angm(i),i=1,ny)
    Close(lun_data)

    ! Deallocate partition function data
    Deallocate (g,angm)

    Return
  End Subroutine format_nuclear_data

End Module nuclear_data

Module reaction_data
  !-------------------------------------------------------------------------------------------------
  ! This module contains the data needed to calculate the cross sections, reaction rates, and to map
  ! to each species those reactions which affect it.
  ! The csect variables are the results, velocity integrated cross section * density dependence.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None
  Integer               :: nreac(4)                                  ! # of reactions with i reactants
  Integer, Allocatable  :: n1i(:,:), n2i(:,:), n3i(:,:), n4i(:,:)    ! List of nuclei affected by each reaction
  Real(dp), Allocatable :: rc1(:,:), rc2(:,:), rc3(:,:), rc4(:,:)    ! REACLIB parameters
  Real(dp), Allocatable :: csect1(:,:), csect2(:,:), csect3(:,:), csect4(:,:)
  Real(dp), Allocatable :: dcsect1dt9(:,:), dcsect2dt9(:,:), dcsect3dt9(:,:), dcsect4dt9(:,:)

  ! Reaction flags to indicate variations in how the rate is calculated
  Integer, Allocatable :: iwk1(:), iwk2(:), iwk3(:), iwk4(:)     ! Weak Reaction
  Integer, Allocatable :: ires1(:), ires2(:), ires3(:), ires4(:) ! Resonant Reaction
  Integer, Allocatable :: irev1(:), irev2(:), irev3(:), irev4(:) ! Reverse Reaction

  ! Additional rate or energy factors
  Real(dp), Allocatable :: q1(:), q2(:), q3(:), q4(:)         ! Reaction Q values

  !-------------------------------------------------------------------------------------------------
  ! Aside from the REACLIB formated data, this dataset includes pointers to sets of external
  ! reaction rates, in the form of indices encoded in rc{1,2,3}(1).
  !-------------------------------------------------------------------------------------------------

  ! FFN are electron and positron capture rates encoded in the data format of Fuller, Fowler & Neuman (1982,1985).
  Integer              :: nffn    ! The number of FFN reactions
  Integer, Allocatable :: iffn(:) ! Pointers to FFN list

  ! NNU are neutrino and antineutrino capture rates from Zinner & Langanke, implemented by Carla Froehlich.
  Integer              :: nnnu    ! The number of NNU reactions
  Integer, Allocatable :: innu(:) ! Pointers to NNU list

  ! Extended reaction mapping arrays
  Integer               :: nan(4)                          ! Size of extended reaction->nuclei arrays
  Integer, Allocatable  :: la(:,:), le(:,:)                ! Extents in extended reaction arrays for each reactant
  Integer, Allocatable  :: mu1(:), mu2(:), mu3(:), mu4(:)  ! Index mapping rates to extended arrays

  ! nij(k) is the jth reactant in i-reactant reactions for reaction k in the extended arrays
  Integer, Allocatable  :: n11(:), n21(:), n22(:), n31(:), n32(:), n33(:), n41(:), n42(:), n43(:), n44(:)
  Integer, Allocatable  :: n10(:), n20(:), n30(:), n40(:)  ! Identifies i-reactant index for nij(k)

  ! Factors to indicate creation/destruction of species and avoid double-counting for identical reactants
  Real(dp), Allocatable :: a1(:), a2(:), a3(:), a4(:)

  ! Reaction rates after folding cross-sections in with counting factors
  Real(dp), Allocatable :: b1(:,:), b2(:,:), b3(:,:), b4(:,:) ! Coefficiencts of the Y terms in Eq. 10 of Hix & Meyer (2006)

Contains

  Subroutine read_reaction_data(data_dir)
    !-----------------------------------------------------------------------------------------------
    ! This routine reads in the necessary reaction data.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny, izmax, nname, zz
    Use xnet_constants, Only: five3rd
    Use xnet_controls, Only: iheat, iscrn, lun_stderr, nzevolve, iweak0, tid
    Use xnet_ffn, Only: read_ffn_data, ffnsum, ffnenu, ffn_ec, ffn_beta, ffn_qval, has_logft, &
      & qkffn, phasei, dphaseidt9, ngrid, rffn, dlnrffndt9
    Use xnet_nnu, Only: read_nnu_data, nnu_match, ntnu, nnuspec, sigmanu, ltnu, fluxnu, rnnu
    Use xnet_parallel, Only: parallel_bcast, parallel_IOProcessor
    Use xnet_types, Only: dp
    Use xnet_util, Only: xnet_terminate
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: data_dir

    ! Local variables
    Integer :: i, j, n, l, ierr
    Integer :: nr1, nr2, nr3, nr4
    Integer :: lun_s3, lun_s4

    ! Read in nuclear set, numbers of reactions, and extents of extended reaction arrays
    Allocate (la(4,ny),le(4,ny))
    If ( parallel_IOProcessor() ) Then
      Open(newunit=lun_s4, file=trim(data_dir)//"/nets4", form='unformatted', status='old', action='read', iostat=ierr)
      If ( ierr /= 0 ) Call xnet_terminate('Failed to open nets4 file',ierr)
      Read(lun_s4) ny
      Read(lun_s4) (nname(i), i=1,ny)
      Read(lun_s4) nffn, nnnu
      Read(lun_s4) (nreac(i), i=1,4)
      Do i = 1, ny
        Read(lun_s4) n, (la(j,i), le(j,i), j=1,4)
        If ( n /= i ) Then
          Write(lun_stderr,*) 'Error in nets4',i,n
          Call xnet_terminate('Error in nets4')
        EndIf
      EndDo
      Close(lun_s4)
    EndIf
    Call parallel_bcast(nffn)
    Call parallel_bcast(nnnu)
    Call parallel_bcast(nreac)
    Call parallel_bcast(la)
    Call parallel_bcast(le)

    ! If there are FFN rates, read in the FFN data and set FFN array sizes
    If ( nffn > 0) Then
      If ( parallel_IOProcessor() ) Then               
        Call read_ffn_data(nffn,data_dir)
      Else
        Allocate (ffnsum(nffn,ngrid),ffnenu(nffn,ngrid))
        Allocate (qkffn(nffn,nzevolve))
        Allocate (has_logft(nffn))
        ffnsum = 0.0
        ffnenu = 0.0
        qkffn = 0.0
        has_logft = 0

        ! Additional data for logft rates
        Allocate (ffn_ec(nffn,ngrid),ffn_beta(nffn,ngrid))
        Allocate (ffn_qval(nffn))
        Allocate (phasei(nffn,nzevolve),dphaseidt9(nffn,nzevolve))
        ffn_ec = 0.0
        ffn_beta = 0.0
        ffn_qval = 0.0
        phasei = 0.0
        dphaseidt9 = 0.0
      EndIf
      Call parallel_bcast(ffnsum)
      Call parallel_bcast(ffnenu)
      Call parallel_bcast(qkffn)
      Call parallel_bcast(has_logft)

      ! Additional data for logft rates
      Call parallel_bcast(ffn_ec)
      Call parallel_bcast(ffn_beta)
      Call parallel_bcast(ffn_qval)
      Call parallel_bcast(phasei)
      Call parallel_bcast(dphaseidt9)
    Else
      Allocate (ffnsum(1,ngrid),ffnenu(1,ngrid))
      Allocate (qkffn(1,nzevolve))
      Allocate (has_logft(1))
      ffnsum = 0.0
      ffnenu = 0.0
      qkffn = 0.0
      has_logft = 0

      ! Additional arrays for logft rates
      Allocate (ffn_ec(1,ngrid),ffn_beta(1,ngrid))
      Allocate (ffn_qval(1))
      Allocate (phasei(1,nzevolve),dphaseidt9(1,nzevolve))
      ffn_ec = 0.0
      ffn_beta = 0.0
      ffn_qval = 0.0
      phasei = 0.0
      dphaseidt9 = 0.0
    EndIf

    ! If there are NNU rates, read in the NNU data and set NNU array sizes
    If ( nnnu > 0 ) Then
      If ( parallel_IOProcessor() ) Then
        Call read_nnu_data(nnnu,data_dir)
      Else
        Allocate (ltnu(nnuspec,nzevolve))
        Allocate (fluxnu(nnuspec,nzevolve))
        Allocate (sigmanu(nnnu,ntnu))
        ltnu = 0.0
        fluxnu = 0.0
        sigmanu = 0.0
      EndIf
      Call parallel_bcast(sigmanu)
    Else
      Allocate (ltnu(nnuspec,nzevolve))
      Allocate (fluxnu(nnuspec,nzevolve))
      Allocate (sigmanu(1,ntnu))
      ltnu = 0.0
      fluxnu = 0.0
      sigmanu = 0.0
    EndIf

    ! Read in reaction arrays for 1 reactant reactions
    nr1 = nreac(1)
    Allocate (n1i(5,nr1),iwk1(nr1),ires1(nr1),irev1(nr1),rc1(7,nr1),q1(nr1))
    If ( parallel_IOProcessor() ) Then
      Open(newunit=lun_s3, file=trim(data_dir)//"/nets3", form='unformatted', status='old', action='read', iostat=ierr)
      If ( ierr /= 0 ) Call xnet_terminate('Failed to open nets3 file',ierr)
      Do j = 1, nr1
        Read(lun_s3) n, (n1i(l,j), l=1,5), iwk1(j), ires1(j), irev1(j), (rc1(l,j), l=1,7), q1(j)
        If ( n /= j ) Then
          Write(lun_stderr,*) 'Error in nets3, 1',j,n
          Call xnet_terminate('Error in nets3, 1')
        EndIf
      EndDo
    EndIf
    Call parallel_bcast(n1i)
    Call parallel_bcast(iwk1)
    Call parallel_bcast(ires1)
    Call parallel_bcast(irev1)
    Call parallel_bcast(rc1)
    Call parallel_bcast(q1)

    ! Read in reaction arrays for 2 reactant reactions
    nr2 = nreac(2)
    Allocate (n2i(6,nr2),iwk2(nr2),ires2(nr2),irev2(nr2),rc2(7,nr2),q2(nr2))
    If ( parallel_IOProcessor() ) Then
      Do j = 1, nr2
        Read(lun_s3) n, (n2i(l,j), l=1,6), iwk2(j), ires2(j), irev2(j), (rc2(l,j), l=1,7), q2(j)
        If ( n /= j ) Then
          Write(lun_stderr,*) 'Error in nets3, 2',j,n
          Call xnet_terminate('Error in nets3, 2')
        EndIf
      EndDo
    EndIf
    Call parallel_bcast(n2i)
    Call parallel_bcast(iwk2)
    Call parallel_bcast(ires2)
    Call parallel_bcast(irev2)
    Call parallel_bcast(rc2)
    Call parallel_bcast(q2)

    ! Read in reaction arrays for 3 reactant reactions
    nr3 = nreac(3)
    Allocate (n3i(6,nr3),iwk3(nr3),ires3(nr3),irev3(nr3),rc3(7,nr3),q3(nr3))
    If ( parallel_IOProcessor() ) Then
      Do j = 1, nr3
        Read(lun_s3) n, (n3i(l,j), l=1,6), iwk3(j), ires3(j), irev3(j), (rc3(l,j), l=1,7), q3(j)
        If ( n /= j ) Then
          Write(lun_stderr,*) 'Error in nets3, 3',j,n
          Call xnet_terminate('Error in nets3, 3')
        EndIf
      EndDo
    EndIf
    Call parallel_bcast(n3i)
    Call parallel_bcast(iwk3)
    Call parallel_bcast(ires3)
    Call parallel_bcast(irev3)
    Call parallel_bcast(rc3)
    Call parallel_bcast(q3)

    ! Read in reaction arrays for 4 reactant reactions
    nr4 = nreac(4)
    Allocate (n4i(6,nr4),iwk4(nr4),ires4(nr4),irev4(nr4),rc4(7,nr4),q4(nr4))
    If ( parallel_IOProcessor() ) Then
      Do j = 1, nr4
        Read(lun_s3) n, (n4i(l,j), l=1,6), iwk4(j), ires4(j), irev4(j), (rc4(l,j), l=1,7), q4(j)
        If ( n /= j ) Then
          Write(lun_stderr,*) 'Error in nets3, 4',j,n
          Call xnet_terminate('Error in nets3, 4')
        EndIf
      EndDo
    EndIf
    Call parallel_bcast(n4i)
    Call parallel_bcast(iwk4)
    Call parallel_bcast(ires4)
    Call parallel_bcast(irev4)
    Call parallel_bcast(rc4)
    Call parallel_bcast(q4)

    ! Calculate pointers to non-REACLIB data
    Allocate (iffn(nr1),innu(nr1))
    Where ( iwk1 == 2 .or. iwk1 == 3 ) ! FFN reaction
      iffn = nint(rc1(1,:))
      innu = 0
    ElseWhere ( iwk1 == 7 .or. iwk1 == 8 ) ! NNU reaction
      iffn = 0
      innu = nint(rc1(1,:))
    ElseWhere
      iffn = 0
      innu = 0
    EndWhere

    ! Calculate reverse pointers for non-REACLIB data
    Call nnu_match(nnnu,nr1,iwk1,innu)

    ! Allocate and read extended reaction->nuclei arrays linking nuclei to the reactions which affect them
    nan = le(:,ny)
    Allocate (mu1(nan(1)),a1(nan(1)),n10(nan(1)),n11(nan(1)))
    Allocate (mu2(nan(2)),a2(nan(2)),n20(nan(2)),n21(nan(2)),n22(nan(2)))
    Allocate (mu3(nan(3)),a3(nan(3)),n30(nan(3)),n31(nan(3)),n32(nan(3)),n33(nan(3)))
    Allocate (mu4(nan(4)),a4(nan(4)),n40(nan(4)),n41(nan(4)),n42(nan(4)),n43(nan(4)),n44(nan(4)))
    If ( parallel_IOProcessor() ) Then
      Do j = 1, nan(1)
        Read(lun_s3) a1(j), mu1(j)
      EndDo
      Do j = 1, nan(2)
        Read(lun_s3) a2(j), mu2(j)
      EndDo
      Do j = 1, nan(3)
        Read(lun_s3) a3(j), mu3(j)
      EndDo
      Do j = 1, nan(4)
        Read(lun_s3) a4(j), mu4(j)
      EndDo
      Close(lun_s3)
    EndIf
    Call parallel_bcast(a1)
    Call parallel_bcast(mu1)
    Call parallel_bcast(a2)
    Call parallel_bcast(mu2)
    Call parallel_bcast(a3)
    Call parallel_bcast(mu3)
    Call parallel_bcast(a4)
    Call parallel_bcast(mu4)

    Do i = 1, ny
      Do j = la(1,i), le(1,i)
        n10(j) = i
        n11(j) = n1i(1,mu1(j))
      EndDo
      Do j = la(2,i), le(2,i)
        n20(j) = i
        n21(j) = n2i(1,mu2(j))
        n22(j) = n2i(2,mu2(j))
      EndDo
      Do j = la(3,i), le(3,i)
        n30(j) = i
        n31(j) = n3i(1,mu3(j))
        n32(j) = n3i(2,mu3(j))
        n33(j) = n3i(3,mu3(j))
      EndDo
      Do j = la(4,i), le(4,i)
        n40(j) = i
        n41(j) = n4i(1,mu4(j))
        n42(j) = n4i(2,mu4(j))
        n43(j) = n4i(3,mu4(j))
        n44(j) = n4i(4,mu4(j))
      EndDo
    EndDo

    Allocate (b1(nan(1),nzevolve))
    Allocate (b2(nan(2),nzevolve))
    Allocate (b3(nan(3),nzevolve))
    Allocate (b4(nan(4),nzevolve))
    Allocate (csect1(nr1,nzevolve))
    Allocate (csect2(nr2,nzevolve))
    Allocate (csect3(nr3,nzevolve))
    Allocate (csect4(nr4,nzevolve))
    Allocate (dcsect1dt9(nr1,nzevolve))
    Allocate (dcsect2dt9(nr2,nzevolve))
    Allocate (dcsect3dt9(nr3,nzevolve))
    Allocate (dcsect4dt9(nr4,nzevolve))

    Allocate (rffn(max(1,nffn),nzevolve))
    Allocate (dlnrffndt9(max(1,nffn),nzevolve))
    Allocate (rnnu(max(1,nnnu),nnuspec,nzevolve))

    !XDIR XENTER_DATA XASYNC(tid) &
    !XDIR XCOPYIN(nreac,nan,la,le,iffn,innu,rc1,rc2,rc3,rc4) &
    !XDIR XCOPYIN(iwk1,iwk2,iwk3,iwk4,irev1,irev2,irev3,irev4) &
    !XDIR XCOPYIN(mu1,mu2,mu3,mu4,a1,a2,a3,a4,n1i,n2i,n3i,n4i) &
    !XDIR XCOPYIN(n10,n11,n20,n21,n22,n30,n31,n32,n33,n40,n41,n42,n43,n44) &
    !XDIR XCOPYIN(iffn,ffnsum,ffnenu,qkffn) &
    !XDIR XCOPYIN(has_logft,ffn_ec,ffn_beta,ffn_qval,phasei,dphaseidt9) &
    !XDIR XCOPYIN(innu,sigmanu,ltnu,fluxnu) &
    !XDIR XCREATE(b1,b2,b3,b4,csect1,csect2,csect3,csect4) &
    !XDIR XCREATE(dcsect1dt9,dcsect2dt9,dcsect3dt9,dcsect4dt9) &
    !XDIR XCREATE(rffn,dlnrffndt9,rnnu)

    Return
  End Subroutine read_reaction_data

  Subroutine enudot(y,sqnu,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the neutrino energy loss rate.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use xnet_constants, Only: epmev, avn
    Use xnet_controls, Only: zb_lo, zb_hi, lzactive, tid
    Use xnet_ffn, Only: qkffn
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: y(ny,zb_lo:zb_hi)

    ! Input/Output variables
    Real(dp), Intent(inout) :: sqnu(zb_lo:zb_hi) ! Neutrino energy loss rate [ergs g^{-1} s^{-1}]

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Real(dp) :: sqnu1
    Integer :: k, izb, nr1
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return
    
    nr1 = nreac(1)

    !XDIR XENTER_DATA XASYNC(tid) &
    !XDIR XCOPYIN(mask,y,sqnu)

    !XDIR XLOOP_OUTER(1) XASYNC(tid) &
    !XDIR XPRESENT(mask,qkffn,iffn,y,sqnu) &
    !XDIR XPRIVATE(sqnu1)
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then

        ! Loop over 1 species weak reactions to calculate neutrino loss
        sqnu1 = 0.0
        !XDIR XLOOP_INNER(1) &
        !XDIR XREDUCTION(+,sqnu1)
        Do k = 1, nr1
          sqnu1 = sqnu1 + y(n1i(1,k),izb) * qkffn(iffn(k),izb)
        EndDo

        ! Change units from MeV/nucleon to erg/g
        sqnu(izb) = epmev * avn * sqnu1
      EndIf
    EndDo

    !XDIR XEXIT_DATA XASYNC(tid) &
    !XDIR XCOPYOUT(sqnu) &
    !XDIR XDELETE(mask,y)

    Return
  End Subroutine enudot

End Module reaction_data
