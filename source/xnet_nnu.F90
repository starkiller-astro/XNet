!***************************************************************************************************
! xnet_nnu.f90 for XNet version 7 3/13/12
! This file contains reaction data structures for neutrino-induced reactions according to Froehlich
! PhD thesis 2007 (using rates from Zinner & Langanke) and routines to read the data and calculate
! the reaction rates.  It also contains the data structures for time history of the neutrino fluxes
! and temperatures.
!
! Credit: Carla Froehlich
!***************************************************************************************************

#include "xnet_macros.fh"

Module xnet_nnu
  !-------------------------------------------------------------------------------------------------
  ! This module contains the data and routines to calculate neutrino-induced reactions.
  !-------------------------------------------------------------------------------------------------
  Use xnet_conditions, Only: nhmx
  Use xnet_types, Only: dp
  Implicit None

  Integer, Parameter :: nnuspec = 4            ! Number of neutrino species in thermo history
  Integer, Parameter :: ntnu = 7               ! Number of neutrino temperature grid points per rate
  Real(dp), Parameter :: tnugrid(ntnu) = &     ! Neutrino temperature grid for NNU rate data
    & (/ 2.8, 3.5, 4.0, 5.0, 6.4, 8.0, 10.0 /)
  Real(dp), Parameter :: ltnugrid(ntnu) = log(tnugrid(:))

  Real(dp), Parameter :: sigmamin = 1.0e-100   ! Floor for neutrino cross-section rate

  Real(dp), Allocatable :: tmevnu(:,:,:)  ! Neutrino temperature [MeV]
  Real(dp), Allocatable :: fluxcms(:,:,:) ! Neutrino fluxes [cm^-2 s^-1]
  Real(dp), Allocatable :: ltnu(:,:)      ! Interpolated neutrino temperature [MeV]
  Real(dp), Allocatable :: fluxnu(:,:)    ! Interpolated neutrino fluxes [cm^-2 s^-1]

  Integer, Allocatable :: irl(:)          ! Reaclib Chapter 1 index of reaction
  Integer, Allocatable :: nuspec(:)       ! The neutrino species involved in the reaction

  Real(dp), Allocatable :: sigmanu(:,:) ! dim(nnnu,ntnu)
  Real(dp), Allocatable :: rnnu(:,:,:)

Contains

  Subroutine read_nnu_data(nnnu,data_dir)
    !-----------------------------------------------------------------------------------------------
    ! This routine reads the neutrino cross sections [in units of 10^-42 cm^2]
    !-----------------------------------------------------------------------------------------------
    Use xnet_controls, Only: nzevolve
    Implicit None

    ! Input variables
    Integer, Intent(in) :: nnnu
    Character(*), Intent(in) :: data_dir

    ! Local variables
    Integer :: i, j, lun_nnu

    Allocate (sigmanu(nnnu,ntnu))
    Allocate (ltnu(nnuspec,nzevolve))
    Allocate (fluxnu(nnuspec,nzevolve))
    sigmanu = 0.0
    ltnu = 0.0
    fluxnu = 0.0

    Open(newunit=lun_nnu, file=trim(data_dir)//"/netneutr", status='old', action='read')
    Do i = 1, nnnu
      Read(lun_nnu,*)
      Read(lun_nnu,*) (sigmanu(i,j), j=1,ntnu)
    EndDo
    Close(lun_nnu)
    sigmanu = max( sigmamin, sigmanu )

    Return
  End Subroutine read_nnu_data

  Subroutine nnu_match(nnnu,nr,iwk,innu)
    !-----------------------------------------------------------------------------------------------
    ! This routine finds placement in the REACLIB list of each neutrino reaction.
    ! From this, the neutrino species involved is determined.
    !-----------------------------------------------------------------------------------------------
    Use, Intrinsic :: iso_fortran_env, Only: lun_stdout=>output_unit
    Use xnet_controls, Only: idiag, lun_diag
    Implicit None

    ! Input variables
    Integer, Intent(in) :: nnnu, nr, iwk(nr), innu(nr)

    ! Local variables
    Integer :: i, j

    Allocate (irl(nnnu),nuspec(nnnu))
    irl = 0
    nuspec = 0

    Do j = 1, nr
      If ( innu(j) > 0 ) Then
        If ( idiag >= 6 ) Then
          Write(lun_diag,"(a,2i5)") "[NU]",j,iwk(j)
        EndIf
        irl(innu(j)) = j
        If ( iwk(j) == 7 ) Then
          nuspec(innu(j)) = 1
        ElseIf (iwk(j) == 8 ) Then
          nuspec(innu(j)) = 2
        ElseIf (iwk(j) == 9 ) Then
          nuspec(innu(j)) = 3
        ElseIf (iwk(j) == 10 ) Then
          nuspec(innu(j)) = 4
        Else
          Write(lun_stdout,"(a,3i5)") 'NNU match error',j,innu(j),iwk(j)
        EndIf
      EndIf
    EndDo

    If ( idiag>= 6 ) Then
      Write(lun_diag,"(a)") 'i,IRL,Nuspec'
      Write(lun_diag,"(3i6)") (i,irl(i),nuspec(i),i=1,nnnu)
    EndIf

    !__dir_enter_data &
    !__dir_async &
    !__dir_copyin(irl,nuspec)

    Return
  End Subroutine nnu_match

  Subroutine nnu_flux(tf,nf,ltnuf,fluxf,ts,ns,tnus,fluxs)
    !__dir_routine_seq
    Use xnet_types, Only: dp
    Use xnet_util, Only: safe_exp
    Implicit None

    ! Input variables
    Integer, Intent(in) :: ns
    Real(dp), Intent(in) :: tf, ts(:), tnus(:,:), fluxs(:,:)

    ! Output variables
    Integer, Intent(out) :: nf
    Real(dp), Intent(out) :: ltnuf(:), fluxf(:)

    ! Local variables
    Real(dp) :: dt, rdt
    Integer :: j, n

    ! For constant conditions (ns = 1), do not interpolate in time
    If ( ns == 1 ) Then
      nf = 1

    ! Otherwise, intepolate the neutrino fluxes and temperature from time history
    Else
      Do n = 1, ns
        If ( tf <= ts(n) ) Exit
      EndDo
      nf = n
      rdt = ( tf - ts(nf-1)) / ( ts(nf) - ts(nf-1) )
    EndIf

    Do j = 1, nnuspec
      If ( nf == 1 ) Then
        ltnuf(j) = log(tnus(1,j))
        fluxf(j) = fluxs(1,j)
      ElseIf ( nf > 1 .and. nf <= ns ) Then
        ltnuf(j) = rdt*log(tnus(nf,j)) + (1.0-rdt)*log(tnus(nf-1,j))
        If ( fluxs(nf,j) == 0.0 .and. fluxs(nf-1,j) == 0.0 ) Then
          fluxf(j) = 0.0
        ElseIf ( fluxs(nf,j) == 0.0 .or. fluxs(nf-1,j) == 0.0 ) Then
          fluxf(j) = rdt*fluxs(nf,j) + (1.0-rdt)*fluxs(nf-1,j)
        Else
          fluxf(j) = safe_exp( rdt*log(fluxs(nf,j)) + (1.0-rdt)*log(fluxs(nf-1,j)) )
        EndIf
      Else
        ltnuf(j) = log(tnus(ns,j))
        fluxf(j) = fluxs(ns,j)
      EndIf
    EndDo

    Return
  End Subroutine nnu_flux

  Subroutine nnu_rate(nnnu,time,rate,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the neutrino nucleus reaction rates [in s^-1] as 
    ! flux [cm^-2 s^-1] * cross section [cm^2].
    ! This flux already includes a factor 10^-42 from the cross sections.
    !
    ! Depending on the output of the prior simulation, the "neutrino flux" may need to be computed
    ! from neutrino luminosities.
    !-----------------------------------------------------------------------------------------------
    Use xnet_conditions, Only: nh, th
    Use xnet_controls, Only: zb_lo, zb_hi, lzactive, ineutrino, idiag, lun_diag, szbatch
    Use xnet_types, Only: dp
    Use xnet_util, Only: safe_exp
    Implicit None

    ! Input variables
    Integer, Intent(in) :: nnnu
    Real(dp), Intent(in) :: time(zb_lo:zb_hi)

    ! Output variables
    Real(dp), Intent(out) :: rate(nnnu,nnuspec,zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Real(dp) :: lsigmanu1, lsigmanu2
    Real(dp) :: rcsnu
    Real(dp) :: rdt, rdltnu, ltnu1, ltnu2, lfl1, lfl2
    Real(dp) :: xrate
    Integer :: izb, it, i, j, k, n, izone
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    !__dir_enter_data &
    !__dir_async &
    !__dir_create(rate) &
    !__dir_copyin(mask,time)

    ! Only interpolate if neutrino reactions are on
    If ( ineutrino == 0 ) Then
      !__dir_loop(3) &
      !__dir_async &
      !__dir_present(mask,rate)
      Do izb = zb_lo, zb_hi
        Do j = 1, nnuspec
          Do k = 1, nnnu
            If ( mask(izb) ) Then
              rate(k,j,izb) = 0.0
            EndIf
          EndDo
        EndDo
      EndDo
    Else

      ! Interpolate flux and neutrino temperature from time history
      !__dir_loop_outer(1) &
      !__dir_async &
      !__dir_present(mask,time,ltnu,fluxnu,th,nh,tmevnu,fluxcms)
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          Call nnu_flux(time(izb),n,ltnu(:,izb),fluxnu(:,izb), &
            th(:,izb),nh(izb),tmevnu(:,:,izb),fluxcms(:,:,izb))
        EndIf
      EndDo

      ! Compute neutrino cross sections
      !__dir_loop_outer(2) &
      !__dir_async &
      !__dir_present(mask,sigmanu,ltnu,fluxnu,nuspec) &
      !__dir_private(it)
      Do izb = zb_lo, zb_hi
        Do j = 1, nnuspec
          If ( mask(izb) ) Then

            Do i = 1, ntnu
              If ( ltnu(j,izb) <= ltnugrid(i) ) Exit
            EndDo
            it = i

            !__dir_loop_inner(1) &
            !__dir_private(ltnu1,ltnu2,lsigmanu1,lsigmanu2,rdltnu,rcsnu,xrate)
            Do k = 1, nnnu

              ! Log interpolation
              If ( it == 1 ) Then
                rcsnu = sigmanu(k,1)
              ElseIf ( it > ntnu ) Then
                rcsnu = sigmanu(k,ntnu)
              Else
                ltnu1 = ltnugrid(it-1)
                ltnu2 = ltnugrid(it)
                lsigmanu1 = log(sigmanu(k,it-1))
                lsigmanu2 = log(sigmanu(k,it))
                rdltnu = (ltnu(j,izb)-ltnu1) / (ltnu2-ltnu1)
                rcsnu = safe_exp( rdltnu*lsigmanu2 + (1.0-rdltnu)*lsigmanu1 )
              EndIf

              ! Calculate rate only for neutrino species involved in reaction
              If (  nuspec(k) == j .or. &
                & ( nuspec(k) == 3 .and. j == 1 ) .or. &
                & ( nuspec(k) == 4 .and. j == 2 ) ) Then
                xrate = fluxnu(j,izb)*rcsnu*1e-42 ! Not a fan of 1e-42 being hard-coded here
                If ( xrate > 1.0e-80 ) Then
                  rate(k,j,izb) = xrate
                Else
                  rate(k,j,izb) = 0.0
                Endif
              Else
                rate(k,j,izb) = 0.0
              EndIf
            EndDo
          EndIf
        EndDo
      EndDo
    EndIf

    If ( idiag >= 6 ) Then
      !__dir_update &
      !__dir_wait &
      !__dir_host(fluxnu,rate)
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          izone = izb + szbatch - zb_lo
          Write(lun_diag,"(a,i5)") 'NNU',izone
          Write(lun_diag,"(a,2i5)") 'Nuspec,Nnnu',nnuspec,nnnu
          Write(lun_diag,"(a,4es23.15)") 'Nu Flux',(fluxnu(j,izb), j=1,nnuspec)
          Write(lun_diag,"(i5,5es13.5)") (k, rcsnu, (rate(k,j,izb),j=1,nnuspec), k=1,nnnu)
        EndIf
      EndDo
    EndIf

    !__dir_exit_data &
    !__dir_async &
    !__dir_copyout(rate) &
    !__dir_delete(mask,time)

    !__dir_wait

    Return
  End Subroutine nnu_rate

End Module xnet_nnu
