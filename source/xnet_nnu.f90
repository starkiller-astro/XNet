!***************************************************************************************************
! xnet_nnu.f90 for XNet version 7 3/13/12
! This file contains reaction data structures for neutrino-induced reactions according to Froehlich
! PhD thesis 2007 (using rates from Zinner & Langanke) and routines to read the data and calculate
! the reaction rates.  It also contains the data structures for time history of the neutrino fluxes
! and temperatures.
!
! Credit: Carla Froehlich
!***************************************************************************************************

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

  Integer, Allocatable :: irl(:)          ! Reaclib Chapter 1 index of reaction
  Integer, Allocatable :: nuspec(:)       ! The neutrino species involved in the reaction

  Real(dp), Allocatable :: sigmanu(:,:) ! dim(nnnu,ntnu)
  Real(dp), Allocatable :: rnnu(:,:,:)

Contains

  Subroutine read_nnu_data(nnnu,data_dir)
    !-----------------------------------------------------------------------------------------------
    ! This routine reads the neutrino cross sections [in units of 10^-42 cm^2]
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Integer, Intent(in) :: nnnu
    Character(*), Intent(in) :: data_dir

    ! Local variables
    Integer :: i, j, lun_nnu

    Allocate (sigmanu(nnnu,ntnu))
    sigmanu = 0.0

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
    Use xnet_controls, Only: idiag, lun_diag
    Implicit None

    ! Input variables
    Integer, Intent(in) :: nnnu, nr, iwk(nr),innu(nr)

    ! Local variables
    Integer :: i, j

    Allocate(irl(nnnu),nuspec(nnnu))
    Do j = 1, nr
      If(innu(j)>0) Then
        If (idiag>5) Then
            write(lun_diag,*) "[NU]",j,iwk(j)
        Endif
        irl(innu(j)) = j
        If(iwk(j).eq.7) Then
          nuspec(innu(j))=1
        ElseIf(iwk(j).eq.8) Then
          nuspec(innu(j))=2
        ElseIf(iwk(j).eq.9 ) Then
          nuspec(innu(j))=3
        ElseIf(iwk(j).eq.10) Then
          nuspec(innu(j))=4
        Else
          Write(6,*) 'nnu match error',j,innu(j),iwk(j)
        EndIf
      EndIf
    EndDo

    If(idiag==6) Then
      Write(lun_diag,"(a)") 'i,IRL,Nuspec'
      Write(lun_diag,"(3i6)") (i,irl(i),nuspec(i),i=1,nnnu)
    EndIf

    Return

  End Subroutine nnu_match

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
    Use xnet_controls, Only: nzevolve, zb_lo, zb_hi, lzactive, ineutrino, idiag, lun_diag, szbatch
    Use xnet_types, Only: dp
    Use xnet_util, Only: safe_exp
    Implicit None

    ! Input variables
    Integer, Intent(in) :: nnnu
    Real(dp), Intent(in) :: time(nzevolve)

    ! Output variables
    Real(dp), Intent(out) :: rate(nnnu,nnuspec,zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Real(dp) :: ltnu, flux(nnuspec)
    Real(dp) :: lsigmanu1, lsigmanu2
    Real(dp) :: rcsnu(nnnu)
    Real(dp) :: rdt, rdltnu, ltnu1, ltnu2, lfl1,lfl2
    Real(dp) :: xrate
    Integer :: izb, it, i, j, k, n, izone
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    ! Only interpolate if neutrino reactions are on
    If ( ineutrino == 0 ) Then
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          Do j = 1, nnuspec
            Do k = 1, nnnu
              rate(k,j,izb) = 0.0
            EndDo
          EndDo
        EndIf
      EndDo
    Else
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then

          ! For constant conditions (nh = 1), do not interpolate in time
          If ( nh(izb) == 1 ) Then
            n = 1

          ! Otherwise, intepolate the neutrino fluxes and temperature from time history
          Else
            Do i = 1, nh(izb)
              If ( time(izb) <= th(i,izb) ) Exit
            EndDo
            n = i
          EndIf

          ! Compute neutrino cross sections
          Do j = 1, nnuspec

            ! Linear interpolation in log-space
            If ( n == 1 ) Then
              ltnu = log(tmevnu(1,j,izb))
              flux(j) = fluxcms(1,j,izb)
            ElseIf ( n > nh(izb) ) Then
              ltnu = log(tmevnu(nh(izb),j,izb))
              flux(j) = fluxcms(nh(izb),j,izb)
            Else
              rdt = (time(izb)-th(n-1,izb)) / (th(n,izb)-th(n-1,izb))
              ltnu = rdt*log(tmevnu(n,j,izb)) + (1.0-rdt)*log(tmevnu(n-1,j,izb))
              
              If ( fluxcms(n,j,izb) ==0. .and. fluxcms(n-1,j,izb)==0. ) then
                      flux(j) = 0.
              Else If ( fluxcms(n,j,izb) ==0. .or. fluxcms(n-1,j,izb)==0. ) then
                      flux(j) = rdt*fluxcms(n,j,izb) + (1.0-rdt)*fluxcms(n-1,j,izb)
              Else
                      flux (j) = safe_exp( rdt*log(fluxcms(n,j,izb)) + (1.0-rdt)*log(fluxcms(n-1,j,izb)) )
              Endif
          EndIf
            Do i = 1, ntnu
              If ( ltnu <= ltnugrid(i) ) Exit
            EndDo
            it = i

            Do k = 1, nnnu

              ! Log interpolation
              If ( it == 1 ) Then
                rcsnu(k) = sigmanu(k,1)
              ElseIf ( it > ntnu ) Then
                rcsnu(k) = sigmanu(k,ntnu)
              Else
                ltnu1 = ltnugrid(it-1)
                ltnu2 = ltnugrid(it)
                lsigmanu1 = log(sigmanu(k,it-1))
                lsigmanu2 = log(sigmanu(k,it))
                rdltnu = (ltnu-ltnu1) / (ltnu2-ltnu1)
                rcsnu(k) = safe_exp( rdltnu*lsigmanu2 + (1.0-rdltnu)*lsigmanu1 )
              EndIf

              ! Calculate rate only for neutrino specie involved in reaction
              If (nuspec(k) == j .or. (nuspec(k)==3 .and. j==1) .or. (nuspec(k)==4 .and. j==2) ) Then
                xrate = flux(j)*rcsnu(k)*1e-42
                If (xrate > 1.d-80) Then
                   rate(k,j,izb) = xrate
                Else
                   rate(k,j,izb) = 0.0d0
                Endif
               
!                write(*,*) "[NU rate]",k,j,nuspec(k),flux(j),rcsnu(k)
!                write(lun_diag,*) "[NU rate]",k,j,flux(j),rcsnu(k)
              Else
                rate(k,j,izb) = 0.0
              EndIf
            EndDo
          EndDo
        EndIf
      EndDo
    EndIf

    If ( idiag >= 6 ) Then
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          izone = izb + szbatch - zb_lo
          Write(lun_diag,"(a,i5)") 'NNU',izone
          Write(lun_diag,"(a,2i5)") 'Nuspec,Nnnu',nnuspec,nnnu
          Write(lun_diag,"(a,4es23.15)") 'Nu Flux',flux
          Write(lun_diag,"(i5,5es13.5)") (k, rcsnu(k), rate(k,:,izb), k=1,nnnu)
        EndIf
      EndDo
    EndIf

    Return
  End Subroutine nnu_rate

End Module xnet_nnu
