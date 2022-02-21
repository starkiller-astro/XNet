!***************************************************************************************************
! xnet_ffn.f90 10/18/17
! This file contains reaction data structures for Fuller, Fowler, Neuman (FFN; 1982,1985) formated
! weak reactions and routines to read the data and calculate the reaction rates.
!***************************************************************************************************

Module xnet_ffn
  !-------------------------------------------------------------------------------------------------
  ! This module contains the data to calculate FFN formatted weak reactions.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None
  Integer, Parameter  :: nt9grid = 13             ! Number of temperature grid points per rate
  Integer, Parameter  :: nenegrid = 11            ! Number of electron density grid points per rate
  Integer, Parameter  :: ngrid = nt9grid*nenegrid ! Total number of grid points per rate
  Real(dp), Parameter :: t9grid(nt9grid) = &      ! Temperature grid for FFN rate data
    & (/ 0.01, 0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0 /)
  Real(dp), Parameter :: enegrid(nenegrid) = &    ! Electron density grid for FFN rate data
    & (/ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0 /)
  Real(dp), Allocatable :: ffn_ec(:,:), ffn_beta(:,:), ffnsum(:,:), ffnenu(:,:) ! dim(nffn,ngrid)
  Real(dp), Allocatable :: ffn_qval(:)
  Integer, Allocatable :: has_logft(:)

Contains

  Subroutine read_ffn_data(nffn,data_dir)
    !-----------------------------------------------------------------------------------------------
    ! This routine allocates and loads the data structures for FFN reaction rates.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Integer, Intent(in)  :: nffn
    Character(*), Intent(in) :: data_dir

    ! Local variables
    Integer :: i, j, lun_ffn
    Character(6) :: nuc1,nuc2
    Character(3) :: desc
    Real(dp) :: qval_in

    Allocate (ffnsum(nffn,ngrid),ffnenu(nffn,ngrid))
    Allocate (ffn_ec(nffn,ngrid),ffn_beta(nffn,ngrid))
    Allocate (ffn_qval(nffn))
    Allocate (has_logft(nffn))
    Open(newunit=lun_ffn, file=trim(data_dir)//"/netweak", status='old')
    Do i = 1, nffn
      Read(lun_ffn,*)  nuc1, nuc2,desc, ffn_qval(i)
      If (desc=='ft+') Then 
           has_logft(i)=2 !positron capture
      Else If (desc=='ft-') Then 
           has_logft(i)=1 !electron capture
      Else
           has_logft(i)=0
      EndIf
      Read(lun_ffn,*) (ffn_beta(i,j),ffn_ec(i,j),ffnsum(i,j), ffnenu(i,j), j=1,ngrid)
     
    ! Read logft, beta-decay rate and sum.
    EndDo
    Close(lun_ffn)

    Where (ffn_ec < -30d0) 
        ffn_ec = 99.d0
    End Where

    Return
  End Subroutine read_ffn_data

  Subroutine ffn_rate(nffn,t9,ene,rf,dlnrfdt9,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the reaction rates for FFN weak rates
    !-----------------------------------------------------------------------------------------------
    Use xnet_constants, Only: ln_10, bok, m_e
    Use xnet_controls, Only: iheat, nzevolve, zb_lo, zb_hi, lzactive
    Use xnet_conditions, Only: etae
    Use xnet_types, Only: dp
   
    Implicit None

    ! Input variables
    Integer, Intent(in)  :: nffn              ! Number of FFN rates
    Real(dp), Intent(in) :: t9(nzevolve)      ! Temperature [GK]
    Real(dp), Intent(in) :: ene(zb_lo:zb_hi)  ! Electron Density [g cm^{-3}]

    ! Output variables
    Real(dp), Intent(out) :: rf(nffn,zb_lo:zb_hi)       ! Temperature and density dependent FFN rates
    Real(dp), Intent(out) :: dlnrfdt9(nffn,zb_lo:zb_hi) ! Temperature and density dependent log FFN rate derivatives

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Real(dp), Parameter :: lrfmin = -30.0
    Real(dp) :: r1, r2, dr1, dr2
    Real(dp) :: rf_beta,rf_ec, phasei, cheme
    Real(dp) :: enel, dt9, dene, rdt9, rdene
    Integer :: i,izb, k, le1, lt1, i1, i2, i3, i4
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then

        ! Find the temperature grid point
        Do i = 1, nt9grid
          If ( t9(izb) <= t9grid(i) ) Exit
        EndDo
        lt1 = max(i-1,1)

        ! Find the density grid point (on log grid)
        enel = log10(ene(izb))
        le1 = max(int(enel),1)

        ! Bi-linear interpolation
        dt9 = t9grid(lt1+1) - t9grid(lt1)
        rdt9 = (t9(izb) - t9grid(lt1)) / dt9
        dene = enegrid(le1+1) - enegrid(le1)
        rdene = (enel - enegrid(le1)) / dene
        i1 = nt9grid*(le1-1) + lt1
        i2 = i1 + 1
        i3 = nt9grid*le1 + lt1
        i4 = i3 + 1

        cheme = etae(izb)*bok*t9(izb)+m_e
        Do k = 1, nffn

          If (has_logft(k) > 0) then

             If (has_logft(k)==1) Then
                Call effphase(t9(izb),cheme,ffn_qval(k),phasei)
             Else If (has_logft(k)==2) Then
                Call effphase(t9(izb),-cheme,ffn_qval(k),phasei)
             Else
                Write(*,*) "Problem in logft option. Should be 1 or 2, but is ",has_logft(k)
                Stop
             Endif 


             dr1 = ffn_ec(k,i2) - ffn_ec(k,i1)
             dr2 = ffn_ec(k,i4) - ffn_ec(k,i3)
             r1 = ffn_ec(k,i1) + rdt9*dr1
             r2 = ffn_ec(k,i3) + rdt9*dr2
             rf_ec = r1 + rdene*(r2 - r1) ! logft


             dr1 = ffn_beta(k,i2) - ffn_beta(k,i1)
             dr2 = ffn_beta(k,i4) - ffn_beta(k,i3)
             r1 = ffn_beta(k,i1) + rdt9*dr1
             r2 = ffn_beta(k,i3) + rdt9*dr2
             rf_beta = r1 + rdene*(r2 - r1)
           !  write(*,*) "rf_beta: ",rf_beta
             If (rf_ec < -30.0d0) then
                     rf_ec=0.
             Else
                     rf_ec = 10.0d0**rf_ec
             EndIf

             If (phasei< 1d-100 .or. rf_ec<1.d-30 ) Then
                     rf_ec = 0.0
             Else
                     rf_ec=dlog(2.0d0)*phasei/rf_ec ! turn into rate
             EndIf

             If (rf_beta > lrfmin) Then
                     rf_beta = 10.0**rf_beta
             Else
                     rf_beta=0.0
             EndIf
             rf(k,izb) = rf_beta + rf_ec
             
             !write(*,'(A6,I5,8ES14.3)') "FFN: ",k,ffn_qval(k),t9(izb),ene(izb),etae(izb),cheme,phasei,rf_beta,rf_ec

             If ( rf_beta < lrfmin .and. rf_ec < lrfmin ) Then
               rf(k,izb) = 0.0
               dlnrfdt9(k,izb) = 0.0
             Else
               ! Ignores temperature dependence of the logft based EC rate at the moment
               dlnrfdt9(k,izb) = ln_10 * ( rdene*dr2 + (1.0-rdene)*dr1 ) / dt9
             EndIf
          Else
             dr1 = ffnsum(k,i2) - ffnsum(k,i1)
             dr2 = ffnsum(k,i4) - ffnsum(k,i3)
             r1 = ffnsum(k,i1) + rdt9*dr1
             r2 = ffnsum(k,i3) + rdt9*dr2
             rf(k,izb) = r1 + rdene*(r2 - r1)
             If ( rf(k,izb) < lrfmin ) Then
               rf(k,izb) = 0.0
               dlnrfdt9(k,izb) = 0.0
             Else
               rf(k,izb) = 10.0**rf(k,izb)
               dlnrfdt9(k,izb) = ln_10 * ( rdene*dr2 + (1.0-rdene)*dr1 ) / dt9
             EndIf

          EndIf
        EndDo
      EndIf
    EndDo

    Return
  End Subroutine ffn_rate

  Subroutine effphase(t9,chem,qec,phasei)
!---------------------------------------------------------------------------
!  Calculates phase space factor for electron capture rates from tab. logft 
!  intent in:
!  temp....temperature in T9
!  chem....electron chemical potential in MeV (including electron rest mass)
!  qec.....Q-value for electron capture (MeV)
!  intent out:
!  phasei..phase space integral
!---------------------------------------------------------------------------
       use xnet_constants, only: m_e, bok
       use xnet_integ_phase, only: t9me,xchem,qcap
       Use xnet_types, Only: dp
       Implicit none
       Real (dp), Intent(in) :: t9,chem,qec
       Real (dp), Intent(out) :: phasei
!---------------------------------------------------------------------------
       Integer :: i
       Real (dp), Parameter :: phase_min=1.d-20
       Real (dp) :: tmev, emin, emax,comparison, conv_check
       Real (dp), external :: xeffc
!---------------------------------------------------------------------------
! Variables for QAGI integration routine
       Real (dp), parameter :: epsabs = 0.0d0, epsrel = 1.0d-08!1.0d-10
       integer, parameter :: inf = 1
       integer, parameter :: lw=1000, liw=lw/4
       Real (dp) ::  abserr
       Real (dp) :: w(lw)
       integer :: iw(liw)
       integer :: ifail,last,neval
!---------------------------------------------------------------------------
       qcap=qec/m_e
       emin=max(-qcap,1.0d0)
       xchem=chem/m_e
       t9me=t9*bok/m_e

       call dqagi(xeffc,emin,INF,epsabs,epsrel,phasei,abserr,neval, &
                 &     ifail ,liw,lw,last,iw,w)
      ! write(*,"(A6,5ES18.4)") "QAGI",emin,qcap,t9me,xchem,phasei
End subroutine effphase

End Module xnet_ffn
