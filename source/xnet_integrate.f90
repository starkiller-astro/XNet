!***************************************************************************************************
! xnet_integrate.f90 10/18/17
! This file contains subroutines shared by the various integrators.
!***************************************************************************************************

Module xnet_integrate
  Implicit None

Contains

  Subroutine timestep(kstep,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the trial timestep.
    ! For tdel >0, this calculation is based on the relative changes of the abundances during the
    !              previous step in the evolution and the integrators estimate of the next timestep.
    !          =0, the timestep is calculated using the time derivatives of the abundances based on
    !              reaction rates.
    !          <0, the timestep is held constant, tdel=-tdel_old.
    ! There is also the provision for limiting the timestep in the event that the thermodynamic
    ! conditions are changing too rapidly.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny, nname
    Use xnet_abundances, Only: y, yo, yt, ydot
    Use xnet_conditions, Only: t, tt, tdel, tdel_next, tdel_old, t9, t9o, t9t, rho, &
      & rhot, t9dot, cv, nt, ntt, ints, intso, tstop, tdelstart, nh, th, t9h, rhoh, t9rhofind
    Use xnet_controls, Only: changemx, changemxt, idiag, iheat, iscrn, iweak, lun_diag, yacc, &
      & ymin, szbatch, zb_lo, zb_hi, lzactive
    Use xnet_types, Only: dp
    Use xnet_util, Only: xnet_terminate
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Real(dp), Parameter :: changeth = 0.1
    Real(dp) :: changest
    Real(dp) :: rtau_y(ny), changey, tdel_dy(zb_lo:zb_hi)
    Real(dp) :: rtau_t9, changet9, tdel_dt9(zb_lo:zb_hi)
    Real(dp) :: tdel_stop(zb_lo:zb_hi), tdel_init
    Real(dp) :: dtherm(zb_lo:zb_hi)
    Integer :: i, j, k, izb, izone
    Logical :: mask_init(zb_lo:zb_hi)
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    ! Retain old values of timestep and thermo and calculate remaining time
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        tdel_old(izb) = tdel(izb)
        mask_init(izb) = ( abs(tdel_old(izb)) < tiny(0.0) )
      Else
        mask_init(izb) = .false.
      EndIf
    EndDo
    Call cross_sect(mask_in = mask_init)
    Call yderiv(mask_in = mask_init)

    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        tdel_stop(izb) = tstop(izb) - t(izb)
        tdel_dy(izb) = 0.0
        tdel_dt9(izb) = 0.0
        If ( mask_init(izb) ) Then
          intso(izb) = 0
          If ( tdelstart(izb) < 0.0 ) Then
            changest = 0.5
          Else
            changest = 0.01
          EndIf

          ! For unevolved initial abundances, as often found in test problems, derivatives may
          ! produce large jumps from zero abundance in the first timestep, though generally 
          ! inconsequential; an additional factor, changest, is used to increase accuracy.
          changey = changest*min(0.1,changemx)
          changet9 = changest*min(0.01,changemxt)
          tdel_init = min(changest*abs(tdelstart(izb)),tdel_stop(izb))
          If ( nh(izb) > 2 ) Then
            tdel_next(izb) = 1.0e-4*tdel_stop(izb)
          Else
            tdel_next(izb) = tdel_stop(izb)
          EndIf
        Else
          changey = changemx
          changet9 = changemxt
          tdel_init = 0.0
        EndIf

        ! Estimate timestep from relevant timescales
        If ( tdel_old(izb) >= 0.0 ) Then
          ! Calculate timescales for abundance changes, Y/(dY/dt)
          Do k = 1, ny
            If ( y(k,izb) > ymin ) Then
              If ( mask_init(izb) ) Then
                ! Calculate timescale directly from derivatives.
                rtau_y(k) = abs(ydot(k,izb)/max(y(k,izb),yacc))
              Else
                ! Calculate timescale from changes in last timestep.
                rtau_y(k) = abs((y(k,izb)-yo(k,izb))/max(y(k,izb),yacc) / tdel_old(izb))
              EndIf
            Else
              rtau_y(k) = 0.0
            EndIf
          EndDo

          ! If derivatives are non-zero, use Y/(dY/dt).
          ints(izb) = maxloc(rtau_y,dim=1)
          If ( abs(rtau_y(ints(izb))) > 0.0 ) Then
            tdel_dy(izb) = changey/rtau_y(ints(izb))
          Else
            tdel_dy(izb) = tdel_next(izb)
          EndIf

          ! Calculate timescale for temperature changes, T/(dT/dt)
          If ( iheat > 0 ) Then
            If ( mask_init(izb) ) Then
              rtau_t9 = abs(t9dot(izb)/t9(izb))
            Else
              rtau_t9 = abs((t9(izb)-t9o(izb))/t9(izb)) / tdel_old(izb)
            EndIf
          Else
            rtau_t9 = 0.0
          EndIf
          If ( rtau_t9 > 0.0 ) Then
            tdel_dt9(izb) = changet9/rtau_t9
          Else
            tdel_dt9(izb) = tdel_next(izb)
          EndIf

          ! Estimate timestep
          tdel(izb) = min( tdel_stop(izb), tdel_next(izb), tdel_dy(izb), tdel_dt9(izb) )

          ! Use the user-defined initial timestep if it is larger than the estimated timestep
          tdel(izb) = max( tdel_init, tdel(izb) )

        Else
          ! Keep timestep constant
          tdel(izb) = -tdel_old(izb)
        EndIf

        ! Update trial time
        tt(izb) = min(t(izb) + tdel(izb), tstop(izb))
      EndIf
    EndDo

    If ( idiag >= 2 ) Then
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          izone = izb + szbatch - zb_lo
          Write(lun_diag,"(a4,2i5,6es12.4,i5)") &
            & 'tdel',kstep,izone,tdel(izb),tdel_old(izb),tdel_stop(izb),tdel_next(izb), &
            & tdel_dy(izb),tdel_dt9(izb),ints(izb)
          !Write(lun_diag,"(a5,i4,2es12.4)") (nname(k),k,y(k,izb),rtau_y(k),k=1,ny)

          ! Retain the index of the species setting the timestep
          If ( ints(izb) /= intso(izb) ) Then
            Write(lun_diag,"(a4,a5,3es23.15)") 'ITC ',nname(ints(izb)),y(ints(izb),izb),t(izb),tdel(izb)
            intso(izb) = ints(izb)
          EndIf
        EndIf
      EndDo
    EndIf

    ! For varying temperature and density, capture thermodynamic features
    If ( iheat == 0 ) Then

      ! Make sure to not skip any features in the temperature or density profiles by checking
      ! for profile monotonicity between t and t+del
      Do izb = zb_lo, zb_hi
        If ( mask(izb) .and. nh(izb) > 1 ) Then
          Call t9rhofind(kstep,tt(izb),ntt(izb),t9t(izb),rhot(izb), &
            & nh(izb),th(:,izb),t9h(:,izb),rhoh(:,izb))
          If ( ntt(izb)-1 > nt(izb) ) Then
            Do j = nt(izb), ntt(izb)-1
              If ( t9h(j,izb) > t9(izb) .and. t9h(j,izb) > t9t(izb) ) Then
                tdel(izb) = th(j,izb) - t(izb)
                Exit
              ElseIf ( t9h(j,izb) < t9(izb) .and. t9h(j,izb) < t9t(izb) ) Then
                tdel(izb) = th(j,izb) - t(izb)
                Exit
              ElseIf ( rhoh(j,izb) > rho(izb) .and. rhoh(j,izb) > rhot(izb) ) Then
                tdel(izb) = th(j,izb) - t(izb)
                Exit
              ElseIf ( rhoh(j,izb) < rho(izb) .and. rhoh(j,izb) < rhot(izb) ) Then
                tdel(izb) = th(j,izb) - t(izb)
                Exit
              EndIf
            EndDo
            tt(izb) = min(t(izb) + tdel(izb), tstop(izb))
          EndIf

          ! Limit timestep if fractional density change is larger than changeth (10% by default)
          ! or fraction temperature change is larger than 0.1*changeth (1% by default)
          dtherm(izb) = 0.0
          Do i = 1, 10
            Call t9rhofind(kstep,tt(izb),ntt(izb),t9t(izb),rhot(izb), &
              & nh(izb),th(:,izb),t9h(:,izb),rhoh(:,izb))
            If ( t9(izb) > 0.0 ) Then
              dtherm(izb) = 10.0*abs(t9t(izb)-t9(izb))/t9(izb) + abs(rhot(izb)-rho(izb))/rho(izb)
            EndIf
            If ( dtherm(izb) < changeth ) Then
              Exit
            Else
              tdel(izb) = 0.5*tdel(izb)
              tt(izb) = min(t(izb) + tdel(izb), tstop(izb))
            EndIf
          EndDo
        EndIf
      EndDo
      If ( idiag >= 2 ) Then
        Do izb = zb_lo, zb_hi
          If ( mask(izb) .and. nh(izb) > 1 ) Then
            izone = izb + szbatch - zb_lo
            If ( dtherm(izb) >= changeth ) Then
              Write(lun_diag,"(a,i2,a,i5,3es24.16)") &
                & 'Error in Thermo variations after ',i,' reductions',izone,tdel(izb),t9t(izb),rhot(izb)
              Call xnet_terminate('Error in Thermo variations')
            EndIf
            Write(lun_diag,"(a5,2i5,2es23.15)") 'T9del',kstep,izone,tdel(izb),dtherm(izb)
          EndIf
        EndDo
      EndIf
    EndIf

    Return
  End Subroutine timestep

  Subroutine update_iweak(t9,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine updates the value of the iweak flag to control the treatment of strong reactions.
    !-----------------------------------------------------------------------------------------------
    Use xnet_controls, Only: iweak0, iweak, lun_stdout, t9min, zb_lo, zb_hi, lzactive
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9(zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer :: izb
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then

        ! Turn off strong interactions if the temperature is less than t9min
        If ( t9(izb) < t9min .and. iweak0 /= 0 ) Then
          iweak(izb) = -1
          Write(lun_stdout,*) 'Warning: Strong reactions ignored for T9 < T9min.'
        Else
          iweak(izb) = iweak0
        EndIf
      EndIf
    EndDo

    Return
  End Subroutine update_iweak

  Subroutine update_eos(mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine updates the dependent thermodynamic state variables by interfacing with the EoS.
    !-----------------------------------------------------------------------------------------------
    Use xnet_abundances, Only: yt, xext, aext, zext
    Use xnet_conditions, Only: t9t, rhot, yet, cv, etae, detaedt9
    Use xnet_controls, Only: zb_lo, zb_hi, lzactive
    Use xnet_eos, Only: eos_interface
    Use xnet_timers, Only: xnet_wtime, start_timer, stop_timer, timer_eos
    Use xnet_types, Only: dp
    Implicit None

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer :: izb
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    start_timer = xnet_wtime()
    timer_eos = timer_eos - start_timer

    Call eos_interface(t9t(zb_lo:zb_hi),rhot(zb_lo:zb_hi),yt(:,zb_lo:zb_hi), &
      & yet(zb_lo:zb_hi),cv(zb_lo:zb_hi),etae(zb_lo:zb_hi),detaedt9(zb_lo:zb_hi), &
      & xext(zb_lo:zb_hi),aext(zb_lo:zb_hi),zext(zb_lo:zb_hi),mask_in = mask)

    stop_timer = xnet_wtime()
    timer_eos = timer_eos + stop_timer

    Return
  End Subroutine update_eos

  Subroutine yderiv(mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates time derivatives for each nuclear species. This calculation is
    ! performed by looping over nuclei and summing the reaction rates for each reaction which
    ! effects that nucleus (see Equation 10 in Hix & Meyer (2006)).
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny, nname, mex
    Use reaction_data, Only: a1, a2, a3, a4, b1, b2, b3, b4, la, le, mu1, mu2, mu3, mu4, n11, n21, &
      & n22, n31, n32, n33, n41, n42, n43, n44, csect1, csect2, csect3, csect4, n1i, n2i, n3i, n4i
    Use xnet_abundances, Only: yt, ydot
    Use xnet_conditions, Only: cv, t9t, t9dot
    Use xnet_controls, Only: idiag, iheat, ktot, lun_diag, nzbatchmx, szbatch, zb_lo, zb_hi, &
      & lzactive
    Use xnet_timers, Only: xnet_wtime, start_timer, stop_timer, timer_deriv
    Use xnet_types, Only: dp
    Implicit None

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer :: i, j, i0, i1, la1, le1, la2, le2, la3, le3, la4, le4, izb, izone
    Real(dp) :: s1, s2, s3, s4, s11, s22, s33, s44, sdot
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    start_timer = xnet_wtime()
    timer_deriv = timer_deriv - start_timer

    ! From the cross sections and the counting array, calculate the reaction rates
    ! Calculate Ydot and T9dot for each nucleus, summing over the reactions which affect it.
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then

        Do i0 = 1, ny
          la1 = la(1,i0)
          la2 = la(2,i0)
          la3 = la(3,i0)
          la4 = la(4,i0)
          le1 = le(1,i0)
          le2 = le(2,i0)
          le3 = le(3,i0)
          le4 = le(4,i0)

          ! Sum over the reactions with 1 reactant
          s1 = 0.0
          Do i1 = la1, le1
            b1(i1,izb) = a1(i1)*csect1(mu1(i1),izb)
            s11 = b1(i1,izb)*yt(n11(i1),izb)
            s1 = s1 + s11
          EndDo

          ! Sum over the reactions with 2 reactants
          s2 = 0.0
          Do i1 = la2, le2
            b2(i1,izb) = a2(i1)*csect2(mu2(i1),izb)
            s22 = b2(i1,izb)*yt(n21(i1),izb)*yt(n22(i1),izb)
            s2 = s2 + s22
          EndDo

          ! Sum over the reactions with 3 reactants
          s3 = 0.0
          Do i1 = la3, le3
            b3(i1,izb) = a3(i1)*csect3(mu3(i1),izb)
            s33 = b3(i1,izb)*yt(n31(i1),izb)*yt(n32(i1),izb)*yt(n33(i1),izb)
            s3 = s3 + s33
          EndDo

          ! Sum over the reactions with 4 reactants
          s4 = 0.0
          Do i1 = la4, le4
            b4(i1,izb) = a4(i1)*csect4(mu4(i1),izb)
            s44 = b4(i1,izb)*yt(n41(i1),izb)*yt(n42(i1),izb)*yt(n43(i1),izb)*yt(n44(i1),izb)
            s4 = s4 + s44
          EndDo

          ! Sum the 4 components of Ydot
          ydot(i0,izb) = s1 + s2 + s3 + s4
        EndDo
      EndIf
    EndDo

    If ( iheat > 0 ) Then

      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          sdot = 0.0
          Do i0 = 1, ny
            sdot = sdot - mex(i0)*ydot(i0,izb) / cv(izb)
          EndDo
          t9dot(izb) = sdot
        EndIf
      EndDo
    EndIf

    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        ktot(4,izb) = ktot(4,izb) + 1
      EndIf
    EndDo

    ! Separate loop for diagnostics so compiler can properly vectorize
    If ( idiag >= 5 ) Then
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          izone = izb + szbatch - zb_lo
          Write(lun_diag,"(a,i5)") 'YDERIV',izone
          Do i0 = 1, ny
            la1 = la(1,i0)
            la2 = la(2,i0)
            la3 = la(3,i0)
            la4 = la(4,i0)
            le1 = le(1,i0)
            le2 = le(2,i0)
            le3 = le(3,i0)
            le4 = le(4,i0)
            s1 = 0.0
            s2 = 0.0
            s3 = 0.0
            s4 = 0.0
            Do i1 = la1, le1
              s11 = b1(i1,izb)*yt(n11(i1),izb)
              s1 = s1 + s11
            EndDo
            Do i1 = la2, le2
              s22 = b2(i1,izb)*yt(n21(i1),izb)*yt(n22(i1),izb)
              s2 = s2 + s22
            EndDo
            Do i1 = la3, le3
              s33 = b3(i1,izb)*yt(n31(i1),izb)*yt(n32(i1),izb)*yt(n33(i1),izb)
              s3 = s3 + s33
            EndDo
            Do i1 = la4, le4
              s44 = b4(i1,izb)*yt(n41(i1),izb)*yt(n42(i1),izb)*yt(n43(i1),izb)*yt(n44(i1),izb)
              s4 = s4 + s44
            EndDo
            If ( idiag >= 6 ) Then
              Write(lun_diag,"(a3,a6,i4)") 'NUC',nname(i0),i0
              Write(lun_diag,"(5a5,'  1  ',4es23.15)") &
                & ((nname(n1i(j,mu1(i1))),j=1,5),b1(i1,izb),yt(n11(i1),izb), &
                & b1(i1,izb)*yt(n11(i1),izb),a1(i1),i1=la1,le1)
              Write(lun_diag,*) '1->',nname(i0),la1,le1,s1
              Write(lun_diag,"(6a5,'  2  ',5es23.15)") &
                & ((nname(n2i(i,mu2(i1))),i=1,6),b2(i1,izb),yt(n21(i1),izb),yt(n22(i1),izb), &
                & b2(i1,izb)*yt(n21(i1),izb)*yt(n22(i1),izb),a2(i1),i1=la2,le2)
              Write(lun_diag,*) '2->',nname(i0),la2,le2,s2
              Write(lun_diag,"(6a5,'  3  ',6es23.15)") &
                & ((nname(n3i(i,mu3(i1))),i=1,6),b3(i1,izb),yt(n31(i1),izb),yt(n32(i1),izb),yt(n33(i1),izb), &
                & b3(i1,izb)*yt(n31(i1),izb)*yt(n32(i1),izb)*yt(n33(i1),izb),a3(i1),i1=la3,le3)
              Write(lun_diag,*) '3->',nname(i0),la3,le3,s3
              Write(lun_diag,"(6a5,'  4  ',7es23.15)") &
                & ((nname(n4i(i,mu4(i1))),i=1,6),b4(i1,izb),yt(n41(i1),izb),yt(n42(i1),izb),yt(n43(i1),izb),yt(n44(i1),izb), &
                & b4(i1,izb)*yt(n41(i1),izb)*yt(n42(i1),izb)*yt(n43(i1),izb)*yt(n44(i1),izb),a4(i1),i1=la4,le4)
              Write(lun_diag,*) '4->',nname(i0),la4,le4,s4
            EndIf
            Write(lun_diag,"(a4,a5,2es23.15,4es11.3)") 'YDOT',nname(i0),yt(i0,izb),ydot(i0,izb),s1,s2,s3,s4
          EndDo
          If ( iheat > 0 ) Write(lun_diag,"(5x,a5,2es23.15)") '   T9',t9t(izb),t9dot(izb)
        EndIf
      EndDo
    EndIf

    stop_timer = xnet_wtime()
    timer_deriv = timer_deriv + stop_timer

    Return
  End Subroutine yderiv

  Subroutine cross_sect(mask_in)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the cross section for each reaction.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: nname, gg, dlngdt9, partf
    Use reaction_data
    Use xnet_abundances, Only: yt
    Use xnet_conditions, Only: tt, rhot, t9t, yet
    Use xnet_constants, Only: third, two3rd, four3rd, five3rd
    Use xnet_controls, Only: idiag, iheat, iscrn, iweak, ktot, lun_diag, nzbatchmx, szbatch, &
      & zb_lo, zb_hi, lzactive
    Use xnet_ffn, Only: ffn_rate, rffn, dlnrffndt9
    Use xnet_nnu, Only: nnu_rate, nnuspec, rnnu
    Use xnet_screening, Only: h1, h2, h3, h4, dh1dt9, dh2dt9, dh3dt9, dh4dt9, screening
    Use xnet_timers, Only: xnet_wtime, start_timer, stop_timer, timer_csect
    Use xnet_types, Only: dp
    Use xnet_util, Only: safe_exp
    Implicit None

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Logical, Parameter :: use_blas = .true.
    Integer, Parameter :: dgemm_nzbatch = 200 ! Min number of zones to use dgemm, otherwise use dgemv
    Real(dp) :: t09(7,zb_lo:zb_hi), dt09(7,zb_lo:zb_hi)
    Real(dp) :: ene(zb_lo:zb_hi), ytot, abar, zbar, z2bar, zibar
    Real(dp) :: ascrn, rhot2, rhot3
    Real(dp) :: rpf1, rpf2, rpf3, rpf4
    Real(dp) :: dlnrpf1dt9, dlnrpf2dt9, dlnrpf3dt9, dlnrpf4dt9
    Real(dp) :: lambda1, lambda2, lambda3, lambda4
    Real(dp) :: dlam1dt9, dlam2dt9, dlam3dt9, dlam4dt9
    Integer :: j, k, izb, izone, nzmask
    Integer :: nr1, nr2, nr3, nr4
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return
    nzmask = count(mask)

    nr1 = nreac(1)
    nr2 = nreac(2)
    nr3 = nreac(3)
    nr4 = nreac(4)

    ! Update thermodynamic state
    Call update_eos(mask_in = mask)

    ! Calculate the screening terms
    If ( iscrn >= 1 ) Then
      ascrn = 1.0
      Call screening(mask_in = mask)
    Else
      ascrn = 0.0
    EndIf

    start_timer = xnet_wtime()
    timer_csect = timer_csect - start_timer

    ! Check for any changes to iweak
    Call update_iweak(t9t(zb_lo:zb_hi),mask_in = mask)

    ! Calculate partition functions for each nucleus at t9t
    Call partf(t9t(zb_lo:zb_hi),mask_in = mask)

    ! Calculate necessary thermodynamic moments
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then
        ene(izb) = yet(izb)*rhot(izb)

        t09(1,izb) = 1.0
        t09(2,izb) = t9t(izb)**(-1)
        t09(3,izb) = t9t(izb)**(-1.0/3.0)
        t09(4,izb) = t9t(izb)**(+1.0/3.0)
        t09(5,izb) = t9t(izb)
        t09(6,izb) = t9t(izb)**(+5.0/3.0)
        t09(7,izb) = log(t9t(izb))

        ! Calculate reaction rate coefficient derivatives
        dt09(1,izb) = 0.0
        dt09(2,izb) = -t9t(izb)**(-2)
        dt09(3,izb) = -t9t(izb)**(-4.0/3.0) / 3.0
        dt09(4,izb) = +t9t(izb)**(-2.0/3.0) / 3.0
        dt09(5,izb) = 1.0
        dt09(6,izb) = +t9t(izb)**(+2.0/3.0) * 5.0/3.0
        dt09(7,izb) = 1.0 / t9t(izb)
      EndIf
    EndDo

    ! If there are any FFN reactions, calculate their rates
    If ( nffn > 0 ) Call ffn_rate(nffn,t9t(zb_lo:zb_hi),ene, &
      & rffn(:,zb_lo:zb_hi),dlnrffndt9(:,zb_lo:zb_hi),mask_in = mask)

    ! If there are any neutrino-nucleus reactions, calculate their rates
    If ( nnnu > 0 ) Call nnu_rate(nnnu,tt(zb_lo:zb_hi),rnnu(:,:,zb_lo:zb_hi),mask_in = mask)

    ! Calculate the REACLIB exponent polynomials and derivatives, adding screening terms
    If ( use_blas ) Then
      If ( nzmask >= dgemm_nzbatch ) Then
        If ( nr1 > 0 ) Call dgemm('T','N',nr1,nzbatchmx,7,1.0,rc1,7,t09,7,ascrn,h1(:,zb_lo:zb_hi),nr1)
        If ( nr2 > 0 ) Call dgemm('T','N',nr2,nzbatchmx,7,1.0,rc2,7,t09,7,ascrn,h2(:,zb_lo:zb_hi),nr2)
        If ( nr3 > 0 ) Call dgemm('T','N',nr3,nzbatchmx,7,1.0,rc3,7,t09,7,ascrn,h3(:,zb_lo:zb_hi),nr3)
        If ( nr4 > 0 ) Call dgemm('T','N',nr4,nzbatchmx,7,1.0,rc4,7,t09,7,ascrn,h4(:,zb_lo:zb_hi),nr4)

        If ( iheat > 0 ) Then
          If ( nr1 > 0 ) Call dgemm('T','N',nr1,nzbatchmx,7,1.0,rc1,7,dt09,7,ascrn,dh1dt9(:,zb_lo:zb_hi),nr1)
          If ( nr2 > 0 ) Call dgemm('T','N',nr2,nzbatchmx,7,1.0,rc2,7,dt09,7,ascrn,dh2dt9(:,zb_lo:zb_hi),nr2)
          If ( nr3 > 0 ) Call dgemm('T','N',nr3,nzbatchmx,7,1.0,rc3,7,dt09,7,ascrn,dh3dt9(:,zb_lo:zb_hi),nr3)
          If ( nr4 > 0 ) Call dgemm('T','N',nr4,nzbatchmx,7,1.0,rc4,7,dt09,7,ascrn,dh4dt9(:,zb_lo:zb_hi),nr4)
        EndIf
      Else
        Do izb = zb_lo, zb_hi
          If ( mask(izb) ) Then
            If ( nr1 > 0 ) Call dgemv('T',7,nr1,1.0,rc1,7,t09(:,izb),1,ascrn,h1(:,izb),1)
            If ( nr2 > 0 ) Call dgemv('T',7,nr2,1.0,rc2,7,t09(:,izb),1,ascrn,h2(:,izb),1)
            If ( nr3 > 0 ) Call dgemv('T',7,nr3,1.0,rc3,7,t09(:,izb),1,ascrn,h3(:,izb),1)
            If ( nr4 > 0 ) Call dgemv('T',7,nr4,1.0,rc4,7,t09(:,izb),1,ascrn,h4(:,izb),1)
            If ( iheat > 0 ) Then
              If ( nr1 > 0 ) Call dgemv('T',7,nr1,1.0,rc1,7,dt09(:,izb),1,ascrn,dh1dt9(:,izb),1)
              If ( nr2 > 0 ) Call dgemv('T',7,nr2,1.0,rc2,7,dt09(:,izb),1,ascrn,dh2dt9(:,izb),1)
              If ( nr3 > 0 ) Call dgemv('T',7,nr3,1.0,rc3,7,dt09(:,izb),1,ascrn,dh3dt9(:,izb),1)
              If ( nr4 > 0 ) Call dgemv('T',7,nr4,1.0,rc4,7,dt09(:,izb),1,ascrn,dh4dt9(:,izb),1)
            EndIf
          EndIf
        EndDo
      EndIf
    Else
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          Do k = 1, nr1
            lambda1 = 0.0
            dlam1dt9 = 0.0
            Do j = 1, 7
              lambda1 = lambda1 + rc1(j,k)*t09(j,izb)
              dlam1dt9 = dlam1dt9 + rc1(j,k)*dt09(j,izb)
            EndDo
            h1(k,izb) = ascrn*h1(k,izb) + lambda1
            dh1dt9(k,izb) = ascrn*dh1dt9(k,izb) + dlam1dt9
          EndDo
          Do k = 1, nr2
            lambda2 = 0.0
            dlam2dt9 = 0.0
            Do j = 1, 7
              lambda2 = lambda2 + rc2(j,k)*t09(j,izb)
              dlam2dt9 = dlam2dt9 + rc2(j,k)*dt09(j,izb)
            EndDo
            h2(k,izb) = ascrn*h2(k,izb) + lambda2
            dh2dt9(k,izb) = ascrn*dh2dt9(k,izb) + dlam2dt9
          EndDo
          Do k = 1, nr3
            lambda3 = 0.0
            dlam3dt9 = 0.0
            Do j = 1, 7
              lambda3 = lambda3 + rc3(j,k)*t09(j,izb)
              dlam3dt9 = dlam3dt9 + rc3(j,k)*dt09(j,izb)
            EndDo
            h3(k,izb) = ascrn*h3(k,izb) + lambda3
            dh3dt9(k,izb) = ascrn*dh3dt9(k,izb) + dlam3dt9
          EndDo
          Do k = 1, nr4
            lambda4 = 0.0
            dlam4dt9 = 0.0
            Do j = 1, 7
              lambda4 = lambda4 + rc4(j,k)*t09(j,izb)
              dlam4dt9 = dlam4dt9 + rc4(j,k)*dt09(j,izb)
            EndDo
            h4(k,izb) = ascrn*h4(k,izb) + lambda4
            dh4dt9(k,izb) = ascrn*dh4dt9(k,izb) + dlam4dt9
          EndDo
        EndIf
      EndDo
    EndIf

    ! Calculate cross sections
    Do izb = zb_lo, zb_hi
      If ( mask(izb) ) Then

        ! Calculate the csect for reactions with 1 reactant
        Do k = 1, nr1
          If ( irev1(k) == 1 ) Then
            rpf1 =   ( gg(n1i(2,k),izb) * gg(n1i(3,k),izb) * gg(n1i(4,k),izb) * gg(n1i(5,k),izb) ) &
              &    / ( gg(n1i(1,k),izb) )
          Else
            rpf1 = 1.0
          EndIf
          If ( iweak(izb) > 0 ) Then
            If ( iwk1(k) == 0 .or. iwk1(k) == 4 ) Then
              csect1(k,izb) = rpf1 * safe_exp(h1(k,izb))
            ElseIf ( iwk1(k) == 1 ) Then
              csect1(k,izb) = rpf1 * safe_exp(h1(k,izb)) * ene(izb)
            ElseIf ( iwk1(k) == 2 .or. iwk1(k) == 3 ) Then ! FFN reaction
              csect1(k,izb) = rffn(iffn(k),izb)
            ElseIf ( iwk1(k) == 7 ) Then ! Electron neutrino capture
              csect1(k,izb) = rpf1 * rnnu(innu(k),1,izb)
            ElseIf ( iwk1(k) == 8 ) Then ! Electron anti-neutrino capture
              csect1(k,izb) = rpf1 * rnnu(innu(k),2,izb)
            EndIf
          ElseIf ( iweak(izb) < 0 ) Then
            If ( iwk1(k) == 0 ) Then
              csect1(k,izb) = 0.0
            ElseIf ( iwk1(k) == 1 ) Then
              csect1(k,izb) = rpf1 * safe_exp(h1(k,izb)) * ene(izb)
            ElseIf ( iwk1(k) == 4 ) Then
              csect1(k,izb) = rpf1 * safe_exp(h1(k,izb))
            ElseIf ( iwk1(k) == 2 .or. iwk1(k) == 3 ) Then ! FFN reaction
              csect1(k,izb) = rffn(iffn(k),izb)
            ElseIf ( iwk1(k) == 7 ) Then ! Electron neutrino capture
              csect1(k,izb) = rpf1 * rnnu(innu(k),1,izb)
            ElseIf ( iwk1(k) == 8 ) Then ! Electron anti-neutrino capture
              csect1(k,izb) = rpf1 * rnnu(innu(k),2,izb)
            EndIf
          Else
            If ( iwk1(k) == 0 ) Then
              csect1(k,izb) = rpf1 * safe_exp(h1(k,izb))
            Else
              csect1(k,izb) = 0.0
            EndIf
          EndIf
        EndDo

        ! Calculate the csect for reactions with 2 reactants
        Do k = 1, nr2
          If ( irev2(k) == 1 ) Then
            rpf2 =   ( gg(n2i(3,k),izb) * gg(n2i(4,k),izb) * gg(n2i(5,k),izb) * gg(n2i(6,k),izb) ) &
              &    / ( gg(n2i(1,k),izb) * gg(n2i(2,k),izb) )
          Else
            rpf2 = 1.0
          EndIf
          If ( iweak(izb) > 0 ) Then
            If ( iwk2(k) == 1 ) Then
              csect2(k,izb) = rhot(izb) * rpf2 * safe_exp(h2(k,izb)) * ene(izb)
            Else
              csect2(k,izb) = rhot(izb) * rpf2 * safe_exp(h2(k,izb))
            EndIf
          ElseIf ( iweak(izb) < 0 ) Then
            If ( iwk2(k) == 0 ) Then
              csect2(k,izb) = 0.0
            ElseIf ( iwk2(k) == 1 ) Then
              csect2(k,izb) = rhot(izb) * rpf2 * safe_exp(h2(k,izb)) * ene(izb)
            Else
              csect2(k,izb) = rhot(izb) * rpf2 * safe_exp(h2(k,izb))
            EndIf
          Else
            If ( iwk2(k) == 0 ) Then
              csect2(k,izb) = rhot(izb) * rpf2 * safe_exp(h2(k,izb))
            Else
              csect2(k,izb) = 0.0
            EndIf
          EndIf
        EndDo

        ! Calculate the csect for reactions with 3 reactants
        Do k = 1, nr3
          If ( irev3(k) == 1 ) Then
            rpf3 =   ( gg(n3i(4,k),izb) * gg(n3i(5,k),izb) * gg(n3i(6,k),izb) ) &
              &    / ( gg(n3i(1,k),izb) * gg(n3i(2,k),izb) * gg(n3i(3,k),izb) )
          Else
            rpf3 = 1.0
          EndIf
          rhot2 = rhot(izb)**2
          If ( iweak(izb) > 0 ) Then
            If ( iwk3(k) == 1 ) Then
              csect3(k,izb) = rhot2 * rpf3 * safe_exp(h3(k,izb)) * ene(izb)
            Else
              csect3(k,izb) = rhot2 * rpf3 * safe_exp(h3(k,izb))
            EndIf
          ElseIf ( iweak(izb) < 0 ) Then
            If ( iwk3(k) == 0 ) Then
              csect3(k,izb) = 0.0
            ElseIf ( iwk3(k) == 1 ) Then
              csect3(k,izb) = rhot2 * rpf3 * safe_exp(h3(k,izb)) * ene(izb)
            Else
              csect3(k,izb) = rhot2 * rpf3 * safe_exp(h3(k,izb))
            EndIf
          Else
            If ( iwk3(k) == 0 ) Then
              csect3(k,izb) = rhot2 * rpf3 * safe_exp(h3(k,izb))
            Else
              csect3(k,izb) = 0.0
            EndIf
          EndIf
          If ( csect3(k,izb) < 1.0e-20 ) csect3(k,izb) = 0.0
        EndDo

        ! Calculate the csect for reactions with 4 reactants
        Do k = 1, nr4
          If ( irev4(k) == 1 ) Then
            rpf4 =   ( gg(n4i(5,k),izb) * gg(n4i(6,k),izb) ) &
              &    / ( gg(n4i(1,k),izb) * gg(n4i(2,k),izb) * gg(n4i(3,k),izb) * gg(n4i(4,k),izb ) )
          Else
            rpf4 = 1.0
          EndIf
          rhot3 = rhot(izb)**3
          If ( iweak(izb) > 0 ) Then
            If ( iwk4(k) == 1 ) Then
              csect4(k,izb) = rhot3 * rpf4 * safe_exp(h4(k,izb)) * ene(izb)
            Else
              csect4(k,izb) = rhot3 * rpf4 * safe_exp(h4(k,izb))
            EndIf
          ElseIf ( iweak(izb) < 0 ) Then
            If ( iwk4(k) == 0 ) Then
              csect4(k,izb) = 0.0
            ElseIf ( iwk4(k) == 1 ) Then
              csect4(k,izb) = rhot3 * rpf4 * safe_exp(h4(k,izb)) * ene(izb)
            Else
              csect4(k,izb) = rhot3 * rpf4 * safe_exp(h4(k,izb))
            EndIf
          Else
            If ( iwk4(k) == 0 ) Then
              csect4(k,izb) = rhot3 * rpf4 * safe_exp(h4(k,izb))
            Else
              csect4(k,izb) = 0.0
            EndIf
          EndIf
        EndDo

        ! Increment counter
        ktot(5,izb) = ktot(5,izb) + 1
      EndIf
    EndDo

    If ( iheat > 0 ) Then


      ! Calculate cross-section derivatives
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then

          ! 1-reactant reactions
          Do k = 1, nr1
            If ( irev1(k) == 1 ) Then
              dlnrpf1dt9 =   (   dlngdt9(n1i(2,k),izb) + dlngdt9(n1i(3,k),izb) &
                &              + dlngdt9(n1i(4,k),izb) + dlngdt9(n1i(5,k),izb) ) &
                &          - (   dlngdt9(n1i(1,k),izb) )
            Else
              dlnrpf1dt9 = 0.0
            EndIf
            If ( iwk1(k) == 2 .or. iwk1(k) == 3 ) Then  ! FFN reaction
              dcsect1dt9(k,izb) = rffn(iffn(k),izb)*dlnrffndt9(iffn(k),izb)
            ElseIf ( iwk1(k) == 7 ) Then  ! Electron neutrino capture
              dcsect1dt9(k,izb) = rnnu(innu(k),1,izb)*dlnrpf1dt9
            ElseIf ( iwk1(k) == 8 ) Then  ! Electron anti-neutrino capture
              dcsect1dt9(k,izb) = rnnu(innu(k),2,izb)*dlnrpf1dt9
            Else
              dcsect1dt9(k,izb) = csect1(k,izb)*(dh1dt9(k,izb)+dlnrpf1dt9)
            EndIf
          EndDo

          ! 2-reactant reactions
          Do k = 1, nr2
            If ( irev2(k) == 1 ) Then
              dlnrpf2dt9 =   (   dlngdt9(n2i(3,k),izb) + dlngdt9(n2i(4,k),izb) &
                               + dlngdt9(n2i(5,k),izb) + dlngdt9(n2i(6,k),izb) ) &
                &          - (   dlngdt9(n2i(1,k),izb) + dlngdt9(n2i(2,k),izb) )
            Else
              dlnrpf2dt9 = 0.0
            EndIf
            dcsect2dt9(k,izb) = csect2(k,izb)*(dh2dt9(k,izb)+dlnrpf2dt9)
          EndDo

          ! 3-reactant reactions
          Do k = 1, nr3
            If ( irev3(k) == 1 ) Then
              dlnrpf3dt9 =   ( dlngdt9(n3i(4,k),izb) + dlngdt9(n3i(5,k),izb) + dlngdt9(n3i(6,k),izb) ) &
                &          - ( dlngdt9(n3i(1,k),izb) + dlngdt9(n3i(2,k),izb) + dlngdt9(n3i(3,k),izb) )
            Else
              dlnrpf3dt9 = 0.0
            EndIf
            dcsect3dt9(k,izb) = csect3(k,izb)*(dh3dt9(k,izb)+dlnrpf3dt9)
          EndDo

          ! 4-reactant reactions
          Do k = 1, nr4
            If ( irev4(k) == 1 ) Then
              dlnrpf4dt9 =   (   dlngdt9(n4i(5,k),izb) + dlngdt9(n4i(6,k),izb) ) &
                &          - (   dlngdt9(n4i(1,k),izb) + dlngdt9(n4i(2,k),izb) &
                &              + dlngdt9(n4i(3,k),izb) + dlngdt9(n4i(4,k),izb) )
            Else
              dlnrpf4dt9 = 0.0
            EndIf
            dcsect4dt9(k,izb) = csect4(k,izb)*(dh4dt9(k,izb)+dlnrpf4dt9)
          EndDo
        EndIf
      EndDo
    EndIf

    If ( idiag >= 6 ) Then
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          izone = izb + szbatch - zb_lo
          Write(lun_diag,"(a,i5)") 'CSect',izone
          Write(lun_diag,"(a,i5)") 'CSect1',nr1
          Write(lun_diag,"(i5,6a5,3i3,2es23.15)") &
            & (k,nname(n1i(1,k)),'-->',(nname(n1i(j,k)),j=2,5), &
            & iwk1(k),irev1(k),ires1(k),h1(k,izb),csect1(k,izb), k=1,nr1)
          Write(lun_diag,"(a,i5)") 'CSect2',nr2
          Write(lun_diag,"(i5,7a5,3i3,2es23.15)") &
            & (k,(nname(n2i(j,k)),j=1,2),'-->',(nname(n2i(j,k)),j=3,6), &
            & iwk2(k),irev2(k),ires2(k),h2(k,izb),csect2(k,izb),k=1,nr2)
          Write(lun_diag,"(a,i5)") 'CSect3',nr3
          Write(lun_diag,"(i5,7a5,3i3,2es23.15)") &
            & (k,(nname(n3i(j,k)),j=1,3),'-->',(nname(n3i(j,k)),j=4,6), &
            & iwk3(k),irev3(k),ires3(k),h3(k,izb),csect3(k,izb),k=1,nr3)
          Write(lun_diag,"(a,i5)") 'CSect4',nr4
          Write(lun_diag,"(i5,7a5,3i3,2es23.15)") &
            & (k,(nname(n4i(j,k)),j=1,4),'-->',(nname(n4i(j,k)),j=5,6), &
            & iwk4(k),irev4(k),ires4(k),h4(k,izb),csect4(k,izb),k=1,nr4)
          If ( iheat > 0 ) Then
            Write(lun_diag,"(a,i5)") 'dCSect1/dT9',nr1
            Write(lun_diag,"(i5,6a5,3i3,2es23.15)") &
              & (k,nname(n1i(1,k)),'-->',(nname(n1i(j,k)),j=2,5), &
              & iwk1(k),irev1(k),ires1(k),dh1dt9(k,izb),dcsect1dt9(k,izb), k=1,nr1)
            Write(lun_diag,"(a,i5)") 'dCSect2/dT9',nr2
            Write(lun_diag,"(i5,7a5,3i3,2es23.15)") &
              & (k,(nname(n2i(j,k)),j=1,2),'-->',(nname(n2i(j,k)),j=3,6), &
              & iwk2(k),irev2(k),ires2(k),dh2dt9(k,izb),dcsect2dt9(k,izb),k=1,nr2)
            Write(lun_diag,"(a,i5)") 'dCSect3/dT9',nr3
            Write(lun_diag,"(i5,7a5,3i3,2es23.15)") &
              & (k,(nname(n3i(j,k)),j=1,3),'-->',(nname(n3i(j,k)),j=4,6), &
              & iwk3(k),irev3(k),ires3(k),dh3dt9(k,izb),dcsect3dt9(k,izb),k=1,nr3)
            Write(lun_diag,"(a,i5)") 'dCSect4/dT9',nr4
            Write(lun_diag,"(i5,7a5,3i3,2es23.15)") &
              & (k,(nname(n4i(j,k)),j=1,4),'-->',(nname(n4i(j,k)),j=5,6), &
              & iwk4(k),irev4(k),ires4(k),dh4dt9(k,izb),dcsect4dt9(k,izb),k=1,nr4)
          EndIf
        EndIf
      EndDo
    EndIf

    stop_timer = xnet_wtime()
    timer_csect = timer_csect + stop_timer

    Return
  End Subroutine cross_sect

End Module xnet_integrate
