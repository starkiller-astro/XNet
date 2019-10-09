!***************************************************************************************************
! xnet_evolve.f90 10/18/17
! The subroutines in this file perform nucleosynthesis for thermodynamic trajectories (usually
! Lagrangian mass particles or finite-volume cells). Multiple trajectories can be involved in step
! with "batching", but for large numbers of trajectories still need to be managed externally.
!***************************************************************************************************

Module xnet_evolve
  !-------------------------------------------------------------------------------------------------
  ! This module contains the data and routines to perform the evolution of the reaction network
  ! along a thermodynamic trajectory.
  !-------------------------------------------------------------------------------------------------
  Implicit None

Contains

  Subroutine full_net(kstep)
    !-----------------------------------------------------------------------------------------------
    ! The abundance evolution is performed over a series of timesteps, with the duration of the
    ! timestep determined by the integration scheme and the changing thermodynamic conditions.
    ! Integration is performed by a choice of methods controlled by the isolv flag.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny, nname
    Use xnet_abundances, Only: yo, y, yt, ydot
    Use xnet_conditions, Only: t, to, tt, tdel, tdel_old, tdel_next, t9, t9o, t9t, t9dot, &
      & rho, rhoo, rhot, yeo, ye, yet, nt, nto, ntt, tstop, t9rhofind
    Use xnet_controls, Only: idiag, iheat, isolv, itsout, kstmx, kmon, ktot, lun_diag, lun_stdout, &
      & lzactive, szbatch, nzbatchmx, nzevolve, zb_lo, zb_hi, tid
    Use xnet_integrate, Only: timestep
    Use xnet_integrate_be, Only: solve_be
    Use xnet_integrate_bdf, Only: solve_bdf
    Use xnet_output, Only: final_output, ts_output
    Use xnet_timers, Only: xnet_wtime, start_timer, stop_timer, timer_xnet
    Use xnet_types, Only: dp
    Use xnet_util, Only: xnet_terminate
    Implicit None

    ! Output variables
    Integer, Intent(out) :: kstep

    ! Local variables
    !Integer, Parameter :: kstep_output = 10
    Real(dp) :: enm, enb, enold(zb_lo:zb_hi), en0(zb_lo:zb_hi)
    Integer :: k, izb, izone, nstep_est, idiag0
    Integer :: its(zb_lo:zb_hi), mykstep(zb_lo:zb_hi)
    Logical :: lzsolve(zb_lo:zb_hi), lzoutput(zb_lo:zb_hi)

    start_timer = xnet_wtime()
    timer_xnet = timer_xnet - start_timer

    ! Set reaction controls not read in from control
    idiag0 = idiag

    ! Initialize trial time step abundances and conditions
    !$acc parallel loop gang async(tid) &
    !$acc present(kmon,ktot,tdel,tdel_old,tdel_next,nt,nto,ntt, &
    !$acc         t9,t9o,t9t,rho,rhoo,rhot,t,to,tt,ye,yeo,yet,y,yo,yt)
    Do izb = zb_lo, zb_hi
      tdel_old(izb) = tdel(izb)
      tdel_next(izb) = tdel(izb)
      nto(izb) = nt(izb)
      ntt(izb) = nt(izb)
      to(izb) = t(izb)
      tt(izb) = t(izb)
      t9o(izb) = t9(izb)
      t9t(izb) = t9(izb)
      rhoo(izb) = rho(izb)
      rhot(izb) = rho(izb)
      yeo(izb) = ye(izb)
      yet(izb) = ye(izb)
      !$acc loop vector
      Do k = 1, ny
        yo(k,izb) = y(k,izb)
        yt(k,izb) = y(k,izb)
      EndDo
      !$acc loop vector
      Do k = 1, 5
        kmon(k,izb) = 0
        ktot(k,izb) = 0
      EndDo
    EndDo

    ! Output initial abundances and conditions
    Call ts_output(0,en0,enold)

    ! Initialize timestep loop flags
    kstep = 0
    mykstep = 0
    lzsolve = lzactive(zb_lo:zb_hi)
    lzoutput = lzsolve
    Do izb = zb_lo, zb_hi
      If ( lzsolve(izb) ) Then
        its(izb) = 0
      Else
        its(izb) = -1
      EndIf
    EndDo

    ! Start evolution
    !$acc enter data async(tid) &
    !$acc copyin(its,mykstep,lzsolve,lzoutput)
    Do kstep = 1, kstmx

      ! Determine if this is an output step
      idiag = idiag0
      !If ( mod(kstep,kstep_output) == 0 ) idiag = 2

      ! Calculate an initial guess for the timestep
      Call timestep(kstep,mask_in = lzsolve)

      ! Take integration step (only worry about solve_be for now)
      Select Case (isolv)
      Case (3)
        Call solve_bdf(kstep,its)
      !Case (2)
      !  Call solve_bd(kstep,its)
      Case Default
        Call solve_be(kstep,its)
      End Select

      ! If convergence is successful, output timestep results
      If ( idiag >= 1 ) Then
        !$acc update wait(tid) &
        !$acc host(its)
        Do izb = zb_lo, zb_hi
          izone = izb + szbatch - zb_lo
          If ( its(izb) == 0 .and. idiag >= 1 ) Then
            !$acc update wait(tid) &
            !$acc host(t(izb),tdel(izb),t9(izb),rho(izb),ye(izb))
            Write(lun_diag,"(2(a,i5),5es14.7)") &
              & 'KStep ',kstep,' Zone ',izone,t(izb),tdel(izb),t9(izb),rho(izb),ye(izb)
            If ( idiag >= 3 ) Then
              !$acc update wait(tid) &
              !$acc host(y(:,izb),yo(:,izb),ydot(:,izb),tdel(izb),&
              !$acc      t9(izb),t9o(izb),t9dot(izb))
              Write(lun_diag,"(a8)") 'delta Y'
              Write(lun_diag,"(8x,a5,4es12.4)") &
                & (nname(k),y(k,izb),yo(k,izb),(y(k,izb)-yo(k,izb)),(tdel(izb)*ydot(k,izb)),k=1,ny)
              If ( iheat > 0 ) Write(lun_diag,"(8x,a5,4es12.4)") &
                & 'T9',t9(izb),t9o(izb),t9(izb)-t9o(izb),tdel(izb)*t9dot(izb)
            EndIf
          ElseIf ( its(izb) == 1 ) Then
            Write(lun_diag,"(a,2(i5,a))") 'KStep ',kstep,' Zone ',izone,' Inter!!!'
          EndIf
        EndDo
        Flush(lun_diag)
      EndIf

      !$acc parallel loop gang async(tid) &
      !$acc present(its,t,tstop,mykstep,lzsolve,lzoutput)
      Do izb = zb_lo, zb_hi
        If ( its(izb) == 0 ) Then

          ! If this zone reaches the stop time, flag it to remove from loop
          If ( t(izb) >= tstop(izb) ) Then
            mykstep(izb) = kstep
            its(izb) = -1
            lzsolve(izb) = .false.
          EndIf

        ! If reduced timesteps fail to successfully integrate, flag it to remove from loop
        ElseIf ( its(izb) == 1 ) Then
          its(izb) = 2
          lzsolve(izb) = .false.
          lzoutput(izb) = .false.
        EndIf
      EndDo

      !$acc update wait(tid) &
      !$acc host(lzoutput,lzsolve)
      Call ts_output(kstep,en0,enold,mask_in = lzoutput)

      ! Test if all zones have stopped
      If ( .not. any( lzsolve ) ) Exit
    EndDo

    !$acc exit data wait(tid) &
    !$acc copyout(its,mykstep) &
    !$acc delete(lzsolve,lzoutput)

    ! Test that the stop time is reached
    Do izb = zb_lo, zb_hi
      If ( lzactive(izb) ) Then
        !$acc update wait(tid) &
        !$acc host(t(izb),tdel(izb))
        If ( t(izb) < tstop(izb) .or. its(izb) > 0 ) Then
          izone = izb + szbatch - zb_lo
          Write(lun_stdout,"(a,i5,a,3(es12.4,a))") &
            & 'Zone ',izone,' Evolution stopped at time= ',t(izb),' with timestep= ',tdel(izb), &
            & ', stop time= ',tstop(izb),' not reached!'
          nstep_est = int( min((tstop(izb)-t(izb))/tdel(izb), real(huge(nstep_est),dp) ) )
          Write(lun_stdout,"(a,i12,a)") &
            & 'Approximately ',nstep_est,' more steps needed'
          Call xnet_terminate('[XNet] Evolution failed to converge')
        EndIf
      EndIf
    EndDo

    kstep = max(1, maxval(mykstep))

    stop_timer = xnet_wtime()
    timer_xnet = timer_xnet + stop_timer

    ! Output final state
    Call final_output(kstep)

    Return
  End Subroutine full_net

End Module xnet_evolve
