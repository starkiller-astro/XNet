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

  Subroutine full_net
    !-----------------------------------------------------------------------------------------------
    ! The abundance evolution is performed over a series of timesteps, with the duration of the
    ! timestep determined by the integration scheme and the changing thermodynamic conditions.
    ! Integration is performed by a choice of methods controlled by the isolv flag.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny, nname, benuc
    Use xnet_abundances, Only: yo, y, yt, ydot
    Use xnet_conditions, Only: t, to, tt, tdel, tdel_old, tdel_next, t9, t9o, t9t, t9dot, rho, rhoo, &
      rhot, yeo, ye, yet, nt, nto, ntt, tstart, tstop, t9rhofind
    Use xnet_controls, Only: idiag, iheat, isolv, itsout, kstmx, kmon, ktot, lun_diag, lun_stdout, &
      lzactive, szbatch, nzbatchmx
    Use xnet_integrate, Only: timestep
    Use xnet_integrate_be, Only: solve_be
    Use xnet_integrate_bdf, Only: solve_bdf
    Use xnet_output, Only: final_output, ts_output
    Use xnet_timers, Only: xnet_wtime, start_timer, stop_timer, timer_xnet
    Use xnet_types, Only: dp
    Use xnet_util, Only: xnet_terminate
    Implicit None

    ! Local variables
    !Integer, Parameter :: kstep_output = 10
    Real(dp) :: enm(nzbatchmx), enb(nzbatchmx), enold(nzbatchmx), en0(nzbatchmx)
    Real(dp) :: delta_en(nzbatchmx), edot(nzbatchmx)
    Real(dp) :: ytot, ztot, atot
    Integer :: idiag0, its(nzbatchmx), mykstep(nzbatchmx)
    Integer :: k, izb, izone, kstep, nstep_est
    Logical :: lzstep(nzbatchmx)

    start_timer = xnet_wtime()
    timer_xnet = timer_xnet - start_timer

    ! Initialize counters
    kstep = 0
    kmon(:,:) = 0
    ktot(:,:) = 0

    ! Set reaction controls not read in from control
    idiag0 = idiag

    ! Initialize trial time step abundances and conditions
    Call t9rhofind(0,t,nt,t9,rho)
    nto(:) = nt(:)
    ntt(:) = nt(:)
    to(:) = t(:)
    tt(:) = t(:)
    yo(:,:) = y(:,:)
    yt(:,:) = y(:,:)
    t9o(:) = t9(:)
    t9t(:) = t9(:)
    rhoo(:) = rho(:)
    rhot(:) = rho(:)
    yet(:) = ye(:)
    yeo(:) = ye(:)
    tdel_old(:) = tdel(:)
    tdel_next(:) = tdel(:)

    edot(:) = 0.0
    en0(:) = 0.0
    enm(:) = 0.0
    delta_en(:) = 0.0
    Do izb = 1, nzbatchmx
      If ( lzactive(izb) ) Then

        ! Calculate the total energy of the nuclei
        Call benuc(yt(:,izb),enb(izb),enm(izb),ytot,ztot,atot)
        en0(izb) = enm(izb)
        edot(izb) = 0.0

        If ( itsout > 0 ) Write(lun_stdout,"(a,i6,a,i2,2(a,es10.3))") &
          & 'Max Step',kstmx,' IDiag=',idiag,' Start Time',tstart(izb),' Stop Time',tstop(izb)
      EndIf
    EndDo

    ! Output initial abundances and conditions
    delta_en(:) = enm(:) - en0(:)
    Call ts_output(kstep,delta_en,edot)

    ! Start evolution
    Where ( lzactive(:) )
      its(:) = 0
    ElseWhere
      its(:) = -1
    EndWhere
    mykstep(:) = 0
    lzstep(:) = ( its(:) < 0 )
    Do kstep = 1, kstmx

      ! Determine if this is an output step
      idiag = idiag0
      !If ( mod(kstep,kstep_output) == 0 ) idiag = 2

      ! Calculate an initial guess for the timestep
      Call timestep(kstep,mask_in = (its(:) == 0))

      ! Take integration step (only worry about solve_be for now)
      Select Case (isolv)
      Case (3)
        Call solve_bdf(kstep,its)
      !Case (2)
      !  Call solve_bd(kstep,its)
      Case Default
        Call solve_be(kstep,its)
      End Select

      Do izb = 1, nzbatchmx
        izone = izb + szbatch - 1

        ! If convergence is successful, output timestep results
        If ( its(izb) == 0 ) Then
          If ( idiag >= 1 ) Write(lun_diag,"(2(a,i5),5es14.7)") &
            & 'KStep ',kstep,' Zone ',izone,t(izb),tdel(izb),t9(izb),rho(izb),ye(izb)
          If ( idiag >= 3 ) Then
            Write(lun_diag,"(a8)") 'delta Y'
            Write(lun_diag,"(8x,a5,4es12.4)") &
              & (nname(k),y(k,izb),yo(k,izb),(y(k,izb)-yo(k,izb)),(tdel(izb)*ydot(k,izb)),k=1,ny)
            If ( iheat > 0 ) Write(lun_diag,"(8x,a5,4es12.4)") &
              & 'T9',t9(izb),t9o(izb),t9(izb)-t9o(izb),tdel(izb)*t9dot(izb)
          EndIf
          enold(izb) = enm(izb)
          Call benuc(yt(:,izb),enb(izb),enm(izb),ytot,ztot,atot)
          edot(izb) = -(enm(izb)-enold(izb)) / tdel(izb)

        ! If reduced timesteps fail to successfully integrate, warn and flag to remove from loop
        ElseIf ( its(izb) == 1 ) Then
          If ( idiag >= 0 ) Write(lun_diag,"(a,2(i5,a))") 'KStep ',kstep,' Zone ',izone,' Inter!!!'
          its(izb) = 2
        EndIf
      EndDo
      delta_en(:) = enm(:) - en0(:)
      Call ts_output(kstep,delta_en,edot,mask_in = (its(:) == 0))

      ! If this zone reaches the stop time, flag it to remove from loop
      Where ( t(:) >= tstop(:) .and. its(:) == 0 )
        mykstep(:) = kstep
        its(:) = -1
      EndWhere

      ! Test if all zones have stopped
      If ( all( its(:) /= 0 ) ) Exit
    EndDo

    ! Test that the stop time is reached
    Do izb = 1, nzbatchmx
      If ( lzactive(izb) ) Then
        If ( t(izb) < tstop(izb) .or. its(izb) > 0 ) Then
          izone = izb + szbatch - 1
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

    stop_timer = xnet_wtime()
    timer_xnet = timer_xnet + stop_timer

    ! Output final state
    Call final_output(kstep)

    Return
  End Subroutine full_net

End Module xnet_evolve
