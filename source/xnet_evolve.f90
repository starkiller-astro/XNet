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
    Use nuclear_data, Only: ny, nname, benuc
    Use xnet_abundances, Only: yo, y, yt, ydot, xext, aext, zext
    Use xnet_conditions, Only: t, to, tt, tdel, tdel_old, tdel_next, t9, t9o, t9t, t9dot, rho, rhoo, &
      & rhot, yeo, ye, yet, nt, nto, ntt, tstart, tstop, t9rhofind, t9start, rhostart, yestart
    Use xnet_controls, Only: idiag, iheat, isolv, itsout, kstmx, kmon, ktot, lun_diag, lun_stdout, &
      & lzactive, szbatch, nzbatchmx, nzevolve, zb_lo, zb_hi, zone_id
    Use xnet_integrate, Only: timestep
    Use xnet_integrate_be, Only: solve_be
    Use xnet_integrate_bdf, Only: solve_bdf
    Use xnet_output, Only: final_output, ts_output, write_xnet_th, write_xnet_inab
    Use xnet_timers, Only: xnet_wtime, start_timer, stop_timer, timer_xnet
    Use xnet_types, Only: dp
    Use xnet_util, Only: xnet_terminate
    Implicit None

    ! Output variables
    Integer, Intent(out) :: kstep

    ! Local variables
    !Integer, Parameter :: kstep_output = 10
    Real(dp) :: enm(zb_lo:zb_hi), enb(zb_lo:zb_hi)
    Real(dp) :: enold(zb_lo:zb_hi), en0(zb_lo:zb_hi)
    Real(dp) :: delta_en(zb_lo:zb_hi), edot(zb_lo:zb_hi)
    Real(dp) :: yout(ny+1)
    Integer :: k, izb, izone, nstep_est, idiag0
    Integer :: its(zb_lo:zb_hi), mykstep(zb_lo:zb_hi)
    Integer :: ierr
    Logical :: lzsolve(zb_lo:zb_hi), lzoutput(zb_lo:zb_hi)
    Character(5) :: nname_out(ny+1)
    Character(128) :: ab_fname, th_fname

    start_timer = xnet_wtime()
    timer_xnet = timer_xnet - start_timer

    ! Initialize counters
    kstep = 0
    Do izb = zb_lo, zb_hi
      kmon(:,izb) = 0
      ktot(:,izb) = 0
    EndDo

    ! Set reaction controls not read in from control
    idiag0 = idiag

    ! Initialize trial time step abundances and conditions
    Call t9rhofind(0,t(zb_lo:zb_hi),nt(zb_lo:zb_hi),t9(zb_lo:zb_hi),rho(zb_lo:zb_hi))
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
      yet(izb) = ye(izb)
      yeo(izb) = ye(izb)
      yo(:,izb) = y(:,izb)
      yt(:,izb) = y(:,izb)
    EndDo

    en0 = 0.0
    enm = 0.0
    delta_en = 0.0
    edot = 0.0
    Do izb = zb_lo, zb_hi
      If ( lzactive(izb) ) Then

        ! Calculate the total energy of the nuclei
        Call benuc(yt(:,izb),enb(izb),enm(izb))
        en0(izb) = enm(izb)
        delta_en(izb) = 0.0
        edot(izb) = 0.0
      EndIf
    EndDo

    ! Output initial abundances and conditions
    Call ts_output(0,delta_en,edot)

    ! Start evolution
    Do izb = zb_lo, zb_hi
      If ( lzactive(izb) ) Then
        its(izb) = 0
      Else
        its(izb) = -1
      EndIf
      lzsolve(izb) = lzactive(izb)
      mykstep(izb) = 0
    EndDo
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

      Do izb = zb_lo, zb_hi
        izone = izb + szbatch - zb_lo

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
          Call benuc(yt(:,izb),enb(izb),enm(izb))
          delta_en(izb) = enm(izb) - en0(izb)
          edot(izb) = -(enm(izb)-enold(izb)) / tdel(izb)

        ! If reduced timesteps fail to successfully integrate, warn and flag to remove from loop
        ElseIf ( its(izb) == 1 ) Then
          If ( idiag >= 0 ) Write(lun_diag,"(a,2(i5,a))") 'KStep ',kstep,' Zone ',izone,' Inter!!!'
          its(izb) = 2
        EndIf
      EndDo

      lzoutput = ( lzsolve .and. its <= 0 )
      Call ts_output(kstep,delta_en,edot,mask_in = lzoutput)

      Do izb = zb_lo, zb_hi
        izone = izb + szbatch - zb_lo

          ! If this zone reaches the stop time, flag it to remove from loop
          If ( t(izb) >= tstop(izb) ) Then
            mykstep(izb) = kstep
            its(izb) = -1
          EndIf
      EndDo

      ! Test if all zones have stopped
      lzsolve = ( its == 0 )
      If ( .not. any( lzsolve ) ) Exit
    EndDo

    ! Test that the stop time is reached
    Do izb = zb_lo, zb_hi
      If ( lzactive(izb) ) Then
        If ( t(izb) < tstop(izb) .or. its(izb) > 0 ) Then
          izone = izb + szbatch - zb_lo
          Write(lun_stdout,"(a,i5,a,3(es12.4,a))") &
            & 'Zone ',izone,' Evolution stopped at time= ',t(izb),' with timestep= ',tdel(izb), &
            & ', stop time= ',tstop(izb),' not reached!'
          nstep_est = int( min((tstop(izb)-t(izb))/tdel(izb), real(huge(nstep_est),dp) ) )
          Write(lun_stdout,"(a,i12,a)") &
            & 'Approximately ',nstep_est,' more steps needed'

          nname_out(1:ny) = nname(1:ny)
          nname_out(ny+1) = '  Aux'
          yout(1:ny) = yo(:,izb) 
          yout(ny+1) = xext(izb)
          Write(ab_fname,"(3(a,i3.3))") 'ab_fail_',zone_id(1,izb),'_',zone_id(2,izb),'_',zone_id(3,izb)
          Call write_xnet_inab(ab_fname,ab_fname,nname_out,yout,ierr)
          Write(th_fname,"(3(a,i3.3))") 'th_fail_',zone_id(1,izb),'_',zone_id(2,izb),'_',zone_id(3,izb)
          Call write_xnet_th(th_fname,th_fname,tstart(izb),tstop(izb),t9start(izb),rhostart(izb),yestart(izb),ierr)

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
