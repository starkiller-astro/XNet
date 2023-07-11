!***************************************************************************************************
! xnet_output.f90 10/18/17
! This file contains the ts_output routine, which controls the information output at the end of each
! timestep, and the final_output routine, which controls the output at the end of the evolution.
!***************************************************************************************************

Module xnet_output
  !-------------------------------------------------------------------------------------------------
  ! This module contains data and routines used for output.
  !-------------------------------------------------------------------------------------------------
  Implicit None

  interface write_xnet_th
    module procedure write_xnet_th_history
    module procedure write_xnet_th_point
  end interface write_xnet_th

Contains

  Subroutine ts_output(kstep,enuc,edot,mask_in)
    !-----------------------------------------------------------------------------------------------
    ! The per timestep output routine. If the flag itsout > 0, full_net calls this routine to handle
    ! stepwise output.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny, aa, nname
    Use xnet_abundances, Only: y, yo, ydot
    Use xnet_conditions, Only: rho, t, t9, tdel, t9o, t9dot
    Use xnet_controls, Only: nnucout, nnucout_string, idiag, inucout, itsout, kmon, lun_diag, &
      & lun_ev, lun_stdout, lun_ts, szbatch, zb_lo, zb_hi, lzactive, itrain, lun_train
    Use xnet_flux, Only: flx, flx_int, flux
    Use xnet_match, Only: iwflx, mflx, nflx
    Use xnet_timers, Only: xnet_wtime, start_timer, stop_timer, timer_output
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep
    Real(dp), Intent(in) :: enuc(zb_lo:zb_hi), edot(zb_lo:zb_hi)

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Character(40) :: ev_format
    Integer :: i, j, k, izb, izone
    Logical, Pointer :: mask(:)

    start_timer = xnet_wtime()
    timer_output = timer_output - start_timer

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf
    If ( .not. any(mask) ) Return

    ! Calculate reaction fluxes
    If ( itsout >= 1 .or. idiag >= 1 ) Then
      If ( kstep > 0 ) Then
        Call flux(mask_in = mask)
      Else
        flx_int = 0.0
      EndIf
    EndIf

    If ( itrain >= 1 ) Then
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          izone = izb + szbatch - zb_lo

          ! Snapshot of potential training data for integrator
          Write(lun_train(izb)) kstep,t(izb),t9(izb),rho(izb),tdel(izb),edot(izb),y(:,izb), &
            & yo(:,izb),ydot(:,izb),t9(izb),t9o(izb),t9dot(izb)

        EndIf
      EndDo
    EndIf

    If ( itsout >= 1 ) Then
      Write(ev_format,"(a)") "(i4,1es15.8,2es10.3,2es10.2,"//trim(nnucout_string)//"es9.2,2i2)"
      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          izone = izb + szbatch - zb_lo

          ! An abundance snapshot is written to the binary file
          Write(lun_ts(izb)) kstep,t(izb),t9(izb),rho(izb),tdel(izb),edot(izb),y(:,izb),flx(:,izb)

          ! For itsout>=2, output important mass fractions to the ASCII file
          If ( itsout >= 2 ) Write(lun_ev(izb),ev_format) &
            & kstep,t(izb),t9(izb),rho(izb),edot(izb),tdel(izb), &
            & (aa(inucout(i))*y(inucout(i),izb),i=1,nnucout),(kmon(j,izb),j=1,2)

          ! For itsout>=3, output time and thermo evolution to the screen
          If ( itsout >= 3 ) Write(lun_stdout,"(3i5,4es12.4,2i3)") &
            & izb, izone, kstep,t(izb),tdel(izb),t9(izb),rho(izb),(kmon(j,izb),j=1,2)

          ! For itsout>=4, output abundances for each timestep to the diagnostic file
          If ( itsout >= 4 .and. idiag >= 0 ) Write(lun_diag,"(4(a5,es14.7,1x))") &
            & (nname(i),aa(i)*y(i,izb),i=1,ny)

          ! For itsout>=5, output fluxes for each timestep to the diagnostic file
          If ( itsout >= 5 .and. idiag >= 0 ) Write(lun_diag,"(i5,10a5,i5,es11.3)") &
            & (k,(nname(nflx(j,k)),j=1,4),' <-> ',(nname(nflx(j,k)),j=5,8),iwflx(k),flx(k,izb),k=1,mflx)
        EndIf
      EndDo
    EndIf

    stop_timer = xnet_wtime()
    timer_output = timer_output + stop_timer

    Return
  End Subroutine ts_output

  Subroutine final_output(kstep,mask_in)
    !-------------------------------------------------------------------------------------------------
    ! full_net calls this routine to handle output at the end of its execution.
    !-------------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny, nname, aa
    Use xnet_abundances, Only: y
    Use xnet_conditions, Only: tt, tstop
    Use xnet_controls, Only: changemx, iconvc, isolv, kitmx, kstmx, ktot, lun_diag, tolc, tolm, yacc, &
      & szbatch, zb_lo, zb_hi, lzactive, idiag
    Use xnet_flux, Only: flx_int
    Use xnet_match, Only: iwflx, mflx, nflx
    Use xnet_timers
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(zb_lo:zb_hi)

    ! Local variables
    Integer, Parameter :: itimer_reset = 0
    Integer :: iflx_mx
    Integer :: i, k, izb, izone
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask(zb_lo:) => mask_in
    Else
      mask(zb_lo:) => lzactive(zb_lo:zb_hi)
    EndIf

    If ( idiag >= 0 ) Then

      start_timer = xnet_wtime()
      timer_output = timer_output - start_timer

      Do izb = zb_lo, zb_hi
        If ( mask(izb) ) Then
          izone = izb + szbatch - zb_lo

          ! Write final abundances to diagnostic output (in ASCII)
          Write(lun_diag,"(a3,3i6,2es14.7)") 'End',izone,kstep,kstmx,tt(izb),tstop(izb)
          Write(lun_diag,"(4(a5,es14.7,1x))") (nname(i),aa(i)*y(i,izb),i=1,ny)
          Write(lun_diag,"(a3,3i6,4es10.3)") 'CNT',isolv,iconvc,kitmx,yacc,tolm,tolc,changemx

          ! Write integrated fluxes to diagnotic output
          If ( idiag >= 1 ) Then
            iflx_mx = maxloc(abs(flx_int(:,izb)),dim=1)
            Write(lun_diag,"(2x,a8,i5,es11.3,9a5)") 'Flux Max',&
              & iflx_mx,flx_int(iflx_mx,izb),(nname(nflx(i,iflx_mx)),i=1,4),' <-> ',(nname(nflx(i,iflx_mx)),i=5,8)
            Write(lun_diag,"(i5,9a5,i5,es11.3)") &
              & (k,(nname(nflx(i,k)),i=1,4),' <-> ',(nname(nflx(i,k)),i=5,8),iwflx(k),flx_int(k,izb),k=1,mflx)
          EndIf

          Write(lun_diag,"(a10,a5,5a10)") 'Counters: ','Zone','TS','NR','Jacobian','Deriv','CrossSect'
          Write(lun_diag,"(10x,i5,5i10)") izone,(ktot(i,izb),i=1,5)
        EndIf
      EndDo

      stop_timer = xnet_wtime()
      timer_output = timer_output + stop_timer

      ! Write timers
      If ( idiag >= 0 ) Then
        Write(lun_diag,"(a8,10a10)") &
          & 'Timers: ','Total','TS','NR','Solver','Jacobian','Deriv','CrossSect','Screening','Setup','Output'
        Write(lun_diag,"(8x,10es10.3)") &
          & timer_xnet,timer_tstep,timer_nraph,timer_solve,timer_jacob,timer_deriv,timer_csect,timer_scrn,timer_setup,timer_output
      EndIf

      ! Use the following line to restart the timers
      If ( itimer_reset > 0 ) Call reset_timers

    EndIf

    Return
  End Subroutine final_output

  Subroutine write_xnet_inab( fname, header, nname, y, ierr )
    Use xnet_controls, Only: lun_stderr
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: fname
    Character(*), Intent(in) :: header
    Character(5), Intent(in) :: nname(:)
    real(dp), Intent(in) :: y(:)

    ! Output variables
    Integer, Intent(out) :: ierr

    ! Local variables
    Integer :: lun_inab
    Integer :: i

    ierr = 0

    ! open file
    open (newunit=lun_inab,file=trim(fname),iostat=ierr)
    if ( ierr /= 0 ) then
      write (lun_stderr,'(3a,i3,a)') '[write_xnet_inab] ERROR: Problem opening file: ',trim(fname),'(',ierr,')'
      return
    end if

    ! write header
    write (lun_inab,'(a)',iostat=ierr) trim(header)
    if ( ierr /= 0 ) then
      write (lun_stderr,'(a,i3,a)') '[write_xnet_inab] ERROR: Problem writing header (',ierr,')'
      close (lun_inab)
      return
    end if

    ! write abundances
    write (lun_inab,'(a5,es23.15)',iostat=ierr) (adjustr(nname(i)), y(i), i=1,size(y))
    if ( ierr /= 0 ) then
      write (lun_stderr,'(a,i3,a)') '[write_xnet_inab] ERROR: Problem writing data (',ierr,')'
      close (lun_inab)
      return
    end if

    ! close file
    close (lun_inab)

    return
  End Subroutine write_xnet_inab

  Subroutine write_xnet_th_history( fname, header, time, t9, rho, ye, ierr )
    Use xnet_controls, Only: lun_stderr
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: fname
    Character(*), Intent(in) :: header
    Real(dp), Intent(in) :: time(:) ! First three values are for header
    Real(dp), Intent(in) :: t9(size(time)-3)
    Real(dp), Intent(in) :: rho(size(time)-3)
    Real(dp), Intent(in) :: ye(size(time)-3)

    ! Output variables
    Integer, Intent(out) :: ierr

    ! Local variables
    Integer :: lun_th
    Integer :: tmp_ierr(4)
    Integer :: i

    ierr = 0
    tmp_ierr = 0

    if ( size(time) < 5 ) then
      write (lun_stderr,'(a)') '[write_xnet_th] ERROR: Not enough data points'
      ierr = 1
      return
    end if

    ! open file
    open (newunit=lun_th,file=trim(fname),iostat=ierr)
    if ( ierr /= 0 ) then
      write (lun_stderr,'(3a,i3,a)') '[write_xnet_th] ERROR: Problem opening file: ',trim(fname),'(',ierr,')'
      return
    end if

    ! write header
    write (lun_th,'(a)',iostat=tmp_ierr(1)) trim(header)
    write (lun_th,'(es23.15,a)',iostat=tmp_ierr(2)) time(1),' start time'
    write (lun_th,'(es23.15,a)',iostat=tmp_ierr(3)) time(2),' stop time'
    write (lun_th,'(es23.15,a)',iostat=tmp_ierr(4)) time(3),' initial timestep'
    if ( any( tmp_ierr /= 0 ) ) then
      ierr = maxval(tmp_ierr)
      write (lun_stderr,'(a,i3,a)') '[write_xnet_th] ERROR: Problem writing header (',ierr,')'
      close (lun_th)
      return
    end if

    write (lun_th,'(4es23.15)',iostat=ierr) (time(i+3),t9(i),rho(i),ye(i),i=1,size(t9))
    if ( ierr /= 0 ) then
      write (lun_stderr,'(a,i3,a)') '[write_xnet_th] ERROR: Problem writing data (',ierr,')'
      close (lun_th)
      return
    end if

    close (lun_th)

    return
  End Subroutine write_xnet_th_history

  Subroutine write_xnet_th_point( fname, header, tstart, tstop, t9, rho, ye, ierr )
    Use xnet_controls, Only: lun_stderr
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: fname
    Character(*), Intent(in) :: header
    Real(dp), Intent(in) :: tstart
    Real(dp), Intent(in) :: tstop
    Real(dp), Intent(in) :: t9
    Real(dp), Intent(in) :: rho
    Real(dp), Intent(in) :: ye

    ! Output variables
    Integer, Intent(out) :: ierr

    ! Local variables
    Integer :: lun_th
    Integer :: tmp_ierr(4)
    Integer :: i

    ierr = 0
    tmp_ierr = 0

    ! open file
    open (newunit=lun_th,file=trim(fname),iostat=ierr)
    if ( ierr /= 0 ) then
      write (lun_stderr,'(3a,i3,a)') '[write_xnet_th] ERROR: Problem opening file: ',trim(fname),'(',ierr,')'
      return
    end if

    ! write header
    write (lun_th,'(a)',iostat=tmp_ierr(1)) trim(header)
    write (lun_th,'(es23.15,a)',iostat=tmp_ierr(2)) tstart,' start time'
    write (lun_th,'(es23.15,a)',iostat=tmp_ierr(3)) tstop,' stop time'
    write (lun_th,'(es23.15,a)',iostat=tmp_ierr(4)) (tstop-tstart)*1.0e-4_dp,' initial timestep'
    if ( any( tmp_ierr /= 0 ) ) then
      ierr = maxval(tmp_ierr)
      write (lun_stderr,'(a,i3,a)') '[write_xnet_th] ERROR: Problem writing header (',ierr,')'
      close (lun_th)
      return
    end if

    write (lun_th,'(4es23.15)',iostat=ierr) tstop, t9, rho, ye
    write (lun_th,'(4es23.15)',iostat=ierr) tstop, t9, rho, ye
    if ( ierr /= 0 ) then
      write (lun_stderr,'(a,i3,a)') '[write_xnet_th] ERROR: Problem writing data (',ierr,')'
      close (lun_th)
      return
    end if

    flush (lun_th)
    close (lun_th)

    return
  End Subroutine write_xnet_th_point

End Module xnet_output
