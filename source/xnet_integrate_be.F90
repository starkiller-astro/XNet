!***************************************************************************************************
! xnet_integrate_be.f90 10/18/17
! Backward Euler solver
!
! The routines in this file perform the Backward Euler time integration for the thermonuclear
! reaction network.
!***************************************************************************************************

#include "xnet_macros.fh"

Module xnet_integrate_be
  Implicit None

Contains

  Subroutine solve_be(kstep,its)
    !-----------------------------------------------------------------------------------------------
    ! This routine handles the potential failure of an individual Backward Euler step. Repeated calls
    ! to the BE step integrator are made with successivley smaller timesteps until success is achieved
    ! or ktsmx trials fail. The abundances are evolved using a Newton-Raphson iteration scheme to
    ! solve the equation (yt(i)-y(i))/tdel = ydot(i), where y(i) is the abundance at the beginning of
    ! the iteration, yt(i) is the trial abundance at the end of the timestep, ydot(i) is the time
    ! derivative of the trial abundance calculated from reaction rates and tdel is the timestep.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny
    Use reaction_data, Only: b1, b2, b3, b4, csect1, csect2, csect3, csect4, &
      & dcsect1dt9, dcsect2dt9, dcsect3dt9, dcsect4dt9
    Use xnet_abundances, Only: y, yo, yt, ydot
    Use xnet_conditions, Only: t, to, tt, tdel, tdel_next, tdelstart, t9, t9o, t9t, t9dot, rho, rhoo, &
      & rhot, yeo, ye, yet, nt, nto, ntt, t9rhofind, cv, etae, detaedt9
    Use xnet_controls, Only: idiag, iheat, kitmx, kmon, ktot, lun_diag, lun_stdout, tdel_maxmult, &
      & szbatch, zb_lo, zb_hi
    Use xnet_integrate, Only: timestep
    Use xnet_screening, Only: h1, h2, h3, h4, dh1dt9, dh2dt9, dh3dt9, dh4dt9
    Use xnet_timers, Only: xnet_wtime, start_timer, stop_timer, timer_tstep
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep

    ! Input/Output variables
    Integer, Intent(inout) :: its(zb_lo:zb_hi) ! On input,  = 0 indicates active zone
                                               ! On output, = 1 if zone fails to converge

    ! Local variables
    Integer, Parameter :: ktsmx = 10
    Integer :: kts, k, izb, izone
    Integer :: inr(zb_lo:zb_hi)
    Integer :: mykts(zb_lo:zb_hi)
    Logical :: lzstep(zb_lo:zb_hi)

    start_timer = xnet_wtime()
    timer_tstep = timer_tstep - start_timer

    !__dir_enter_data &
    !__dir_async &
    !__dir_create(b1,b2,b3,b4,csect1,csect2,csect3,csect4) &
    !__dir_create(dcsect1dt9,dcsect2dt9,dcsect3dt9,dcsect4dt9) &
    !__dir_create(h1,h2,h3,h4,dh1dt9,dh2dt9,dh3dt9,dh4dt9) &
    !__dir_copyin(t,to,tt,tdel,tdel_next,tdelstart,nto,nt,ntt) &
    !__dir_copyin(t9o,t9,t9t,t9dot,rhoo,rho,rhot,yo,y,yt,ydot) &
    !__dir_copyin(yeo,ye,yet,cv,etae,detaedt9) &
    !__dir_copyin(kmon,ktot) &
    !__dir_create(inr,mykts,lzstep) &
    !__dir_copyin(its)

    ! If the zone has previously converged or failed, do not iterate
    !__dir_loop_outer(1) &
    !__dir_async &
    !__dir_present(its,inr,lzstep,mykts)
    Do izb = zb_lo, zb_hi
      If ( its(izb) /= 0 ) Then
        inr(izb) = -1
        lzstep(izb) = .false.
      Else
        inr(izb) = 0
        lzstep(izb) = .true.
      EndIf
      mykts(izb) = 0
    EndDo

    !-----------------------------------------------------------------------------------------------
    ! For each trial timestep, tdel, the NR iteration is attempted.
    ! The possible results for a timestep are a converged set of yt or a failure to converge.
    ! If convergence fails, retry with the trial timestep reduced by tdel_maxmult.
    !-----------------------------------------------------------------------------------------------
    Do kts = 1, ktsmx

      ! Attempt Backward Euler integration over desired timestep
      Call step_be(kstep,inr)

      !__dir_loop_outer(1) &
      !__dir_async &
      !__dir_present(its,inr,tdel,tt,t,yet,ye,yt,y,mykts,kmon,ktot)
      Do izb = zb_lo, zb_hi

        ! If integration fails, reset abundances, reduce timestep and retry.
        If ( inr(izb) == 0 ) Then
          tdel(izb) = tdel(izb) / tdel_maxmult
          tt(izb) = t(izb) + tdel(izb)
          yet(izb) = ye(izb)
          mykts(izb) = kts+1

          !__dir_loop_inner(1)
          Do k = 1, ny
            yt(k,izb) = y(k,izb)
          EndDo

          ! Record number of NR iterations
          If ( its(izb) == 0 ) Then
            kmon(2,izb) = kitmx + 1
            ktot(2,izb) = ktot(2,izb) + kitmx + 1
          EndIf

        ! If integration is successfull, flag the zone for removal from zone loop
        ElseIf ( inr(izb) > 0 ) Then
          mykts(izb) = kts

          ! Record number of NR iterations
          If ( its(izb) == 0 ) Then
            kmon(2,izb) = inr(izb)
            ktot(2,izb) = ktot(2,izb) + inr(izb)
          EndIf
        EndIf
        lzstep(izb) = ( inr(izb) == 0 )
      EndDo
      !__dir_update &
      !__dir_wait &
      !__dir_host(lzstep)

      ! Reset temperature and density for failed integrations
      Call t9rhofind(kstep,tt(zb_lo:zb_hi),ntt(zb_lo:zb_hi), &
        & t9t(zb_lo:zb_hi),rhot(zb_lo:zb_hi),mask_in = lzstep)
      If ( iheat > 0 ) Then
        !__dir_loop_outer(1) &
        !__dir_async &
        !__dir_present(lzstep,t9t,t9)
        Do izb = zb_lo, zb_hi
          If ( lzstep(izb) ) Then
            t9t(izb) = t9(izb)
          EndIf
        EndDo
      EndIf

      ! For the last attempt, re-calculate timestep based on derivatives
      ! as is done for the first timestep.
      If ( kts == ktsmx-1 ) Then
        !__dir_loop_outer(1) &
        !__dir_async &
        !__dir_present(lzstep,tdel,tdelstart)
        Do izb = zb_lo, zb_hi
          If ( lzstep(izb) ) Then
            tdel(izb) = 0.0
            tdelstart(izb) = 0.0
          EndIf
        EndDo
        Call timestep(kstep,mask_in = lzstep)
      EndIf

      ! Log the failed integration attempts
      If ( idiag >= 2 ) Then
        !__dir_update &
        !__dir_wait &
        !__dir_host(inr,tt,tdel)
        Do izb = zb_lo, zb_hi
          If ( inr(izb) == 0 ) Then
            izone = izb + szbatch - zb_lo
            Write(lun_diag,"(a,i5,i3,3es12.4)") 'BE TS Reduce',izone,kts,tt(izb),tdel(izb),tdel_maxmult
          EndIf
        EndDo
      EndIf

      ! Test if all zones have converged
      If ( .not. any( lzstep ) ) Exit
    EndDo

    ! Mark TS convergence only for zones which haven't previously failed or converged
    !__dir_loop_outer(1) &
    !__dir_async &
    !__dir_present(its,inr,kmon,ktot,mykts,tdel,tdel_next,nt,nto,ntt,t,to,tt) &
    !__dir_present(t9,t9o,t9t,rho,rhoo,rhot,ye,yeo,yet,y,yo,yt)
    Do izb = zb_lo, zb_hi
      If ( inr(izb) >= 0 ) Then
        kmon(1,izb) = mykts(izb)
        ktot(1,izb) = ktot(1,izb) + mykts(izb)
        If ( inr(izb) > 0 ) Then
          tdel_next(izb) = tdel(izb) * tdel_maxmult
          nto(izb) = nt(izb)
          nt(izb) = ntt(izb)
          to(izb) = t(izb)
          t(izb) = tt(izb)
          t9o(izb) = t9(izb)
          t9(izb) = t9t(izb)
          rhoo(izb) = rho(izb)
          rho(izb) = rhot(izb)
          yeo(izb) = ye(izb)
          ye(izb) = yet(izb)
          Do k = 1, ny
            yo(k,izb) = y(k,izb)
            y(k,izb) = yt(k,izb)
          EndDo
        Else
          its(izb) = 1
        EndIf
      EndIf
    EndDo

    ! Log TS success/failure
    If ( idiag >= 0 ) Then
      !__dir_update &
      !__dir_wait &
      !__dir_host(inr,mykts,t,tdel,t9t,rhot)
      Do izb = zb_lo, zb_hi
        izone = izb + szbatch - zb_lo
        If ( inr(izb) > 0 .and. idiag >= 2 ) Then
          Write(lun_diag,"(a,2i5,2i3)") &
            & 'BE TS Success',kstep,izone,mykts(izb),inr(izb)
        ElseIf ( inr(izb) == 0 ) Then
          Write(lun_diag,"(a,2i5,4es12.4,2i3)") &
            & 'BE TS Fail',kstep,izone,t(izb),tdel(izb),t9t(izb),rhot(izb),inr(izb),mykts(izb)
          Write(lun_stdout,*) 'Timestep retrys fail after ',mykts(izb),' attempts'
        EndIf
      EndDo
    EndIf

    !__dir_exit_data &
    !__dir_async &
    !__dir_delete(b1,b2,b3,b4,csect1,csect2,csect3,csect4) &
    !__dir_delete(dcsect1dt9,dcsect2dt9,dcsect3dt9,dcsect4dt9) &
    !__dir_delete(h1,h2,h3,h4,dh1dt9,dh2dt9,dh3dt9,dh4dt9) &
    !__dir_copyout(t,to,tt,tdel,tdel_next,tdelstart,nto,nt,ntt) &
    !__dir_copyout(t9o,t9,t9t,t9dot,rhoo,rho,rhot,yo,y,yt,ydot) &
    !__dir_copyout(yeo,ye,yet,cv,etae,detaedt9) &
    !__dir_copyout(kmon,ktot) &
    !__dir_delete(inr,mykts,lzstep) &
    !__dir_copyout(its)

    !__dir_wait

    stop_timer = xnet_wtime()
    timer_tstep = timer_tstep + stop_timer

    Return
  End Subroutine solve_be

  Subroutine step_be(kstep,inr)
    !-----------------------------------------------------------------------------------------------
    ! This routine attempts to integrate a single Backward Euler step for the timestep tdel.
    ! If successful, inr = 1
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: ny, aa, nname
    Use xnet_abundances, Only: y, ydot, yt, xext
    Use xnet_conditions, Only: t9, t9dot, t9t, tdel, nh
    Use xnet_controls, Only: iconvc, idiag, iheat, ijac, kitmx, lun_diag, tolc, tolm, tolt9, ymin, &
      & szbatch, zb_lo, zb_hi
    Use xnet_integrate, Only: cross_sect, yderiv
    Use xnet_jacobian, Only: jacobian_bksub, jacobian_decomp, jacobian_build, jacobian_solve
    Use xnet_timers, Only: xnet_wtime, start_timer, stop_timer, timer_nraph
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep

    ! Input/Output variables
    Integer, Intent(inout) :: inr(zb_lo:zb_hi) ! On input,  = 0 indicates active zone
                                               !            =-1 indicates inactive zone
                                               ! On output, > 0 indicates # NR iterations if converged

    ! Local variables
    Integer :: irdymx, idymx
    Integer :: i, k, kit, izb, izone
    Real(dp) :: s1, s2, s3, ytmp
    Real(dp) :: testc(zb_lo:zb_hi), testc2(zb_lo:zb_hi), testm(zb_lo:zb_hi)
    Real(dp) :: testn(zb_lo:zb_hi), toln(zb_lo:zb_hi)
    Real(dp) :: xtot(zb_lo:zb_hi), xtot_init(zb_lo:zb_hi), rdt(zb_lo:zb_hi), mult(zb_lo:zb_hi)
    Real(dp) :: yrhs(ny,zb_lo:zb_hi), dy(ny,zb_lo:zb_hi), reldy(ny,zb_lo:zb_hi)
    Real(dp) :: t9rhs(zb_lo:zb_hi), dt9(zb_lo:zb_hi), relt9(zb_lo:zb_hi)
    Logical :: iterate(zb_lo:zb_hi), eval_rates(zb_lo:zb_hi), rebuild(zb_lo:zb_hi)

    start_timer = xnet_wtime()
    timer_nraph = timer_nraph - start_timer

    !__dir_enter_data &
    !__dir_async &
    !__dir_create(iterate,eval_rates,rebuild,testc,testc2,testm,testn,toln) &
    !__dir_create(xtot,xtot_init,rdt,mult,yrhs,dy,reldy,t9rhs,dt9,relt9)

    !__dir_loop_outer(1) &
    !__dir_async &
    !__dir_present(inr,iterate,xtot_init,rdt,mult,aa,y,tdel,toln,xext) &
    !__dir_private(s1)
    Do izb = zb_lo, zb_hi
      If ( inr(izb) == 0 ) Then

        ! Create mask for active zone integrations
        iterate(izb) = .true.

        ! Calculate initial total mass fraction
        s1 = 0.0
        !__dir_loop_inner(1) &
        !__dir_reduction(+,s1)
        Do k = 1, ny
          s1 = s1 + aa(k)*y(k,izb)
        EndDo
        xtot_init(izb) = s1 + xext(izb) - 1.0

        ! Jacobian scaling factors
        rdt(izb) = 1.0 / tdel(izb)
        mult(izb) = -1.0
      Else
        iterate(izb) = .false.
        xtot_init(izb) = 0.0
        rdt(izb) = 0.0
        mult(izb) = 0.0
      EndIf
      If ( iconvc /= 0 ) Then
        toln(izb) = tolc
      Else
        toln(izb) = tolm
      EndIf
    EndDo
    !__dir_update &
    !__dir_wait &
    !__dir_host(iterate)

    ! The Newton-Raphson iteration occurs for at most kitmx iterations.
    Do kit = 1, kitmx
      Do izb = zb_lo, zb_hi

        ! Rebuild and update LU factorization of the jacobian every ijac iterations
        rebuild(izb) = ( iterate(izb) .and. mod(kit-1,ijac) == 0 )

        ! If thermodynamic conditions are changing, rates need to be udpated
        If ( iheat > 0 ) Then
          eval_rates(izb) = iterate(izb)
        ElseIf ( kit == 1 ) Then
          eval_rates(izb) = ( iterate(izb) .and. nh(izb) > 1 )
        Else
          eval_rates(izb) = .false.
        EndIf
      EndDo
      !__dir_update &
      !__dir_async &
      !__dir_device(rebuild,eval_rates)

      ! Calculate the reaction rates and abundance time derivatives
      Call cross_sect(mask_in = eval_rates)
      Call yderiv(mask_in = iterate)
      Call jacobian_build(diag_in = rdt,mult_in = mult,mask_in = rebuild)
      Call jacobian_decomp(kstep,mask_in = rebuild)

      ! Calculate equation to zero
      !__dir_loop_outer(1) &
      !__dir_async &
      !__dir_present(iterate,yrhs,y,yt,rdt,ydot,t9rhs,t9,t9t,t9dot)
      Do izb = zb_lo, zb_hi
        If ( iterate(izb) ) Then
          If ( kit > 1 ) Then
            !__dir_loop_inner(1)
            Do k = 1, ny
              yrhs(k,izb) = (y(k,izb)-yt(k,izb))*rdt(izb) + ydot(k,izb)
            EndDo
            If ( iheat > 0 ) t9rhs(izb) = (t9(izb)-t9t(izb))*rdt(izb) + t9dot(izb)
          Else
            !__dir_loop_inner(1)
            Do k = 1, ny
              yrhs(k,izb) = ydot(k,izb)
            EndDo
            If ( iheat > 0 ) t9rhs(izb) = t9dot(izb)
          EndIf
        Else
          !__dir_loop_inner(1)
          Do k = 1, ny
            yrhs(k,izb) = 0.0
          EndDo
          If ( iheat > 0 ) t9rhs(izb) = 0.0
        EndIf
      EndDo
      If ( idiag >= 4 ) Then
        !__dir_update &
        !__dir_wait &
        !__dir_host(yrhs,ydot,yt,t9rhs,t9dot,t9t,rdt)
        Do izb = zb_lo, zb_hi
          If ( iterate(izb) ) Then
            izone = izb + szbatch - zb_lo
            Write(lun_diag,"(a3,2i5,es14.7)") 'RHS',kstep,izone,rdt(izb)
            Write(lun_diag,"(a5,4es23.15)") (nname(i),yrhs(i,izb),ydot(i,izb),yt(i,izb),y(i,izb),i=1,ny)
            If ( iheat > 0 ) Write(lun_diag,"(a5,4es23.15)") 'T9',t9rhs(izb),t9dot(izb),t9t(izb),t9(izb)
          EndIf
        EndDo
      EndIf

      ! Solve the jacobian and calculate the changes in abundances, dy
      !Call jacobian_solve(kstep,yrhs,dy,t9rhs,dt9,mask_in = iterate)
      Call jacobian_bksub(kstep,yrhs,dy,t9rhs,dt9,mask_in = iterate)

      ! Evolve the abundances and calculate convergence tests
      !-----------------------------------------------------------------------------------------
      ! There are 3 included convergence tests:
      ! testc which measures relative changes,
      ! testc2 which measures total abundance changes, and
      ! testm which tests mass conservation.
      !
      ! testc is the most stringent test, and requires the most iterations.
      ! testm is the most lax, and therefore the fastest, often requiring only one iteration.
      ! Considering the uncertainties of the reaction rates, it is doubtful that the
      ! increased precision of testc is truly increased accuracy.
      !-----------------------------------------------------------------------------------------
      !__dir_loop_outer(1) &
      !__dir_async &
      !__dir_present(iterate,yt,dy,reldy,t9t,dt9,relt9,xtot,xtot_init,xext) &
      !__dir_present(aa,inr,testm,testc,testc2,testn,toln) &
      !__dir_private(s1,s2,s3)
      Do izb = zb_lo, zb_hi
        If ( iterate(izb) ) Then
          s1 = 0.0
          s2 = 0.0
          s3 = 0.0
          !__dir_loop_inner(1) &
          !__dir_private(ytmp) &
          !__dir_reduction(+,s1,s2,s3)
          Do k = 1, ny
            ytmp = yt(k,izb) + dy(k,izb)
            If ( ytmp < ymin ) Then
              dy(k,izb) = -yt(k,izb)
              yt(k,izb) = 0.0
              reldy(k,izb) = 0.0
            Else
              yt(k,izb) = ytmp
              reldy(k,izb) = abs(dy(k,izb) / yt(k,izb))
            EndIf
            s1 = s1 + reldy(k,izb)
            s2 = s2 + aa(k)*dy(k,izb)
            s3 = s3 + aa(k)*yt(k,izb)
          EndDo
          testc(izb)  = s1
          testc2(izb) = s2
          xtot(izb)   = s3 + xext(izb) - 1.0
          testm(izb)  = xtot(izb) - xtot_init(izb)

          If ( iheat > 0 ) Then
            t9t(izb) = t9t(izb) + dt9(izb)
            If ( abs(t9t(izb)) > tiny(0.0) ) Then
              relt9(izb) = abs(dt9(izb) / t9t(izb))
            Else
              relt9(izb) = 0.0
            EndIf
          Else
            relt9(izb) = 0.0
          EndIf

          ! Ordinarily, test for true convergence
          If ( iconvc /= 0 ) Then
            testn(izb) = testc(izb)

          ! Otherwise, use mass conservation for convergence condition
          Else
            testn(izb) = testm(izb)
          EndIf

          ! If converged, exit NR loop
          If ( abs(testn(izb)) <= toln(izb) .and. relt9(izb) < tolt9 ) Then
            inr(izb) = kit
            iterate(izb) = .false.
          Else
            iterate(izb) = .true.
          EndIf
        EndIf
      EndDo
      !__dir_update &
      !__dir_wait &
      !__dir_host(iterate)

      If ( idiag >= 3 ) Then
        !__dir_update &
        !__dir_wait &
        !__dir_host(inr,testm,testc,testc2,yt,dy,reldy,t9t,dt9,relt9)
        Do izb = zb_lo, zb_hi
          If ( inr(izb) >= 0 ) Then
            izone = izb + szbatch - zb_lo
            If ( idiag >= 4 ) Then
              irdymx = maxloc(reldy(:,izb),dim=1)
              idymx = maxloc(dy(:,izb),dim=1)
              Write(lun_diag,"(a3,2i5,i3,2(a5,2es23.15))") &
                & ' dY',kstep,izone,kit,nname(idymx),dy(idymx,izb),y(idymx,izb), &
                & nname(irdymx),reldy(irdymx,izb),y(irdymx,izb)
              If ( idiag >= 5 ) Write(lun_diag,"(a5,2es23.15,es12.4,2es23.15)") &
                & (nname(k),yt(k,izb),dy(k,izb),reldy(k,izb), &
                &  (aa(k)*dy(k,izb)),(aa(k)*yt(k,izb)),k=1,ny)
              If ( iheat > 0 ) Write(lun_diag,"(a3,2i5,i3,5x,2es23.15,5x,es12.4)") &
                & 'dT9',kstep,izone,kit,dt9(izb),t9t(izb),relt9(izb)
            EndIf
            Write(lun_diag,"(a,3i5,3es14.6)") 'NR',kstep,izone,kit,testm(izb),testc(izb),testc2(izb)
          EndIf
        EndDo
      EndIf

      ! Check that all zones are converged
      If ( .not. any(iterate) ) Exit
    EndDo

    If ( idiag >= 2 ) Then
      !__dir_update &
      !__dir_wait &
      !__dir_host(inr,xtot,testn,toln)
      Do izb = zb_lo, zb_hi
        If ( inr(izb) >= 0 ) Then
          izone = izb + szbatch - zb_lo
          If ( inr(izb) > 0 ) Then
            Write(lun_diag,"(a,3i5,3es12.4)") 'Conv',kstep,izone,inr(izb),xtot(izb),testn(izb),toln(izb)
          ElseIf ( inr(izb) == 0 ) Then
            Write(lun_diag,"(a,3i5,2es10.2)") 'BE Failure',izone,inr(izb),kitmx,xtot(izb),testn(izb)
          EndIf
        EndIf
      EndDo
    EndIf

    !__dir_exit_data &
    !__dir_async &
    !__dir_delete(iterate,eval_rates,rebuild,testc,testc2,testm,testn,toln) &
    !__dir_delete(xtot,xtot_init,rdt,mult,yrhs,dy,reldy,t9rhs,dt9,relt9)

    stop_timer = xnet_wtime()
    timer_nraph = timer_nraph + stop_timer

    Return
  End Subroutine step_be

End Module xnet_integrate_be
