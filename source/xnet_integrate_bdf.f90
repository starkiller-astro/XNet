!***************************************************************************************************
! xnet_integrate_bdf.f90 10/18/17
! Backward differentiation formula ("Gear" methods; BDF) solver
!
! The routines in this file perform the BDF time integration for the thermonuclear reaction network.
! This particular implementation was adapated from:
!   https://github.com/starkiller-astro/Microphysics/tree/master/integration/VBDF
!     Original author: Matthew Emmett
!     Modifications by: Adam Jacobs and Mike Zingale
!
! References:
!    Longland, R., Martin, D. and José, J. (2014) Performance improvements for nuclear reaction
!      network integration, 563, A67. doi.org/10.1051/0004-6361/201321958
!    Brown, P., Byrne, G. & Hindmarsh, A. (1989) VODE: A Variable-Coefficient ODE Solver, SIAM J. 
!      Sci. Stat. Comput., 10, 1038. doi.org/10.1137/0910062
!    Jackson, K. R. & Sacks-Davis, R. (1980) An Alternative Implementation of Variable Step-Size
!      Multistep Formulas for Stiff ODEs, ACM Trans. Math. Softw., 6, 295.
!      doi.org/10.1145/355900.355903
!    Byrne, G. D. & Hindmarsh, A. C. (1975) A Polyalgorithm for the Numerical Solution of Ordinary
!      Differential Equations, ACM Trans. Math. Soft., 1, 71. doi.org/10.1145/355626.355636
!    Gear, C. W. (1971) The Automatic Integration of Ordinary Differential Equations, Commun. ACM,
!      14, 176. doi.org/10.1145/362566.362571
!    Nordsieck, A. (1962) On numerical integration of ordinary differential equations, Mathematics
!      of Computation, 16, 22. doi.org/10.1090/S0025-5718-1962-0136519-5
! 
!***************************************************************************************************

Module xnet_integrate_bdf
  Use nuclear_data, Only: ny
  Use xnet_controls, Only: idiag, iheat, iscrn, nzbatchmx, szbatch, lzactive, lun_diag, lun_stdout
  Use xnet_types, Only: dp
  Implicit None
  Private

  Logical, Parameter :: check_stability = .false.

  Integer, Parameter :: max_order = 5     ! Maximum BDF order
  Integer, Parameter :: max_it_ts = 20    ! Maximum number of attempts to advance step
  Integer, Parameter :: max_it_nr = 3     ! Maximum number of NR iterations with same Jacobian
  Integer, Parameter :: max_p_age = -1    ! Maximum number of steps between Newton matrix evaluations
  Integer, Parameter :: max_j_age = -1    ! Maximum number of steps between Jacobian evaluations
  Integer            :: max_it_nrslv      ! Maximum number of NR iterations per solve attempt (=max(kitmx,max_it_nr))

  Integer, Parameter :: max_nerr_ts = 10     ! Maximum number of error test failures before giving up
  Integer, Parameter :: max_ncf_ts = 10      ! Maximum number of NR convergence failures before giving up
  Integer, Parameter :: nerr_ts_small = 2    ! Number of solver iteration failures before forcing eta = min(eta,eta_reset)
  Integer, Parameter :: nerr_ts_qreduce = 3  ! Number of solver iteration failures before forcing an order reduction
  Integer, Parameter :: nstep_qwait = 10     ! Number of steps to delay change in order if q==1 and nerr_ts > nerr_ts_qreduce

  Real(dp), Parameter :: dt_min = 1.0e-24_dp      ! Minimum step-size
  Real(dp), Parameter :: eta_min = 0.1_dp         ! Minimum step-size shrink factor
  Real(dp), Parameter :: eta_small = 0.2_dp       ! Step-size reduction factor if nerr_ts > nerr_ts_small
  Real(dp), Parameter :: eta_reduce = 0.25_dp     ! Factor to reduce step-size if NR iteration fails to converge
  Real(dp), Parameter :: eta_max = 10.0_dp        ! Maximum step-size growth factor
  Real(dp), Parameter :: eta_thresh = 1.5_dp      ! Step-size growth threshold
  Real(dp), Parameter :: dgmaxp = 0.3_dp          ! Maximum change in gamma between Newton matrix evaluations
  Real(dp), Parameter :: dgmaxj = 0.2_dp          ! Maximum change in gamma between Jacobian evaluations
  Real(dp), Parameter :: tol_nr = 0.1_dp          ! Used to set convergence criteria for NR iteration
  Real(dp), Parameter :: cr_down = 0.3_dp         ! Used in estimation of the convergence rate
  Real(dp), Parameter :: cr_diverge = 2.0_dp      ! Stop NR iteration if convergence rate exceeds this
  Real(dp), Parameter :: small = 1.0e-10_dp       ! Used in stability detection

  Integer, Allocatable :: A_pascal(:,:)    ! Pascal triangle matrix
  Integer, Allocatable :: Ainv_pascal(:,:) ! Pascal triangle matrix inverse

  Logical, Allocatable :: bdf_active(:) ! Logical mask for zones being integrated
  Logical, Allocatable :: refactor(:)   ! Logical mask for zones being re-factored
  Logical, Allocatable :: retry_ts(:)   ! Logical mask for zones to retry non-linear solve
  !$omp threadprivate(bdf_active,refactor,retry_ts)

  Integer, Allocatable :: ierr_nr(:)    ! Convergence status flag for NR iterations
  Integer, Allocatable :: ierr_ts(:)    ! Convergence status flag for solve
  Integer, Allocatable :: nit_nrslv(:)  ! Number of NR iterations in current solve
  Integer, Allocatable :: nit_nr(:)     ! Number of NR iterations in current solve with same Jacobian
  Integer, Allocatable :: nit_ts(:)     ! Number of attempts to advance step
  Integer, Allocatable :: nerr_ts(:)    ! Number of solver error test failures during step
  Integer, Allocatable :: ncf_ts(:)     ! Number of solver NR convergence failures during step
  !$omp threadprivate(ierr_nr,ierr_ts,nit_nrslv,nit_nr,nit_ts,nerr_ts,ncf_ts)

  Integer, Allocatable :: q(:)          ! Current order
  Integer, Allocatable :: deltaq(:)     ! Change in order for next step
  Integer, Allocatable :: q_age(:)      ! Number of consecutive steps at current order
  Integer, Allocatable :: j_age(:)      ! Number of times Jacobian has been reused 
  Integer, Allocatable :: p_age(:)      ! Number of times Newton matrix has been reused
  Integer, Allocatable :: nscon(:)      ! Number of steps since last stability limit detection
  !$omp threadprivate(q,deltaq,q_age,j_age,p_age,nscon)

  Real(dp), Allocatable :: sddat(:,:,:) ! Scaled derivative data used to detect stability limits
  !$omp threadprivate(sddat)

  Real(dp), Allocatable :: z(:,:,:)     ! Nordsieck vector at current time t
  Real(dp), Allocatable :: z0(:,:,:)    ! Nordsieck vector at beginning of solve for current time t
  Real(dp), Allocatable :: zt0(:,:,:)   ! Predicted Nordsieck vector at trial time tt
  !$omp threadprivate(z,z0,zt0)

  Real(dp), Allocatable :: hvec(:,:)    ! Step sizes
  Real(dp), Allocatable :: lvec(:,:)    ! Coefficients for updating Nordsieck vector
  Real(dp), Allocatable :: dt_scale(:)  ! Last stepsize used to rescale Nordsieck vector
  Real(dp), Allocatable :: etaq(:,:)    ! Change in timestep dictated by change in order
  Real(dp), Allocatable :: eta(:)       ! Change in timestep between iterations
  Real(dp), Allocatable :: eta_next(:)  ! Change in timestep between for next step
  Real(dp), Allocatable :: gam(:)       ! Current gamma = hvec(0)/lvec(1)
  Real(dp), Allocatable :: gamhat(:)    ! Gamma last used to build Newton matrix
  Real(dp), Allocatable :: gamratio(:)  ! gam/gamhat
  !$omp threadprivate(hvec,lvec,dt_scale,etaq,eta,eta_next,gam,gamhat,gamratio)

  Real(dp), Allocatable :: ewt(:,:)     ! Error weights
  Real(dp), Allocatable :: acor(:,:)    ! Accumulated correction vector
  Real(dp), Allocatable :: acorp(:,:)   ! Accumulated correction vector, previous step
  Real(dp), Allocatable :: acnrm(:)     ! Weighted norm of accumulated correction vector
  Real(dp), Allocatable :: crate(:)     ! Estimated corrector convergence rate
  Real(dp), Allocatable :: tq(:,:)      ! Error estimation coefficients
  Real(dp), Allocatable :: tq2save(:)   ! Saved value of tq(2)
  Real(dp), Allocatable :: rtol(:)      ! Relative tolerances
  Real(dp), Allocatable :: atol(:)      ! Absolute tolerances
  !$omp threadprivate(ewt,acor,acorp,acnrm,crate,tq,tq2save)

  Real(dp), Allocatable :: ydot0(:,:)   ! Saved abundance derivatives from beginning of timestep
  Real(dp), Allocatable :: t9dot0(:)    ! Saved temperature derivatives from beginning of timestep
  !$omp threadprivate(ydot0,t9dot0)

  Integer :: neq

  Integer, Parameter :: BDF_STATUS_SKIP    = -1
  Integer, Parameter :: BDF_STATUS_SUCCESS =  0
  Integer, Parameter :: BDF_STATUS_FAIL    =  1

  Public :: solve_bdf, bdf_init

  Private :: bdf_reset, bdf_adjust, bdf_update, bdf_predict, bdf_step, bdf_check
  Private :: bdf_correct, bdf_eta, bdf_stab
  Private :: bdf_stab_detect, bdf_adjust_order, bdf_rescale, bdf_restore, pascal_build, normw

Contains

  Subroutine solve_bdf(kstep,ierr)
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    Use xnet_abundances, Only: y, yo, yt, ydot
    Use xnet_conditions, Only: t, to, tt, tdel, tdel_next, t9, t9o, t9t, t9dot, rho, rhoo, rhot, &
      yeo, ye, yet, nt, nto, ntt
    Use xnet_controls, Only: kitmx, kmon, ktot
    Use xnet_timers, Only: xnet_wtime, start_timer, stop_timer, timer_tstep
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep

    ! Input/Output variables
    Integer, Intent(inout) :: ierr(:) ! On input,  = BDF_STATUS_SUCCESS indicates active zone
                                      ! On output, = BDF_STATUS_FAIL    if zone fails to converge

    ! Local variables
    Integer :: kts, izb, izone

    ! Initial setup and allocations
    If ( kstep == 1 ) Call bdf_reset

    start_timer = xnet_wtime()
    timer_tstep = timer_tstep - start_timer

    ! Copy status flag to local copy
    ierr_ts(:) = ierr(:)

    ! If the zone has previously converged or failed, do not iterate
    Where ( ierr_ts(:) /= BDF_STATUS_SUCCESS )
      ierr_nr(:) = BDF_STATUS_SKIP
    ElseWhere
      ierr_nr(:) = BDF_STATUS_FAIL
    EndWhere

    ! Make adjustments to account for any change in stepsize or order from previous step
    If ( kstep /= 1 ) Call bdf_adjust(kstep)
    ydot0(:,:) = ydot(:,:)
    If ( iheat > 0 ) t9dot0(:) = t9dot(:)

    ! Reset maximum stepsize for next step
    eta_next(:) = eta_max

    ! Reset step counters
    nit_ts(:) = 0
    nerr_ts(:) = 0
    ncf_ts(:) = 0
    Do kts = 1, max_it_ts

      bdf_active(:) = ( ierr_nr(:) == BDF_STATUS_FAIL )

      Call bdf_update(kstep)
      Call bdf_predict(kstep)
      Call bdf_step(kstep)

      ! Check for errors and prepare for next attempt if necessary
      Call bdf_check(kstep)

      If ( .not. any( ierr_nr(:) == BDF_STATUS_FAIL ) ) Exit

    EndDo

    Call bdf_correct(kstep)
    Call bdf_eta(kstep)
    If ( check_stability ) Call bdf_stab(kstep)

    ! Set the factor for the next timestep
    tdel_next(:) = eta(:) * tdel(:)

    Do izb = 1, nzbatchmx
      If ( ierr_nr(izb) == BDF_STATUS_SUCCESS ) Then
        yt(:,izb) = z(1:ny,0,izb)
        If ( iheat > 0 ) t9t(izb) = z(neq,0,izb)
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
        yo(:,izb) = y(:,izb)
        y(:,izb) = yt(:,izb)
      EndIf
    EndDo

    ! Record iteration counters
    kmon(1,:) = nit_ts(:)
    ktot(1,:) = ktot(1,:) + nit_ts(:)
    kmon(2,:) = nit_nrslv(:)

    ! Log TS success/failure
    Do izb = 1, nzbatchmx
      izone = izb + szbatch - 1
      If ( ierr_nr(izb) == BDF_STATUS_SUCCESS ) Then
        ierr_ts(izb) = BDF_STATUS_SUCCESS
        If ( idiag >= 2 ) Write(lun_diag,"(a,2i5,2i3)") &
          'BDF TS Success',kstep,izone,nit_ts(izb),nit_nr(izb)
      ElseIf ( ierr_nr(izb) == BDF_STATUS_FAIL .or. ierr_ts(izb) == BDF_STATUS_FAIL ) Then
        ierr_ts(izb) = BDF_STATUS_FAIL
        If ( idiag >= 0 ) Write(lun_diag,"(a,2i5,4es12.4,2i3)") &
          'BDF TS Fail',kstep,izone,t(izb),tdel(izb),t9t(izb),rhot(izb),nit_nr(izb),nit_ts(izb)
      EndIf
    EndDo

    ! Copy local status flag to output
    ierr(:) = ierr_ts(:)

    stop_timer = xnet_wtime()
    timer_tstep = timer_tstep + stop_timer

    Return
  End Subroutine solve_bdf

  Subroutine bdf_init
    !-----------------------------------------------------------------------------------------------
    ! This routine initializes the BDF data structures and prepares for the first integration step.
    !-----------------------------------------------------------------------------------------------
    Use xnet_controls, Only: kitmx, yacc, ymin, tolc, tolt9
    Implicit None

    ! Local variables
    Integer, Parameter :: nheat = 1 ! Currently only one self-heating equation

    ! Set maximum number of NR iterations per solve attempt
    max_it_nrslv = max( kitmx, max_it_nr )

    ! Initialize Pascal matrix
    If ( .not. allocated(A_pascal) ) Then
      Allocate (A_pascal(0:max_order,0:max_order))
      Allocate (Ainv_pascal(0:max_order,0:max_order))
      Call pascal_build
    EndIf

    ! Set number of RHS equations (ny species + possible self-heating)
    If ( iheat > 0 ) Then
      neq = ny + nheat
    Else
      neq = ny
    EndIf

    ! Set tolerances
    If ( .not. allocated(rtol) ) Allocate (rtol(neq))
    If ( .not. allocated(atol) ) Allocate (atol(neq))
    rtol(1:ny) = tolc
    atol(1:ny) = yacc
    If ( iheat > 0 ) Then
      rtol(neq) = tolt9
      atol(neq) = 1.0e-99_dp
    EndIf

    !$omp parallel default(shared)

    If ( .not. allocated(bdf_active) ) Allocate (bdf_active(nzbatchmx))
    If ( .not. allocated(refactor) ) Allocate (refactor(nzbatchmx))
    If ( .not. allocated(retry_ts) ) Allocate (retry_ts(nzbatchmx))
    bdf_active(:) = .false.
    refactor(:) = .false.
    retry_ts(:) = .false.

    If ( .not. allocated(nit_nrslv) ) Allocate (nit_nrslv(nzbatchmx))
    If ( .not. allocated(nit_nr) ) Allocate (nit_nr(nzbatchmx))
    If ( .not. allocated(nit_ts) ) Allocate (nit_ts(nzbatchmx))
    If ( .not. allocated(ierr_nr) ) Allocate (ierr_nr(nzbatchmx))
    If ( .not. allocated(ierr_ts) ) Allocate (ierr_ts(nzbatchmx))
    If ( .not. allocated(nerr_ts) ) Allocate (nerr_ts(nzbatchmx))
    If ( .not. allocated(ncf_ts) ) Allocate (ncf_ts(nzbatchmx))
    nit_nrslv(:) = 0
    nit_nr(:) = 0
    nit_ts(:) = 0
    ierr_nr(:) = 0
    ierr_ts(:) = 0
    nerr_ts(:) = 0
    ncf_ts(:) = 0

    If ( .not. allocated(q) ) Allocate (q(nzbatchmx))
    If ( .not. allocated(deltaq) ) Allocate (deltaq(nzbatchmx))
    If ( .not. allocated(j_age) ) Allocate (j_age(nzbatchmx))
    If ( .not. allocated(p_age) ) Allocate (p_age(nzbatchmx))
    If ( .not. allocated(q_age) ) Allocate (q_age(nzbatchmx))
    If ( .not. allocated(nscon) ) Allocate (nscon(nzbatchmx))
    q(:) = 1
    deltaq(:) = 0
    j_age(:) = max_j_age + 1 ! This forces a rebuild/refactor on first step
    p_age(:) = max_p_age + 1 ! This forces a rebuild/refactor on first step
    q_age(:) = 0
    nscon(:) = 0

    If ( .not. allocated(sddat) ) Allocate (sddat(5,3,nzbatchmx))
    sddat(:,:,:) = 0.0_dp

    If ( .not. allocated(z)) Allocate (z(neq,0:max_order,nzbatchmx))
    If ( .not. allocated(z0)) Allocate (z0(neq,0:max_order,nzbatchmx))
    If ( .not. allocated(zt0) ) Allocate (zt0(neq,0:max_order,nzbatchmx))
    z(:,:,:) = 0.0_dp
    z0(:,:,:) = 0.0_dp
    zt0(:,:,:) = 0.0_dp

    If ( .not. allocated(hvec) ) Allocate (hvec(0:max_order,nzbatchmx))
    If ( .not. allocated(lvec) ) Allocate (lvec(0:max_order,nzbatchmx))
    If ( .not. allocated(dt_scale) ) Allocate (dt_scale(nzbatchmx))
    If ( .not. allocated(etaq) ) Allocate (etaq(-1:1,nzbatchmx))
    If ( .not. allocated(eta) ) Allocate (eta(nzbatchmx))
    If ( .not. allocated(eta_next) ) Allocate (eta_next(nzbatchmx))
    If ( .not. allocated(gam) ) Allocate (gam(nzbatchmx))
    If ( .not. allocated(gamhat) ) Allocate (gamhat(nzbatchmx))
    If ( .not. allocated(gamratio) ) Allocate (gamratio(nzbatchmx))
    hvec(:,:) = 0.0_dp
    lvec(:,:) = 0.0_dp
    dt_scale(:) = 0.0_dp
    etaq(:,:) = 0.0_dp
    eta(:) = 0.0_dp
    eta_next(:) = 0.0_dp
    gam(:) = 0.0_dp
    gamhat(:) = 0.0_dp
    gamratio(:) = 0.0_dp

    If ( .not. allocated(ewt) ) Allocate (ewt(neq,nzbatchmx))
    If ( .not. allocated(acor) ) Allocate (acor(neq,nzbatchmx))
    If ( .not. allocated(acorp) ) Allocate (acorp(neq,nzbatchmx))
    If ( .not. allocated(acnrm) ) Allocate (acnrm(nzbatchmx))
    If ( .not. allocated(crate) ) Allocate (crate(nzbatchmx))
    If ( .not. allocated(tq) ) Allocate (tq(-1:2,nzbatchmx))
    If ( .not. allocated(tq2save) ) Allocate (tq2save(nzbatchmx))
    ewt(:,:) = 0.0_dp
    acor(:,:) = 0.0_dp
    acorp(:,:) = 0.0_dp
    acnrm(:) = 0.0_dp
    crate(:) = 0.0_dp
    tq(:,:) = 0.0_dp
    tq2save(:) = 0.0_dp

    If ( .not. allocated(ydot0) ) Allocate (ydot0(ny,nzbatchmx))
    If ( .not. allocated(t9dot0) ) Allocate (t9dot0(nzbatchmx))
    ydot0(:,:) = 0.0_dp
    t9dot0(:) = 0.0_dp

    !$omp end parallel

    Return
  End Subroutine bdf_init

  Subroutine bdf_reset
    !-----------------------------------------------------------------------------------------------
    ! This routine initializes the BDF data structures and prepares for the first integration step.
    !-----------------------------------------------------------------------------------------------
    Use xnet_abundances, Only: yt, ydot
    Use xnet_conditions, Only: t9t, t9dot, tdel
    Implicit None

    ! Local variables
    Integer :: izb

    bdf_active(:) = .false.
    refactor(:) = .false.
    retry_ts(:) = .false.

    nit_nrslv(:) = 0
    nit_nr(:) = 0
    nit_ts(:) = 0
    ierr_nr(:) = 0
    ierr_ts(:) = 0
    nerr_ts(:) = 0
    ncf_ts(:) = 0

    q(:) = 1
    deltaq(:) = 0
    j_age(:) = max_j_age + 1 ! This forces a rebuild/refactor on first step
    p_age(:) = max_p_age + 1 ! This forces a rebuild/refactor on first step
    q_age(:) = 0
    nscon(:) = 0

    sddat(:,:,:) = 0.0_dp

    z(:,:,:) = 0.0_dp
    Do izb = 1, nzbatchmx
      z(1:ny,0,izb) = yt(:,izb)
      z(1:ny,1,izb) = ydot(:,izb) * tdel(izb)
      If ( iheat > 0 ) Then
        z(neq,0,izb) = t9t(izb)
        z(neq,1,izb) = t9dot(izb) * tdel(izb)
      EndIf
    EndDo
    z0(:,:,:) = 0.0_dp
    zt0(:,:,:) = 0.0_dp

    Do izb = 1, nzbatchmx
      hvec(:,izb) = tdel(izb)
    EndDo
    lvec(:,:) = 0.0_dp
    dt_scale(:) = tdel(:)
    etaq(:,:) = 0.0_dp
    eta(:) = 0.0_dp
    eta_next(:) = eta_max
    gam(:) = tdel(:)
    gamhat(:) = tdel(:)
    gamratio(:) = 1.0_dp

    ewt(:,:) = 0.0_dp
    acor(:,:) = 0.0_dp
    acorp(:,:) = 0.0_dp
    acnrm(:) = 0.0_dp
    crate(:) = 0.0_dp
    tq(:,:) = 0.0_dp
    tq2save(:) = 0.0_dp

    ydot0(:,:) = ydot(:,:)
    t9dot0(:) = t9dot(:)

    Return
  End Subroutine bdf_reset

  Subroutine bdf_adjust(kstep)
    !-----------------------------------------------------------------------------------------------
    ! This routine handles any required adjustments to the Nordsieck vector if the current stepsize
    ! is different than the previous. Also, any change in order is applied here.
    !-----------------------------------------------------------------------------------------------
    Use xnet_conditions, Only: tdel, tdel_old
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep

    ! Local variables
    Integer :: i, j, izb
    Logical :: rescale(nzbatchmx)

    ! See if stepsize has changed from last step
    rescale(:) = .false.
    Where ( ierr_ts(:) == BDF_STATUS_SUCCESS )
      Where ( abs(tdel(:)/tdel_old(:) - 1.0_dp) > epsilon(0.0_dp) )
        eta(:) = tdel(:) / tdel_old(:)
        rescale(:) = .true.
      EndWhere
    EndWhere

    ! Apply any change in order determined during previous step
    Call bdf_adjust_order(kstep)

    ! Scale Nordsieck vector to account for change in stepsize
    Call bdf_rescale(rescale)

    Return
  End Subroutine bdf_adjust

  Subroutine bdf_predict(kstep)
    !-----------------------------------------------------------------------------------------------
    ! This subroutine computes the predictor step (i.e. apply Pascal matrix):
    !   $z^{(0)}_{n+1} = A(q) z_n$
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: nname, ny
    Use xnet_abundances, Only: yt, ydot
    Use xnet_conditions, Only: tdel, t9t, t9dot
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep

    ! Local variables
    Integer :: i, j, izb, izone

    Do izb = 1, nzbatchmx
      If ( bdf_active(izb) ) Then
        z0(:,:,izb) = z(:,:,izb)
        zt0(:,:,izb) = z(:,:,izb)
        Do i = 0, q(izb)
          Do j = i+1, q(izb)
            zt0(:,i,izb) = zt0(:,i,izb) + z(:,j,izb) * A_pascal(j,i)
          EndDo
        EndDo
      EndIf
    EndDo

    If ( idiag >= 3 ) Then
      Do izb = 1, nzbatchmx
        izone = izb + szbatch - 1
        Write(lun_diag,"(a,2i5,i3,1es12.4)") 'BDF Predict',kstep,izone,q(izb),tdel(izb)
        Write(lun_diag,"(2x,a5,6es15.7)") (nname(i), &
          zt0(i,0,izb),z(i,0,izb),yt(i,izb), &
          zt0(i,1,izb),z(i,1,izb),ydot(i,izb)*tdel(izb),i=1,ny)
        If ( iheat > 0 ) Write(lun_diag,"(2x,a5,6es23.15)") 'T9', &
          zt0(neq,0,izb),z(neq,0,izb),t9t(izb), &
          zt0(neq,1,izb),z(neq,1,izb),t9dot(izb)*tdel(izb)
      EndDo
    EndIf

    Return
  End Subroutine bdf_predict

  Subroutine bdf_update(kstep)
    !-----------------------------------------------------------------------------------------------
    ! This subroutine computes the l coeffiecients for updating the Nordsieck vector (see Section 5
    ! of Jackson & Sacks-Davis (1980)) and error coefficients tq.
    !
    ! The components of the array l are the coefficients of a polynomial:
    !   Lambda(x) = l_0 + l_1 x + ... + l_q x^q, given by
    !
    !                                  q-1
    !   Lambda(x) = (1 + x / xi*_q) * PRODUCT (1 + x / xi_i) , where
    !                                  i=1
    !   xi_i = [t_n - t_(n-i)] / h.
    !
    ! The array tq is set to test quantities used in the convergence
    ! test, the error test, and the selection of h at a new order.
    !   tq(-1) : error coeff. for order q-1
    !   tq(0)  : error coeff. for order q
    !   tq(1)  : error coeff. for order q+1
    !   tq(2)  : error coeff. for order q+1 used to get the order q+2 derivative vector
    !
    ! This routine was largely adapted from CVODE
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: nname
    Use xnet_conditions, Only: tdel, tt
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep

    ! Local variables
    Real(dp) :: hsum, xi_inv, xistar_inv, alpha0, alpha0_hat, c, cinv
    Real(dp) :: a1, a2, a3, a4, a5, a6
    Character(5) :: ewtname(2)
    Integer :: iewt(2), i, j, izb, izone

    Do izb = 1, nzbatchmx
      If ( bdf_active(izb) ) Then

        ! Compute lvec and tq
        lvec(0,izb) = 1.0_dp
        lvec(1,izb) = 1.0_dp
        lvec(2:max_order,izb) = 0.0_dp
        xi_inv = 1.0_dp
        xistar_inv = 1.0_dp
        alpha0 = -1.0_dp
        alpha0_hat = -1.0_dp
        hsum = tdel(izb)
        If ( q(izb) > 1 ) Then
          Do j = 2, q(izb)-1
            hsum = hsum + hvec(j-1,izb)
            xi_inv = tdel(izb) / hsum
            alpha0 = alpha0 - 1.0_dp / real(j,dp)

            ! lvec(i) are coefficients of product(1 to j) (1 + x/xi_i)
            Do i = j, 1, -1
              lvec(i,izb) = lvec(i,izb) + lvec(i-1,izb) * xi_inv
            EndDo
          EndDo

          ! j = q
          alpha0 = alpha0 - 1.0_dp / real(q(izb),dp)
          xistar_inv = -lvec(1,izb) - alpha0
          hsum = hsum + hvec(q(izb)-1,izb)
          xi_inv = tdel(izb) / hsum
          alpha0_hat = -lvec(1,izb) - xi_inv
          Do i = q(izb), 1, -1
            lvec(i,izb) = lvec(i,izb) + lvec(i-1,izb) * xistar_inv
          EndDo
        EndIf

        ! Compute tq
        a1 = 1.0_dp - alpha0_hat + alpha0
        a2 = 1.0_dp + q(izb) * a1
        tq(0,izb) = abs(a1 / (alpha0*a2))
        tq(2,izb) = abs(a2*xistar_inv / (lvec(q(izb),izb)*xi_inv))
        If ( q(izb) > 1 ) Then
          c = xistar_inv / lvec(q(izb),izb)
          a3 = alpha0 + 1.0_dp / q(izb)
          a4 = alpha0_hat + xi_inv
          tq(-1,izb) = abs(c*(1.0_dp - a4 + a3)/a3)
        Else
          tq(-1,izb) = 1.0_dp
        EndIf

        hsum = hsum + hvec(q(izb),izb)
        xi_inv = tdel(izb) / hsum
        a5 = alpha0 - 1.0_dp / real(q(izb)+1,dp)
        a6 = alpha0_hat - xi_inv
        cinv = (1.0_dp - a6 + a5) / a2
        tq(1,izb) = abs(cinv / (xi_inv * real(q(izb)+2,dp) * a5))

        ! Pre-compute error weights
        Call error_weights( z(:,0,izb), ewt(:,izb) )

        gam(izb) = tdel(izb) / lvec(1,izb)
        If ( kstep > 1 ) Then
          gamratio(izb) = gam(izb) / gamhat(izb)
        Else
          gamhat(izb) = gam(izb)
          gamratio(izb) = 1.0_dp
        EndIf

        If ( idiag >= 3 ) Then
          izone = izb + szbatch - 1
          Write(lun_diag,"(a,2i5,i3,7es12.4)") 'BDF Update',kstep,izone,q(izb),tdel(izb),tt(izb)
          Write(lun_diag,"(2x,a5,6es12.4)") ' H',(hvec(i,izb),i=0,max_order)
          Write(lun_diag,"(2x,a5,4es12.4)") 'TQ',(tq(i,izb),i=-1,2)
          iewt(1) = minloc(ewt(:,izb),dim=1)
          iewt(2) = maxloc(ewt(:,izb),dim=1)
          Do i = 1, 2
            If ( iewt(i) <= ny ) Then
              ewtname(i) = nname(iewt(i))
            Else
              ewtname(i) = 'T9'
            EndIf
          EndDo
          Write(lun_diag,"(2x,a5,2(a5,es14.7))") 'EWT',(ewtname(i),ewt(iewt(i),izb),i=1,2)
        EndIf
      EndIf
    EndDo

    Return
  End Subroutine bdf_update

  Subroutine bdf_step(kstep)
    !-----------------------------------------------------------------------------------------------
    ! This routine iteratively calculates the correction vector necessary to adjust each corrector
    ! step m until convergence.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: aa, nname
    Use xnet_abundances, Only: y, yt, ydot
    Use xnet_conditions, Only: tt, tdel, ntt, t9t, t9dot, rhot, ntt, cv, nh, t9rhofind
    Use xnet_controls, Only: ymin
    Use xnet_integrate, Only: cross_sect, yderiv
    Use xnet_jacobian, Only: jacobian_build, jacobian_scale, jacobian_decomp, jacobian_bksub
    Use xnet_timers, Only: xnet_wtime, start_timer, stop_timer, timer_nraph
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep

    ! Local variables
    Real(dp) :: diag(nzbatchmx), cstar
    Real(dp) :: yrhs(ny,nzbatchmx), dy(ny,nzbatchmx)
    Real(dp) :: t9rhs(nzbatchmx), dt9(nzbatchmx)
    Real(dp) :: dvec(neq), del(nzbatchmx), delp(nzbatchmx), dcon(nzbatchmx)
    Integer :: idymx, i, k, kit, izb, izone
    Logical :: rebuild(nzbatchmx), iterate(nzbatchmx), converged(nzbatchmx), eval_rates(nzbatchmx)

    start_timer = xnet_wtime()
    timer_nraph = timer_nraph - start_timer

    ! Load prediction into abundance vector
    Do izb = 1, nzbatchmx
      If ( bdf_active(izb) ) Then
        acor(:,izb) = 0.0_dp
        yt(:,izb) = zt0(1:ny,0,izb)
        If ( iheat > 0 ) t9t(izb) = zt0(neq,0,izb)
      EndIf
    EndDo
    If ( iheat == 0 ) Call t9rhofind(kstep,tt,ntt,t9t,rhot,mask_in = retry_ts(:))

    ! Determine which zones to rebuild the Newton matrix and factor
    refactor(:) = ( bdf_active(:) .and. &
      ( retry_ts(:) .or. p_age(:) > max_p_age .or. abs(gamratio(:)-1.0_dp) > dgmaxp ) )

    ! Determine which zones' Jacobians need to be rebuilt
    rebuild(:) = ( refactor(:) .and. j_age(:) > max_j_age )

    ! Do the NR iteration
    Where ( bdf_active(:) )
      nit_nr(:) = 0
      nit_nrslv(:) = 0
    EndWhere
    iterate(:) = bdf_active(:)
    converged(:) = .false.
    crate(:) = 1.0_dp
    delp(:) = 0.0_dp
    diag(:) = 1.0_dp
    Do kit = 1, max_it_nrslv

      ! Increment NR counter
      Where ( iterate(:) )
        nit_nr(:) = nit_nr(:) + 1
        nit_nrslv(:) = nit_nrslv(:) + 1
      EndWhere

      ! Only re-evaluate rates if first iteration or if need to re-evaluate Jacobian
      If ( iheat > 0 ) Then
        eval_rates(:) = iterate(:)
      ElseIf ( kit == 1 ) Then
        eval_rates(:) = ( iterate(:) .and. nh(:) > 1 )
      Else
        eval_rates(:) = .false.
      EndIf

      ! Update abundance derivatives and build Jacobian
      Call cross_sect(mask_in = eval_rates(:))
      Call yderiv(mask_in = iterate(:))
      Call jacobian_build(mask_in = rebuild(:))
      Call jacobian_scale(diag,-gam,mask_in = refactor(:))
      Call jacobian_decomp(kstep,mask_in = refactor(:))
      If ( idiag >= 4 ) Then
        Do izb = 1, nzbatchmx
          izone = izb + szbatch - 1
          If ( refactor(izb) ) Write(lun_diag,"(a,3i5,2es14.7,i5)") &
            'BDF Refactor',kstep,izone,nit_nr(izb),gam(izb),abs(gamratio(izb)-1.0_dp),p_age(izb)
          If ( rebuild(izb) ) Write(lun_diag,"(a,3i5,2es14.7,i5)") &
            'BDF Rebuild',kstep,izone,nit_nr(izb),gam(izb),abs(gamratio(izb)-1.0_dp),j_age(izb)
        EndDo
      EndIf
      Where ( refactor(:) )
        gamhat(:) = gam(:)
        gamratio(:) = 1.0_dp
        crate(:) = 1.0_dp
        p_age(:) = 0
        refactor(:) = .false.
      EndWhere
      Where ( rebuild(:) )
        j_age(:) = 0
        rebuild(:) = .false.
      EndWhere

      ! Calculate RHS
      Do izb = 1, nzbatchmx
        If ( iterate(izb) ) Then
          cstar = 2.0_dp / ( 1.0_dp + gamratio(izb) )
          yrhs(:,izb) =  cstar * ( gam(izb)*ydot(:,izb) - zt0(1:ny,1,izb)/lvec(1,izb) - acor(1:ny,izb) )
          If ( iheat > 0 ) t9rhs(izb) = cstar * ( gam(izb)*t9dot(izb) - zt0(neq,1,izb)/lvec(1,izb) - acor(neq,izb) )
        EndIf
      EndDo
      If ( idiag >= 4 ) Then
        Do izb = 1, nzbatchmx
          If ( iterate(izb) ) Then
            izone = izb + szbatch - 1
            Write(lun_diag,"(a,3i5,3es14.7)") &
              'BDF RHS',kstep,izone,nit_nr(izb),gam(izb),gamratio(izb),tdel(izb)
            Write(lun_diag,"(2x,a5,4es23.15)") &
              (nname(i),yrhs(i,izb),ydot(i,izb),yt(i,izb),zt0(i,0,izb),i=1,ny)
            If ( iheat > 0 ) Write(lun_diag,"(2x,a5,4es23.15)") &
              'T9',t9rhs(izb),t9dot(izb),t9t(izb),zt0(neq,0,izb)
          EndIf
        EndDo
      EndIf

      ! Solve using factorized Newton matrix
      Call jacobian_bksub(kstep,yrhs,dy,t9rhs,dt9,mask_in = iterate(:))

      ! Apply the iterative correction
      Do izb = 1, nzbatchmx
        If ( iterate(izb) ) Then

          ! Prevent abundance from becoming too small
          yt(:,izb) = zt0(1:ny,0,izb) + acor(1:ny,izb) + dy(:,izb)
          Where ( yt(1:ny,izb) < ymin )
            yt(1:ny,izb) = 0.0_dp
            dy(1:ny,izb) = -zt0(1:ny,0,izb) - acor(1:ny,izb)
          ElseWhere ( yt(1:ny,izb) * aa(1:ny) > 1.0_dp )
            yt(1:ny,izb) = 1.0_dp / aa(1:ny)
            dy(1:ny,izb) = 1.0_dp / aa(1:ny) - zt0(1:ny,0,izb) - acor(1:ny,izb)
          EndWhere
          acor(1:ny,izb) = acor(1:ny,izb) + dy(:,izb)
          dvec(1:ny) = dy(:,izb)
          If ( iheat > 0 ) Then
            dvec(neq) = dt9(izb)
            acor(neq,izb) = acor(neq,izb) + dt9(izb)
            t9t(izb) = zt0(neq,0,izb) + acor(neq,izb)
          EndIf
          del(izb) = normw( dvec, ewt(:,izb) )
          If ( nit_nr(izb) > 1 ) crate(izb) = max( cr_down*crate(izb), del(izb)/delp(izb) )
          dcon(izb) = del(izb) * min(1.0_dp,crate(izb)) * tq(0,izb)
        EndIf
      EndDo
      If ( idiag >= 3 ) Then
        Do izb = 1, nzbatchmx
          If ( iterate(izb) ) Then
            izone = izb + szbatch - 1
            If ( idiag >= 4 ) Then
              idymx = maxloc(dy(:,izb),dim=1)
              Write(lun_diag,"(a10,3i5,a5,2es23.15)") &
                'dY',kstep,izone,nit_nr(izb),nname(idymx),dy(idymx,izb),y(idymx,izb)
              If ( iheat > 0 ) Write(lun_diag,"(a10,3i5,5x,2es23.15)") &
                'dT9',kstep,izone,nit_nr(izb),dt9(izb),t9t(izb)
              If ( idiag >= 5 ) Write(lun_diag,"(a5,4es12.4)") &
                (nname(k),yt(k,izb),dy(k,izb),(aa(k)*dy(k,izb)),(aa(k)*yt(k,izb)),k=1,ny)
            EndIf
            Write(lun_diag,"(a,3i5,3es14.6)") &
              'BDF NR',kstep,izone,nit_nr(izb),del(izb),dcon(izb),crate(izb)
          EndIf
        EndDo
      EndIf

      ! Test for convergence
      Do izb = 1, nzbatchmx
        If ( iterate(izb) ) Then
          converged(izb) = ( dcon(izb) < tol_nr )

          ! If converged, remove from NR iteration
          If ( converged(izb) ) Then
            iterate(izb) = .false.
            If ( nit_nr(izb) > 1 ) Then
              acnrm(izb) = normw( acor(:,izb), ewt(:,izb) )
            Else
              acnrm(izb) = del(izb)
            EndIf

          ElseIf ( nit_nrslv(izb) == max_it_nrslv ) Then
            iterate(izb) = .false.
            nit_nr(izb) = max_it_nr + 1

          ! If at max number of NR iterations or if diverging....
          ElseIf ( nit_nr(izb) == max_it_nr .or. &
              ( nit_nr(izb) > 1 .and. del(izb) > delp(izb) * cr_diverge ) ) Then

            ! Reset the iteration and try again with a new Jacobian
            If ( j_age(izb) > 0 ) Then
              yt(:,izb) = zt0(1:ny,0,izb)
              If ( iheat > 0 ) t9t(izb) = zt0(neq,0,izb)
              refactor(izb) = .true.
              rebuild(izb) = ( abs(gamratio(izb) - 1.0_dp) < dgmaxj )
              nit_nr(izb) = 0
              acor(:,izb) = 0.0_dp
            Else
              iterate(izb) = .false.
              nit_nr(izb) = max_it_nr + 1
            EndIf
          EndIf

          ! Save for future iterations
          delp(izb) = del(izb)
        EndIf
      EndDo

      ! Check that all zones are converged
      If ( .not. any(iterate(:)) ) Exit

    EndDo
    
    If ( idiag >= 2 ) Then
      Do izb = 1, nzbatchmx
        izone = izb + szbatch - 1
        If ( converged(izb) ) Then
          Write(lun_diag,"(a,2i5,2i3,3es12.4,2i3)") &
            'BDF NR Conv',kstep,izone,nit_nr(izb),nit_nrslv(izb),dcon(izb),gamratio(izb),crate(izb),j_age(izb),p_age(izb)
        ElseIf ( bdf_active(izb) ) Then
          Write(lun_diag,"(a,2i5,2i3,3es12.4,2i3)") &
            'BDF NR Failure',kstep,izone,nit_nr(izb),nit_nrslv(izb),dcon(izb),gamratio(izb),crate(izb),j_age(izb),p_age(izb)
        EndIf
      EndDo
    EndIf

    stop_timer = xnet_wtime()
    timer_nraph = timer_nraph + stop_timer

    Return
  End Subroutine bdf_step

  Subroutine bdf_check(kstep)
    !-----------------------------------------------------------------------------------------------
    ! This routine checks for errors in the NR iteration and rescales the timestep as necessary.
    !-----------------------------------------------------------------------------------------------
    Use xnet_abundances, Only: y
    Use xnet_conditions, Only: t9, t, tdel
    Use xnet_controls, Only: kitmx, kmon, ktot
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep

    ! Local variables
    Real(dp), Parameter :: biasq = 6.0_dp
    Real(dp), Parameter :: addon = 1.0e-6_dp
    Real(dp) :: error
    Integer :: i, j, izb, izone
    Logical :: restore(nzbatchmx), rescale(nzbatchmx)

    retry_ts(:) = .false.
    restore(:) = .false.
    rescale(:) = .false.
    Do izb = 1, nzbatchmx
      If ( bdf_active(izb) ) Then
        izone = izb + szbatch - 1

        ! NR iteration failed to converge...
        If ( nit_nr(izb) > max_it_nr ) Then

          ! Increment convergence failure counter
          ncf_ts(izb) = ncf_ts(izb) + 1

          ! Prevent stepsize from growing after this step
          eta_next(izb) = 1.0_dp

          ! Solver failed too many times, give up
          If ( nit_ts(izb) >= max_it_ts ) Then
            ierr_nr(izb) = BDF_STATUS_SKIP
            ierr_ts(izb) = BDF_STATUS_FAIL
            Write(lun_stdout,"(a,2i5,1x,a)") 'BDF TS Fail',kstep,izone,'(Too many TS failures)'

          ! Solver failed to converge with minimum timestep, give up
          ElseIf ( tdel(izb) <= spacing(t(izb)) ) Then
            ierr_nr(izb) = BDF_STATUS_SKIP
            ierr_ts(izb) = BDF_STATUS_FAIL
            Write(lun_stdout,"(a,2i5,1x,a)") 'BDF TS Fail',kstep,izone,'(Timestep too small)'

          ! Solver failed to converge too many times, give up
          ElseIf ( ncf_ts(izb) == max_ncf_ts ) Then
            ierr_nr(izb) = BDF_STATUS_SKIP
            ierr_ts(izb) = BDF_STATUS_FAIL
            Write(lun_stdout,"(a,2i5,1x,a)") 'BDF TS Fail',kstep,izone,'(Too many NR failures)'

          ! Try again with smaller timestep
          Else
            ierr_nr(izb) = BDF_STATUS_FAIL
            ierr_ts(izb) = BDF_STATUS_SUCCESS
            retry_ts(izb) = .true.
            restore(izb) = .true.

            ! Reduce timestep
            eta(izb) = max( eta_reduce, spacing(t(izb))/tdel(izb) )
            rescale(izb) = .true.

            ! Force sparse solvers to redo the analysis phase of the decomposition
            kmon(2,izb) = kitmx + 1
          EndIf

        ! NR iteration converged, but error may still be too large
        Else ! nit_nr(izb) <= max_it_nr

          ! Check if we need to shrink dt
          error = tq(0,izb) * acnrm(izb)

          ! Error sufficiently small, successful step complete
          If ( error <= 1.0_dp ) Then

            ! Flag for removal from step loop
            ierr_nr(izb) = BDF_STATUS_SUCCESS
            ierr_ts(izb) = BDF_STATUS_SUCCESS

          ! Error too large
          Else ! error > 1.0_dp

            ! Increment error test failure counter
            nerr_ts(izb) = nerr_ts(izb) + 1

            ! Prevent stepsize from growing after this step
            eta_next(izb) = 1.0_dp

            If ( idiag >= 2 ) Write(lun_diag,"(a,2i5,2i3,es12.4)") &
              'BDF TS Error',kstep,izone,nit_ts(izb),nerr_ts(izb),error

            ! Solver failed error test too many times, give up
            If ( nerr_ts(izb) == max_nerr_ts ) Then
              ierr_nr(izb) = BDF_STATUS_SKIP
              ierr_ts(izb) = BDF_STATUS_FAIL
              Write(lun_stdout,"(a,2i5,1x,a)") 'BDF TS Fail',kstep,izone,'(Too many failed error tests)'

            ! Solver failed error test with minimum timestep, give up
            ElseIf ( tdel(izb) <= spacing(t(izb)) ) Then
              ierr_nr(izb) = BDF_STATUS_SKIP
              ierr_ts(izb) = BDF_STATUS_FAIL
              Write(lun_stdout,"(a,2i5,1x,a)") 'BDF TS Fail',kstep,izone,'(Timestep too small)'

            ! Solver failed error test after multiple attempts, reduce timestep, reduce order, and retry
            ElseIf ( nerr_ts(izb) > nerr_ts_qreduce ) Then
              ierr_nr(izb) = BDF_STATUS_FAIL
              ierr_ts(izb) = BDF_STATUS_SUCCESS
              retry_ts(izb) = .true.

              ! Reduce timestep
              eta(izb) = max( eta_min, spacing(t(izb))/tdel(izb) )
              rescale(izb) = .true.

              ! If already at order 1, restart...
              If ( q(izb) == 1 ) Then

                ! Wait a while before considering a change in order
                q_age(izb) = q(izb) - nstep_qwait

                ! Reload Nordsieck vector from scratch instead of restoring to state prior to prediction
                z(1:ny,0,izb) = y(:,izb)
                z(1:ny,1,izb) = ydot0(:,izb) * tdel(izb)
                If ( iheat > 0 ) Then
                  z(neq,0,izb) = t9(izb)
                  z(neq,1,izb) = t9dot0(izb) * tdel(izb)
                EndIf
                restore(izb) = .false.

              ! Otherwise, decrease order
              Else ! q(izb) > 1
                restore(izb) = .true.
                deltaq(izb) = -1
              EndIf

            ! Solver failed error test reduce timestep and retry
            Else ! nerr_ts(izb) <= nerr_ts_qreduce
              ierr_nr(izb) = BDF_STATUS_FAIL
              ierr_ts(izb) = BDF_STATUS_SUCCESS
              retry_ts(izb) = .true.
              restore(izb) = .true.

              ! Reduce timestep
              eta(izb) = 1.0_dp / ( (biasq * error) ** (1.0_dp/(q(izb)+1)) + addon )
              eta(izb) = min( eta_max, max( eta_min, spacing(t(izb))/tdel(izb), eta(izb) ) )
              If ( nerr_ts(izb) >= nerr_ts_small ) eta(izb) = min( eta_small, eta(izb) )
              rescale(izb) = .true.
            EndIf
          EndIf
        EndIf
      EndIf
    EndDo

    ! Record iteration counters
    Where ( bdf_active(:) )
      ktot(2,:) = ktot(2,:) + nit_nrslv(:)
      nit_ts(:) = nit_ts(:) + 1
    EndWhere

    ! Log the failed integration attempts
    If ( idiag >= 2 ) Then
      Do izb = 1, nzbatchmx
        If ( retry_ts(izb) ) Then
          izone = izb + szbatch - 1
          Write(lun_diag,"(a,i5,3i3,3es12.4,i3)") &
            'BDF TS Reduce',izone,nit_ts(izb),ncf_ts(izb),nerr_ts(izb),t(izb),tdel(izb),eta(izb),q(izb)
        EndIf
      EndDo
    EndIf

    ! Undo previous step
    Call bdf_restore(restore)

    ! Reduce order
    Call bdf_adjust_order(kstep)

    ! Scale Nordsieck vector to account for change in stepsize
    Call bdf_rescale(rescale)

    Return
  End Subroutine bdf_check

  Subroutine bdf_correct(kstep)
    !-----------------------------------------------------------------------------------------------
    ! This routine performs various update operations if the solver converges properly;
    ! the corrections are applied to the Nordsieck vector, history are shifted, and previous values
    ! are saved for the next iteration.
    !-----------------------------------------------------------------------------------------------
    Use nuclear_data, Only: aa, nname
    Use xnet_abundances, Only: yt, ydot
    Use xnet_conditions, Only: tdel, t9t, t9dot
    Use xnet_controls, Only: ymin
    Use xnet_integrate, Only: yderiv
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep

    ! Local variables
    Integer :: i, j, izb, izone

    Do izb = 1, nzbatchmx
      If ( ierr_nr(izb) == BDF_STATUS_SUCCESS ) Then

        ! Apply corrections to Nordsieck vector
        z(1:ny,0,izb) = yt(:,izb)
        If ( iheat > 0 ) z(neq,0,izb) = t9t(izb)
        Do j = 1, q(izb)
          z(:,j,izb) = zt0(:,j,izb) + acor(:,izb) * lvec(j,izb)
        EndDo
        !Call dgemm('N','T',neq,q(izb)+1,1,1.0_dp,acor(:,izb),neq,lvec(:,izb),max_order+1,zt0(:,:,izb),neq)
        !Call dger(neq,q(izb)+1,1.0_dp,acor(:,izb),1,lvec(:,izb),1,zt0(:,:,izb),neq)
        !z(:,:,izb) = zt0(:,:,izb)

        ! Shift timestep history
        Do j = max_order, 1, -1
          hvec(j,izb) = hvec(j-1,izb)
        EndDo

        ! Increment step counters
        q_age(izb) = q_age(izb) + 1
        j_age(izb) = j_age(izb) + 1
        p_age(izb) = p_age(izb) + 1
        nscon(izb) = nscon(izb) + 1
      EndIf
    EndDo

    If ( idiag >= 3 ) Then
      Do izb = 1, nzbatchmx
        If ( ierr_nr(izb) == BDF_STATUS_SUCCESS ) Then
          izone = izb + szbatch - 1
          Write(lun_diag,"(a,2i5,i3,1es12.4)") 'BDF Correct',kstep,izone,q(izb),acnrm(izb)
          Write(lun_diag,"(2x,a5,8es15.7)") (nname(i), &
            acor(i,izb),zt0(i,0,izb),z(i,0,izb),yt(i,izb), &
            acor(i,izb)*lvec(1,izb),zt0(i,1,izb),z(i,1,izb),ydot(i,izb)*tdel(izb),i=1,ny)
          If ( iheat > 0 ) Write(lun_diag,"(2x,a5,8es15.7)") 'T9', &
            acor(neq,izb),zt0(neq,0,izb),z(neq,0,izb),t9t(izb), &
            acor(neq,izb)*lvec(1,izb),zt0(neq,1,izb),z(neq,1,izb),t9dot(izb)*tdel(izb)
        EndIf
      EndDo
    EndIf

    Return
  End Subroutine bdf_correct

  Subroutine bdf_eta(kstep)
    !-----------------------------------------------------------------------------------------------
    ! This routine determines the maximum increase in stepsize for the next step by potentially
    ! increasing or decreasing the order. The changes to state variables are made at the beginning
    ! of the next step.
    !-----------------------------------------------------------------------------------------------
    Use xnet_conditions, Only: tdel
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep

    ! Local variables
    Real(dp), Parameter :: biasqm1 = 6.0_dp
    Real(dp), Parameter :: biasq = 6.0_dp
    Real(dp), Parameter :: biasqp1 = 10.0_dp
    Real(dp), Parameter :: addon = 1.0e-6_dp
    Real(dp) :: c, error, acorhat(neq)
    Real(dp) :: alpha0, alpha1, prod, xi, xiold, hsum, a1
    Integer :: i, j, izb, izone

    Do izb = 1, nzbatchmx
      If ( ierr_nr(izb) == BDF_STATUS_SUCCESS ) Then
        izone = izb + szbatch - 1

        ! Only allow change in order of there were no problems in the solve
        If ( eta_next(izb) > 1.0_dp ) Then

          ! Compute eta for each potential change in q
          etaq(:,izb) = 1.0_dp
          error = tq(0,izb) * acnrm(izb)
          etaq(0,izb) = 1.0_dp / ( (biasq * error) ** (1.0_dp/(q(izb)+1)) + addon )
          If ( idiag >= 3 ) Write(lun_diag,"(a,2es23.15)") 'BDF Eta(0)',error,etaq(0,izb)
          If ( q_age(izb) > q(izb) ) Then

            ! eta for q-1
            If ( q(izb) > 1 ) Then
              error = tq(-1,izb) * normw( z(:,q(izb),izb), ewt(:,izb) )
              etaq(-1,izb) = 1.0_dp / ( (biasqm1 * error) ** (1.0_dp/q(izb)) + addon )
              If ( idiag >= 3 ) Write(lun_diag,"(a,2es23.15)") 'BDF Eta(-1)',error,etaq(-1,izb)
            EndIf

            ! eta for q+1
            If ( q(izb) < max_order ) Then
              c = (tq(2,izb) / tq2save(izb)) * (hvec(1,izb) / hvec(2,izb)) ** (q(izb)+1)
              acorhat(:) = acor(:,izb) - c * acorp(:,izb)
              error = tq(1,izb) * normw( acorhat(:), ewt(:,izb) )
              etaq(1,izb) = 1.0_dp / ( (biasqp1 * error) ** (1.0_dp/(q(izb)+2)) + addon )
              If ( idiag >= 3 ) Write(lun_diag,"(a,2es23.15)") 'BDF Eta(+1)',error,etaq(+1,izb)
            EndIf

            ! Wait a bit before considering another order change
            q_age(izb) = q(izb)-1
          EndIf

          ! Choose maximum eta, but prioritize q, then q-1, then q+1
          eta(izb) = 1.0_dp
          If ( etaq(0,izb) > eta(izb) ) Then
            eta(izb) = etaq(0,izb)
            deltaq(izb) = 0
          EndIf
          If ( etaq(-1,izb) > eta(izb) ) Then
            eta(izb) = etaq(-1,izb)
            deltaq(izb) = -1
          EndIf
          If ( etaq(1,izb) > eta(izb) ) Then
            eta(izb) = etaq(1,izb)
            deltaq(izb) = +1
          EndIf

          ! Don't increase stepsize if eta < eta_thresh
          If ( eta(izb) <= eta_thresh ) Then
            eta(izb) = 1.0_dp
            deltaq(izb) = 0
          EndIf
          If ( deltaq(izb) < 0 ) nscon(izb) = 0

          ! Limit stepsize to maximum
          eta(izb) = min( eta(izb), eta_max )

          If ( idiag >= 2 ) Write(lun_diag,"(a,2i5,2i3,sp,i3,ss,4es12.4)") &
            'BDF Eta',kstep,izone,q_age(izb),q(izb),deltaq(izb),eta(izb),(etaq(i,izb),i=-1,1)
        Else

          ! Defer order change if stepsize is stuck
          q_age(izb) = min( q_age(izb), q(izb)-1 )
        EndIf

        ! Save these for the next timestep
        acorp(:,izb) = acor(:,izb)
        tq2save(izb) = tq(2,izb)
      EndIf
    EndDo

    Return
  End Subroutine bdf_eta

  Subroutine bdf_stab(kstep)
    !-----------------------------------------------------------------------------------------------
    ! This routine handles BDF stability limit detection (STALD). If a stability limit violation
    ! is detected, the order is reduced and the step size adjusted accordingly.
    !-----------------------------------------------------------------------------------------------
    Use xnet_util, Only: ifactorial
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep

    ! Local variables
    Real(dp) :: sq, sqm1, sqm2
    Integer :: ldflag(nzbatchmx), fact, i, k, izb, izone
    Logical :: detect(nzbatchmx)

    Do izb = 1, nzbatchmx
      If ( ierr_nr(izb) == BDF_STATUS_SUCCESS ) Then
        If ( q(izb) >= 3 ) Then
          Do k = 1, 3
            Do i = 5, 2, -1
              sddat(i,k,izb) = sddat(i-1,k,izb)
            EndDo
          EndDo
          fact = ifactorial(q(izb)-1)
          sq = fact * q(izb) * (q(izb)+1) * acnrm(izb) / max( tq(2,izb), small )
          sqm1 = fact * q(izb) * normw( z(:,q(izb),izb), ewt(:,izb) )
          sqm2 = fact * normw( z(:,q(izb)-1,izb), ewt(:,izb) )
          sddat(1,1,izb) = sqm2*sqm2
          sddat(1,2,izb) = sqm1*sqm1
          sddat(1,3,izb) = sq*sq
        EndIf
      EndIf
    EndDo

    detect(:) = ( ierr_nr(:) == BDF_STATUS_SUCCESS .and. &
      q(:) >= 3 .and. deltaq(:) >= 0 .and. nscon(:) >= q(:)+5 )
    Call bdf_stab_detect(detect,ldflag)
    Where ( detect(:) )
      Where ( ldflag(:) > 3 )
        deltaq(:) = -1
        eta(:) = min( eta_next(:), etaq(-1,:) )
      ElseWhere ( deltaq(:) /= 0 )
        nscon(:) = 0
      EndWhere
    EndWhere
    If ( idiag >= 2 ) Then
      Do izb = 1, nzbatchmx
        If ( ldflag(izb) > 3 ) Then
          izone = izb + szbatch - 1
          Write(lun_diag,"(a,2i5,i3,es12.4,sp,i3)") &
            'BDF Stab',kstep,izone,q(izb),eta(izb),ldflag(izb)
        EndIf
      EndDo
    EndIf

    Return
  End Subroutine bdf_stab

  Subroutine bdf_stab_detect(mask,kflag)
    !-----------------------------------------------------------------------------------------------
    ! This routine detects stability limit violations from scaled derivative history data. If a
    ! violation is detected (rr > rrcut), kflag > 3 is returned.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Logical, Intent(in) :: mask(:)

    ! Output variables
    Integer, Intent(out) :: kflag(size(mask))

    ! Local variables
    Real(dp), Parameter :: rrcut = 0.98_dp
    Real(dp), Parameter :: vrrtol = 1.0e-4_dp
    Real(dp), Parameter :: vrrt2 = 5.0e-4_dp
    Real(dp), Parameter :: sqtol = 1.0e-3_dp
    Real(dp), Parameter :: rrtol = 1.0e-2_dp
    Real(dp) :: rat(4,3), rav(3), qkr(3), sigsq(3), smin(3), smax(3), ssmax(3), ssdat(5,3)
    Real(dp) :: qp(3), rrc(3), sqmx(3), qjk(3,3), vrat(4), qc(5,3), qco(5,3)
    Real(dp) :: rr, smink, smaxk, sumrat, sumrsq, vmin, vmax, drrmax, adrr
    Real(dp) :: tem, sqmax, saqk, s, sqmaxk, saqj, sqmin
    Real(dp) :: rsa, rsb, rsc, rsd, rd1a, rd1b, rd1c
    Real(dp) :: rd2a, rd2b, rd3a, cest1, corr1
    Real(dp) :: ratp, ratm, qfac1, qfac2, bb, rrb
    Integer :: i, j, k, kmin, it, izb, izone

    kflag(:) = 0
    loop0: Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
        ssdat(:,:) = sddat(:,:,izb)
        rr = 0.0_dp
        Do k = 1, 3
          smink = minval(ssdat(:,k))
          smaxk = maxval(ssdat(:,k))
          If ( smink < small*smaxk ) Then
            kflag(izb) = -1
            Exit
          EndIf
          smin(k) = smink
          smax(k) = smaxk
          ssmax(k) = smaxk*smaxk
          sumrat = 0.0_dp
          sumrsq = 0.0_dp
          Do i = 1, 4
            rat(i,k) = ssdat(i,k) / ssdat(i+1,k)
            sumrat = sumrat + rat(i,k)
            sumrsq = sumrsq + rat(i,k)*rat(i,k)
          EndDo
          rav(k) = sumrat*0.25_dp
          vrat(k) = abs(sumrsq*0.25_dp - rav(k)*rav(k))
          qc(5,k) = ssdat(1,k)*ssdat(3,k) - ssdat(2,k)*ssdat(2,k)
          qc(4,k) = ssdat(2,k)*ssdat(3,k) - ssdat(1,k)*ssdat(4,k)
          qc(3,k) = 0.0_dp
          qc(2,k) = ssdat(2,k)*ssdat(5,k) - ssdat(3,k)*ssdat(4,k)
          qc(1,k) = ssdat(4,k)*ssdat(4,k) - ssdat(3,k)*ssdat(5,k)
          qco(:,k) = qc(:,k)
        EndDo
        vmin = minval(vrat(:))
        vmax = maxval(vrat(:))
        If ( vmin < vrrtol*vrrtol ) Then
          If ( vmax > vrrt2*vrrt2 ) Then
            kflag(izb) = -2
            Cycle loop0
          Else
            rr = sum(rav(:)) / 3.0_dp
            drrmax = maxval( abs(rav(:) - rr) )
            If ( drrmax > vrrt2 ) Then
              kflag(izb) = -3
              Cycle loop0
            Else
              kflag(izb) = 1
            EndIf
          EndIf
        Else
          If ( abs(qco(1,1)) < small*ssmax(1) ) Then
            kflag(izb) = -4
            Cycle loop0
          EndIf
          qco(2:5,2) = qco(2:5,2) - (qco(1,2)/qco(1,1))*qco(2:5,1)
          qco(1,2) = 0.0_dp
          qco(2:5,3) = qco(2:5,3) - (qco(1,3)/qco(1,1))*qco(2:5,1)
          qco(1,3) = 0.0_dp
          If ( abs(qco(2,2)) < small*ssmax(2) ) Then
            kflag(izb) = -4
            Cycle loop0
          EndIf
          qco(3:5,3) = qco(3:5,3) - (qco(2,3)/qco(2,2))*qco(3:5,2)
          If ( abs(qco(4,3)) < small*ssmax(3) ) Then
            kflag(izb) = -4
            Cycle loop0
          EndIf
          rr = -qco(5,3)/qco(4,3)
          If ( rr < small .or. rr > 100.0_dp ) Then
            kflag(izb) = -5
            Cycle loop0
          EndIf
          qkr(:) = qc(5,:) + rr*(qc(4,:) + rr*rr*(qc(2,:) + rr*qc(1,:)))
          sqmax = maxval( abs(qkr(:))/ssmax(:) )
          If ( sqmax < sqtol ) Then
            kflag(izb) = 2
          Else
            Do it = 1, 3
              qp(:) = qc(4,:) + rr*rr*(3.0_dp*qc(2,:) + 4.0_dp*rr*qc(1,:))
              Where ( abs(qp(:)) > small*ssmax(:) )
                rrc(:) = rr - qkr(:) / qp(:)
              ElseWhere
                rrc(:) = rr
              EndWhere
              Do k = 1, 3
                s = rrc(k)
                qjk(:,k) = qc(5,:) + s*(qc(4,:) + s*s*(qc(2,:) + s*qc(1,:)))
                sqmx(k) = maxval( abs(qjk(:,k))/ssmax(:) )
              EndDo
              kmin = minloc(sqmx(:), dim=1)
              sqmin = sqmx(kmin)
              rr = rrc(kmin)
              If ( sqmin < sqtol ) Then
                kflag(izb) = 3
                Exit
              Else
                qkr(:) = qjk(:,kmin)
              EndIf
            EndDo
            If ( sqmin > sqtol ) Then
              kflag(izb) = -6
              Cycle loop0
            EndIf
          EndIf
        EndIf
        Do k = 1, 3
          rsa = ssdat(1,k)
          rsb = ssdat(2,k)*rr
          rsc = ssdat(3,k)*rr*rr
          rsd = ssdat(4,k)*rr*rr*rr
          rd1a = rsa - rsb
          rd1b = rsb - rsc
          rd1c = rsc - rsd
          rd2a = rd1a - rd1b
          rd2b = rd1b - rd1c
          rd3a = rd2a - rd2b
          If ( abs(rd1b) < small*smax(k) ) Then
            kflag(izb) = -7
            Cycle loop0
          EndIf
          cest1 = -rd3a/rd1b
          If ( cest1 < small .or. cest1 > 4.0_dp ) Then
            kflag(izb) = -7
            Cycle loop0
          EndIf
          corr1 = (rd2b/cest1)/(rr*rr)
          sigsq(k) = ssdat(3,k) + corr1
        EndDo
        If ( sigsq(2) < small ) Then
          kflag(izb) = -8
          Cycle loop0
        EndIf
        ratp = sigsq(3) / sigsq(2)
        ratm = sigsq(1) / sigsq(2)
        qfac1 = 0.25_dp * real(q(izb)*q(izb) - 1,dp)
        qfac2 = 2.0_dp / real(q(izb) - 1,dp)
        bb = ratp*ratm - 1.0_dp - qfac1*ratp
        tem = 1.0_dp - qfac2*bb
        If ( abs(tem) < small ) Then
          kflag(izb) = -8
          Cycle loop0
        EndIf
        rrb = 1.0_dp / tem
        If ( abs(rrb - rr) > rrtol ) Then
          kflag(izb) = -9
          Cycle loop0
        EndIf
        If ( rr > rrcut ) Then
          If ( kflag(izb) == 1 ) kflag(izb) = 4
          If ( kflag(izb) == 2 ) kflag(izb) = 5
          If ( kflag(izb) == 3 ) kflag(izb) = 6
        EndIf
      EndIf
    EndDo loop0

    Return
  End Subroutine bdf_stab_detect

  Subroutine bdf_adjust_order(kstep)
    !-----------------------------------------------------------------------------------------------
    ! This subroutine restores the Nordsieck vector back to its state before the prediction.
    !-----------------------------------------------------------------------------------------------
    Use xnet_conditions, Only: tdel
    Implicit None

    ! Input variables
    Integer, Intent(in) :: kstep

    ! Local variables
    Real(dp) :: alpha0, alpha1, prod, xi, xiold, hsum, a1
    Integer :: i, j, izb, izone

    Do izb = 1, nzbatchmx

      If ( deltaq(izb) == -1 ) Then

        ! Decrease order
        lvec(:,izb) = 0.0_dp
        lvec(2,izb) = 1.0_dp
        hsum = 0.0_dp
        Do j = 1, q(izb)-2
          hsum = hsum + hvec(j,izb)
          xi = hsum / tdel(izb)
          Do i = j+2, 2, -1
            lvec(i,izb) = lvec(i,izb)*xi + lvec(i-1,izb)
          EndDo
        EndDo
        Do j = 2, q(izb)-1
          z(:,j,izb) = z(:,j,izb) - z(:,q(izb),izb) * lvec(j,izb)
        EndDo
        !Call dgemm('N','T',neq,q(izb)-2,1,-1.0_dp,z(:,q(izb),izb),neq,lvec(2,izb),max_order+1,z(1,2,izb),neq)
        !Call dger(neq,q(izb)-2,-1.0_dp,z(:,q(izb),izb),1,lvec(2,izb),1,z(1,2,izb),neq)
        z(:,q(izb),izb) = 0.0_dp
        q(izb) = q(izb) - 1
        q_age(izb) = 0

      ElseIf ( deltaq(izb) == +1 ) Then

        ! Increase order
        lvec(:,izb) = 0.0_dp
        lvec(2,izb) = 1.0_dp
        alpha0 = -1.0_dp
        alpha1 = 1.0_dp
        prod = 1.0_dp
        xiold = 1.0_dp
        hsum = tdel(izb)
        Do j = 1, q(izb)-1
          hsum = hsum + hvec(j+1,izb)
          xi = hsum / tdel(izb)
          prod = prod * xi
          alpha0 = alpha0 - 1.0_dp / real(j+1,dp)
          alpha1 = alpha1 + 1.0_dp / xi
          Do i = j+2, 2, -1
            lvec(i,izb) = lvec(i,izb)*xiold + lvec(i-1,izb)
          EndDo
          xiold = xi
        EndDo
        a1 = -(alpha0 + alpha1) / prod
        z(:,q(izb)+1,izb) = a1 * acor(:,izb)
        Do j = 2, q(izb)
          z(:,j,izb) = z(:,j,izb) + a1 * acor(:,izb) * lvec(j,izb)
        EndDo
        !Call dgemm('N','T',neq,q(izb),1,a1,acor(:,izb),neq,lvec(2,izb),max_order+1,z(1,2,izb),neq)
        !Call dger(neq,q(izb),a1,acor(:,izb),1,lvec(2,izb),1,z(1,2,izb),neq)
        q(izb) = q(izb) + 1
        q_age(izb) = 0
      EndIf
    EndDo
    If ( idiag >= 2 ) Then
      Do izb = 1, nzbatchmx
        izone = izb + szbatch - 1
        If ( deltaq(izb) /= 0 ) Write(lun_diag,"(a,2i5,i3,sp,i3)") &
          'BDF Order Change',kstep,izone,q(izb),deltaq(izb)
      EndDo
    EndIf
    deltaq(:) = 0

    Return
  End Subroutine bdf_adjust_order

  Subroutine bdf_rescale(mask)
    !-----------------------------------------------------------------------------------------------
    ! This subroutine rescales the Nordsieck vector and timestep for a new eta.
    !-----------------------------------------------------------------------------------------------
    Use xnet_conditions, Only: tdel, t, tt
    Implicit None

    ! Input variables
    Logical, Intent(in) :: mask(:)

    ! Local variables
    Integer :: i, j, izb

    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
        tdel(izb) = eta(izb) * dt_scale(izb)
        hvec(0,izb) = tdel(izb)
        tt(izb) = t(izb) + tdel(izb)
        ! Scale Nordsieck vector to account for change in stepsize
        Do j = 1, q(izb)
          z(:,j,izb) = z(:,j,izb) * eta(izb)**j
        EndDo
        dt_scale(izb) = tdel(izb)
        nscon(izb) = 0
      EndIf
    EndDo

    Return
  End Subroutine bdf_rescale

  Subroutine bdf_restore(mask)
    !-----------------------------------------------------------------------------------------------
    ! This subroutine restores the Nordsieck vector back to its state before the prediction.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Logical, Intent(in) :: mask(:)

    ! Local variables
    Integer :: i, j, izb

    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
        z(:,:,izb) = z0(:,:,izb)
      EndIf
    EndDo

    Return
  End Subroutine bdf_restore

  Subroutine pascal_build
    !-----------------------------------------------------------------------------------------------
    ! Construct the triangular Pascal matrix:
    !   $A^{ij}(q) = \binom{i}{j} = \frac{i!}{j! (i-j)!), i \ge j; i, j = 0, 1, 2, ..., q$
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Local variables
    Integer :: i, j

    ! Initialize everything to zero
    A_pascal(:,:) = 0
    Ainv_pascal(:,:) = 0

    ! Set the first column
    A_pascal(0:max_order,0) = 1

    ! Build the matrix
    Do i = 1, max_order
      Do j = 1, i
        A_pascal(i,j) = A_pascal(i-1,j) + A_pascal(i-1,j-1)
        If ( mod(i-j,2) == 0 ) Then
          Ainv_pascal(i,j) = A_pascal(i,j)
        Else
          Ainv_pascal(i,j) = -A_pascal(i,j)
        EndIf
      EndDo
    EndDo

    Return
  End Subroutine pascal_build

  Subroutine error_weights( y, wt )
    !-----------------------------------------------------------------------------------------------
    ! This routine returns the error weights for comparing errors in y
    !-----------------------------------------------------------------------------------------------
    Use xnet_controls, Only: iconvc, ymin
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: y(:)

    ! Output variables
    Real(dp), Intent(out) :: wt(size(y))

    ! Local variables
    Integer :: i, n

    n = size(y)
    If ( iconvc == 0 .or. iconvc == 1 ) Then
      Where ( y(:) < ymin )
        wt(:) = 0.0_dp
      ElseWhere ( y(:) < atol(:) )
        wt(:) = 1.0_dp / (rtol(:) * atol(:))
      ElseWhere
        wt(:) = 1.0_dp / (rtol(:) * abs(y(:)))
      EndWhere
    Else
      wt(:) = 1.0_dp / ( rtol(:) * max(abs(y(:)),ymin) + atol(:) )
    EndIf

    Return
  End Subroutine error_weights

  Function normw( x, wt ) Result( xnrm )
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the weighted l-norm of vector x, with l deteremined from the value of
    ! the iconvc control:
    !   iconvc = 0 : infinity norm
    !   iconvc = 1 : 1-norm
    !   iconvc = 2 : 2-norm
    !   iconvc = 3 : RMS norm
    !-----------------------------------------------------------------------------------------------
    Use xnet_controls, Only: iconvc
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: x(1:), wt(1:)

    ! Function variable
    Real(dp) :: xnrm

    ! Local variables
    Real(dp) :: xw(size(x))
    Integer :: n, incx

    n = size(x)
    xw(:) = x(:) * wt(:)
    If ( iconvc == 0 ) Then
      xnrm = maxval( abs(xw(:)) )
    ElseIf ( iconvc == 1 ) Then
      xnrm = sum( abs(xw(:)) )
    ElseIf ( iconvc == 2 ) Then
      xnrm = sqrt( sum(xw(:)**2) )
    ElseIf ( iconvc == 3 ) Then
      xnrm = sqrt( sum(xw(:)**2) / real(n,dp) )
    Else
      xnrm = 0.0_dp
    EndIf

    Return
  End Function normw

End Module xnet_integrate_bdf