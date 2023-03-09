Module xnet_nse
  !-------------------------------------------------------------------------------------------------
  ! Fast, thread-safe solver for composition in nuclear statistical equilibrium (NSE). Solution is
  ! obtained by solving the Saha equation using globally convergent Newton method with backtracking
  ! line search. Corrections from Coulomb interactions (screening) and nuclear partition functions
  ! are included. For theoretical review, see Seitenzahl et al., ADNT 95 (2009) 96. For numerical
  ! review, see Dennis & Schnabel (1996) Numerical Methods for Unconstrained Optimization and
  ! Nonlinear Equations (doi.org/10.1137/1.9781611971200).
  !-------------------------------------------------------------------------------------------------
  Use nuclear_data, Only: ny, aa, nn, zz, ia, iz, izmax, be, mex, mm, nname, angm, g, ng, t9i, &
    & zseq, zseq53, zseqi
  Use xnet_constants, Only: avn, bi, bip1, bok, cds, clt, e2, epmev, hbar, kbi, m_n, m_p, &
    & m_e, pi, third
  Use xnet_controls, Only: idiag, itsout, iscrn, lun_diag, lun_stderr, lun_stdout
  Use xnet_timers
  Use xnet_types, Only: dp
  Use xnet_util, Only: xnet_terminate, safe_exp
  Implicit None
  Private

  ! NSE solver parameters
  Integer, Parameter  :: lnorm    = 0       ! Which lnorm to use for determining convergence
  Real(dp), Parameter :: tolu     = 1.0e-14_dp ! Tolerance for testu convergence test
  Real(dp), Parameter :: tolf     = 1.0e-8_dp  ! Tolerance for testf convergence test
  Real(dp), Parameter :: tolmin   = 1.0e-14_dp ! Tolerance for minimum change per iteration
  Real(dp), Parameter :: tolbeta  = 3.0e-15_dp ! Tolerance for minimum allowable "step" length
  Real(dp), Parameter :: maxstep0 = 2.0_dp     ! Scaling factor for maximum relative "step" length
  Real(dp), Parameter :: typf     = 1.0_dp

  Real(dp), Parameter :: alpha  = 1.0e-8_dp
  Real(dp), Parameter :: gamma1 = 0.1_dp
  Real(dp), Parameter :: gamma2 = 0.5_dp

  Integer, Parameter :: nritmax = 200 ! Maximum number of Newton-Raphson iterations
  Integer, Parameter :: lsitmax = 30  ! Maximum number of line-search iterations

  ! NSE solution variables
  Real(dp), Allocatable, Public :: unse(:) ! Chemical potentials for each species
  Real(dp), Allocatable, Public :: xnse(:) ! Mass fractions
  Real(dp), Allocatable, Public :: ynse(:) ! Molar fractions
  !$omp threadprivate(xnse,ynse,unse)

  ! Other NSE state variables
  Real(dp) :: rhonse, t9nse, yense
  Real(dp), Allocatable :: mm52(:), ggnse(:)
  !$omp threadprivate(rhonse,t9nse,yense,ggnse)

  ! These contain the solution vector and Jacobian of the NSE linear system
  Real(dp) :: fvec(2), fjac(2,2)
  !$omp threadprivate(fvec,fjac)

  ! Scaling variables
  Real(dp) :: typu(2), typfvec(2)
  Real(dp) :: scaleu(2), scalefvec(2)
  !$omp threadprivate(typu,typfvec,scaleu,scalefvec)

  ! Counters
  Integer :: knrtot(3), knr(3)
  Logical :: use_CP98
  !$omp threadprivate(knrtot,knr,use_CP98)

  ! Screening variables
  Integer, Allocatable :: intz(:)
  Real(dp), Allocatable :: hnse(:)
  !$omp threadprivate(hnse)

  ! Network indices for some characteristic nuclei
  Integer :: ibe                           ! Maximally bound nucleus
  Integer :: iza_min, iza_max              ! min/max Z/A (A>1)
  Integer :: iye_za, iye_za_l, iye_za_u    ! closest Z/A to Ye from above and below (A>1)
  Integer :: ihvy_min, ihvy_max            ! min/max Z/A (Z>=20)
  Integer :: iye_hvy, iye_hvy_l, iye_hvy_u ! closest Z/A to Ye from abbove and below (Z>=20)
  Integer :: ife_min, ife_max              ! min/max Z/A (common Fe/Ni isotopes)
  Integer :: iye_fe, iye_fe_l, iye_fe_u    ! closest Z/A to Ye from above and below (common Fe/Ni isotopes)
  Integer :: i_nn, i_pp, i_he4, i_si28, i_ni56
  Real(dp) :: zatst, za_min, za_max
  !$omp threadprivate(iye_za,iye_za_l,iye_za_u,iye_hvy,iye_hvy_l,iye_hvy_u,iye_fe,iye_fe_l,iye_fe_u,zatst)

  Public :: nse_initialize, nse_solve

Contains

  Function l2norm( x ) Result( y )
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the l2norm of vector x
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: x(:)

    ! Function variable
    Real(dp) :: y

    ! Local variables
    Integer :: n, incx

    ! External function
    Real(dp), External :: dnrm2

    n = size(x)
    incx = 1

    y = dnrm2( n, x, incx )

    Return
  End Function l2norm

  Function dot_prod( x, y ) Result( z )
    !-----------------------------------------------------------------------------------------------
    ! This routines calculates the dot product of vectors x and y.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: x(:)
    Real(dp), Intent(in) :: y(:)

    ! Function variable
    Real(dp) :: z

    ! Local variables
    Integer :: n, incx, incy

    ! External function
    Real(dp), External :: ddot

    n = size(x)
    incx = 1
    incy = 1

    z = ddot( n, x, incx, y, incy )

    Return
  End Function dot_prod

  Subroutine nse_initialize
    !-----------------------------------------------------------------------------------------------
    ! This routine allocates and initializes NSE arrays.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Local variables
    Integer :: i

    start_timer = xnet_wtime()
    timer_nse = timer_nse - start_timer

    !$omp parallel default(shared)
    ! Allocate NSE composition arrays
    Allocate (xnse(ny),ynse(ny),unse(ny))
    ! Allocate NSE partition function array
    Allocate (ggnse(ny))
    ! Allocate NSE screening arrays
    Allocate (hnse(izmax))

    ! Typical relative values for u and f(u) (used for scaling)
    typu(:)    = (/ 1.0_dp, 1.0_dp /)
    typfvec(:) = (/ 1.0_dp, 1.0_dp /)
    scaleu(:)    = 1.0_dp / typu(:)
    scalefvec(:) = 1.0_dp / typfvec(:)
    !$omp end parallel

    ! Pre-calculate m_i^{5/2}, which is used throughout the solver
    Allocate (mm52(ny))
    mm52(:) = mm(:) * mm(:) * sqrt(mm(:))

    Allocate (intz(ny))
    Where ( zz(:) > 1.0_dp )
      intz(:) = nint(zz(:))
    ElseWhere
      intz(:) = 1
    EndWhere

    ! Find maximally bound nucleus
    ibe = maxloc( be(:)/aa(:), 1 )

    ! Find min/max Z/A in network (A>1)
    iza_min = minloc( zz(:)/aa(:), 1, mask = ia(:) > 1 )
    iza_max = maxloc( zz(:)/aa(:), 1, mask = ia(:) > 1 )
    za_min = zz(iza_min)/aa(iza_min)
    za_max = zz(iza_max)/aa(iza_max)

    ! Find min/max Z/A for Z>=20 nuclei
    ihvy_min = minloc( zz(:)/aa(:), 1, mask = iz(:) >= 20 )
    ihvy_max = ny + 1 - maxloc( zz(ny:1:-1)/aa(ny:1:-1), 1, mask = iz(ny:1:-1) >= 20 )

    ! Find min/max Z/A for common Fe/Ni nuclei
    ife_min = minloc( zz(:)/aa(:), 1, ( iz(:) == 26 .or. iz(:) == 28 ) .and. &
      &                               ( ia(:) == 52 .or. ia(:) == 54 .or. &
      &                                 ia(:) == 56 .or. ia(:) == 58 ) )
    ife_max = ny + 1 - maxloc( zz(ny:1:-1)/aa(ny:1:-1), 1, ( iz(ny:1:-1) == 26 .or. iz(ny:1:-1) == 28) .and. &
      &                                                    ( ia(ny:1:-1) == 52 .or. ia(ny:1:-1) == 54 .or. &
      &                                                      ia(ny:1:-1) == 56 .or. ia(ny:1:-1) == 58) )

    i_nn = 0 ; i_pp = 0 ; i_he4 = 0 ; i_si28 = 0 ; i_ni56 = 0
    Do i = 1, ny
      If ( iz(i) == 0  .and. ia(i) == 1  ) i_nn = i
      If ( iz(i) == 1  .and. ia(i) == 1  ) i_pp = i
      If ( iz(i) == 2  .and. ia(i) == 4  ) i_he4 = i
      If ( iz(i) == 14 .and. ia(i) == 28 ) i_si28 = i
      If ( iz(i) == 28 .and. ia(i) == 56 ) i_ni56 = i
    EndDo
    If ( i_nn == 0 ) Then
      Write(lun_stdout,'(2(a,a5))') 'Could not find ','n','; setting i_nn to min Z/A: ',nname(iza_min)
      i_nn = iza_min
    EndIf
    If ( i_pp == 0 ) Then
      Write(lun_stdout,'(2(a,a5))') 'Could not find ','p','; setting i_pp to max Z/A: ',nname(iza_max)
      i_pp = iza_max
    EndIf
    If ( i_he4 == 0 ) Then
      Write(lun_stdout,'(2(a,a5))') 'Could not find ','he4','; setting i_he4 to closest Z/A to Ye: ',nname(iye_za)
      i_he4 = iye_za
    EndIf
    If ( i_si28 == 0 ) Then
      Write(lun_stdout,'(2(a,a5))') 'Could not find ','si28','; setting i_si28 to closest heavy Z/A to Ye: ',nname(iye_hvy)
      i_si28 = iye_hvy
    EndIf
    If ( i_ni56 == 0 ) Then
      Write(lun_stdout,'(2(a,a5))') 'Could not find ','ni56','; setting i_ni56 to max B/A: ',nname(ibe)
      i_ni56 = ibe
    EndIf

    stop_timer = xnet_wtime()
    timer_nse = timer_nse + stop_timer

    Return
  End Subroutine nse_initialize

  Subroutine nse_inuc
    !-----------------------------------------------------------------------------------------------
    ! This routine determines the species indiecs for some characteristic nuclei to be used in the
    ! nse_guess routine.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Local variables
    Integer :: i

    ! Find A>1 nuclei with Z/A closest to Ye from above and below
    iye_za_u = ny + 1 - minloc( abs(zz(ny:1:-1)/aa(ny:1:-1) - yense), 1, &
      &                         mask = zz(ny:1:-1)/aa(ny:1:-1) >= yense .and. &
      &                                ia(ny:1:-1) > 1 )
    If ( iye_za_u > ny ) Then
      iye_za_u = iza_max
    EndIf
    iye_za_l = minloc( abs(zz(:)/aa(:) - yense), 1, &
      &                mask = zz(:)/aa(:) <= yense .and. &
      &                       ia(:) > 1 )
    If ( iye_za_l < 1 ) Then
      iye_za_l = iza_min
    EndIf
    zatst = zz(iye_za_l)/aa(iye_za_l)
    If ( abs(zatst - yense) < abs(zz(iye_za_u)/aa(iye_za_u) - yense) ) Then
      iye_za = iye_za_l
    Else
      iye_za = iye_za_u
    EndIf

    ! Find Z>=20 nuclei with Z/A closest to Ye from above and below
    iye_hvy_u = ny + 1 - minloc( abs(zz(ny:1:-1)/aa(ny:1:-1) - yense), 1, &
      &                         mask = zz(ny:1:-1)/aa(ny:1:-1) >= yense .and. &
      &                                iz(ny:1:-1) >= 20 )
    If ( iye_hvy_u > ny ) Then
      iye_hvy_u = ihvy_max
    EndIf
    iye_hvy_l = minloc( abs(zz(:)/aa(:) - yense), 1, &
      &                mask = zz(:)/aa(:) <= yense .and. &
      &                       iz(:) >= 20 )
    If ( iye_hvy_l < 1 ) Then
      iye_hvy_l = ihvy_min
    EndIf
    zatst = zz(iye_hvy_l)/aa(iye_hvy_l)
    If ( abs(zatst - yense) < abs(zz(iye_hvy_u)/aa(iye_hvy_u) - yense) ) Then
      iye_hvy = iye_hvy_l
    Else
      iye_hvy = iye_hvy_u
    EndIf

    ! Find common Fe/Ni nuclei with Z/A closest to Ye from above and below
    iye_fe_u = ny + 1 - minloc( abs(zz(ny:1:-1)/aa(ny:1:-1) - yense), 1, &
      &                         mask = zz(ny:1:-1)/aa(ny:1:-1) >= yense .and. &
      &                                ( iz(ny:1:-1) == 26 .or. iz(ny:1:-1) == 28 ) .and. &
      &                                ( ia(ny:1:-1) == 52 .or. ia(ny:1:-1) == 54 .or. &
      &                                  ia(ny:1:-1) == 56 .or. ia(ny:1:-1) == 58 ) )
    If ( iye_fe_u > ny ) Then
      iye_fe_u = ife_max
    EndIf
    iye_fe_l = minloc( abs(zz(:)/aa(:) - yense), 1, &
    &                  mask = zz(:)/aa(:) <= yense .and. &
    &                         ( iz(:) == 26 .or. iz(:) == 28 ) .and. &
    &                         ( ia(:) == 52 .or. ia(:) == 54 .or. &
    &                           ia(:) == 56 .or. ia(:) == 58 ) )
    If ( iye_fe_l < 1 ) Then
      iye_fe_l = ihvy_min
    EndIf
    zatst = zz(iye_fe_l)/aa(iye_fe_l)
    If ( abs(zatst - yense) < abs(zz(iye_fe_u)/aa(iye_fe_u) - yense) ) Then
      iye_fe = iye_fe_l
    Else
      iye_fe = iye_fe_u
    EndIf

    Return
  End Subroutine nse_inuc

  Subroutine nse_solve(rho,t9,ye)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the NSE composition by solving for the neutron and proton chemical
    ! potentials which satisfy the Saha equation under the constraints of mass and charge
    ! conservation. The numerical method employed is largely adapated from Algorithm D6.1.3 of:
    !     Dennis, J. & Schnabel, R. (1996) Numerical Methods for Unconstrained Optimization and
    !       Nonlinear Equations, doi.org/10.1137/1.9781611971200
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: rho, t9, ye

    ! Local variables
    Real(dp) :: uvec(2), uvec0(2)
    Integer :: i, ii
    Logical :: check
    Integer :: info, info0, iguess, iguess0

    start_timer = xnet_wtime()
    timer_nse = timer_nse - start_timer

    ! Reset counters
    knrtot(:) = 0

    ! Initialize NSE variables
    rhonse = rho
    t9nse = t9
    yense = ye
    Call nse_pf

    ! Identify species in network that are close in Z/A to Ye
    Call nse_inuc

    ! Start with easy screening to get a good guess
    use_CP98 = .true.
    Call nse_screen_CP98
    Do iguess = 1, 6
      unse(:) = 0.0_dp
      xnse(:) = 0.0_dp
      ynse(:) = 0.0_dp
      Call nse_guess(iguess,uvec)
      Call nse_nr(uvec,check,info)
      If ( info > 0 ) Exit
    EndDo
    If ( itsout >= 3 ) Write(lun_stdout,'(a5,2i2,8es23.15)') 'CP98:', &
    & info,iguess,uvec(1),uvec(2),xnse(i_nn),xnse(i_pp),xnse(i_he4),xnse(i_ni56),fvec(1),fvec(2)

    ! Save the CP98 result
    iguess0 = iguess
    info0 = info
    uvec0(:) = uvec(:)

    ! Now use XNet's screening to match network
    If ( iscrn > 0 ) Then
      use_CP98 = .false.
      If ( info > 0 ) Then
        iguess = 0

        ! Iterate a few times to make sure screening and composition are consistent
        Do i = 1, 3
          Call nse_screen
          Call nse_nr(uvec,check,info)
          If ( itsout >= 3 ) Write(lun_stdout,'(a5,2i2,8es23.15)') 'XNET:', &
          & info,iguess,uvec(1),uvec(2),xnse(i_nn),xnse(i_pp),xnse(i_he4),xnse(i_ni56),fvec(1),fvec(2)
        EndDo
      EndIf

      ! Try different guesses if both screening approaches fail
      If ( info <= 0 ) Then
        Do iguess = 1, 6
          unse(:) = 0.0_dp
          xnse(:) = 0.0_dp
          ynse(:) = 0.0_dp
          hnse(:) = 0.0_dp
          Call nse_guess(iguess,uvec)
          Call nse_screen
          Call nse_nr(uvec,check,info)
          If ( info > 0 ) Exit
        EndDo
        If ( itsout >= 3 ) Write(lun_stdout,'(a5,2i2,8es23.15)') 'XNET:', &
        & info,iguess,uvec(1),uvec(2),xnse(i_nn),xnse(i_pp),xnse(i_he4),xnse(i_ni56),fvec(1),fvec(2)
      EndIf
    EndIf

    If ( info <= 0 ) Then
      If ( itsout >= 3 ) Then
        If ( info == 0 ) Then
          Write(lun_stdout,'(a)') 'Exceeded max NR iterations'
        ElseIf ( info == -1 ) Then
          Write(lun_stdout,'(a)') 'Exceeded max LS iterations'
        ElseIf ( info == -2 ) Then
          Write(lun_stdout,'(a)') 'Failed to find sufficient decrease in LS'
        ElseIf ( info == -3 ) Then
          Write(lun_stdout,'(a)') 'Could not converge to solution'
        ElseIf ( info == -4 ) Then
          Write(lun_stdout,'(a)') 'Non-negative slope in LS'
        ElseIf ( info == -6 ) Then
          Write(lun_stdout,'(a)') 'Minimizer is not a root'
        ElseIf ( info == -7 ) Then
          Write(lun_stdout,'(a)') 'Singular matrix'
        EndIf
      EndIf
      If ( info0 > 0 ) Then
        If ( itsout >= 3 ) Write(lun_stdout,'(a)') 'XNet screening failed; using CP98 screening result'
      Else
        Write(lun_stdout,'(a,i2,a,3es23.15)') &
          & 'Check convergence of root finder (',info,'), rho,T,Ye= ',rhonse,t9nse,yense
        Call xnet_terminate('NSE ERROR: NR Failed',info)
      EndIf
    EndIf

    ! Evaluate NSE composition with converged solution
    If ( info > 0 ) Then
      use_CP98 = .false.
      Call nse_eval(uvec)
      knrtot(3) = knrtot(3) + 1

    ! If XNet screening failed to converge, use the CP98 result
    ElseIf ( info0 > 0 ) Then
      use_CP98 = .true.
      uvec(:) = uvec0(:)
      Call nse_screen_CP98
      Call nse_eval(uvec)
      knrtot(3) = knrtot(3) + 1
    EndIf

    If ( itsout >= 2 ) Then
      Write(lun_stdout,'(a)') 'NSE solved'
      Write(lun_stdout,'(a,1x,3i15)')     'NR iters, LS iters, F evals =',knrtot(1),knrtot(2),knrtot(3)
      Write(lun_stdout,'(a,2es23.15)')    'roots                       =',uvec(1),uvec(2)
      Write(lun_stdout,'(a,4es23.15)')    'xneut, xprot, xalpha, xni56 =',xnse(i_nn),xnse(i_pp),xnse(i_he4),xnse(i_ni56)
      Write(lun_stdout,'(a,2es23.15)')    'residuals                   =',fvec(1),fvec(2)
      ii = maxloc(xnse,1)
      Write(lun_stdout,'(a,a23,es23.15)') 'largest mass fraction       =',nname(ii),xnse(ii)
    EndIf

    stop_timer = xnet_wtime()
    timer_nse = timer_nse + stop_timer

    Return
  End Subroutine nse_solve

  Subroutine nse_guess(iguess,uvec)
    !-----------------------------------------------------------------------------------------------
    ! This routine provides the initial guess for the neutron and chemical potentials. The guess is
    ! determined by assuming a composition of only two nuclei in the network and solving the inverse
    ! Saha equation, with mass fractions chosen such that they reproduce the electron fraction. A
    ! good guess can drastically reduce the number of iterations required and, in some cases, is
    ! necessary to converge at all. Different initial guesses are returned depending on the value
    ! of iguess. The order of these guesses reflects their predicted likelihood to quickly converge:
    !
    ! iguess = 1 and low density and high temperature, assume free nucleons:
    !               iye(1) = n                          ; iye(2) = p
    ! iguess > 6,
    !               iye(1) = n                          ; iye(2) = p
    ! otherwise...
    ! if Ye < min(Z/A),
    !   iguess = 1: iye(1) = heavy w/ closest Z/A <= Ye ; iye(2) = n
    !   iguess = 2: iye(1) = ni56                       ; iye(2) = n
    !   iguess = 3: iye(1) = ni56
    !   iguess = 4: iye(1) = heavy w/ closest Z/A <= Ye
    !   iguess = 5: iye(1) = max(B/A)                   ; iye(2) = n
    !   iguess = 6: iye(1) = max(B/A)
    ! if min(Z/A) <= Ye <= ~0.5,
    !   iguess = 1: iye(1) = Fe/Ni w/ closest Z/A <= Ye ; iye(2) = Fe/Ni w/ closest Z/A >  Ye
    !   iguess = 2: iye(1) = ni56                       ; iye(2) = Fe/Ni w/ closest Z/A <= Ye
    !   iguess = 3: iye(1) = he4                        ; iye(2) = Fe/Ni w/ closest Z/A <= Ye
    !   iguess = 4: iye(1) = ni56
    !   iguess = 5: iye(1) = heavy w/ closest Z/A <= Ye ; iye(2) = heavy w/ closest Z/A >  Ye
    !   iguess = 6: iye(1) = max(B/A)
    ! if Ye ~= 0.5
    !   iguess = 1: iye(1) = ni56
    !   iguess = 2: iye(1) = Fe/Ni w/ closest Z/A <= Ye
    !   iguess = 3: iye(1) = max(B/A)
    !   iguess = 4: iye(1) = heavy w/ closest Z/A >  Ye
    !   iguess = 5: iye(1) = si28
    !   iguess = 6: iye(1) = he4
    ! if ~0.5 < Ye
    !   iguess = 1: iye(1) = ni56                       ; iye(2) = p
    !   iguess = 2: iye(1) = ni56
    !   iguess = 3: iye(1) = he4                        ; iye(2) = closest Z/A >  Ye
    !   iguess = 4: iye(1) = ni56                       ; iye(2) = closest Z/A >  Ye
    !   iguess = 5: iye(1) = max(B/A)
    !   iguess = 6: iye(1) = max(B/A)                   ; iye(2) = p
    ! based on prior experience.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Integer, Intent(in) :: iguess

    ! Output variables
    Real(dp), Intent(out) :: uvec(2)

    ! Local variables
    Real(dp), Parameter :: xmax = 1.0_dp, xmin = 0.0_dp
    Real(dp) :: bkt, c1, c2, det, lhs(2,2), rhs(2)
    Real(dp) :: za_ye(2)
    Integer :: i, ii, iye(2)

    If ( iguess <= 1 .and. ( (rhonse <= 1.0e5_dp  .and. t9nse >= 7.0_dp ) .or. &
      &                      (rhonse <= 1.0e6_dp  .and. t9nse >= 8.0_dp ) .or. &
      &                      (rhonse <= 1.0e7_dp  .and. t9nse >= 9.0_dp ) .or. &
      &                      (rhonse <= 1.0e8_dp  .and. t9nse >= 11.0_dp) .or. &
      &                      (rhonse <= 1.0e9_dp  .and. t9nse >= 14.0_dp) .or. &
      &                      (rhonse <= 1.0e10_dp .and. t9nse >= 17.0_dp) .or. &
      &                      (                          t9nse >= 20.0_dp) ) ) Then
      iye(1) = i_nn
      iye(2) = i_pp
    ElseIf ( yense < za_min ) Then
      If ( iguess <= 1 ) Then
        iye(1) = iye_hvy_l
        iye(2) = i_nn
      ElseIf ( iguess == 2 ) Then
        iye(1) = i_ni56
        iye(2) = i_nn
      ElseIf ( iguess == 3 ) Then
        iye(:) = i_ni56
      ElseIf ( iguess == 4 ) Then
        iye(:) = iye_hvy_l
      ElseIf ( iguess == 5 ) Then
        iye(1) = ibe
        iye(2) = i_nn
      ElseIf ( iguess == 6 ) Then
        iye(:) = ibe
      Else
        iye(1) = i_nn
        iye(2) = i_pp
      EndIf
    ElseIf ( yense >= za_min .and. yense < 0.5_dp ) Then
      If ( iguess <= 1 ) Then
        iye(1) = iye_fe_l
        iye(2) = iye_fe_u
      ElseIf ( iguess == 2 ) Then
        iye(1) = i_ni56
        iye(2) = iye_fe_l
      ElseIf ( iguess == 3 ) Then
        iye(1) = i_he4
        iye(2) = iye_fe_l
      ElseIf ( iguess == 4 ) Then
        iye(:) = i_ni56
      ElseIf ( iguess == 5 ) Then
        iye(1) = iye_hvy_l
        iye(2) = iye_hvy_u
      ElseIf ( iguess == 6 ) Then
        iye(:) = ibe
      Else
        iye(1) = i_nn
        iye(2) = i_pp
      EndIf
    ElseIf ( yense >= 0.501_dp ) Then
      If ( iguess <= 1 ) Then
        iye(1) = i_ni56
        iye(2) = i_pp
      ElseIf ( iguess == 2 ) Then
        iye(:) = i_ni56
      ElseIf ( iguess == 3 ) Then
        iye(1) = i_he4
        iye(2) = i_pp
      ElseIf ( iguess == 4 ) Then
        iye(1) = i_ni56
        iye(2) = iye_za_u
      ElseIf ( iguess == 5 ) Then
        iye(:) = ibe
      ElseIf ( iguess == 6 ) Then
        iye(1) = ibe
        iye(2) = i_pp
      Else
        iye(1) = i_nn
        iye(2) = i_pp
      EndIf
    Else
      If ( iguess <= 1 ) Then
        iye(:) = i_ni56
      ElseIf ( iguess == 2 ) Then
        iye(:) = iye_fe_l
      ElseIf ( iguess == 3 ) Then
        iye(:) = ibe
      ElseIf ( iguess == 4 ) Then
        iye(:) = iye_hvy_u
      ElseIf ( iguess == 5 ) Then
        iye(:) = i_si28
      ElseIf ( iguess == 6 ) Then
        iye(:) = i_he4
      Else
        iye(1) = i_nn
        iye(2) = i_pp
      EndIf
    EndIf
    za_ye(:) = zz(iye(:))/aa(iye(:))

    ! Calculate mass fractions consistent with Ye and nuclei above
    If ( iguess /= 0 ) Then
      If ( iye(1) == iye(2) ) Then
        xnse(iye(1)) = xmax
      ElseIf ( abs(za_ye(1) - za_ye(2)) < tiny(0.0_dp) ) Then
        xnse(iye(1)) = xmax
        xnse(iye(2)) = xmin
      Else
        xnse(iye(1)) = max( min( ( yense - za_ye(2) ) / ( za_ye(1) - za_ye(2) ), xmax ), xmin )
        xnse(iye(2)) = max( min( 1.0_dp - xnse(iye(1)), xmax ), xmin )
      EndIf
      ynse(iye) = xnse(iye) / (mm(iye)*avn)

    ! If iguess = 0, assume we have a full composition, so update the guess by picking the
    ! most abundant species with Z/A closest to Ye (i.e. the least wrong nuclei)
    Else
      iye(:) = minloc( (zz(:)*xnse(:)/aa(:) - yense)**2 + (xnse(:)-1.0_dp)**2, 1 )
!     iye(1) = i_nn
!     iye(2) = i_pp
    EndIf

    ! Update screening corrections
    Call nse_screen

    ! Calculate the initial guess (un and up)
    bkt = t9nse*bok*epmev
    c1 = bkt / (2.0_dp*pi*hbar*hbar*epmev*epmev)

    ! Assume un = up if only using 1 species
    If ( iye(1) == iye(2) .or. any( xnse(iye(:)) <= 0.0_dp ) ) Then

      ii = maxloc( xnse(iye(:)), 1 )
      i = iye(ii)
      c2 = mm52(i) * angm(i)*ggnse(i) * c1 * sqrt(c1) / rhonse
      uvec(1) = (bkt*(log(xnse(i)/c2) - hnse(intz(i))) - be(i)*epmev) / aa(i)
      uvec(2) = uvec(1)

    ! Solve for un and up
    Else
      Do ii = 1, 2
        i = iye(ii)
        c2 = mm52(i) * angm(i)*ggnse(i) * c1 * sqrt(c1) / rhonse
        lhs(ii,1) = nn(i)
        lhs(ii,2) = zz(i)
        rhs(ii) = bkt*(log(xnse(i)/c2) - hnse(intz(i))) - be(i)*epmev
      EndDo
      det = lhs(1,1)*lhs(2,2) - lhs(1,2)*lhs(2,1)
      uvec(1) = rhs(1)*lhs(2,2) - rhs(2)*lhs(1,2)
      uvec(2) = rhs(2)*lhs(1,1) - rhs(1)*lhs(2,1)
      uvec(:) = uvec(:) / det
    EndIf

    If ( itsout >= 2 ) Then
      Write(lun_stdout,'(a,2es23.15)') 'Initial guess               =', uvec(1), uvec(2)
      Write(lun_stdout,'(a5,f11.7,es23.15)') nname(iye(1)), za_ye(1), xnse(iye(1))
      Write(lun_stdout,'(a5,f11.7,es23.15)') nname(iye(2)), za_ye(2), xnse(iye(2))
    EndIf
    If ( itsout >= 3 ) Then
      Write(lun_stdout,'(a,2a5)') 'iza_min,   iza_max   =', nname(iza_min),   nname(iza_max)
      Write(lun_stdout,'(a,2a5)') 'iye_za_l,  iye_za_u  =', nname(iye_za_l),  nname(iye_za_u)
      Write(lun_stdout,'(a,2a5)') 'iye_hvy_l, iye_hvy_u =', nname(iye_hvy_l), nname(iye_hvy_u)
      Write(lun_stdout,'(a,2a5)') 'iye_fe_l,  iye_fe_u  =', nname(iye_fe_l),  nname(iye_fe_u)
    EndIf

    Return
  End Subroutine nse_guess

  Subroutine nse_composition(uvec)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the full composition vector from the n and p chemical potentials.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: uvec(2)

    ! Local variables
    Real(dp) :: c1, temp(ny)
    Real(dp) :: bkt, bktinv, rhoinv
    Integer :: ii, jj

    ! Useful scalars
    bkt = t9nse*bok*epmev
    bktinv = 1.0_dp / bkt
    rhoinv = 1.0_dp / rhonse
    c1 = bkt / (2.0_dp*pi*hbar*hbar*epmev*epmev)
    c1 = c1*sqrt(c1)

    ! Evaluate chemical potentials and mass fractions
    unse(:) = nn(:)*uvec(1) + zz(:)*uvec(2)
    xnse(:) = c1 * rhoinv * angm(1:ny)*ggnse(1:ny) * mm52(:)
    temp(:) = (unse(:) + be(:)*epmev)*bktinv + hnse(intz(:))
    If ( itsout >= 5 ) Then
      ii = maxloc( temp, 1 )
      jj = minloc( temp, 1 )
      Write(lun_stdout,'(a,a5,4es13.5)') 'max exponent: ', &
        & nname(ii), temp(ii), unse(ii)*bktinv, be(ii)*epmev*bktinv, hnse(intz(ii))
      Write(lun_stdout,'(a,a5,4es13.5)') 'min exponent: ', &
        & nname(jj), temp(jj), unse(jj)*bktinv, be(jj)*epmev*bktinv, hnse(intz(jj))
    EndIf
    xnse(:) = xnse(:) * safe_exp( temp(:) )
    xnse(:) = max( 0.0_dp, min( 1.0_dp, xnse(:) ) )
    ynse(:) = xnse(:) / (mm(:)*avn)

    Return
  End Subroutine nse_composition

  Subroutine nse_eval(uvec)
    !-----------------------------------------------------------------------------------------------
    ! This routine updates the RHS with the updated full composition
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: uvec(2)

    ! Local variables
    Real(dp) :: temp(ny)

    ! Update the NSE composition variables, unse, xnse, ynse
    Call nse_composition(uvec)

    ! Mass and charge conservation
    fvec(1) = sum(xnse(:)) - 1.0_dp

    temp(:) = (zz(:)-aa(:)*yense)/(avn*mm(:))
    fvec(2) = sum( temp(:)*xnse(:) )
!   temp(:) = zz(:)/(avn*mm(:))
!   fvec(2) = sum( temp(:)*xnse(:) ) - yense

    knr(3) = knr(3) + 1

    Return
  End Subroutine nse_eval

  Subroutine nse_jac(uvec,feval)
    !-----------------------------------------------------------------------------------------------
    ! This routine updates the Jacobian (LHS) with updated full composition. If feval = .true., then
    ! also update the RHS.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: uvec(2)
    Logical, Intent(in) :: feval

    ! Local variables
    Real(dp) :: dxdn(ny), dxdp(ny)
    Real(dp) :: bkt, bktinv, temp(ny)

    ! Useful scalars
    bkt = t9nse*bok*epmev
    bktinv = 1.0_dp / bkt

    If ( feval ) Then
      ! Update RHS
      Call nse_eval(uvec)
    Else
      ! Update the NSE composition variables, unse, xnse, ynse
      Call nse_composition(uvec)
    EndIf

    temp(:) = (zz(:) - aa(:)*yense) / (avn*mm(:))
!   temp(:) = zz(:)/(avn*mm(:))

    ! Jacobian
    dxdn(:)   = bktinv * nn(:) * xnse(:)
    fjac(1,1) = sum( dxdn(:) )
    fjac(2,1) = sum( temp(:)*dxdn(:) )

    dxdp(:)   = bktinv * zz(:) * xnse(:)
    fjac(1,2) = sum( dxdp(:) )
    fjac(2,2) = sum( temp(:)*dxdp(:) )

    Return
  End Subroutine nse_jac

  Subroutine nse_nr(uvec,check,info)
    !-----------------------------------------------------------------------------------------------
    ! This routine uses the Newton-Raphson method up to nritmax iterations to try to find a solution
    ! for a supplied initial guess.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input/Output variables
    Real(dp), Intent(inout) :: uvec(2)

    ! Output variables
    Logical, Intent(out) :: check
    Integer, Intent(out) :: info

    ! Local variables
    Real(dp) :: uvec0(2), gvec(2), pvec(2)
    Real(dp) :: det, jcond, fjaci(2,2)
    Real(dp) :: f, f0, maxstep, unorm
    Real(dp) :: testf, testu, testg

    Integer :: nrit, ipiv(2)
    Integer :: neval
    Integer :: i, n, info2

    ! Initialize return values
    info = 0
    info2 = -1
    check = .false.
    knr(:) = 0

    n = size(uvec)
    testu = 0.0_dp
    testg = 0.0_dp

    ! Check for initial guess being correct, but with stricter tolerance
    f = nse_func(uvec)
    knrtot(:) = knrtot(:) + knr(:)
    testf = maxval( abs(scalefvec(:)*fvec(:)), 1 )
    If ( testf < 0.01_dp*tolf ) Then
      check = .false.
      info = 1
      Return
    EndIf

    ! Calculate maximum relative step length
    unorm = l2norm( scaleu(:)*uvec(:) )
    maxstep = maxstep0 * max( unorm, l2norm( scaleu(:) ) )

    If ( itsout >= 3 ) Write(lun_stdout,'(3i3,7es23.15)') &
      & 0, knr(2), knr(3), uvec(1), uvec(2), fvec(1), fvec(2), f, testf, testu
    If ( itsout >= 4 ) Write(lun_stdout,'(9x,2es23.15)') xnse(i_nn), xnse(i_pp)

    ! Do the NR iterations
    Do nrit = 1, nritmax

      ! Reset per iteration counters
      knr(:) = 0

      ! Build jacobian (build in nse_func instead)
!     Call nse_jac(uvec,.true.)

      ! Get estimate of condition number
      det = fjac(1,1)*fjac(2,2) - fjac(1,2)*fjac(2,1)
      If ( det > 2.0_dp*epsilon(det) ) Then
        fjaci(1,1) = fjac(2,2)
        fjaci(2,1) = -fjac(2,1)
        fjaci(1,2) = -fjac(1,2)
        fjaci(2,2) = fjac(1,1)
        fjaci(:,:) = fjaci(:,:) / det
        jcond = maxval( sum( abs(fjac(:,:)), 1 ) ) * maxval( sum( abs(fjaci(:,:)), 1 ) )
      Else
        jcond = 0.0_dp
      EndIf

      ! Check for singular matrix
      If ( det <= 2.0_dp*epsilon(det) .or. jcond*tolmin > 1.0_dp ) Then
        If ( itsout >=1 ) Write(lun_stdout,'(a,2es23.15)') 'Singular Matrix: det,jcond = ',det,jcond
        info = -7
        Exit
      EndIf

      ! Calculate search direction, pvec
      pvec(1) = -fvec(1)*fjaci(1,1) - fvec(2)*fjaci(1,2)
      pvec(2) = -fvec(1)*fjaci(2,1) - fvec(2)*fjaci(2,2)

      ! Calculate gradient of fvec, gvec
      Call dgemv('T',n,n,1.0_dp,fjac,n,scalefvec(:)*scalefvec(:)*fvec(:),1,0.0_dp,gvec,1)

!     pvec(:) = -fvec(:)
!     Call dgesv(n,1,fjac,n,ipiv,pvec,n,info)

      ! Save old values
      uvec0(:) = uvec(:)
      f0 = f

      ! Backtracking line search to get step length and update uvec, f, fvec
      Call nse_lnsrch(uvec0,f0,uvec,f,gvec,pvec,maxstep,check,info2)

      knrtot(:) = knrtot(:) + knr(:)

      ! Test for convergence
      If ( lnorm == 0 ) Then
        testf = maxval( abs(scalefvec(:)*fvec(:)), 1 )
        testu = maxval( abs(uvec(:)-uvec0(:)) / merge( abs(uvec(:)), typu(:), abs(uvec(:))>typu(:) ), 1 )
        testg = maxval( abs(gvec(:)) * merge( abs(uvec(:)), typu(:), abs(uvec(:))>typu(:) ) / max(f, 0.5_dp*n), 1 )
      ElseIf ( lnorm == 1 ) Then
        testf = sum( abs(scalefvec(:)*fvec(:)) )
        testu = sum( abs(uvec(:)-uvec0(:)) / merge( abs(uvec(:)), typu(:), abs(uvec(:))>typu(:) ) )
        testg = sum( abs(gvec(:)) * merge( abs(uvec(:)), typu(:), abs(uvec(:))>typu(:) ) / max(f, 0.5_dp*n) )
      ElseIf ( lnorm == 2 ) Then
        testf = l2norm( abs(scalefvec(:)*fvec(:)) )
        testu = l2norm( abs(uvec(:)-uvec0(:)) / merge( abs(uvec(:)), typu(:), abs(uvec(:))>typu(:) ) )
        testg = l2norm( abs(gvec(:)) * merge( abs(uvec(:)), typu(:), abs(uvec(:))>typu(:) ) / max(f, 0.5_dp*n) )
      EndIf

      If ( itsout >= 3 ) Write(lun_stdout,'(3i3,8es23.15)') &
        & nrit, knr(2), knr(3), uvec(1), uvec(2), fvec(1), fvec(2), f, testf, testu, testg
      If ( itsout >= 4 ) Write(lun_stdout,'(7x,i2,6es23.15)') &
        & info2, xnse(1), xnse(2), pvec(1), pvec(2), gvec(1), gvec(2)

      If ( info2 < 0 ) Then
        info = info2
        Exit
      ElseIf ( testf < tolf ) Then
        info = 1
        Exit
      ElseIf ( check ) Then
        If ( testg < tolmin ) Then
          info = -6
          Exit
        Else
          info = -3
          Exit
        EndIf
      ElseIf ( testu < tolu ) Then
        info = 2
        Exit
      Else
        info = 0
      EndIf
    EndDo

    If ( info <= 0 ) Then
      check = .true.
    EndIf
    knr(1) = nrit
    knrtot(1) = knrtot(1) + min( nrit, nritmax )

    Return
  End Subroutine nse_nr

  Subroutine nse_lnsrch(uvec0,f0,uvec,f,gvec,pvec,maxstep,check,info)
    !-----------------------------------------------------------------------------------------------
    ! This routine performs a backtracking line search to get step length and update uvec, f, fvec
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: uvec0(2), f0, maxstep

    ! Input/Output variables
    Real(dp), Intent(inout) :: pvec(2), gvec(2)

    ! Output variables
    Real(dp), Intent(out) :: uvec(2), f
    Logical, Intent(out) :: check
    Integer, Intent(out) :: info

    ! Local variables
    Integer :: i, lsit, ipiv(4), ierr, j

    Real(dp) :: gveck(2)
    Real(dp) :: betak(0:lsitmax+1), fk(0:lsitmax+1), slopek(0:lsitmax+1)
    Real(dp) :: newtlen, slope, plength, beta_min
    Real(dp) :: rbsq, rb0sq, rbeta, rbeta0, rdbeta
    Real(dp) :: rhs(4), c1, c2, c3, disc, vmat(4,4)

    check = .false.
    info = -1

    ! Scale step length if necessary
    newtlen = l2norm( scaleu(:)*pvec(:) )
    If ( newtlen > maxstep ) Then
      pvec(:) = pvec(:) * maxstep/newtlen
    EndIf

    ! Calculate "slope" of step direction
    slope = dot_prod( gvec(:), pvec(:) )
!   slope = -dot_prod( fvec(:), fvec(:) )
    If ( slope >= 0.0_dp ) Then
      info = -4
      Return
    EndIf

    ! Determine minimum allowable step length
    plength = maxval( abs(pvec(:)) / merge( abs(uvec0(:)), typu(:), abs(uvec0(:))>typu(:) ), 1 )
    beta_min = tolbeta / plength

    ! Start with full step length
    fk(:) = f0
    betak(0) = 0.0_dp
    betak(1) = 1.0_dp
    slopek(:) = slope
    If ( itsout >= 4 ) Write(lun_stdout,'(3x,i3,3x,4es23.15)') 0, beta_min, f0, f0+alpha*betak(1)*slope, slope

    Do lsit = 1, lsitmax

      ! Update uvec, f, fvec
      uvec(:) = uvec0(:) + betak(lsit)*pvec(:)
      fk(lsit) = nse_func(uvec)

      If ( itsout >= 4 ) Write(lun_stdout,'(3x,i3,3x,5es23.15)') &
        & lsit, betak(lsit), fk(lsit), fk(lsit-1) + alpha*betak(lsit)*slopek(lsit-1), slopek(lsit-1), &
        & dot_prod( (gveck(:)-gvec(:)), pvec )

      ! Test if sufficient decrease condition has been met
      If ( fk(lsit) <= f0 + alpha*betak(lsit)*slope ) Then
        check = .false.
        info = 0
        f = fk(lsit)
        Exit
      ElseIf ( betak(lsit) < beta_min ) Then
        check = .true.
        info = 3
        uvec(:) = uvec0(:)
        f = f0
        Exit
      ElseIf ( lsit > 200 ) Then
        Call dgemv('T',2,2,1.0_dp,fjac,2,scalefvec(:)*scalefvec(:)*fvec(:),1,0.0_dp,gveck,1)
        slopek(lsit) = dot_prod( gveck(:), pvec(:) )
        rhs(:) = fk(lsit-3:lsit)
        vmat(:,1) = betak(lsit-3:lsit)**3
        vmat(:,2) = betak(lsit-3:lsit)**2
        vmat(:,3) = betak(lsit-3:lsit)
        vmat(:,4) = 1.0_dp
      ElseIf ( lsit > 100 ) Then
        Call dgemv('T',2,2,1.0_dp,fjac,2,scalefvec(:)*scalefvec(:)*fvec(:),1,0.0_dp,gveck,1)
        slopek(lsit) = dot_prod( gveck(:), pvec(:) )
        rhs(1:2) = fk(lsit-1:lsit)
        rhs(3:4) = slopek(lsit-1:lsit)
        vmat(1:2,1) = betak(lsit-1:lsit)**3
        vmat(3:4,1) = 3.0_dp*betak(lsit-1:lsit)**2
        vmat(1:2,2) = betak(lsit-1:lsit)**2
        vmat(3:4,2) = 2.0_dp*betak(lsit-1:lsit)
        vmat(1:2,3) = betak(lsit-1:lsit)
        vmat(3:4,3) = 1.0_dp
        vmat(1:2,4) = 1.0_dp
        vmat(3:4,4) = 0.0_dp
      ElseIf ( lsit >= 1 ) Then
!       ! Choose next step length that minimizes the quadratic interpolating polynomial:
!       !     P = c1*beta(k+1)^2 + c2*beta(k+1) + c3
!       c2 = (fk(lsit) - f0 - slope*betak(lsit))*rbsq
!       c3 = slope
!       betak(lsit+1) = -0.5 * c2 / c1

        ! We know fk(0), slope(0), fk(k), and slope(k), so...
        ! Find betak(k+1) that minimizes the cubic interpolating polynomial:
        !     P = c1*beta(k+1)^3 + c2*beta(k+1)^2 + c3*beta(k+1) + c4
        !     Get the coefficients from cn = V^-1 * y_n, for Vandermonde matrix V
        Call dgemv('T',2,2,1.0_dp,fjac,2,scalefvec(:)*scalefvec(:)*fvec(:),1,0.0_dp,gveck,1)
        slopek(lsit) = dot_prod( gveck(:), pvec(:) )
        rhs(1:2) = fk(lsit-1:lsit)
        rhs(3:4) = slopek(lsit-1:lsit)
        vmat(1:2,1) = betak(lsit-1:lsit)**3
        vmat(3:4,1) = 3.0_dp*betak(lsit-1:lsit)**2
        vmat(1:2,2) = betak(lsit-1:lsit)**2
        vmat(3:4,2) = 2.0_dp*betak(lsit-1:lsit)
        vmat(1:2,3) = betak(lsit-1:lsit)
        vmat(3:4,3) = 1.0_dp
        vmat(1:2,4) = 1.0_dp
        vmat(3:4,4) = 0.0_dp
      Else
        Call dgemv('T',2,2,1.0_dp,fjac,2,scalefvec(:)*scalefvec(:)*fvec(:),1,0.0_dp,gveck,1)
        slopek(lsit) = dot_prod( gveck(:), pvec(:) )

        ! We know fk(0), slope(0), fk(k-1), and fk(k), so...
        rhs(1:3) = fk(lsit-2:lsit)
        rhs(4) = slopek(lsit-2)
        vmat(1:3,1) = betak(lsit-2:lsit)**3
        vmat(4,1) = 3.0_dp*betak(lsit-2)**2
        vmat(1:3,2) = betak(lsit-2:lsit)**2
        vmat(4,2) = 2.0_dp*betak(lsit-2)
        vmat(1:3,3) = betak(lsit-2:lsit)
        vmat(4,3) = 1.0_dp
        vmat(1:3,4) = 1.0_dp
        vmat(4,4) = 0.0_dp
      EndIf

      Call dgesv(4,1,vmat,4,ipiv,rhs,4,ierr)
      c1 = rhs(1)
      c2 = rhs(2)
      c3 = rhs(3)
      If ( abs(c1) < tiny(0.0_dp) ) Then
        ! Minimum of quadratic interpolating polynomial
        betak(lsit+1) = -0.5_dp*c3/c2
      Else
        ! Choose the root which corresponds to the local minimum
        disc = c2*c2 - 3.0_dp*c1*c3
        If ( disc < 0.0_dp ) Then
          ! No roots exist, so set next step length to its upper-bound
          betak(lsit+1) = gamma2*betak(lsit)
        ElseIf ( c2 < 0.0_dp ) Then
          betak(lsit+1) = (-c2 + sqrt(disc)) / (3.0_dp*c1)
        Else
          betak(lsit+1) = -c3 / (c2 + sqrt(disc))
        EndIf
      EndIf

      ! Limit change in step length
      betak(lsit+1) = min( betak(lsit+1), gamma2*betak(lsit) )
      betak(lsit+1) = max( betak(lsit+1), gamma1*betak(lsit) )
    EndDo

    If ( lsit > lsitmax ) Then
      check = .true.
      f = fk(lsitmax)
    EndIf
    knr(2) = lsit

    Return
  End Subroutine nse_lnsrch

  Function nse_func(uvec)
    !-----------------------------------------------------------------------------------------------
    ! This is the actual function we are trying to minimize ( 0.5 * f dot f )
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: uvec(2)

    ! Function variable
    Real(dp) :: nse_func

    Call nse_jac(uvec,.true.)

    nse_func = 0.5_dp * dot_prod( scalefvec(:)*fvec(:), scalefvec(:)*fvec(:) )

    Return
  End Function nse_func

  Subroutine nse_pf
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates the nuclear partition functions as a function of t9nse.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Local variables
    Real(dp) :: rdt9
    Integer :: i, ii

    Do i = 1, ng
      If ( t9i(i) >= t9nse ) Exit
    EndDo
    ii = i

    Select Case (ii)
    Case (1)
      ggnse(1:ny) = g(1,1:ny)
    Case (ng+1)
      ggnse(1:ny) = g(ng,1:ny)
    Case Default
      rdt9 = (t9nse-t9i(ii-1)) / (t9i(ii)-t9i(ii-1))
      ggnse(1:ny) = safe_exp( rdt9*log(g(ii,1:ny)) + (1-rdt9)*log(g(ii-1,1:ny)) )
    End Select

    Return
  End Subroutine nse_pf

  Subroutine nse_screen
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates screening factors (i.e. Coulomb corrections) for NSE in a manner that
    ! is consistent with XNet reaction network screening calculation.
    !-----------------------------------------------------------------------------------------------
    Use xnet_eos, Only: eos_interface, eos_screen
    Implicit None

    ! Local variables
    Integer, Parameter :: iz1 = 1
    Real(dp), Parameter :: z1 = 1.0_dp

    Real(dp), Parameter :: lam_1 = 0.1_dp
    Real(dp), Parameter :: lam_2 = 0.125_dp
    Real(dp), Parameter :: lam_3 = 2.0_dp
    Real(dp), Parameter :: lam_4 = 2.15_dp
    Real(dp), Parameter :: lam_5 = 4.85_dp
    Real(dp), Parameter :: lam_6 = 5.0_dp

    Integer :: i, j, iz2(izmax-1), izc(izmax-1)
    Real(dp) :: z2(izmax-1), zetaw(izmax-1), zetai(izmax-1)
    Real(dp) :: h0(izmax-1), hw(izmax-1), hi(izmax-1), hs(izmax-1), lambda12(izmax-1)
    Real(dp) :: theta12iw(izmax-1), theta12is(izmax-1), theta12si(izmax-1)
    Real(dp) :: fhs(0:izmax+1), fhi(0:izmax+1), gammaz(0:izmax+1)
    Real(dp) :: cv, etae, detaedt9, ztot, ztilde, zinter, lambda0, gammae, dztildedt9, s0

    If ( iscrn <= 0 ) Then
      hnse = 0.0_dp
    ElseIf ( .not. use_CP98 ) Then

      ! Call EOS to get plasma quantities
      call eos_interface(t9nse,rhonse,ynse,ztot,cv,etae,detaedt9)
      Call eos_screen(t9nse,rhonse,ynse,etae,detaedt9,ztilde,zinter,lambda0,gammae,dztildedt9)

      !---------------------------------------------------------------------------------------------
      ! Calculate screening energies as a function of Z, for prescriptions that follow this approach
      !---------------------------------------------------------------------------------------------
      gammaz(0) = 0.0_dp
      gammaz(1:izmax+1) = gammae * zseq53(1:izmax+1)
      fhi(0) = 0.0_dp
      fhi(1:izmax+1) = kbi * zinter * lambda0**bi * zseqi(1:izmax+1)
      fhs(0) = 0.0_dp
      fhs(1:izmax+1) = + cds(1) * gammaz(1:izmax+1) &
        &              + cds(2) * gammaz(1:izmax+1)**cds(5) / cds(5) &
        &              + cds(3) * log(gammaz(1:izmax+1)) &
        &              + cds(4)

      z2 = zseq(1:izmax-1)
      iz2 = nint(z2)
      izc = iz1 + iz2

      ! Weak and intermediate screening factors, Table 4 of Graboske et al. (1973)
      zetaw = z1 * z2
      lambda12 = zetaw * ztilde * lambda0
      hw = lambda12
      hi = fhi(izc) - fhi(iz1) - fhi(iz2)
!     zetai = (z1 + z2)**bip1 - z1**bip1 - z2**bip1
!     hi = kbi * zinter * lambda0**bi * zetai

      ! Strong screening from Dewitt & Slattery (2003) using linear mixing.
      hs = fhs(iz1) + fhs(iz2) - fhs(izc)

      ! Select screening factor
      theta12iw = max( 0.0, min( 1.0, (lambda12 - lam_1) / (lam_2 - lam_1) ) )
      theta12is = max( 0.0, min( 1.0, (lambda12 - lam_5) / (lam_6 - lam_5) ) )
      theta12si = max( 0.0, min( 1.0, (lambda12 - lam_3) / (lam_4 - lam_3) ) )
      Where ( iz2 == 0 )
        h0 = 0.0_dp
      ElseWhere ( lambda12 < lam_1 )
        h0 = hw
      ElseWhere ( lambda12 < lam_3 )
        h0 = theta12iw*hi + (1.0_dp - theta12iw)*hw
      ElseWhere ( lambda12 > lam_6 )
        h0 = hs
      ElseWhere ( hi < hs )
        h0 = theta12is*hi + (1.0_dp - theta12is)*hs
      ElseWhere
        h0 = theta12si*hs + (1.0_dp - theta12si)*hi
      EndWhere

      ! Add succeeding screening factors
      hnse(1) = 0.0_dp
      s0 = 0.0_dp
      Do j = 2, izmax
        s0 = s0 + h0(j-1)
        hnse(j) = s0
      EndDo

      If ( idiag >= 5 ) Write(lun_diag,'(a5,3i6,6es23.15)') &
        & ('HNSE',iz1,iz2(j-1),j,lambda12(j-1),h0(j-1),hw(j-1),hi(j-1),hs(j-1),hnse(j),j=2,izmax)

    EndIf

    Return
  End Subroutine nse_screen

  Subroutine nse_screen_CP98
    !-----------------------------------------------------------------------------------------------
    ! This routine calculates a simpler Coulomb correction term as a backup if XNet screening fails.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Local variables
    Integer, Parameter :: iz1 = 1
    Real(dp), Parameter :: z1 = 1.0_dp
    Real(dp), Parameter :: a1 = -0.9052_dp, a2 = 0.6322_dp
    Real(dp), Parameter :: a2inv = 1.0_dp / a2
    Real(dp), Parameter :: a3 = -sqrt(3.0_dp)*2.0_dp - a1*sqrt(a2inv)

    Integer :: j, iz2(izmax), izc(izmax)
    Real(dp) :: z2(izmax), h0(izmax)
    Real(dp) :: fh0(izmax+1)
    Real(dp) :: ztot, ztilde, zinter, lambda0, gammae, dztildedt9
    Real(dp) :: gz, z, sqrtgz, ae

    If ( iscrn <= 0 ) Then
      hnse = 0.0_dp
    Else
      ae = (3.0_dp / (4.0_dp*pi*avn*rhonse*yense))**third ! electron-sphere radius
      gammae = e2 / (ae*t9nse*bok) ! electron Coulomb coupling parameter

      !---------------------------------------------------------------------------------------------
      ! Calculate screening energies as a function of Z, for prescriptions that follow this approach
      !---------------------------------------------------------------------------------------------
      Do j = 1, izmax+1
        gz = gammae * zseq53(j)
        sqrtgz = sqrt(gz)
        fh0(j) = a1*( sqrt(gz*(a2 + gz)) - a2*log(sqrt(gz*a2inv) + sqrt(1.0_dp+gz*a2inv)) ) &
          &      + 2.0_dp*a3*( sqrtgz - atan(sqrtgz) )
      EndDo

      z2 = zseq(1:izmax)
      iz2 = nint(z2)
      izc = iz1 + iz2

!     h0 = z2*fh0(1) - fh0
      h0 = fh0(iz1) + fh0(iz2) - fh0(iz1+iz2)

      hnse(1) = 0.0_dp
      hnse(2:izmax) = z2(2:izmax)*fh0(1) - fh0(2:izmax)

      If ( idiag >= 5 ) Write(lun_diag,"(a5,3i6,23x,es23.15,69x,es23.15)") &
        & ('HNSE',iz1,iz2(j-1),j,sum(h0(1:(j-1))),hnse(j),j=2,izmax)

    EndIf

    Return
  End Subroutine nse_screen_CP98

End Module xnet_nse
