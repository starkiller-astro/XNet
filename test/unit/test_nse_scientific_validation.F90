Module test_nse_scientific_validation
  Use testdrive, Only: check, error_type, new_unittest, unittest_type
  Use xnet_types, Only: dp
  Use nse_validation_support, Only: nse_fail_charge, nse_fail_dominant, nse_fail_duplicate, &
    & nse_fail_finite, nse_fail_identity, nse_fail_l1, nse_fail_linf, nse_fail_mass, &
    & nse_fail_nonnegative, nse_fail_ye, nse_validation_metrics, nse_validation_tolerances, &
    & apply_nse_validation_gates, evaluate_nse_candidate
  Implicit None
  Private

  Character(32), Allocatable :: state_id(:)
  Character(5), Allocatable :: expected_names(:,:)
  Integer, Allocatable :: expected_a(:,:), expected_n(:,:), expected_z(:,:)
  Real(dp), Allocatable :: expected_x(:,:), state_rho(:), state_t9(:), state_ye(:)
  Type(nse_validation_tolerances), Allocatable :: state_tolerances(:)
  Integer :: state_count = 0

  Public :: collect_nse_scientific_validation, load_validation_fixture

Contains

  Subroutine collect_nse_scientific_validation(testsuite)
    Implicit None

    Type(unittest_type), Allocatable, Intent(out) :: testsuite(:)

    testsuite = [ &
      & new_unittest('three independent unscreened NSE states',test_scientific_states), &
      & new_unittest('dominant and complete composition mutations',test_composition_mutations), &
      & new_unittest('Ye and normalization mutations',test_conservation_mutations), &
      & new_unittest('species identity and invalid-value mutations',test_identity_mutations), &
      & new_unittest('binding-energy mutation',test_binding_energy_mutation), &
      & new_unittest('every numerical tolerance boundary',test_tolerance_boundaries), &
      & new_unittest('scientific-state supplied guess repeatability', &
      & test_scientific_guess_repeatability) ]

    Return
  End Subroutine collect_nse_scientific_validation

  Subroutine load_validation_fixture(network_directory,reference_path)
    Implicit None

    Character(*), Intent(in) :: network_directory, reference_path

    Call load_network(network_directory)
    Call load_reference(reference_path)

    Return
  End Subroutine load_validation_fixture

  Subroutine load_network(network_directory)
    Use nuclear_data, Only: aa, angm, be, g, ia, iz, izmax, mex, mm, ng, nname, nn, ny, t9i, &
      & zz, zz2, zzi, zseq, zseq53, zseqi
    Use xnet_constants, Only: avn, bip1, five3rd, thbim1
    Use xnet_controls, Only: idiag, iscrn, itsout
    Use xnet_nse, Only: nse_initialize
    Implicit None

    Character(*), Intent(in) :: network_directory

    Character(512) :: filename
    Character(5) :: record_name
    Integer :: ierr, inuc, it9i(ng), j, lun, neutron_index, proton_index
    Integer, Allocatable :: neutron_number(:)
    Real(dp) :: spin

    filename = trim(network_directory)//'/netwinv'
    Open(newunit=lun,file=trim(filename),status='old',action='read',iostat=ierr)
    If ( ierr /= 0 ) Call fail_fixture('cannot open NSE validation netwinv fixture')
    Read(lun,'(i5)',iostat=ierr) ny
    If ( ierr /= 0 .OR. ny /= 489 ) Call fail_fixture('invalid NSE validation species count')
    Read(lun,'(24i3)',iostat=ierr) it9i
    If ( ierr /= 0 ) Call fail_fixture('invalid NSE validation temperature grid')

    Allocate(aa(ny),angm(ny),be(ny),g(ng,ny),ia(ny),iz(ny),mex(ny),mm(ny),nname(ny), &
      & nn(ny),t9i(ng),zz(ny),zz2(ny),zzi(ny),neutron_number(ny))
    Do inuc = 1, ny
      Read(lun,'(a5)',iostat=ierr) nname(inuc)
      If ( ierr /= 0 ) Call fail_fixture('invalid NSE validation species list')
    EndDo
    Do inuc = 1, ny
      Read(lun,*,iostat=ierr) record_name, aa(inuc), iz(inuc), neutron_number(inuc), &
        & spin, mex(inuc)
      If ( ierr /= 0 .OR. adjustl(record_name) /= adjustl(nname(inuc)) ) &
        & Call fail_fixture('invalid NSE validation nuclear record')
      Read(lun,*,iostat=ierr) (g(j,inuc),j=1,ng)
      If ( ierr /= 0 ) Call fail_fixture('invalid NSE validation partition factors')
      angm(inuc) = 2.0_dp*spin + 1.0_dp
    EndDo
    Close(lun)

    ia = nint(aa)
    zz = real(iz,dp)
    nn = real(neutron_number,dp)
    neutron_index = find_species_index(0,1)
    proton_index = find_species_index(1,1)
    If ( neutron_index == 0 .OR. proton_index == 0 ) &
      & Call fail_fixture('NSE validation network requires free n and p')
    be = nn*mex(neutron_index) + zz*mex(proton_index) - mex
    mm = aa/avn
    t9i = real(it9i,dp)
    t9i(1:ng-1) = 0.01_dp*t9i(1:ng-1)
    t9i(ng) = 0.1_dp*t9i(ng)
    zz2 = zz*zz
    zzi = zz**thbim1
    izmax = maxval(iz)
    Allocate(zseq(0:izmax+2),zseq53(0:izmax+2),zseqi(0:izmax+2))
    zseq = (/ (real(j,dp),j=0,izmax+2) /)
    zseq53 = zseq**five3rd
    zseqi = zseq**bip1
    Deallocate(neutron_number)

    idiag = 0
    iscrn = 0
    itsout = 0
    Call nse_initialize

    Return
  End Subroutine load_network

  Subroutine load_reference(reference_path)
    Use nuclear_data, Only: ny
    Implicit None

    Character(*), Intent(in) :: reference_path

    Character(512) :: line
    Character(64) :: order_hash
    Integer :: ierr, inuc, lun, reference_species_count, state

    Open(newunit=lun,file=trim(reference_path),status='old',action='read',iostat=ierr)
    If ( ierr /= 0 ) Call fail_fixture('cannot open NSE validation reference data')
    Read(lun,'(a)',iostat=ierr) line
    If ( ierr /= 0 .OR. trim(line) /= 'XNET_NSE_REFERENCE_V2' ) &
      & Call fail_fixture('invalid NSE validation reference schema')
    Read(lun,*,iostat=ierr) reference_species_count, state_count
    If ( ierr /= 0 .OR. reference_species_count /= ny .OR. state_count /= 3 ) &
      & Call fail_fixture('invalid NSE validation reference dimensions')
    Read(lun,'(a64)',iostat=ierr) order_hash
    If ( ierr /= 0 .OR. len_trim(order_hash) /= 64 ) &
      & Call fail_fixture('invalid NSE validation order hash')

    Allocate(state_id(state_count),state_rho(state_count),state_t9(state_count), &
      & state_ye(state_count),state_tolerances(state_count))
    Allocate(expected_names(ny,state_count),expected_a(ny,state_count), &
      & expected_z(ny,state_count),expected_n(ny,state_count),expected_x(ny,state_count))
    Do state = 1, state_count
      Read(lun,'(a)',iostat=ierr) line
      If ( ierr /= 0 .OR. line(1:6) /= 'STATE ' ) &
        & Call fail_fixture('invalid NSE validation state header')
      state_id(state) = adjustl(line(7:))
      Read(lun,*,iostat=ierr) state_rho(state),state_t9(state),state_ye(state)
      If ( ierr /= 0 ) Call fail_fixture('invalid NSE validation state inputs')
      Read(lun,*,iostat=ierr) state_tolerances(state)%mass, &
        & state_tolerances(state)%charge,state_tolerances(state)%ye, &
        & state_tolerances(state)%l1,state_tolerances(state)%linf
      If ( ierr /= 0 ) Call fail_fixture('invalid NSE validation state tolerances')
      Do inuc = 1, ny
        Read(lun,'(a)',iostat=ierr) line
        If ( ierr /= 0 ) Call fail_fixture('invalid NSE validation composition row')
        expected_names(inuc,state) = line(1:5)
        Read(line(7:),*,iostat=ierr) expected_a(inuc,state),expected_z(inuc,state), &
          & expected_n(inuc,state),expected_x(inuc,state)
        If ( ierr /= 0 ) Call fail_fixture('invalid NSE validation mass fraction')
      EndDo
    EndDo
    Read(lun,'(a)',iostat=ierr) line
    If ( ierr == 0 ) Call fail_fixture('unexpected trailing NSE validation reference data')
    Close(lun)

    Return
  End Subroutine load_reference

  Integer Function find_species_index(proton_number,mass_number) Result(index)
    Use nuclear_data, Only: aa, iz, ny
    Implicit None

    Integer, Intent(in) :: mass_number, proton_number

    Integer :: inuc

    index = 0
    Do inuc = 1, ny
      If ( iz(inuc) == proton_number .AND. nint(aa(inuc)) == mass_number ) Then
        index = inuc
        Exit
      EndIf
    EndDo

    Return
  End Function find_species_index

  Subroutine solve_state(state)
    Use xnet_controls, Only: iscrn
    Use xnet_nse, Only: nse_solve
    Implicit None

    Integer, Intent(in) :: state

    iscrn = 0
    Call nse_solve(state_rho(state),state_t9(state),state_ye(state))

    Return
  End Subroutine solve_state

  Subroutine evaluate_state(state,candidate,candidate_names,metrics,failures)
    Use nuclear_data, Only: aa, nname, nn, zz
    Implicit None

    Integer, Intent(in) :: state
    Character(5), Intent(in) :: candidate_names(:)
    Real(dp), Intent(in) :: candidate(:)
    Type(nse_validation_metrics), Intent(out) :: metrics
    Integer, Intent(out) :: failures

    Call evaluate_nse_candidate(expected_names(:,state),candidate_names, &
      & expected_a(:,state),expected_z(:,state),expected_n(:,state),aa,zz,nn, &
      & expected_x(:,state),candidate,state_ye(state),state_tolerances(state),metrics,failures)

    Return
  End Subroutine evaluate_state

  Subroutine test_scientific_states(error)
    Use, Intrinsic :: iso_fortran_env, Only: error_unit
    Use nuclear_data, Only: nname
    Use xnet_nse, Only: knrtot, xnse
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error

    Integer :: failures, state
    Type(nse_validation_metrics) :: metrics

    Do state = 1, state_count
      Call solve_state(state)
      Call evaluate_state(state,xnse,nname,metrics,failures)
      Write(error_unit,'(2a,5(a,es12.4),a,3(i0,1x))') '# NSE metrics ', &
        & trim(state_id(state)),' mass=',metrics%mass_residual, &
        & ' charge=',metrics%charge_residual,' Ye=',metrics%ye_error, &
        & ' L1=',metrics%l1,' Linf=',metrics%linf,' counters=',knrtot
      If ( failures /= 0 ) Then
        Write(error_unit,'(2a,i0,5(a,es12.4))') 'NSE validation state ',trim(state_id(state)), &
          & failures,' mass=',metrics%mass_residual,' charge=',metrics%charge_residual, &
          & ' Ye=',metrics%ye_error,' L1=',metrics%l1,' Linf=',metrics%linf
      EndIf
      Call check(error,failures,0)
      If ( allocated(error) ) Return
      Call check(error,knrtot(3) > 0)
      If ( allocated(error) ) Return
    EndDo

    Return
  End Subroutine test_scientific_states

  Subroutine test_composition_mutations(error)
    Use nuclear_data, Only: nname, ny
    Use xnet_nse, Only: xnse
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error

    Integer :: dominant, failures
    Real(dp), Allocatable :: mutated_expected(:)
    Type(nse_validation_metrics) :: metrics

    Call solve_state(1)
    Allocate(mutated_expected(ny))
    dominant = maxloc(expected_x(:,1),dim=1)
    mutated_expected = expected_x(:,1)
    mutated_expected(dominant) = mutated_expected(dominant) &
      & + 2.0_dp*state_tolerances(1)%linf
    Call evaluate_nse_candidate(expected_names(:,1),nname, &
      & expected_a(:,1),expected_z(:,1),expected_n(:,1),get_mass_numbers(), &
      & get_proton_numbers(),get_neutron_numbers(),mutated_expected,xnse, &
      & state_ye(1),state_tolerances(1),metrics,failures)
    Call require_flag(error,failures,nse_fail_dominant)
    If ( allocated(error) ) Return
    Call require_flag(error,failures,nse_fail_linf)
    If ( allocated(error) ) Return

    mutated_expected = 0.99_dp*expected_x(:,1)
    Call evaluate_nse_candidate(expected_names(:,1),nname, &
      & expected_a(:,1),expected_z(:,1),expected_n(:,1),get_mass_numbers(), &
      & get_proton_numbers(),get_neutron_numbers(),mutated_expected,xnse, &
      & state_ye(1),state_tolerances(1),metrics,failures)
    Call require_flag(error,failures,nse_fail_l1)
    Deallocate(mutated_expected)

    Return
  End Subroutine test_composition_mutations

  Subroutine test_conservation_mutations(error)
    Use nuclear_data, Only: nname, ny
    Use xnet_nse, Only: i_nn, i_pp, xnse
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error

    Integer :: failures
    Real(dp) :: delta
    Real(dp), Allocatable :: candidate(:)
    Type(nse_validation_metrics) :: metrics

    Call solve_state(1)
    Allocate(candidate(ny))
    candidate = xnse*(1.0_dp + 2.0_dp*state_tolerances(1)%mass)
    Call evaluate_state(1,candidate,nname,metrics,failures)
    Call require_flag(error,failures,nse_fail_mass)
    If ( allocated(error) ) Return

    candidate = xnse
    delta = 2.0_dp*max(state_tolerances(1)%charge,state_tolerances(1)%ye)
    candidate(i_nn) = candidate(i_nn) + delta
    candidate(i_pp) = candidate(i_pp) - delta
    Call evaluate_state(1,candidate,nname,metrics,failures)
    Call require_flag(error,failures,nse_fail_charge)
    If ( allocated(error) ) Return
    Call require_flag(error,failures,nse_fail_ye)
    Deallocate(candidate)

    Return
  End Subroutine test_conservation_mutations

  Subroutine test_identity_mutations(error)
    Use, Intrinsic :: ieee_arithmetic, Only: ieee_quiet_nan, ieee_value
    Use nuclear_data, Only: aa, nname, nn, ny, zz
    Use xnet_nse, Only: xnse
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error

    Character(5) :: saved_name
    Character(5), Allocatable :: candidate_names(:)
    Integer :: failures
    Real(dp), Allocatable :: candidate(:), candidate_a(:)
    Type(nse_validation_metrics) :: metrics

    Call solve_state(1)
    Allocate(candidate(ny),candidate_a(ny),candidate_names(ny))
    candidate = xnse
    candidate_names = nname
    saved_name = candidate_names(1)
    candidate_names(1) = candidate_names(2)
    candidate_names(2) = saved_name
    Call evaluate_state(1,candidate,candidate_names,metrics,failures)
    Call require_flag(error,failures,nse_fail_identity)
    If ( allocated(error) ) Return

    candidate_names = nname
    candidate_names(2) = candidate_names(1)
    Call evaluate_state(1,candidate,candidate_names,metrics,failures)
    Call require_flag(error,failures,nse_fail_duplicate)
    If ( allocated(error) ) Return

    Call evaluate_nse_candidate(expected_names(:,1),nname(:ny-1),expected_a(:,1), &
      & expected_z(:,1),expected_n(:,1),aa(:ny-1),zz(:ny-1),nn(:ny-1), &
      & expected_x(:,1),candidate(:ny-1),state_ye(1),state_tolerances(1),metrics,failures)
    Call require_flag(error,failures,nse_fail_identity)
    If ( allocated(error) ) Return

    candidate_a = aa
    candidate_a(ny) = candidate_a(ny) + 1.0_dp
    Call evaluate_nse_candidate(expected_names(:,1),nname,expected_a(:,1), &
      & expected_z(:,1),expected_n(:,1),candidate_a,zz,nn,expected_x(:,1), &
      & candidate,state_ye(1),state_tolerances(1),metrics,failures)
    Call require_flag(error,failures,nse_fail_identity)
    If ( allocated(error) ) Return

    candidate = xnse
    candidate(1) = ieee_value(candidate(1),ieee_quiet_nan)
    Call evaluate_state(1,candidate,nname,metrics,failures)
    Call require_flag(error,failures,nse_fail_finite)
    If ( allocated(error) ) Return

    candidate = xnse
    candidate(1) = -tiny(1.0_dp)
    Call evaluate_state(1,candidate,nname,metrics,failures)
    Call require_flag(error,failures,nse_fail_nonnegative)
    Deallocate(candidate,candidate_a,candidate_names)

    Return
  End Subroutine test_identity_mutations

  Subroutine test_binding_energy_mutation(error)
    Use nuclear_data, Only: be, nname
    Use xnet_nse, Only: xnse
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error

    Integer :: failures, target
    Real(dp) :: saved_binding
    Type(nse_validation_metrics) :: metrics

    target = find_name_index(' co55')
    If ( target == 0 ) Call fail_fixture('binding mutation species is missing')
    saved_binding = be(target)
    be(target) = saved_binding + 0.010_dp
    Call solve_state(1)
    Call evaluate_state(1,xnse,nname,metrics,failures)
    be(target) = saved_binding
    Call require_flag(error,failures,nse_fail_l1)
    If ( allocated(error) ) Return
    Call require_flag(error,failures,nse_fail_linf)

    Return
  End Subroutine test_binding_energy_mutation

  Subroutine test_tolerance_boundaries(error)
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error

    Integer :: failures, state
    Type(nse_validation_metrics) :: metrics
    Type(nse_validation_tolerances) :: tolerances

    Do state = 1, state_count
      tolerances = state_tolerances(state)
      metrics = nse_validation_metrics()
      failures = 0
      metrics%mass_residual = nearest(tolerances%mass,1.0_dp)
      Call apply_nse_validation_gates(metrics,tolerances,failures)
      Call require_flag(error,failures,nse_fail_mass)
      If ( allocated(error) ) Return

      metrics = nse_validation_metrics()
      failures = 0
      metrics%charge_residual = nearest(tolerances%charge,1.0_dp)
      Call apply_nse_validation_gates(metrics,tolerances,failures)
      Call require_flag(error,failures,nse_fail_charge)
      If ( allocated(error) ) Return

      metrics = nse_validation_metrics()
      failures = 0
      metrics%ye_error = nearest(tolerances%ye,1.0_dp)
      Call apply_nse_validation_gates(metrics,tolerances,failures)
      Call require_flag(error,failures,nse_fail_ye)
      If ( allocated(error) ) Return

      metrics = nse_validation_metrics()
      failures = 0
      metrics%l1 = nearest(tolerances%l1,1.0_dp)
      Call apply_nse_validation_gates(metrics,tolerances,failures)
      Call require_flag(error,failures,nse_fail_l1)
      If ( allocated(error) ) Return

      metrics = nse_validation_metrics()
      failures = 0
      metrics%linf = nearest(tolerances%linf,1.0_dp)
      Call apply_nse_validation_gates(metrics,tolerances,failures)
      Call require_flag(error,failures,nse_fail_linf)
      If ( allocated(error) ) Return

      metrics = nse_validation_metrics()
      failures = 0
      metrics%dominant_error = nearest(tolerances%linf,1.0_dp)
      Call apply_nse_validation_gates(metrics,tolerances,failures)
      Call require_flag(error,failures,nse_fail_dominant)
      If ( allocated(error) ) Return
    EndDo

    Return
  End Subroutine test_tolerance_boundaries

  Subroutine test_scientific_guess_repeatability(error)
    Use nuclear_data, Only: nname
    Use xnet_nse, Only: knrtot, nse_solve, unse, xnse
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error

    Integer :: failures, state
    Real(dp) :: supplied_guess(2)
    Type(nse_validation_metrics) :: metrics

    Do state = 1, state_count
      Call solve_state(state)
      supplied_guess = unse(1:2) + (/ 2.0_dp, -2.0_dp /)
      Call nse_solve(state_rho(state),state_t9(state),state_ye(state),supplied_guess)
      Call evaluate_state(state,xnse,nname,metrics,failures)
      Call check(error,failures,0)
      If ( allocated(error) ) Return
      Call check(error,knrtot(3) > 0)
      If ( allocated(error) ) Return
    EndDo

    Return
  End Subroutine test_scientific_guess_repeatability

  Function get_mass_numbers() Result(values)
    Use nuclear_data, Only: aa
    Implicit None

    Real(dp) :: values(size(aa))

    values = aa

    Return
  End Function get_mass_numbers

  Function get_proton_numbers() Result(values)
    Use nuclear_data, Only: zz
    Implicit None

    Real(dp) :: values(size(zz))

    values = zz

    Return
  End Function get_proton_numbers

  Function get_neutron_numbers() Result(values)
    Use nuclear_data, Only: nn
    Implicit None

    Real(dp) :: values(size(nn))

    values = nn

    Return
  End Function get_neutron_numbers

  Integer Function find_name_index(name) Result(index)
    Use nuclear_data, Only: nname, ny
    Implicit None

    Character(5), Intent(in) :: name

    Integer :: inuc

    index = 0
    Do inuc = 1, ny
      If ( nname(inuc) == name ) Then
        index = inuc
        Exit
      EndIf
    EndDo

    Return
  End Function find_name_index

  Subroutine require_flag(error,failures,flag)
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error
    Integer, Intent(in) :: failures, flag

    Call check(error,iand(failures,flag) /= 0)

    Return
  End Subroutine require_flag

  Subroutine fail_fixture(message)
    Use, Intrinsic :: iso_fortran_env, Only: error_unit
    Implicit None

    Character(*), Intent(in) :: message

    Write(error_unit,'(a)') trim(message)
    Stop 1
  End Subroutine fail_fixture

End Module test_nse_scientific_validation

Program nse_scientific_validation_tests
  Use, Intrinsic :: iso_fortran_env, Only: error_unit
  Use test_nse_scientific_validation, Only: collect_nse_scientific_validation, &
    & load_validation_fixture
  Use testdrive, Only: new_testsuite, run_testsuite, testsuite_type
  Implicit None

  Character(512) :: network_directory, reference_path
  Integer :: stat
  Type(testsuite_type), Allocatable :: testsuites(:)

  Call get_command_argument(1,network_directory)
  Call get_command_argument(2,reference_path)
  If ( len_trim(network_directory) == 0 .OR. len_trim(reference_path) == 0 ) Then
    Write(error_unit,'(a)') &
      & 'usage: nse_scientific_validation_tests NETWORK_DIR REFERENCE_DATA'
    Stop 1
  EndIf
  Call load_validation_fixture(trim(network_directory),trim(reference_path))

  testsuites = [ new_testsuite('XNet independent NSE validation', &
    & collect_nse_scientific_validation) ]
  stat = 0
  Call run_testsuite(testsuites(1)%collect,error_unit,stat,parallel=.False.)
  If ( stat > 0 ) Stop 1

End Program nse_scientific_validation_tests
