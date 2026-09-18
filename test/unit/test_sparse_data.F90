Module sparse_data_fixture
  Use xnet_sparse, Only: sparse_data, read_sparse_ind
  Implicit None
  Private

  Integer, Parameter :: map_sizes(4) = (/ 3, 2, 1, 1 /)
  Integer, Parameter :: n10(3) = (/ 1, 2, 3 /)
  Integer, Parameter :: n11(3) = (/ 2, 3, 4 /)
  Integer, Parameter :: n20(2) = (/ 1, 2 /)
  Integer, Parameter :: n21(2) = (/ 1, 1 /)
  Integer, Parameter :: n22(2) = (/ 2, 3 /)
  Integer, Parameter :: n30(1) = 3
  Integer, Parameter :: n31(1) = 1
  Integer, Parameter :: n32(1) = 2
  Integer, Parameter :: n33(1) = 3
  Integer, Parameter :: n40(1) = 4
  Integer, Parameter :: n41(1) = 1
  Integer, Parameter :: n42(1) = 2
  Integer, Parameter :: n43(1) = 3
  Integer, Parameter :: n44(1) = 4
  Integer, Parameter :: ny = 4

  Character(256) :: work_directory = ''

  Public :: map_sizes
  Public :: n10, n11, n20, n21, n22, n30, n31, n32, n33
  Public :: n40, n41, n42, n43, n44
  Public :: ny
  Public :: read_fixture
  Public :: set_work_directory
  Public :: sparse_file_name
  Public :: write_fixture

Contains

  Subroutine set_work_directory(path)
    Implicit None

    Character(*), Intent(in) :: path

    work_directory = trim(path)

    Return
  End Subroutine set_work_directory

  Character(512) Function sparse_file_name(mutation)
    Implicit None

    Character(*), Intent(in) :: mutation

    sparse_file_name = trim(work_directory)//'/sparse_ind-'//trim(mutation)

    Return
  End Function sparse_file_name

  Subroutine read_fixture(mutation,data,status,io_status,message,requested_map_sizes)
    Implicit None

    Character(*), Intent(in) :: mutation
    Type(sparse_data), Intent(out) :: data
    Integer, Intent(out) :: status, io_status
    Character(*), Intent(out) :: message
    Integer, Intent(in), Optional :: requested_map_sizes(4)

    Integer :: effective_map_sizes(4)

    effective_map_sizes = map_sizes
    If ( present(requested_map_sizes) ) effective_map_sizes = requested_map_sizes
    Call read_sparse_ind(sparse_file_name(mutation),ny,effective_map_sizes,n10,n11,n20,n21, &
      & n22,n30,n31,n32,n33,n40,n41,n42,n43,n44,data,status,io_status,message)

    Return
  End Subroutine read_fixture

  Subroutine write_fixture(mutation)
    Implicit None

    Character(*), Intent(in) :: mutation

    Integer, Parameter :: lval = 13
    Integer :: cidx(lval), lun_sparse, lval_out, pb(ny+1), ridx(lval)
    Integer :: ns11(3), ns21(2), ns22(2), ns31(1), ns32(1), ns33(1)
    Integer :: ns41(1), ns42(1), ns43(1), ns44(1)

    ridx = (/ 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /)
    cidx = (/ 1, 2, 1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4 /)
    pb = (/ 1, 3, 6, 10, 14 /)
    ns11 = (/ 2, 5, 9 /)
    ns21 = (/ 1, 3 /)
    ns22 = (/ 2, 5 /)
    ns31 = 6
    ns32 = 7
    ns33 = 8
    ns41 = 10
    ns42 = 11
    ns43 = 12
    ns44 = 13
    lval_out = lval

    Open(newunit=lun_sparse,file=sparse_file_name(mutation),status='replace',form='unformatted')
    Select Case (trim(mutation))
    Case ('truncated-header')
      Close(lun_sparse)
      Return
    Case ('malformed-header')
      Write(lun_sparse) 'x'
      Close(lun_sparse)
      Return
    Case ('invalid-count')
      lval_out = ny*ny + 1
    End Select
    Write(lun_sparse) lval_out
    If ( trim(mutation) == 'invalid-count' ) Then
      Close(lun_sparse)
      Return
    EndIf

    Select Case (trim(mutation))
    Case ('truncated-topology')
      Write(lun_sparse) ridx(1:lval-1), cidx, pb
      Close(lun_sparse)
      Return
    Case ('pointer-start')
      pb(1) = 0
    Case ('pointer-terminal')
      pb(ny+1) = lval
    Case ('pointer-order')
      pb(3) = pb(2)
    Case ('coordinate-column')
      cidx(2) = ny + 1
    Case ('coordinate-row')
      ridx(2) = 2
    Case ('missing-diagonal')
      cidx(1:2) = (/ 2, 3 /)
    Case ('unordered-columns')
      cidx(1:2) = (/ 2, 1 /)
    End Select
    Write(lun_sparse) ridx, cidx, pb

    If ( trim(mutation) == 'map-dimensions' ) Then
      Write(lun_sparse) 3, 2, 1, 0
      Close(lun_sparse)
      Return
    EndIf
    Write(lun_sparse) map_sizes

    Select Case (trim(mutation))
    Case ('map-index')
      ns11(1) = 0
    Case ('map-coordinate-ns11')
      ns11(1) = 1
    Case ('map-coordinate-ns21')
      ns21(1) = 2
    Case ('map-coordinate-ns22')
      ns22(1) = 1
    Case ('map-coordinate-ns31')
      ns31(1) = 7
    Case ('map-coordinate-ns32')
      ns32(1) = 6
    Case ('map-coordinate-ns33')
      ns33(1) = 9
    Case ('map-coordinate-ns41')
      ns41(1) = 11
    Case ('map-coordinate-ns42')
      ns42(1) = 10
    Case ('map-coordinate-ns43')
      ns43(1) = 11
    Case ('map-coordinate-ns44')
      ns44(1) = 12
    Case ('truncated-map')
      Write(lun_sparse) ns11
      Close(lun_sparse)
      Return
    End Select
    Write(lun_sparse) ns11, ns21, ns22
    Write(lun_sparse) ns31
    Write(lun_sparse) ns32
    Write(lun_sparse) ns33
    Write(lun_sparse) ns41
    Write(lun_sparse) ns42
    Write(lun_sparse) ns43
    Write(lun_sparse) ns44
    Close(lun_sparse)

    Return
  End Subroutine write_fixture


End Module sparse_data_fixture

Module test_sparse_data
  Use sparse_data_fixture
  Use testdrive, Only: check, error_type, new_unittest, unittest_type
  Use xnet_sparse, Only: augment_crs_heat, read_sparse_ind, sparse_data, sparse_ind_invalid, &
    & sparse_ind_ok, sparse_ind_read_error
  Implicit None
  Private

  Public :: collect_sparse_data_tests

Contains

  Subroutine collect_sparse_data_tests(testsuite)
    Implicit None

    Type(unittest_type), Allocatable, Intent(out) :: testsuite(:)

    testsuite = [ &
      & new_unittest('valid sparse_ind schema',test_valid_sparse_ind), &
      & new_unittest('zero-length reaction maps',test_zero_length_maps), &
      & new_unittest('malformed and truncated header',test_invalid_header), &
      & new_unittest('incompatible dimensions',test_invalid_dimensions), &
      & new_unittest('truncated topology and maps',test_truncated_records), &
      & new_unittest('row pointer invariants',test_invalid_pointers), &
      & new_unittest('coordinate and topology invariants',test_invalid_topology), &
      & new_unittest('reaction-map invariants',test_invalid_maps), &
      & new_unittest('self-heating CRS augmentation and remapping',test_crs_heat), &
      & new_unittest('self-heating CRS rejects invalid input',test_crs_heat_rejection) ]

    Return
  End Subroutine collect_sparse_data_tests

  Subroutine test_valid_sparse_ind(error)
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error
    Type(sparse_data) :: data
    Character(256) :: message
    Integer :: io_status, status

    Call write_fixture('valid')
    Call read_fixture('valid',data,status,io_status,message)
    Call check(error,status,sparse_ind_ok)
    If ( allocated(error) ) Return
    Call check(error,io_status,0)
    If ( allocated(error) ) Return
    Call check(error,data%lval == 13 .and. &
      & all((/ data%l1s, data%l2s, data%l3s, data%l4s /) == map_sizes))
    If ( allocated(error) ) Return
    Call check(error,all(data%ridx == (/ 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /)) .and. &
      & all(data%cidx == (/ 1, 2, 1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4 /)))
    If ( allocated(error) ) Return
    Call check(error,all(data%pb == (/ 1, 3, 6, 10, 14 /)))
    If ( allocated(error) ) Return
    Call check(error,all(data%ns11 == (/ 2, 5, 9 /)) .and. &
      & all(data%ns21 == (/ 1, 3 /)) .and. all(data%ns22 == (/ 2, 5 /)))
    If ( allocated(error) ) Return
    Call check(error,data%ns31(1) == 6 .and. data%ns32(1) == 7 .and. data%ns33(1) == 8)
    If ( allocated(error) ) Return
    Call check(error,data%ns41(1) == 10 .and. data%ns42(1) == 11 .and. &
      & data%ns43(1) == 12 .and. data%ns44(1) == 13)

    Return
  End Subroutine test_valid_sparse_ind

  Subroutine test_zero_length_maps(error)
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error
    Type(sparse_data) :: data
    Character(512) :: file_name
    Character(256) :: message
    Integer, Allocatable :: empty(:)
    Integer :: io_status, lun_sparse, status

    Allocate (empty(0))
    file_name = sparse_file_name('zero-maps')
    Open(newunit=lun_sparse,file=trim(file_name),status='replace',form='unformatted')
    Write(lun_sparse) 2
    Write(lun_sparse) (/ 1, 2 /), (/ 1, 2 /), (/ 1, 2, 3 /)
    Write(lun_sparse) 0, 0, 0, 0
    Write(lun_sparse) empty, empty, empty
    Write(lun_sparse) empty
    Write(lun_sparse) empty
    Write(lun_sparse) empty
    Write(lun_sparse) empty
    Write(lun_sparse) empty
    Write(lun_sparse) empty
    Write(lun_sparse) empty
    Close(lun_sparse)
    Call read_sparse_ind(trim(file_name),2,(/ 0, 0, 0, 0 /),empty,empty,empty,empty,empty, &
      & empty,empty,empty,empty,empty,empty,empty,empty,empty,data,status,io_status,message)
    Call check(error,status,sparse_ind_ok)
    If ( allocated(error) ) Return
    Call check(error,data%lval == 2 .and. all(data%pb == (/ 1, 2, 3 /)))
    If ( allocated(error) ) Return
    Call check(error,size(data%ns11) == 0 .and. size(data%ns44) == 0)
    Deallocate (empty)

    Return
  End Subroutine test_zero_length_maps

  Subroutine test_invalid_header(error)
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error

    Call expect_failure(error,'truncated-header',sparse_ind_read_error,'header record')
    If ( allocated(error) ) Return
    Call expect_failure(error,'malformed-header',sparse_ind_read_error,'header record')

    Return
  End Subroutine test_invalid_header

  Subroutine test_invalid_dimensions(error)
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error
    Type(sparse_data) :: data
    Character(256) :: message
    Character(3), Parameter :: targets(15) = (/ 'ny ', 'n10', 'n11', 'n20', 'n21', 'n22', &
      & 'n30', 'n31', 'n32', 'n33', 'n40', 'n41', 'n42', 'n43', 'n44' /)
    Integer :: io_status, status, target

    Call expect_failure(error,'invalid-count',sparse_ind_invalid,'nonzero count')
    If ( allocated(error) ) Return
    Call expect_failure(error,'map-dimensions',sparse_ind_invalid,'map sizes')
    If ( allocated(error) ) Return

    Call read_fixture('unused',data,status,io_status,message, &
      & requested_map_sizes=(/ 3, 2, 1, -1 /))
    Call check(error,status,sparse_ind_invalid)
    If ( allocated(error) ) Return
    Call check(error,index(message,'dimensions are incompatible') > 0)
    If ( allocated(error) ) Return

    Do target = 1, size(targets)
      Call expect_caller_dimension_failure(error,trim(targets(target)))
      If ( allocated(error) ) Return
    EndDo

    Return
  End Subroutine test_invalid_dimensions

  Subroutine expect_caller_dimension_failure(error,mismatch)
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error
    Character(*), Intent(in) :: mismatch

    Type(sparse_data) :: data
    Character(256) :: message
    Integer, Allocatable :: n10_arg(:), n11_arg(:), n20_arg(:), n21_arg(:), n22_arg(:)
    Integer, Allocatable :: n30_arg(:), n31_arg(:), n32_arg(:), n33_arg(:)
    Integer, Allocatable :: n40_arg(:), n41_arg(:), n42_arg(:), n43_arg(:), n44_arg(:)
    Integer :: io_status, requested_ny, status

    requested_ny = ny
    If ( mismatch == 'ny' ) requested_ny = 0
    Allocate (n10_arg(size(n10)+merge(1,0,mismatch == 'n10')))
    Allocate (n11_arg(size(n11)+merge(1,0,mismatch == 'n11')))
    Allocate (n20_arg(size(n20)+merge(1,0,mismatch == 'n20')))
    Allocate (n21_arg(size(n21)+merge(1,0,mismatch == 'n21')))
    Allocate (n22_arg(size(n22)+merge(1,0,mismatch == 'n22')))
    Allocate (n30_arg(size(n30)+merge(1,0,mismatch == 'n30')))
    Allocate (n31_arg(size(n31)+merge(1,0,mismatch == 'n31')))
    Allocate (n32_arg(size(n32)+merge(1,0,mismatch == 'n32')))
    Allocate (n33_arg(size(n33)+merge(1,0,mismatch == 'n33')))
    Allocate (n40_arg(size(n40)+merge(1,0,mismatch == 'n40')))
    Allocate (n41_arg(size(n41)+merge(1,0,mismatch == 'n41')))
    Allocate (n42_arg(size(n42)+merge(1,0,mismatch == 'n42')))
    Allocate (n43_arg(size(n43)+merge(1,0,mismatch == 'n43')))
    Allocate (n44_arg(size(n44)+merge(1,0,mismatch == 'n44')))

    Call read_sparse_ind(sparse_file_name('unused'),requested_ny,map_sizes,n10_arg,n11_arg, &
      & n20_arg,n21_arg,n22_arg,n30_arg,n31_arg,n32_arg,n33_arg,n40_arg,n41_arg,n42_arg, &
      & n43_arg,n44_arg,data,status,io_status,message)
    Call check(error,status,sparse_ind_invalid)
    If ( allocated(error) ) Return
    If ( mismatch == 'ny' ) Then
      Call check(error,index(message,'network dimension must be positive') > 0)
    Else
      Call check(error,index(message,'dimensions are incompatible') > 0)
    EndIf

    Return
  End Subroutine expect_caller_dimension_failure

  Subroutine test_truncated_records(error)
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error

    Call expect_failure(error,'truncated-topology',sparse_ind_read_error,'topology record')
    If ( allocated(error) ) Return
    Call expect_failure(error,'truncated-map',sparse_ind_read_error,'reactant map record')

    Return
  End Subroutine test_truncated_records

  Subroutine test_invalid_pointers(error)
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error

    Call expect_failure(error,'pointer-start',sparse_ind_invalid,'pointer')
    If ( allocated(error) ) Return
    Call expect_failure(error,'pointer-terminal',sparse_ind_invalid,'pointer')
    If ( allocated(error) ) Return
    Call expect_failure(error,'pointer-order',sparse_ind_invalid, &
      & 'pointers are not strictly ordered')

    Return
  End Subroutine test_invalid_pointers

  Subroutine test_invalid_topology(error)
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error

    Call expect_failure(error,'coordinate-column',sparse_ind_invalid,'coordinate')
    If ( allocated(error) ) Return
    Call expect_failure(error,'coordinate-row',sparse_ind_invalid,'coordinate')
    If ( allocated(error) ) Return
    Call expect_failure(error,'missing-diagonal',sparse_ind_invalid,'diagonal')
    If ( allocated(error) ) Return
    Call expect_failure(error,'unordered-columns',sparse_ind_invalid, &
      & 'columns are not strictly ordered')

    Return
  End Subroutine test_invalid_topology

  Subroutine test_invalid_maps(error)
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error

    Call expect_failure(error,'map-index',sparse_ind_invalid,'index is out of range')
    If ( allocated(error) ) Return
    Call expect_failure(error,'map-coordinate-ns11',sparse_ind_invalid, &
      & 'ns11 does not resolve to its reaction coordinate')
    If ( allocated(error) ) Return
    Call expect_failure(error,'map-coordinate-ns21',sparse_ind_invalid, &
      & 'ns21 does not resolve to its reaction coordinate')
    If ( allocated(error) ) Return
    Call expect_failure(error,'map-coordinate-ns22',sparse_ind_invalid, &
      & 'ns22 does not resolve to its reaction coordinate')
    If ( allocated(error) ) Return
    Call expect_failure(error,'map-coordinate-ns31',sparse_ind_invalid, &
      & 'ns31 does not resolve to its reaction coordinate')
    If ( allocated(error) ) Return
    Call expect_failure(error,'map-coordinate-ns32',sparse_ind_invalid, &
      & 'ns32 does not resolve to its reaction coordinate')
    If ( allocated(error) ) Return
    Call expect_failure(error,'map-coordinate-ns33',sparse_ind_invalid, &
      & 'ns33 does not resolve to its reaction coordinate')
    If ( allocated(error) ) Return
    Call expect_failure(error,'map-coordinate-ns41',sparse_ind_invalid, &
      & 'ns41 does not resolve to its reaction coordinate')
    If ( allocated(error) ) Return
    Call expect_failure(error,'map-coordinate-ns42',sparse_ind_invalid, &
      & 'ns42 does not resolve to its reaction coordinate')
    If ( allocated(error) ) Return
    Call expect_failure(error,'map-coordinate-ns43',sparse_ind_invalid, &
      & 'ns43 does not resolve to its reaction coordinate')
    If ( allocated(error) ) Return
    Call expect_failure(error,'map-coordinate-ns44',sparse_ind_invalid, &
      & 'ns44 does not resolve to its reaction coordinate')

    Return
  End Subroutine test_invalid_maps


  Subroutine test_crs_heat(error)
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error
    Type(sparse_data) :: sparse_ind
    Type(sparse_data) :: sparse_ind_heat
    Character(256) :: message, mutation
    Integer :: io_status, length, status, variable_status

    Call write_fixture('augmentation')
    Call read_fixture('augmentation',sparse_ind,status,io_status,message)
    Call check(error,status,sparse_ind_ok)
    If ( allocated(error) ) Return
    Call augment_crs_heat(sparse_ind,ny,sparse_ind_heat,status,message)
    Call check(error,status,sparse_ind_ok)
    If ( allocated(error) ) Return
    Call check(error,sparse_ind_heat%lval == 22 .and. all((/ sparse_ind_heat%l1s, &
      & sparse_ind_heat%l2s, sparse_ind_heat%l3s, sparse_ind_heat%l4s /) == map_sizes))
    If ( allocated(error) ) Return

    mutation = ''
    Call get_environment_variable('XNET_CRS_AUGMENTATION_MUTATION',mutation, &
      & length=length,status=variable_status)
    If ( variable_status == 0 ) Then
      If ( mutation(1:length) == 'missing-temperature-entry' ) &
        & sparse_ind_heat%cidx(sparse_ind_heat%pb(2)-1) = 1
    EndIf

    Call check(error,all(sparse_ind_heat%ridx == &
      & (/ 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5 /)))
    If ( allocated(error) ) Return
    Call check(error,all(sparse_ind_heat%cidx == &
      & (/ 1, 2, 5, 1, 2, 3, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5 /)))
    If ( allocated(error) ) Return
    Call check(error,all(sparse_ind_heat%pb == (/ 1, 4, 8, 13, 18, 23 /)))
    If ( allocated(error) ) Return
    Call check(error,all(sparse_ind_heat%ns11 == (/ 2, 6, 11 /)))
    If ( allocated(error) ) Return
    Call check(error,all(sparse_ind_heat%ns21 == (/ 1, 4 /)) .and. &
      & all(sparse_ind_heat%ns22 == (/ 2, 6 /)))
    If ( allocated(error) ) Return
    Call check(error,sparse_ind_heat%ns31(1) == 8 .and. sparse_ind_heat%ns32(1) == 9 .and. &
      & sparse_ind_heat%ns33(1) == 10)
    If ( allocated(error) ) Return
    Call check(error,sparse_ind_heat%ns41(1) == 13 .and. sparse_ind_heat%ns42(1) == 14 .and. &
      & sparse_ind_heat%ns43(1) == 15 .and. sparse_ind_heat%ns44(1) == 16)

    Return
  End Subroutine test_crs_heat

  Subroutine test_crs_heat_rejection(error)
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error
    Type(sparse_data) :: sparse_ind
    Type(sparse_data) :: sparse_ind_heat
    Character(256) :: message
    Integer :: io_status, status

    Call write_fixture('augmentation-rejection')
    Call read_fixture('augmentation-rejection',sparse_ind,status,io_status,message)
    Call check(error,status,sparse_ind_ok)
    If ( allocated(error) ) Return
    sparse_ind%pb(2) = sparse_ind%pb(2) + 1
    Call augment_crs_heat(sparse_ind,ny,sparse_ind_heat,status,message)
    Call check(error,status,sparse_ind_invalid)
    If ( allocated(error) ) Return
    Call check(error,index(message,'declared row') > 0 .or. index(message,'strictly ordered') > 0)

    Return
  End Subroutine test_crs_heat_rejection

  Subroutine expect_failure(error,mutation,expected_status,diagnostic)
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error
    Character(*), Intent(in) :: mutation, diagnostic
    Integer, Intent(in) :: expected_status

    Type(sparse_data) :: data
    Character(256) :: message
    Integer :: io_status, status

    Call write_fixture(mutation)
    Call read_fixture(mutation,data,status,io_status,message)
    Call check(error,status,expected_status)
    If ( allocated(error) ) Return
    Call check(error,index(message,trim(diagnostic)) > 0)

    Return
  End Subroutine expect_failure

End Module test_sparse_data

Program sparse_data_test_runner
  Use, Intrinsic :: iso_fortran_env, Only: error_unit
  Use sparse_data_fixture, Only: set_work_directory
  Use test_sparse_data, Only: collect_sparse_data_tests
  Use testdrive, Only: new_testsuite, run_testsuite, testsuite_type
  Implicit None

  Character(256) :: work_directory
  Integer :: stat
  Type(testsuite_type), Allocatable :: testsuites(:)

  If ( command_argument_count() /= 1 ) Then
    Write(error_unit,*) 'usage: sparse data test WORK_DIRECTORY'
    Stop 1
  EndIf
  Call get_command_argument(1,work_directory)
  Call set_work_directory(trim(work_directory))

  stat = 0
  testsuites = [ new_testsuite('sparse data',collect_sparse_data_tests) ]
  Write(error_unit,'("# Testing: ",a)') testsuites(1)%name
  Call run_testsuite(testsuites(1)%collect,error_unit,stat,parallel=.False.)
  If ( stat > 0 ) Then
    Write(error_unit,'(i0,1x,a)') stat,'test(s) failed'
    Stop 1
  EndIf
End Program sparse_data_test_runner
