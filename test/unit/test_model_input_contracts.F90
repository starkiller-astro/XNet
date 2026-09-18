Module test_model_input_contracts
  Use testdrive, Only: check, error_type, new_unittest, unittest_type
  Use xnet_types, Only: dp
  Implicit None
  Private

  Public :: collect_model_input_contracts

Contains

  Subroutine collect_model_input_contracts(testsuite)
    Implicit None

    Type(unittest_type), Allocatable, Intent(out) :: testsuite(:)

    testsuite = [ &
      & new_unittest('initial abundance reader missing file outputs', test_missing_file_outputs), &
      & new_unittest('initial abundance reader valid file', test_valid_file) ]

    Return
  End Subroutine collect_model_input_contracts

  Subroutine test_missing_file_outputs(error)
    Use model_input_ascii, Only: read_inab_file
    Use nuclear_data, Only: initialize_nuclear_fixture
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error

    Character(80) :: abund_desc
    Integer :: i, ierr
    Real(dp) :: aext, xext, xnet, yein, yin(2), zext

    Call initialize_nuclear_fixture
    abund_desc = 'unset'
    yein = -1.0_dp
    yin = -1.0_dp
    xext = -1.0_dp
    aext = -1.0_dp
    zext = -1.0_dp
    xnet = -1.0_dp
    ierr = 0
    Call read_inab_file('missing-initial-abundance-file',abund_desc,yein,yin,xext,aext,zext,xnet,ierr)

    Call check(error,ierr /= 0)
    If ( allocated(error) ) Return
    Call check(error,len_trim(abund_desc),0)
    If ( allocated(error) ) Return
    Call check(error,yein,0.0_dp)
    If ( allocated(error) ) Return
    Do i = 1, size(yin)
      Call check(error,yin(i),0.0_dp)
      If ( allocated(error) ) Return
    EndDo
    Call check(error,xext,0.0_dp)
    If ( allocated(error) ) Return
    Call check(error,aext,1.0_dp)
    If ( allocated(error) ) Return
    Call check(error,zext,0.0_dp)
    If ( allocated(error) ) Return
    Call check(error,xnet,0.0_dp)

    Return
  End Subroutine test_missing_file_outputs

  Subroutine test_valid_file(error)
    Use model_input_ascii, Only: read_inab_file
    Use nuclear_data, Only: initialize_nuclear_fixture
    Implicit None

    Type(error_type), Allocatable, Intent(out) :: error

    Character(80) :: abund_desc
    Integer :: i, ierr, lun
    Real(dp) :: aext, expected_yin(2), xext, xnet, yein, yin(2), zext

    Call initialize_nuclear_fixture
    Open(newunit=lun,file='model-input-valid-abundances',action='write',status='replace')
    Write(lun,"(a)") 'valid initial abundances'
    Write(lun,"(a)") 'h1 0.5 he4 0.125 ye 0.5 xx 0.0'
    Close(lun)

    Call read_inab_file('model-input-valid-abundances',abund_desc,yein,yin,xext,aext,zext,xnet,ierr)

    Call check(error,ierr,0)
    If ( allocated(error) ) Return
    Call check(error,trim(abund_desc),'valid initial abundances')
    If ( allocated(error) ) Return
    Call check(error,yein,0.5_dp)
    If ( allocated(error) ) Return
    expected_yin = (/ 0.5_dp, 0.125_dp /)
    Do i = 1, size(yin)
      Call check(error,yin(i),expected_yin(i))
      If ( allocated(error) ) Return
    EndDo
    Call check(error,xext,0.0_dp)
    If ( allocated(error) ) Return
    Call check(error,aext,1.0_dp)
    If ( allocated(error) ) Return
    Call check(error,zext,0.0_dp)
    If ( allocated(error) ) Return
    Call check(error,xnet,1.0_dp)

    Open(newunit=lun,file='model-input-valid-abundances',status='old')
    Close(lun,status='delete')

    Return
  End Subroutine test_valid_file

End Module test_model_input_contracts

Program model_input_contract_tests
  Use, Intrinsic :: iso_fortran_env, Only: error_unit
  Use test_model_input_contracts, Only: collect_model_input_contracts
  Use testdrive, Only: new_testsuite, run_testsuite, testsuite_type
  Implicit None

  Integer :: stat
  Type(testsuite_type), Allocatable :: testsuites(:)

  stat = 0
  testsuites = [ new_testsuite('model input contracts',collect_model_input_contracts) ]
  Call run_testsuite(testsuites(1)%collect,error_unit,stat,parallel=.False.)
  If ( stat > 0 ) Stop 1
End Program model_input_contract_tests
