Program unit_test_main
  Use, Intrinsic :: iso_fortran_env, Only: error_unit
  Use testdrive, Only: new_testsuite, run_testsuite, testsuite_type
  Use test_input_parsing, Only: collect_input_parsing
  Use test_name_ordered, Only: collect_name_ordered
  Use test_safe_exp, Only: collect_safe_exp
  Use test_types, Only: collect_types
  Implicit None

  Integer :: i, stat
  Type(testsuite_type), Allocatable :: suites(:)

  suites = [ &
    new_testsuite("types", collect_types), &
    new_testsuite("safe_exp", collect_safe_exp), &
    new_testsuite("input parsing", collect_input_parsing), &
    new_testsuite("generated filename suffixes", collect_name_ordered) &
  ]

  stat = 0
  Do i = 1, Size(suites)
    Write(error_unit, '("# Testing:", 1x, a)') suites(i)%name
    Call run_testsuite(suites(i)%collect, error_unit, stat)
  EndDo

  If ( stat > 0 ) Then
    Write(error_unit, '(i0, 1x, a)') stat, "test(s) failed"
    Error Stop 1
  EndIf
End Program unit_test_main
