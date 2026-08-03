Module test_name_ordered
  Use testdrive, Only: check, error_type, new_unittest, unittest_type
  Use xnet_util, Only: name_ordered
  Implicit None
  Private

  Public :: collect_name_ordered

Contains

  Subroutine collect_name_ordered(tests)
    Type(unittest_type), Allocatable, Intent(out) :: tests(:)

    tests = [ &
      new_unittest("generated filename suffixes", test_generated_filename_suffixes) &
    ]
  End Subroutine collect_name_ordered

  Subroutine test_generated_filename_suffixes(error)
    Type(error_type), Allocatable, Intent(out) :: error
    Character(16) :: name

    name = "zone"
    Call name_ordered(name, 1, 1)
    Call check(error, Trim(name) == "zone1", &
      "nmax=1 must append the literal one-digit suffix")
    If ( Allocated(error) ) Return

    name = "diag"
    Call name_ordered(name, 3, 12)
    Call check(error, Trim(name) == "diag03", &
      "nmax=12 must append the literal two-digit padded suffix")
  End Subroutine test_generated_filename_suffixes

End Module test_name_ordered
