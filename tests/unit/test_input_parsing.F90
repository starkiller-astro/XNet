Module test_input_parsing
  Use testdrive, Only: check, error_type, new_unittest, unittest_type
  Use xnet_types, Only: dp
  Use xnet_util, Only: readnext, replace_tabs, string_lc
  Implicit None
  Private

  Public :: collect_input_parsing

Contains

  Subroutine collect_input_parsing(tests)
    Type(unittest_type), Allocatable, Intent(out) :: tests(:)

    tests = [ &
      new_unittest("tab and blank separated values", test_tab_and_blank_values), &
      new_unittest("ASCII case folding", test_ascii_case_folding) &
    ]
  End Subroutine collect_input_parsing

  Subroutine test_ascii_case_folding(error)
    Type(error_type), Allocatable, Intent(out) :: error
    Character(4) :: name

    name = "Fe56"
    Call string_lc(name)

    Call check(error, name == "fe56", "ASCII case folding must lowercase F and preserve e56")
  End Subroutine test_ascii_case_folding

  Subroutine test_tab_and_blank_values(error)
    Type(error_type), Allocatable, Intent(out) :: error
    Character(64) :: line
    Integer :: pos
    Real(dp) :: value

    line = "  1.25" // Achar(9) // " -2.5" // Achar(9) // "3.75"
    Call replace_tabs(line)
    pos = 1

    Call readnext(line, pos, value)
    Call check(error, value == 1.25_dp, "first token must be read after leading blanks")
    If ( Allocated(error) ) Return

    Call readnext(line, pos, value)
    Call check(error, value == -2.5_dp, "second token must be read after a replaced tab")
    If ( Allocated(error) ) Return

    Call readnext(line, pos, value)
    Call check(error, value == 3.75_dp, "third token must be read after a replaced tab")
    If ( Allocated(error) ) Return

    Call readnext(line, pos, value)
    Call check(error, value == 0.0_dp, "missing token must return zero")
    If ( Allocated(error) ) Return
    Call check(error, pos == 0, "missing token must set position to zero")
  End Subroutine test_tab_and_blank_values

End Module test_input_parsing
