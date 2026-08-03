Module test_types
  Use, Intrinsic :: iso_fortran_env, Only: int32, int64, real32, real64
  Use testdrive, Only: check, error_type, new_unittest, unittest_type
  Use xnet_types, Only: dp, i4, i8, sp
  Implicit None
  Private

  Public :: collect_types

Contains

  Subroutine collect_types(tests)
    Type(unittest_type), Allocatable, Intent(out) :: tests(:)

    tests = [ &
      new_unittest("kind values", test_kind_values), &
      new_unittest("storage sizes", test_storage_sizes) &
    ]
  End Subroutine collect_types

  Subroutine test_kind_values(error)
    Type(error_type), Allocatable, Intent(out) :: error

    Call check(error, dp == real64, "dp must equal iso_fortran_env real64")
    If ( Allocated(error) ) Return
    Call check(error, sp == real32, "sp must equal iso_fortran_env real32")
    If ( Allocated(error) ) Return
    Call check(error, i4 == int32, "i4 must equal iso_fortran_env int32")
    If ( Allocated(error) ) Return
    Call check(error, i8 == int64, "i8 must equal iso_fortran_env int64")
  End Subroutine test_kind_values

  Subroutine test_storage_sizes(error)
    Type(error_type), Allocatable, Intent(out) :: error

    Call check(error, Storage_Size(0.0_dp) == 64, "real(dp) must occupy 64 bits")
    If ( Allocated(error) ) Return
    Call check(error, Storage_Size(0.0_sp) == 32, "real(sp) must occupy 32 bits")
    If ( Allocated(error) ) Return
    Call check(error, Storage_Size(0_i4) == 32, "integer(i4) must occupy 32 bits")
    If ( Allocated(error) ) Return
    Call check(error, Storage_Size(0_i8) == 64, "integer(i8) must occupy 64 bits")
  End Subroutine test_storage_sizes

End Module test_types
