Module xnet_gpu
  Use, Intrinsic :: iso_c_binding, Only: c_loc, c_ptr
  Use xnet_types, Only: dp
  Implicit None

  Interface dev_ptr
    Module Procedure dev_ptr_int
    Module Procedure dev_ptr_dp
  End Interface dev_ptr

Contains

  Function dev_ptr_int(value) Result(pointer)
    Implicit None

    Integer, Intent(in), Target :: value
    Type(c_ptr) :: pointer

    pointer = c_loc(value)

    Return
  End Function dev_ptr_int

  Function dev_ptr_dp(value) Result(pointer)
    Implicit None

    Real(dp), Intent(in), Target :: value
    Type(c_ptr) :: pointer

    pointer = c_loc(value)

    Return
  End Function dev_ptr_dp

End Module xnet_gpu

Module xnet_linalg
  Use xnet_types, Only: dp
  Implicit None

Contains

  Subroutine LinearSolve_CPU(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
    Implicit None

    Character, Intent(in) :: trans
    Integer, Intent(in) :: n, nrhs, lda, ldb
    Real(dp), Intent(in) :: a(lda,*)
    Integer, Intent(inout) :: ipiv(*)
    Real(dp), Intent(inout) :: b(ldb,*)
    Integer, Intent(out) :: info

    Real(dp) :: work(n,n)

    If ( trans /= 'N' .or. nrhs /= 1 ) Then
      info = -1
      Return
    EndIf
    work = a(1:n,1:n)
    Call dgesv(n,nrhs,work,n,ipiv,b,ldb,info)

    Return
  End Subroutine LinearSolve_CPU

  Subroutine LUDecomp_CPU(m,n,a,lda,ipiv,info)
    Implicit None

    Integer, Intent(in) :: m, n, lda
    Real(dp), Intent(inout) :: a(lda,*)
    Integer, Intent(out) :: ipiv(*), info

    Integer :: i

    Do i = 1, min(m,n)
      ipiv(i) = i
    EndDo
    info = 0

    Return
  End Subroutine LUDecomp_CPU

  Subroutine LUBksub_CPU(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
    Implicit None

    Character, Intent(in) :: trans
    Integer, Intent(in) :: n, nrhs, lda, ldb
    Real(dp), Intent(in) :: a(lda,*)
    Integer, Intent(inout) :: ipiv(*)
    Real(dp), Intent(inout) :: b(ldb,*)
    Integer, Intent(out) :: info

    Call LinearSolve_CPU(trans,n,nrhs,a,lda,ipiv,b,ldb,info)

    Return
  End Subroutine LUBksub_CPU

  Subroutine LinearSolveBatched_GPU
    Implicit None

    Return
  End Subroutine LinearSolveBatched_GPU

  Subroutine LUDecompBatched_GPU
    Implicit None

    Return
  End Subroutine LUDecompBatched_GPU

  Subroutine LUBksubBatched_GPU
    Implicit None

    Return
  End Subroutine LUBksubBatched_GPU

End Module xnet_linalg
