module matrix_util
  use, intrinsic :: iso_fortran_env, only: dp=>real64
  implicit none
  private

  public :: matrix_cond

contains

  subroutine matrix_cond(n, A, cond)
    implicit none
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: A(n,n)
    real(dp), intent(out) :: cond

    real(dp) :: work(4*n)
    integer  :: ipiv(n), iwork(n)
    real(dp) :: anorm, rcond
    integer  :: info
    real(dp) :: dlange
    external :: dgetrf, dgecon, dlange

    ! Get inf-norm of A
    anorm = dlange('I', n, n, A, n, work)

    ! LU factorization
    call dgetrf(n, n, A, n, ipiv, info)
    if (info /= 0) then
      print *, "DGETRF failed with info = ", info
      cond = -1.0
      return
    end if

    ! Estimate condition number
    call dgecon('I', n, A, n, anorm, rcond, work, iwork, info)
    if (info /= 0) then
      print *, "DGECON failed with info = ", info
      cond = -1.0
      return
    end if

    cond = 1.0_dp / rcond

    return
  end subroutine matrix_cond

end module matrix_util