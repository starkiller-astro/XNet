!***************************************************************************************************
! xnet_util.f90 10/18/17
! This file contains various utility routines that are common throughout XNet.
!***************************************************************************************************

Module xnet_util
  !-------------------------------------------------------------------------------------------------
  ! This module contains the data and routines for an assortment of (mostly independent) utility
  ! routines.
  !-------------------------------------------------------------------------------------------------
  Use xnet_constants, Only: ln_2
  Use xnet_types, Only: dp
  Implicit None

  Real(dp), Parameter :: exp_max = maxexponent(1.0_dp)*ln_2*0.99_dp
  Real(dp), Parameter :: exp_min = minexponent(1.0_dp)*ln_2*0.99_dp

  Interface safe_exp
    Module Procedure safe_exp_scalar
    Module Procedure safe_exp_vector
  End Interface safe_exp

  Interface readnext
    Module Procedure readnext_i
    Module Procedure readnext_r
    Module Procedure readnext_d
  End Interface readnext

Contains

  Integer Function ifactorial(n)
    !-----------------------------------------------------------------------------------------------
    ! This function computes factorials.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Integer, Intent(in) :: n

    ! Local variable
    Integer :: i

    ! Compute factorial of n
    ifactorial = 1
    Do i = 1, n
      ifactorial = ifactorial*i
    EndDo

    Return
  End Function ifactorial

  Function safe_exp_scalar( x ) Result( y )
    !-----------------------------------------------------------------------------------------------
    ! This routine safely calculates e^{x} with x constrained to prevent overflow and underflow.
    !-----------------------------------------------------------------------------------------------
    Use xnet_types, Only: dp
    Implicit None
    !$acc routine seq

    ! Input variables
    Real(dp), Intent(in) :: x

    ! Function variable
    Real(dp) :: y

    y = exp( min( exp_max, max( exp_min, x ) ) )

    Return
  End Function safe_exp_scalar

  Function safe_exp_vector( x ) Result( y )
    !-----------------------------------------------------------------------------------------------
    ! This routine safely calculates e^{x} with x constrained to prevent overflow and underflow.
    !-----------------------------------------------------------------------------------------------
    Use xnet_types, Only: dp
    Implicit None
    !$acc routine seq

    ! Input variables
    Real(dp), Intent(in) :: x(:)

    ! Function variable
    Real(dp) :: y(size(x))

    y = exp( min( exp_max, max( exp_min, x ) ) )

    Return
  End Function safe_exp_vector

  Integer Function getNewUnit(unit)
    !-----------------------------------------------------------------------------------------------
    ! Get a free unit number within range 7-999.
    !-----------------------------------------------------------------------------------------------
    Implicit None
    Integer, Intent(out), Optional :: unit
    Logical :: connected
    Integer :: number
    getNewUnit = 0
    Do number = 7,999
      Inquire(unit=number, opened=connected)
      If ( .not. connected ) Then
        getNewUnit = number
        Exit
      EndIf
    EndDo
    If ( present(unit) ) unit = getNewUnit
    Return
  End Function getNewUnit

  Subroutine replace_tabs(string)
    !-----------------------------------------------------------------------------------------------
    ! Replace all tabs in a string with spaces
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input/Output variables
    Character(*), Intent(inout) :: string

    ! Local variables
    Integer, Parameter :: iachar_tab = 9
    Integer :: i, iachar_buffer

    Do i = 1, len_trim(string)
      iachar_buffer = iachar(string(i:i))
      If ( iachar_buffer == iachar_tab ) string(i:i) = ' '
    EndDo

    Return
  End Subroutine replace_tabs

  Subroutine readnext_i(string, pos, xx)
    !-----------------------------------------------------------------------------------------------
    ! Given a starting position in a string, skip any whitespace and read the next number as a int.
    ! Return the value (zero if non-existant) and the new position in the string.
    !-----------------------------------------------------------------------------------------------
    Implicit None
    Character(*), Intent(in) :: string
    Integer, Intent(inout) :: pos
    Integer, Intent(out) :: xx

    Integer :: i1, i2

    i2 = len_trim(string)

    ! initial values
    If ( pos > i2 ) Then
      pos = 0
      xx = 0
      Return
    EndIf

    ! skip blanks
    i1 = pos
    Do
      If ( string(i1:i1) /= ' ' ) Exit
      i1 = i1 + 1
    EndDo

    ! read real value and set pos
    Read(string(i1:i2), *) xx
    pos = scan(string(i1:i2), ' ')
    If ( pos == 0 ) Then
      pos = i2 + 1
    Else
      pos = pos + i1 - 1
    EndIf

    Return
  End Subroutine readnext_i

  Subroutine readnext_r(string, pos, xx)
    !-----------------------------------------------------------------------------------------------
    ! Given a starting position in a string, skip any whitespace and read the next number as a real.
    ! Return the value (zero if non-existant) and the new position in the string.
    !-----------------------------------------------------------------------------------------------
    Use xnet_types, Only: sp
    Implicit None
    Character(*), Intent(in) :: string
    Integer, Intent(inout) :: pos
    Real(sp), Intent(out) :: xx

    Integer :: i1, i2

    i2 = len_trim(string)

    ! initial values
    If ( pos < 1 .or. pos > i2 ) Then
      pos = 0
      xx = 0.0
      Return
    EndIf

    ! skip blanks
    i1 = pos
    Do
      If ( string(i1:i1) /= ' ' ) Exit
      i1 = i1 + 1
    EndDo

    ! read real value and set pos
    Read(string(i1:i2), *) xx
    pos = scan(string(i1:i2), ' ')
    If ( pos == 0 ) Then
      pos = i2 + 1
    Else
      pos = pos + i1 - 1
    EndIf

    Return
  End Subroutine readnext_r

  Subroutine readnext_d(string, pos, xx)
    !-----------------------------------------------------------------------------------------------
    ! Given a starting position in a string, skip any whitespace and read the next number as a dbl.
    ! Return the value (zero if non-existant) and the new position in the string.
    !-----------------------------------------------------------------------------------------------
    Use xnet_types, Only: dp
    Implicit None
    Character(*), Intent(in) :: string
    Integer, Intent(inout) :: pos
    Real(dp), Intent(out) :: xx

    Integer :: i1, i2

    i2 = len_trim(string)

    ! initial values
    If ( pos > i2 ) Then
      pos = 0
      xx = 0.0
      Return
    EndIf

    ! skip blanks
    i1 = pos
    Do
      If ( string(i1:i1) /= ' ' ) Exit
      i1 = i1 + 1
    EndDo

    ! read real value and set pos
    Read(string(i1:i2), *) xx
    pos = scan(string(i1:i2), ' ')
    If ( pos == 0 ) Then
      pos = i2 + 1
    Else
      pos = pos + i1 - 1
    EndIf

    Return
  End Subroutine readnext_d

  Subroutine norm(yy,aa)
    !-----------------------------------------------------------------------------------------------
    ! This routine renormalizes the abundances to guarantee mass conservation.
    !-----------------------------------------------------------------------------------------------
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: aa(:)

    ! Input/Output variables
    Real(dp), Intent(inout) :: yy(:)

    ! Local variables
    Real(dp) :: xtot, rxt

    ! Renormalize total mass fraction to 1
    xtot = sum(yy * aa)
    rxt  = 1.0 / xtot
    yy   = yy * rxt

    Return
  End Subroutine norm

  Subroutine ye_norm(yy,ye,zz,nn,aa)
    !-----------------------------------------------------------------------------------------------
    ! This routine renormalizes the abundances to guarantee mass and charge conservation.
    !-----------------------------------------------------------------------------------------------
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: ye
    Real(dp), Intent(in) :: zz(:)
    Real(dp), Intent(in) :: nn(:)
    Real(dp), Intent(in) :: aa(:)

    ! Input/Output variables
    Real(dp), Intent(inout) :: yy(:)

    ! Local variables
    Real(dp) :: zy, nny, zny, zzy, beta, alph

    ! Calculate moments needed to calculate normalization coefficients
    nny  = sum(nn * yy)
    zy   = sum(zz * yy)
    zny  = sum(nn * zz * yy / aa)
    zzy  = sum(zz * zz * yy / aa)
    beta = (ye*nny - zny) / (nny*zzy - zy*zny)
    alph = (1.0 - beta*zy) / nny
    yy   = yy * (alph*nn + beta*zz) / aa

    Return
  End Subroutine ye_norm

  Subroutine plasma(t9,rho,ytot,ye,zbar,zibar,ztilde,zinter,lambda0,gammae)
    !-----------------------------------------------------------------------------------------------
    ! This routine calculate various plasma quantities
    !-----------------------------------------------------------------------------------------------
    Use xnet_constants, Only: avn, bok, e2, pi, third, thbim2, twm2bi
    Use xnet_types, Only: dp
    Implicit None
    !$acc routine seq

    ! Input variables
    Real(dp), Intent(in) :: t9, rho, ytot, ye, zbar, zibar, ztilde

    ! Output variables
    Real(dp), Intent(out) :: zinter, lambda0, gammae

    ! Local variables
    Real(dp) :: ae, bkt, nb, ni, ne

    bkt = bok*t9
    nb = avn*rho
    ni = nb*ytot
    ne = nb*ye
    lambda0 = sqrt(4.0*pi*ni) * (e2/bkt)**1.5 ! DGC, Eq. 3
    ae = (3.0 / (4.0*pi*ne))**third ! electron-sphere radius
    gammae = e2 / (ae*bkt) ! electron Coulomb coupling parameter
    zinter = zibar / (ztilde**thbim2 * zbar**twm2bi) ! GDC, Table 4

    Return
  End Subroutine plasma

  Subroutine name_ordered(base_string,n,nmax)
    !-----------------------------------------------------------------------------------------------
    ! This routine appends the integer n (padded with "0"s up to the size of the integer nmax) to
    ! the character variable base_string.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input variables
    Integer, Intent(in) :: n, nmax

    ! Input/Output variables
    Character(*), Intent(inout) :: base_string

    ! Local variables
    Character(1) :: nmax_digits_string
    Character(6) :: n_format
    Character(9) :: n_string

    ! Find character length of imax
    Write(nmax_digits_string,"(i1)") int(log10(real(nmax))) + 1

    ! Construct format spec and write n as zero padded string
    n_format = "(i"//nmax_digits_string//"."//nmax_digits_string//")"
    Write(n_string,n_format) n

    ! Append n_string to base_string
    base_string = trim(base_string)//trim(n_string)

    Return
  End Subroutine name_ordered

  Subroutine string_lc(string)
    !-----------------------------------------------------------------------------------------------
    ! This routine converts an ASCII string to all lower case.
    !-----------------------------------------------------------------------------------------------
    Implicit None

    ! Input/Output variables
    Character(*), Intent(inout) :: string

    ! Local variables
    Integer, Parameter :: lc_a_ascii=iachar('a')
    Integer, Parameter :: uc_a_ascii=iachar('A')
    Integer, Parameter :: uc_z_ascii=iachar('Z')
    Integer :: i, x

    Do i = 1, len_trim(string)
      x = iachar(string(i:i))
      If ( x >= uc_a_ascii .and. x <= uc_z_ascii ) Then
        x = x + (lc_a_ascii - uc_a_ascii)
        string(i:i) = achar(x)
      EndIf
    EndDo

    Return
  End Subroutine string_lc

  Subroutine xnet_terminate(c_diagnostic,i_diagnostic)
    !-----------------------------------------------------------------------------------------------
    ! This routine gracefully exits XNet with a diagnostic statement in the event of an error.
    !-----------------------------------------------------------------------------------------------
    Use, Intrinsic :: iso_fortran_env, Only: lun_stdout=>output_unit
    Use xnet_parallel, Only: parallel_abort
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: c_diagnostic
    Integer, Intent(in), Optional :: i_diagnostic

    ! Print the diagnostic statement
    If ( present(i_diagnostic) ) Then
      Call parallel_abort(c_diagnostic, i_diagnostic)
    Else
      Call parallel_abort(c_diagnostic)
    EndIf

    Stop
  End Subroutine xnet_terminate

End Module xnet_util
