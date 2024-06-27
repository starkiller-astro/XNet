module eos_type_module

  use xnet_types, only: dp
  use xnet_util, only: xnet_terminate

  implicit none

  real(dp), parameter, private :: ONE = 1.0_dp
  real(dp), parameter, private :: ZERO = 0.0_dp
  real(dp), parameter, private :: small_x = 0.0_dp

  integer, parameter :: eos_input_rt = 1  ! rho, T are inputs
  integer, parameter :: eos_input_rh = 2  ! rho, h are inputs
  integer, parameter :: eos_input_tp = 3  ! T, p are inputs
  integer, parameter :: eos_input_rp = 4  ! rho, p are inputs
  integer, parameter :: eos_input_re = 5  ! rho, e are inputs
  integer, parameter :: eos_input_ps = 6  ! p, s are inputs
  integer, parameter :: eos_input_ph = 7  ! p, h are inputs
  integer, parameter :: eos_input_th = 8  ! T, h are inputs

  ! these are used to allow for a generic interface to the
  ! root finding
  integer, parameter :: itemp = 1
  integer, parameter :: idens = 2
  integer, parameter :: iener = 3
  integer, parameter :: ienth = 4
  integer, parameter :: ientr = 5
  integer, parameter :: ipres = 6

  ! error codes
  integer, parameter :: ierr_general         = 1
  integer, parameter :: ierr_input           = 2
  integer, parameter :: ierr_iter_conv       = 3
  integer, parameter :: ierr_neg_e           = 4
  integer, parameter :: ierr_neg_p           = 5
  integer, parameter :: ierr_neg_h           = 6
  integer, parameter :: ierr_neg_s           = 7
  integer, parameter :: ierr_iter_var        = 8
  integer, parameter :: ierr_init            = 9
  integer, parameter :: ierr_init_xn         = 10
  integer, parameter :: ierr_out_of_bounds   = 11
  integer, parameter :: ierr_not_implemented = 12

  ! Minimum and maximum thermodynamic quantities permitted by the EOS.

  real(dp), allocatable :: mintemp
  real(dp), allocatable :: maxtemp
  real(dp), allocatable :: mindens
  real(dp), allocatable :: maxdens
  real(dp), allocatable :: minx
  real(dp), allocatable :: maxx
  real(dp), allocatable :: minye
  real(dp), allocatable :: maxye
  real(dp), allocatable :: mine
  real(dp), allocatable :: maxe
  real(dp), allocatable :: minp
  real(dp), allocatable :: maxp
  real(dp), allocatable :: mins
  real(dp), allocatable :: maxs
  real(dp), allocatable :: minh
  real(dp), allocatable :: maxh

  !$acc declare &
  !$acc create(mintemp, maxtemp, mindens, maxdens, minx, maxx, minye, maxye) &
  !$acc create(mine, maxe, minp, maxp, mins, maxs, minh, maxh)

  ! A generic structure holding thermodynamic quantities and their derivatives,
  ! plus some other quantities of interest.

  ! rho      -- mass density (g/cm**3)
  ! T        -- temperature (K)
  ! p        -- the pressure (dyn/cm**2)
  ! h        -- the enthalpy (erg/g)
  ! e        -- the internal energy (erg/g)
  ! s        -- the entropy (erg/g/K)
  ! c_v      -- specific heat at constant volume
  ! c_p      -- specific heat at constant pressure
  ! ne       -- number density of electrons + positrons
  ! np       -- number density of positrons only
  ! eta      -- degeneracy parameter
  ! pele     -- electron pressure + positron pressure
  ! ppos     -- position pressure only
  ! mu       -- mean molecular weight
  ! mu_e     -- mean number of nucleons per electron
  ! y_e      -- electron fraction == 1 / mu_e
  ! dPdT     -- d pressure/ d temperature
  ! dPdr     -- d pressure/ d density
  ! dedT     -- d energy/ d temperature
  ! dedr     -- d energy/ d density
  ! dsdT     -- d entropy/ d temperature
  ! dsdr     -- d entropy/ d density
  ! dhdT     -- d enthalpy/ d temperature
  ! dhdr     -- d enthalpy/ d density
  ! dPdX     -- d pressure / d xmass
  ! dhdX     -- d enthalpy / d xmass at constant pressure
  ! gam1     -- first adiabatic index (d log P/ d log rho) |_s
  ! cs       -- sound speed
  ! abar     -- average atomic number ( sum_k {X_k} ) / ( sum_k {X_k/A_k} )
  ! zbar     -- average proton number ( sum_k {Z_k X_k/ A_k} ) / ( sum_k {X_k/A_k} )
  ! dpdA     -- d pressure/ d abar
  ! dpdZ     -- d pressure/ d zbar
  ! dedA     -- d energy/ d abar
  ! dedZ     -- d energy/ d zbar
  ! dpde     -- d pressure / d energy |_rho
  ! dpdr_e   -- d pressure / d rho |_energy
  ! conductivity -- thermal conductivity (in erg/cm/K/sec)

  type :: eos_t

    real(dp) :: rho
    real(dp) :: T
    real(dp) :: p
    real(dp) :: e
    real(dp) :: h
    real(dp) :: s

    real(dp) :: dpdT
    real(dp) :: dpdr
    real(dp) :: dedT
    real(dp) :: dedr
    real(dp) :: dhdT
    real(dp) :: dhdr
    real(dp) :: dsdT
    real(dp) :: dsdr
    real(dp) :: dpde
    real(dp) :: dpdr_e

    real(dp) :: cv
    real(dp) :: cp
    real(dp) :: xne
    real(dp) :: xnp
    real(dp) :: eta
    real(dp) :: detadt
    real(dp) :: pele
    real(dp) :: ppos
    real(dp) :: mu
    real(dp) :: mu_e
    real(dp) :: y_e
    real(dp) :: gam1
    real(dp) :: cs

    real(dp) :: abar
    real(dp) :: zbar

    real(dp) :: conductivity

  end type eos_t

contains

  ! Provides a copy subroutine for the eos_t type to
  ! avoid derived type assignment (OpenACC and CUDA can't handle that)
  subroutine copy_eos_t(to_eos, from_eos)

    implicit none

    type(eos_t) :: to_eos, from_eos

    to_eos % rho = from_eos % rho
    to_eos % T = from_eos % T
    to_eos % p = from_eos % p
    to_eos % e = from_eos % e
    to_eos % h = from_eos % h
    to_eos % s = from_eos % s

    to_eos % dpdT = from_eos % dpdT
    to_eos % dpdr = from_eos % dpdr
    to_eos % dedT = from_eos % dedT
    to_eos % dedr = from_eos % dedr
    to_eos % dhdT = from_eos % dhdT
    to_eos % dhdr = from_eos % dhdr
    to_eos % dsdT = from_eos % dsdT
    to_eos % dsdr = from_eos % dsdr
    to_eos % dpde = from_eos % dpde
    to_eos % dpdr_e = from_eos % dpdr_e

    to_eos % cv = from_eos % cv
    to_eos % cp = from_eos % cp
    to_eos % xne = from_eos % xne
    to_eos % xnp = from_eos % xnp
    to_eos % eta = from_eos % eta
    to_eos % pele = from_eos % pele
    to_eos % ppos = from_eos % ppos
    to_eos % mu = from_eos % mu
    to_eos % mu_e = from_eos % mu_e
    to_eos % y_e = from_eos % y_e

    to_eos % gam1 = from_eos % gam1
    to_eos % cs = from_eos % cs

    to_eos % abar = from_eos % abar
    to_eos % zbar = from_eos % zbar

    to_eos % conductivity = from_eos % conductivity

  end subroutine copy_eos_t


  ! Ensure that inputs are within reasonable limits.

  subroutine clean_state(state)

    implicit none

    type (eos_t), intent(inout) :: state

    state % T = min(maxtemp, max(mintemp, state % T))
    state % rho = min(maxdens, max(mindens, state % rho))

  end subroutine clean_state



  ! Print out details of the state.

  subroutine print_state(state)

    implicit none

    type (eos_t), intent(in) :: state

    print *, 'DENS = ', state % rho
    print *, 'TEMP = ', state % T
    print *, 'Y_E  = ', state % y_e
    print *, 'ABAR  = ', state % abar
    print *, 'ZBAR  = ', state % zbar

  end subroutine print_state


  subroutine eos_get_small_temp(small_temp_out)

    !$acc routine seq

    implicit none

    real(dp), intent(out) :: small_temp_out

    small_temp_out = mintemp

  end subroutine eos_get_small_temp



  subroutine eos_get_small_dens(small_dens_out)

    !$acc routine seq

    implicit none

    real(dp), intent(out) :: small_dens_out

    small_dens_out = mindens

  end subroutine eos_get_small_dens



  subroutine eos_get_max_temp(max_temp_out)

    !$acc routine seq

    implicit none

    real(dp), intent(out) :: max_temp_out

    max_temp_out = maxtemp

  end subroutine eos_get_max_temp



  subroutine eos_get_max_dens(max_dens_out)

    !$acc routine seq

    implicit none

    real(dp), intent(out) :: max_dens_out

    max_dens_out = maxdens

  end subroutine eos_get_max_dens


  ! Check to see if variable ivar is a valid
  ! independent variable for the given input
  function eos_input_has_var(input, ivar) result(has)

    implicit none

    integer, intent(in) :: input, ivar
    logical :: has

    has = .false.
    
    select case (ivar)

    case (itemp)

       if (input == eos_input_rt .or. &
           input == eos_input_tp .or. &
           input == eos_input_th) then

          has = .true.

       endif

    case (idens)

       if (input == eos_input_rt .or. &
           input == eos_input_rh .or. &
           input == eos_input_rp .or. &
           input == eos_input_re) then

          has = .true.

       endif

    case (iener)

       if (input == eos_input_re) then

          has = .true.

       endif
       
    case (ienth)

       if (input == eos_input_rh .or. &
           input == eos_input_ph .or. &
           input == eos_input_th) then

          has = .true.

       endif

    case (ientr)

       if (input == eos_input_ps) then

          has = .true.

       endif

    case (ipres)

       if (input == eos_input_tp .or. &
           input == eos_input_rp .or. &
           input == eos_input_ps .or. &
           input == eos_input_ph) then

          has = .true.

       endif

    case default

       call xnet_terminate("EOS: invalid independent variable")

    end select

  end function eos_input_has_var

end module eos_type_module
