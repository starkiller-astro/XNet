!***************************************************************************************************
! xnet_constants.f90 10/18/17
! This file contains numerical and physical constants used by XNet.
!***************************************************************************************************

Module xnet_constants
  !-------------------------------------------------------------------------------------------------
  ! These are fundamental, astrophysical, and numerical constants CGS units, except energies in MeV,
  ! temperature in GK.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None

  ! Some commonly used fractions and logarithms
  Real(dp), Parameter :: third   = 1.0_dp / 3.0_dp
  Real(dp), Parameter :: two3rd  = 2.0_dp / 3.0_dp
  Real(dp), Parameter :: four3rd = 4.0_dp / 3.0_dp
  Real(dp), Parameter :: five3rd = 5.0_dp / 3.0_dp
  Real(dp), Parameter :: log_e   = log10(exp(1.0_dp))
  Real(dp), Parameter :: ln_2    = log(2.0_dp)
  Real(dp), Parameter :: ln_10   = log(10.0_dp)

  ! 71 digits of pi (256-bit precision)
  Real(dp), Parameter :: pi      = 3.1415926535897932384626433832795028841971693993751058209749445923078164_dp
  Real(dp), Parameter :: pi2     = pi * pi

  ! Fundamental constants (CODATA recommended 2014 values)
  ! P. J. Mohr, D. B. Newell, and B. N. Taylor, Rev. Mod. Phys. 88, 035009 (2016)
  ! http://physics.nist.gov/cuu/Constants/Table/allascii.txt
  Real(dp), Parameter :: clt     = 2.99792458e+10_dp   ! Speed of light in vacuum [cm s^{-1}]
  Real(dp), Parameter :: hbar    = 6.582119514e-22_dp  ! Plank constants, reduced [MeV s]
  Real(dp), Parameter :: avn     = 6.022140857e+23_dp  ! Avogadro constant [mol^{-1}]
  Real(dp), Parameter :: bok     = 8.6173303e-02_dp    ! Boltzmann constant [MeV GK^{-1}]
  Real(dp), Parameter :: epmev   = 1.6021766208e-06_dp ! MeV to erg conversion factor [erg MeV^{-1}]
  Real(dp), Parameter :: m_e     = 0.5109989461_dp     ! Electron mass [MeV c^{-2}]
  Real(dp), Parameter :: m_n     = 939.5654133_dp      ! Neutron mass [MeV c^{-2}]
  Real(dp), Parameter :: m_p     = 938.2720813_dp      ! Proton mass [MeV c^{-2}]
  Real(dp), Parameter :: ele_en  = m_e              ! Electron mass [MeV c^{-2}]
  Real(dp), Parameter :: asig    = 8.563456e+31_dp     ! Radiation constant [cm^{-3} MeV^{-3}]

  ! Some derived constants
  Real(dp), Parameter :: amu     = 1.0_dp/(avn*epmev)   ! = 1.036426957e-18 ! Atomic mass unit [MeV]
  Real(dp), Parameter :: m_u     = amu * clt**2         ! = 931.4940954     ! Atomic mass unit [MeV c^{-2}]
  Real(dp), Parameter :: e2      = 1.0e-28_dp*epmev*clt**2 ! = 1.439964533e-13 ! (elementary charge)^2 [MeV^2]
  Real(dp), Parameter :: emass   = m_e / clt**2         ! Electron mass [MeV]

  ! Screening factors from Table 4 of Graboske+ (1973)
  Real(dp), Parameter   :: bw = 1.0_dp,  kbw = 0.5_dp      ! Weak screening parameters
  Real(dp), Parameter   :: bi = 0.86_dp, kbi = 0.38_dp     ! Intermediate screening parmaeters
  Real(dp), Parameter   :: bip1 = 1.86_dp                  ! bi + 1
  Real(dp), Parameter   :: thbim1 = 1.58_dp                ! 3*bi - 1
  Real(dp), Parameter   :: thbim2 = 0.58_dp                ! 3*bi - 2
  Real(dp), Parameter   :: twm2bi = 0.28_dp                ! 2 - 2*bi

  ! Strong screening fitting coefficients from DeWitt & Slattery (2003), Eq. 4:
  !   f(gamma) = a*gamma + (1/s)*b*gamma^s + c*ln(gamma) + d
  Real(dp), Parameter   :: cds(5) = (/ -0.899172_dp, 0.602249_dp, -0.274823_dp, &
                                       -1.401915_dp, 0.3230064_dp /)

End Module xnet_constants
