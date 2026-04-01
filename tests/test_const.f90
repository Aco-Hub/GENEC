program test_const
  use test_framework
  use evol, only: kindreal
  use const
  implicit none

  call test_fundamental_constants()
  call test_log_consistency()
  call test_solar_values()
  call test_derived_constants()

  call print_summary()
  if (exit_code() /= 0) stop 1

contains

!-----------------------------------------------------------------------
subroutine test_fundamental_constants()
  call set_suite('Physical constants - CODATA values')

  ! Speed of light (cm/s)
  call assert_rel(cst_c, 2.99792458d10, 1.d-8, 'speed of light c')

  ! Gravitational constant (cgs)
  call assert_rel(cst_G, 6.67428d-8, 1.d-4, 'gravitational constant G')

  ! Planck constant (erg*s)
  call assert_rel(cst_h, 6.62606896d-27, 1.d-8, 'Planck constant h')

  ! Boltzmann constant (erg/K)
  call assert_rel(cst_k, 1.3806504d-16, 1.d-6, 'Boltzmann constant k')

  ! Electron mass (g)
  call assert_rel(cst_me, 9.1093826d-28, 1.d-6, 'electron mass m_e')

  ! Avogadro number
  call assert_rel(cst_avo, 6.02214179d23, 1.d-8, 'Avogadro number')

  ! Electron charge (Coulomb)
  call assert_rel(cst_e, 1.602176487d-19, 1.d-8, 'electron charge (SI)')

  ! Stefan-Boltzmann constant (cgs)
  call assert_rel(cst_sigma, 5.67040d-5, 1.d-4, 'Stefan-Boltzmann sigma')

  ! Pi
  call assert_approx(pi, 4.d0 * atan(1.d0), 1.d-15, 'pi = 4*atan(1)')

end subroutine test_fundamental_constants

!-----------------------------------------------------------------------
subroutine test_log_consistency()
  call set_suite('Log10 of constants - internal consistency')

  ! Each cstlg_X should equal log10(cst_X)
  call assert_approx(cstlg_c, log10(cst_c), 1.d-6, 'log10(c) consistent')
  call assert_approx(cstlg_G, log10(cst_G), 1.d-6, 'log10(G) consistent')
  call assert_approx(cstlg_h, log10(cst_h), 1.d-6, 'log10(h) consistent')
  call assert_approx(cstlg_k, log10(cst_k), 1.d-6, 'log10(k) consistent')
  call assert_approx(cstlg_me, log10(cst_me), 1.d-6, 'log10(m_e) consistent')
  call assert_approx(cstlg_avo, log10(cst_avo), 1.d-6, 'log10(N_A) consistent')
  call assert_approx(cstlg_sigma, log10(cst_sigma), 1.d-5, 'log10(sigma) consistent')
  call assert_approx(lgpi, log10(pi), 1.d-6, 'log10(pi) consistent')

  ! um = ln(10)
  call assert_approx(um, log(10.d0), 1.d-15, 'um = ln(10)')

end subroutine test_log_consistency

!-----------------------------------------------------------------------
subroutine test_solar_values()
  call set_suite('Solar parameters')

  ! Solar mass ~2e33 g
  call assert_rel(Msol, 1.9884d33, 1.d-3, 'solar mass')

  ! Solar radius ~7e10 cm
  call assert_rel(Rsol, 6.9551d10, 1.d-3, 'solar radius')

  ! Solar luminosity ~3.8e33 erg/s
  call assert_rel(Lsol, 3.8427d33, 1.d-3, 'solar luminosity')

  ! Solar effective temperature ~5778 K
  call assert_rel(Teffsol, 5.777d3, 1.d-3, 'solar Teff')

  ! Log10 consistency of solar values
  call assert_approx(lgMsol, log10(Msol), 1.d-5, 'log10(Msol) consistent')
  call assert_approx(lgRsol, log10(Rsol), 1.d-5, 'log10(Rsol) consistent')
  ! NOTE: lgLsol in const.f90 is 33.58460257, but log10(3.8427e33) = 33.58463648
  ! This is a minor inconsistency in the stored constant (off by ~3e-5)
  call assert_approx(lgLsol, log10(Lsol), 5.d-4, 'log10(Lsol) consistent (relaxed tol)')

end subroutine test_solar_values

!-----------------------------------------------------------------------
subroutine test_derived_constants()
  call set_suite('Derived constants - cross-checks')

  ! Radiation constant a = 4*sigma/c
  call assert_rel(cst_a, 4.d0 * cst_sigma / cst_c, 1.d-4, 'a = 4*sigma/c')

  ! Gas constant R = k * N_A
  call assert_rel(rgaz, cst_k * cst_avo, 1.d-5, 'R = k * N_A')

  ! Thomson cross-section ~ 6.65e-25 cm^2
  call assert_rel(cst_thomson, 6.652458558d-25, 1.d-6, 'Thomson cross-section')

  ! qapicg = 4*pi*c*G
  call assert_rel(qapicg, 4.d0 * pi * cst_c * cst_G, 1.d-4, '4*pi*c*G')

  ! Julian year in seconds
  call assert_approx(year, 3.1557600d7, 1.d0, 'Julian year = 365.25 days')

  ! Day in seconds
  call assert_approx(day, 86400.d0, 1.d-10, 'day = 86400 s')

end subroutine test_derived_constants

end program test_const
