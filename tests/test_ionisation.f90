program test_ionisation
  use test_framework
  use ionisation
  implicit none

  call test_atomic_data()
  call test_module_dimensions()

  call print_summary()
  if (exit_code() /= 0) stop 1

contains

!-----------------------------------------------------------------------
subroutine test_atomic_data()
  call set_suite('Ionisation module - atomic data')

  ! Atomic weights (public: a_ion)
  call assert_approx(a_ion(1), 1.00794d0, 1.d-4, 'H atomic weight = 1.00794')
  call assert_approx(a_ion(2), 4.002602d0, 1.d-4, 'He atomic weight = 4.002602')
  call assert_approx(a_ion(3), 12.0107d0, 1.d-3, 'C atomic weight = 12.0107')
  call assert_approx(a_ion(4), 15.9994d0, 1.d-3, 'O atomic weight = 15.9994')
  call assert_approx(a_ion(5), 20.1797d0, 1.d-3, 'Ne atomic weight = 20.1797')
  call assert_approx(a_ion(6), 24.305d0, 1.d-2, 'Mg atomic weight = 24.305')

  ! Atomic numbers (public: iz)
  call assert_equal_int(iz(1), 1, 'H atomic number = 1')
  call assert_equal_int(iz(2), 2, 'He atomic number = 2')
  call assert_equal_int(iz(3), 6, 'C atomic number = 6')
  call assert_equal_int(iz(4), 8, 'O atomic number = 8')
  call assert_equal_int(iz(5), 10, 'Ne atomic number = 10')
  call assert_equal_int(iz(6), 12, 'Mg atomic number = 12')

  ! Number of atoms included
  call assert_equal_int(iatoms, 6, '6 elements included')

  ! Atomic weights must be positive
  call assert_true(all(a_ion > 0.d0), 'all atomic weights positive')

  ! Atomic numbers must be positive and increasing
  call assert_true(iz(1) < iz(2), 'iz(H) < iz(He)')
  call assert_true(iz(2) < iz(3), 'iz(He) < iz(C)')
  call assert_true(iz(3) < iz(4), 'iz(C) < iz(O)')
  call assert_true(iz(4) < iz(5), 'iz(O) < iz(Ne)')
  call assert_true(iz(5) < iz(6), 'iz(Ne) < iz(Mg)')

end subroutine test_atomic_data

!-----------------------------------------------------------------------
subroutine test_module_dimensions()
  call set_suite('Ionisation module - dimensions and state')

  ! Module state variables should be initialisable
  abond = 0.d0
  abond(1) = 0.70d0
  abond(2) = 0.28d0
  call assert_approx(abond(1), 0.70d0, 1.d-15, 'abond(H) settable to 0.70')
  call assert_approx(abond(2), 0.28d0, 1.d-15, 'abond(He) settable to 0.28')

  ! vnu should be settable
  vnu = 0.d0
  vnu(1) = abond(1) / a_ion(1)
  vnu(2) = abond(2) / a_ion(2)
  call assert_true(vnu(1) > 0.d0, 'vnu(H) > 0 after setting')
  call assert_true(vnu(2) > 0.d0, 'vnu(He) > 0 after setting')

  ! H should have more atoms per gram than He
  call assert_true(vnu(1) > vnu(2), 'vnu(H) > vnu(He) for solar composition')

  ! list should be settable
  list = 0
  list(1) = 1
  list(2) = 2
  call assert_equal_int(list(1), 1, 'list(1) = H')
  call assert_equal_int(list(2), 2, 'list(2) = He')

  ! ionized flag
  ionized = 0
  call assert_equal_int(ionized, 0, 'ionized = 0 (not fully ionised)')
  ionized = 1
  call assert_equal_int(ionized, 1, 'ionized = 1 (fully ionised)')

end subroutine test_module_dimensions

end program test_ionisation
