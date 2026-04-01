module test_framework
  implicit none

  integer, save :: tests_run = 0
  integer, save :: tests_passed = 0
  integer, save :: tests_failed = 0
  character(len=256), save :: current_suite = ''

contains

  subroutine set_suite(name)
    character(len=*), intent(in) :: name
    current_suite = name
    write(*,'(a)') ''
    write(*,'(a,a)') '=== ', trim(name)
  end subroutine set_suite

  subroutine assert_true(condition, test_name)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: test_name
    tests_run = tests_run + 1
    if (condition) then
      tests_passed = tests_passed + 1
      write(*,'(a,a)') '  PASS: ', trim(test_name)
    else
      tests_failed = tests_failed + 1
      write(*,'(a,a)') '  FAIL: ', trim(test_name)
    end if
  end subroutine assert_true

  subroutine assert_equal_int(actual, expected, test_name)
    integer, intent(in) :: actual, expected
    character(len=*), intent(in) :: test_name
    tests_run = tests_run + 1
    if (actual == expected) then
      tests_passed = tests_passed + 1
      write(*,'(a,a)') '  PASS: ', trim(test_name)
    else
      tests_failed = tests_failed + 1
      write(*,'(a,a,a,i0,a,i0)') '  FAIL: ', trim(test_name), ' (got ', actual, ', expected ', expected
    end if
  end subroutine assert_equal_int

  subroutine assert_approx(actual, expected, tol, test_name)
    real(8), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: test_name
    tests_run = tests_run + 1
    if (abs(actual - expected) <= tol) then
      tests_passed = tests_passed + 1
      write(*,'(a,a)') '  PASS: ', trim(test_name)
    else
      tests_failed = tests_failed + 1
      write(*,'(a,a)') '  FAIL: ', trim(test_name)
      write(*,'(a,es20.12,a,es20.12,a,es10.2)') '        got=', actual, ' expected=', expected, ' tol=', tol
    end if
  end subroutine assert_approx

  subroutine assert_rel(actual, expected, reltol, test_name)
    real(8), intent(in) :: actual, expected, reltol
    character(len=*), intent(in) :: test_name
    real(8) :: denom
    tests_run = tests_run + 1
    denom = max(abs(expected), 1.0d-30)
    if (abs(actual - expected) / denom <= reltol) then
      tests_passed = tests_passed + 1
      write(*,'(a,a)') '  PASS: ', trim(test_name)
    else
      tests_failed = tests_failed + 1
      write(*,'(a,a)') '  FAIL: ', trim(test_name)
      write(*,'(a,es20.12,a,es20.12,a,es10.2)') '        got=', actual, ' expected=', expected, ' reltol=', reltol
    end if
  end subroutine assert_rel

  subroutine print_summary()
    write(*,'(a)') ''
    write(*,'(a)') '==============================='
    write(*,'(a,i0,a,i0,a)') 'Results: ', tests_passed, '/', tests_run, ' passed'
    if (tests_failed > 0) then
      write(*,'(a,i0,a)') 'FAILURES: ', tests_failed, ' test(s) failed'
    else
      write(*,'(a)') 'All tests passed.'
    end if
    write(*,'(a)') '==============================='
  end subroutine print_summary

  integer function exit_code()
    if (tests_failed > 0) then
      exit_code = 1
    else
      exit_code = 0
    end if
  end function exit_code

end module test_framework
