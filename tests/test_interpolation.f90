program test_interpolation
  use test_framework
  use interpolation
  implicit none

  call test_indic()
  call test_fipoi()
  call test_fipoi1()
  call test_fipoi2()
  call test_quint()
  call test_qua()
  call test_quad_gg()
  call test_flin()
  call test_spline_splint()
  call test_getd()

  call print_summary()
  if (exit_code() /= 0) stop 1

contains

!-----------------------------------------------------------------------
subroutine test_indic()
  integer :: k
  real(8) :: x_asc(5), x_desc(5)

  call set_suite('indic - binary search in monotonic table')

  x_asc = (/1.d0, 2.d0, 3.d0, 4.d0, 5.d0/)
  x_desc = (/5.d0, 4.d0, 3.d0, 2.d0, 1.d0/)

  ! Value in the middle of ascending table
  k = indic(2.5d0, x_asc, 5)
  call assert_equal_int(k, 2, 'ascending table, value between x(2) and x(3)')

  ! Value at a table point
  k = indic(3.0d0, x_asc, 5)
  call assert_equal_int(k, 3, 'ascending table, exact match at x(3)')

  ! Value below table range
  k = indic(0.5d0, x_asc, 5)
  call assert_equal_int(k, 1, 'ascending table, below range returns 1')

  ! Value above table range
  k = indic(5.5d0, x_asc, 5)
  call assert_equal_int(k, 4, 'ascending table, above range returns m-1')

  ! Descending table
  k = indic(3.5d0, x_desc, 5)
  call assert_equal_int(k, 2, 'descending table, value between x(2) and x(3)')

  ! Descending below range
  k = indic(0.5d0, x_desc, 5)
  call assert_equal_int(k, 4, 'descending table, below range returns m-1')

end subroutine test_indic

!-----------------------------------------------------------------------
subroutine test_fipoi()
  real(8) :: a(4), b(4), result

  call set_suite('fipoi - linear interpolation')

  a = (/1.d0, 2.d0, 3.d0, 4.d0/)
  b = (/10.d0, 20.d0, 30.d0, 40.d0/)

  ! Midpoint interpolation
  result = fipoi(2.5d0, 4, a, b)
  call assert_approx(result, 25.d0, 1.d-12, 'midpoint gives 25')

  ! At a table point
  result = fipoi(3.0d0, 4, a, b)
  call assert_approx(result, 30.d0, 1.d-12, 'exact point gives 30')

  ! Quarter point
  result = fipoi(1.25d0, 4, a, b)
  call assert_approx(result, 12.5d0, 1.d-12, 'quarter point gives 12.5')

  ! Non-uniform y-values: y = x^2
  b = (/1.d0, 4.d0, 9.d0, 16.d0/)
  result = fipoi(1.5d0, 4, a, b)
  call assert_approx(result, 2.5d0, 1.d-12, 'linear interp of x^2 at 1.5 gives 2.5')

end subroutine test_fipoi

!-----------------------------------------------------------------------
subroutine test_fipoi1()
  real(8) :: a(3), b(3), r1

  call set_suite('fipoi1 - linear interpolation subroutine')

  a = (/0.d0, 1.d0, 2.d0/)
  b = (/0.d0, 5.d0, 10.d0/)

  call fipoi1(0.5d0, 3, a, b, r1)
  call assert_approx(r1, 2.5d0, 1.d-12, 'midpoint of [0,5] gives 2.5')

  ! Consistency with fipoi
  call assert_approx(r1, fipoi(0.5d0, 3, a, b), 1.d-15, 'fipoi1 equals fipoi')

end subroutine test_fipoi1

!-----------------------------------------------------------------------
subroutine test_fipoi2()
  real(8) :: a(4), b(4), r1, r2

  call set_suite('fipoi2 - linear interpolation with derivative')

  a = (/1.d0, 2.d0, 3.d0, 4.d0/)
  b = (/10.d0, 20.d0, 30.d0, 40.d0/)

  call fipoi2(2.5d0, 4, a, b, r1, r2)
  call assert_approx(r1, 25.d0, 1.d-12, 'value at 2.5 is 25')
  call assert_approx(r2, 10.d0, 1.d-12, 'derivative is 10 (constant slope)')

  ! Non-uniform spacing
  a = (/0.d0, 1.d0, 3.d0, 6.d0/)
  b = (/0.d0, 2.d0, 8.d0, 20.d0/)

  call fipoi2(2.0d0, 4, a, b, r1, r2)
  ! Between a(2)=1 and a(3)=3, slope = (8-2)/(3-1) = 3
  call assert_approx(r1, 5.d0, 1.d-12, 'non-uniform: value at 2.0')
  call assert_approx(r2, 3.d0, 1.d-12, 'non-uniform: derivative at 2.0')

end subroutine test_fipoi2

!-----------------------------------------------------------------------
subroutine test_quint()
  real(8) :: y

  call set_suite('quint - quadratic interpolation (equidistant)')

  ! Parabola y = x^2, points at x=0,1,2 with h=1
  call quint(0.5d0, 0.d0, 1.d0, 0.d0, 1.d0, 4.d0, y)
  call assert_approx(y, 0.25d0, 1.d-12, 'x^2 at 0.5 gives 0.25')

  call quint(1.5d0, 0.d0, 1.d0, 0.d0, 1.d0, 4.d0, y)
  call assert_approx(y, 2.25d0, 1.d-12, 'x^2 at 1.5 gives 2.25')

  ! Linear function y = 3x+1, points at x=0,1,2
  call quint(0.7d0, 0.d0, 1.d0, 1.d0, 4.d0, 7.d0, y)
  call assert_approx(y, 3.1d0, 1.d-12, '3x+1 at 0.7 gives 3.1')

  ! Passes through the data points
  call quint(0.d0, 0.d0, 1.d0, 0.d0, 1.d0, 4.d0, y)
  call assert_approx(y, 0.d0, 1.d-12, 'passes through (0, 0)')

  call quint(1.d0, 0.d0, 1.d0, 0.d0, 1.d0, 4.d0, y)
  call assert_approx(y, 1.d0, 1.d-12, 'passes through (1, 1)')

  call quint(2.d0, 0.d0, 1.d0, 0.d0, 1.d0, 4.d0, y)
  call assert_approx(y, 4.d0, 1.d-12, 'passes through (2, 4)')

end subroutine test_quint

!-----------------------------------------------------------------------
subroutine test_qua()
  real(8) :: result

  call set_suite('qua - Lagrange quadratic interpolation')

  ! Parabola y = x^2, through (1,1), (2,4), (3,9)
  result = qua(1.d0, 2.d0, 3.d0, 1.d0, 4.d0, 9.d0, 1.5d0)
  call assert_approx(result, 2.25d0, 1.d-12, 'x^2 at 1.5 gives 2.25')

  result = qua(1.d0, 2.d0, 3.d0, 1.d0, 4.d0, 9.d0, 2.5d0)
  call assert_approx(result, 6.25d0, 1.d-12, 'x^2 at 2.5 gives 6.25')

  ! Passes through the three data points
  result = qua(1.d0, 2.d0, 3.d0, 1.d0, 4.d0, 9.d0, 1.d0)
  call assert_approx(result, 1.d0, 1.d-12, 'passes through (1,1)')

  result = qua(1.d0, 2.d0, 3.d0, 1.d0, 4.d0, 9.d0, 2.d0)
  call assert_approx(result, 4.d0, 1.d-12, 'passes through (2,4)')

  result = qua(1.d0, 2.d0, 3.d0, 1.d0, 4.d0, 9.d0, 3.d0)
  call assert_approx(result, 9.d0, 1.d-12, 'passes through (3,9)')

  ! Linear function y = 2x+1, through (0,1), (1,3), (2,5)
  result = qua(0.d0, 1.d0, 2.d0, 1.d0, 3.d0, 5.d0, 0.75d0)
  call assert_approx(result, 2.5d0, 1.d-12, 'linear 2x+1 at 0.75 gives 2.5')

end subroutine test_qua

!-----------------------------------------------------------------------
subroutine test_quad_gg()
  real(8) :: deriv

  call set_suite('quad_gg - Lagrange quadratic derivative')

  ! y = x^2 through (1,1), (2,4), (3,9) => dy/dx = 2x
  deriv = quad_gg(1.d0, 2.d0, 3.d0, 1.d0, 4.d0, 9.d0, 2.d0)
  call assert_approx(deriv, 4.d0, 1.d-12, 'dy/dx of x^2 at x=2 is 4')

  deriv = quad_gg(1.d0, 2.d0, 3.d0, 1.d0, 4.d0, 9.d0, 1.5d0)
  call assert_approx(deriv, 3.d0, 1.d-12, 'dy/dx of x^2 at x=1.5 is 3')

  ! Linear function y = 5x => dy/dx = 5 everywhere
  deriv = quad_gg(0.d0, 1.d0, 2.d0, 0.d0, 5.d0, 10.d0, 0.5d0)
  call assert_approx(deriv, 5.d0, 1.d-12, 'dy/dx of 5x is 5')

  ! Consistency: numerical derivative of qua should match quad_gg
  block
    real(8) :: h, f_plus, f_minus, num_deriv
    h = 1.d-6
    f_plus  = qua(1.d0, 2.d0, 3.d0, 1.d0, 4.d0, 9.d0, 2.d0 + h)
    f_minus = qua(1.d0, 2.d0, 3.d0, 1.d0, 4.d0, 9.d0, 2.d0 - h)
    num_deriv = (f_plus - f_minus) / (2.d0 * h)
    deriv = quad_gg(1.d0, 2.d0, 3.d0, 1.d0, 4.d0, 9.d0, 2.d0)
    call assert_approx(deriv, num_deriv, 1.d-6, 'quad_gg matches numerical derivative of qua')
  end block

end subroutine test_quad_gg

!-----------------------------------------------------------------------
subroutine test_flin()
  real(8) :: result

  call set_suite('flin - linear interpolation/extrapolation')

  ! Interpolation: (1,10) to (3,30), at x=2
  result = flin(1.d0, 3.d0, 10.d0, 30.d0, 2.d0)
  call assert_approx(result, 20.d0, 1.d-12, 'interpolation at midpoint')

  ! At endpoint
  result = flin(1.d0, 3.d0, 10.d0, 30.d0, 1.d0)
  call assert_approx(result, 10.d0, 1.d-12, 'value at x1')

  result = flin(1.d0, 3.d0, 10.d0, 30.d0, 3.d0)
  call assert_approx(result, 30.d0, 1.d-12, 'value at x2')

  ! Extrapolation beyond x2
  result = flin(1.d0, 3.d0, 10.d0, 30.d0, 5.d0)
  call assert_approx(result, 50.d0, 1.d-12, 'extrapolation beyond x2')

  ! Extrapolation below x1
  result = flin(1.d0, 3.d0, 10.d0, 30.d0, -1.d0)
  call assert_approx(result, -10.d0, 1.d-12, 'extrapolation below x1')

end subroutine test_flin

!-----------------------------------------------------------------------
subroutine test_spline_splint()
  integer, parameter :: n = 5
  real(8) :: x(n), y(n), y2(n)
  real(8) :: yval, ypval

  call set_suite('spline/splint - cubic spline interpolation')

  ! Test with y = x^3 on [0, 4]
  x = (/0.d0, 1.d0, 2.d0, 3.d0, 4.d0/)
  y = (/0.d0, 1.d0, 8.d0, 27.d0, 64.d0/)

  call spline(x, y, n, y2)

  ! Spline should pass through data points
  call splint(x, y, n, y2, 1.d0, yval, ypval)
  call assert_approx(yval, 1.d0, 1.d-6, 'spline passes through (1, 1)')

  call splint(x, y, n, y2, 2.d0, yval, ypval)
  call assert_approx(yval, 8.d0, 1.d-6, 'spline passes through (2, 8)')

  call splint(x, y, n, y2, 3.d0, yval, ypval)
  call assert_approx(yval, 27.d0, 1.d-6, 'spline passes through (3, 27)')

  ! Interpolated value should be reasonable
  call splint(x, y, n, y2, 1.5d0, yval, ypval)
  ! x^3 at 1.5 = 3.375, cubic spline should be close
  call assert_approx(yval, 3.375d0, 0.5d0, 'spline at 1.5 close to 3.375')

  ! Second derivatives should be computed (non-zero for cubic)
  call assert_true(any(y2 /= 0.d0), 'y2 array is non-zero for cubic data')

end subroutine test_spline_splint

!-----------------------------------------------------------------------
subroutine test_getd()
  integer, parameter :: n = 10
  real(8) :: f(n), d(n), fp1, fpn
  integer :: i

  call set_suite('getd - spline coefficients (unit intervals)')

  ! Linear function f(i) = 2*i + 1 (on unit intervals)
  ! Second derivative should be ~0 for linear data
  do i = 1, n
    f(i) = 2.d0 * dble(i) + 1.d0
  end do

  call getd(f, n, d, fp1, fpn)

  ! For linear data, first derivative at endpoints should be 2
  call assert_approx(fp1, 2.d0, 1.d-10, 'endpoint derivative fp1 = 2 for linear')
  call assert_approx(fpn, 2.d0, 1.d-10, 'endpoint derivative fpn = 2 for linear')

  ! Second derivatives should be ~0 for linear data
  do i = 1, n
    call assert_approx(d(i), 0.d0, 1.d-10, 'second deriv ~0 for linear data')
  end do

end subroutine test_getd

end program test_interpolation
