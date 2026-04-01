program test_smallfunc
  use test_framework
  use SmallFunc
  implicit none

  call test_neg_root()
  call test_expf10()
  call test_exphi()
  call test_primus_drimus()
  call test_secun_dsecun()
  call test_tertiu_dterti()
  call test_pos()
  call test_flin_valinterp()
  call test_girl()

  call print_summary()
  if (exit_code() /= 0) stop 1

contains

!-----------------------------------------------------------------------
subroutine test_neg_root()
  call set_suite('neg_root - signed power function')

  ! Positive base
  call assert_approx(neg_root(4.d0, 0.5d0), 2.d0, 1.d-12, '4^0.5 = 2')
  call assert_approx(neg_root(8.d0, 1.d0/3.d0), 2.d0, 1.d-12, '8^(1/3) = 2')
  call assert_approx(neg_root(1.d0, 5.d0), 1.d0, 1.d-12, '1^5 = 1')

  ! Negative base: sign preserved
  call assert_approx(neg_root(-8.d0, 1.d0/3.d0), -2.d0, 1.d-12, '-8^(1/3) = -2')
  call assert_approx(neg_root(-4.d0, 0.5d0), -2.d0, 1.d-12, '-4^0.5 = -2 (sign preserved)')
  call assert_approx(neg_root(-1.d0, 2.d0), -1.d0, 1.d-12, '-1^2 = -1 (sign preserved)')

  ! Zero
  call assert_approx(neg_root(0.d0, 2.d0), 0.d0, 1.d-12, '0^2 = 0')

end subroutine test_neg_root

!-----------------------------------------------------------------------
subroutine test_expf10()
  call set_suite('expf10 - 10^x function')

  call assert_approx(expf10(0.d0), 1.d0, 1.d-12, '10^0 = 1')
  call assert_approx(expf10(1.d0), 10.d0, 1.d-10, '10^1 = 10')
  call assert_approx(expf10(2.d0), 100.d0, 1.d-8, '10^2 = 100')
  call assert_approx(expf10(-1.d0), 0.1d0, 1.d-12, '10^-1 = 0.1')
  call assert_approx(expf10(0.5d0), sqrt(10.d0), 1.d-10, '10^0.5 = sqrt(10)')

end subroutine test_expf10

!-----------------------------------------------------------------------
subroutine test_exphi()
  real(8) :: x, exact

  call set_suite('exphi - numerically stable 1-exp(x)')

  ! Large |x|: uses standard formula
  call assert_approx(exphi(1.d0), 1.d0 - exp(1.d0), 1.d-12, 'x=1: 1-exp(1)')
  call assert_approx(exphi(-1.d0), 1.d0 - exp(-1.d0), 1.d-12, 'x=-1: 1-exp(-1)')

  ! Small |x|: uses Taylor expansion for stability
  x = 1.d-4
  exact = 1.d0 - exp(x)
  call assert_approx(exphi(x), exact, 1.d-10, 'x=1e-4: Taylor series accurate')

  x = -1.d-4
  exact = 1.d0 - exp(x)
  call assert_approx(exphi(x), exact, 1.d-10, 'x=-1e-4: Taylor series accurate')

  ! Zero
  call assert_approx(exphi(0.d0), 0.d0, 1.d-15, 'x=0: exphi = 0')

end subroutine test_exphi

!-----------------------------------------------------------------------
subroutine test_primus_drimus()
  real(8) :: f, df, t9

  call set_suite('primus/drimus - reaction rate and derivative')

  t9 = 1.d0

  ! primus(a, b, e1, e2, t9) = a / t9^e1 * exp(-b / t9^e2)
  ! With a=1, b=0, e1=0, e2=1: result = 1
  f = primus(1.d0, 0.d0, 0.d0, 1.d0, t9)
  call assert_approx(f, 1.d0, 1.d-12, 'primus with b=0 gives a=1')

  ! With a=2, b=1, e1=0, e2=1, t9=1: result = 2*exp(-1)
  f = primus(2.d0, 1.d0, 0.d0, 1.d0, 1.d0)
  call assert_approx(f, 2.d0 * exp(-1.d0), 1.d-12, 'primus = 2*exp(-1)')

  ! Underflow protection: very large negative exponent
  f = primus(1.d0, 1.d10, 0.d0, 1.d0, 1.d0)
  call assert_approx(f, 0.d0, 1.d-30, 'primus returns 0 for extreme underflow')

  ! Derivative consistency: drimus = d(log(primus))/d(log(t9))
  ! drimus(b, e1, e2, t9) = -e1 + e2*b/t9^e2
  df = drimus(1.d0, 2.d0, 1.d0, 1.d0)
  ! = -2 + 1*1/1^1 = -1
  call assert_approx(df, -1.d0, 1.d-12, 'drimus(-e1 + e2*b/t9^e2)')

end subroutine test_primus_drimus

!-----------------------------------------------------------------------
subroutine test_secun_dsecun()
  real(8) :: f

  call set_suite('secun/dsecun - extended reaction rate')

  ! secun(a, b, c, e1, e2, e3, t9) = a / t9^e1 * exp(-b/t9^e2 - (t9/c)^e3)
  ! With b=0, c=1e30 (negligible second term): reduces to a/t9^e1
  f = secun(5.d0, 0.d0, 1.d30, 2.d0, 1.d0, 1.d0, 2.d0)
  call assert_approx(f, 5.d0 / 4.d0, 1.d-8, 'secun reduces to a/t9^e1 when b=0, c>>1')

  ! Full calculation: a=1, b=1, c=1, e1=0, e2=1, e3=1, t9=1
  ! = 1 * exp(-1 - 1) = exp(-2)
  f = secun(1.d0, 1.d0, 1.d0, 0.d0, 1.d0, 1.d0, 1.d0)
  call assert_approx(f, exp(-2.d0), 1.d-12, 'secun = exp(-2)')

  ! dsecun(b, c, e1, e2, e3, t9) = -e1 + e2*b/t9^e2 - e3*(t9/c)^e3
  ! b=1, c=2, e1=1, e2=1, e3=1, t9=1 => -1 + 1 - 0.5 = -0.5
  call assert_approx(dsecun(1.d0, 2.d0, 1.d0, 1.d0, 1.d0, 1.d0), -0.5d0, 1.d-12, &
       'dsecun = -e1 + e2*b/t9^e2 - e3*(t9/c)^e3')

end subroutine test_secun_dsecun

!-----------------------------------------------------------------------
subroutine test_tertiu_dterti()
  real(8) :: f, df

  call set_suite('tertiu/dterti - polynomial reaction rate')

  ! tertiu = a + b*t9^e1 + c*t9^e2 + d*t9^e3 + e*t9^e4 + f*t9^e5
  ! Simple case: all zero except a
  f = tertiu(7.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 1.d0)
  call assert_approx(f, 7.d0, 1.d-12, 'constant term only')

  ! Linear: a=0, b=3, e1=1, rest zero, t9=2 => 0 + 3*2 = 6
  f = tertiu(0.d0, 3.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 2.d0)
  call assert_approx(f, 6.d0, 1.d-12, 'linear: 3*t9')

  ! dterti = e1*b*t9^e1 + e2*c*t9^e2 + ...
  ! For b=3, e1=1, rest zero, t9=2 => 1*3*2 = 6
  df = dterti(3.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 2.d0)
  call assert_approx(df, 6.d0, 1.d-12, 'dterti linear: e1*b*t9^e1 = 6')

  ! Quadratic: a=1, b=2, c=3, e1=1, e2=2, t9=2
  ! f = 1 + 2*2 + 3*4 = 1 + 4 + 12 = 17
  f = tertiu(1.d0, 2.d0, 3.d0, 0.d0, 0.d0, 0.d0, 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 2.d0)
  call assert_approx(f, 17.d0, 1.d-12, 'quadratic: 1 + 2*t9 + 3*t9^2')

  ! df = 1*2*2 + 2*3*4 = 4 + 24 = 28
  df = dterti(2.d0, 3.d0, 0.d0, 0.d0, 0.d0, 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 2.d0)
  call assert_approx(df, 28.d0, 1.d-12, 'dterti quadratic derivative')

end subroutine test_tertiu_dterti

!-----------------------------------------------------------------------
subroutine test_pos()
  integer :: k
  real(8) :: x(5)

  call set_suite('pos - binary search (same as indic)')

  x = (/1.d0, 3.d0, 5.d0, 7.d0, 9.d0/)

  k = pos(4.d0, x, 5)
  call assert_equal_int(k, 2, 'value 4 between x(2)=3 and x(3)=5')

  k = pos(6.d0, x, 5)
  call assert_equal_int(k, 3, 'value 6 between x(3)=5 and x(4)=7')

  k = pos(1.d0, x, 5)
  call assert_equal_int(k, 1, 'value at x(1)=1')

  k = pos(10.d0, x, 5)
  call assert_equal_int(k, 4, 'value beyond table returns m-1')

end subroutine test_pos

!-----------------------------------------------------------------------
subroutine test_flin_valinterp()
  real(8) :: result

  call set_suite('ValInterp/OmInterp - 2D interpolation')

  ! ValInterp(x1, x2, y1, y2, y3) = (x1*(y2-y3) + x2*(y3-y1)) / (y2-y1)
  ! Simple case: x1=0, x2=10, y1=0, y2=1, y3=0.5 => (0*0.5 + 10*0.5)/1 = 5
  result = ValInterp(0.d0, 10.d0, 0.d0, 1.d0, 0.5d0)
  call assert_approx(result, 5.d0, 1.d-12, 'ValInterp midpoint')

  ! At y3=y1: result = x1
  result = ValInterp(3.d0, 7.d0, 2.d0, 5.d0, 2.d0)
  call assert_approx(result, 3.d0, 1.d-12, 'ValInterp at y3=y1 gives x1')

  ! At y3=y2: result = x2
  result = ValInterp(3.d0, 7.d0, 2.d0, 5.d0, 5.d0)
  call assert_approx(result, 7.d0, 1.d-12, 'ValInterp at y3=y2 gives x2')

end subroutine test_flin_valinterp

!-----------------------------------------------------------------------
subroutine test_girl()
  ! girl inverts a matrix using Gauss elimination
  ! Input: a is a 1D array of size n*(n+m) representing [A | I]
  ! Output: b is a 1D array of size n*m representing A^-1
  real(8) :: prod

  call set_suite('girl - matrix inversion (Gauss elimination)')

  ! 2x2 identity matrix: A = I, A^-1 = I
  ! Layout: column-major, a = [A | I] stored as n*(n+m) = 2*4 = 8... wait
  ! Actually n=2, m=2, so a has size n*(n+m) = 2*4 = 8
  ! But the columns are interleaved: a(1:n) is column 1, a(n+1:2n) is column 2, etc.

  ! 2x2 matrix A = [[2, 1], [1, 3]], augmented with I
  ! Column-major: col1=[2,1], col2=[1,3], col3=[1,0], col4=[0,1]
  ! a = [2, 1, 1, 3, 1, 0, 0, 1]
  block
    real(8) :: aa(8), bb(4)
    integer :: fl

    aa = (/2.d0, 1.d0, 1.d0, 3.d0, 1.d0, 0.d0, 0.d0, 1.d0/)
    call girl(aa, bb, 2, 2, fl)

    call assert_equal_int(fl, 0, '2x2 inversion succeeds (flag=0)')

    ! A^-1 = (1/5) * [[3, -1], [-1, 2]]
    ! Column-major: col1 = [3/5, -1/5], col2 = [-1/5, 2/5]
    call assert_approx(bb(1), 0.6d0, 1.d-12, 'A^-1(1,1) = 3/5')
    call assert_approx(bb(2), -0.2d0, 1.d-12, 'A^-1(2,1) = -1/5')
    call assert_approx(bb(3), -0.2d0, 1.d-12, 'A^-1(1,2) = -1/5')
    call assert_approx(bb(4), 0.4d0, 1.d-12, 'A^-1(2,2) = 2/5')

    ! Verify A * A^-1 = I by manual multiplication
    ! (A * A^-1)(1,1) = 2*(3/5) + 1*(-1/5) = 6/5 - 1/5 = 1
    prod = 2.d0 * bb(1) + 1.d0 * bb(2)  ! A(1,:) * A^-1(:,1)
    call assert_approx(prod, 1.d0, 1.d-12, 'A*A^-1 (1,1) = 1')

    prod = 2.d0 * bb(3) + 1.d0 * bb(4)  ! A(1,:) * A^-1(:,2)
    call assert_approx(prod, 0.d0, 1.d-12, 'A*A^-1 (1,2) = 0')

    prod = 1.d0 * bb(1) + 3.d0 * bb(2)  ! A(2,:) * A^-1(:,1)
    call assert_approx(prod, 0.d0, 1.d-12, 'A*A^-1 (2,1) = 0')

    prod = 1.d0 * bb(3) + 3.d0 * bb(4)  ! A(2,:) * A^-1(:,2)
    call assert_approx(prod, 1.d0, 1.d-12, 'A*A^-1 (2,2) = 1')
  end block

  ! 1x1 matrix: A = [5], A^-1 = [0.2]
  block
    real(8) :: aa1(2), bb1(1)
    integer :: fl1
    aa1 = (/5.d0, 1.d0/)
    call girl(aa1, bb1, 1, 1, fl1)
    call assert_equal_int(fl1, 0, '1x1 inversion succeeds')
    call assert_approx(bb1(1), 0.2d0, 1.d-12, '1/5 = 0.2')
  end block

end subroutine test_girl

end program test_smallfunc
