!> \file interpolmod.f95
!> \brief Interpolation module
!>
!> Linear interpolation
module interpolation

  implicit none

  private
  public :: indic,fipoi,fipoi1,fipoi2,quint
  public:: qua,quad,quad_gg,flin,spline,splint,getd

contains

!=======================================================================
  integer function indic(x0,x,m)
!> recherche rapide de la position d une valeur x0 dans une table
!! monotone croissante ou decroissante de m nombres x(i)
! si  k = indice(x0,x,m)  on aura   x0 compris entre x(k) et x(k+1)
!                              ou   x0 = x(k).
! si x0 exterieur a la table indice =   1 si x0 du cote de x(1)
!                            indice = m-1 si x0 du cote de x(m)
!-----------------------------------------------------------------------
    integer,intent(in):: m
    real(8),intent(in):: x0
    real(8),dimension(m),intent(in):: x

    integer:: i=0,k,n
!-----------------------------------------------------------------------
    n = m
    k = 1

    do while (n-k-1 /= 0)
      i=(k+n)/2
      if((x(i)-x0)*(x(n)-x(1)) <= 0) then
        k = i
      else
        n = i
      endif
    enddo

    indic = k
    return

  end function indic
!=======================================================================
  real(8) function fipoi(x,n,a,b)
!-----------------------------------------------------------------------
    real(8),intent(in)::x
    integer,intent(in)::n
    integer::k
    real(8), dimension(n),intent(in):: a,b
!-----------------------------------------------------------------------
    k=indic(x,a,n)
    fipoi=b(k)+(b(k+1)-b(k))*(x-a(k))/(a(k+1)-a(k))

    return

  end function fipoi
!=======================================================================
  subroutine fipoi1(x,n,a,b,r1)
!-----------------------------------------------------------------------
! interpolation lineaire pour trouver r1 a partir de x en interpolant
! entre a et b

! utilise le sous-programme search, sous-programme d interpolation.
!-----------------------------------------------------------------------
    implicit none

    integer,intent(in)::n
    real(8),intent(in)::x
    real(8),intent(out)::r1

    integer::k
    real(8), dimension(n),intent(in):: a,b
!-----------------------------------------------------------------------
    k=indic(x,a,n)
! interpolation lineaire.
    r1=b(k)+(b(k+1)-b(k))*(x-a(k))/(a(k+1)-a(k))

    return

  end subroutine fipoi1
!=======================================================================
  subroutine fipoi2(x,n,a,b,r1,r2)
!-----------------------------------------------------------------------
! interpolation lineaire pour trouver r1 et sa derivee r2 a partir de x
! en interpolant entre a et b

! utilise le sous-programme search, sous-programme d interpolation.
!-----------------------------------------------------------------------
    implicit none

    integer,intent(in)::n
    real(8),intent(in)::x
    real(8),intent(out)::r1,r2

    integer::k
    real(8), dimension(n):: a,b
!-----------------------------------------------------------------------
    k=indic(x,a,n)
! interpolation lineaire.
    r1=b(k)+(b(k+1)-b(k))*(x-a(k))/(a(k+1)-a(k))
    r2=(b(k+1)-b(k))/(a(k+1)-a(k))

    return

  end subroutine fipoi2
!=======================================================================
  subroutine quint(x,x0,h,y0,y1,y2,y)
!-----------------------------------------------------------------------
!...... Quadratic interpolation for equidistant points
!...... y0=y(x0),y1=y(x1),y2=y(x2); h=x1-x0=x2-x1;
!...... Computes y=y(x)
!-----------------------------------------------------------------------
    implicit none

    real(8),intent(in)::x0,x,h,y0,y1,y2
    real(8),intent(out)::y

    real(8)::d1,d2,t
!-----------------------------------------------------------------------
    d1 = y1 - y0
    d2 = y2 - 2.d0*y1 + y0
    t  = (x - x0)/h

    y  = y0 + t*d1 + 0.5d0*t*(t-1)*d2

    return
  end subroutine quint
!=======================================================================
  real(kind=8) function qua(x1,x2,x3,y1,y2,y3,x0)
!-----------------------------------------------------------------------
! --  Interpolation quadratique par son polynome de Lagrange
!     cf. Numerical Recipes p.80
!     x et y sont les coordonnees des trois points par lesquels
!     passe la parabole
!     x0 est la valeur pour laquelle on veut connaitre y
!     qua est la valeur de y en x0
!-----------------------------------------------------------------------
  implicit none

  real(kind=8),intent(in):: x1,x2,x3,y1,y2,y3,x0
  real(kind=8):: a23,a12,a13,a01,a02,a03
!-----------------------------------------------------------------------
  a23  = x2 - x3
  a12  = x1 - x2
  a13  = x1 - x3
  a01  = x0 - x1
  a02  = x0 - x2
  a03  = x0 - x3
  qua = a02 * a03 / a12 / a13 * y1 - a01 * a03 / a12 / a23 * y2 + a01 * a02 / a13 / a23 * y3

  return

  end function qua
!=======================================================================
  real(kind=8) function quad(ic,i,x_quad,y1,y2,y3,x1,x2,x3,dkap_quad)
!-----------------------------------------------------------------------
! this function performs a quadratic interpolation.
!-----------------------------------------------------------------------
  implicit none

  save

  integer,intent(in):: ic,i
  real(kind=8),intent(in):: x1,x2,x3,y1,y2,y3
  real(kind=8):: c1,c2,c3,x_quad,dkap_quad
  real(kind=8),dimension(3)::  xx_quad,yy_quad
  real(kind=8),dimension(30):: xx12,xx13,xx23,xx1sq,xx1pxx2
!-----------------------------------------------------------------------
  xx_quad(1)=x1
  xx_quad(2)=x2
  xx_quad(3)=x3
  yy_quad(1)=y1
  yy_quad(2)=y2
  yy_quad(3)=y3
  if(ic == 0) then
    xx12(i)=1.d0/(xx_quad(1)-xx_quad(2))
    xx13(i)=1.d0/(xx_quad(1)-xx_quad(3))
    xx23(i)=1.d0/(xx_quad(2)-xx_quad(3))
    xx1sq(i)=xx_quad(1)*xx_quad(1)
    xx1pxx2(i)=xx_quad(1)+xx_quad(2)
  endif
  c3=(yy_quad(1)-yy_quad(2))*xx12(i)
  c3=c3-(yy_quad(2)-yy_quad(3))*xx23(i)
  c3=c3*xx13(i)
  c2=(yy_quad(1)-yy_quad(2))*xx12(i)-(xx1pxx2(i))*c3
  c1=yy_quad(1)-xx_quad(1)*c2-xx1sq(i)*c3
  dkap_quad=c2+(x_quad+x_quad)*c3
  quad=c1+x_quad*(c2+x_quad*c3)

  return

  end function quad
!=======================================================================
  real(kind=8) function quad_gg(x1,x2,x3,y1,y2,y3,x0)
!-----------------------------------------------------------------------
! --  Interpolation quadratique par son polynome de Lagrange
!     cf. Numerical Recipes p.80
!     x et y sont les coordonnees des trois points par lesquels
!     passe la parabole
!     x0 est la valeur de x en laquelle on veut connaitre la valeur de la derivee
!     quad_gg est la derivee en x0
!-----------------------------------------------------------------------
  implicit none

  real(kind=8),intent(in):: x1,x2,x3,y1,y2,y3,x0
  real(kind=8):: a01,a02,a03,a12,a13,a23
!-----------------------------------------------------------------------
  a23  = x2 - x3
  a12  = x1 - x2
  a13  = x1 - x3
  a01  = x0 - x1
  a02  = x0 - x2
  a03  = x0 - x3

  quad_gg = ((a02 + a03)/a12/a13*y1) - ((a01 + a03)/a12/a23*y2) + ((a01 + a02)/a13/a23*y3)

  return

  end function quad_gg
!=======================================================================
  real(kind=8) function flin(x1,x2,y1,y2,x0)
!-----------------------------------------------------------------------
! --  Interpolation lineaire pouvant aussi extrapoler en x0
!     (x1,y1) et (x2,y2) sont les deux points par lesquels passe la droite.
!     on donne x0 et flin renvoie la valeur de y correspondante.
!-----------------------------------------------------------------------
  implicit none

  real(kind=8),intent(in):: x1,x2,y1,y2,x0
  real(kind=8):: a01,a02,a12
!-----------------------------------------------------------------------
  a12  = x1 - x2
  a01  = x0 - x1
  a02  = x0 - x2

  flin =  (a02 * y1 - a01 * y2)/a12

  return

  end function flin
!=======================================================================
  subroutine spline(X,Y,N,Y2)
!-----------------------------------------------------------------------
  implicit none

  integer, parameter:: nmax=100
  integer:: i,k
  integer,intent(in):: n
  real(kind=8):: yp1,ypn,sig,p,qn,un
  real(kind=8), dimension(n),intent(in):: x,y
  real(kind=8), dimension(n),intent(out):: y2
  real(kind=8), dimension(nmax):: u
!-----------------------------------------------------------------------
!     FIRST DERIVATIVES AT END POINTS USING CUBIC FIT
  yp1=((y(3)-y(1))*(x(2)-x(1))**2.d0-(y(2)-y(1))*(x(3)-x(1))**2.d0)/ &
      ((x(3)-x(1))*(x(2)-x(1))*(x(2)-x(3)))
  ypn=((y(n-2)-y(n))*(x(n-1)-x(n))**2.d0-(y(n-1)-y(n))*(x(n-2)-x(n))**2.d0)/ &
      ((x(n-2)-x(n))*(x(n-1)-x(n))*(x(n-1)-x(n-2)))

  y2(1)=-0.5d0
  u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  do i=2,n-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2.d0
    y2(i)=(sig-1.d0)/p
    u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  enddo
  qn=0.5d0
  un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
  do k=n-1,1,-1
    y2(k)=y2(k)*y2(k+1)+u(k)
  enddo

  return

  end subroutine spline
!=======================================================================
  subroutine splint(XA,YA,N,Y2A,X,Y,YP)
!-----------------------------------------------------------------------
  implicit none

  integer,intent(in):: n
  integer:: klo,khi,k
  real(kind=8):: h,a,b
  real(kind=8),intent(in):: x
  real(kind=8),intent(out):: y,yp
  real(kind=8),dimension(n):: xa,ya,y2a
!-----------------------------------------------------------------------
  klo=1
  khi=n
  do while (khi-klo>1)
    k=(khi+klo)/2
    if(xa(k) > x) then
      khi=k
    else
      klo=k
    endif
  enddo
  h=xa(khi)-xa(klo)
  if (h == 0.d0) then
    stop 'Bad XA input.'
  endif
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)+((a**3.d0-a)*y2a(klo)+(b**3.d0-b)*y2a(khi))*(h**2.d0)/6.d0
  yp=0.05d0*((-ya(klo)+ya(khi))/h+( -(3.d0*a**2.d0-1.d0)*y2a(klo)+(3.d0*b**2.d0-1.d0)*y2a(khi) )*h/6.d0 )

  return

  end subroutine splint
!=======================================================================
  subroutine getd(f,n,d,fp1,fpn)
!-----------------------------------------------------------------------
!  SIMPLIFIED CODE FOR SPLINE COEFFICIENTS, FOR CASE OF INTERVALS
!  OF UNITY.
!-----------------------------------------------------------------------
  implicit none
  integer,intent(in):: n
  integer:: j
  real(kind=8), intent(out):: fp1,fpn
  real(kind=8),dimension(n):: f,d
  real(kind=8),dimension(85):: t
!-----------------------------------------------------------------------
  fp1=(-11.d0*f(1)+18.d0*f(2)-9.d0*f(3)+2.d0*f(4))/6.d0
  fpn=(11.d0*f(n)-18.d0*f(n-1)+9.d0*f(n-2)-2.d0*f(n-3))/6.d0

  d(1)=-.5d0
  t(1)=.5d0*(-f(1)+f(2)-fp1)

  do j=2,n-1
    d(j)=-1.d0/(4.d0+d(j-1))
    t(j)=-d(j)*(f(j-1)-2.d0*f(j)+f(j+1)-t(j-1))
  enddo

  d(n)=(fpn+f(n-1)-f(n)-t(n-1))/(2.d0+d(n-1))

  do j=n-1,1,-1
    d(j)=d(j)*d(j+1)+t(j)
  enddo

  return

  end subroutine getd
!=======================================================================
end module interpolation
