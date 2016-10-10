!> \file interpolmod.f95
!> \brief Interpolation module
!>
!> Linear interpolation
module interpolmod

  private
  public :: fipoi

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  real(kind=8) function fipoi(x,n,a,b)

    real(kind=8),intent(in)::x
    integer,intent(in)::n
    integer::k
    real(kind=8), dimension(n),intent(in):: a,b

    k=indic(x,a,n)
    fipoi=b(k)+(b(k+1)-b(k))*(x-a(k))/(a(k+1)-a(k))

    return

  end function fipoi
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> recherche rapide de la position d une valeur x0 dans une table
!! monotone croissante ou decroissante de m nombres x(i)
  integer function indic(x0,x,m)

    real(kind=8),intent(in):: x0
    integer,intent(in):: m
    integer:: i=0,k,n
    real(kind=8), intent(in),dimension(m):: x

! si  k = indice(x0,x,m)  on aura   x0 compris entre x(k) et x(k+1)
!                              ou   x0 = x(k).
! si x0 exterieur a la table indice=1 si x0 du cote de x(1)
!                            indice = m-1 si x0 du cote de x(m)

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
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end module interpolmod
