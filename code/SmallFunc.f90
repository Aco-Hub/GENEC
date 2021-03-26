module SmallFunc

implicit none

private
public:: girl
public:: CheckProfile,SmoothProfile
public:: tridiago
public:: pos,neg_root
public:: primus,drimus,secun,dsecun,tertiu,dterti,expf10,exphi,ValInterp,OmInterp
public:: weighed_smoothing, threshold_smoothing

contains
!=======================================================================
  subroutine girl(a,b,n,m,flag)
!-----------------------------------------------------------------------
! matrix inversion by Gauss elimination
!-----------------------------------------------------------------------
  implicit none

  integer,intent(in)::n,m
  integer,intent(inout):: flag
  real(8),dimension(n*(n+m)),intent(inout):: a
  real(8),dimension(n*m),intent(out):: b

  integer::j,npm,nj,jj,j1,jm,i,ij,i1,i2,ik,k,jk
  real(8):: amax,zwi,faktor
!-----------------------------------------------------------------------
  flag=0
  npm = n+m
  do  j = 1,n
   nj = (j-1)*n
   jj = nj + j
   j1 = j + 1
   amax = abs (a(jj))
   jm = j
   if ( (j1 - n)  <= 0 ) then
     do  i = j1,n
      ij = nj + i
      if ((abs(a(ij))-amax) > 0.d0) then
        amax = abs (a(ij))
        jm = i
      endif
     enddo
     if ( (jm - j)  <  0 ) then
       flag=1
       write(*,*) ' JM=',jm,'J=',j
       write(*,*) ' APPEL GIRL AVEC n=',n,' et m=',m
       return
     else
       ! partial pivoting if pivot position is =0
       if ( (jm - j)  >  0 ) then
         i1 = jm + nj
         i2 = jj
         do i = j,npm
          zwi = a(i1)
          a(i1) = a(i2)
          a(i2) = zwi
          i1 = i1 + n
          i2 = i2 + n
         enddo
       endif
     endif
   endif
   if ( a(jj) == 0.d0 )then
     flag=2
     write(*,*) 'GIRL:', jj,j
     write(*,*) ' JM=',jm,'J=',j
     write(*,*) ' APPEL GIRL AVEC n=',n,' et m=',m
     return
   endif
   do i = 1,n
     ! elimination loop
    if ( i  /=  j ) then
      ij = nj + i
      ik = nj + i
      jk = jj
      faktor = - a(ij)/a(jj)
      do  k = j1,npm
       jk = jk + n
       ik = ik + n
       a(ik) = a(ik) + faktor * a(jk)
      enddo
    endif
   enddo
   ! division of pivot row
   jk = jj
   faktor = 1.d0/a(jj)
   do k = j1,npm
    jk = jk + n
    a(jk) = a(jk) * faktor
   enddo
  enddo
  i1 = n*m
  i2 = n*n
  do i = 1,i1
   i2 = i2 + 1
   b(i) = a(i2)
  enddo

  return

  end subroutine girl
!=======================================================================
subroutine CheckProfile(Prof,Profout,mesh,m)
!-----------------------------------------------------------------------
! Checks that no jumps occur in a curve, and if there is, replaces the point by the mean value of its neighbours
!-----------------------------------------------------------------------
  use evol,only: ldi,kindreal

  implicit none

  integer,intent(in):: m
  real(kindreal),dimension(ldi),intent(in):: Prof,mesh
  real(kindreal),dimension(ldi),intent(out):: Profout

  integer, parameter:: win=2
  integer:: i,j
  real(kindreal),parameter:: Tol=1.10d0
  real(kindreal):: mean
!-----------------------------------------------------------------------
  Profout(:) = Prof(:)

  do i=win+1,m-(win+1)
    mean=0.d0
    do j=1,win
      mean=mean + Prof(i-j)+Prof(i+j)
    enddo
    mean=(mean+Prof(i))/real(2*win+1)

    if (Prof(i)/mean > Tol .or. Prof(i)/mean < 1.d0/Tol) then
      Profout(i) = (Profout(i+1)-Profout(i-1))/(mesh(i+1)-mesh(i-1))*(mesh(i)-mesh(i-1)) + Profout(i-1)
    endif
  enddo

  return

end subroutine CheckProfile
!=======================================================================
subroutine SmoothProfile(Prof,Profout,Derivout,mesh,min,max,WinSize,m)
!-----------------------------------------------------------------------
! Smoothing of a curve by the locally weighted scatterplot smoothing technique
!-----------------------------------------------------------------------
  use evol,only: ldi,kindreal
  use inputparam,only: idebug,verbose
  use nagmod,only: e02acf

  implicit none

  integer,intent(in):: min,max,WinSize,m
  real(kindreal),dimension(ldi),intent(in):: Prof,mesh
  real(kindreal),dimension(ldi),intent(out):: Profout,Derivout

  integer:: i,j
  integer,parameter:: ordre=1
  real(kindreal):: deviation
  real(kindreal),dimension(2*WinSize+1):: Xdata,Ydata
  real(kindreal),dimension(ordre+1):: aaa
!-----------------------------------------------------------------------
  Xdata(:) = 0.d0
  Ydata(:) = 0.d0
  aaa(:) = 0.d0
  deviation = 0.d0

  if (max-min < 2*WinSize+1) then
    if(verbose) write(*,*)'Zone too small to fit, only',max-min,' layers'
    return
  endif

  if (idebug > 0) then
    write(*,*) 'in SmoothProfile, WinSize,min,max:',WinSize,min,max
  endif
  do i=min,max
   if (i - Winsize < 1) then
     Xdata(WinSize+2-i:2*WinSize+1)=mesh(1:i+WinSize)
     Ydata(WinSize+2-i:2*WinSize+1)=Prof(1:i+WinSize)
   elseif (i + Winsize > m) then
     Xdata(1:m-i+WinSize+1)=mesh(i-WinSize:m)
     Ydata(1:m-i+WinSize+1)=Prof(i-WinSize:m)
   else
     Xdata=mesh(i-WinSize:i+WinSize)
     Ydata=Prof(i-WinSize:i+WinSize)
   endif
   if (i < min + WinSize) then
     call e02acf(Xdata(min-i+WinSize+1:2*WinSize+1),Ydata(min-i+WinSize+1:2*WinSize+1),WinSize+i-min+1,aaa,ordre+1,deviation)
   elseif (i > max - WinSize) then
     call e02acf(Xdata(1:max-i+WinSize+1),Ydata(1:max-i+WinSize+1),WinSize-i+max+1,aaa,ordre+1,deviation)
   else
     call e02acf(Xdata,Ydata,2*WinSize+1,aaa,ordre+1,deviation)
   endif
   Profout(i)=0.d0
   Derivout(i)=0.d0
   do j=1,ordre+1
     Profout(i)=Profout(i)+aaa(j)*mesh(i)**real(j-1)
   enddo
   do j=2,ordre+1
     Derivout(i)=Derivout(i)+real(j-1)*aaa(j)*mesh(i)**real(j-2)
   enddo
  enddo

  return

end subroutine SmoothProfile
!=======================================================================
  subroutine tridiago(at,bt,ct,vvx,k)
!-----------------------------------------------------------------------
    use evol,only: ldi,kindreal

    implicit none

    integer,intent(in):: k
    real(kindreal),dimension(ldi),intent(in):: at,bt,ct
    real(kindreal),dimension(ldi),intent(inout):: vvx

    integer:: i
    real(kindreal):: bet
    real(kindreal),dimension(ldi):: rt,gam
!-----------------------------------------------------------------------
    rt(1:k)=vvx(1:k)

    bet=bt(1)
    vvx(1)=rt(1)/bet
    do i=2,k
    if (vvx(i-1) < 1.d-75) then
      vvx(i-1) = 0.d0
    endif
     gam(i)=ct(i-1)/bet
     bet=bt(i)-at(i)*gam(i)
     if (bet == 0.0d0) then
        !        print*,max(abs(at(:)),1.d0),max(abs(bt(:)),1.d0),max(abs(ct(:)),1.d0)
        print*,"at(i),bt(i),ct(i),i"
         print*,at(i),bt(i),ct(i),i
       stop ' tridag failed'
     endif
     vvx(i)=(rt(i)-at(i)*vvx(i-1))/bet
    enddo
    do i=k-1,1,-1
     vvx(i)=vvx(i)-gam(i+1)*vvx(i+1)
    enddo
    if (isnan(at(i))) stop '"at" is a NaN'
    if (isnan(bt(i))) stop '"bt" is a NaN'
    if (isnan(ct(i))) stop '"ct" is a NaN'
   
    return

  end subroutine tridiago
!=======================================================================
  subroutine weighed_smoothing(dd, n, ns, preserve_sign, ddold)
    use evol,only: ldi,kindreal
!     based on routine written by S.-C. Yoon, 18 Sept. 2002
!     for smoothing  any variable (dd) with size n over 2*ns+1 cells.
  real(kindreal), intent(inout) :: dd(:) ! (n)
  integer, intent(in) :: n, ns
  logical, intent(in) :: preserve_sign
  real(kindreal), intent(inout) :: ddold(:) ! (n) work array
  integer :: nweight, mweight, i, j, k
  real(kindreal) :: weight(2*ns+1), sweight, v0

  do i = 1,n
     ddold(i) = dd(i)
  end do

  !--preparation for smoothing --------
  nweight = ns
  mweight = 2*nweight+1
  do i = 1,mweight
     weight(i) = 0d0
  end do
  weight(1) = 1d0
  do i = 1,mweight-1
     do j = i+1,2,-1
        weight(j) = weight(j) + weight(j-1)
     end do
  end do

  !--smoothing ------------------------
  do i=2,n-1
     sweight=0d0
     dd(i)=0d0
     v0 = ddold(i)
     do j = i, max(1,i-nweight), -1
        k=j-i+nweight+1
        if (preserve_sign .and. v0*ddold(j) <= 0) exit
        sweight = sweight+weight(k)
        dd(i) = dd(i)+ddold(j)*weight(k)
     end do
     do j = i+1, min(n,i+nweight)
        k=j-i+nweight+1
        if (preserve_sign .and. v0*ddold(j) <= 0) exit
        sweight = sweight+weight(k)
        dd(i) = dd(i)+ddold(j)*weight(k)
     end do
     if (sweight > 0) then
        sweight = 1d0/sweight
        dd(i) = dd(i)*sweight
     end if
  end do

end subroutine weighed_smoothing
!================================================================================
subroutine threshold_smoothing (dd, dd_thresh, n, ns, preserve_sign, ddold)
    use evol,only: ldi,kindreal
  ! Same as weighed_smoothing, but only smooth contiguous regions where |dd| >= dd_thresh
  ! NOTE: this can be adapted to any smoothing algorithm
  real(kindreal), intent(inout) :: dd(:)    ! (n)
  real(kindreal), intent(in)    :: dd_thresh
  integer, intent(in)     :: n
  integer, intent(in)     :: ns
  logical, intent(in)     :: preserve_sign
  real(kindreal), intent(inout) :: ddold(:) ! (n) work array

  logical :: in_region
  integer :: i
  integer :: i_a
  integer :: i_b

  ! Process regions

  in_region = .FALSE.

  i_a = 1
  do i = 1, n

     if (in_region) then

        if (ABS(dd(i)) < dd_thresh) then
           i_b = i-1
           if (i_b > i_a) call weighed_smoothing(dd(i_a:i_b), i_b-i_a+1, ns, preserve_sign, ddold(i_a:i_b))
           in_region = .FALSE.
        endif

     else
        if (ABS(dd(i)) >= dd_thresh) then
           i_a = i
           in_region = .TRUE.
        endif

     end if

  end do

  ! Handle the final region

  if (in_region) then

     i_b = n
     if (i_b > i_a) call weighed_smoothing(dd(i_a:i_b), i_b-i_a+1, ns, preserve_sign, ddold(i_a:i_b))

  endif

  ! Finish

  return
end subroutine threshold_smoothing
!=======================================================================
!> This function protects the possibly negative value of the variable raised to a real power
!! @parameter[in]  a: variable; b: power
  real(8) function neg_root(a,b)
!-----------------------------------------------------------------------
    implicit none

    real(8),intent(in):: a,b
    real(8):: c
!-----------------------------------------------------------------------
    c = abs(a)**b
    if (a < 0.d0) then
      neg_root = -c
    else
      neg_root = c
    endif

end function neg_root
!=======================================================================
  real(8) function primus(a,b,e1,e2,t9)
!-----------------------------------------------------------------------
    implicit none

    real(8),intent(in):: a,b,e1,e2,t9
!-----------------------------------------------------------------------
    if (-b/t9**e2 > -706.d0) then
      primus=a/t9**e1*exp(-b/t9**e2)
    else
      primus=0.d0
    endif

    return

  end function primus
!=======================================================================
  real(8) function drimus(b,e1,e2,t9)
!-----------------------------------------------------------------------
    implicit none

    real(8),intent(in):: b,e1,e2,t9
!-----------------------------------------------------------------------
    drimus=-e1+e2*b/t9**e2

    return

  end function drimus
!=======================================================================
  real(8) function secun(a,b,c,e1,e2,e3,t9)
!-----------------------------------------------------------------------
    implicit none

    real(8),intent(in):: a,b,c,e1,e2,e3,t9
!-----------------------------------------------------------------------
    secun=a/t9**e1*exp(-b/t9**e2-(t9/c)**e3)

    return

  end function secun
!=======================================================================
  real(8) function dsecun(b,c,e1,e2,e3,t9)
!-----------------------------------------------------------------------
    implicit none

    real(8),intent(in):: b,c,e1,e2,e3,t9
!-----------------------------------------------------------------------
    dsecun=-e1+e2*b/t9**e2-e3*(t9/c)**e3

    return

  end function dsecun
!=======================================================================
  real(8) function tertiu(a,b,c,d,e,f,e1,e2,e3,e4,e5,t9)
!-----------------------------------------------------------------------
    implicit none

    real(8),intent(in):: a,b,c,d,e,f,e1,e2,e3,e4,e5,t9
!-----------------------------------------------------------------------
    tertiu=a+b*t9**e1+c*t9**e2+d*t9**e3+e*t9**e4+f*t9**e5

    return

  end function tertiu
!=======================================================================
  real(8) function dterti(b,c,d,e,f,e1,e2,e3,e4,e5,t9)
!-----------------------------------------------------------------------
    implicit none

    real(8),intent(in):: b,c,d,e,f,e1,e2,e3,e4,e5,t9
!-----------------------------------------------------------------------
    dterti=e1*b*t9**e1+e2*c*t9**e2+e3*d*t9**e3+e4*e*t9**e4+e5*f*t9**e5

    return

  end function dterti
!=======================================================================
  real(8) function expf10(x)
!-----------------------------------------------------------------------
    implicit none

    real(8),intent(in):: x

    real(8):: vloge
!-----------------------------------------------------------------------
    vloge = log10(exp(1.d0))
    expf10 = exp(x/vloge)

    return

  end function expf10
!=======================================================================
  real(8) function exphi(x)
!-----------------------------------------------------------------------
    implicit none

    real(8),intent(in):: x
!-----------------------------------------------------------------------
    if ((abs(x) - 3.d-2) <= 0.d0) then
      exphi = -((x*(1.d0/3.d0)+1.d0)*x*0.5d0+1.d0)*x
    else
      exphi = 1.d0 - exp(x)
    endif

    return

  end function exphi
!=======================================================================
  integer function pos(x0,x,m)
!-----------------------------------------------------------------------
!  -- recherche de la position d'une valeur x0 dans une table monotone
!     croissante ou decroissante de m nombre x(i)

!     si k = pos(x0,x,m) on aura x0 compris entre x(k) et x(k+1)
!     ou x0 = x(k)

!     si x0 est en dehors de  la table pos = 1 si x0 du cote de x(1)
!     pos = m-1 si x0 du cote de x(m)
!-----------------------------------------------------------------------
    implicit none

    integer,intent(in):: m
    real(8),intent(in):: x0

    integer:: n,k,i
    real(8), dimension(m):: x
!-----------------------------------------------------------------------
    n=m
    k=1
    do while (n-k-1  /=  0)
     i = (k+n)/2
     if ((x(i)-x0)*(x(n)-x(1))  <= 0) then
       k=i
     else
       n=i
     endif
    enddo

    pos = k

    return

  end function pos
!=======================================================================
  real(8) function ValInterp(x1,x2,y1,y2,y3)
!-----------------------------------------------------------------------
      implicit none

      real(8), intent(in)::x1,x2,y1,y2,y3
!-----------------------------------------------------------------------
      ValInterp = (x1*(y2-y3)+x2*(y3-y1))/(y2-y1)

      return

      end function ValInterp
!=======================================================================
  real(8) function OmInterp(x1,x2,y1,y2,y3,z1,z2,z3)
!-----------------------------------------------------------------------
  implicit none

  real(8), intent(in)::x1,x2,y1,y2,y3,z1,z2,z3
!-----------------------------------------------------------------------
  OmInterp = (x1*(y2-y3)*z1+x2*(y3-y1)*z2)/((y2-y1)*z3)

  return

  end function OmInterp
!=======================================================================
end module SmallFunc
