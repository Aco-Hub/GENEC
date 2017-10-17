module convection

  use evol,only: kindreal
  use const,only: cst_G,Msol
  use caramodele,only: gms,nwmd
  use inputparam,only: idifcon,iover,iunder,dovhp,dunder,verbose
  use strucmod,only: m,adgrad,zensi,q,r,p,rho

  implicit none

  integer,parameter:: iirad=300,ixzc=100

  integer,save:: iconra,nzrad,nzcon,muce,mlce,muci,mlci,jzint,jwint,iidraw
  integer,dimension(iirad),save:: nxzrad,nxzcon
  integer,dimension(ixzc),save:: izc

  real(kindreal),save:: xzce,xunder,xz,xover,r_core,qmnc,qbc,BaseZC
  real(kindreal),dimension(40):: drawcon
  real(kindreal),dimension(iirad),save:: xmomin,xmocin
  real(kindreal),dimension(ixzc),save:: xzc

private
public :: bordn,over1,unders,fconva,rechzco,konvek
public :: qmnc,qbc,xover,r_core,xmomin,xmocin
public :: nzrad,nxzrad,nzcon,nxzcon,iconra,muce,mlce,muci,mlci,jzint,jwint,ixzc,izc,xzc
public :: CZdraw,iidraw,drawcon,BaseZC

contains
!=======================================================================
subroutine bordn
!----------------------------------------------------------------------
! RECHERCHE LES BORDS DES ZONES CONVECTIVES
! SI IOVER=1: DETERMINATION DU BORD DU COEUR AVEC OVERSHOOTING
! SI IUNDER=1: DETERMINATION DE LA BASE DE LA ZC EXTERNE AVEC UNDERSHOOTING
!----------------------------------------------------------------------
  use interpolation,only: fipoi

  implicit none

  integer:: k,int,nt,n,i1,i2,mwz,mwy
  real(kindreal):: qm1,qm2,xlfer,binf,xlover,bord
!----------------------------------------------------------------------
  k=m-1
  int=1
  jwint=0
  do nt=1,ixzc
   xzc(nt)=0.d0
  enddo
  do while (k /= 0)
   n=0
   do while (adgrad(k) <= 0.d0)
     n=n+1
     if (k <= 1) then
       k=0
       exit
     endif
     k=k-1
   enddo
   if (n > 0) then
     i1=k
     i2=k+n
     if (i2 == m-1) then
       xzc(int)=10000.d0
     else
       qm1=(1.d0-exp(q(i2)))
       qm2=(1.d0-exp(q(i2+1)))
       xzc(int)=(qm2*adgrad(i2)-qm1*adgrad(i2+1))/(adgrad(i2)-adgrad(i2+1))
     endif
     if (i1 == 0) then
       xzc(int+1)=1.d0
     else
       qm1=(1.d0-exp(q(i1)))
       qm2=(1.d0-exp(q(i1+1)))
       xzc(int+1)=(qm1*adgrad(i1+1)-qm2*adgrad(i1))/(adgrad(i1+1)-adgrad(i1))
     endif
     int=int+2
     if (int > ixzc) then
       print*,' ATTENTION XZC SS-DIM: int= ',int
       rewind(222)
       write (222,*) nwmd,': xzc sous-dim ==> STOP'
       stop
     endif
     jwint=(int-1)/2
   endif
   if (k <= 1) then
     exit
   endif
   k=k-1
  enddo

  if (xzc(1) == 10000.d0) then
    mwz=m
    if (xover /= 0.d0) then
      xlover=log(xover)
      bord=fipoi(xlover,mwz,r,q)
      qmnc=1.d0-exp(bord)
      r_core = xover
    else
      qmnc=xzc(2)
      r_core=exp(fipoi(log(1.0d0-qmnc),mwz,q,r))
    endif
  else
    qmnc=0.d0
  endif
  mwy=m
  if (iunder == 0 .or. xunder == 0.d0) then
    qbc=1.0d0
  else
    xlfer=log(xunder)
    binf=fipoi(xlfer,mwy,r,q)
    qbc=1.d0-exp(binf)
  endif

  return

end subroutine bordn
!=======================================================================
subroutine over1
!----------------------------------------------------------------------
! CALCUL DE L'OVERSHOOTING
!----------------------------------------------------------------------
  use inputparam,only: itminc

  implicit none

  integer:: l
  real(kindreal):: wy2,wy1,wx2,wx1,hpp
!----------------------------------------------------------------------
  if (adgrad(m-1) > 0.d0) then
    xover=0.d0
    return
  endif

  l=2
  do while (l <= m-1)
   if (adgrad(m-l) > 0.d0) then
     exit
   endif
   l=l+1
  enddo

  wy2=adgrad(m-l)
  wy1=adgrad(m-l+1)
  wx2=exp(r(m-l))
  wx1=exp(r(m-l+1))
  xz=wx1-wy1*((wx1-wx2)/(wy1-wy2))
  hpp=exp(p(m-l))*wx2
  hpp=hpp/(cst_G*Msol*exp(rho(m-l)))
  hpp=hpp*wx2/(gms*(1.d0-exp(q(m-l))))
  if (hpp < xz) then
    xover=xz+dovhp*hpp
    write(3,'(1x,a,3(3x,1pe11.4))') ' OVERSHOOTING XZ,HPP,XOVER= ',xz,hpp,xover
  else
    xover=xz*(1.d0+dovhp)
    if (verbose .or. itminc == 1) then
      print*,' HP >= RCONV'
    endif
    if (itminc == 1) then
      write(3,*)' HP >= RCONV'
      write(3,'(1x,a,3(3x,1pe11.4))')' OVERSHOOTING XZ,HPP,XOVER= ',xz,hpp,xover
    endif
  endif

  return

end subroutine over1
!=======================================================================
subroutine unders
!----------------------------------------------------------------------
! CALCUL DE L'UNDERSHOOTING
!----------------------------------------------------------------------
  implicit none

  integer:: l
  real(kindreal):: wx1,wx2,wy1,wy2,hpp
!----------------------------------------------------------------------
  xunder=0.d0
  if (adgrad(1) > 0.d0) then
    return
  endif

  l=2
  do while (l <= m)
   if (adgrad(l) > 0.d0) then
     exit
   endif
   l=l+1
  enddo
  wy2=adgrad(l)
  wy1=adgrad(l-1)
  wx2=exp(r(l))
  wx1=exp(r(l-1))
  xzce=wx1-wy1*((wx2-wx1)/(wy2-wy1))
! Hpp : echelle de pression
  hpp=exp(p(l))*wx2
  hpp=hpp/(Msol*cst_G*exp(rho(l)))
  hpp=hpp*wx2/(gms*(1.d0-exp(q(l))))
  xunder=xzce-dunder*hpp
  write(3,'(2x,a,1x,3(3x,1pe11.4))') 'UNDERSHOOTING XZ,HPP,XUNDER=',xzce,hpp,xunder

  return

end subroutine unders
!=======================================================================
subroutine fconva(dwdo2,f,fz,x,a,b,c)
!----------------------------------------------------------------------
! FONCTION DE LA CONVECTION POUR L  =  ALPHA * HRHO
!----------------------------------------------------------------------
  implicit none

  real(kindreal),intent(in):: dwdo2,a,b,c
  real(kindreal),intent(out):: f,fz,x

  real(kindreal), parameter:: e=1.d-20
  real(kindreal):: z2,y,w,yz,xz
!----------------------------------------------------------------------
  z2 = dwdo2*dwdo2
  y = dwdo2*(z2 + (8.d0/9.d0)*(dwdo2 + 2.d0))

  if (z2 < e) then
    w = dwdo2*sqrt((y/dwdo2)**2.d0 + a*z2*z2*z2)
  else
    w = z2*sqrt((y/z2)**2.d0 + a*z2*z2)
  end if

  x = ((9.d0/16.d0)*b*(y + w))**0.25d0
  f = c*x**3.d0*(x - 1.d0) - dwdo2*(dwdo2 + 2.d0)
  yz = 3.d0*z2 + (16.d0/9.d0)*(dwdo2 + 1.d0)
  xz = (9.d0/64.d0)*b*((1.d0+y/w)*yz+4.d0*a*(z2/w)*z2*z2*dwdo2)/x**3.d0
  fz = c*x*x*(4.d0*x - 3.d0)*xz - 2.d0*(dwdo2 + 1.d0)

  return

end subroutine fconva
!=======================================================================
subroutine rechzco
!----------------------------------------------------------------------
! Recherche des limites des zones radiatives

! repere le type de configuration
! iconra=0 * entierement convective
! iconra=1 * entierement radiative
! iconra=2 rad----rad
! iconra=3 conv---conv
! iconra=4 rad----conv
! iconra=5 conv---rad
! repere les limites des zones convectives
! calcule le moment d inertie des zones convectives
!----------------------------------------------------------------------
  use inputparam, only: itminc

  implicit none

  integer:: n,jrad,nrad,izin,izex,jord,kn,jinte,jexte,last
!----------------------------------------------------------------------
! Detection des zones radiatives
! On va du centre vers le bord

! Compteur de couches
  n=m
! nombre de zones radiatives
  jrad=0
! compteur du nombre de couches dans la zone radiative
  nrad=0

  do while (n  >=  1)
   if (zensi(n) <= 0.d0) then
     nrad=nrad+1
     if (n /= 1) then
       n=n-1
     else
       if (nrad /= 0) then
         izin=n+nrad-1
         izex=n
         jrad=jrad+1
         jord=2*jrad-1
         nxzrad(jord)=izin
         nxzrad(jord+1)=izex
       endif
       exit
     endif
   else
     if (nrad /= 0) then
       izin=n+nrad
       izex=n+1
       jrad=jrad+1
       jord=2*jrad-1
       nxzrad(jord)=izin
       nxzrad(jord+1)=izex
       if (n /= 1) then
         n=n-1
         nrad=0
       else
         exit
       endif
     else
       if (n /= 1) then
         n=n-1
       else
         exit
       endif
     endif
   endif
  enddo

! impression pour controle
  nzrad=jrad
  if (verbose .or. itminc == 1) then
    write(*,*) ' nombre de zones radiatives=', nzrad,m
  endif
  do kn=1,jrad
   jinte=2*kn-1
   jexte=2*kn
   if (verbose .or. itminc == 1) then
     write(*,*) ' ', kn,' ', nxzrad(jinte),' ', nxzrad(jexte)
   endif
   write(3,*) 'zones rad: ', kn,' ', 1.d0-exp(q(nxzrad(jinte))),' ',1.d0-exp(q(nxzrad(jexte)))
   write(3,*) nxzrad(jinte),nxzrad(jexte)
  enddo

  if (nzrad == 0) then
    nzcon=1
    nxzcon(1)=m
    nxzcon(2)=1
    iconra=0
    return
  endif

  last=2*nzrad
  if (nxzrad(1) == m .and. nxzrad(last) == 1) then
    if (nzrad == 1) then
      nzcon=0
      iconra=1
      return
    endif
    nzcon=nzrad-1
    iconra=2
    do n=1,nzcon
     nxzcon(2*n-1)=nxzrad(2*n)-1
     nxzcon(2*n)=nxzrad(2*n+1)+1
    enddo
    muci=nxzcon(2*nzcon)
    mlci=nxzcon(2*nzcon-1)
  else if (nxzrad(1) /= m .and. nxzrad(last) /= 1) then
    iconra=3
    nzcon=nzrad+1
    nxzcon(1)=m
    nxzcon(2)=nxzrad(1)+1
    nxzcon(2*nzcon-1)=nxzrad(2*nzcon-2)-1
    nxzcon(2*nzcon)=1
    if (nzrad > 1) then
      do n=2,nzcon-1
       nxzcon(2*n-1)=nxzrad(2*n-2)-1
       nxzcon(2*n)=nxzrad(2*n-1)+1
      enddo
    endif
    muci=nxzcon(2*(nzcon-1))
    mlci=nxzcon(2*(nzcon-1)-1)
    muce=nxzcon(2*nzcon)
    mlce=nxzcon(2*nzcon-1)
  else if (nxzrad(1) == m .and. nxzrad(last) /= 1) then
    iconra=4
    nzcon=nzrad
    nxzcon(2*nzcon-1)=nxzrad(2*nzcon)-1
    nxzcon(2*nzcon)=1
    if (nzcon > 1) then
      do n=1,nzcon-1
       nxzcon(2*n-1)=nxzrad(2*n)-1
       nxzcon(2*n)=nxzrad(2*n+1)+1
      enddo
    endif
    if (nzcon == 1) then
      muce=nxzcon(2)
      mlce=nxzcon(1)
    else
      muci=nxzcon(2*(nzcon-1))
      mlci=nxzcon(2*(nzcon-1)-1)
      muce=nxzcon(2*nzcon)
      mlce=nxzcon(2*nzcon-1)
    endif
  else if (nxzrad(1) /= m .and. nxzrad(last) == 1) then
    iconra=5
    nzcon=nzrad
    nxzcon(1)=m
    nxzcon(2)=nxzrad(1)+1
    if (nzcon > 1) then
      do n=2,nzcon
       nxzcon(2*n-1)=nxzrad(2*n-2)-1
       nxzcon(2*n)=nxzrad(2*n-1)+1
      enddo
    endif
    muci=nxzcon(2*nzcon)
    mlci=nxzcon(2*nzcon-1)
  endif

  return

end subroutine rechzco
!=======================================================================
subroutine konvek
!----------------------------------------------------------------------
  use strucmod,only: nr,konvv,konv,kzone,it2

  implicit none
!----------------------------------------------------------------------
  nr=nr+1

  if ((konvv-konv) < 0.d0) then

    kzone=kzone+1
    if (it2 == 6) then
      write(3,'(47x,a,1x,i1,a/)') 'BEGINNING OF THE',kzone,'th CONVECTIVE ZONE'
    endif
  else
    if ((konvv-konv) > 0.d0) then
      if (it2 == 6) then
        write(3,'(49x,a,1x,i1,a/)') 'end of the',kzone,'th convective zone'
      endif
    endif

  endif

  konvv=konv

  return

end subroutine konvek
!=======================================================================
subroutine CZdraw
!-----------------------------------------------------------------------
  implicit none

  integer:: ii
  real(kindreal):: convyn
!-----------------------------------------------------------------------
  drawcon(:)=0.d0
  iidraw=2
  convyn=1.d0

  if (zensi(m) < 0.d0) then
    iidraw=iidraw-1
    convyn=-1.d0
  endif
  do ii = m-1,1,-1
   if ((convyn > 0.d0 .and. zensi(ii) < 0.d0) .or. (convyn < 0.d0 .and. zensi(ii) > 0.d0)) then
     drawcon(iidraw)=1.d0-exp(q(ii))
     iidraw=iidraw+1
     convyn=-1.d0*convyn
   endif
   if (iidraw == 41) then
     if (verbose) then
       print*,'too many conv. zones for drawing'
     endif
     exit
   endif
  enddo
  if (mod(iidraw-1,2) == 1) then
    BaseZC = drawcon(iidraw-1)
  else
    BaseZC = -1.d0
  endif

  return

end subroutine CZdraw
!=======================================================================
end module convection
