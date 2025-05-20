module advection

use io_definitions
use evol,only: ldi,kindreal
use inputparam,only: iadvec,idebug,verbose
use rotmod,only: omegi
use safestop,only: safe_stop

implicit none

integer,parameter:: nminim=3
integer:: ncdiff

private
public :: advect,gbar,gtilgm

contains
!======================================================================
subroutine badv
! ---------------------------------------------------------------------
! Calcul des coefficients necessaires dans la methode de henyey
!   lors du calcul de l'evolution du moment cinetique ;
!   condition limite a la surface de l'etoile

! B_OMEGA est appele par HEN_OMEGA

! mtu est le numero de la premiere couche de la zone radiative
! condition au bord uniquement lorsque ibext=0
! ---------------------------------------------------------------------
use const, only: pi
use strucmod,only: rb,rho
use rotmod, only: ur,xldoex,vvomeg,vsuminenv,Flux_remaining
use diffadvmod,only: mtu,ibext,xbext,xcext
use equadiffmod,only: b1,b1o,b1u
use timestep,only: dzeit

implicit none

real(kindreal):: fact,xbextv,xcextv,meridflux_part,meridflux_tot
! ---------------------------------------------------------------------
  if (ibext == 0) then ! external layer is radiative, xbext and xcext come from the envelope
    xbextv = vsuminenv
    xcextv = vsuminenv*vvomeg(mtu)
  else ! external layer is convective, xbext and xcext were determined in confi
    xbextv = xbext
    xcextv = xcext
  endif

  fact = 4.d0*pi*2.d0/(3.d0*xbextv)
  meridflux_part = -0.2d0*exp(4.d0*rb(mtu)+rho(mtu))
  meridflux_tot = meridflux_part*omegi(mtu)*ur(mtu)
  b1=fact*(meridflux_tot+(Flux_remaining + xldoex)*3.d0/(8.d0*pi))*dzeit-(omegi(mtu)-xcextv/xbextv)
  b1o=ur(mtu)*fact*meridflux_part*dzeit-1.d0
  b1u=omegi(mtu)*fact*meridflux_part*dzeit

  return

end subroutine badv
!======================================================================
subroutine confi
!----------------------------------------------------------------------
! determine le type de configuration, ainsi que le type de conditions au bord a appliquer pour la convection
!    inoadv=1: pas d'application de l'advection
!    inoadv=0: application de l'advection dans une zone
!
! mtu numero de la couche superieure de la region avec advection
! npasr numero de la couche inferieure de la region avec advection
!
! ibext=0: condition bord ext ur=0
! ibext=1: condition bord ext chgmt moment cine zone convective
! ibint=0: condition bord int ur=0
! ibint=1: condition bord int chgmt moment cine zone convective
!-----------------------------------------------------------------------
use inputparam, only: phase
use strucmod,only: m
use convection,only: nzrad,nxzrad,nzcon,xmomin,xmocin,iconra
use diffadvmod,only: inoadv,ibext,ibint,mtu,npasr,xbint,xcint,xomint,xbext,xcext,xomext

implicit none

integer:: numer
!-----------------------------------------------------------------------
  inoadv=0
  xbext = 0.d0
  xcext = 0.d0
  xbint = 0.d0
  xcint = 0.d0

  if (idebug > 0) then
    write(*,*) 'iconra=',iconra
  endif

  select case(iconra)

    case(0)
! star is entirely convective
      inoadv = 1
    case(1)    ! star is entirely radiative
      ibext = 0
      ibint = 0
      mtu = 1
      npasr = m
    case(2)
! star is radiative at centre and surface, and some convective zones inbetween
      ibext = 0
      ibint = 1
      mtu = 1
      numer = 2*nzrad-1
      npasr = nxzrad(numer)
      xbint = xmomin(nzcon)
      xcint = xmocin(nzcon)
      xomint = xcint/xbint
    case(3)
! star is convective at centre and surface, and some other convective zones inbetween
      ibext = 1
      ibint = 1
      numer = 2*nzrad-1
      npasr = nxzrad(1)
      mtu = nxzrad(numer+1)
      xbint = xmomin(1)
      xbext = xmomin(nzcon)
      xcint = xmocin(1)
      xcext = xmocin(nzcon)
      xomint = xcint / xbint
      xomext = xcext / xbext
    case(4)
! star is radiative at centre and convective at surface
      if (nzrad == 1) then
        npasr = m
        mtu = nxzrad(2)
        ibint = 0
        ibext = 1
        xbext = xmomin(1)
        xcext = xmocin(1)
        xomext = xcext / xbext
      else
        numer = 2*nzrad-1
        npasr = nxzrad(numer)
        mtu = nxzrad(numer+1)
        ibext = 1
        ibint = 1
        xbext = xmomin(nzcon)
        xbint = xmomin(nzcon-1)
        xcext = xmocin(nzcon)
        xcint = xmocin(nzcon-1)
        xomext = xcext / xbext
        xomint = xcint / xbint
      endif
    case(5)
! star is convective at centre and radiative at surface
      ibext = 0
      ibint = 1
      numer = 1
      npasr = nxzrad(numer)
      mtu = 1
      xbint = xmomin(1)
      xcint = xmocin(1)
      xomint = xcint / xbint
    case default
      call safe_stop('Pb iconra > 5')
  end select

  ncdiff = mtu-npasr
! At the end of the MS of some models (high rot and/or low Z), some
! additional convectives zones appear above the core. This makes the
! computation of the advection difficult and not accurate. In this case,
! we don't apply the advection, and pass by the local angular momentum
! conservation.
!  if (nzrad > 1 .and. phase == 1) then
!    inoadv = 2
!    write(io_logs,*) 'Too many rad zones, advection not applied'
!    if (verbose) then
!      write(*,*) 'Too many rad zones, advection not applied'
!    endif
!  endif

  return

end subroutine confi
!======================================================================
subroutine coradv
!-----------------------------------------------------------------------
use rotmod,only: deladv,vvomeg
use convection,only: muci,mlci,muce,mlce
use diffadvmod,only: mtu,npasr,ibint,xomint,ibext,xomext

implicit none
!-----------------------------------------------------------------------
  deladv(mtu:npasr) = omegi(mtu:npasr)-vvomeg(mtu:npasr)

  if (ibint /= 0) then ! internal limit is convective
    deladv(muci:mlci) = omegi(npasr)-xomint
  endif
  if (ibext /= 0) then ! external limit is convective
    deladv(muce:mlce) = omegi(mtu)-xomext
  endif

  return

end subroutine coradv
!======================================================================
subroutine ebems
!-----------------------------------------------------------------------
use const,only: Lsol
use caramodele,only: glm,gls
use abundmod,only: eps,epsy,epsc
use strucmod,only: q,sb,m
use rotmod,only: ebem

implicit none

integer::n
real(kindreal)::zwi1,epsbar,epsmoy
!-----------------------------------------------------------------------
  zwi1 = 1.d0 / (exp(sb(1))-1.d0)

  do n=1,m-1
   epsbar = eps(n) + epsy(n) + epsc(n)
   epsmoy = ((exp(sb(n))-1.d0) * zwi1 * gls*Lsol) / (exp(glm)*(1.d0-exp(q(n))))
   ebem(n) = epsbar / epsmoy
  enddo

  ebem(m) = 1.d0

  return

end subroutine ebems
!======================================================================
subroutine echtem
!-----------------------------------------------------------------------
use const,only: Msol
use caramodele,only: glm
use strucmod,only: m,tb,rb,q,amub
use rotmod,only: ht,xmu

implicit none

integer::i
real(kindreal)::yyim,yyi,yyip,xxim,xxi,xxip,wpena,wpenb,wfa,wfb,xmso
!-----------------------------------------------------------------------
  do i=2,m-1
   yyim = tb(i-1)
   yyi = tb(i)
   yyip = tb(i+1)
   xxim = rb(i-1)
   xxi = rb(i)
   xxip = rb(i+1)
   wpena = (yyi-yyip) / (xxi-xxip)
   wpenb = (yyim-yyi) / (xxim-xxi)
   wfa= (xxi-xxim) / (xxip-xxim)
   wfb= (xxip-xxi) / (xxip-xxim)
   ht(i) = abs(exp(xxi) / (wpena*wfa+wpenb*wfb))
   xmso = exp(glm) * (1.d0-exp(q(i))) / Msol
  enddo
  ht(m) = ht(m-1)
  ht(1) = ht(2)

  xmu(1:m) = amub(1:m)

  return

end subroutine echtem
!======================================================================
subroutine gbar
!-----------------------------------------------------------------------
! calcul de la valeur moyenne de la gravite sur une isobare
!-----------------------------------------------------------------------
use const,only: cst_G
use caramodele ,only: glm
use equadiffmod,only: ccg1
use strucmod,only: m,rb,q
use rotmod,only: xmeg
use geomod, only: rpsi_min,rpsi_max,geomeang

implicit none

integer:: n
real(kindreal):: xpsin,xmegn
!-----------------------------------------------------------------------
  do n=1,m-1
   xpsin = exp(rb(n)) / ((cst_G*exp(ccg1)*(1.d0-exp(q(n)))) / (omegi(n)*omegi(n)))**(1.d0/3.d0)
   if (xpsin >= rpsi_min) then
       if (xpsin > rpsi_max) then
         xpsin = 0.9999999999d0*rpsi_max
       endif
     call geomeang(xpsin,xmegn)
     xmeg(n) = (cst_G * (1.d0-exp(q(n))) * exp(glm) * omegi(n)*omegi(n)*omegi(n)*omegi(n))**(1.d0/3.d0)*xmegn
   else
     xmeg(n) = cst_G * (1.d0-exp(q(n)))*exp(glm) / (exp(rb(n))*exp(rb(n)))
   endif
  enddo

  xmeg(m) = 0.d0

  return

end subroutine gbar
!======================================================================
subroutine giadv(j,j1)
!-----------------------------------------------------------------------
! calcule les elements de matrice correspondant aux
! 4 equations permettant de resoudre la partie advective
! de l'equation de transport du moment angulaire.
!-----------------------------------------------------------------------
use const,only: pi,cst_G
use caramodele,only: glm
use strucmod,only: rho,q,rb,xb,delt,opac,opact,epsit,H_P,Nabla_mu
use rotmod,only: xmeg,theta,ht,ur,xmu,aux,ur1,ebem,gtilde,vomegi
use diffadvmod,only: D_eff,xueff,D_h,K_ther
use equadiffmod,only: ag1,ag1t,ag1t1,ag1a,ag1a1,ag1u,ag1u1,ag1o,ag1o1,ag2,ag2t,ag2t1,ag2a,ag2a1,ag2u,ag2u1,ag2o,ag2o1,&
                      ag3,ag3t,ag3t1,ag3a,ag3a1,ag3u,ag3u1,ag3o,ag3o1,ag4,ag4t,ag4t1,ag4a,ag4a1,ag4u,ag4u1,ag4o,ag4o1
use timestep,only: dzeit

implicit none

integer,intent(in):: j,j1

real(kindreal):: xurdt1,xure,aux3t1,aurj1,anitj,anitj1,aurj,auro1,auro,aux1t,aux1j2,aux3t,aux2j2,aux1t1,aux3j2,&
  aux5t1,aux5t,aux5j2,aux4j2,bur4t1,bb22,bb11,bibo,bur3j1,bur2a,bur1,bur2,bur1j,bur1j1,bur3,bur2a1,bur3j,&
  bur4t,bur4,bur3t,bur3t1,cc11,cap_gi,ca,cap1_gi,capt_gi,capt1_gi,cb,cc33,cc22,cq,xurbt,tlno,cur3j,cur2,cur1,&
  cur2a,cur3,cur2a1,cur4t,cur4,cur3t,cur3j1,cur3t1,cur5,cur4t1,cur5j,cur5j1,epst_gi,epst1_gi,fac,kmuj1,kmuj,&
  qlno,qlnd,xktj1,xeffj1,qur,rhomj1,rh1,rh,rhomj,rht1,rht,ur1o,thedj,thedj1,tqrmr,tlnoj,tlnoj1,ur1o1,xdq,xdodq,&
  xdom,xeffj,xktj,xlamj1,xlamj,xuret1,xurej1,xurca1,xurb,xthe,xur,xuetj1,xuetj,xueuj,xueuj1,xura,xurao1,xurao,&
  xurba,xurba1,xurc,xurbt1,xurca,xurej,xurdo1,xurd,xurct1,xurct,xurdo,xurdt,xuret,ylamj,xureu1,xureu,ccc
!-----------------------------------------------------------------------
  rh=rho(j)
  rh1=rho(j1)
  rht=-delt(j)
  rht1=-delt(j1)
  cap_gi=opac(j)
  cap1_gi=opac(j1)
  capt_gi=opact(j)
  capt1_gi=opact(j1)
  epst_gi=epsit(j)
  epst1_gi=epsit(j1)

! equation (1) definition de teta
  xdq=q(j1)-q(j)
  xdodq=(log(omegi(j1))-log(omegi(j)))/xdq
  ca=exp(glm+0.5d0*(q(j)+q(j1))-log(8.d0*pi/3.d0)-2.d0*(rb(j)+rb(j1))-0.5d0*(rh+rh1))*sqrt(xmeg(j)*xmeg(j1))
  xdom=exp(-(log(omegi(j))+log(omegi(j1))))
  xthe=(theta(j)+theta(j1))/2.d0
  ag1=xdodq/(xdom)+ca*xthe
  ag1t=0.5d0*ca
  ag1t1=0.5d0*ca
  ag1a=0.d0
  ag1a1=0.d0
  ag1u=0.d0
  ag1u1=0.d0
  ag1o=1.d0/(xdom*omegi(j))*(xdodq-1.d0/xdq)
  ag1o1=1.d0/(xdom*omegi(j1))*(xdodq+1.d0/xdq)

  ccc=100.0d0
  ag1=ag1*ccc
  ag1t=ag1t*ccc
  ag1t1=ag1t1*ccc
  ag1a=0.d0
  ag1a1=0.d0
  ag1u=0.d0
  ag1u1=0.d0
  ag1o=ag1o*ccc
  ag1o1=ag1o1*ccc

! equation (2) definition de aux
! cb=M/(4pi r2 rho)
  cb=exp(glm-log(4.d0*pi)-(rb(j)+rb(j1))-0.5d0*(rh+rh1))
  cq=exp(0.5d0*(q(j)+q(j1)))
! tqrmr=-e^Q * M/(4pi r2 rho)
  tqrmr=-cb*cq

  thedj=-theta(j)/rht
  thedj1=-theta(j1)/rht1
  aux1j2=-(ht(j)+ht(j1))/2.d0*(thedj1-thedj)/xdq
  aux1j2=aux1j2/tqrmr
  aux1t=-(ht(j)+ht(j1))/2.d0*(-1.d0)/(xdq*tqrmr)*(-1.d0/rht)
  aux1t1=-(ht(j)+ht(j1))/2.d0*(1.d0)/(xdq*tqrmr)*(-1.d0/rht1)

  if (ur(j) /= 0.0d0) then
    xlamj=exp(rb(j))/(60.d0*H_P(j))*Nabla_mu(j)*ur(j)/abs(ur(j))
  else
    xlamj=0.0d0
  endif
  if (ur(j1) /= 0.0d0) then
    xlamj1=exp(rb(j1))/(60.d0*H_P(j1))*Nabla_mu(j1)*ur(j1)/abs(ur(j1))
  else
    xlamj1=0.0d0
  endif

  if (D_eff(j) /= 0.d0 .and. ur(j) /= 0.0d0) then
    ylamj=xlamj*10.d0*6.d0/exp(rb(j))*5.d0*D_eff(j)/ur(j)
  endif

  xlamj=0.0d0
  xlamj1=0.0d0
  aux2j2=(ht(j)+ht(j1))/2.d0*(xlamj1/(-rht1)-xlamj/(-rht))/xdq
  aux2j2=aux2j2/tqrmr

  aux3j2=0.5d0*(thedj+thedj1)
  aux3t=0.5d0*(-1.d0)/rht
  aux3t1=0.5d0*(-1.d0)/rht1

  kmuj=0.16d0/(exp(cap_gi)*xmu(j))-1.d0
  kmuj1=0.16d0/(exp(cap1_gi)*xmu(j1))-1.d0
  xktj=3.d0-capt_gi-rht
  xktj1=3.d0-capt1_gi-rht1
  cc11=0.25d0*(kmuj+kmuj1)*(xlamj+xlamj1)
  cc22=0.25d0*(xktj/(-rht)+xktj1/(-rht1))*(xlamj+xlamj1)
  cc33=0.25d0*(1.d0/(-rht)+1.d0/(-rht1))*(xlamj+xlamj1)
  aux4j2=-cc11-cc22-cc33

  bibo   = 0.5d0*(xktj+xktj1)-1.d0
  aux5j2 = bibo*0.5d0*(theta(j)+theta(j1))
  aux5t  = bibo*0.5d0
  aux5t1 = bibo*0.5d0

  bb11=0.5d0*(aux(j)+aux(j1))
  bb22=aux1j2+aux2j2+aux3j2+aux4j2+aux5j2
  ag2=bb11-bb22
  ag2t=-(aux1t+aux3t+aux5t)
  ag2t1=-(aux1t1+aux3t1+aux5t1)
  ag2a=0.5d0
  ag2a1=0.5d0
  ag2u=0.d0
  ag2u1=0.d0
  ag2o=0.d0
  ag2o1=0.d0

  ccc=1.0d-2
  ag2=ag2*ccc
  ag2t=ag2t*ccc
  ag2t1=ag2t1*ccc
  ag2a=ag2a*ccc
  ag2a1=ag2a1*ccc
  ag2u=0.d0
  ag2u1=0.d0
  ag2o=0.d0
  ag2o1=0.d0

! equation (3) definition de ur
  rhomj=exp(glm)*(1.d0-exp(q(j))) /(4.d0*pi/3.d0*exp(3.d0*rb(j)))
  rhomj1=exp(glm)*(1.d0-exp(q(j1)))/(4.d0*pi/3.d0*exp(3.d0*rb(j1)))
  ur1o=ur1(j)/(pi*cst_G*rhomj/omegi(j)-omegi(j)/2.d0)
  ur1o1=ur1(j1)/(pi*cst_G*rhomj1/omegi(j1)-omegi(j1)/2.d0)

  aurj=2.d0*(1.d0-(omegi(j)*omegi(j))/(2.d0*pi*cst_G*exp(rh))-ebem(j))*gtilde(j)
  auro=-4.d0*gtilde(j)*omegi(j)/(2.d0*pi*cst_G*exp(rh))+aurj*2.d0/omegi(j)
  aurj1=2.d0*(1.d0-(omegi(j1)*omegi(j1))/(2.d0*pi*cst_G*exp(rh1))-ebem(j1))*gtilde(j1)
  auro1=-4.d0*gtilde(j1)*omegi(j1)/(2.d0*pi*cst_G*exp(rh1))+aurj1*2.d0/omegi(j1)
  xura=0.5d0*(aurj+aurj1)
  xurao=0.5d0*(auro)
  xurao1=0.5d0*(auro1)

  bur1j=rhomj/exp(rh)
  bur1j1=rhomj1/exp(rh1)
  bur1=0.5d0*(bur1j+bur1j1)
  bur2=exp(0.5d0*(rb(j)+rb(j1)))/3.d0*(aux(j1)-aux(j))/(q(j1)-q(j))*1.d0/tqrmr
  bur2a=exp(0.5d0*(rb(j)+rb(j1)))/3.d0*(-1.d0)/(q(j1)-q(j))*1.d0/tqrmr
  bur2a1=exp(0.5d0*(rb(j)+rb(j1)))/3.d0*1.d0/(q(j1)-q(j))*1.d0/tqrmr
  bur3j=2.d0*ht(j)/exp(rb(j))*(theta(j)/(-rht)-xlamj/(-rht))
  bur3t=2.d0*ht(j)/(exp(rb(j))*(-rht))
  bur3j1=2.d0*ht(j1)/exp(rb(j1))*(theta(j1)/(-rht1)-xlamj1/(-rht1))
  bur3t1=2.d0*ht(j1)/(exp(rb(j1))*(-rht1))
  bur3=0.5d0*(bur3j+bur3j1)
  bur4=-2.d0/3.d0*0.5d0*(theta(j)+theta(j1))
  bur4t=-2.d0/3.d0*0.5d0
  bur4t1=-2.d0/3.d0*0.5d0
  xurb=bur1*(bur2+bur3+bur4)
  xurbt=bur1*(0.5d0*bur3t+bur4t)
  xurbt1=bur1*(0.5d0*bur3t1+bur4t1)
  xurba=bur1*(bur2a)
  xurba1=bur1*(bur2a1)

  cur1=0.5d0*(ebem(j)+ebem(j1))
  cur2=0.5d0*(aux(j)+aux(j1))
  cur2a=0.5d0
  cur2a1=0.5d0
  cur3j=-theta(j)*(epst_gi/(-rht)-xktj/(-rht))
  cur3t=-(epst_gi/(-rht)-xktj/(-rht))
  cur3j1=-theta(j1)*(epst1_gi/(-rht1)-xktj1/(-rht1))
  cur3t1=-(epst1_gi/(-rht1)-xktj1/(-rht1))
  cur3=0.5d0*(cur3j+cur3j1)
  cur4=-0.5d0*(xktj*theta(j)+xktj1*theta(j1))
  cur4t=-0.5d0*xktj
  cur4t1=-0.5d0*xktj1
! attention epsi mu correct seulement pour H burn
  cur5j=xlamj*(-4.d0/(5.d0*xb(j)*xmu(j))+epst_gi/(-rht))
  cur5j1=xlamj1*(-4.d0/(5.d0*xb(j1)*xmu(j1))+epst1_gi/(-rht1))
  cur5=0.5d0*(cur5j+cur5j1)
  xurc=cur1*(cur2+cur3+cur4+cur5)
  xurct=cur1*(0.5d0*cur3t+cur4t)
  xurct1=cur1*(0.5d0*cur3t1+cur4t1)
  xurca=cur1*(cur2a)
  xurca1=cur1*(cur2a1)

  xurd=0.0d0
  xurdt=0.0d0
  xurdt1=0.0d0
  xurdo=0.0d0
  xurdo1=0.0d0

  anitj  = 2.d0*ht(j)*bur1j
  anitj1 = 2.d0*ht(j1)*bur1j1
  xurej  = anitj/exp(rb(j))*D_h(j)/K_ther(j)*thedj
  xurej1 = anitj1/exp(rb(j1))*D_h(j1)/K_ther(j1)*thedj1
  xure   = 0.5d0*(xurej+xurej1)
  xuetj  =-anitj/exp(rb(j))*D_h(j)/(K_ther(j)*rht)
  xuetj1 =-anitj1/exp(rb(j1))*D_h(j1)/(K_ther(j1)*rht1)
  xueuj  = 0.d0
  xueuj1 = 0.d0
  xuret  = 0.5d0*xuetj
  xuret1 = 0.5d0*xuetj1
  xureu  = 0.5d0*xueuj
  xureu1 = 0.5d0*xueuj1

  ccc=1.0d-4
  xur=0.5d0*(ur1(j)+ur1(j1))*(xura+xurb+xurc+xurd+xure)*xueff
  ag3=0.5d0*(ur(j)+ur(j1))-xur
  ag3t=-0.5d0*(ur1(j)+ur1(j1))*(xurbt+xurct+xurdt+xuret)*xueff
  ag3t1=-0.5d0*(ur1(j)+ur1(j1))*(xurbt1+xurct1+xurdt1+xuret1)*xueff
  ag3a=-0.5d0*(ur1(j)+ur1(j1))*(xurba+xurca)*xueff
  ag3a1=-0.5d0*(ur1(j)+ur1(j1))*(xurba1+xurca1)*xueff
  ag3u=0.5d0-0.5d0*(ur1(j)+ur1(j1))*xureu*xueff
  ag3u1=0.5d0-0.5d0*(ur1(j)+ur1(j1))*xureu1*xueff
  ag3o=-(0.5d0*ur1o*(xura+xurb+xurc+xurd+xure)+0.5d0*(ur1(j)+ur1(j1))*(xurao+xurdo))*xueff
  ag3o1=-(0.5d0*ur1o1*(xura+xurb+xurc+xurd+xure)+0.5d0*(ur1(j)+ur1(j1))*(xurao1+xurdo1))*xueff
  ag3=ccc*ag3
  ag3t=ccc*ag3t
  ag3t1=ccc*ag3t1
  ag3a=ccc*ag3a
  ag3a1=ccc*ag3a1
  ag3u=ccc*ag3u
  ag3u1=ccc*ag3u1
  ag3o=ccc*ag3o
  ag3o1=ccc*ag3o1

! equation (4) equation differentielle pour omega
  xeffj=1.d0
  xeffj1=1.d0
  tlnoj=(log(omegi(j))-log(vomegi(j)))/dzeit
  tlnoj1=(log(omegi(j1))-log(vomegi(j1)))/dzeit
  tlno=0.5d0*(tlnoj+tlnoj1)
  qlnd=(rh1-rh)/xdq
  qlno=(log(omegi(j1))-log(omegi(j)))/xdq
  qur=(xeffj1*ur(j1)-xeffj*ur(j))/xdq

  fac=4.d0*pi/5.d0*exp(rb(j)+rb(j1)+0.5d0*(rh+rh1))/(exp(0.5d0*(q(j)+q(j1)))*exp(glm))

  ag4=tlno+fac*(0.5d0*(xeffj*ur(j)+xeffj1*ur(j1))*(qlnd+qlno)+qur)-2.d0/5.d0*(xeffj*ur(j)+xeffj1*ur(j1))/exp(0.5d0*(rb(j)+rb(j1)))
  ag4t=0.d0
  ag4t1=0.d0
  ag4a=0.d0
  ag4a1=0.d0
  ag4u=fac*(0.5d0*xeffj*(qlnd+qlno)-xeffj/xdq)-2.d0*xeffj/(5.d0*exp(0.5d0*(rb(j)+rb(j1))))
  ag4u1=fac*(0.5d0*xeffj1*(qlnd+qlno)+xeffj1/xdq)-2.d0*xeffj1/(5.d0*exp(0.5d0*(rb(j)+rb(j1))))
  ag4o=0.5d0*(1.d0/dzeit*1.d0/omegi(j))+fac*(0.5d0*(xeffj*ur(j)+xeffj1*ur(j1))*(-1.d0/xdq*1.d0/omegi(j)))
  ag4o1=0.5d0*(1.d0/dzeit*1.d0/omegi(j1))+fac*(0.5d0*(xeffj*ur(j)+xeffj1*ur(j1))*(1.d0/xdq*1.d0/omegi(j1)))
  ccc=1.d+5
  ag4=ag4*ccc
  ag4t=0.d0
  ag4t1=0.d0
  ag4a=0.d0
  ag4a1=0.d0
  ag4u=ag4u*ccc
  ag4u1=ag4u1*ccc
  ag4o=ag4o*ccc
  ag4o1=ag4o1*ccc

  return

end subroutine giadv
!======================================================================
subroutine gtilgm
!-----------------------------------------------------------------------
!calcul de g_tilde/g_moyen: fluctuation de g sur l'isobare par rapport au g moyen de la couche
!-----------------------------------------------------------------------
use const,only: cst_G,pi
use caramodele ,only: glm
use strucmod,only: q,rb,rho,m
use rotmod,only: gtilde

implicit none

integer:: n
real(kindreal):: acentri,gravm
!-----------------------------------------------------------------------
  do n=1,m-1
   acentri = omegi(n)*omegi(n) * exp(rb(n))
   gravm = cst_G*exp(glm)*(1.d0-exp(q(n))) / exp(2.d0*rb(n))
   gtilde(n) = 4.d0/3.d0 * acentri / gravm
  enddo

  gtilde(m) = omegi(m)*omegi(m) / (cst_G*pi*exp(rho(m)))

  return

end subroutine gtilgm
!======================================================================
subroutine henadv(alph1,flag_girl)
!-----------------------------------------------------------------------
! Resolution des equations aux derivees partielles du premier ordre
!   decrivant l'advection du moment cinetique par la methode de relaxation
!   de Henyey.

! HENADV est appelee par ADVEC
!                et appelle B_OMEGA, GI_OMEGA, Z_OMEGA et GIRL

! Derniere version : 13 novembre 1996
!-----------------------------------------------------------------------
use caramodele,only: nwmd
use rotmod,only: theta,ur,aux
use strucmod,only: m
use diffadvmod,only: mtu,npasr
use equadiffmod,only:ag1,ag1t,ag1t1,ag1a,ag1a1,ag1u,ag1u1,ag1o,ag1o1,ag2,ag2t,ag2t1,ag2a,ag2a1,ag2u,ag2u1,ag2o,ag2o1,ag3,ag3t, &
  ag3t1,ag3a,ag3a1,ag3u,ag3u1,ag3o,ag3o1,ag4,ag4t,ag4t1,ag4a,ag4a1,ag4u,ag4u1,ag4o,ag4o1,jterma,b1O,b1U,b1,az1O1,az1U1,az1
use SmallFunc,only: girl

implicit none

integer,intent(inout):: flag_girl
real(kindreal),intent(inout):: alph1

integer::j,j1,jgg1,jgg2,jgg3,jgg4,mtu1,jgdu,jgda,jdga,jgdo,jgdt,iterad,iSE,jSE
real(kindreal), parameter:: itmax=40
real(kindreal):: agdu,agdth,agdo,agda,agg,agmax,gg1,gg2,gg3,gg4,dur,dteta,dom1,daux1,dur1,gdu,gda,gdo,gdt,dom,daux,gdutamp, &
  gdttamp,gdotamp,gdatamp,dteta1
real(kindreal), dimension(5):: zu
real(kindreal), dimension(ldi):: uu,vu,wu,ua,va,wa,uo,vo,wo,ut,vt,wt
real(kindreal), dimension(4,3):: hu
real(kindreal), dimension(4,7):: ha
real(kindreal), dimension(5,3):: u
real(kindreal), dimension(5,6):: za
real(kindreal), dimension(5,8):: a

logical:: endIter
!-----------------------------------------------------------------------
  if (idebug > 0) then
    write(*,*) 'HENADV: mtu,npasr,m',mtu,npasr,m
  endif
  agdu = 0.1d-7
  agdth = agdu
  agdo = agdu
  agda = agdu
  agg = 1.d-4
  agmax = 1.d6
  endIter = .false.

  alph1 = 1.d0/8.d0
  if (mtu /= 1) then
    theta(mtu) = 0.d0
  else
! thext1 condition au bord avec correction de omega(1)
! thext2 condition au bord sans correction de omega(1)
!    theta(mtu) = thext2
    theta(mtu) = 0.d0
  endif
  theta(npasr) = 0.d0

  iterad = 0
  jterma = 0
  do while (.not. endIter)
   flag_girl=0
   iterad = iterad+1 ! iterad est le compteur d'iterations effectuees.
   if (iterad  <=  3) then
     alph1 = 2.d0*alph1
   else
     alph1 = 0.8d0
   endif

! Calcul du bloc 1. Debut de la zone radiative.
! Condition limite a la surface (flux de moment cinetique).
!-----------------------------------------------------------
   if (idebug > 0) then
     write(*,*) 'call badv'
   endif
   call badv

   a(1,1) = b1o
   a(1,2) = 0.d0
   a(1,3) = b1u
   a(1,4) = 0.d0
   a(1,5) = 0.d0
   a(1,6) = 0.d0
   a(1,7) = 0.d0
   a(1,8) = -b1

! Calcul de Gi, dGi/domj, dGi/dauxj, dGi/durj, dGi/dtetaj.
!----------------------------------------------------------
! i=1,2,3,4, j=mtu,mtu+1.
   j = mtu
   j1 = mtu+1
   if (idebug > 0) then
     write(*,*) 'call giadv'
   endif
   call giadv(j,j1)

   gg1 = ag1     ! erreur maximale
   gg2 = ag2
   gg3 = ag3
   gg4 = ag4
   jgg1 = mtu   ! # de la couche ou se trouve l'erreur maximale
   jgg2 = mtu
   jgg3 = mtu
   jgg4 = mtu

   a(2,1) = ag1o
   a(2,2) = ag1a
   a(2,3) = ag1u
   a(2,4) = ag1o1
   a(2,5) = ag1a1
   a(2,6) = -ag1u1
   a(2,7) = -ag1t1
   a(2,8) = -ag1
   a(3,1) = ag2o
   a(3,2) = ag2a
   a(3,3) = ag2u
   a(3,4) = ag2o1
   a(3,5) = ag2a1
   a(3,6) = -ag2u1
   a(3,7) = -ag2t1
   a(3,8) = -ag2
   a(4,1) = ag3o
   a(4,2) = ag3a
   a(4,3) = ag3u
   a(4,4) = ag3o1
   a(4,5) = ag3a1
   a(4,6) = -ag3u1
   a(4,7) = -ag3t1
   a(4,8) = -ag3
   a(5,1) = ag4o
   a(5,2) = ag4a
   a(5,3) = ag4u
   a(5,4) = ag4o1
   a(5,5) = ag4a1
   a(5,6) = -ag4u1
   a(5,7) = -ag4t1
   a(5,8) = -ag4

! On exprime les variables par combinaisons lineaires de dur2
! et dteta2 :
!             dom1=uo(1)*dur2+vo(1)*dteta2+wo(1)
!             ...
!             daux2=ua(2)*dur2+va(2)*dteta2+wa(2)

   if (idebug > 0) then
     write(*,*) 'call girl (5,3)'
   endif
   flag_girl = 0
   call girl(a,u,5,3,flag_girl)
   if (flag_girl /= 0) then
     if (idebug>0) then
       write(*,*) 'henadv - matrix a(5,8),flag:',flag_girl
       do iSE=1,5
        do jSE=1,8
         write(*,'("a(",i1,",",i1,") :",d22.12)') iSE,jSE,a(iSE,jSE)
       enddo
      enddo
     endif
     rewind(io_runfile)
     write(io_runfile,*) nwmd,':girl crashes in henadv with matrix a(5,8)'
     call safe_stop('girl crashes in henadv with matrix a(5,8)')
   endif

! Stockage des coefficients
   mtu1 = mtu+1
   uo(mtu) = u(1,1)
   vo(mtu) = u(1,2)
   wo(mtu) = u(1,3)
   ua(mtu) = u(2,1)
   va(mtu) = u(2,2)
   wa(mtu) = u(2,3)
   uu(mtu) = u(3,1)
   vu(mtu) = u(3,2)
   wu(mtu) = u(3,3)
   uo(mtu1) = u(4,1)
   vo(mtu1) = u(4,2)
   wo(mtu1) = u(4,3)
   ua(mtu1) = u(5,1)
   va(mtu1) = u(5,2)
   wa(mtu1) = u(5,3)

! Calcul du bloc j. Couches mtu+1 a npasr-2 de l'interieur.
!---------------------------------------------------------------
! Les variables sont exprimees sous forme de combinaisons lineaires
! de durj et dtetaj

   if (idebug > 0) then
     write(*,*) 'call giadv'
   endif
   do
    j = j+1     ! Dernier appel pour j=npasr-1
    j1 = j+1
    call giadv(j,j1)

    if (abs(gg1) < abs(ag1)) then
      gg1 = ag1
      jgg1 = j
    endif
    if (abs(gg2) < abs(ag2)) then
      gg2 = ag2
      jgg2 = j
    endif
    if (abs(gg3) < abs(ag3)) then
      gg3 = ag3
      jgg3 = j
    endif
    if (abs(gg4) < abs(ag4)) then
      gg4 = ag4
      jgg4 = j
    endif
    if (j1 >= npasr .or. j1 >= m) then
      exit
    endif
! Si j1 < npasr, on passe le calcul de la couche noralement.
! Si j1 = npasr, on passe au dernier bloc, limite au centre.

! Remplissage du tableau ha, qui contient les elements de la matrice
!   a inverser  (p.152)
    ha(1,1) = ag1u+uo(j)*ag1o+ua(j)*ag1a
    ha(1,2) = ag1t+vo(j)*ag1o+va(j)*ag1a
    ha(1,3) = ag1o1
    ha(1,4) = ag1a1
    ha(1,5) = -ag1u1
    ha(1,6) = -ag1t1
    ha(1,7) = -ag1-wo(j)*ag1o-wa(j)*ag1a
    ha(2,1) = ag2u+uo(j)*ag2o+ua(j)*ag2a
    ha(2,2) = ag2t+vo(j)*ag2o+va(j)*ag2a
    ha(2,3) = ag2o1
    ha(2,4) = ag2a1
    ha(2,5) = -ag2u1
    ha(2,6) = -ag2t1
    ha(2,7) = -ag2-wo(j)*ag2o-wa(j)*ag2a
    ha(3,1) = ag3u+uo(j)*ag3o+ua(j)*ag3a
    ha(3,2) = ag3t+vo(j)*ag3o+va(j)*ag3a
    ha(3,3) = ag3o1
    ha(3,4) = ag3a1
    ha(3,5) = -ag3u1
    ha(3,6) = -ag3t1
    ha(3,7) = -ag3-wo(j)*ag3o-wa(j)*ag3a
    ha(4,1) = ag4u+uo(j)*ag4o+ua(j)*ag4a
    ha(4,2) = ag4t+vo(j)*ag4o+va(j)*ag4a
    ha(4,3) = ag4o1
    ha(4,4) = ag4a1
    ha(4,5) = -ag4u1
    ha(4,6) = -ag4t1
    ha(4,7) = -ag4-wo(j)*ag4o-wa(j)*ag4a

! Calcul de (Ui,Vi,Wi) par inversion
! Le tableau hu(m,n) contient les coefficients U,V,W

    if (idebug > 0) then
      write(*,*) 'call girl (4,3)'
    endif
    flag_girl = 0
    call girl(ha,hu,4,3,flag_girl)
    if (flag_girl /= 0) then
      if (idebug>0) then
        write(*,*) 'henadv - matrix ha(4,7),flag:',flag_girl
        do iSE=1,4
         do jSE=1,7
          write(*,'("ha(",i1,",",i1,") :",d22.12)') iSE,jSE,ha(iSE,jSE)
        enddo
       enddo
      endif
      rewind(io_runfile)
      write(io_runfile,*) nwmd,':girl crashes in henadv with matrix ha(4,7)'
      call safe_stop('girl crashes in henadv with matrix ha(4,7)')
    endif


    uu(j) = hu(1,1)
    vu(j) = hu(1,2)
    wu(j) = hu(1,3)
    ut(j) = hu(2,1)
    vt(j) = hu(2,2)
    wt(j) = hu(2,3)
    uo(j1) = hu(3,1)
    vo(j1) = hu(3,2)
    wo(j1) = hu(3,3)
    ua(j1) = hu(4,1)
    va(j1) = hu(4,2)
    wa(j1) = hu(4,3)
   enddo

! Remplissage du dernier bloc j=npasr-1, j1=npasr.
!--------------------------------------------------
! On a cinq equations (4Gi + condition limite au centre)
! Les Gi sont deja calculees

   za(1,1) = ag1u+uo(j)*ag1o+ua(j)*ag1a
   za(1,2) = ag1t+vo(j)*ag1o+va(j)*ag1a
   za(1,3) = ag1o1
   za(1,4) = ag1a1
   za(1,5) = ag1u1
   za(1,6) = -ag1-wo(j)*ag1o-wa(j)*ag1a
   za(2,1) = ag2u+uo(j)*ag2o+ua(j)*ag2a
   za(2,2) = ag2t+vo(j)*ag2o+va(j)*ag2a
   za(2,3) = ag2o1
   za(2,4) = ag2a1
   za(2,5) = ag2u1
   za(2,6) = -ag2-wo(j)*ag2o-wa(j)*ag2a
   za(3,1) = ag3u+uo(j)*ag3o+ua(j)*ag3a
   za(3,2) = ag3t+vo(j)*ag3o+va(j)*ag3a
   za(3,3) = ag3o1
   za(3,4) = ag3a1
   za(3,5) = ag3u1
   za(3,6) = -ag3-wo(j)*ag3o-wa(j)*ag3a
   za(4,1) = ag4u+uo(j)*ag4o+ua(j)*ag4a
   za(4,2) = ag4t+vo(j)*ag4o+va(j)*ag4a
   za(4,3) = ag4o1
   za(4,4) = ag4a1
   za(4,5) = ag4u1
   za(4,6) = -ag4-wo(j)*ag4o-wa(j)*ag4a

   if (idebug > 0) then
     write(*,*) 'call ziadv'
   endif
   call ziadv

   za(5,1) = 0.d0
   za(5,2) = 0.d0
   za(5,3) = az1o1
   za(5,4) = 0.d0
   za(5,5) = az1u1
   za(5,6) = -az1

   if (idebug > 0) then
     write(*,*) 'call girl (5,1)'
   endif
   flag_girl = 0
   call girl(za,zu,5,1,flag_girl)
   if (flag_girl /= 0) then
     if (idebug>0) then
       write(*,*) 'henadv - matrix za(5,6),flag:',flag_girl
       do iSE=1,5
        do jSE=1,6
         write(*,'("za(",i1,",",i1,") :",d22.12)') iSE,jSE,za(iSE,jSE)
       enddo
      enddo
     endif
     rewind(io_runfile)
     write(io_runfile,*) nwmd,':girl crashes in henadv with matrix za(5,6)'
     call safe_stop('girl crashes in henadv with matrix za(5,6)')
   endif

   dur = zu(1)      ! npasr-1
   dteta = zu(2)
   dom1 = zu(3)     ! npasr
   daux1 = zu(4)
   dur1 = zu(5)
! Introduction des variables qui notent la correction maximale
   if (abs(ur(npasr)) > 1.d-15) then
     gdu = dur1/ur(npasr)
     jgdu = npasr
   else
     gdu = 0.d0
     jgdu = npasr
   endif
   if (abs(aux(npasr)) > 1.d-15) then
     gda = daux1/aux(npasr)
     jgda = npasr
   else
     gda = 0.d0
     jdga = npasr
   endif
   if (abs(omegi(npasr)) > 1.d-15) then
     gdo = dom1/omegi(npasr)
     jgdo = npasr
   else
     gdo = 0.d0
     jgdo = npasr
   endif
   if (abs(theta(npasr-1)) > 1.d-15) then
     gdt = dteta/theta(npasr-1)
     jgdt = npasr-1
   else
     gdt = 0.d0
     jgdt = npasr-1
   endif

   omegi(npasr) = abs(omegi(npasr)+alph1*dom1)
   aux(npasr) = aux(npasr)+alph1*daux1
   ur(npasr) = ur(npasr)+alph1*dur1

   dom = uo(j)*dur+vo(j)*dteta+wo(j)
   daux = ua(j)*dur+va(j)*dteta+wa(j)

! Remontee des corrections
!--------------------------
   do j=npasr-1,mtu,-1

! Recherche de l'erreur relative maximale
    if (abs(ur(j)) > 1.d-15) then
      gdutamp = dur/ur(j)
      if (abs(gdutamp) > abs(gdu)) then
        gdu = gdutamp
        jgdu = j
      endif
    endif
    if (abs(theta(j)) > 1.d-15) then
      gdttamp = dteta/theta(j)
      if (abs(gdttamp) > abs(gdt)) then
        gdt = gdttamp
        jgdt = j
      endif
    endif
    if (abs(omegi(j)) > 1.d-15) then
      gdotamp = dom/omegi(j)
      if (abs(gdotamp) > abs(gdo)) then
        gdo = gdotamp
        jgdo = j
      endif
    endif
    if (abs(aux(j)) > 1.d-15) then
      gdatamp = daux/aux(j)
      if (abs(gdatamp) > abs(gda)) then
        gda = gdatamp
        jgda = j
      endif
    endif

! Correction des valeurs
    omegi(j) = abs(omegi(j) + alph1*dom)
    aux(j) = aux(j) + alph1*daux
    ur(j) = ur(j) + alph1*dur
    theta(j) = theta(j) + alph1*dteta

! Corrections pour le point suivant
    j1=j-1
    if (j1 > mtu) then
      dur1 = dur
      dteta1 = dteta
      dur = uu(j1)*dur1 + vu(j1)*dteta1 + wu(j1)
      dteta = ut(j1)*dur1 + vt(j1)*dteta1 + wt(j1)
      dom = uo(j1)*dur + vo(j1)*dteta + wo(j1)
      daux = ua(j1)*dur + va(j1)*dteta + wa(j1)
    else
      if (j1 == mtu) then
        dur1 = dur
        dteta1 = dteta
        dur = uu(j1)*dur1 + vu(j1)*dteta1 + wu(j1)
        dteta = 0.d0
        dom = uo(j1)*dur1 + vo(j1)*dteta1 + wo(j1)
        daux = ua(j1)*dur1 + va(j1)*dteta1 + wa(j1)
      endif
    endif
   enddo

   write(io_logs,'(a,i3/,11x,a,e8.2)') 'MOMENT CINETIQUE : iteration ',iterad,'alpha : ',alph1
   write(io_logs,'(a,12x,4(i5,1x,e10.4,1x))') ' Plus grand Gi',jgg1,gg1,jgg2,gg2,jgg3,gg3,jgg4,gg4
   write(io_logs,'(a,3x,4(i5,1x,e8.2,1x)/)') ' Plus grande correction',jgdu,gdu,jgdt,gdt,jgda,gda,jgdo,gdo
   write(io_logs,'(1x,a,e9.2,1x,a,e9.2/)') 'Valeur de Z1 :',az1,'Valeur de B1 :',b1

! Si une des corrections est superieure a la correction maximale
!   toleree, il faut effectuer une nouvelle iteration.
   if (idebug > 0) then
     write(*,*) 'call ur1s'
   endif
   call ur1s
   if (idebug > 0) then
     write(*,*) 'call gtilgm'
   endif
   call gtilgm
   if (abs(gdu)>=agdu .or. abs(gda)>=agda .or. abs(gdo)>=agdo .or. abs(gdt)>=agdth .or. abs(gg1)>=agg .or. abs(gg2)>=agg &
       .or. abs(gg3)>=agg .or. abs(gg4)>=agg) then
! On n'a pas obtenue la convergence
     if (iterad < itmax) then
       if (abs(gdu)>=agmax .or. abs(gda)>=agmax .or. abs(gdo)>=agmax .or. abs(gdt)>=agmax .or. abs(gg1)>=agmax &
           .or. abs(gg2)>=agmax .or. abs(gg3)>=agmax .or. abs(gg4)>=agmax) then
         write(io_logs,*)' attention correction trop grande'
         write(io_logs,'(a,i3/,11x,a,e8.2)')'MOMENT CINETIQUE : iteration ',iterad,'alpha : ',alph1
         write(io_logs,'(a,12x,4(i5,1x,e10.4,1x))') ' Plus grand Gi',jgg1,gg1,jgg2,gg2,jgg3,gg3,jgg4,gg4
         write(io_logs,'(a,3x,4(i5,1x,e8.2,1x)/)')' Plus grande correction',jgdu,gdu,jgdt,gdt,jgda,gda,jgdo,gdo
         write(io_logs,'(1x,a,e9.2,1x,a,e9.2/)') 'Valeur de Z1 :',az1,'Valeur de B1 :',b1
       endif
     else
       write(io_logs,*)' itmax atteint sans convergence'
       jterma=1
       endIter = .true.
     endif
   else

! Ici, on a obtenu convergence
     write(io_logs,'(a,i3/,11x,a,e8.2)') 'MOMENT CINETIQUE : iteration ',iterad,'alpha : ',alph1
     write(io_logs,'(a,12x,4(i5,1x,e10.4,1x))') ' Plus grand Gi',jgg1,gg1,jgg2,gg2,jgg3,gg3,jgg4,gg4
     write(io_logs,'(a,3x,4(i5,1x,e8.2,1x)/)') ' Plus grande correction',jgdu,gdu,jgdt,gdt,jgda,gda,jgdo,gdo
     write(io_logs,'(1x,a,e9.2,1x,a,e9.2/)') 'Valeur de Z1 :',az1,'Valeur de B1 :',b1
     endIter = .true.
   endif

  enddo

  return

end subroutine henadv
!======================================================================
subroutine initia
!-----------------------------------------------------------------------
! Sous-routine d'initialisation des valeurs inconnues
!-----------------------------------------------------------------------
use const,only: pi,cst_G,Msol
use caramodele,only: glm
use strucmod,only: rb,H_P,Nabla_mu,rho,q,delt,opac,opact,xb,epsit
use rotmod,only: xlo,xmeg,theta,dlodlr,ht,xmu,aux,gtilde,ebem,ur
use diffadvmod,only: mtu,npasr,D_shear
use abundmod,only: epst

implicit none

integer:: n
real(kindreal):: cb,cq,tqrmr,yx1,yx2,yx3,yy1,yy2,yy3,wpena,wpenb,wfa,wfb,xdtedq,aux1,aux2,aux3,aux4,xdladq,chimu,chit,xktj, &
  rhom,bur1,bur2,bur3,bur4,xura,xurb,xurc,xurd,cur1,cur2,cur3,cur4,xepsmu,xojo,vm
real(kindreal), dimension(ldi):: xlam,xdaux
!-----------------------------------------------------------------------
  if (ncdiff < nminim) then
    if (verbose) then
      write(*,'(a,i2)') ' nombre de couches radiatives inf. a ',nminim
    endif
    return
  endif

  do n=mtu+1,npasr-1
   xlo(n) = log(omegi(n))
   if (xmeg(n) /= 0.d0) then
! Maeder 2009: Eq. (11.41)
     theta(n) = 2.d0/3.d0 * exp(rb(n)) * omegi(n)*omegi(n) / xmeg(n) * dlodlr(n)
   else
     theta(n) = 0.d0
   endif
  enddo
  theta(mtu) = 0.d0
  theta(npasr) = 0.d0

! initialisation de aux : cf. Maeder 2009, Eq. (11.33)
  xlam(mtu:npasr) = exp(rb(mtu:npasr)) / (6.d0*H_P(mtu:npasr)) * Nabla_mu(mtu:npasr)

  do n=mtu+1, npasr-1
   cb = exp(glm-log(4.d0)-log(pi)-2.d0*rb(n)-rho(n))
   cq = exp(q(n))
! tqrmr=-e^Q * M/(4pi r2 rho) : passage de dq a dr
   tqrmr = -cb*cq

   yx1 = q(n)-q(n-1)
   yx3 = q(n)-q(n+1)
   yx2 = q(n+1)-q(n-1)
   yy1 = theta(n-1) / delt(n-1)
   yy2 = theta(n) / delt(n)
   yy3 = theta(n+1) / delt(n+1)
   wpena = (yy2-yy3) / yx3
   wpenb = (yy2-yy1) / yx1
   wfa = yx1 / yx2
   wfb = -yx3 / yx2
   xdtedq = wpena*wfa + wpenb*wfb
! aux1: H_T * d(Theta/delta)/dr   entre dans E_Omega
   aux1 = -ht(n) * xdtedq / tqrmr

   yy1 = xlam(n-1) / delt(n-1)
   yy2 = xlam(n) / delt(n)
   yy3 = xlam(n+1) / delt(n+1)
   wpena = (yy2-yy3) / yx3
   wpenb = (yy2-yy1) / yx1
   wfa = yx1 / yx2
   wfb = -yx3 / yx2
   xdladq = wpena*wfa + wpenb*wfb
! aux2: H_T * d(phi Lambda / delta)/dr   entre dans E_mu
   aux2 = ht(n) * xdladq / tqrmr

   aux3 = theta(n) / delt(n)


   chimu = 0.16d0 / (exp(opac(n)) * xmu(n))-1.d0
   chit = 3.d0 - opact(n) + delt(n)
! aux4:    entre dans E_mu
   aux4 = -xlam(n) * (chimu + 1.d0 / delt(n) * (chit + 1.d0))

   aux(n) = aux1 + aux2 + aux3 + aux4
  enddo
  aux(mtu) = aux(mtu+1)
  aux(npasr) = aux(npasr-1)

! calcul de dAux/dr
  do n=mtu+1,npasr-1
   yx1 = exp(rb(n))-exp(rb(n-1))
   yx3 = exp(rb(n))-exp(rb(n+1))
   yx2 = exp(rb(n+1))-exp(rb(n-1))
   yy1 = aux(n-1)
   yy2 = aux(n)
   yy3 = aux(n+1)
   wpena = (yy2-yy3) / yx3
   wpenb = (yy2-yy1) / yx1
   wfa = yx1 / yx2
   wfb = -yx3 / yx2
   xdaux(n) = wpena*wfa + wpenb*wfb
  enddo
  xdaux(mtu) = xdaux(mtu+1)
  xdaux(npasr) = xdaux(npasr-1)

  do n=mtu,npasr
! initialisation de Ur
   xura = 2.d0*gtilde(n) * (1.d0-omegi(n)*omegi(n) / (2.d0*pi*cst_G*exp(rho(n)))-ebem(n))

   rhom = exp(glm)*(1.d0-exp(q(n))) * 3.d0/(4.d0*pi * exp(3.d0*rb(n)))
! bur1: rho_m/\bar{rho}
   bur1 = rhom / exp(rho(n))
   bur2 = exp(rb(n)) / 3.d0*xdaux(n)
! bur3: partiellement dans 11.67 et 11.73
   bur3 = 2.d0*ht(n) / exp(rb(n)) * (theta(n)/delt(n) - xlam(n)/delt(n))
   bur4 = -2.d0/3.d0 * theta(n)
   xurb = bur1 * (bur2+bur3+bur4)

! cur1: \bar{eps}/eps_m
   cur1 = ebem(n)
   cur2 = aux(n)
! xktj = chi_T plus haut
   xktj = 3.d0-opact(n) + delt(n)
   cur3 = -theta(n) * (epst/delt(n) - xktj/delt(n))
   xepsmu = -4.d0/5.d0 * 1.d0/(xmu(n)*xb(n))
! cur4:
   cur4 = -theta(n)+xlam(n)*(xepsmu+epsit(n))
   xurc = cur1*(cur2+cur3+cur4)

   xurd = -theta(n) + theta(n)*(1.d0-omegi(n)*omegi(n) / (2.d0*pi*cst_G*exp(rho(n))))

   xojo = -5.d0 * D_shear(n) * dlodlr(n) / exp(rb(n))
   ur(n) = xojo
   vm = exp(glm)*(1.d0-exp(q(n)))/Msol
  enddo

  return

end subroutine initia
!======================================================================
subroutine momcon
!-----------------------------------------------------------------------
! repere le type de configuration
!    iconra=0 * entierement convective
!    iconra=1 * entierement radiative
!    iconra=2 rad----rad
!    iconra=3 conv---conv
!    iconra=4 rad----conv
!    iconra=5 conv---rad
! repere les limites des zones convectives
! calcule le moment d inertie des zones convectives
!-----------------------------------------------------------------------
use const,only: pi
use caramodele ,only: glm
use strucmod,only: m,q,vr,rho
use rotmod,only: vsuminenv,vvomeg
use convection,only: nzcon,nxzcon,xmomin,xmocin

implicit none

integer:: n,nc,izin,izex,im
real(kindreal):: dxmomi,dmo,xmi,xme
!-----------------------------------------------------------------------
  if (nzcon == 0 .or. (nxzcon(1) == m .and. nxzcon(2) == 1)) then
    return
  endif

  do n=1,nzcon
   xmomin(n) = 0.d0
   xmocin(n) = 0.d0
   nc = 2*n-1
   izin = nxzcon(nc)
   izex = nxzcon(nc+1)

   im = izex
! Au moment d'inertie de la couche 1 il
! il faudrait ajouter le moment d'inertie de l'enveloppe
   do while (im <= izin)
    if (im /= 1 .and. im /= m) then
      xme = (1.d0-exp(q(im-1)))*exp(glm)
      xmi = (1.d0-exp(q(im+1)))*exp(glm)
      dmo = (xme-xmi)/2.d0
      xmomin(n) = xmomin(n) + exp(2.d0*vr(im)) * dmo * 2.d0/3.d0
!      xmocin(n) = xmocin(n) + vomegi(im) * exp(2.d0*vr(im)) * dmo * 2.d0/3.d0
      xmocin(n) = xmocin(n) + vvomeg(im) * exp(2.d0*vr(im)) * dmo * 2.d0/3.d0
    endif
    if (im == 1) then
      xme = (1.d0-exp(q(im)))*exp(glm)
      xmi = (1.d0-exp(q(im+1)))*exp(glm)
      dmo = (xme-xmi)/2.d0
! On ajoute la contribution de l'enveloppe dans la couche 1.
      xmomin(n) = xmomin(n) + exp(2.d0*vr(im)) * dmo * 2.d0/3.d0 + vsuminenv
      xmocin(n) = xmocin(n) + vvomeg(im) * exp(2.d0*vr(im)) * dmo * 2.d0/3.d0 + vsuminenv*vvomeg(im)
    endif
    if (im == m) then
      dmo = (1.d0-exp(q(im-1)))*exp(glm)/2.d0
      dxmomi = dmo * (3.d0*dmo/(4.d0*pi*exp(rho(im))))**(2.d0/3.d0)*2.d0/5.d0
      xmomin(n) = xmomin(n) + dxmomi
      xmocin(n) = xmocin(n) + vvomeg(im)*dxmomi
      exit
    endif
    im = im+1
   enddo
  enddo

  return

end subroutine momcon
!======================================================================
subroutine ur1s
!-----------------------------------------------------------------------
! calcul de la premiere partie de U2
! xnabyy est le xnabj calcule dans nabgam:
! nabla de l'eq. (6.39) de Maeder 1997, A&A 321, 134
!-----------------------------------------------------------------------
use const,only: pi,cst_G,Lsol
use caramodele,only: glm,gls
use strucmod,only: sb,m,q,rb,Nabla_ad,delt,Nabla_mu
use rotmod,only: xmeg,ur1
use diffadvmod,only: xnabyy

implicit none

integer:: n
real(kindreal):: ura,urb,urc,adramu,xlumi,xmst,rhom,zwi1
!-----------------------------------------------------------------------
  zwi1 = 1.d0 / (exp(sb(1))-1.d0)

  do n=1,m-1
! ATTENTION PROVISOIREMENT NABLA RAD EST UTILISE
! calcul de rhom
   rhom = exp(glm)*(1.d0-exp(q(n)))*3.d0 /(4.d0*pi*exp(3.d0*rb(n)))
   xmst = exp(glm)*(1.d0-exp(q(n)))*(1.d0-omegi(n)*omegi(n)/(2.d0*pi*cst_G*rhom))
! L/(M* X g)
   xlumi = (exp(sb(n))-1.d0) * zwi1*gls*Lsol
   ura = xlumi / (xmst*xmeg(n))
! P/(rho Cp T)
   urb = Nabla_ad(n) / delt(n)
! 1./(Ad-nab)
   adramu = Nabla_ad(n) - xnabyy(n) + Nabla_mu(n)/delt(n)
   if (adramu > 0.0d0) then
     urc = 1.d0 / (adramu)
   else
     urc = 100.d0
   endif
   ur1(n) = ura*urb*urc
  enddo

! MISE A ZERO AU CENTRE PROVISOIREMENT !
  ur1(m) = 0.d0

  return

end subroutine ur1s
!======================================================================
subroutine vrcirc
!-----------------------------------------------------------------------
! Sous-routine de calcul de la composante tangentielle
! de la vitesse de circulation meridienne
!-----------------------------------------------------------------------
use strucmod,only: m,q,rb,rho
use rotmod,only: vcirc,ur
use diffadvmod,only: npasr,mtu,xueff,ursmooth
use SmallFunc,only: SmoothProfile

implicit none

integer, parameter:: WS = 5
integer:: n,ilim1,ilim2
real(kindreal),parameter:: Mrlim1 = 0.65d0,Mrlim2 =0.80d0
real(kindreal):: yx1,yx2,yx3,yy1,yy2,yy3,wpena,wpenb,wfa,wfb,qlim1,qlim2,pond
real(kindreal),dimension(ldi):: xdulnr,xdulnrsmooth,xxdulnr,uursmooth,xxdulnrsmooth,rrb,uur,xldlnr,xddulnr,ursmooth1,&
                                xdulnrsmooth1,ursmooth2,xdulnrsmooth2,ursmooth3,xdulnrsmooth3
!-----------------------------------------------------------------------
! initialisation
  if (idebug > 0) then
    write(*,*) 'VRCIRC, initialisations'
  endif
  vcirc(:)=0.d0

  ursmooth(:) = ur(:)
  ursmooth1(:) = ur(:)
  ursmooth2(:) = ur(:)
  xdulnrsmooth(:) = 0.d0
  xdulnrsmooth1(:) = 0.d0
  xdulnrsmooth2(:) = 0.d0
  xdulnr(:) = 0.d0
  xddulnr(:) = 0.d0
  uursmooth(:) = ur(:)
  xxdulnrsmooth(:) = 0.d0
  xxdulnr(:) = 0.d0
  uur(:) = 0.d0
  rrb(:) = 0.d0

  qlim1 = log(1.d0-Mrlim1)
  qlim2 = log(1.d0-Mrlim2)
  ilim1 = 1
  ilim2 = 1
  if (idebug > 0) then
    write(*,*) 'first do loop'
  endif
  do n=1,m
   if (q(n)  <=  qlim2) then
     ilim2 = n
   elseif (q(n)  <=  qlim1) then
     ilim1 = n
   endif
  enddo

  if (idebug > 0) then
    write(*,*) 'check npasr-mtu > 3'
  endif

  if (ncdiff < nminim .and. verbose) then
    write(*,'(a,i2)') 'return: nombre de couches radiatives inf. a ',nminim
    return
  endif

! cas intermediaires
  if (idebug > 0) then
    write(*,*) 'second do loop'
  endif
  do n=mtu+1, npasr-1
   yx1=rb(n)-rb(n-1)
   yx3=rb(n)-rb(n+1)
   yx2=rb(n+1)-rb(n-1)
   yy1=rho(n-1)
   yy2=rho(n)
   yy3=rho(n+1)
   wpena=(yy2-yy3)/yx3
   wpenb=(yy2-yy1)/yx1
   wfa=yx1/yx2
   wfb=-yx3/yx2
   xldlnr(n)=wpena*wfa+wpenb*wfb
   yy1=ur(n-1)
   yy2=ur(n)
   yy3=ur(n+1)
   wpena=(yy2-yy3)/yx3
   wpenb=(yy2-yy1)/yx1
   wfa=yx1/yx2
   wfb=-yx3/yx2
   xddulnr(n)=wpena*wfa+wpenb*wfb
  enddo

  if (idebug > 0) then
    write(*,*) 'ext border'
  endif

! bords ext
  if (mtu /= 1) then
    yx1=rb(mtu)-rb(mtu-1)
    yx3=rb(mtu)-rb(mtu+1)
    yx2=rb(mtu+1)-rb(mtu-1)
    yy1=rho(mtu-1)
    yy2=rho(mtu)
    yy3=rho(mtu+1)
    wpena=(yy2-yy3)/yx3
    wpenb=(yy2-yy1)/yx1
    wfa=yx1/yx2
    wfb=-yx3/yx2
    xldlnr(mtu)=wpena*wfa+wpenb*wfb
  else
    xldlnr(mtu)=(rho(1)-rho(2))/(rb(1)-rb(2))
  endif
  xddulnr(mtu)=(ur(mtu)-ur(mtu+1))/(rb(mtu)-rb(mtu+1))
  if (idebug > 0) then
    write(*,*) 'int. border'
  endif
! bords int
  if (npasr /= m) then
    yx1=rb(npasr)-rb(npasr-1)
    yx3=rb(npasr)-rb(npasr+1)
    yx2=rb(npasr+1)-rb(npasr-1)
    yy1=rho(npasr-1)
    yy2=rho(npasr)
    yy3=rho(npasr+1)
    wpena=(yy2-yy3)/yx3
    wpenb=(yy2-yy1)/yx1
    wfa=yx1/yx2
    wfb=-yx3/yx2
    xldlnr(npasr)=wpena*wfa+wpenb*wfb
  else
    xldlnr(npasr)=(rho(m-1)-rho(m))/(rb(m-1)-rb(m))
  endif
  xddulnr(npasr)=(ur(npasr-1)-ur(npasr))/(rb(npasr-1)-rb(npasr))

  rrb(1:m) = rb(m:1:-1)
  uur(1:m) = ur(m:1:-1)
  xxdulnr(1:m) = xddulnr(m:1:-1)

  if (idebug > 0) then
    write(*,*) 'call SmoothProfile 1'
  endif
  uursmooth = uur
  xxdulnrsmooth = xxdulnr
  call SmoothProfile(uur,uursmooth,xxdulnrsmooth,rrb,m-npasr+1,m-mtu+1,WS,m)

  ursmooth1(1:m) = uursmooth(m:1:-1)
  xdulnrsmooth1(1:m) = xxdulnrsmooth(m:1:-1)/xueff
  xdulnrsmooth1(mtu) = xddulnr(mtu)/xueff

  if (idebug > 0) then
    write(*,*) 'call SmoothProfile 2'
  endif
  uursmooth = uur
  xxdulnrsmooth = xxdulnr
  call SmoothProfile(uur,uursmooth,xxdulnrsmooth,rrb,m-npasr+1,m-mtu+1,6*WS,m)

  ursmooth2(1:m) = uursmooth(m:1:-1)
  xdulnrsmooth2(1:m) = xxdulnrsmooth(m:1:-1)/xueff
  xdulnrsmooth2(mtu) = xddulnr(mtu)/xueff

  if (idebug > 0) then
    write(*,*) 'call SmoothProfile 3'
  endif
  uursmooth = uur
  xxdulnrsmooth = xxdulnr
  call SmoothProfile(uur,uursmooth,xxdulnrsmooth,rrb,m-npasr+1,m-mtu+1,20*WS,m)

  ursmooth3(1:m) = uursmooth(m:1:-1)
  xdulnrsmooth3(1:m) = xxdulnrsmooth(m:1:-1)/xueff
  xdulnrsmooth3(mtu) = xddulnr(mtu)/xueff

  if (idebug > 0) then
    write(*,*) 'third do loop'
  endif
  do n=1,m
   if (n <= ilim2) then
     if (n >= ilim2 - 25) then
       pond = 0.04*real(n-ilim2+25)
       ursmooth(n) = pond*ursmooth2(n) + (1.d0-pond)*ursmooth3(n)
       xdulnrsmooth(n) = pond*xdulnrsmooth2(n) + (1.d0-pond)*xdulnrsmooth3(n)
     else
       ursmooth(n) = ursmooth3(n)
       xdulnrsmooth(n) = xdulnrsmooth3(n)
     endif
   elseif (n <= ilim1) then
     if (n >= ilim1 - 25) then
       pond = 0.04*real(n-ilim1+25)
       ursmooth(n) = pond*ursmooth1(n) + (1.d0-pond)*ursmooth2(n)
       xdulnrsmooth(n) = pond*xdulnrsmooth1(n) + (1.d0-pond)*xdulnrsmooth2(n)
     else
       ursmooth(n) = ursmooth2(n)
       xdulnrsmooth(n) = xdulnrsmooth2(n)
     endif
   else
     ursmooth(n) = ursmooth1(n)
     xdulnrsmooth(n) = xdulnrsmooth1(n)
   endif
  enddo

  xxdulnr(1:m)=xxdulnr(1:m)/xueff
  xddulnr(1:m)=xddulnr(1:m)/xueff

  if (idebug > 0) then
    write(*,*) 'fourth do loop, vcirc'
  endif
  do n=mtu,npasr
   vcirc(n)=ursmooth(n)/(3.d0*xueff)+ursmooth(n)/(6.d0*xueff)*xldlnr(n)+xdulnrsmooth(n)/6.d0
  enddo

  return

end subroutine vrcirc
!======================================================================
subroutine ziadv
!-----------------------------------------------------------------------
! Calcul des coefficients necessaires dans la methode de henyey
!   lors du calcul de l'evolution du moment cinetique ;
!   condition limite au centre de l'etoile

! Z_OMEGA est appele par HEN_OMEGA

! npasr est le numero de la derniere couche de la zone radiative
! ATTENTION APPORTER CHANGEMENT POUR CONDITIONS AU BORD INTERIEUR
!-----------------------------------------------------------------------
use const, only: pi
use strucmod, only: H_P,rb,rho
use rotmod, only: ur
use diffadvmod,only:npasr,ibint,xbint,xcint
use equadiffmod, only: az1,az1o1,az1u1
use timestep, only: dzeit

implicit none

real(kindreal):: gamml,xeffl,tamp1,vava
!-----------------------------------------------------------------------
  if (ibint == 0) then
    az1 = ur(npasr)
    az1o1 = 0.d0
    az1u1 = 1.d0
  else
! calcul de la vitesse des cellules convectives en npasr+1
    gamml = H_P(npasr) / (6.d0*exp(rb(npasr)))
    xeffl = 1.0d0
    tamp1 = 4.d0*pi*0.2d0*exp(4.d0*rb(npasr)+rho(npasr))*2.d0/3.d0
    vava = dzeit * tamp1 / xbint
    az1 = omegi(npasr) * xeffl * ur(npasr) * vava - (omegi(npasr)-xcint/xbint)
    az1o1 = xeffl*ur(npasr)*vava-1.d0
    az1u1 = omegi(npasr)*vava*xeffl
  endif

  return

end subroutine ziadv
!======================================================================
subroutine advect
!-----------------------------------------------------------------------
use inputparam,only: phase,xcn,itminc,idebug
use caramodele,only: nwmd,glm,firstmods,inum,xLstarbefHen,gms
use abundmod,only: x
use equadiffmod,only: iadnok,jterma
use strucmod,only: m,q,rb,r
use rotmod,only: deladv,vvomeg,omegp,omegd,vsuminenv,xldoex,BTotal_EndAdvect,Flux_remaining,timestep_control
use diffadvmod,only: inoadv,xueff
use timestep,only: dzeit

implicit none

logical:: AdvecTest
real(kindreal), dimension(ldi):: xmr,oommgg
real(kindreal):: alph1,btota,xmocin,dbrmr,dbrmr1,dbrmrm,btoto,xdmax,xdibb,btota1,btota2, &
                 max_tolerance = 1.d-3,Moment_inertie,xMoCinScale
integer:: inzr,npair,n,flag_girl=0
!-----------------------------------------------------------------------
  if (inum == 0) then
    xueff=(1.d0+xcn)
  else if (inum > 0) then
    xueff=2.d0
  endif

! Initialisation de AdvecTest
  AdvecTest = .true.
  max_tolerance = 1.d-3

! initialisation de deladv et de omegi
  deladv(:)=0.d0
  omegi(1:m)=vvomeg(1:m)
  iadnok=0

  if (iadvec == 1) then
! calcul du moment d'inertie des zones convectives
    if (idebug > 0) then
      write(*,*) 'in ADVEC, call momcon'
    endif
    call momcon
! determine si il existe une region ou l'advection
! est appliquee. Si oui, determine le type de condition au bord
    if (idebug > 0) then
      write(*,*) 'call confi'
    endif
    call confi
    if (inoadv /= 1) then
! calcul de l'effet de l'advection entre npasr et mtu
! calcul de certaines quantites, dans toute l'etoile
! gbar appelle geomeang
      if (idebug > 0) then
        write(*,*) 'call gbar'
      endif
      call gbar
      if (idebug > 0) then
        write(*,*) 'call echtem'
      endif
      call echtem
      if (idebug > 0) then
        write(*,*) 'call ur1s'
      endif
      call ur1s
      if (idebug > 0) then
        write(*,*) 'call gtilgm'
      endif
      call gtilgm
      if (idebug > 0) then
        write(*,*) 'call ebems'
      endif
      call ebems
! initialisation des inconnues
      inzr=1
      if (idebug > 0) then
        write(*,*) 'call initia'
      endif
      call initia
      alph1=1.0d0
! resolution du systeme d'equations
! henadv appelle badv, giadv et ziadv
      if (idebug > 0) then
        write(*,*) 'call henadv'
      endif
      call henadv(alph1,flag_girl)
      if (idebug > 0) then
        write(*,*) 'call coradv'
      endif
      call coradv
! calcule de la composante tangentielle de la
! vitesse de circulation meridienne
      if (idebug > 0) then
        write(*,*) 'call vrcirc'
      endif
      call vrcirc
      if (idebug > 0) then
        write(*,*) 'end vrcirc'
      endif
      npair=mod(nwmd,2)
      if (idebug > 0) then
        write(*,*) 'npair=',npair
      endif
! 18 octobre 2005
! Si imagn /= 0 alors on n'applique pas la correction due a l'advection.
! On ne passe que dans la partie diffusive del'equation de transport
!      if (imagn > 0) npair=0
      if (npair == 1) then
!-----------------------------------------------------------------------
! cas ou seulement l'advection est calculee
        write(io_logs,*) 'PASSAGE PAR ADVECTION'
        if (idebug > 0) then
          write(*,*) 'PASSAGE PAR ADVECTION'
        endif
! Ici, on pre-corrige le profile de rotation obtenu avec la conservation locale du
! moment cinetique des effets de la diffusion, approximes comme omegd - vvomeg.
! vvomeg est le profile de omega du pas de temps precedent, et est utilise en entree
! dans diffom. omegd est le profile sorti par diffom.
! ATTENTION: ceci est un peu artificiel, et ne permet EN AUCUN CAS de garantir
! la conservation du moment cinetique. Pour cette raison, on peut avoir besoin de le
! corriger (fait plus loin).

        oommgg(1:m)=vvomeg(1:m)+deladv(1:m)
        xmr(1:m)=exp(glm)*(1.d0-exp(q(1:m)))

! [Modif CG]
! Ici, on s'assure que la diffusion conserve bien le moment cinetique. C'est-a-dire
! que l'on verifie que le moment cinetique apres diffusion btoto est egal au moment
! cinetique avant diffusion btota +/- le flux entrant/sortant de moment cinetique.
! Il s'agit ici d'un test purement numerique de la qualite de la diffusion.

        btota=0.d0
        do n=2,m-1
         xmocin=exp(2.d0*rb(n))*vvomeg(n)*2.d0/3.d0
         dbrmr=xmocin*(xmr(n-1)-xmr(n+1))/2.d0
         btota=btota+dbrmr
        enddo
        xmocin=2.d0/3.d0*exp(2.d0*rb(1))*vvomeg(1)
        dbrmr1=xmocin*(xmr(1)-xmr(2))/2.d0
        dbrmrm=2.d0/5.d0*xmr(m-1)/2.d0*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*vvomeg(m)
        btota=btota+dbrmr1+dbrmrm+vsuminenv*vvomeg(1) + (xldoex+Flux_remaining)*dzeit

        btoto=0.d0
        do n=2,m-1
         xmocin=exp(2.d0*rb(n))*oommgg(n)*2.d0/3.d0
         dbrmr=xmocin*(xmr(n-1)-xmr(n+1))/2.d0
         btoto=btoto+dbrmr
        enddo
        xmocin=2.d0/3.d0*exp(2.d0*rb(1))*oommgg(1)
        dbrmr1=xmocin*(xmr(1)-xmr(2))/2.d0
        dbrmrm=2.d0/5.d0*xmr(m-1)/2.d0*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*oommgg(m)
        btoto=btoto+dbrmr1+dbrmrm+vsuminenv*oommgg(1)

        timestep_control = abs((xldoex+Flux_remaining)*dzeit/(vsuminenv*oommgg(1)))

        if (idebug > 0 .or. itminc == 1) then
          write(*,*) 'CONSERVATION ADVECTION: ', abs(btoto/btota),timestep_control,Flux_remaining*dzeit/btota
        endif

! [Modif CG]
! On ajout ici un test qui verifie la conservation du moment cinetique
! total de l'etoile durant le processus d'advection.
! Interruption si la variation est superieure a 1d-5.
        if (abs(btoto/btota -1.d0) > 1.d-2) then
          if (itminc == 1) then
            write(io_logs,'(a,a)') 'Total angular momentum variation during advection greater than 10^-2.'
            write(io_logs,'(2(a,d14.8))') 'Linitial = ',btota,'      Lfinal = ',btoto
            write(*,*) 'Advection not applied in this model.'
            write(io_logs,*) 'ADVECTION : not applied in this model.'
            rewind(io_runfile)
            write(io_runfile,*) nwmd,': Problem during advection ==> STOP'
            call safe_stop('problem during advection')
          endif
        else if (abs(btoto/btota -1.d0) > max_tolerance) then
          if (verbose) then
            write(*,*) 'Total angular momentum variation during advection greater than',max_tolerance
            write(*,*) 'L avant : ', btota
            write(*,*) 'L apres : ', btoto
          endif
          if (.not. firstmods) then
            AdvecTest = .false.
            write(io_logs,'(a,1pe7.1)') 'Total angular momentum variation during advection greater than',max_tolerance
            write(io_logs,'(2(a,d14.8))') 'Linitial = ',btota,'      Lfinal = ',btoto
            write(*,*) 'Advection not applied in this model.'
            write(io_logs,*) 'ADVECTION : not applied in this model.'
            if (itminc == 1) then
              write(*,*) "Problem with conservation of angular momentum during advection."
              rewind(io_runfile)
              write(io_runfile,*) nwmd,': Problem during advection ==> STOP'
              call safe_stop('problem during advection')
            endif
          else
            write(*,*) 'Advection applied nevertheless.'
            write(io_logs,*) 'ADVECTION : applied in this model.'
            write(io_logs,'(2(a,d14.8))') 'Linitial = ',btota,'      Lfinal = ',btoto
          endif
        endif
! [/Modif]

        xdmax=max_tolerance
        xdibb=abs(btoto-btota)/btota
        if (xdibb > xdmax .or. jterma == 1 .or. (.not.AdvecTest) .or. inoadv == 2) then
          iadnok=1
          if (itminc == 1) then
            write(io_sfile,'(1x,i5,a)') nwmd,'!!!!! NOT OK ADVECTION !!!!!!'
          endif
          if (verbose) then
            write(*,*) ' PAS D ADVECTION'
          endif
! Si l'advection est mal calculee, alors on utilise la conservation locale du moment
! cinetique seulement.
          omegi(1:m) = omegp(1:m)
        else
! Si l'advection est prise en compte malgre une trop grande variation,
! on stoppe l'execution.
          if (itminc == 1) then
            write(io_sfile,'(1x,i5,a)') nwmd,' OK ADVECTION'
          endif
! Si on applique l'advection, on applique la meme procedure que lors de la
! pre-correction ci-dessus, mais cette fois-ci sur omegi.
          omegi(1:m) = omegp(1:m) + deladv(1:m)
        endif

! [Modif CG]
! Ici, on s'assure que le moment cinetique final soit bien celui qu'on attend
! (voir remarque precedente). On compare donc le nouveau moment cinetique calcule
! avec le profile de omegi et la nouvelle structure rb, que l'on compare avec le
! moment cinetique initial du debut du calcul +/- le flux gagne/perdu.
        btota=0.d0
        do n=2,m-1
         xmocin=exp(2.d0*rb(n))*omegi(n)*2.d0/3.d0
         dbrmr=xmocin*(xmr(n-1)-xmr(n+1))/2.d0
         btota=btota+dbrmr
        enddo
        xmocin=2.d0/3.d0*exp(2.d0*rb(1))*omegi(1)
        dbrmr1=xmocin*(xmr(1)-xmr(2))/2.d0
        dbrmrm=2.d0/5.d0*xmr(m-1)/2.d0*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*omegi(m)
        btota=btota+dbrmr1+dbrmrm+vsuminenv*omegi(1)!+(xldoex+Flux_remaining)*dzeit

! Verifie que l'on n'est pas trop eloigne du moment initial.
        if (abs((xLstarbefHen+(xldoex+Flux_remaining)*dzeit)/btota - 1.d0) > max_tolerance) then
          if (verbose) then
            write(*,*) 'Problem during advection.'
            write(*,*) 'Old angular momentum: ', xLstarbefHen,' New angular momentum: ', btota
          endif
          write(io_logs,*) 'Problem during advection.'
          write(io_logs,'(2(a,d14.8))') 'Old angular momentum: ',xLstarbefHen,' New angular momentum: ', btota
          if (x(m) < 7.d-1) then
            rewind(io_runfile)
            write(io_runfile,*) nwmd,': Problem during advection ==> STOP'
            call safe_stop('problem during advection')
          endif
        endif
! Corrige le profil de rotation afin de garantir la conservation du moment angulaire.
        omegi(1:m) = omegi(1:m)*(xLstarbefHen+(xldoex+Flux_remaining)*dzeit)/btota
        btota=0.d0
        do n=2,m-1
         xmocin=exp(2.d0*rb(n))*omegi(n)*2.d0/3.d0
         dbrmr=xmocin*(xmr(n-1)-xmr(n+1))/2.d0
         btota=btota+dbrmr
        enddo
        xmocin=2.d0/3.d0*exp(2.d0*rb(1))*omegi(1)
        dbrmr1=xmocin*(xmr(1)-xmr(2))/2.d0
        dbrmrm=2.d0/5.d0*xmr(m-1)/2.d0*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*omegi(m)
        btota=btota+dbrmr1+dbrmrm+vsuminenv*omegi(1)!+(xldoex+Flux_remaining)*dzeit

! Sauvegarde des moments cinetiques de chaque coquille, afin de corriger le
! moment cinetique calcule avec r et omegi (actuellement vr et omegi).
        do n=2,m-1
         xmocin=exp(2.d0*rb(n))*omegi(n)*2.d0/3.d0
         xMoCinScale=xmocin*(xmr(n-1)-xmr(n+1))/2.d0
         Moment_inertie = exp(2.d0*r(n))*2.d0/3.d0*(xmr(n-1)-xmr(n+1))/2.d0
         omegi(n) = xMoCinScale/Moment_inertie
        enddo
        xmocin=2.d0/3.d0*exp(2.d0*rb(1))*omegi(1)
        xMoCinScale=xmocin*(xmr(1)-xmr(2))/2.d0 + vsuminenv*omegi(1)
        Moment_inertie = 2.d0/3.d0*exp(2.d0*r(1))*(xmr(1)-xmr(2))/2.d0 + vsuminenv
        omegi(1) = xMoCinScale/Moment_inertie
        xMoCinScale=2.d0/5.d0*xmr(m-1)/2.d0*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*omegi(m)
        Moment_inertie = 2.d0/5.d0*xmr(m-1)/2.d0*exp(r(m-1))/(2.d0)**(1.d0/3.d0)*exp(r(m-1))/(2.d0)**(1.d0/3.d0)
        omegi(m) = xMoCinScale/Moment_inertie
! [/Modif]

        BTotal_EndAdvect = btota - Flux_remaining*dzeit
        return
      endif
    endif
  endif
!-----------------------------------------------------------------------
! cas ou seulement la diffusion est calculee
  write(io_logs,*) 'PASSAGE PAR DIFFUSION'
  if (idebug > 0) then
    write(*,*) 'PASSAGE PAR DIFFUSION'
  endif
  xmr(1:m)=exp(glm)*(1.d0-exp(q(1:m)))

  if (phase == 1 .and. .not.firstmods .and. gms > 1.4d0) then
    max_tolerance = 1.d-5
  else
    max_tolerance = 1.d-3
  endif

! [Modif CG]
! Ici, on s'assure que la diffusion conserve bien le moment cinetique. C'est-a-dire
! que l'on verifie que le moment cinetique apres diffusion btota2 est egal au moment
! cinetique avant diffusion btota1 +/- le flux entrant/sortant de moment cinetique.
! Il s'agit ici d'un test purement numerique de la qualite de la diffusion.

  btota=0.d0
  do n=2,m-1
   xmocin=exp(2.d0*rb(n))*vvomeg(n)*2.d0/3.d0
   dbrmr=xmocin*(xmr(n-1)-xmr(n+1))/2.d0
   btota=btota+dbrmr
  enddo
  xmocin=2.d0/3.d0*exp(2.d0*rb(1))*vvomeg(1)
  dbrmr1=xmocin*(xmr(1)-xmr(2))/2.d0
  dbrmrm=2.d0/5.d0*xmr(m-1)/2.d0*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*vvomeg(m)
  btota1=btota+dbrmr1+dbrmrm +vsuminenv*vvomeg(1)+(xldoex+Flux_remaining)*dzeit

  btota=0.d0
  do n=2,m-1
   xmocin=exp(2.d0*rb(n))*omegd(n)*2.d0/3.d0
   dbrmr=xmocin*(xmr(n-1)-xmr(n+1))/2.d0
   btota=btota+dbrmr
  enddo
  xmocin=2.d0/3.d0*exp(2.d0*rb(1))*omegd(1)
  dbrmr1=xmocin*(xmr(1)-xmr(2))/2.d0
  dbrmrm=2.d0/5.d0*xmr(m-1)/2.d0*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*omegd(m)
  btota2=btota+dbrmr1+dbrmrm+vsuminenv*omegd(1)

  timestep_control = abs((xldoex+Flux_remaining)*dzeit/(vsuminenv*omegd(1)))

! [Modif CG]
! On ajoute ici un test qui verifie la conservation du moment cinetique
! total de l'etoile durant le processus de diffusion.
! Interruption si la variation est superieure a 1d-5.
  if (idebug > 0 .or. itminc == 1) then
    write(*,*) 'CONSERVATION DIFFUSION: ',abs(btota2/btota1),(xldoex+Flux_remaining)*dzeit
  endif

  if (abs(btota2/btota1 -1.d0) > max_tolerance) then
    if (verbose .or. itminc == 1) then
      write(io_logs,*) 'Angular momentum variation during diffusion: ', abs(btota2/btota1 -1.d0)
      write(*,*) 'Angular momentum variation during diffusion: ', abs(btota2/btota1 -1.d0)
    endif
    if (itminc == 1) then
      rewind(io_runfile)
      write(io_runfile,*) nwmd,': Ang. mom. variation too large during diffusion ==> STOP'
      write(*,'(a,es7.1,a)') 'Total angular momentum variation during diffusion greater than ',max_tolerance,'. Aborting...'
      call safe_stop('Total angular momentum variation during diffusion too large')
    endif
  endif
! [/Modif]

! Ici, on corrige le profile de rotation obtenu avec la conservation locale du
! moment cinetique des effets de la diffusion, approximes comme omegd - vvomeg.
! vvomeg est le profile de omega du pas de temps precedent, et est utilise en entree
! dans diffom. omegd est le profile sorti par diffom.
! ATTENTION: ceci est un peu artificiel, et ne permet EN AUCUN CAS de garantir
! la conservation du moment cinetique. Pour cette raison, on peut avoir besoin de le
! corriger (fait plus loin).

  omegi(1:m)= omegp(1:m) + (omegd(1:m)-vvomeg(1:m))

! [Modif CG]
! Ici, on s'assure que le moment cinetique final soit bien celui qu'on attend
! (voir remarque precedente). On compare donc le nouveau moment cinetique calcule
! avec le profile de omegi et la nouvelle structure rb, que l'on compare avec le
! moment cinetique initial du debut du calcul +/- le flux gagne/perdu.
  btota=0.d0
  do n=2,m-1
   xmocin=exp(2.d0*rb(n))*omegi(n)*2.d0/3.d0
   dbrmr=xmocin*(xmr(n-1)-xmr(n+1))/2.d0
   btota=btota+dbrmr
  enddo
  xmocin=2.d0/3.d0*exp(2.d0*rb(1))*omegi(1)
  dbrmr1=xmocin*(xmr(1)-xmr(2))/2.d0
  dbrmrm=2.d0/5.d0*xmr(m-1)/2.d0*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*omegi(m)
  btota=btota+dbrmr1+dbrmrm+vsuminenv*omegi(1)

! Verifie que l'on n'est pas trop eloigne du moment initial.
  if (abs((xLstarbefHen+(xldoex+Flux_remaining)*dzeit)/btota - 1.d0) > max_tolerance &
     .and. itminc == 1) then
    write(*,*) 'Old angular momentum: ', xLstarbefHen,' New angular momentum: ', btota
    rewind(io_runfile)
    write(io_runfile,*) nwmd,': Problem during diffusion ==> STOP'
    call safe_stop('Problem during diffusion.')
  endif

! Corrige le profil de rotation afin de garantir la conservation du moment angulaire.
  omegi(1:m) = omegi(1:m)*(xLstarbefHen+(xldoex+Flux_remaining)*dzeit)/btota
! Sauvegarde des moments cinetiques de chaque coquille, afin de corriger le
! moment cinetique calcule avec r et omegi (actuellement vr et omegi).
  btota=0.d0
  do n=2,m-1
   xmocin=exp(2.d0*rb(n))*omegi(n)*2.d0/3.d0
   dbrmr=xmocin*(xmr(n-1)-xmr(n+1))/2.d0
   btota=btota+dbrmr
  enddo
  xmocin=2.d0/3.d0*exp(2.d0*rb(1))*omegi(1)
  dbrmr1=xmocin*(xmr(1)-xmr(2))/2.d0
  dbrmrm=2.d0/5.d0*xmr(m-1)/2.d0*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*omegi(m)
  btota=btota+dbrmr1+dbrmrm+vsuminenv*omegi(1)!+(xldoex+Flux_remaining)*dzeit

  do n=2,m-1
   xmocin=exp(2.d0*rb(n))*omegi(n)*2.d0/3.d0
   xMoCinScale=xmocin*(xmr(n-1)-xmr(n+1))/2.d0
   Moment_inertie = exp(2.d0*r(n))*2.d0/3.d0*(xmr(n-1)-xmr(n+1))/2.d0
   omegi(n) = xMoCinScale/Moment_inertie
  enddo
  xmocin=2.d0/3.d0*exp(2.d0*rb(1))*omegi(1)
  xMoCinScale=xmocin*(xmr(1)-xmr(2))/2.d0 + vsuminenv*omegi(1)
  Moment_inertie = 2.d0/3.d0*exp(2.d0*r(1))*(xmr(1)-xmr(2))/2.d0 + vsuminenv
  omegi(1) = xMoCinScale/Moment_inertie
  xMoCinScale=2.d0/5.d0*xmr(m-1)/2.d0*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*omegi(m)
  Moment_inertie = 2.d0/5.d0*xmr(m-1)/2.d0*exp(r(m-1))/(2.d0)**(1.d0/3.d0)*exp(r(m-1))/(2.d0)**(1.d0/3.d0)
  omegi(m) = xMoCinScale/Moment_inertie

  BTotal_EndAdvect = btota - Flux_remaining*dzeit
! [/Modif]

  return

end subroutine advect
!======================================================================
end module advection
