module nablas

use evol, only: kindreal,ldi

implicit none

private
public:: grapmui
public:: nabla,nabgam

contains
!***********************************************************************
!=======================================================================
subroutine nabla
!-----------------------------------------------------------------------
! Derniere version : 28 septembre 1992
!-----------------------------------------------------------------------
  use const,only: cst_G,um,cst_k,cst_u
  use inputparam,only: irot
  use caramodele ,only: hh6
  use equadiffmod,only: ccg1
  use EOS,only: rh1,toni,rhe,pl,rht1,uta,num
  use strucmod,only: r,j1,p,t,m,j,q,s,zrad1,ccrad1,cap1,rad1,zradm,zrad,radm,cap,vmye,beta1,vmyo,adi1,adim,adi,adip1,adit1
  use rotmod,only: omegi
  use geomod, only: rpsi_max,geom
  use SmallFunc,only: exphi

  implicit none

  real(kindreal):: pt368,xpsij,xpsij1,xfpj,xraj,dxfpj,dxraj,xfpj1,xraj1,dxfpj1,dxraj1,xra,pfak,urt,cp_nablar,hfak
!-----------------------------------------------------------------------
! cf Patenaude 74 eq. B.34: critere de degenerescence
  pt368=p(j1)-2.5d0*t(j1)+3.68d0

  if (irot == 1) then
    if (j1 /= m) then
      xpsij=exp(r(j))/((cst_G*exp(ccg1)*(1.d0-exp(q(j))))/(omegi(j)*omegi(j)))**(1.d0/3.d0)
      xpsij1=exp(r(j1))/((cst_G*exp(ccg1)*(1.d0-exp(q(j1))))/(omegi(j1)*omegi(j1)))**(1.d0/3.d0)
      if (xpsij > rpsi_max) xpsij=0.9999999999d0*rpsi_max
      if (xpsij1 > rpsi_max) xpsij1=0.9999999999d0*rpsi_max
      call geom(xpsij,xfpj,xraj,dxfpj,dxraj)
      call geom(xpsij1,xfpj1,xraj1,dxfpj1,dxraj1)
      xra=sqrt(xraj*xraj1)
    else
      xpsij=exp(r(j))/((cst_G*exp(ccg1)*(1.d0-exp(q(j))))/(omegi(j)*omegi(j)))**(1.d0/3.d0)
      if (xpsij > rpsi_max) xpsij=0.9999999999d0*rpsi_max
      call geom(xpsij,xfpj,xraj,dxfpj,dxraj)
      xraj1=1.d0
      xra=sqrt(xraj*xraj1)
    endif
  else   !   irot=0
    xraj1=1.d0
    xra=1.d0
  endif   !   irot

  if ((j1 - m) < 0) then
! gradient radiatif dans la premiere coquille, sans correction rotation.
    zrad1=-exphi(s(j1))*exp(ccrad1+hh6+cap1+p(j1)-4.d0*t(j1))/exphi(q(j1))
! correction rotation
    rad1=xraj1*zrad1
! moyenne des termes.
    zradm = 0.5d0*(zrad1+zrad)
    radm=xra*zradm
  else
! gradient radiatif de l'interieur.
    zradm =-exphi(s(m-1))*exp(ccrad1+hh6+0.5d0*(cap1+cap+p(m)+p(m-1))-2.d0*(t(m)+t(m-1)))/exphi(q(m-1))
    radm=xra*zradm
  endif

  if (rh1 > 4.d0*um .and. t(j1) > 7.d0*um .and. pt368 > 0.d0) then
! pfak=rho T / P
    pfak=toni*rhe*vmye/pl
    urt=((3.d0*beta1-4.d0)*rht1+12.d0*(1.d0-beta1))/pfak
! cf Patenaude 74 eq. B.60: cp dans le cas degenere
    cp_nablar=3.d0*cst_k/(vmyo*2.d0*cst_u)+uta+urt
    adi1=(-rht1/pfak)/cp_nablar
    hfak = cp_nablar*pfak
  else
! Denominateur de B.63, Patenaude 74
    hfak=-(4.d0-1.5d0*beta1)*rht1+6.d0*(1.d0-beta1)
    if (num /= 0) then
      adi1=-rht1/hfak
    endif
    if (num == -1000) then
      num = 0
    endif
  endif

  adim=0.5d0*(adi1+adi)
  adip1=-adi1*4.d0*(1.d0-beta1)*(1.d0/adi1-4.d0+1.5d0*beta1-(1.5d0-0.375d0*rht1)*beta1**2.d0)/(hfak*beta1**2.d0)
  adit1=-4.d0*adip1

  return

end subroutine nabla

!=======================================================================
subroutine nabgam
!-----------------------------------------------------------------------
! Cette sous-routine calcule le nabla dans toute l'etoile
! Eq. (6.39)  du Papier II (Maeder 1997, A&A 321, 134)

! Sous-routine appelee par NABLA, pour le calcul du nabla
! a utiliser dans GI et ZI pour le transfert du rayonnement
!-----------------------------------------------------------------------
  use const,only: Msol,cst_G,cst_a,cst_c,pi
  use inputparam,only: iledou
  use caramodele ,only: gms
  use EOS,only: rh1,rh,rht1
  use strucmod,only: r,xnabj1,rad1,j1,m,q,p,t,adi1,cap1,Nabla_mu,xnabm,xnabj
  use rotmod,only: omegi,dlodlr
  use geomod, only: rpsi_min,rpsi_max,geocalc
  use nagmod,only: c02agf

  implicit none

  integer::jpos,ifail,ndegre,nroot
  real(kindreal):: xgamj1,gmsu,xgpj1,gravj1,echpj1,diftj1,grazj1,admu,delna,delmu,croch1,delsh,aa0,aa1,aa2,aa3,xgam, &
                   g21gj1,xpsij1

  real(kindreal), dimension(6):: www2
  real(kindreal), dimension(8):: www3
  real(kindreal), dimension(0:2):: apol2
  real(kindreal), dimension(0:3):: apol3
  real(kindreal), dimension(2,2):: zero2
  real(kindreal), dimension(2,3):: zero3

  logical, parameter:: scale=.true.
!-----------------------------------------------------------------------
  jpos=0
  xgamj1=0.d0
  xnabj1=rad1
  gmsu=gms*Msol

! calcul de la gravite
! j1
  if (j1 /= m) then
    xpsij1=exp(r(j1))/((cst_G*gmsu*(1.d0-exp(q(j1))))/(omegi(j1)*omegi(j1)))**(1.d0/3.d0)

    if (xpsij1 >= rpsi_min) then
      if (xpsij1 > rpsi_max) xpsij1=0.9999999999d0*rpsi_max
! ancien appel: call geograv(xpsij1,xgpj1)
      call geocalc(xpsij1,xgpj1,1)
      gravj1=cst_G*gmsu*(1.d0-exp(q(j1)))*omegi(j1)**4.d0
      gravj1=gravj1**(1.d0/3.d0)*xgpj1
    else
      gravj1=cst_G*gmsu*(1.d0-exp(q(j1)))/exp(2.d0*r(j1))
    endif
  else
    gravj1=0.d0
  endif

! Calcule  de l'echelle de pression
  if (j1 /= m) then
    echpj1=exp(p(j1))/(exp(rh1)*gravj1)
  else
    echpj1=exp(p(m-1))/(exp(rh)*cst_G*gmsu*(1.d0-q(m-1))/exp(2.d0*r(m-1)))
  endif
! du coefficient de diffusion radiative
  diftj1=4.d0*cst_a*cst_c*exp(4.d0*t(j1))*adi1/(3.d0*exp(cap1+rh1+p(j1))*(-rht1))
! du gradient de vitesse verticale
! 9 pi/32: resultat de int_0^pi (sin^4(theta)dtheta)
!                     /int_0^pi (sin^3(theta)dtheta)
  grazj1=abs((9.d0*pi/32.d0)*omegi(j1)*dlodlr(j1))

! critere de LEDOUX ou de SCHWARZSCHILD
  if (iledou == 1) then
    admu=adi1+Nabla_mu(j1)/(-rht1)
  else
    admu=adi1
  endif

! ZONES RADIATIVES
  if (rad1 < admu) then
    delna=adi1-rad1
    delmu=Nabla_mu(j1)/(-rht1)
    croch1=grazj1*grazj1
    if (gravj1 /= 0.0d0) then
      delsh=echpj1/(gravj1*(-rht1))*croch1
    else
      delsh=0.d0
    endif
    aa0=12.d0*delmu
    aa1=2.d0*(delna+delmu)-6.d0*delsh
    aa2=2.d0*(delna+delmu)+4.d0*delna-delsh
    aa3=-delsh
    ifail=0
    if (aa0 /= 0.d0) then
      ndegre=3
      apol3(0)=aa0
      apol3(1)=aa1
      apol3(2)=aa2
      apol3(3)=aa3
      call c02agf(apol3,ndegre,scale,zero3,www3,ifail)
      nroot=1
      do while (nroot <= ndegre)
       if (zero3(2,nroot) == 0.d0) then
         if (zero3(1,nroot) > 0.d0) then
           xgam=zero3(1,nroot)
           xgamj1=max(xgamj1,xgam)
           jpos=jpos+1
           nroot=nroot+1
         else
           nroot=nroot+1
         endif
       else
         nroot=nroot+1
       endif
      enddo
    else
      if (aa1 /= 0.d0) then
        ndegre=2
        apol2(0)=aa1
        apol2(1)=aa2
        apol2(2)=aa3
        call c02agf(apol2,ndegre,scale,zero2,www2,ifail)
        nroot=1
        do while (nroot <= ndegre)
         if (zero2(2,nroot) == 0.d0) then
           if (zero2(1,nroot) > 0.d0) then
             xgam=zero2(1,nroot)
             xgamj1=max(xgamj1,xgam)
             jpos=jpos+1
             nroot=nroot+1
           else
             nroot=nroot+1
           endif
         else
           nroot=nroot+1
         endif
        enddo
      else
        if (aa2 /= 0.d0) then
          xgamj1=delsh/aa2
          if (xgamj1 >= 0.d0) jpos=jpos+1
        else
          jpos=-1
        endif
      endif
    endif
  endif

  if (xgamj1 > 0.d0) then
    g21gj1=(6.d0*xgamj1*xgamj1)/(1.d0+xgamj1)
    xnabj1=(rad1+g21gj1*adi1)/(1.d0+g21gj1)
  endif
  xnabm=0.5d0*(xnabj+xnabj1)

  return

end subroutine nabgam

!=======================================================================
subroutine grapmui
!-----------------------------------------------------------------------
! calcule du gradient de mu par rapport a la pression
!-----------------------------------------------------------------------
  use inputparam,only: z,ialflu,itminc,idialu,nwseq,verbose
  use abundmod,only: x,y3,y,xc12,xc13,xc14,xn14,xn15,xo16,xo17,xo18,xf18,xf19,xne20,xne21,xne22,xna23,xmg24,xmg25,xmg26,xal26, &
    xal27,xsi28,xprot,xneut,xbid,xbid1,wx,wy3,wy,wxc12,wxc13,wxc14,wxn14,wxn15,wxo16,wxo17,wxo18,wxf18,wxf19,wxne20,wxne21, &
    wxne22,wxna23,wxmg24,wxmg25,wxmg26,wxal26g,wxal27,wxsi28,wxprot,wxneut,wxbid,wxbid1,zabelx,nbelx,wabelx,abelx,nbzel,nbael
  use strucmod,only: m,Nabla_mu,p,t,r,s,vr,zensi,pb,tb,rb,sb,xb,amu,amub,xmufit
  use SmallFunc,only: CheckProfile,SmoothProfile

  implicit none

  integer:: n,nm,ii,i,jmu,jnul,nnul,izin,izex,jord,jmu0,kn,jinte,ii11,jexte,ii22,jor0,jmu1=0,ideb,ifin,ncou,l,refit
  integer, dimension(300):: gmu0,gmu1,xzmu
  real(kindreal):: xmoin1,xplus1,ymoin1,yiiiii,yplus1,xten,cconv,econv,chimu2,distmax,chimutol,amumax
  real(kindreal), dimension(10):: aaa
  real(kindreal), dimension(ldi):: press,xmul,ray,xmu_dl,wpmu,rln,rdif,xmu2,amu1,amu2,amu3,wpent,wpentini,amulisse,wpent_inv
!-----------------------------------------------------------------------
  chimu2 = 0.d0
  distmax = 0.d0
  amumax = 0.d0
  chimutol = 0.d0
  refit = 1

  xzmu(:)=0
  gmu0(:)=0
  gmu1(:)=0
  aaa(:)=0.0d0
  amu1(:)=0.d0
  amu2(:)=0.d0
  amu3(:)=0.d0
  amulisse(:)=0.d0

! On sauve dans les variables les grandeurs
! utiles au calcul de l'advection. On doit par coherence
! calculer l'advection en utilisant les variables de
! l'iteration precedente
  pb(:)=p(:)
  tb(:)=t(:)
  rb(:)=r(:)
  sb(:)=s(:)
  xb(:)=x(:)

!*******************************************************
! renversement de la numerotation et calcul de mu
  do n=1,m
   nm=m-n+1
   press(n)=exp(p(nm))
   wx(n)=x(nm)
   wy3(n)=y3(nm)
   wy(n)=y(nm)
   wxc12(n)=xc12(nm)
   wxc13(n)=xc13(nm)
   wxn14(n)=xn14(nm)
   wxn15(n)=xn15(nm)
   wxo16(n)=xo16(nm)
   wxo17(n)=xo17(nm)
   wxo18(n)=xo18(nm)
   wxne20(n)=xne20(nm)
   wxne22(n)=xne22(nm)
   wxmg24(n)=xmg24(nm)
   wxmg25(n)=xmg25(nm)
   wxmg26(n)=xmg26(nm)
   amu1(n)=2.d0*wx(n)+wy3(n)+0.75d0*wy(n)+7.d0*(wxc12(n)/12.d0+wxc13(n)/13.d0)+8.d0*(wxn14(n)/14.d0+wxn15(n)/15.d0) &
            +9.d0*(wxo16(n)/16.d0+wxo17(n)/17.d0)+13.d0*(wxmg24(n)/24.d0+wxmg25(n)/25.d0)+11.d0/20.d0*wxne20(n) &
            +0.5d0*(wxo18(n)+wxne22(n)+wxmg26(n))
   if (ialflu == 1) then
     wxf19(n)=xf19(nm)
     wxne21(n)=xne21(nm)
     wxna23(n)=xna23(nm)
     wxal26g(n)=xal26(nm)
     wxal27(n)=xal27(nm)
     wxsi28(n)=xsi28(nm)
     wxneut(n)=xneut(nm)
     wxprot(n)=xprot(nm)
     wxc14(n)=xc14(nm)
     wxf18(n)=xf18(nm)
     wxbid(n)=xbid(nm)
     wxbid1(n)=xbid1(nm)
     amu2(n)=10.d0/19.d0*wxf19(n)+11.d0/21.d0*wxne21(n)+14.d0/26.d0*wxal26g(n)+14.d0/27.d0*wxal27(n)+15.d0/28.d0*wxsi28(n) &
              +12.d0/23.d0*wxna23(n)+10.d0/18.d0*wxf18(n)+0.50d0*(wxc14(n)+wxbid(n)+wxbid1(n))
   endif
   amu3(n)=0.5d0*zabelx
   do ii=1,nbelx
    wabelx(ii,n)=abelx(ii,nm)
    amu3(n)=amu3(n)+wabelx(ii,n)*(nbzel(ii)+1.)/nbael(ii)
   enddo
   amu(n)=amu1(n)+amu2(n)+amu3(n)
  enddo
  call CheckProfile(amu,amulisse,press,m)

! sauvetage du poinds moleculaire moyen pour le calcul
! de l'advection
  amub(1:m) = 1.d0/amulisse(m:1:-1)

! premiere estimation du gradient de mu et detection des zones ou
! gradient de mu est nul
  do i=2,m-1
   nm=m-i+1
! on ne met un gradient non nul que lorsque zensi negatif !
   if (zensi(nm) <= 0.d0) then
     xmoin1=log(press(i-1))-log(press(i))
     xplus1=log(press(i+1))-log(press(i))
     ymoin1=-log(amulisse(i-1))
     yiiiii=-log(amulisse(i))
     yplus1=-log(amulisse(i+1))
     wpent(i)=(xplus1*xplus1*(ymoin1-yiiiii)+xmoin1*xmoin1*(yiiiii-yplus1))/(xmoin1*xplus1*(xplus1-xmoin1))
! on prend ici la valeur absolue car on veut eliminer
! l'effet de gradients negatifs qui peuvent apparaitre
! marginalement au cours de l'evolution, qui ne prennent
! semble-t-il jamais une amplitude telle qu'il faille en
! tenir compte
     if (wpent(i) < 0.d0) wpent(i)=abs(wpent(i))
   else
     wpent(i)=0.d0
   endif
   if (abs(wpent(i)) <= 1.0d-8) wpent(i)=0.d0
  enddo
  wpent(1)=wpent(2)
  wpent(m)=wpent(m-1)
  wpentini=wpent
! On comparera le chi^2 a la valeur maximale de mu, c-a-d au centre
  amumax=abs(-log(amulisse(1)))
!*******************************************************
! impression pour controle
  rdif (1)=0.d0
  do i=1,m
   xmul(i)=-log(amulisse(i))
   rln(i)=log(press(i))
   if (i > 1) rdif (i)=rln(i)-rln(i-1)
   if (i > 1 .and. rdif (i) >= 0) then
     if (verbose .or. itminc == 1) then
       print*,' grapmui, rdif',i
     endif
   endif
  enddo

! fite le poids moleculaire moyen
! avec un polynome si xdial different de zero
  if (idialu == 0) then
    do n=1,m
     if (wpent(n) < 0.d0) wpent(n)=abs(wpent(n))
     l=m-n+1
     Nabla_mu(l)=wpent(n)
    enddo
    return   !   <== EXIT
  else ! idialu=1
! detection des zones ou le gradient de mu est nul
! on va du centre vers le bord
    jmu=1
! nombre de zones avec gradients de mu=zero
    jnul=0
! compteur du nombre de couches avec gradient de mu=zero
! mis egal a 1 lors de premiere couche a mu nul.
 3  nnul=0
 1  if (wpent(jmu) <= 1.0d-8) then
      nnul=nnul+1
      if (jmu /= m) then
        jmu=jmu+1
        go to 1
      else
        go to 2
      endif
    else
      if (nnul == 0) then
        if (jmu /= m) then
          jmu=jmu+1
          go to 1
        else
          go to 2
        endif
      else
        izin=jmu-nnul
        izex=jmu
        jnul=jnul+1
        jord=2*jnul-1
        xzmu(jord)=izin
        xzmu(jord+1)=izex
        if (jmu /= m) then
          jmu=jmu+1
          go to 3
        else
          go to 2
        endif
      endif
    endif
2   if (nnul /= 0) then
      izin=jmu-nnul+1
      izex=jmu
      jnul=jnul+1
      jord=2*jnul-1
      xzmu(jord)=izin
      xzmu(jord+1)=izex
    endif

! le vecteur xzmu contient le numero de la premiere couche avec
! grad mu zero et le numero de la premiere couche ou apres une serie
! de grad mu zero on a grad mu non zero. Cas zone convective exterieure
! numero de la couche au bord --> 1.

! Les zones convectives qui couvrent une extension qui est inferieure
! a un pourcent du rayon en fitm sont consideree comme non existante
! dans le calcul du gradient de mu.

! calcul de l'etendue de la zone convective. Les zones convectives
! dont l'etendue est inferieure au pourcent du rayon de la couche 1
! ne sont pas consideree comme des zones convectives pour le calcul
! du gradient de mu
    jmu0=0
    do kn=1,jnul
     jinte=2*kn-1
     ii11=xzmu(jinte)
     jexte=2*kn
     ii22=xzmu(jexte)
     xten=(exp(vr(m-ii11+1))-exp(vr(m-ii22+1)))/(-exp(vr(1)))

! on met les bornes de cette zone convective
! dans le vecteur gmu0
     if (xten >= 0.0) then
       jmu0=jmu0+1
       jor0=2*jmu0-1
       gmu0(jor0)=ii11
       gmu0(jor0+1)=ii22
     endif
    enddo

! determination des zones ou le gradient de mu doit etre recalcule
! 1) combien a-t-on de zones ou gradient de mu doit
! etre recalcule ? Si jmu0 est le nombre de zones
! ou le gradient de mu est nul, alors le nombre de
! zones ou le gradient de mu est non nul est compris
! entre jmu0-1 et jmu0+1 (nous suppososons jmu0 sup ou egal a 1
    cconv=0.0d0
    econv=0.0d0
    if (gmu0(1) == 1) cconv=1
    if (jmu0 /= 0) then
      if (gmu0(2*jmu0) == m) econv=1.d0
    endif

    if (cconv == 1.d0 .and. econv == 1.d0) then

      jmu1=jmu0-1
! etoile avec gradient de mu=0 partout
      if (jmu1 == 0) then
        do n=1,m
         if (wpent(n) < 0.d0) wpent(n)=abs(wpent(n))
         l=m-n+1
         Nabla_mu(l)=wpent(n)
        enddo
        return   !   <== EXIT
      endif
      do kn=1,jmu1
       jinte=2*kn-1
       jexte=2*kn
       gmu1(jinte)=gmu0(jexte)
       gmu1(jexte)=gmu0(jexte+1)-1
      enddo

    endif

    if (cconv == 1.d0 .and. econv == 0.d0) then

      jmu1=jmu0
      do kn=1,jmu1
       jinte=2*kn-1
       jexte=2*kn
       gmu1(jinte)=gmu0(jexte)
       if (kn /= jmu1) then
         gmu1(jexte)=gmu0(jexte+1)-1
       else
         gmu1(jexte)=m
       endif
      enddo

    endif

    if (cconv == 0.d0 .and. econv == 1.d0) then

      jmu1=jmu0
      do kn=1,jmu1
       jinte=2*kn-1
       jexte=2*kn
       if (kn == 1) then
         gmu1(jinte)=1
       else
         gmu1(jinte)=gmu0(jexte-2)
       endif
       gmu1(jexte)=gmu0(jinte)-1
      enddo

    endif

    if (cconv == 0.d0 .and. econv == 0.d0) then
      jmu1=jmu0+1
      do kn=1,jmu1
       jinte=2*kn-1
       jexte=2*kn
       if (kn == 1) then
         gmu1(jinte)=1
       else
         gmu1(jinte)=gmu0(jexte-2)
       endif
       if (kn /= jmu1) then
         gmu1(jexte)=gmu0(jinte)-1
       else
         gmu1(jexte)=m
       endif
      enddo

    endif

    ray(1:m) = rln(m:1:-1)
    xmu_dl(1:m) = xmul(m:1:-1)
    wpent_inv(1:m) = wpent(m:1:-1)
    xmu2 = xmu_dl
    wpmu = wpent_inv

    ncou = 0
    do kn=1,jmu1
     jinte=2*kn-1
     jexte=2*kn
! mise des valeurs du rayon dans le vecteur ray
! mise des valeurs de mu dans le vecteur xmu
     ideb=gmu1(jinte)
     ifin=gmu1(jexte)

     if (jmu0 == 0) then
       ideb=1
       ifin=m
     endif
     ncou = ncou+ifin-ideb
     if (verbose) write(*,*)'Fit on layers',m-ifin+2,m-ideb
     call SmoothProfile(xmu_dl,wpent_inv,xmu2,wpmu,ray,m-ifin+2,m-ideb,5,m)
    enddo
    xmufit = xmu2
    Nabla_mu = abs(wpmu)
    wpent(1:m) = abs(wpmu(m:1:-1))

    chimu2=0.d0
    distmax=0.d0
    do l=1,m
! recherche de la distance max entre fit et courbe + chi^2
     chimu2=chimu2+(xmufit(l)+log(amulisse(m-l+1)))**2.d0
     if (abs(xmufit(l)+log(amulisse(m-l+1)))  >  distmax) then
       distmax=abs(xmufit(l)+log(amulisse(m-l+1)))
     endif
    enddo
    if (ncou  /=  0) then
      chimu2=sqrt(chimu2)/real(ncou)
    else
      chimu2=sqrt(chimu2)
    endif
! on accepte un ecart de 0.005*mu(surface) et un chi^2 < 1.d-5:
    chimutol=1.d-2*abs(-log(amulisse(m)))
    write(3,'(a,1p,2(1x,d23.16))')'GRAPMUI - chi^2,chi^2 accepted:',chimu2,1.d-4*amumax
    write(3,'(a,1p,2(1x,d23.16))')'          distmax,chimutol:',distmax,chimutol

    return
  endif

  return

end subroutine grapmui

!=======================================================================
!***********************************************************************
end module nablas
