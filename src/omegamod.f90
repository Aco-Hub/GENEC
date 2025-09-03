module omegamod

use io_definitions
use evol,only: ldi,kindreal
use const,only: Msol
use caramodele,only: glm,gms
use strucmod,only: m,q,r,rb,zensi
use rotmod,only: vxmoci,xinert,omegi,vsuminenv

implicit none

real(kindreal):: xjspe1,xjspe2

private
public:: VcritCalc
public:: omenew,omescale
public:: dlonew,omconv,vomcon
public:: momevo,momspe
public:: omenex,om2old,omesta
public:: xjspe1,xjspe2

contains
!***********************************************************************
subroutine VcritCalc(ivcalc,vcrit1,vcrit2,vequat)
!------------------------------------------------------------------------
! Calcul des vitesses critiques et des rapports V/Vc et O/Oc
! Calcul de la deformation
!------------------------------------------------------------------------
  use const, only: pi,lgLsol,Lsol,cstlg_sigma,cst_sigma,cstlg_G,lgMsol,lgRsol,cst_G,Msol,Rsol
  use inputparam, only: irot
  use caramodele, only: gls,teff,eddesc
  use rotmod, only: rapcri,rapom2,rapvco,rapvc2,xobla,vpsi
  use geomod, only: sund_max,GammaEddmax_min,geocalc,geomedd

  implicit none

  logical,intent(in):: ivcalc
  real(kindreal),intent(out):: vcrit1,vcrit2,vequat

  real(kindreal):: fffff,ffff3,yrequa,requa2,omegaMax,xoblaMax,rpole,rstar
!------------------------------------------------------------------------
  if (irot == 0) then
    rapcri = 0.0d0
    rapom2 = 0.0d0
    vequat = 0.0d0
! Computation of the stellar radius, and the first critical velocity. As the star is not rotating and is
! spherical, the critical equatorial radius is 1.5 times the radius.
    rstar = sqrt(gls*Lsol/(4.d0*pi*cst_sigma))/teff**2.d0
    vcrit1 = sqrt(cst_G*gms*Msol/(1.5d0*rstar))
    vcrit1=vcrit1/1.0d+05
    vcrit2 = 0.0d0
    rapvco = 0.0d0
    rapvc2 = 0.0d0
    xobla = 1.0d0
    fffff=1.d0
    return
  endif

  if (ivcalc) then
!!! Vpsi est le Spsi' dans la doc !
    vpsi=log10(gls)+lgLsol-4.d0*log10(teff)-cstlg_sigma
    vpsi=10.d0**(vpsi-2.d0/3.d0*(cstlg_G+log10(gms)+lgMsol-2.d0*log10(omegi(1))))
    if (vpsi > sund_max) then
      print*,' omega superieur a omega critique'
      vpsi=0.999999d0*sund_max
    endif
! ancien appel: call geocrit(vpsi,xobla)
    call geocalc(vpsi,xobla,4)
  else
    xobla = 1.0d0
  endif

  if (.not.ivcalc .or. abs(1.d0-xobla) < 1.0d-10) then
    rstar = sqrt(gls*Lsol/(4.d0*pi*cst_sigma))/teff**2.d0
    xobla=1.d0
    if (omegi(1) >= 0.0d0) then
      rapcri=omegi(1)/sqrt(8.d0*cst_G*gms*Msol/(27.d0*rstar**3.d0))
      vequat=omegi(1)*rstar/1.0d+05
    else
      rapcri=0.0d0
      vequat=0.0d0
    endif
    rapom2=rapcri

! Case of small rotation. We assume that the star is spherical, and thus, that the critical equatorial radius
! is 1.5 times the spherical radius.
    vcrit1=sqrt(cst_G*gms*Msol/(1.5d0*rstar))
    vcrit1=vcrit1/1.0d+05
    rapvco=vequat/vcrit1
    vcrit2=0.d0
    rapvc2=rapvco
    fffff=1.d0
    write(*,*) 'vcritcalc, low rot: vequat=',vequat
  else
    fffff=1.d0/xobla
    ffff3=fffff*fffff*fffff
    rapcri=sqrt((3.d0/2.d0)**3.d0*2.d0*(fffff-1.d0)/ffff3)
    if (xobla < 1.d0) then
      vequat=(2.d0*(1.d0/xobla-1.d0))**(1.d0/3.d0)
      vequat=(cst_G*gms*Msol*omegi(1))**(1.d0/3.d0)*vequat
      vequat=vequat/1d+05
      rpole= (2.d0*(1.d0/xobla-1.d0))**(1.d0/3.d0)/fffff
      rpole= rpole*(cst_G*gms*Msol/(omegi(1)*omegi(1)))**(1.d0/3.d0)
    else
      yrequa=1.d0/2.d0*(lgLsol-log10(4.d0*pi)-cstlg_sigma)-lgRsol
      yrequa=10.d0**(yrequa+0.5d0*log10(gls)-2.d0*log10(teff))
      rpole = yrequa
      vequat=(yrequa*Rsol*omegi(1))/1.0d+05
    endif
    rapvco=rapcri*2.d0/(3.d0*xobla)
! calcul de vcrit1
    vcrit1=vequat/rapvco
! calcul de vcrit2
    if (eddesc > GammaEddmax_min) then
      call geomEdd(eddesc,xoblaMax,omegaMax)

! calcul de Omega/Omegacrit2
! rapcri: Omega/Omegacrit1 et omegaMax: Omegacrit2/Omegacrit1
      if (omegaMax /= 0.0d0) then
        rapom2 = rapcri/omegaMax
      else
        write(*,*)'Beware: omegaMax=0... Is there a problem ??'
        rapom2 = 1.0d0
      endif
! calcul de requa a Vcrit2
      requa2=rpole/xoblaMax
      vcrit2=requa2*omegi(1)/(rapom2*1.0d+05)
      rapvc2=vequat/vcrit2
      write(*,*)'2nd critical limit computed'
      write(io_logs,*)'2nd critical limit computed'
    else
      vcrit2=0.0d0
      rapvc2=rapvco
      rapom2=rapcri

      write(*,*)'1st critical limit computed'
      write(io_logs,*)'1st critical limit computed'
    endif
    write(*,*)'         rapport V/Vc= ', rapvc2
    write(io_logs,*)'         rapport V/Vc= ', rapvc2
    write(*,*)'         rapport O/Oc= ', rapom2
    write(io_logs,*)'         rapport O/Oc= ', rapom2
    write(io_logs,'(3(a,f9.3))') '    vequat= ',vequat,' vcrit1= ',vcrit1,' vcrit2= ',vcrit2
    write(*,'(3(a,f9.3))') '    vequat= ',vequat,' vcrit1= ',vcrit1,' vcrit2= ',vcrit2

  endif   !   ivcalc

  return

end subroutine vcritcalc
!=======================================================================
subroutine dlonew
!------------------------------------------------------------------------
! Calcul la derivee par rapport a lnr de ln omega
! pour le nouveau profil de rotation en cours d'iteration
!------------------------------------------------------------------------
  use inputparam,only: itminc,idialo,rapcrilim,xdial,phase,nwseq,verbose
  use caramodele,only: nwmd,printerror
  use strucmod,only: xomegafit,npcoucheeff
  use rotmod,only: dlodlr
  use SmallFunc,only: CheckProfile,SmoothProfile
  use safestop,only: safe_stop

  implicit none

  integer, parameter:: ncoutr=50
  integer:: couchemin,WinSize,refit
  integer::i,k,n,nm,jmu,jnul,nnul,izin,izex,jord,jmu0,kn,ncou,jinte,ii11,jexte,ii22,jor0,jmu1=0,ideb,ifin,l
  integer, dimension(300):: gmu0, gmu1,xzmu
  real(kindreal):: xlolo,yx1,yx3,yx2,yy1,yy2,yy3,wpena,wpenb,wfa,wfb,y1,y2,x1,x2,yyaa,yybb,xten,cconv,econv
  real(kindreal), parameter:: xminex=0.0d0
  real(kindreal):: chiom2,distmax,chiomtol,omegamax
  real(kindreal), dimension(10):: aaa
  real(kindreal), dimension(ldi):: xmul,rln
  real(kindreal), dimension(ldi):: dlokn,xmasol,rayon,xmu2
  real(kindreal), dimension(ldi):: dlodlrini,omegilisse,dlokn_fit

! rln:   ln rayon numerotation centre bord
! xmul:  ln omegi numerotation centre-bord
! dlokn: dln om/d ln r numerot centre-bord
! ray:   ln rayon dans zone de fit
! xmu_dl:   ln omegi dans zone de fit
! wpmu:  dln om/d ln r dans zone de fit

  logical::checkMesh = .true.
!------------------------------------------------------------------------
  chiom2 = 0.d0
  distmax = 0.d0
  omegamax = 0.d0
  chiomtol = 0.d0
  refit = 1

  xzmu=0
  gmu0=0
  gmu1=0
  aaa=0.0d0

  k=m
  checkMesh = .true.

! normalement omega est toujours positif mais il y a des cas ou sur une
! couche ou deux il peut devenir negatif
  do i=2,m
   if (omegi(i) < 0.d0) then
     omegi(i)=max(omegi(i-1),1.d-12)
     if ((verbose .or. itminc == 1) .and. PrintError) then
       print*,' ATTENTION OMEGI NEGATIF, COUCHE=',i
       print*, 'CHANGED to ',omegi(i)
       PrintError = .false.
     endif
   endif
   if (omegi(i) > omegamax) then

     omegamax = omegi(i)
   endif
  enddo
  call CheckProfile(omegi,omegilisse,r,m)

  do i=2,m-1
   xlolo=r(i)-r(i+1)
   if (xlolo <= 0.d0) then

     write(*,'(3(a,i0),a)')' i=',i,': r(',i,') - r(',i+1,') < 0'
     if (itminc == 1) then
       rewind(io_runfile)
       write(io_runfile,'(i7,a,i0)') nwmd,': Problem with radius inversion in layer ',i
       call safe_stop('Problem with radius inversion')
     endif
   endif
  enddo

  do i=2,m-1
   if (zensi(i) <= 0.d0) then
     yx1=r(i)-r(i-1)
     yx3=r(i)-r(i+1)
     yx2=r(i+1)-r(i-1)
     if (omegilisse(i-1) /= 0.0d0) then
       yy1=log(omegilisse(i-1))
     else
       yy1=-50.0d0
     endif
     if (omegilisse(i) /= 0.0d0) then
       yy2=log(omegilisse(i))
     else
       yy2=-50.0d0
     endif
     if (omegilisse(i+1) /= 0.0d0) then
       yy3=log(omegilisse(i+1))
     else
       yy3=-50.0d0
     endif
     wpena=(yy2-yy3)/yx3
     if (yx1 /= 0.d0) then
       wpenb=(yy2-yy1)/yx1
     else
       wpenb=0.d0
     endif
     wfa=yx1/yx2
     wfb=-yx3/yx2
     dlodlr(i)=wpena*wfa+wpenb*wfb
     if (abs(dlodlr(i)) <= 1.0d-08) then
       dlodlr(i)=0.0d0
     endif
   else
     dlodlr(i)=0.0d0
   endif
  enddo

  dlodlr(1)=dlodlr(2)
  y1=dlodlr(k-1)
  y2=dlodlr(k-2)
  x1=r(k-1)
  x2=r(k-2)
  yyaa=(y1-y2)/(x1-x2)
  yybb=(x1*y2-x2*y1)/(x1-x2)
  dlodlr(k)=yyaa*r(k)+yybb
  if (abs(dlodlr(k)) <= 1.0d-08) then
    dlodlr(i)=0.0d0
  endif
  dlodlrini=dlodlr

!-----------------------------------------------------------------------
! renversement de la numerotation
  do n=1,k
   nm=m-n+1
   rln(n)=r(nm)
   rayon(n)=exp(r(nm))
   xmul(n)=log(omegilisse(nm))
   dlokn(n)=dlodlr(nm)
   xmasol(n)=(1.d0-exp(q(nm)))*gms
  enddo
! fite   le profil de omega
! avec un polynome si xdial different de zero
  if (idialo == 0 .and. rapcrilim < 1.d-5) then
    do n=1,k
     l=m-n+1
     dlodlr(l)=dlokn(n)
    enddo
    return   !   <== EXIT
  endif

  if (idialo == 0 .and. rapcrilim > 1.d-5) then
    couchemin = m-(NPcoucheEff+50)
  else
    couchemin = 1
  endif
! detection des zones ou le gradient de mu est nul
! on va du centre vers le bord
  jmu = couchemin
! nombre de zones avec gradients de mu=zero
  jnul=0
! compteur du nombre de couches avec gradient de mu=zero
! mis egal a 1 lors de premiere couche a mu nul.
 3 nnul=0
 1 if (abs(dlokn(jmu)) <= 1.0d-8) then
      nnul=nnul+1
      if (jmu /= k) then
        jmu=jmu+1
        go to 1
      else
        go to 2
      endif
    else
      if (nnul == 0) then
        if (jmu /= k) then
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
        if (jmu /= k) then
          jmu=jmu+1
          go to 3
        else
          go to 2
        endif
      endif
    endif
 2 if (nnul /= 0) then
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

! calcul de l'etendue de la zone convective. Les zones convectives
! dont l'etendue est inferieure au pourmille du rayon de la couche 1
! ne sont pas consideree comme des zones convectives pour le calcul
! du gradient de omega, sauf s'il s'agit du noyau convectif
  jmu0=0
  do kn=1,jnul
   jinte=2*kn-1
   ii11=xzmu(jinte)
   jexte=2*kn
   ii22=xzmu(jexte)
   xten=(rayon(ii11)-rayon(ii22))/(-rayon(m))

! xminex est la fraction minimale du rayon que doit
! couvrir une zone a gradient de omega nul pour qu'elle
! soit consideree comme une region a gradient de omega nul
   if (xten >= xminex .or. ii11 == 1) then
! on met les bornes de cette zone convective
! dans le vecteur gmu0
     jmu0=jmu0+1
     jor0=2*jmu0-1
     gmu0(jor0)=ii11
     gmu0(jor0+1)=ii22
   endif
  enddo

! determination des zones ou le gradient de mu
! doit etre recalcule
! 1) combien a-t-on de zones ou gradient de mu doit
! etre recalcule ? Si jmu0 est le nombre de zones
! ou le gradient de mu est nul, alors le nombre de
! zones ou le gradient de mu est non nul est compris
! entre jmu0-1 et jmu0+1 (nous suppososons jmu0 sup ou egal a 1
  cconv=0.0d0
  econv=0.0d0
  if (gmu0(1) == couchemin) then
    cconv=1.d0
  endif
  if (jmu0 /= 0) then
    if (gmu0(2*jmu0) == k) then
      econv=1.d0
    endif
  endif

  if (cconv == 1.d0 .and. econv == 1.d0) then

    jmu1=jmu0-1
! etoile avec gradient de mu=0 partout
    if (jmu1 == 0) then
      do n=1,k
       l=m-n+1
       dlodlr(l)=dlokn(n)
      enddo
      return   !   <== EXIT
    endif
    do kn=1,jmu1
     jinte=2*kn-1
     jexte=2*kn
     gmu1(jinte)=gmu0(jexte)-1
     gmu1(jexte)=gmu0(jexte+1)
    enddo

  endif

  if (cconv == 1.d0 .and. econv == 0.d0) then

    jmu1=jmu0
    do kn=1,jmu1
     jinte=2*kn-1
     jexte=2*kn
     gmu1(jinte)=gmu0(jexte)-1
     if (kn /= jmu1) then
       gmu1(jexte)=gmu0(jexte+1)
     else
       gmu1(jexte)=k
     endif
    enddo

  endif

  if (cconv == 0.d0 .and. econv == 1.d0) then

   jmu1=jmu0
   do kn=1,jmu1
    jinte=2*kn-1
    jexte=2*kn
    if (kn == 1) then
      gmu1(jinte)=couchemin
    else
      gmu1(jinte)=gmu0(jexte-2)-1
    endif
    gmu1(jexte)=gmu0(jinte)
   enddo

  endif

  if (cconv == 0.d0 .and. econv == 0.d0) then

  jmu1=jmu0+1
   do kn=1,jmu1
   jinte=2*kn-1
   jexte=2*kn
    if (kn == 1) then
     gmu1(jinte)=couchemin
    else
     gmu1(jinte)=gmu0(jexte-2)-1
    endif

    if (kn /= jmu1) then
     gmu1(jexte)=gmu0(jinte)
    else
     gmu1(jexte)=k
    endif
   enddo

  endif

! recherche du fit polynomial pour les differentes regions
  if (xdial == 0.0d0 .and. rapcrilim < 1.d-5) then
   jmu1=1
  endif

  xmu2 = xmul
  dlokn_fit = dlokn
  ncou=0
  do kn=1,jmu1
   jinte=2*kn-1
   jexte=2*kn
   ideb=gmu1(jinte)
   ifin=gmu1(jexte)

   if (jmu0 == 0) then
     ideb = couchemin
     ifin = k
   endif

   ncou=ncou+ifin-ideb
   if (verbose) then
     write(*,*)'Fit on layers',m-ifin+2,m-ideb
   endif
!   if (phase == 1) then
!     WinSize = 5
!   else
!     WinSize = 12
!   endif
! WinSize is defined as a function of m
   WinSize = max(5,ceiling(m / 120.0d0))

! Second condition because if (ifin-ideb-2) < 2*WinSize+1 the smoothing is not performed (voir smallfunc)
   if ((ifin-ideb-2) < 4*WinSize+1 .and. (ifin-ideb-2) >= 2*WinSize+1) then
     write(io_logs,*) 'omega profile smoothed on a relatively small zone: window = ',2*WinSize+1,' zone to fit : ',&
     ifin-ideb-2,' in shells.'
   endif
   do n=ideb+1,ifin-1
    if (rln(n+1)  <  rln(n)) then
      checkMesh = .false.
    endif
   enddo
   if (checkMesh .or. itminc  ==  1) then
     call SmoothProfile(xmul,xmu2,dlokn_fit,rln,ideb+1,ifin-1,WinSize,m)
   endif
  enddo
  xomegafit(1:m) = xmu2(m:1:-1)
  dlokn = dlokn_fit
  dlodlr(1:m)=dlokn_fit(m:1:-1)

  chiom2=0.d0
  distmax=0.d0
  do n=1,k
! recherche de la distance max entre fit et courbe + chi^2
   chiom2=chiom2+(xomegafit(n)-log(omegilisse(n)))**2.d0
   if (abs(xomegafit(n)-log(omegilisse(n)))  >  distmax) then
     distmax=abs(xomegafit(n)-log(omegilisse(n)))
   endif
  enddo
  if (ncou  /=  0) then
    chiom2=sqrt(chiom2)/real(ncou)
  else
    chiom2=sqrt(chiom2)
  endif
  chiomtol=1.d-3*abs((log(omegilisse(m-couchemin+1))+log(omegilisse(1)))/2.d0)
  write(io_logs,'(a,1p,2(1x,d23.16))')'DLOGRA - chi^2,chi^2 accepted:',chiom2,1.d-4*omegamax
  write(io_logs,'(a,1p,2(1x,d23.16))')'         distmax,chiomtol:',distmax,chiomtol

  return

end subroutine dlonew

!=======================================================================
subroutine omenew
! Calcul le nouveau profil de rotation
!------------------------------------------------------------------------
  use caramodele,only: xLstarbefHen
  use rotmod,only: omegp

  implicit none

  integer:: jo
  real(kindreal), dimension(ldi):: xmr
!------------------------------------------------------------------------
  xmr=0.0d0
  do jo=1,m
   xmr(jo)=exp(glm)*(1.d0-exp(q(jo)))
  enddo

! calcul du moment d'inertie de chaque coquille du modele courant
! vxmoci contient le moment cinetique du modele precedent
  do jo=2,m-1
   xinert(jo)=2.d0/3.d0*exp(2.0d0*rb(jo))*(xmr(jo-1)-xmr(jo+1))/2.d0
  enddo
  xinert(1)=2.d0/3.d0*exp(2.0d0*rb(1))*(xmr(1)-xmr(2))/2.d0+vsuminenv
  xinert(m)=2.d0/5.d0*exp(2.0d0*rb(m-1))/2.d0**(2.d0/3.d0)*xmr(m-1)/2.d0

! calcul du nouvel omega
  do jo=1,m
   omegp(jo)=vxmoci(jo)/xinert(jo)
  enddo

! calcul du moment angulaire total
  xLstarbefHen=0.d0
  do jo=1,m
   xLstarbefHen=xLstarbefHen+xinert(jo)*omegp(jo)
  enddo

  return

end subroutine omenew

!=======================================================================
subroutine omescale
!------------------------------------------------------------------------
!  La conservation du moment cinetique avec la perte de masse est assuree
! avec les valeur rb et omegd, mais pas avec r et omegi. Afin de corriger
! les petites erreurs, on utilise les valeurs de moments cinetiques de chaque
! coquille calculees dans advect, et on les utilise ici pour recalculer
! le profil de rotation omegi.
!------------------------------------------------------------------------
  use evol,only: npondcouche
  use inputparam,only: verbose,rapcrilim,idifcon
  use caramodele,only: xLtotbeg
  use strucmod,only: NPcoucheEff
  use rotmod,only: dlelex,CorrOmega

  implicit none

  integer::jo=0,k=0,n=0,i1=0,i2=0,i=0
  real(kindreal):: btota=0.0d0,sumine=0.0d0,sumcin=0.0d0,dmine=0.0d0,dmcine=0.0d0,dmin1=0.0d0,dmcin1=0.0d0,omegm=0.0d0, &
                   xLdiff=0.0d0,StillToLose=0.0d0

! poids: tableau contenant les poids de chaque couche (ici, lineaire,
! avec poids de la (nponcouche + 1)eme couche = 0).
  integer, parameter:: NZero = 0
  integer:: ProfZC=0,ProfCor=0
  real(kindreal):: qdenom=0.0d0,qcorr=0.0d0,Lcorr=0.d0
  real(kindreal), dimension (npondcouche):: Li,poids
  real(kindreal), dimension(ldi):: omegacorr,xmr
!------------------------------------------------------------------------
! Initialisation des poids:
  do jo = 1,NPcoucheEff
   poids(jo) = 1.d0/NPcoucheEff
  enddo

  do jo=1,m
   xmr(jo)=exp(glm)*(1.d0-exp(q(jo)))
  enddo

! calcul du moment d'inertie de chaque coquille du modele courant
  btota=0.d0
  do jo=2,m-1
   xinert(jo)=2.d0/3.d0*exp(2.d0*r(jo))*(xmr(jo-1)-xmr(jo+1))/2.d0
   btota=btota+xinert(jo)*omegi(jo)
  enddo
  xinert(1)=2.d0/3.d0*exp(2.d0*r(1))*(xmr(1)-xmr(2))/2.d0+vsuminenv
  xinert(m)=2.d0/5.d0*exp(2.d0*r(m-1))/2.d0**(2.d0/3.d0)*xmr(m-1)/2.d0
  btota=btota+xinert(1)*omegi(1)+xinert(m)*omegi(m)

! Calcul du moment cinetique total perdu par l'etoile depuis le debut du
! calcul du modele en cours, a savoir: pertes intrinseque liees au
! changement de la masse de l'etoile.
  write(io_logs,*)
  write(io_logs,*) 'Dans omescale : '
  write(io_logs,*) 'L initial = ', xltotbeg,'    L apres conv =',btota
  xLdiff = xltotbeg - btota
  write(io_logs,*) 'Perdu "naturellement" par le code = ', xLdiff
  write(io_logs,*) 'A perdre effectivement = ', dlelex
! Encore a perdre / regagner:
! [Modif Conserve]
! Pour ancienne conservation (tout sur NPcoucheEFF):
  StillToLose = dlelex - xLdiff
! Pour nouvelle conservation (variation purement numerique sur
!       l'ensemble de l'etoile et la variation due a la perte de masse
!       sur NPcoucheEFF):
!      StillToLose = dlelex
! [/Modif Conserve]
  if (verbose) then
    write(*,*) 'Lini, Lact: ', xltotbeg,btota
    write(*,*) 'Actuellement perdu: ', xLdiff
    write(*,*) 'A perdre: ', dlelex
    write(*,*) 'Encore a perdre: ', StillToLose
    write(*,*) 'A atteindre: ', xltotbeg-dlelex
  endif

! [Modif Conserve]
! On corrige sur l'ensemble de l'etoile la variation purement numerique
! Pour ancienne conservation: enlever cette boucle
!      do jo=1,m
!        omegi(jo) = omegi(jo)*xltotbeg/btota
!      enddo
! [/Modif Conserve]

! Recherche de la profondeur de la zone convective superficielle.
  ProfZC=0
  do jo=1,m
    if (zensi(jo) > 0.d0) then
      ProfZC=jo
    else
      exit
    endif
  enddo

  if (zensi(1) > 0.d0 .and. ProfZC > 25) then
    ProfCor=min(ProfZC,NPcoucheEff)
  else
    ProfCor=NPcoucheEff
  endif

! Calcul du facteur global de correction q (cf. documentation) donne par
! q = -StillToLose / (sum_i Li pi + xLatmos p1)
  qdenom = 0.d0
  Lcorr = 0.d0
  do jo=1,ProfCor
   Li(jo) = xinert(jo)*omegi(jo)
   qdenom = qdenom + Li(jo)*poids(jo)
   Lcorr = Lcorr + Li(jo)
  enddo

! Si le contenu en moment cinetique des NPcoucheEff couches est suffisant,
! la correction est calculee de maniere standard (avec correction
! sur ces couches, et application au debut du modele suivant.
!  if (Lcorr >= abs(StillToLose)) then
    qcorr = -StillToLose / qdenom

! Pour chaque couche, la correction est donnee par:
! CorrOmega(i) = omega(i)*p(i)*q
    do jo=1,ProfCor
     CorrOmega(jo) = omegi(jo)*poids(jo)*qcorr
    enddo
    do jo=ProfCor+1,NPcoucheEff
     CorrOmega(jo) = 0.d0
    enddo
! Si rapcrilim est nul, la correction n'est plus appliquee.
! On la fixe donc a 0.
    if (rapcrilim <= 1.d-5) then
      CorrOmega=0.d0
    endif
    if (verbose) then
      write(*,'(a,2(1x,d14.8))') 'omega, correction sur Omega(1):',omegi(1),CorrOmega(1)
    endif

! Si il y a des zones convectives a la surface, il faut faire attention
! a la correction. Ici, on applique la moyenne sur les zones convectives
! afin de corriger la correction dans celles-ci.
    do jo=1,NPcoucheEff
     omegacorr(jo) = omegi(jo) + CorrOmega(jo)
    enddo
    do jo = NPcoucheEff + 1,m
     omegacorr(jo) = omegi(jo)
    enddo

! calcul du moment angulaire total
    btota=0.d0
    do jo=1,m
     btota=btota+xinert(jo)*omegacorr(jo)
    enddo

! La correction prend en compte l'enveloppe. Il faut donc maintenant aussi en prendre compte.
!    xinert(1) = xinert(1) + xIatmos

    if (idifcon /= 1) then
      k=m-1
      couche2 : do while (k >= 0)
       n=0
       sumine=0.d0
       sumcin=0.d0

       do while (zensi(k) > 0.d0)
        n=n+1
        if (n == 1 .and. k < m-1) then
          dmine=xinert(k+1)
          dmcine=omegacorr(k+1)*dmine
          sumine=sumine+dmine
          sumcin=sumcin+dmcine
        endif
        if (k == m-1) then
          dmine=xinert(k+1)
          dmcine=omegacorr(k+1)*dmine
          sumine=sumine+dmine
          sumcin=sumcin+dmcine
        endif
        dmin1=xinert(k)
        dmcin1=omegacorr(k)*dmin1
        sumine=sumine+dmin1
        sumcin=sumcin+dmcin1

        if (k-1 < 0) exit couche2
        if (k-1 == 0) then
          k=0
          exit
        endif
        k=k-1
       enddo
       if (n > 0) then
         omegm=sumcin/sumine
         i1=k+1
         i2=k+n+1
         do i=i1,i2
          omegacorr(i)=omegm
         enddo
       endif
       if (k-1 <= 0) exit
       k=k-1
      enddo couche2
    endif ! idifcon

! Une fois les valeurs moyennees sur la zone convective, on calcule la correction finale.
    do jo = 1,NPcoucheEff
     CorrOmega(jo) = omegacorr(jo) - omegi(jo)
    enddo
    do jo=NPcoucheEff+1,m
      omegi(jo) = omegacorr(jo)
    enddo

  write(io_logs,'(a,d14.8)')'Omega(1) apres melange convectif = ',omegi(1)
  write(io_logs,'(a,d14.8)')'Correction calculee pour Omega(1) = ',CorrOmega(1)

  return

end subroutine omescale
!=======================================================================
subroutine omconv
!------------------------------------------------------------------------
! Calcul le nouveau profil de rotation dans les zones convectives
!------------------------------------------------------------------------
  implicit none

  integer:: jo,k,n,i1,i2,i
  real(kindreal):: sumine,sumcin,dmine,dmcine,dmin1,dmcin1,omegm

  real(kindreal), dimension(ldi):: xmr
!------------------------------------------------------------------------
  xmr=0.0d0

  do jo=1,m
   xmr(jo)=exp(glm)*(1.d0-exp(q(jo)))
  enddo

! calcul du moment d'inertie de chaque coquille du modele courant
! vxmoci contient le moment cinetique du modele precedent
  do jo=2,m-1
   xinert(jo)=2.d0/3.d0*exp(2.d0*r(jo))*(xmr(jo-1)-xmr(jo+1))/2.d0
  enddo
  xinert(1)=2.d0/3.d0*exp(2.d0*r(1))*(xmr(1)-xmr(2))/2.d0+vsuminenv
  xinert(m)=2.d0/5.d0*exp(2.d0*r(m-1))/(2.d0**(2.d0/3.d0))*xmr(m-1)/2.d0

  k=m-1
  do while (k > 0)
   n=0
   sumine=0.d0
   sumcin=0.d0
   do while (zensi(k) > 0.d0)
    n=n+1
    if (n == 1 .and. k < m-1) then
      dmine=xinert(k+1)
      dmcine=omegi(k+1)*dmine
      sumine=sumine+dmine
      sumcin=sumcin+dmcine
    endif
    if (k == m-1) then
      dmine=xinert(k+1)
      dmcine=omegi(k+1)*dmine
      sumine=sumine+dmine
      sumcin=sumcin+dmcine
    endif
    dmin1=xinert(k)
    dmcin1=omegi(k)*dmin1
    sumine=sumine+dmin1
    sumcin=sumcin+dmcin1
    if (k < 1) then
      return
    else if (k == 1) then
      k=0
      exit
    endif
    k=k-1
   enddo   ! do while (zensi(k))

   if (n > 0) then

! c'est ici qu'il faut rajouter l'enveloppe
     omegm=sumcin/sumine
     i1=k+1
     i2=k+n+1
     do i=i1,i2
      omegi(i)=omegm
     enddo
   endif   ! if (n)
   k=k-1
  enddo   ! do while (k)

  return

end subroutine omconv

!=======================================================================
subroutine momevo(r,xoread,xLtot,Corr,mode)
!------------------------------------------------------------------------
! sous-routine de calcul du moment angulaire
!------------------------------------------------------------------------
  use inputparam,only: verbose
  use const,only: Msol
  use rotmod,only: bmomin,btotq,bmomit,btot,btotatm
  use evol, only: npondcouche
  use strucmod, only: NPcoucheEff

  implicit none

  real(kindreal),dimension(npondcouche),intent(in):: Corr
  real(kindreal),dimension(ldi),intent(in):: r,xoread
  logical,intent(in):: mode
  real(kindreal),intent(out):: xLtot

  integer:: i,nm,imb,j
  real(kindreal):: xdm1,btoti,btoto,xLatmCG
  real(kindreal), dimension(ldi):: xq,xr,xo,btotr,bmomr
!------------------------------------------------------------------------
! Ajout de la correction a la vitesse pour les NPcoucheEff premieres couches.
  xo = xoread
  do i=1,NPcoucheEff
   xo(i) = xo(i) + Corr(i)
  enddo

  nm=m
  do i=1,m
   xq(i)=(1.d0-exp(q(i)))*gms*Msol
   if (i  /=  m) then
     xr(i)=exp(2.0d0*r(i))
   else
     xr(i)=exp(2.d0*r(i-1))/2.d0**(2.d0/3.d0)
   endif
  enddo

! calcul du moment d'inertie de chaque coquille
  bmomin(1)=2.d0/3.d0*xr(1)*(xq(1)-xq(2))/2.d0
  do imb=2,nm-1
   xdm1=(xq(imb-1)-xq(imb+1))/2.d0
   bmomin(imb)=2.d0/3.d0*xr(imb)*xdm1
  enddo
  bmomin(nm)=2.d0/5.d0*xr(nm)*xq(nm-1)/2.d0

! calcul du moment angulaire de chaque coquille
  btotq(1)=bmomin(1)*xo(1)
  do imb=2,nm-1
   btotq(imb)=bmomin(imb)*xo(imb)
  enddo
  btotq(nm)=bmomin(nm)*xo(nm)
! calcul du moment d inertie de chaque sphere de rayon r
  btoti=0.d0
  do j=1,nm
   btoti=btoti+bmomin(nm-j+1)
   bmomr(nm-j+1)=btoti
  enddo
! calcul du moment d inertie total
  bmomit=bmomr(1)
! calcul du moment angulaire de chaque sphere de rayon r
  btoto=0.d0
  do j=1,nm
   btoto=btoto+btotq(nm-j+1)
   btotr(nm-j+1)=btoto
  enddo
! calcul du moment angulaire total
  btot=btotr(1)

  xLatmCG = vsuminenv*xo(1)
  btotatm=btot+xLatmCG

  if (.not.mode) then
    xltot = btotatm
  else
    xltot = btot/1.d53
  endif

  if (verbose) then
    write(*,'(a,2(d14.8,2x))') 'Momevo: omega, correction sur Omega(1): ',xo(1),Corr(1)
  endif
  return

end subroutine momevo
!=======================================================================
subroutine vomcon
!------------------------------------------------------------------------
! On homogeneise vomegi dans les nouvelles zones convectives avant
! d'appliquer conserv.locale, diffusion et advection
!------------------------------------------------------------------------
  use const,only: Rsol
  use strucmod,only: vr
  use inputparam,only: idifcon
  use rotmod,only: vvomeg,vomegi,vvsuminenv

  implicit none

  integer:: jo,k,n,i,i1,i2
  real(kindreal):: sumine,sumcin,dmine,dmcine,dmin1,dmcin1,omegm
  real(kindreal), dimension(ldi):: xmr
!------------------------------------------------------------------------
  do jo=1,m
   xmr(jo)=exp(glm)*(1.d0-exp(q(jo)))
   vvomeg(jo)=vomegi(jo)
  enddo

! calcul du moment d'inertie de chaque coquille du modele courant
  do jo=2,m-1
   xinert(jo)=2.d0/3.d0*exp(2.d0*vr(jo))*(xmr(jo-1)-xmr(jo+1))/2.d0
  enddo
  xinert(1)=2.d0/3.d0*exp(2.d0*vr(1))*(xmr(1)-xmr(2))/2.d0  + vvsuminenv
  xinert(m)=2.d0/5.d0*exp(2.d0*vr(m-1))/(2.d0**(2.d0/3.d0))*xmr(m-1)/2.d0

  if (idifcon /= 1) then
    k=m-1
    couche : do while (k >= 0)
     n=0
     sumine=0.d0
     sumcin=0.d0

     do while (zensi(k) > 0)
      n=n+1
      if (n == 1 .and. k < m-1) then
        dmine=xinert(k+1)
        dmcine=vomegi(k+1)*dmine
        sumine=sumine+dmine
        sumcin=sumcin+dmcine
      endif
      if (k == m-1) then
        dmine=xinert(k+1)
        dmcine=vomegi(k+1)*dmine
        sumine=sumine+dmine
        sumcin=sumcin+dmcine
      endif
      dmin1=xinert(k)
      dmcin1=vomegi(k)*dmin1
      sumine=sumine+dmin1
      sumcin=sumcin+dmcin1

      if (k-1 < 0) exit couche
      if (k-1 == 0) then
        k=0
        exit
      endif
      k=k-1
     enddo
     if (n > 0) then

       omegm=sumcin/sumine
       i1=k+1
       i2=k+n+1
       do i=i1,i2
        vvomeg(i)=omegm
       enddo
     endif
     if (k-1 <= 0) exit
     k=k-1
    enddo couche
  endif ! idifcon
! vxmoci contient le moment cinetique du modele precedent
  do jo=1,m
   vxmoci(jo)=xinert(jo)*vvomeg(jo)
  enddo

  return

end subroutine vomcon

!=======================================================================
subroutine momspe(o,xjspe1,xjspe2,gms)
!------------------------------------------------------------------------
! sous-routine de calcul du moment specifique (omega r2)
! aux masses xma1 et xma2
!------------------------------------------------------------------------
  use const,only: Rsol
  use interpolation,only: fipoi

  implicit none

  real(kindreal),intent(in):: gms
  real(kindreal),intent(in), dimension(ldi):: o
  real(kindreal),intent(out):: xjspe1,xjspe2

  integer:: i
  real(kindreal), parameter:: xma1=3.d0,xma2=5.d0
  real(kindreal):: xmlog,xrm,xom
  real(kindreal), dimension(ldi):: xq,xr,xo

  do i=1,m
   if (i == m) then
     xq(i) = -15.d0
   else
     xq(i)=log10((1.d0-exp(q(i)))*gms)
   endif
   xr(i)=log10(exp(r(i))/Rsol)
   xo(i)=log10(max(o(i), 1.d-50))
  enddo
  xmlog=log10(xma1)
  xrm=10.d0**fipoi(xmlog,m,xq,xr)
  xom=10.d0**fipoi(xmlog,m,xq,xo)

  xjspe1=(xom*xrm*xrm)*Rsol *Rsol/1.d16
  xmlog=log10(xma2)
  xrm=10.d0**fipoi(xmlog,m,xq,xr)
  xom=10.d0**fipoi(xmlog,m,xq,xo)

  xjspe2=(xom*xrm*xrm)*Rsol*Rsol/1.d16

  return

end subroutine momspe

!=======================================================================
subroutine omenex
!------------------------------------------------------------------------
! Calcul le nouveau profil de rotation
!------------------------------------------------------------------------
  use strucmod,only: rprov

  implicit none

  integer::jo
  real(kindreal), dimension(ldi):: xmr
!------------------------------------------------------------------------
  xmr=0.0d0

  do jo=1,m
   xmr(jo)=exp(glm)*(1.d0-exp(q(jo)))
  enddo

! calcul du moment d'inertie de chaque coquille du modele courant
! vxmoci is the angular momentum of the current model computed in om2old
  do jo=2,m-1
   xinert(jo)=2.d0/3.d0*exp(2.d0*rprov(jo))*(xmr(jo-1)-xmr(jo+1))/2.d0
  enddo
  xinert(1)=2.d0/3.d0*exp(2.d0*rprov(1))*(xmr(1)-xmr(2))/2.d0+vsuminenv
  xinert(m)=2.d0/5.d0*exp(2.d0*rprov(m-1))/2.d0**(2.d0/3.d0)*xmr(m-1)/2.d0

! calcul du nouvel omega
  do jo=1,m
   omegi(jo)=vxmoci(jo)/xinert(jo)
  enddo

  return

end subroutine omenex

!=======================================================================
subroutine om2old
!------------------------------------------------------------------------
! Computes the angular momentum of the current model,
! with r coming from the last iteration.
! vomegi has been identified to omegi just before the call.
!------------------------------------------------------------------------
  use rotmod,only: vomegi

  implicit none

  integer:: jo
  real(kindreal), dimension(ldi):: xmr
!------------------------------------------------------------------------
  xmr(:)=0.0d0

  do jo=1,m
   xmr(jo)=exp(glm)*(1.d0-exp(q(jo)))
  enddo

! vxmoci contient le moment cinetique du modele courant
  do jo=2,m-1
   vxmoci(jo)=2.d0/3.d0*exp(2.d0*r(jo))*(xmr(jo-1)-xmr(jo+1))/2.d0*vomegi(jo)
  enddo
  vxmoci(1)=2.d0/3.d0*exp(2.d0*r(1))*(xmr(1)-xmr(2))/2.d0*vomegi(1)+vsuminenv*vomegi(1)
  vxmoci(m)=2.d0/5.d0*exp(2.d0*r(m-1))/(2.d0**(2.d0/3.d0))*xmr(m-1)/2.d0*vomegi(m)

  return

end subroutine om2old

!=======================================================================
subroutine omesta
! Calcul le nouveau profil de rotation
!------------------------------------------------------------------------
  use const,only: Rsol
  implicit none

  integer::jo
  real(kindreal), dimension(ldi):: xmr
!------------------------------------------------------------------------
  xmr=0.0d0
  do jo=1,m
   xmr(jo)=exp(glm)*(1.d0-exp(q(jo)))
  enddo

! calcul du moment d'inertie de chaque coquille du modele courant
! vxmoci contient le moment cinetique du modele precedent
  do jo=2,m-1
   xinert(jo)=2.d0/3.d0*exp(2.0d0*r(jo))*(xmr(jo-1)-xmr(jo+1))/2.d0
  enddo
  xinert(1)=2.d0/3.d0*exp(2.0d0*r(1))*(xmr(1)-xmr(2))/2.d0 + vsuminenv
  xinert(m)=2.d0/5.d0*exp(2.0d0*r(m-1))/2.d0**(2.d0/3.d0)*xmr(m-1)/2.d0

! calcul du nouvel omega
  do jo=1,m
   omegi(jo)=vxmoci(jo)/xinert(jo)
  enddo

  return

end subroutine omesta

!***********************************************************************
end module omegamod
