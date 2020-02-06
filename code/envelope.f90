module envelope
!-----------------------------------------------------------------------
! The entry points in the envelope subroutines are DRECK and DRECKF
! DRECK computes the values at the tips of the triangle
! DRECKF calls GGW
! GGW calls ATMOS, ANFITG, and DIFF3
!   ANFITG calls RSGL or RSGL1, and PRINT1
!   DIFF3 calls RSGL or RSGL1
!-----------------------------------------------------------------------

use evol,only: kindreal
use const,only: um,Msol,Lsol,lgLsol,lgqapicg
use inputparam,only: irot,omega,my,elph,verbose
use caramodele,only: nwmd,gls,glm,teff
use strucmod,only: id1,id2,it2,ih,ihv,f,g,h,u,rtp,rtt,rtc,kk,dk,drl,drte,drp,drt,drr,rlp,rlt,rlc,rrp,rrt,rrc,neudr,Eddmax,profIon,&
  fitmIon,vlm,vll,vlr,vlt,vlte,vlro,vlmm,vlrr,vlll,vlgr,rap1_atm,xft_atm,rap2_atm,xgmoym_atm,xpsi_atm,chem,ychem,vmol,vny,izsa, &
  nr,x_env,konv,konvv,vlro,vlt,vlhp,vmion,vmionp,vmiont,vna,vnr,vlka,vkap,vkat,y4,y5int,uvlpm,uvlp,vltm,vlrm,nenvel,envel, &
  vlro_thermo,rhp_thermo,rht_thermo,xmol_thermo,pradlim,xgmoym_atm,tatmos,patmos,vlttau,taulim,beta_env,dwdo2,cp,vlmix,kzone

implicit none

logical:: notFullyIonised,supraEdd

private
public :: dreckf,dreck
public :: notFullyIonised,supraEdd

contains
!======================================================================
subroutine dreckf
!-----------------------------------------------------------------------
! Derniere version : 1 septembre 1992
  use strucmod,only: q,p,t

  implicit none

  integer:: i
  real(kindreal):: vmkrit,nggw,vlnl,vlnte,pggw,tggw,rggw,det,x5dr,x6,x7,x8,x9
!-----------------------------------------------------------------------
  vmkrit= -log(1.d0-exp(q(1)))/um
!           : log fitm : log M/Mr

  if (id2 /= 5) then

! Dans le cas (id2=5), il faut recalculer certains sommets du triangle.
! La nouvelle integration sera a tester pour s'assurer de la correction
! du nouveau modele. L'atmosphere ainsi obtenue dans ggw ne sera pas
! imprimee.
! Dans le cas (id2=6, en entree), l'appel de ggw se fait pour des valeurs
! testees. L'atmosphere obtenue lors de cet appel de ggw est celle du
! nouveau modele, et peut etre imprimee.
    if (id2  ==  1) then
       nggw = 5
    else
       nggw = 6
    endif
    id2=5
    vlnl= log(gls) + log(Lsol)
!     : Ln L*
    vlnte= rtp*p(1) + rtt*t(1) + rtc
!     : Ln Teff

    call ggw(glm,vlnl,vlnte,vmkrit,6,pggw,tggw,rggw)

  else   ! id2 = 5
! Calcul de (Pf,Tf,rf)(i) pour les sommets (i) modifies du triangle
    write(3,*) 'atmos: P       tau       T_tau       Teff       mu'
    do i=1,3
     if (kk(i) /= 1) cycle
     vlnl= drl(i)
     vlnte=drte(i)
!      write(3,'(a,i2,2(1x,d24.18))')'i,drl(i),drte(i):',i,drl(i),drte(i)

! Calcul de P=Ln Pf(i)=drp(i), drt(i) et drr(i)
     call ggw(glm,vlnl,vlnte,vmkrit,5,drp(i),drt(i),drr(i))
    enddo

! Calcul de (alpha(i),beta(i),gama(i)) (i=1,2,3)
! intervenant dans les equations (27,28,29)
    det=(drp(1)-drp(2))*(drt(1)-drt(3))-(drp(1)-drp(3))*(drt(1)-drt(2))
! alpha2:
    rlp=((drl(1)-drl(2))*(drt(1)-drt(3))-(drl(1)-drl(3))*(drt (1)-drt(2)))/det
! beta2:
    rlt=-((drl(1)-drl(2))*(drp(1)-drp(3))-(drl(1)-drl(3))*(drp(1)-drp(2)))/det
! gama2:
    rlc=drl(1)-rlp*drp(1)-rlt*drt(1)
! alpha1:
    rrp=((drr(1)-drr(2))*(drt (1)-drt (3))-(drr(1)-drr(3))*(drt(1)-drt(2)))/det
! beta1:
    rrt=-((drr(1)-drr(2))*(drp(1)-drp(3))-(drr(1)-drr(3))*(drp(1)-drp(2)))/det
! gama1:
    rrc=drr(1)-rrp*drp(1)-rrt*drt(1)
! alpha3:
    rtp=((drte(2)-drte(1))*(drt(3)-drt(1))-(drte(3)-drte(1))*(drt(2)-drt(1)))/det
! beta3:
    rtt=((drp(2)-drp(1))*(drte(3)-drte(1))-(drte(2)-drte(1))*(drp(3)-drp(1)))/det
! gama3:
    rtc=drte(1)-rtp*drp(1)-rtt*drt(1)
    write(3,*)'Triangle:      log l  log te    log p   log t   log r'

! Passage en log(base 10) pour l'impression
    do i=1,3
     x5dr=(drl(i)-log(Lsol))/um
     x6=drte(i)/um
     x7=drp(i)/um
     x8=drt(i)/um
     x9=drr(i)/um
     write(3,'(12x,2f8.4,f9.4,f8.4,f9.4)') x5dr,x6,x7,x8,x9
    enddo

  endif

  return

end subroutine dreckf
!======================================================================
subroutine dreck(nndr_in)
!-----------------------------------------------------------------------
! using nndr_in here since main calls dreck with
! either 0 or nndr as input
!-----------------------------------------------------------------------
  use inputparam,only: deltat,deltal

  implicit none

  integer,intent(in):: nndr_in
  real(kindreal):: deltat2,deltal2,drel,dret,x3,x4
!-----------------------------------------------------------------------
  if (id1 /= 1) then
    if (id1  ==  2) then
! change of both triangle location and size
      drl(1) = log(Lsol)*1.27d0
      drte(1) = log(300.d0)*1.27d0
      deltat = deltat**(1.05d0-real(mod(nwmd,2))/10.d0)
      deltal = deltal**(1.05d0-real(mod(nwmd,2))/10.d0)
      write(997,'(i7.7,2(a,f8.5))')nwmd,': DELTAL=',deltal,' DELTAT=',deltat
      write(*,'(2(a,f8.5))')'    DELTAL=',deltal,' DELTAT=',deltat
      id1 = 0
    else
      drl(1)=log(Lsol)
      drte(1)=log(300.d0)
    endif
    if (nndr_in  ==  2) then
      deltat2 = 1.2d0*deltat
      deltal2 = 1.2d0*deltal
    else
      deltat2 = deltat
      deltal2 = deltal
    endif

    drl(2)=drl(1)
    drte(2)=drte(1)-deltat2*um
    drl(3)=drl(1)+deltal2*um
    drte(3)=drte(1)
    dk=-1.d0
  endif

  drel=log(gls)+log(Lsol)
  dret=log(teff)
  neudr=0

  kk=0

  do
   if (dk*((drl(2)-drl(3))*(dret-drte(2))+(drte(3)-drte(2))*(drel-drl(2))) < 0.d0) then
     drl(1)=drl(3)+drl(2)-drl(1)
     drte(1)=drte(3)+drte(2)-drte(1)
     dk=-dk
     kk(1)=1
     cycle
   endif
   if (dk*((drl(3)-drl(1))*(dret-drte(3))+(drte(1)-drte(3))*(drel-drl(3))) < 0.d0) then
     drl(2)=2.d0*drl(1)-drl(2)
     drte(2)=2.d0*drte(1)-drte(2)
     dk=-dk
     kk(2)=1
     cycle
   endif
   if (dk*((drl(1)-drl(2))*(dret-drte(1))+(drte(2)-drte(1))*(drel-drl(1))) < 0.d0) then
     drl(3)=2.d0*drl(1)-drl(3)
     drte(3)=2.d0*drte(1)-drte(3)
     dk=-dk
     kk(3)=1
     cycle
   else
     exit
   endif
  enddo
  if ((kk(1)+kk(2)+kk(3)) /= 0) then
    neudr=1
    x3=drel/um-lgLsol
    x4=dret/um
    write(3,'(////,a,3x,a,f8.4,3x,a,f8.4)')'Estimated values','log l=',x3,'log te=',x4
  endif

  return

end subroutine dreck
!======================================================================
subroutine ggw(vlnm,vlnl,vlnte,vmkrit,it,p,t,r)
!-----------------------------------------------------------------------
! Calcule, par interpolation lineaire :
!          P = Ln Pf(i) = DRP(i),
!          T = Ln Tf(i) = DRT(i),
!          R = Ln rf(i) = DRR(i), pour le sommet i du triangle.

! GGW est appellee par DRECKF,
!      et appelle ATMOS, ANFITG, DIFF3, et KONVEK.

! Derniere version : 8 septembre 1993
!---------------------------------------------------------------------
! Entrees : vlnm   : Ln M
!           vlnl   : Ln L
!           vlnte  : Ln Teff
!           vmkrit : log(M/Mr) = -log(FITM)
!           it     : nsix (impression de l'atmosphere) ou ncinq

! Sorties : p : Ln Pf
!           t : Ln Tf
!           r : Ln rf
!----------------------------------------------------------------------
  use const,only: cstlg_sigma,lgpi,cstlg_G,cst_G
  use rotmod,only: suminenv
  use geomod, only: sund_max,geomat
  use ionisation,only: iatoms,abond,a_ion,list,vnu,iz,vmyion,ionized
  use convection,only: konvek

  implicit none

  integer,intent(in):: it
  real(kindreal),intent(in):: vlnm,vlnl,vlnte,vmkrit
  real(kindreal),intent(out):: p,t,r

  integer:: n,i
  real(kindreal),parameter:: vmkrit_tol=0.001d0
  real(kindreal):: vmms,vlls,vpsi,vsum,e,vlmg,xllEdd,ff,FITM
!----------------------------------------------------------------------
! [Modif CG]
! Initialisation de la variable contenant la premiere couche completement ionisee de l'enveloppe.
  ProfIon = 1
  fitmIon = 1.d0
! [\Modif]

  Eddmax = 5.d0

  it2 = it
! Si it = nsix : Impression de l'atmosphere
  vmms= exp(vlnm-log(Msol)) ! M/Msoleil
  vlls=exp(vlnl-log(Lsol))  ! L/Lsoleil
  vlm=vlnm/um               ! log M
  vll=vlnl/um               ! log L
  vlte=vlnte/um             ! log Teff

  if (irot == 0 .or. omega<=0.d0) then
! r (L=4.pi.sigma.r.r.Teff^4):
    vlr= 0.5d0*(vll-4.d0*vlte - cstlg_sigma - log10(4.d0) - lgpi)
    vlmm=vlm
    vlrr=vlr
    vlll=vll
    g=10.d0**(cstlg_G + vlm - 2.d0*vlr)
  else
    vpsi=vll-4.d0*vlte-cstlg_sigma
    vpsi=10.d0**(vpsi-2.d0/3.d0*(cstlg_G+vlm-2.d0*log10(omega)))
    if (vpsi < 0.04298317498d0) then
! r (L=4.pi.sigma.r.r.Teff^4):
      vlr= 0.5d0*(vll-4.d0*vlte - cstlg_sigma - log10(4.d0) - lgpi)
      vlmm=vlm
      vlrr=vlr
      vlll=vll
      g=1.d0/10.d0**(cstlg_G + vlm - 2.d0*vlr)
      rap1_atm=1.0d0
      xft_atm=1.0d0
      rap2_atm=1.0d0
      xgmoym_atm=g
    else
      if (vpsi > sund_max) then
        vpsi=0.999999d0* sund_max
      endif
      call geomat(vpsi,xpsi_atm,rap1_atm,xft_atm,rap2_atm,xgmoym_atm)
      if (xpsi_atm <= 0.d0) then
        stop ' (xpsi <= 0) dans ggwr.f'
      endif
      vlr=1.d0/3.d0*(cstlg_G+vlm-2.d0*log10(omega))+log10(xpsi_atm)
      xgmoym_atm=xgmoym_atm/(cst_G*10.d0**(vlm)*omega**4.d0)**(1.d0/3.d0)
      vlmm=vlm
      vlrr=vlr
      vlll=vll
      g=xgmoym_atm
    endif
  endif

  if (it2 == 6) write(3,'(3x,a//10x,a,f10.4,2x,a,f13.4,2x,a,e11.3/10x,a,f9.4,2x,a,f10.4,2x,a,f11.3,2x,a,f9.3//)') &
                 'Input data for outer layers integration','mass= ',vmms,'luminosity= ',vlls,'chem= ',chem,'vlm=',vlm, &
                 'vll=',vll,'vlte=',vlte,'vlr=',vlr
!---------------------------------------------------------------------
!  MODIFICATIONS MAI 1990, D.SCHAERER:
!  IONISATION PARTIELLE POUR PLUSIEURS ELEMENTS...
! vmol : poids moleculaire moyen de la matiere non ionisee.
! vny(i) : fraction du nombre d'atomes de l'element i
!          dans le melange chimique.
  vmol=0.d0
  vsum=0.d0
  do n=1,iatoms
   vmol=vmol+abond(n)/a_ion(n)
   vsum=vsum+abond(n)
  enddo
  vmol=vmol+0.5d0*(1.d0-vsum)
  vmol=1.d0/vmol
!  Mettre a jour anciennes variables VNY() (NE PAS UTILISER):
  vny(1)=vmol*abond(1)
  vny(2)=vmol*abond(2)/a_ion(2)
  vny(3)=vny(2)
  e=0.d0
  do n=1,iatoms
   if (list(n) == 0) exit
   vnu(list(n))=vmol*abond(list(n))/a_ion(list(n))
   e=e+iz(list(n))*vnu(list(n))
  enddo
  vmyion=vmol/(1.d0+e)
  ionized=0
!-----------------------FIN PREMIERE PARTIE...-------------------------
  izsa=0
  vlmg=vlm
  call atmos
  call anfitg
! anfitg initialise ne a zero
!-----------------------------------------------------------------------
!     MODIFICATIONS D.SCHAERER:
!     Pour diminuer le temps de calcul passe dans IONPART:
!     Des que dans l'integration de l'enveloppe le poids moleculaire (VM
!     atteint la valeur qui correspond a l'ionisation totale (VMYION)
!     ceci sera signale a IONPART en mettant IONIZED=1.

! Integration de l'enveloppe (ou "sub-atmosphere")

! Resolution des equations (21,22,23) entre tautilde et fitm.
! Boucle sur les couches de l'enveloppe.
  FITM = 10.d0**(-vmkrit)
!  do while ((vlmg-vlm-vmkrit) <= 0.d0)
  do while ((vlmg-vlm+log10(FITM+vmkrit_tol*(1.d0-FITM))) <= 0.d0)
   call diff3(vlmg,vmkrit)
   if (nr > 500) then
      rewind(222)
      write (222,*) nwmd,': nr greater than 500 in GGW'
     stop 'NR GREATER THAN 500 IN GGW.'
   endif

! Situation par rapport a une zone convective
   call konvek
! konvek incremente nr, numero de couche

   if (x_env(1) < 0.9990d0 .or. x_env(3) < 0.9990d0) then
     notFullyIonised = .true.
   endif

   if (x_env(3)  <  1.0d0 .or. x_env(1)  <  1.0d0) then
     ProfIon = nr
     fitmIon = 10.d0**(vlm-vlmm)
   endif

   if (it2 == 6) call print1

   nenvel=nr
   envel(nr,1)=real(nr)
   envel(nr,2)=real(konv)
   envel(nr,3)=vlr
   envel(nr,4)=vlro
   envel(nr,5)=vlm
   envel(nr,6)=vlt
   envel(nr,7)=vmion
   envel(nr,8)=vna
   envel(nr,9)=vnr
   envel(nr,11)=vll
   envel(nr,15)=x_env(2)
   envel(nr,16)=x_env(3)
   envel(nr,17)=vlka
   envel(nr,18)=vkap
   envel(nr,19)=vkat

! L / L_Eddington
   xllEdd = 10.d0**(vll-vlm+vlka-lgqapicg)
     if (xllEdd  >  Eddmax) then
       supraEdd = .true.
     endif

   if (vmion-vmyion < 1.d-3) ionized=1

  enddo
! Dans ce cas, on n'a pas atteint fitm. Il faut continuer.

  ionized=0
!----------------------FIN DES MODIFICATIONS----------------------------
! Interpolation lineaire pour recalculer log Pf, log Tf, log rf
  ff=(vlmg-vmkrit-y4(3))/(vlm-y4(3))
  uvlpm=uvlp-h+ff*h
  vltm=y4(1)+ff*(vlt-y4(1))
  vlrm=y4(2)+ff*(vlr-y4(2))
  if (it2 == 6) write(3,'(/////,3(a,f9.4),///)') '  uvlpm=',uvlpm,'  vltm=',vltm,'  vlrm=',vlrm
  p=uvlpm*um   ! Ln Pf
  t=vltm*um    ! Ln Tf
  r=vlrm*um    ! Ln rf

  suminenv = 2.d0/3.d0*10.d0**(2.0d0*envel(1,3))*(10.d0**(envel(1,5))-10.d0**(envel(2,5)))/2.d0
  do i=2,nr-1
! calcul du moment d'inertie de chaque coquille
   suminenv = suminenv + 2.d0/3.d0*10.d0**(2.0d0*envel(i,3))*(10.d0**(envel(i-1,5))-10.d0**(envel(i+1,5)))/2.d0
  enddo
  suminenv = suminenv + 2.d0/3.d0*10.d0**(2.0d0*vlrm)*(10.d0**(envel(nr-1,5))-exp(vlnm)*10.d0**(-vmkrit))/2.d0
  if (isnan(suminenv)) then
    write(*,*) 'suminenv=NaN - nr,vlrm,envel(nr-1,5),vlnm,vmkrit:',nr,vlrm,envel(nr-1,5),vlnm,vmkrit
    stop
  endif
  return

end subroutine ggw
!======================================================================
subroutine atmos
!-----------------------------------------------------------------------
! Derniere version : 28 juillet 1993
!-----------------------------------------------------------------------
  use const,only: cstlg_a,rgazlg,um
  use opacity,only: kappa
  use ionisation,only: ionpart
  use equadiffmod,only: iprc
  use PrintAll, only:StoreStructure_atm

  implicit none

  integer::l,lq,j_kap
  integer::Loop26,Loop22,Atmos_Layer

  real(kindreal),parameter:: h1m=8.0d-2
  real(kindreal):: arg1,arg2,hf,rhoatmos,taum,xcmp
  real(kindreal), dimension(2), save:: ed1=(/1.0d-4,1.0d-2/)
  real(kindreal), dimension(3):: hp,taust,vltta
!-----------------------------------------------------------------------
  taum=2.d0/3.d0

  Atmos_Layer = 1

  if (irot == 1) then
    taum=(4.d0/3.d0-(2.d0/3.d0)*rap2_atm)/(rap1_atm*xft_atm)
  endif

  if (it2 == 6) then
    write(3,'(a,///,a,/,a)') ' Atmosphere','integration step','    uvlp    tau     vlt    beta     vlro      vlka'
  endif
  h=2.0d-03
  j_kap = 1

  taust(1)=0.d0
! vltta(1): T0 = Teff/(2**0.25) = Teff a tau=0 (Eq. 17 rapport Patenaude)
  vltta(1)=vlte+(1.d0/4.d0)*log10(1.d0/2.d0)

  if (irot == 1) then
   vltta(1)=vlte+(1.d0/4.d0)*log10(1.d0/2.d0)+0.25d0*log10(rap2_atm)
  endif

! hp(1): a/3 * T0**4
  hp(1)= cstlg_a - log10(3.d0) + 4.d0*vltta(1)
!--------------------------------------------------------------------
  do l=2,4
   xcmp=chem
   call ionpart(hp(l-1),vltta(l-1))

   vlro=log10(vmion*beta_env)+hp(l-1)-vltta(l-1)-rgazlg
   rhp_thermo= 1.d0/beta_env + vmionp
   rht_thermo = -4.d0/beta_env + 3.d0 + vmiont

   arg1=vlro*um
   arg2=vltta(l-1)*um
   call kappa(arg1,arg2,rhp_thermo,rht_thermo,chem,ychem,vlka,vkap,vkat,j_kap)

   vlka=vlka/um
! f(3-5,1) = kappa * P / g = d tau / d ln P (Eq. 19 rapport Patenaude)
   f(l+1,1)=10.d0**(vlka+hp(l-1))*um/g
   if (irot == 1) then
     f(l+1,1)=10.d0**(vlka+hp(l-1))*um*rap1_atm*xgmoym_atm
   endif
   if (l == 4) cycle
   taust(l)=taust(l-1) + h*f(l+1,1)
! Calcul de la temperature a la profondeur optique tau:
! vlte: Teff, taust: tau
   vltta(l)=vlte+0.25d0*log10(0.75d0*(taust(l)+(2.d0/3.d0)))
   if (irot == 1) then
     vltta(l)=vlte+0.25d0*log10(0.75d0*(taust(l)*rap1_atm*xft_atm+(2.d0/3.d0)*rap2_atm))
   endif
   hp(l)=hp(l-1) + h
  enddo
!--------------------------------------------------------------------
  if (it2 == 6) then
    do lq=1,3
     write(3,'(1x,f7.3,f9.4,f8.3,f8.3,f10.3,f9.3,f10.4)') hp(lq),taust(lq),vltta(lq),beta_env,vlro,vlka
    enddo
    lq=3
    if (iprc == 1) then
      if (hp(lq) < 1.e-30) then
        hp(lq) = 0.0d0
      endif
      call StoreStructure_atm(Atmos_Layer,vltta(lq),vlro,hp(lq),rhp_thermo,-rht_thermo,vna,cp,vlka,vkap,vkat, &
           taust(lq),vmion,vmol,x_env(:))
    endif
    Atmos_Layer = Atmos_Layer + 1
! Calcul de la valeur  Prad en Tau=0
    pradlim = (1.d0 - beta_env) * 10.d0**hp(1)
  endif
  f(2,1)=0.d0
  y5int(1)=taust(3)
  uvlp=hp(3)
  ih=1
  Loop22 = 0
22 select case (ih)
  case (1,2)
    uvlp = uvlp + h
    y4(1) = y5int(1)
    do l = 1,4
     f(l,1) = f(l+1,1)
    enddo
  case (3)
    h=2.d0*h
    uvlp = uvlp + h
    y4(1) = y5int(1)
    f(2,1) = f(1,1)
    f(4,1) = f(5,1)
  case default
    stop 'problem with ih in atmos'
  end select
  Loop26 = 0
! Methode Adams Bashforth predictor (cf. numerical recipe)
! f(2,1): derivee en n-2; f(3,1): derivee en n-1; f(4.1): derivee en n
26 vlt=y4(1)+h*(5.d0*f(2,1)-16.d0*f(3,1)+23.d0*f(4,1))/12.d0
! Calcul de la temperature a la profondeur optique tau:
! vlte: Teff, vlt: tau
  vlttau=vlte+0.25d0*log10(0.75d0*(vlt+2.d0/3.d0))

! Correction due a la rotation, cf. Meynet & Maeder 1997, annexe 3
  if (irot == 1) then
    vlttau=vlte+0.25d0*log10(0.75d0*(vlt*rap1_atm*xft_atm+(2.d0/3.d0)*rap2_atm))
  endif

  xcmp=chem
!---------------- Modifications de Schaerer 1990 -----------------------
  call ionpart(uvlp,vlttau)
  vlro=log10(vmion*beta_env) + uvlp - vlttau - rgazlg
  rhp_thermo=1.d0/beta_env + vmionp
  rht_thermo=-4.d0/beta_env+3.d0+vmiont
!-----------------------Fin des modifications---------------------------
  arg1=vlro*um
  arg2=vlttau*um

  call kappa(arg1,arg2,rhp_thermo,rht_thermo,chem,ychem,vlka,vkap,vkat,j_kap)

  vlka=vlka/um
  f(5,1)=10.d0**(vlka+uvlp)*um/g

  if (irot == 1) then
    f(5,1)=10.d0**(vlka+uvlp)*um*rap1_atm*xgmoym_atm
  endif

  y5int(1)=y4(1)+h*((1.d0/24.d0)*f(2,1)-(5.d0/24.d0)*f(3,1) + (19.d0/24.d0)*f(4,1)+(3.d0/8.d0)*f(5,1))
  if (abs(y5int(1)-vlt) > ed1(2)) then
    ih = 1
    h = 0.5d0 * h
    uvlp = uvlp - h
    hf = f(2,1)
    f(2,1) = f(3,1)
    f(3,1) = (3.d0/8.d0)*f(4,1)+(3.d0/4.d0)*f(3,1)-(1.d0/8.d0)*hf
    Loop26 = Loop26 + 1
    if (Loop26 >= 500) then
      rewind(222)
      write (222,*) nwmd,': Loop 26 problem in atmos'
      stop 'Loop 26 problem in atmos'
    endif
    go to 26
  endif
  if (abs(y5int(1)-vlt)>=ed1(1) .or. ih/=2 .or. 2.d0*h>h1m) then
    ih = 2
  else
    ih = 3
  endif
  if (it2 == 6) then
    write(3,'(1x,f7.3,f9.4,f8.3,f8.3,f10.3,f9.3,f10.4)') uvlp,y5int(1),vlttau,beta_env,vlro,vlka,xmol_thermo
    if (iprc == 1) then
      call StoreStructure_atm(Atmos_Layer,vlttau,vlro,uvlp,rhp_thermo,-rht_thermo,vna,cp,vlka,vkap, &
           vkat,y5int(1),vmion,vmol,x_env(:))
    endif
    Atmos_Layer = Atmos_Layer + 1
  endif
!-----------------------------------------------------------------------
  if (y5int(1) < taum) then
    Loop22 = Loop22 + 1
    if (Loop22 >= 500) then
      rewind(222)
      write (222,*) nwmd,': Loop 22 problem in atmos'
      stop 'Loop 22 problem in atmos'
    endif
    go to 22
  endif
  uvlp=uvlp-h+h*(taum-y4(1))/(y5int(1)-y4(1))
! Calcul de la temperature a la profondeur optique tau:
! vlte: Teff, taum: tau
  vlttau=vlte+0.25d0*log10(0.75d0*(taum+(2.d0/3.d0)))

  if (irot == 1) then
    vlttau=vlte+0.25d0*log10(0.75d0*(taum*rap1_atm*xft_atm+(2.d0/3.d0)*rap2_atm))
  endif
  if (it2 == 6) then
    write(3,'(1x,f7.3,f9.4,f8.3,f8.3,f10.3,f9.3,f10.4)') uvlp,taum,vlttau,vlte,vlte,vlte,xmol_thermo
    if (iprc == 1) then
      call StoreStructure_atm(Atmos_Layer,vlttau,vlro,uvlp,rhp_thermo,-rht_thermo,vna,cp,vlka,vkap, &
           vkat,y5int(1),vmion,vmol,x_env(:))
    endif
    Atmos_Layer = Atmos_Layer + 1
  endif
!...................POUR CORRWIND...........................
  patmos=uvlp
  tatmos=vlttau
  rhoatmos=log10(vmion*beta_env) + uvlp - vlttau - rgazlg
!.......................................................................
  if (it2 == 6) write(3,'(a)')'End of the atmosphere'
!---------------------------------------------------------------------------------
  return

end subroutine atmos
!======================================================================
subroutine anfitg
!----------------------------------------------------------------------
! Integration de l'enveloppe au voisinage de (tau tilde).
!  Lr=L*.

! Integration des equations (21), (22), (23) et (24).

! Convection traitee de maniere non-adiabatique.

! ANFITG est appelee par GGW
!        et appelle RSGL ou RSGL1, et KONVEK.

! Derniere version : 28 septembre 1992
!----------------------------------------------------------------------
  use convection,only: konvek

  implicit none

  integer,parameter:: iterm=10
  real(kindreal),dimension(2),save:: ed3=(/1.0d-4,1.0d-2/)

  integer:: iter,j,n
  real(kindreal):: hvlt3,huvlp,vlt1,vlt3
!----------------------------------------------------------------------
! log P(tautilde) : uvlp
  huvlp=uvlp
  konvv=0
  kzone=0
  ihv=2
  j=1

12 vlt1 = vlttau
! log T(tautilde):
  vlt3 = vlttau
  iter = 0

! Preparation des coefficients de l'interpolation de Adams.
13 do n = 1,3
   select case (n)
   case(1)
     uvlp = huvlp-h
     vlt = vlt1
   case(2)
     uvlp = huvlp
     vlt = vlttau
   case(3)
     uvlp = huvlp+h
     vlt = vlt3
   case default
     write(*,*) ' problem in anfitg...'
     rewind(222)
     write (222,*) nwmd,': problem in anfitg ==> STOP'
     stop
   end select

   if (my >= 1) then
     call rsgl1
   else
     call rsgl
   endif

   f(n+2,1) = f(5,1)
   f(n+2,2) = 0.d0
   f(n+2,3) = 0.d0
   f(5,2) = 0.d0
   f(5,3) = 0.d0
   if (j == 1) then
     cycle
   endif
   call konvek

! Impression des parametres de l'enveloppe
   if (it2 == 6) call print1

   nenvel = nr
   envel(nr,1) = real(nr)
   envel(nr,2) = real(konv)
   envel(nr,3) = vlr
   envel(nr,4) = vlro
   envel(nr,5) = vlm
   envel(nr,6) = vlt
   envel(nr,7) = vmion
   envel(nr,8) = vna
   envel(nr,9) = vnr
   envel(nr,11) = vll
   envel(nr,17) = vlka
   envel(nr,18) = vkap
   envel(nr,19) = vkat

   if (n == 3) then
     ih=1
     if (it2 == 6) write(3,'(a,3x,a,i2,2x,a,e12.5)') 'End of first iteration','iter=',iter,'h=',h
     return
   endif
  enddo
  continue

! Interpolation de Adams
! y(n-1):
  vlt1=vlttau +h*(-(5.d0/12.d0)*f(3,1)-(2.d0/3.d0)*f(4,1)+(1.d0/12.d0)*f(5,1))
! y(n+1):
  hvlt3=vlttau +h*(-(1.d0/12.d0)*f(3,1)+(2.d0/3.d0)*f(4,1)+(5.d0/12.d0)*f(5,1))
  iter = iter+1
! check this:
!      print*,'if the code runs through here, please check
!     & the next statements: IF is not necessary'
  if ( abs(vlt3-hvlt3)- sqrt(ed3(1)*ed3(2))) 22,22,23
22 j = 2
  nr = 0
  go to 24
23 if (iter-iterm) 24,24,25
24 vlt3 = hvlt3
  go to 13
25 h=0.5d0*h
  ihv = 1
  go to 12

  return

end subroutine anfitg
!======================================================================
subroutine diff3(vlmg,vmkrit)
!----------------------------------------------------------------------
! Derniere version : 25 septembre 1992
! vmkrit = -log(FITM)
! vlmg = log(M_total) total mass of the star
! vlm = log(M_r) mass at the current level
!----------------------------------------------------------------------
  implicit none

  real(kindreal),parameter:: h3m=1.0d-1

  integer:: k,l,n
  real(kindreal), intent(in):: vmkrit,vlmg
  real(kindreal):: hf
  real(kindreal),dimension(2),save:: ed3=(/1.0d-4,1.0d-2/)
!----------------------------------------------------------------------
  if (ih <= 2) then
    uvlp=uvlp+h
    y4(1) = vlt
    y4(2) = vlr
    y4(3) = vlm
    do  k=1,3
     do  l=1,4
      f(l,k)=f(l+1,k)
     enddo
    enddo
  else
    h=2.d0*h
    uvlp = uvlp + h
    y4(1) = vlt
    y4(2) = vlr
    y4(3) = vlm
    do  k = 1,3
      f(2,k) = f(1,k)
      f(4,k) = f(5,k)
    enddo
  endif

  do
    vlt=y4(1)+h*((5.d0/12.d0)*f(2,1)-(4.d0/3.d0)*f(3,1)+(23.d0/12.d0)*f(4,1))
    vlr=y4(2)+h*((5.d0/12.d0)*f(2,2)-(4.d0/3.d0)*f(3,2)+(23.d0/12.d0)*f(4,2))
    vlm=y4(3)+h*((5.d0/12.d0)*f(2,3)-(4.d0/3.d0)*f(3,3)+(23.d0/12.d0)*f(4,3))
    if (my < 1) then
      call rsgl
    else
      call rsgl1
    endif
    do n=1,3
     y5int(n)=y4(n)+h*((1.d0/24.d0)*f(2,n)-(5.d0/24.d0)*f(3,n)+(19.d0/24.d0)*f(4,n)+(9.d0/24.d0)*f(5,n))
    enddo
    if (y5int(3) - vlmg + vmkrit > 0.d0) then
      if(abs(y5int(1)-vlt)-ed3(2)<= 0.d0) exit
    endif
    ih = 1
    h = 0.5d0 * h
    uvlp = uvlp - h
    do k = 1,3
     hf=f(2,k)
     f(2,k) = f(3,k)
     f(3,k) = 0.375d0 * f(4,k)+0.75d0*f(3,k)-0.125d0*hf
    enddo
  enddo

  ihv = ih
  if (abs(y5int(1)-vlt)-ed3(1) < 0.d0) then
    if (ih == 2) then
      if (2.d0*h-h3m <= 0.d0) then
        ih = 3
      else
        ih = 2
      endif
    else
      ih = 2
    endif
  else
    ih = 2
  endif
  vlt = y5int(1)
  vlr = y5int(2)
  vlm = y5int(3)
  if (my < 1) then
    call rsgl
  else
    call rsgl1
  endif

  return

end subroutine diff3
!======================================================================
subroutine rsgl
!----------------------------------------------------------------------
! Membres de droite des equations de structure interne.

! La convection est traite a la boehm-vitense.

! Derniere version : 28 juillet 1993
!----------------------------------------------------------------------
  use const,only: cst_G,rgazlg,cstlg_a,cstlg_G,cstlg_sigma,lgpi
  use geomod,only: rpsi_min,rpsi_max,geom
  use opacity,only: kappa
  use ionisation,only: ionpart
  use convection,only: fconva

! for ifort compiler, uncomment the next line:
!  use, INTRINSIC:: IEEE_ARITHMETIC, only: isnan => IEEE_IS_NAN

  implicit none

  integer:: iter_number,j_kap
  real(kindreal):: grat,grar,gram,xcmp,xpsi,xfp,xratp,dxfp,dxratp,arg1,arg2,y_rsg,d_rsg
!-----------------------------------------------------------------------
  iter_number = 1
  j_kap = 1

  grat = f(5,1)
  grar = f(5,2)
  gram = f(5,3)

! EOS de Dappen 1980 et 1992
  xcmp=chem

  if (irot == 1) then
    xpsi=((omega*omega)/(cst_G*10.d0**vlm))**(1.d0/3.d0)*10.d0**vlr
    if (xpsi >= rpsi_min) then
      if (xpsi > rpsi_max) xpsi=0.9999999999d0*rpsi_max
      call geom(xpsi,xfp,xratp,dxfp,dxratp)
    else
      xfp=1.d0
      xratp=1.d0
    endif
  else   !   irot=0
    xfp=1.d0
    xratp=1.d0
  endif   !   irot

!---------------- Modifications de Schaerer 1990 -----------------------
  call ionpart(uvlp,vlt)
  vlro = log10(vmion*beta_env)+uvlp-vlt-rgazlg
  rhp_thermo = 1.d0/beta_env+vmionp
  rht_thermo = 3.d0-4.d0/beta_env+vmiont
!---------------- Fin des modifications --------------------------------
  arg1 = vlro*um
  arg2 = vlt*um
  call kappa(arg1,arg2,rhp_thermo,rht_thermo,chem,ychem,vlka,vkap,vkat,j_kap)
  vlka = vlka/um
  vnr  = 10.d0**(log10(3.d0)-log10(4.d0)-lgqapicg-cstlg_a+vlka+vll+uvlp-vlm-4.d0*vlt)*xratp

! CRITERE DE SCHWARZSCHILD
  if (vnr <= vna) then
!  EQUILIBRE RADIATIF
    konv = 0
    grat = vnr
  else
!  EQUILIBRE CONVECTIF
    konv = 1
    vlgr = vlm-2.d0*vlr+cstlg_G
    vlhp = uvlp-vlro-vlgr
    vlmix = log10(elph)+vlhp
! Maeder (2009) Eq. (5.71)
    u = 10.d0**(3.d0*vlt-2.d0*(vlro+vlmix)-vlka+0.5d0*(vlhp-vlgr)+cstlg_sigma+log10(24.d0)+log10(2.d0)/2.d0- &
        rgazlg) * vmol/(cp*sqrt(-rht_thermo))
    y_rsg = 8.d0*(vnr-vna)/(u*u)
    if (y_rsg < 64.d0/3.d0) then
      dwdo2 = y_rsg/16.d0
    else
      dwdo2 = (y_rsg/9.d0)**(1.d0/3.d0)
    end if

    do
     d_rsg = (((9.d0*dwdo2+8.d0)*dwdo2+16.d0)*dwdo2-y_rsg)/((27.d0*dwdo2+16.d0)*dwdo2+16.d0)
     dwdo2 = dwdo2-d_rsg
     if (isnan(dwdo2) .or. isnan(d_rsg)) then
       write(222,*) nwmd,': Nan in rsgl.'
       stop 'Nan in rsgl. '
     endif
     if (abs(d_rsg/dwdo2) <= 1.d-7) then
       exit
     endif
     iter_number = iter_number + 1
     if (iter_number > 10000) then
       write(222,*) nwmd,': convergence problem in rsgl.'
       stop 'Convergence problem in rsgl.'
     endif
    enddo

    grat = vna+u*u*dwdo2*(dwdo2+2.d0)

  endif

! dln(r)/dln(P) = -P / (rho g) *1/fp
  grar = -10.d0**(vlr+uvlp-cstlg_G-vlro-vlm)*1.d0/xfp
! dln(M)/dln(P) = -4pir^4 P / (GM^2) *1/fp
  gram = -10.d0**(log10(4.d0)+lgpi-cstlg_G+4.d0*vlr+uvlp-vlm-vlm)*1.d0/xfp

  f(5,1) = grat
  f(5,2) = grar
  f(5,3) = gram

  return

end subroutine rsgl
!======================================================================
subroutine rsgl1
!-----------------------------------------------------------------------
!     MEMBRE DE DROITE DU SYSTEME D'EQUATION DE STRUCTURE DE L'ENVELOPPE

!     LA CONVECTION EST TRAITEE DE LA MANIERE SUIVANTE:
!       - LA LONGUEUR DE MELANGE EST PROPORTIONNELLE A HRHO :
!                  L = ALPHA * HRHO
!       - LE FLUX ACOUSTIQUE EST TRAITE AVEC LA FORMULE DE LIGHTHILL :
!                  FAC = ETA * V**8 / VSON**5
!       - LA PRESSION TURBULENTE EST PRISE EN CONSIDERATION :
!                  PTOT = PTURB + PGAZ + PRAD
!                  PTURB = ZETA * RHO * VCONV**2

! Derniere version : 28 juillet 1993
!-----------------------------------------------------------------------
  use const, only: cst_G,rgazlg,cstlg_a,cstlg_G,rgaz,cstlg_c,lgpi
  use geomod, only: rpsi_min,rpsi_max,geom
  use opacity,only: kappa
  use ionisation,only: ionpart
  use convection,only: fconva

  implicit none

  integer:: i,j,j_kap
  real(kindreal), parameter:: zeta=1.d0/3.d0,eta=1000.d0
  real(kindreal):: grat,grar,gram,xcmp,xpsi,xfp,xratp,dxfp,dxratp,arg1,arg2,vlvs,cp2,vnro,vlmxa,vlua1,vlua2,ua,xlamb, &
                   ca,cb,cc,ff,fz,xx,d,vlv,vlhp1
!-----------------------------------------------------------------------
  j_kap = 1

  grat = f(5,1)
  grar = f(5,2)
  gram = f(5,3)

! EOS de Dappen 1980 et 1992
  xcmp=chem

  if (irot == 1) then
    xpsi=((omega*omega)/(cst_G*10.d0**vlm))**(1.d0/3.d0)*10.d0**vlr
    if (xpsi >= rpsi_min) then
      if (xpsi > rpsi_max) xpsi=0.9999999999d0*rpsi_max
      call geom(xpsi,xfp,xratp,dxfp,dxratp)
    else
      xfp=1.d0
      xratp=1.d0
    endif
  else
    xfp=1.d0
    xratp=1.d0
  endif

!---------------- Modifications de Schaerer 1990 -----------------------
  call ionpart(uvlp,vlt)
  vlro = log10(vmion*beta_env)+uvlp-vlt-rgazlg
  rhp_thermo  = 1.d0/beta_env+vmionp
  rht_thermo  = 3.d0-4.d0/beta_env+vmiont
!----------------- Fin des modifications -------------------------------

  arg1 = vlro*um
  arg2 = vlt *um
  call kappa(arg1,arg2,rhp_thermo,rht_thermo,chem,ychem,vlka,vkap,vkat,j_kap)
  vlka = vlka/um
  vnr  = 10.d0**(log10(3.d0)-log10(4.d0)-lgqapicg-cstlg_a+vlka+vll+uvlp-vlm-4.d0*vlt)*xratp
  vlgr = vlm-2.d0*vlr+cstlg_G
  vlhp = uvlp-vlro-vlgr

! CRITERE DE K. SCHWARZSCHILD
  if (vnr <= vna) then
! EQUILIBRE RADIATIF
    konv = 0
    grat = vnr
  else
! EQUILIBRE CONVECTIF
    konv = 1
! vlvs: c_sound
    vlvs = 0.5d0*(uvlp-vlro-log10(rhp_thermo+rht_thermo*vna))
    cp2  = rgaz*cp/vmol
    vnro = -rhp_thermo/rht_thermo
    do i = 1,50
! Voir Maeder 2009 eq. 5.71 pp. 98
! vlmxa= l = elph*H_rho
! Gamma1 = 1/(alpha-delta*Nabla_ad) : Maeder (2009) Eq. 7.57, 7.58, 7.66
     vlmxa = log10(elph/(rhp_thermo+rht_thermo*vna))+vlhp
     vlua1 = 3.d0*vlt-2.d0*vlro-vlka+log10(3.d0)+cstlg_a+cstlg_c-log10(cp2)-vlmxa
     vlua2 = log10(-8.d0/rht_thermo)+vlhp-2.d0*vlmxa-vlgr
     ua = 10.d0**(vlua1+0.5d0*vlua2)
     dwdo2  = 1.d0/ua
     xlamb = eta/9.d0*(rht_thermo*rht_thermo/cp2)*10.d0**(3.d0*(vlua1+vlmxa)-5.d0*vlvs-vlt+2.d0*vlgr-vlhp)
     ca = (vnr-vna)*xlamb
     cb = ua*ua/(vnr-vna)
     cc = (vnro-vna)/(ua*ua)
     do
      call fconva(dwdo2,ff,fz,xx,ca,cb,cc)
      if (fz > 0.d0) exit
      dwdo2 = dwdo2+dwdo2
     enddo

     do j = 1,100
      d = ff/fz
      dwdo2 = dwdo2-d
      if (abs(d/dwdo2) < 1.d-6) exit
      call fconva(dwdo2,ff,fz,xx,ca,cb,cc)
      if (j == 100 .and. verbose) then
        print*,' CONVERGENCE LENTE DANS RSGL : D,dwdo2 :',d,dwdo2
      endif
     enddo

     u  = ua/(xx*xx)
     vlmix= log10(xx)+vlmxa
     grat = vna+u*u*dwdo2*(dwdo2+2.d0)
     vlv  = 0.5d0*(vlgr+2.d0*vlmix-vlhp+log10(-rht_thermo*u*u*dwdo2*dwdo2/8.d0))
     vlhp1= vlhp
     vlhp = log10(zeta*10.d0**(vlv+vlv-vlgr)*(-rht_thermo*(vnro-grat))+10.d0**(uvlp-vlro-vlgr))
     if (abs(vlhp-vlhp1) < 1.d-4 ) exit
     if (i == 50 .and. verbose) then
       print *,' PROBLEME DE CONVERGENCE DANS RSGL: VLHP1,VLHP :',vlhp1,vlhp
     endif
    enddo
  endif

! grar = dln(r)/dln(rho)
! gram = dln(m)/dln(rho) = dln(M)/dln(P) * dln(P)/dln(r) * dln(r)/dln(rho)
!      = 4pi r^2 rho Hp / M
  grar = -10.d0**(vlhp-vlr)*1.d0/xfp
  gram = -10.d0**(log10(4.d0)+lgpi+vlr+vlr+vlro-vlm+vlhp)*1.d0/xfp

  f(5,1) = grat
  f(5,2) = grar
  f(5,3) = gram

  return

end subroutine rsgl1
!======================================================================
subroutine print1
!-----------------------------------------------------------------------
! Impression des parametres de l'enveloppe.
! Les log sont des logarithmes decimaux
!     NR   := numero de la couche (n=1 a la surface)
!     UVLP := log(P)
!     UVLPT:= log(P) avec correction pour la turbulence
!     BETA := P_gas/P_tot
!     VLRO := log(rho [g/cm^3])
!     VLR  := log(r)
!     VLM  := log(mr [g])
!     X(1) := fraction de H+
!     X(2) := fraction de He+
!     X(3) =: fraction of He++
!     VLKA := log(kappa), l'opacite dans la couche nr
!     GA1  := Gamma1
!     GA2  := Gamma2
!     GA3  := Gamma3
!     VLT  := log(T)
!     CP   := Cp*mu0/R ou Cp=chaleur specifique a P=cst.[erg/K/g]
!     VMOL := mu0 poids moleculaire moyen
!     RHP  := dln(rho)/dln(P) (elph)
!     RHT  := dln(rho)/dln(T) (-delta)
!     VMION:= poids moleculaire
!     VMIONP:= dln(mu)/dln(P)
!     VMIONT:= dln(mu)/dln(T)
!     VKAP := dln(KAPPA)/dln(P) a T=cst.
!     BETA := BETA = Pg/Ptot
!     GRAM := dln(M)/dln(P)
!     GRAR := dln(r)/dln(P)
!     GRAT := gradient ext (dlnT/dlnP) ext
!     TM   := turnover time
!     LM   := longueur de melange
!     LM/HP:= id. en unite de H_P
!     VM   := vitesse de convection
!     V/VS := v_conv/v_son
!     VNA  := gradient adiabatique (dlnT/dlnP) adia.
!     VNE  := gradient thermique
!     VNR  := gradient radiatif
!     FR   := flux radiat./flux total
!     FC   := flux conv./flux total
!     FA   := 1-f_rad-f_conv
!     VLMIX:= log(l) ou l= longueur de melange
!     VKAT := (dln(kappa)/dln(T)) a P=cte
!     VKAP := (dln(kappa)/dln(P)) a T=cte
!     R/RTOT:= R tot.
!     M/MTOT:= Masse totale
!     VLL  := Luminosite totale
!     L/LEDD:= Luminosite sur L_Eddington
!     U    := variable (5.71) dans Maeder 2009
!     Z    := variable de transformation pour resoudre U
!     IHV  := compteur dans l'integration de l'enveloppe
!--------------------------------------------------------------------
! Derniere version : 29 juillet 1993
!-----------------------------------------------------------------------
  use equadiffmod,only: iprc
  use PrintAll,only: StoreStructure_env

  implicit none

  real(kindreal), parameter:: eta=1000.d0
  real(kindreal):: zeta,grat,grar,gram,r,xm,rhp,rht,ga1,ga2,ga3,d,xl,vm,tm,vlv,uvlpt,vne,vsvs,al,fr,fc,fa,p_turb
  character(16):: Var
!-----------------------------------------------------------------------
  zeta=1.d0/3.d0
  grat = f(5,1)
  grar = f(5,2)
  gram = f(5,3)

  if (nr-1 == 14*((nr-1)/14)) write (3,'(2x,"NR",6x,"UVLPT",5x,"VLRO",7x,"VLR",8x,"VLM",3x,"X(1)",4x,"VMION",1x,"VLKA",1x,&
     & "GAMMA1",8x,"CP",10x,"VNR",6x,"FR",9x,"U",8x,"LM",5x,"LM/HP"/1x,"IHV",7x,"VLT",7x,"RHP",6x,"GRAR",7x,"GRAM",3x,"X(2)",3x,&
     & "VMIONP",1x,"VKAP",1x,"GAMMA2",5x,"BETA",9x,"GRAT",6x,"FC",9x,"Z",8x,"VM",6x,"V/VS"/10x,"UVLP",7x,"RHT",4x,"R/RTOT",5x,&
     & "M/MTOT",3x,"X(3)",3x,"VMIONT",1x,"VKAT",1x,"GAMMA3",6x,"VNA",10x,"VNE",6x,"FA",18x,"TM",4x,"L/LEDD"/1x,131("-"))')

  r = 10.d0**(vlr-vlrr)
  xm = 10.d0**(vlm-vlmm)
  xm=1.d0-xm

! Modifications de Schaerer 1990
  rhp = 1.d0/beta_env+vmionp
  rht = vmiont+3.d0-4.d0/beta_env

  ga1 = 1.d0/(rhp+rht*vna)
  ga2 = 1.d0/(1.d0-vna)
  ga3 = 1.d0+vna*ga1
! L / L_Eddington
  d = 10.d0**(vll-vlm+vlka-lgqapicg)

  if (konv == 0) then
    write (3,'(1x,i3,2f10.4,f10.5,f11.5,f7.4,3f7.3,e10.2,f13.5,f8.5,2e10.3,f10.4)') &
      nr,uvlp,vlro,vlr,vlm,x_env(1),vmion,vlka,ga1,cp
    write (3,'(1x,i3,2f10.4,f10.5,f11.5,f7.4,3f7.3,e10.2,f13.5,f8.5,2e10.3,f10.4)') &
      ihv,vlt,rhp,grar,gram,x_env(2),vmionp,vkap,ga2,beta_env,grat
    write (3,'(14x,f10.4,f10.5,e11.4,f7.4,3f7.3,f10.4,41x,f10.4/1x,131("-"))') rht,r,xm,x_env(3),vmiont,vkat,ga3,vna,d
    vne = vnr
    uvlpt = uvlp
    p_turb = 0.0d0
    vm = 0.0d0
    tm = 0.0d0
  else
    xl = 10.d0**vlmix
    vm = 10.d0**(.5d0*(vlgr-vlhp)+vlmix)*sqrt(-rht/8.d0)*u*dwdo2
    tm = xl/vm
    vlv=log10(vm)
    uvlpt=log10(10.d0**uvlp+zeta*10.d0**(vlro+vlv+vlv))
    p_turb = zeta*10.d0**(vlro+vlv+vlv)
          ! Convection: P=P_gaz+P_turb (Maeder 2008, Eq. 5.90)
    vne = grat-(u*dwdo2)**2.d0
    vsvs = vm/(sqrt(ga1*10.d0**(uvlp-vlro)))
    al = xl/10.d0**vlhp
    fr = grat/vnr
    fc=(9.d0/8.d0)*dwdo2*(u*dwdo2)**2.d0/vnr
    fa=1.d0-(fr+fc)
    write (3,'(1x,i3,2f10.4,f10.5,f11.5,f7.4,3f7.3,e10.2,f13.5,f8.5,2e10.3,f10.4)') &
      nr,uvlpt,vlro,vlr,vlm,x_env(1),vmion,vlka,ga1,cp,vnr,fr,u,xl,al
    write (3,'(1x,i3,2f10.4,f10.5,f11.5,f7.4,3f7.3,e10.2,f13.5,f8.5,2e10.3,f10.4)') &
      ihv,vlt,rhp,grar,gram,x_env(2),vmionp,vkap,ga2,beta_env,grat,fc,dwdo2,vm,vsvs
    write (3,'(4x,2f10.4,f10.5,e10.4,f7.4,3f7.3,f10.4,f13.5,f8.5,10x,e10.3,f10.4/1x,131("-"))') &
      uvlp,rht,r,xm,x_env(3),vmiont,vkat,ga3,vna,vne,fa,tm,d
  endif
  if (iprc == 1) then
    if (my >= 1) then
      Var = "logrho"
    else
      Var = "logP"
    endif
    call StoreStructure_env(nr,vlr,vlm,vlt,vlro,p_turb,uvlp,rhp,-rht,vna,ga1,cp,grat,vlka,vkap,vkat, &
         gram,grar,vmion,vmol,vm,tm,x_env(:),Var)
  endif

  return

end subroutine print1
!======================================================================
end module envelope
