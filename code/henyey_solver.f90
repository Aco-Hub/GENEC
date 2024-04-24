module henyey_solver

use evol, only: kindreal
use const, only: um
use inputparam, only: ialflu,ibasnet,irot,itminc,isugi,verbose,EOS,inetwork
use caramodele, only: gms,nwmd
use abundmod,only: x,y3,y,xc12,xc13,xc14,xn14,xn15,xo16,xo17,xo18,xf18,xf19,xne20,xne21,xne22,xna23,xmg24,xmg25,xmg26, &
                   xal26,xal27,xsi28,xprot,xneut,xbid,xbid1,nbelx,nbael,nbzel,abelx,eps,epsy,epsc,epsn,epsyy,epsyc,epsyo, &
                   eg,en

implicit none

integer,save:: nsugi
character(256):: correction_message
logical,save:: henyey_last

private
public:: henyey
public:: nsugi,correction_message,henyey_last

contains
!-----------------------------------------------------------------------
subroutine printhenyey(log_rho,x8,x10,x11,x12,x13,x14,x15,x16,zwi1)
! Printing all the data in henyey.
  use evol,only: ldi
  use const,only: Msol,cst_G,Rsol
  use caramodele, only: gls
  use inputparam,only: idifcon,idiff
  use abundmod,only: snube7,snub8,snu,b11,b33,b34,b112,b113,b114,b115a,b115g,b116,b117a,b117g,b118a,b118g,epcne,eps20,c144, &
    c184,c224,c134,b119a,b119g,b120,b121,b122,b123g,b123a,b124,b125g,b125m,b1mg26,b1al26,b127g,b127a, &
    e24ag,e17ag,e21ag,e18an,e21na,e25an,e20ng,e21ng,e22ng,e23ng,e24ng,e25ng,e26ng,e27ng,e28ng,a26ga,a26gp,e14np,ec14pg,ec14ag, &
    ef18na,e15ag,ef18np,e18pa,ec14ng,e19ap,e14be,e18be,e26be,e18ng
  use EOS, only: psi,gamma1_Timmes,gamma1_dichte,entropy_timmes
  use strucmod, only: p,j,q,t,r,s,vr,radm,zensi,adim,Nabla_mu,m,gravi,H_P,rho,vmyhelio,vmye,xomegafit,xmufit,amu,vmyo
  use rotmod, only: omegi,dlodlr,omegp,vomegi,btotq,omegd,deladv,theta,aux,ur,vcirc,xoblaj
  use magmod,only:D_magx,D_mago,etask,Nmag,bphi,alven,D_circh,qmin
  use diffadvmod,only:D_conv,D_shear,D_eff,D_h,xnabyy,Richardson,K_ther,ucicoe,vcicoe,mtu,npasr
  use timestep, only: alter,dzeit
  use energy,only: nucal
  use PrintAll, only:StoreStructure_int

  implicit none


  integer::ii
  real(kindreal),intent(in):: zwi1,x14,x15,log_rho,x10,x11,x12,x13,x8,x16
  real(kindreal):: vm,logP,logT,logR,vl,vmasse,gmsu,rrsol,gamma1,entropy
  real(kindreal),dimension(ldi):: qv

  character(*),parameter:: headvf='#j   xmr       p           t         r                lr            X              Y&
    &              C12            O16              eps         epsy        epsc          Nabrad       rho       zensi&
    &         epsnu         dkdP        dkdT          dEdP         dEdT         drhodP       delta        psi       eps3a&
    &      epsCO       epsONe     egrav         Nabad       kappa         beta              Y3             C13            N14&
    &            N15              O17            O18            Ne20           Ne22             Mg24             Mg25&
    &             Mg26             mu            omega          Nablamu        Ri             Dconv          Dshear&
    &         Deff          Mr      dlnOmega/dr      K_ther          U               V               D_circ          HP&
    &              g               Dh              Omegp           vr              vomegi          Dmago           Dmagx&
    &           eta             N^2             B_phi           Alfven          q_min           mu_e      F19            Ne21&
    &           Na23           Al26           Al27           Si28alu        C14            F18            nalu           palu&
    &           xbid           neut           Si28           S32            Ar36           Ca40           Ti44           Cr48&           
    &           Cr56           Fe52           Fe53           Fe54           Fe55           Fe56           Co55           Co57&           
    &           Ni56           Btotq          xomegafit      xmufit         vmu           xobla           Gamma'

  character(*),parameter:: headvfgenet48='#j   xmr       p           t         r                lr            X              Y&
    &              C12            O16              eps         epsy        epsc          Nabrad       rho       zensi&
    &         epsnu         dkdP        dkdT          dEdP         dEdT         drhodP       delta        psi       eps3a&
    &      epsCO       epsONe     egrav         Nabad       kappa         beta              Y3             C13            N14&
    &            N15              O17            O18            Ne20           Ne22             Mg24             Mg25&
    &             Mg26             mu            omega          Nablamu        Ri             Dconv          Dshear&
    &         Deff          Mr      dlnOmega/dr      K_ther          U               V               D_circ          HP&
    &              g               Dh              Omegp           vr              vomegi          Dmago           Dmagx&
    &           eta             N^2             B_phi           Alfven          q_min           mu_e      F19            Ne21&
    &           Na23           Al26           Al27           Si28alu        C14            F18            nalu           palu&
    &           xbid           neut           Si28           Si30           P31            S32            S34            Cl35&           
    &           Ar36           Ar38            K39           Ca40           Ca42           Ti44           Ti46           Cr48&
    &           Cr50           Cr56           Fe52           Fe53           Fe54           Fe55           Fe56           Co55& 
    &           Co56           Co57           Ni56           Btotq          xomegafit      xmufit         vmu           xobla&                     
    &           Gamma          entropy'


  vm=1.d0- exp(q(j))             ! Mr/M
  logP=p(j)/um                     ! log P(j)
  logT=t(j)/um                     ! log T(j)
  logR=r(j)/um                     ! log r(j)
  vl=( exp(s(j))-1.d0)*zwi1      ! Lr/L

  call Calcvmyhelio

  !ADAM having some bugs with this
  if (x(j) < 1d-75) then
    x(j) = 0d0
  end if

  if (verbose .or. j <= 1) then
    write(3,'(1x,i3,f9.6,4f8.4,f8.5,1x,f8.5,1x,1pe8.2,1x,1pe8.2,1x,1pe10.2,2e11.2,0p,2f8.3,f7.1/4x,1pe9.2,0p,&
      &4f8.4,2f7.4,1pe9.1,e9.2,e10.2,2e11.2,0p,f8.3,f8.2,f7.3/4x,1p,e9.3,1x,e10.4,1x,e8.2,1x,e8.2,3x,1pe10.2,1p,&
      &1x,e8.2,1x,e8.2,1x,e8.2,2x,e8.2,3x,e8.2,3x,e8.2,10x,e10.4/1x,6(1x,e12.5),/1x,7(1x,e12.5))') &
      j,vm,logP,logT,logR,vl,x(j),y(j),xc12(j),xo16(j),eps(j),epsy(j),epsc(j),radm,log_rho,zensi(j),epsn ,x10,x11,&
      x12,x13,x14,x15,psi,epsyy(j),epsyc(j),epsyo(j),eg,adim,x8,x16,y3(j),xc13(j),xn14(j),xn15(j),xo17(j),xo18(j),&
      xne20(j),xne22(j),xmg24(j),xmg25(j),xmg26(j),omegi(j),Nabla_mu(j),D_h(j),xnabyy(j),D_conv(j),D_shear(j),D_eff(j),D_mago(j), &
      D_magx(j),etask(j),Nmag(j),bphi(j),alven(j),qmin(j)

    if (ialflu == 1) then
      write(3,'(11(1x,e11.5))') xf19(j),xne21(j),xna23(j),xal26(j),xal27(j),xsi28(j),xc14(j),xf18(j),xneut(j),xprot(j),xbid(j)
    endif

    write(3,'(17x,78(i4,")",e9.2))') (ii,abelx(ii,j),ii=1,nbelx)
  endif

  vmasse=vm*gms
 
  if (j == 1) then
    write(29,'(a53)') '# modnb   age                   mtot  nbshell  deltat'
    write(29,'(i6,1x,1pe20.13,0p,1x,f10.5,i7,1pe20.13)') nwmd,alter,gms,m,dzeit
    write(29,'(a)')trim(headvfgenet48)

  endif

  if ((irot == 0.and.idifcon == 0) .or. (irot==1.and.idiff==0)) then
    gmsu=gms*Msol
    qv(j)=(1.d0-exp(q(j)))*gmsu
    if (j == m) then
      gravi(j)=0.0d0
    else
      gravi(j)=cst_G*qv(j)/(exp(2.d0*r(j)))
    endif
    if (gravi(j) /= 0.d0) then
      H_P(j)=exp(p(j))/(exp(rho(j))*gravi(j))
    endif
  endif


  if ( ( (log_rho) .lt. 2.8d0) .or. (logT .lt. 7.55d0) )  then
    entropy= 0.0d0
    gamma1=gamma1_dichte
  ELSE
    ! write(*,*) "ENTROPY", entropy_timmes,log_rho,logT
    entropy=entropy_timmes
    gamma1=gamma1_Timmes
  endif

  if (EOS ==0) then
    gamma1 = gamma1_dichte
  endif
!23 --> 15 if lower network, to automize
  write(29,'(i4,3(f10.7,1x),f14.11,1x,e14.6,4(1x,e14.7),3x,1p,3(e11.4,1x),2x,e11.4,1x,0pf11.6,1x,1pe12.5,1x,e11.4,&
    &3x,6(e12.5,1x),e9.2,1x,e9.2,1x,e10.2,1x,e11.2,3x,4(e12.5,1x),5x,0p,4(e14.7,1x),2x,4(e14.7,1x),2x,3(e14.7,3x),&
    &f9.6,2x,1p,6(3x,e12.5),1x,0p,f9.4,18(1x,e15.8),1x,f9.6,1p,11(1x,e14.7),26(1x,e14.7),5(1x,e14.7),1x,0pf9.6,1x,e14.7)') & 
    j,vm,logP,logT,logR,vl,x(j),y(j),xc12(j),xo16(j),eps(j),epsy(j),epsc(j),radm,log_rho,zensi(j),epsn ,x10,x11,x12,x13,x14, &
    x15,psi,epsyy(j),epsyc(j),epsyo(j),eg,adim,x8,x16,y3(j),xc13(j),xn14(j),xn15(j),xo17(j),xo18(j),xne20(j),xne22(j), &
    xmg24(j),xmg25(j),xmg26(j),vmyhelio(j),omegi(j),Nabla_mu(j),Richardson(j),D_conv(j),D_shear(j),D_eff(j),vmasse, &
    dlodlr(j),K_ther(j),ucicoe(j),vcicoe(j),D_circh(j),H_P(j),gravi(j),D_h(j),omegp(j),vr(j),vomegi(j),D_mago(j), &
    D_magx(j),etask(j),Nmag(j),bphi(j),alven(j),qmin(j),vmye,xf19(j),xne21(j),xna23(j), &
    xal26(j),xal27(j),xsi28(j),xc14(j),xf18(j),xneut(j),xprot(j),xbid(j),(abelx(ii,j),ii=1,nbelx),btotq(j), &
    exp(xomegafit(j)),exp(xmufit(j)),1.d0/amu(m-j+1),xoblaj,gamma1,entropy

  call StoreStructure_int(j,logR,vm*gms,logT,log_rho,logP,x14,x15,adim,radm,x8,vl*gls,x10,x11,en,x12,x13, &
                          x(j),y(j),omegi(j),vmyhelio(j),vmyo)

  if (j == mtu) then
    write(39,'(1x," masse=",f10.6," age=",e20.6)') gms,alter
  endif
  if (j >= mtu.and.j <= npasr) then
    rrsol=10.d0**logR/Rsol
    write(39,'(i4,12(1x,1pe15.8))') j,vmasse,rrsol,vomegi(j),omegp(j),omegd(j),omegi(j),deladv(j),theta(j),aux(j),ur(j), &
                                    vcirc(j),Nabla_mu(j)
  endif

! Calcul du flux de neutrinos
  if (j == m) then
    write(3,'(/2x,a,2x,e8.2,2x,a/20x,a,6x,e8.2,2x,a/20x,a/20x,a,5x,e8.2,2x,a)') 'Flux de neutrinos : BERYLIUM',snube7, &
                                            'SNU',': BORE',snub8,'SNU','-------------------------',': TOTAL',snu,'SNU'
  endif

  if (j == m) then
    write(3,'(/1x,10e13.4/9e13.4)') b11(j),b33(j),b34(j),b112(j),b113(j),b114(j),b115a(j),b115g(j),b116(j),b117a(j),b117g(j), &
                                    b118a(j),b118g(j),epcne(j),eps20(j),c144(j),c184(j),c224(j),c134(j)

    if (ialflu == 1) then
      write(3,'(9e13.4,/9e13.4/9e13.4/9e13.4/8e13.4)') b119a(j),b119g(j),b120(j),b121(j),b122(j),b123g(j),b123a(j),b124(j), &
        b125g(j),b125m(j),b1mg26(j),b1al26(j),b127g(j),b127a(j),e24ag(j),e17ag(j),e21ag(j),ec14ag(j),e15ag(j),e18an(j), &
        e25an(j),e21na(j),ef18na(j),e20ng(j),e21ng(j),e22ng(j),e23ng(j),e24ng(j),e25ng(j),e26ng(j),e27ng(j),e28ng(j), &
        e18ng(j),ec14ng(j),e14np(j),ef18np(j),a26ga(j),a26gp(j),ec14pg(j),e18pa(j),e19ap(j),e14be(j),e18be(j),e26be(j)
    endif

! Calcul des neutrinos solaires
    if (gms == 1.d0) call nucal
  endif

  return

  end subroutine printhenyey

!***********************************************************************
subroutine Calcvmyhelio

  use const,only: cst_avo,pi,uastr
  use caramodele,only: glm
  use inputparam,only: z
  use abundmod,only: b34neu,ybe7,tauxbe7el,tauxbe7pg,yb8,xnube7,xnub8,snube7,snub8,xnbrbe7,xnbrb8,snu,zabelx
  use strucmod,only: j,m,q,vmyhelio

  implicit none

  integer::k=0,ii=0

! Calcul de 7Be et 8B
!***  xnu(j) en nbre de neutrinos par seconde et par gramme
  if (b34neu(j) /= 0.d0) then
    ybe7(j)=y3(j)/3.d0*y(j)/4.d0*b34neu(j)/(0.5d0*(x(j)+1)*tauxbe7el(j)+x(j)*tauxbe7pg(j))
    yb8(j)=1.11d0*x(j)*ybe7(j)*tauxbe7pg(j)
    xnube7(j)=cst_avo*0.5d0*(x(j)+1)*ybe7(j)*tauxbe7el(j)
    xnub8(j)=cst_avo*ybe7(j)*x(j)*tauxbe7pg(j)
  endif
  if (j == m .and. itminc == 1) then
    snube7=0.d0
    snub8=0.d0
    do k=2,m-1
     xnbrbe7(k)=0.125d0/pi/uastr**2.d0*xnube7(k)*exp(glm)*(exp(q(k+1))-exp(q(k-1)))
     xnbrb8(k)=0.125d0/pi/uastr**2.d0*xnub8(k)*exp(glm)*(exp(q(k+1))-exp(q(k-1)))
     snube7=snube7+xnbrbe7(k)*2.38d-10
     snub8=snub8+xnbrb8(k)*1.08d-6
    enddo

    xnbrbe7(1)=0.125d0/pi/uastr**2.d0*xnube7(1)*exp(glm)*(exp(q(2))+exp(q(1)))
    xnbrb8(1)=0.125d0/pi/uastr**2.d0*xnub8(1)*exp(glm)*(exp(q(2))+exp(q(1)))
    xnbrbe7(m)=0.125d0/pi/uastr**2.d0*xnube7(m)*exp(glm)*(1.d0-exp(q(m-1)))
    xnbrb8(m)=0.125d0/pi/uastr**2.d0*xnub8(m)*exp(glm)*(1.d0-exp(q(m-1)))
    snube7=snube7+(xnbrbe7(1)+xnbrbe7(m))*2.38d-10
    snub8=snub8+(xnbrb8(1)+xnbrb8(m))*1.08d-6
    snu = snube7 +snub8
  endif

  vmyhelio(j)=2.d0*x(j)+y3(j)+3.d0/4.d0*y(j)+7.d0/12.d0*xc12(j)+7.d0/13.d0*xc13(j)+8.d0/14.d0*xn14(j)+8.d0/15.d0*xn15(j)+ &
              9.d0/16.d0*xo16(j)+9.d0/17.d0*xo17(j)+9.d0/18.d0*xo18(j)+11.d0/20.d0*xne20(j)+11.d0/22.d0*xne22(j)+13.d0/24.d0* &
              xmg24(j)+13.d0/25.d0*xmg25(j)+13.d0/26.d0*xmg26(j)
  if (ialflu == 1) then
    vmyhelio(j)=vmyhelio(j)+10.d0/19.d0*xf19(j)+11.d0/21.d0*xne21(j)+12.d0/23.d0*xna23(j)+14.d0/26.d0*xal26(j)+ &
                   14.d0/27.d0*xal27(j)+15.d0/28.d0*xsi28(j)+0.5d0*xc14(j)+10.d0/18.d0*xf18(j)
  endif

!  z elements taken ~ as Ni56: A=56 but (nbzel(ii)+1.)/nbael(ii)-->0.5
  do ii=1,nbelx
   vmyhelio(j)= vmyhelio(j)+abelx(ii,j)*(nbzel(ii)+1.d0)/nbael(ii)
  enddo
  vmyhelio(j)= vmyhelio(j)+0.5d0*zabelx

  vmyhelio(j)=1.d0/vmyhelio(j)

  return

  end subroutine Calcvmyhelio

!***********************************************************************
subroutine gisu

  use const,only: year,cst_G,Msol,pi
  use inputparam,only: ipop3,iledou,iover
  use caramodele,only: hh6
  use abundmod,only: egp,egt,epsn1,epsp1,enuep1,epsp,enuep,epst1,enuet1,epst,enuet,enue,egp1,egt1
  use equadiffmod,only: g1,g1s1,g1s,g1p1,g1p,g1t1,g1t,ccg2,ccg3,g2,g2r1,g2r,g2p1,g2p,g3,g3p1,g3p,g3t1,g3t,g3r1,g3r, &
                        g4,g4r1,g4r,g4p1,ccg1,g4p,g4t1,g4t,g4s1,g4s
  use EOS,only: rh,rht,rhp,rh1,rht1,rhp1
  use strucmod,only: m,j,j1,beta,p,vp,vt,adi,adip,adit,beta1,adi1,adip1,adit1,r,q,s,t,e,adim,radm,xbruj1,rad1,rad,adgrad, &
                     Nabla_mu,zradm,capp1,capp,capt1,capt,zensi,gradp,gradt
  use rotmod,only: omegi,xoblaj
  use timestep,only: alter,dzeit
  use convection,only: xover
  use geomod, only: rpsi_min,rpsi_max,geocalc,geom
  use SmallFunc,only: exphi

  implicit none

  real(kindreal):: egx,egy,egxp=0.d0,egyp=0.d0,egxt=0.d0,egyt=0.d0,egx1=0.d0,egy1=0.d0,egx1p1=0.d0,egy1p1=0.d0,egx1t1=0.d0, &
    egy1t1=0.d0,xpsij,xpsij1,xfpj,xraj,dxfpj,dxraj,xfpj1,xraj1,dxfpj1,dxraj1,xfp=0.d0,xra,djxfp=0.d0,dj1xfp=0.d0,djxra,dj1xra, &
    xgpsi,xggj,xgpsi1,xggj1,xggjm,xhpjm,hnenn,d1,ff1,f1,hfak,hfakn,d2,f2,hfak2,d3,f3,xnbrun,xbruj,admu,admu1,wrm,d4,f4,xhpj1, &
    tnorm,g2a,g2ap1,g2ap,g2ar1,g2ar,ff10,f10,f41,f40,f410,dtat,dtap,ega,egb,dtat1,dtap1,ega1,egb1,eg1,eg1p1=0.d0,eg1t1=0.d0, &
    en1=0.d0,sugib1,sugib4

! coef. beta1/4 in Sugimoto's method
  select case (isugi)
    case (1)
      sugib1=1.d0
    case (2)
      sugib1=0.d0
    case (3)
      sugib1=0.75d0
    case (4)
      sugib1=0.25d0
    case (5)
      sugib1=0.50d0
    case default
      stop 'problem in girsu: isugi not well defined'
  end select

  if (j == m-2 .and. nsugi >= 0) then
    if (epsc(j) == 0.0d0 .or. isugi == 5) then
      nsugi=m
    else
      nsugi =nsugi-1
    endif
    if (verbose) then
      write(*,*) 'nsugi= ',nsugi
    endif
  endif
  if (j < nsugi) sugib1=0.50d0

  sugib4=1.d0-sugib1
  if (alter <= 15.d0*dzeit/year) then
    if (j1 == m-10 .and. verbose) write(*,*)'isugi,sugib1, sugib4= ', isugi,sugib1, sugib4
  endif

  if (sugib1 /= 0.50d0) then
    dtap = -4.d0/beta/beta*(1.d0-beta)
    dtat = -4.d0*dtap
    ega  = exp(p(j)-rh)/dzeit
    egb  = vp(j)-vt(j)/adi
    eg   = ega*(-rht)*egb
    egp=(1.d0-rhp)*eg+ega*(dtap*egb-rht*(1.d0+vt(j)*adip/adi/adi))
    egt=-rht*eg+ega*(dtat*egb+rht*(1.d0/adi-vt(j)*adit/adi/adi))
    egx  = ega*(-rht)*vp(j)
    egy  = ega*rht*vt(j)/adi
    egxp  = (1.d0-rhp)*egx+ ega*(-rht+dtap*vp(j))
    egyp  = (1.d0-rhp)*egy+ ega*(-rht*vt(j)*adip/adi/adi-dtap*vt(j)/adi)
    egxt  = -rht*egx+ ega*dtat*vp(j)
    egyt  = -rht*egy + ega*(rht/adi-rht*vt(j)*adit/adi/adi -dtat*vt(j)/adi)
    dtap1 = -4.d0/beta1/beta1*(1.d0-beta1)
    dtat1 = -4.d0*dtap1
    ega1  = exp(p(j1)-rh1)/dzeit
    egb1  = vp(j1)-vt(j1)/adi1
    eg1   = ega1*(-rht1)*egb1
    eg1p1  = (1.d0-rhp1)*eg1+ega1*(dtap1*egb1-rht1*(1.d0+vt(j1)*adip1/adi1/adi1))
    eg1t1  = -rht1*eg1+ega1*(dtat1*egb1+rht1*(1.d0/adi1-vt(j1)*adit1/adi1/adi1))
    egx1  = ega1*(-rht1)*vp(j1)
    egy1  = ega1*rht1*vt(j1)/adi1
    egx1p1  = (1.d0-rhp1)*egx1+ ega1*(-rht1+dtap1*vp(j1))
    egy1p1  = (1.d0-rhp1)*egy1 + ega1*(-rht1*vt(j1)*adip1/adi1/adi1-dtap1*vt(j1)/adi1)
    egx1t1  = -rht1*egx1+ ega1*dtat1*vp(j1)
    egy1t1  = -rht1*egy1+ ega1*(rht1/adi1-rht1*vt(j1)*adit1/adi1/adi1 -dtat1*vt(j1)/adi1)
    en   = (eps(j)+epsy(j)+epsc(j))
    en1   = (eps(j1)+epsy(j1)+epsc(j1))
    if (ipop3 == 0) then
      if ((eps(j) /= 0.0d0 .and. epsy(j) /= 0.0d0) .or. (eps(j) /= 0.0d0 .and. epsc(j) /= 0.0d0) &
             .or. (epsy(j) /= 0.0d0 .and. epsc(j) /= 0.0d0)) then
        stop 'stop in gisu, line 306'
      endif
    else   ! ipop3=1
      if ((eps(j) /= 0.0d0 .and. epsc(j) /= 0.0d0) .or. (epsy(j) /= 0.0d0 .and. epsc(j) /= 0.0d0)) then
        stop 'stop in gisu, line 310'
      endif
    endif
  endif

  if (irot==1 .and. omegi(j)>1.d-20 .and. omegi(j1)>1.d-20) then
    if (j1 /= m) then
      xpsij=exp(r(j))/((cst_G*exp(ccg1)*(1.d0-exp(q(j))))/(omegi(j)*omegi(j)))**(1.d0/3.d0)
      xpsij1=exp(r(j1))/((cst_G*exp(ccg1)*(1.d0-exp(q(j1))))/(omegi(j1)*omegi(j1)))**(1.d0/3.d0)
      if (xpsij > rpsi_max) xpsij=0.9999999999d0*rpsi_max
      if (xpsij1 > rpsi_max) xpsij1=0.9999999999d0*rpsi_max
      call geom(xpsij,xfpj,xraj,dxfpj,dxraj)
      call geom(xpsij1,xfpj1,xraj1,dxfpj1,dxraj1)
! ancien appel: call geocrit(xpsij,xoblaj)
      call geocalc(xpsij,xoblaj,3)
      xfp=sqrt(xfpj*xfpj1)
      xra=sqrt(xraj*xraj1)
      djxfp=dxfpj*xpsij*xfpj1/(2.d0*xfp)
      dj1xfp=dxfpj1*xpsij1*xfpj/(2.d0*xfp)
      djxra=dxraj*xpsij*xraj1/(2.d0*xra)
      dj1xra=dxraj1*xpsij1*xraj/(2.d0*xra)
    else
      xpsij=exp(r(j))/((cst_G*exp(ccg1)*(1.d0-exp(q(j))))/(omegi(j)*omegi(j)))**(1.d0/3.d0)
      if (xpsij > rpsi_max) xpsij=0.9999999999d0*rpsi_max
      call geom(xpsij,xfpj,xraj,dxfpj,dxraj)
! ancien appel: call geocrit(xpsij,xoblaj)
      call geocalc(xpsij,xoblaj,3)
      xfpj1=1.0d0
      xraj1=1.0d0
      dxfpj1=0.0d0
      dxraj1=0.0d0
      xfp=sqrt(xfpj*xfpj1)
      xra=sqrt(xraj*xraj1)
      djxfp=dxfpj*xpsij*xfpj1/(2.d0*xfp)
      dj1xfp=dxfpj1*xpsij1*xfpj/(2.d0*xfp)
      djxra=dxraj*xpsij*xraj1/(2.d0*xra)
      dj1xra=dxraj1*xpsij1*xraj/(2.d0*xra)
    endif

! calcul de la gravite. Intervient dans la frequence
! de Brunt-Vaisala avec effet de rotation inclu
    if (j1 /= m) then
      xpsij=exp(r(j))/((cst_G*(1.d0-exp(q(j)))*gms*Msol)/(omegi(j)*omegi(j)))**(1.d0/3.d0)
      xpsij1=exp(r(j1))/((cst_G*(1.d0-exp(q(j1)))*gms*Msol)/(omegi(j1)*omegi(j1)))**(1.d0/3.d0)

      if (xpsij >= rpsi_min .and. xpsij1 >= rpsi_min) then
        if (xpsij > rpsi_max) xpsij=0.9999999999d0*rpsi_max
! ancien appel: call geograv(xpsij,xgpsi)
        call geocalc(xpsij,xgpsi,1)
        xggj=cst_G*(1.d0-exp(q(j)))*gms*Msol*omegi(j)*omegi(j)*omegi(j)*omegi(j)
        xggj=xggj**(1.d0/3.d0)*xgpsi
        if (xpsij1 > rpsi_max) xpsij1=0.9999999999d0*rpsi_max
! ancien appel: call geograv(xpsij1,xgpsi1)
        call geocalc(xpsij1,xgpsi1,1)
        xggj1=cst_G*(1.d0-exp(q(j1)))*gms*Msol*omegi(j1)*omegi(j1)*omegi(j1)*omegi(j1)
        xggj1=xggj1**(1.d0/3.d0)*xgpsi1
        xggjm=sqrt(xggj*xggj1)
      else
        xggj=cst_G*(1.d0-exp(q(j)))*gms*Msol/exp(2.d0*r(j))
        xggj1=cst_G*(1.d0-exp(q(j1)))*gms*Msol/exp(2.d0*r(j1))
        xggjm=sqrt(xggj*xggj1)
      endif
    else
! on utilise ci-dessous: r=3/(4pi rhoc)^(1/3)*m_m^(1/3) et m_m=M_(m-1)/2
! Kippenhahn & Weigert, p. 68, Eq. (10.3)
      xggjm=cst_G*(4.d0*pi/3.d0)**(2.d0/3.d0)*2.d0**(-1.d0/3.d0)*exp(2.d0*rh1/3.d0)*((1.d0-exp(q(j)))*gms*Msol)**(1.d0/3.d0)
    endif

! calcul de l'echelle de pression
    if (j1 /= m) then
      xhpjm=exp(0.5d0*(p(j)+p(j1)))/(exp(0.5d0*(rh+rh1))*xggjm)
      xhpj1=exp(p(j1))/(exp(rh1)*xggj1)
    else
      xhpjm=exp(0.5d0*(p(j)+p(j1)))/(exp(0.5d0*(rh+rh1))*xggjm)
! au centre l'echelle de pression diverge, on prend xhpjm
      xhpj1=xhpjm
    endif
  else
    xoblaj = 1.0d0
  endif   !   irot

  hnenn=1.d0/(q(j1)-q(j))
  d1=hnenn*(s(j1)-s(j))

  if (sugib1 /= 0.50d0) then
    ff1=exp(ccg1-hh6+q(j)-s(j))   *sugib1
    ff10=exp(ccg1-hh6+q(j1)-s(j1))*(1-sugib1)
    f1=(en+egx+egy+epsn)*ff1
    f10=(en1+egx1+egy1+epsn1)*ff10
    g1=d1+f1+f10
    if (s(j) < -1.d-3) write(*,*) 's neg.',j,s(j),en,egx,egy,epsn
    g1s1=hnenn-f10
    g1s=-hnenn-f1
    g1p1=ff10*(en1*epsp1+eg1p1+epsn1*enuep1)
    g1p=ff1*(en*epsp+egp+epsn*enuep)
    g1t1=ff10*(en1*epst1+eg1t1+epsn1*enuet1)
    g1t=ff1*(en*epst+egt+epsn*enuet)
    g1p1=ff10*(en1*epsp1+egx1p1+egy1p1+epsn1*enuep1)
    g1p=ff1*(en*epsp+egxp+egyp+epsn*enuep)
    g1t1=ff10*(en1*epst1+egx1t1+egy1t1+epsn1*enuet1)
    g1t=ff1*(en*epst+egxt+egyt+epsn*enuet)
    if (abs(g1p)>HUGE(g1p) .or. abs(g1s)>HUGE(g1s)) then
      write(*,*)'j1,g1s,hnenn,f10,g1p,ff10,en,epsp,egp,epsn,enuep:',&
         j1,g1s,hnenn,f10,g1p,ff10,en,epsp,egp,epsn,enuep
    endif
  else
    ff1=exp(ccg1-hh6+0.5d0*(q(j)+q(j1)-s(j)-s(j1)))
    f1=(en+eg-enue)*ff1
    g1=d1+f1
    g1s1=hnenn-0.5d0*f1
    g1s=-hnenn-0.5d0*f1
    hfak=en/2.d0
    hfakn=-enue*0.5d0
    g1p1=ff1*(hfak*epsp1+egp1+hfakn*enuep1)
    g1p=ff1*(hfak*epsp+egp+hfakn*enuep)
    g1t1=ff1*(hfak*epst1+egt1+hfakn*enuet1)
    g1t=ff1*(hfak*epst+egt+hfakn*enuet)
    if (abs(g1p)>HUGE(g1p) .or. abs(g1s)>HUGE(g1s)) then
      write(*,*)'j1,g1s,hnenn,f1,g1p,ff1,hfak,epsp,egp,hfakn,enuep,enue:',&
         j1,g1s,hnenn,f1,g1p,ff1,hfak,epsp,egp,hfakn,enuep,enue
    endif

  endif

  tnorm=dzeit/1.0d+14
  g1  =tnorm*g1
  g1s =tnorm*g1s
  g1s1=tnorm*g1s1
  g1p =tnorm*g1p
  g1p1=tnorm*g1p1
  g1t =tnorm*g1t
  g1t1=tnorm*g1t1
  d2=hnenn*(p(j1)-p(j))
  f2=-e(j)*exp(ccg2-2.d0*(r(j)+r(j1))+0.5d0*(q(j)+q(j1)-p(j)-p(j1)))
  g2a=0.0d0
  g2ap1= 0.0d0
  g2ap = 0.0d0
  g2ar1= 0.0d0
  g2ar = 0.0d0

  if (irot==0 .or. omegi(j)<=1.d-20 .or. omegi(j1)<=1.d-20) then
    g2=d2+f2+g2a
    g2r1=-2.d0*f2 +g2ar1
    g2r=g2r1 +g2ar
    hfak2=0.5d0*f2
    g2p1=hnenn-hfak2 + g2ap1
    g2p=-hnenn-hfak2 + g2ap
  else
    g2=d2+f2*xfp+g2a
    g2r1=-2.d0*f2*xfp+f2*dj1xfp +g2ar1
    g2r=-2.d0*f2*xfp+f2*djxfp +g2ar
    hfak2=0.5d0*f2*xfp
    g2p1=hnenn-hfak2 + g2ap1
    g2p=-hnenn-hfak2 + g2ap
  endif

  d3=hnenn*(r(j1)-r(j))
  f3=exp(ccg3-1.5d0*(r(j)+r(j1))+0.5d0*(q(j)+q(j1)-rh1-rh))
  g3=d3+f3
  hfak=-0.5d0*f3
  g3p1=hfak*rhp1
  g3p=hfak*rhp
  g3t1=hfak*rht1
  g3t=hfak*rht
  g3r1=hnenn-1.5d0*f3
  g3r=-hnenn-1.5d0*f3

! critere de LEDOUX ou de SCHWARZSCHILD
! Si l'on veut remettre le critere de Hoiland-Solberg
! il faut enlever les commentaires du type H-S
!H-S      if (irot == 0) then
  if (iledou == 0) then
    xnbrun=adim-radm
    xbruj1=adi1-rad1
    xbruj =adi-rad
    xnbrun=xbruj*sugib4+xbruj1*(1.d0-sugib4)
    adgrad(j)=adi-rad
  else
    admu=0.5d0*(adi+Nabla_mu(j)/(-rht)+adi1+Nabla_mu(j1)/(-rht1))
    admu1=adi1+Nabla_mu(j1)/(-rht1)
    xnbrun=admu-radm
    xbruj1=admu1-rad1
    xbruj=adi +Nabla_mu(j)/(-rht)-rad
    xnbrun=xbruj*sugib4+xbruj1*(1.d0-sugib4)
    adgrad(j)=adi +Nabla_mu(j)/(-rht)-rad
  endif
!H-S      endif
!H-S      if (irot == 1.and.iledou == 0) then
!H-S      deltam=-0.5*(rht+rht1)
!H-S      xommoy=0.5*(omegi(j)+omegi(j1))
!H-S      dxommo=0.5*(dlodlr(j)+dlodlr(j1))
!H-S      xnbrun=xggjm*deltam/xhpjm*(adim-radm)
!H-S     &       +1.633*xommoy*xommoy*(2.+dxommo)
!H-S      xbruj1=xggj1*(-rht1)/xhpj1*(adi1-rad1)
!H-S     &       +1.633*omegi(j1)*omegi(j1)*(2.+dlodlr(j1))
!H-S      endif
!H-S      if (irot == 1.and.iledou == 1) then
!H-S      admu=adim+(xdmudp(j)+xdmudp(j1))/(-rht-rht1)
!H-S      admu1=adi1+xdmudp(j1)/(-rht1)
!H-S      deltam=-0.5*(rht+rht1)
!H-S      xommoy=0.5*(omegi(j)+omegi(j1))
!H-S      dxommo=0.5*(dlodlr(j)+dlodlr(j1))
!H-S      xnbrun=xggjm*deltam/xhpjm*(admu-radm)
!H-S     &       +1.633*xommoy*xommoy*(2.+dxommo)
!H-S      xbruj1=xggj1*(-rht1)/xhpj1*(admu1-rad1)
!H-S     &       +1.633*omegi(j1)*omegi(j1)*(2.+dlodlr(j1))
!H-S      endif
  wrm=0.5d0*(exp(r(j1))+exp(r(j)))
  if (xnbrun < 0.d0 .or. (iover == 1 .and. xover > wrm)) then
    if (irot==0 .or. omegi(j)<=1.d-20 .or. omegi(j1)<=1.d-20) then

! convective zone
      d4=hnenn*(t(j1)-t(j))
      f41=-(1.d0-exp(q(j1)))*exp(ccg2+q(j1)-p(j1)-4.d0*r(j1))*(1.d0-sugib4)
      f410=-(1.d0-exp(q(j)))*exp(ccg2+q(j)-p(j)-4.d0*r(j))*sugib4
      f4=adi1*f41
      f40=adi*f410
      g4=d4+f4+f40
      g4r1=-4.d0*f4
      g4r=-4.d0*f40
      g4p1=f41*(adip1-adi1)
      g4p=f410*(adip-adi)
      g4t1=hnenn+f41*adit1
      g4t=-hnenn+f410*adit
      g4s1=0.0d0
      g4s=0.0d0
    else
      d4=hnenn*(t(j1)-t(j))
! f41(rotation)= f41(no rotation)*xfpj1
      f41=-(1.d0-exp(q(j1)))*exp(ccg2+q(j1)-p(j1)-4.d0*r(j1))*xfpj1*(1.d0-sugib4)
      f410=-(1.d0-exp(q(j)))*exp(ccg2+q(j)-p(j)-4.d0*r(j))*xfpj*sugib4
      f4=adi1*f41
      f40=adi*f410
      g4=d4+f4+f40
      g4r1=-4.d0*f4+ f4*dxfpj1*xpsij1/xfpj1
      g4r= -4.d0*f40+f40*dxfpj*xpsij/xfpj
      g4p1=f41*(adip1-adi1)
      g4p=f410*(adip-adi)
      g4t1=hnenn+f41*adit1
      g4t=-hnenn+f410*adit
      g4s1=0.d0
      g4s=0.d0
    endif
  else
    d4=hnenn*(t(j1)-t(j))
    f4=zradm*f2
    if (irot==0 .or. omegi(j)<=1.d-20 .or. omegi(j1)<=1.d-20) then

! radiative zone
      f41=-(1.d0-exp(q(j1)))*exp(ccg2+q(j1)-p(j1)-4.d0*r(j1))*(1.d0-sugib4)
      f410=-(1.d0-exp(q(j)))*exp(ccg2+q(j)-p(j)-4.d0*r(j))*sugib4
      f4=rad1*f41
      f40=rad*f410
      g4=d4+f4+f40
      g4r1=-4.d0*f4
      g4r=-4.d0*f40
      g4s1=f4/exphi(-s(j1))
      g4s=f40/exphi(-s(j))
      g4p1=capp1*f4
      g4p=capp*f40
      g4t1=hnenn+f4*(capt1-4.d0)
      g4t=-hnenn+f40*(capt-4.d0)
    else
      f41=-(1.d0-exp(q(j1)))*exp(ccg2+q(j1)-p(j1)-4.d0*r(j1))*xfpj1*xraj1*(1.d0-sugib4)
      f410=-(1.d0-exp(q(j)))*exp(ccg2+q(j)-p(j)-4.d0*r(j))*xfpj*xraj*sugib4
      f4=rad1*f41
      f40=rad*f410
      g4=d4+f4+f40
      g4r1=-4.d0*f4 +f4*xpsij1*(dxraj1/xraj1+dxfpj1/xfpj1)
      g4r=-4.d0*f40 +f40*xpsij*(dxraj/xraj+dxfpj/xfpj)
      g4s1=f4/exphi(-s(j1))
      g4s=f40/exphi(-s(j))
      g4p1=capp1*f4
      g4p= capp*f40
      g4t1=hnenn+f4*(capt1-4.d0)
      g4t=-hnenn+f40*(capt-4.d0)
    endif

    if (xbruj > 0.d0) then
      zensi(j)=-zensi(j)
    endif
  endif

  if (irot==0 .or. omegi(j)<=1.d-20 .or. omegi(j1)<=1.d-20) then
    gradp(j)=-f2
    gradt(j)=-f4
  else
    gradp(j)=-f2*xfp
    gradt(j)=-f4
  endif

  return

end subroutine gisu
!***********************************************************************
subroutine gi

  use const,only: year,cst_G,Msol,pi,cst_a,cst_c
  use inputparam,only: ipop3,iledou,iover,idebug
  use caramodele,only: hh6
  use abundmod,only: egp,egt,epsp1,enuep1,epsp,enuep,epst1,enuet1,epst,enuet,enue,egp1,egt1
  use equadiffmod,only: g1,g1s1,g1s,g1p1,g1p,g1t1,g1t,ccg2,ccg3,g2,g2r1,g2r,g2p1,g2p,g3,g3p1,g3p,g3t1,g3t,g3r1,g3r, &
                        g4,g4r1,g4r,g4p1,ccg1,g4p,g4t1,g4t,g4s1,g4s
  use EOS,only: rh,rht,rhp,rh1,rht1,rhp1
  use strucmod,only: m,j,j1,p,adi,adip,adit,adi1,adip1,adit1,r,q,s,t,e,adim,radm,xbruj1,rad1,rad, &
                     adgrad,Nabla_mu,zradm,capp1,capp,capt1,capt,zensi,gradp,gradt,cap,zrad1,zrad
  use rotmod,only: omegi,xoblaj
  use timestep,only: alter,dzeit
  use convection,only: xover
  use geomod, only: rpsi_min,rpsi_max,geocalc,geom
  use SmallFunc,only: exphi

  implicit none

  real(kindreal):: xpsij,xpsij1,xfpj,xraj,dxfpj,dxraj,xfpj1,xraj1,dxfpj1,dxraj1,xfp=0.d0,xra=0.d0,djxfp=0.d0,dj1xfp=0.d0, &
                   djxra=0.d0,dj1xra=0.d0,xgpsi,xggj,xgpsi1,xggj1,xggjm,xhpjm,hnenn,d1,ff1,f1,hfak,hfakn,d2,f2,hfak2,d3,f3, &
                   xnbrun,xbruj,admu,admu1,wrm,d4,f4,xhpj1,tnorm,f44=0.d0,f44rj,f44rj1,tdifth,tdifth2

! calcul du temps dans diffusion thermique dans une couche et dans
! toute l'etoile.
! On a : tau = r^2 / K, K = 4act^3 / (3kappa rho c_p)
  tdifth=-rht/adi*exp(log(3.d0)-log(4.d0)-log(cst_a)-log(cst_c) + cap+rh+p(j)-4.d0*t(j))*(exp(r(j))-exp(r(j1)))**2.d0
  tdifth2=-rht/adi*exp(log(3.d0)-log(4.d0)-log(cst_a)-log(cst_c) + cap+rh+p(j)-4.d0*t(j)) *(exp(r(j)))**2.d0

  if (dzeit < tdifth2) then
    if (isugi >= 1) then
      if (alter >= 9.d0*dzeit/year) then
        if (idebug>1) then
          write(*,*) 'call gisu'
        endif
        call gisu
        return
      endif
    endif
  endif

  if (irot == 1) then
    if (j1 /= m) then
      xpsij=exp(r(j))/((cst_G*exp(ccg1)*(1.d0-exp(q(j))))/(omegi(j)*omegi(j)))**(1.d0/3.d0)
      xpsij1=exp(r(j1))/((cst_G*exp(ccg1)*(1.d0-exp(q(j1))))/(omegi(j1)*omegi(j1)))**(1.d0/3.d0)
      if (xpsij > rpsi_max) xpsij=0.9999999999d0*rpsi_max
      if (xpsij1 > rpsi_max) xpsij1=0.9999999999d0*rpsi_max
      call geom(xpsij,xfpj,xraj,dxfpj,dxraj)
      call geom(xpsij1,xfpj1,xraj1,dxfpj1,dxraj1)
! ancien appel: call geocrit(xpsij,xoblaj)
      call geocalc(xpsij,xoblaj,3)
      xfp=sqrt(xfpj*xfpj1)
      xra=sqrt(xraj*xraj1)
      djxfp=dxfpj*xpsij*xfpj1/(2.d0*xfp)
      dj1xfp=dxfpj1*xpsij1*xfpj/(2.d0*xfp)
      djxra=dxraj*xpsij*xraj1/(2.d0*xra)
      dj1xra=dxraj1*xpsij1*xraj/(2.d0*xra)
    else
      xpsij=exp(r(j))/((cst_G*exp(ccg1)*(1.d0-exp(q(j))))/(omegi(j)*omegi(j)))**(1.d0/3.d0)
      if (xpsij > rpsi_max) xpsij=0.9999999999d0*rpsi_max
      call geom(xpsij,xfpj,xraj,dxfpj,dxraj)
! ancien appel: call geocrit(xpsij,xoblaj)
      call geocalc(xpsij,xoblaj,3)
      xfpj1=1.0d0
      xraj1=1.0d0
      dxfpj1=0.0d0
      dxraj1=0.0d0
      xfp=sqrt(xfpj*xfpj1)
      xra=sqrt(xraj*xraj1)
      djxfp=dxfpj*xpsij*xfpj1/(2.d0*xfp)
      dj1xfp=dxfpj1*xpsij1*xfpj/(2.d0*xfp)
      djxra=dxraj*xpsij*xraj1/(2.d0*xra)
      dj1xra=dxraj1*xpsij1*xraj/(2.d0*xra)
    endif

! calcul de la gravite. Intervient dans la frequence
! de Brunt-Vaisala avec effet de rotation inclu
    if (j1 /= m) then
      xpsij=exp(r(j))/((cst_G*(1.d0-exp(q(j)))*gms*Msol)/(omegi(j)*omegi(j)))**(1.d0/3.d0)
      xpsij1=exp(r(j1))/((cst_G*(1.d0-exp(q(j1)))*gms*Msol)/(omegi(j1)*omegi(j1)))**(1.d0/3.d0)
      if (xpsij >= rpsi_min .and. xpsij1 >= rpsi_min) then
        if (xpsij > rpsi_max) xpsij=0.9999999999d0*rpsi_max
! ancien appel: call geograv(xpsij,xgpsi)
        call geocalc(xpsij,xgpsi,1)
        xggj=cst_G*(1.d0-exp(q(j)))*gms*Msol*omegi(j)*omegi(j)*omegi(j)*omegi(j)
        xggj=xggj**(1.d0/3.d0)*xgpsi
        if (xpsij1 > rpsi_max) xpsij1=0.9999999999d0*rpsi_max
! ancien appel: call geograv(xpsij1,xgpsi1)
        call geocalc(xpsij1,xgpsi1,1)
        xggj1=cst_G*(1.d0-exp(q(j1)))*gms*Msol*omegi(j1)*omegi(j1)*omegi(j1)*omegi(j1)
        xggj1=xggj1**(1.d0/3.d0)*xgpsi1
        xggjm=sqrt(xggj*xggj1)
      else
        xggj=cst_G*(1.d0-exp(q(j)))*gms*Msol/exp(2.d0*r(j))
        xggj1=cst_G*(1.d0-exp(q(j1)))*gms*Msol/exp(2.d0*r(j1))
        xggjm=sqrt(xggj*xggj1)
      endif
    else
! on utilise ci-dessous: r=3/(4pi rhoc)^(1/3)*m_m^(1/3) et m_m=M_(m-1)/2
! Kippenhahn & Weigert, p. 68, Eq. (10.3)
      xggjm=cst_G*(4.d0*pi/3.d0)**(2.d0/3.d0)*2.d0**(-1.d0/3.d0)*exp(2.d0*rh1/3.d0)*((1.d0-exp(q(j)))*gms*Msol)**(1.d0/3.d0)
    endif
! calcul de l'echelle de pression
    if (j1 /= m) then
      xhpjm=exp(0.5d0*(p(j)+p(j1)))/(exp(0.5d0*(rh+rh1))*xggjm)
      xhpj1=exp(p(j1))/(exp(rh1)*xggj1)
    else
      xhpjm=exp(0.5d0*(p(j)+p(j1)))/(exp(0.5d0*(rh+rh1))*xggjm)
! au centre l'echelle de pression diverge, on prend xhpjm
      xhpj1=xhpjm
    endif
  else
    xoblaj = 1.0d0
  endif   !   irot

! G1: structure equation 3, dL/dM
  hnenn=1.d0/(q(j1)-q(j))
  d1=hnenn*(s(j1)-s(j))
  ff1=exp(ccg1-hh6+0.5d0*(q(j)+q(j1)-s(j)-s(j1)))
  f1=(en+eg-enue)*ff1
  g1=d1+f1
  g1s1=hnenn-0.5d0*f1
  g1s=-hnenn-0.5d0*f1
  hfak=en/2.d0
  hfakn=-enue*0.5d0
  g1p1=ff1*(hfak*epsp1+egp1+hfakn*enuep1)
  g1p=ff1*(hfak*epsp+egp+hfakn*enuep)
  g1t1=ff1*(hfak*epst1+egt1+hfakn*enuet1)
  g1t=ff1*(hfak*epst+egt+hfakn*enuet)
  tnorm=dzeit/1.0d+14
  g1  =tnorm*g1
  g1s =tnorm*g1s
  g1s1=tnorm*g1s1
  g1p =tnorm*g1p
  g1p1=tnorm*g1p1
  g1t =tnorm*g1t
  g1t1=tnorm*g1t1

! G2: structure equation 1, dP/dM
  d2=hnenn*(p(j1)-p(j))
  f2=-e(j)*exp(ccg2-2.d0*(r(j)+r(j1))+0.5d0*(q(j)+q(j1)-p(j)-p(j1)))

  if (irot==0 .or. omegi(j)<=1.d-20 .or. omegi(j1)<=1.d-20) then
    g2=d2+f2
    g2r1=-2.d0*f2
    g2r=g2r1
    hfak2=0.5d0*f2
    g2p1=hnenn-hfak2
    g2p=-hnenn-hfak2
  else
    g2=d2+f2*xfp
    g2r1=-2.d0*f2*xfp+f2*dj1xfp
    g2r=-2.d0*f2*xfp+f2*djxfp
    hfak2=0.5d0*f2*xfp
    g2p1=hnenn-hfak2
    g2p=-hnenn-hfak2
  endif

  d3=hnenn*(r(j1)-r(j))
  f3=exp(ccg3-1.5d0*(r(j)+r(j1))+0.5d0*(q(j)+q(j1)-rh1-rh))
  g3=d3+f3
  hfak=-0.5d0*f3
  g3p1=hfak*rhp1
  g3p=hfak*rhp
  g3t1=hfak*rht1
  g3t=hfak*rht
  g3r1=hnenn-1.5d0*f3
  g3r=-hnenn-1.5d0*f3

! critere de LEDOUX ou de SCHWARZSCHILD
! Si l'on veut remettre le critere de Hoiland-Solberg
! il faut enlever les commentaires du type H-S
!H-S      if (irot == 0) then
  if (iledou == 0) then
    xnbrun=adim-radm
    xbruj1=adi1-rad1
    xbruj =adi-rad
    adgrad(j)=adi-rad
  else
    admu=0.5d0*(adi +Nabla_mu(j)/(-rht) +adi1+Nabla_mu(j1)/(-rht1))
    admu1=adi1+Nabla_mu(j1)/(-rht1)
    xnbrun=admu-radm
    xbruj1=admu1-rad1
    xbruj=adi +Nabla_mu(j)/(-rht)-rad
    adgrad(j)=adi +Nabla_mu(j)/(-rht)-rad
  endif
!H-S      endif
!H-S      if (irot == 1.and.iledou == 0) then
!H-S      deltam=-0.5*(rht+rht1)
!H-S      xommoy=0.5*(omegi(j)+omegi(j1))
!H-S      dxommo=0.5*(dlodlr(j)+dlodlr(j1))
!H-S      xnbrun=xggjm*deltam/xhpjm*(adim-radm)
!H-S     &       +1.633*xommoy*xommoy*(2.+dxommo)
!H-S      xbruj1=xggj1*(-rht1)/xhpj1*(adi1-rad1)
!H-S     &       +1.633*omegi(j1)*omegi(j1)*(2.+dlodlr(j1))
!H-S      endif
!H-S      if (irot == 1.and.iledou == 1) then
!H-S      admu=adim+(xdmudp(j)+xdmudp(j1))/(-rht-rht1)
!H-S      admu1=adi1+xdmudp(j1)/(-rht1)
!H-S      deltam=-0.5*(rht+rht1)
!H-S      xommoy=0.5*(omegi(j)+omegi(j1))
!H-S      dxommo=0.5*(dlodlr(j)+dlodlr(j1))
!H-S      xnbrun=xggjm*deltam/xhpjm*(admu-radm)
!H-S     &       +1.633*xommoy*xommoy*(2.+dxommo)
!H-S      xbruj1=xggj1*(-rht1)/xhpj1*(admu1-rad1)
!H-S     &       +1.633*omegi(j1)*omegi(j1)*(2.+dlodlr(j1))
!H-S      endif
  wrm=0.5d0*(exp(r(j1))+exp(r(j)))
  if (xnbrun < 0.d0 .or. (iover == 1 .and. xover > wrm)) then
! Cas convectif
    if (irot==0 .or. omegi(j)<=1.d-20 .or. omegi(j1)<=1.d-20) then
      d4=hnenn*(t(j1)-t(j))
      f4=adim*f2
      g4=d4+f4
      g4r1=-2.d0*f4
      g4r=g4r1
      g4p1=hfak2*(adip1-adim)
      g4p=hfak2*(adip-adim)
      g4t1=hnenn+hfak2*adit1
      g4t=-hnenn+hfak2*adit
      g4s1=0.d0
      g4s=0.d0
    else
      d4=hnenn*(t(j1)-t(j))
      f4=adim*f2*xfp
      f44=1.d0
      g4=d4+f4
      g4r1=-2.d0*f4+adim*f2*dj1xfp
      g4r=-2.d0*f4+adim*f2*djxfp
      g4p1=hfak2*(adip1-adim)
      g4p=hfak2*(adip-adim)
      g4t1=hnenn+hfak2*adit1
      g4t=-hnenn+hfak2*adit
      g4s1=0.d0
      g4s=0.d0
    endif
  else
! Cas radiatif
    d4=hnenn*(t(j1)-t(j))
    f4=zradm*f2
    if (irot==0 .or. omegi(j)<=1.d-20 .or. omegi(j1)<=1.d-20) then
      g4=d4+f4
      g4r1=-2.d0*f4
      g4r=g4r1
      g4s1=hfak2*zrad1
      g4s=f4-g4s1
      g4p1=g4s1*(capp1+1.d0)-0.5d0*f4
      g4p=g4s*(capp+1.d0)-0.5d0*f4
      g4t1=hnenn+g4s1*(capt1-4.d0)
      g4t=-hnenn+g4s*(capt-4.d0)
      g4s1=g4s1/exphi(-s(j1))
      g4s=g4s/exphi(-s(j))
    else
      f44=xra*xfp
      f44rj=djxra*xfp+xra*djxfp
      f44rj1=dj1xra*xfp+xra*dj1xfp
      g4=d4+f4*f44
      g4r1=-2.d0*f4*f44+f4*f44rj1
      g4r=-2.d0*f4*f44+f4*f44rj
      g4s1=hfak2*zrad1*xra
      g4s=hfak2*zrad*xra
      g4p1=g4s1*(capp1+1.d0)-0.5d0*f4*f44
      g4p=g4s*(capp+1.d0)-0.5d0*f4*f44
      g4t1=hnenn+g4s1*(capt1-4.d0)
      g4t=-hnenn+g4s*(capt-4.d0)
      g4s1=g4s1/exphi(-s(j1))
      g4s=g4s/exphi(-s(j))
    endif

    if (xbruj > 0.d0) then
      zensi(j)=-zensi(j)
    endif
  endif

  if (irot==0 .or. omegi(j)<=1.d-20 .or. omegi(j1)<=1.d-20) then
    gradp(j)=-f2
    gradt(j)=-f4
  else
    gradp(j)=-f2*xfp
    gradt(j)=-f4*f44
  endif

  return

end subroutine gi
!***********************************************************************
subroutine zi

  use evol,only: kindreal
  use const,only: cst_G,pi,Msol
  use inputparam,only: irot,iledou
  use caramodele,only: gms,glm,hh6
  use abundmod,only: en,eg,enue,epsp1,enuep1,egp1,epsp,epst1,enuet1,epst,enuet,egt,enuep,egp,egt1
  use equadiffmod,only: ccg1,z1,z1p1,z1p,z1t1,z1t,ccz2,z2,z2p,z2p1,z2t1,z3,ccz3,z3p1,z3t1,z4,z4p1,z4p,z4t1,z4t,z4s
  use EOS,only: rh1,rh,rhp1,rht1,rht
  use strucmod,only: r,q,m,p,s,t,radm,xbruj1,adi1,rad1,adi,j,j1,adgrad,Nabla_mu,rad,adip1,adip,adit1,adit,capp,capt1,capt,zensi,&
                     capp1,adim
  use rotmod,only: omegi
  use geomod,only: rpsi_max,geom
  use SmallFunc,only: exphi

  implicit none
  real(kindreal):: xpsij,xfpj,xraj,dxfpj,xggjm,xggj1,xhpjm,dxraj,xhpj1,ff1,f1,fh1,fh,hfak,hfakn,d2,f2,z2rot,xnbrun=0.d0,admu,&
                   admu1,f4

  if (irot == 1) then
    xpsij=exp(r(m-1))/((cst_G*exp(ccg1)*(1.d0-exp(q(m-1))))/(omegi(m-1)*omegi(m-1)))**(1.d0/3.d0)
    if (xpsij > rpsi_max) xpsij=0.9999999999d0*rpsi_max
    call geom(xpsij,xfpj,xraj,dxfpj,dxraj)
! calcul de la gravite et de l'echelle de pression
! qui interviennent dans la frequence
! de Brunt-Vaisala avec effet de rotation inclu
! on utilise ci-dessous: r=3/(4pi rhoc)^(1/3)*m_m^(1/3) et m_m=M_(m-1)/2
! Kippenhahn & Weigert, p. 68, Eq. (10.3)
    xggjm=cst_G*(4.d0*pi/3.d0)**(2.d0/3.d0)*2.d0**(-1.d0/3.d0)*exp(2.d0*rh1/3.d0)*((1.d0-exp(q(j)))*gms*Msol)**(1.d0/3.d0)
    xggj1=xggjm
    xhpjm=exp(0.5d0*(p(j)+p(j1)))/(exp(0.5d0*(rh+rh1))*xggjm)
    xhpj1=xhpjm
  endif

  ff1=exphi(q(m-1))
  f1=log(ff1)
  fh1=exp(glm-hh6)*ff1
  fh=(en+eg-enue)*fh1
  z1=s(m-1)-log(1.d0+fh)
  if (isnan(z1)) then
  write (*,*)"hh6,exp(glm-hh6),ff1,enue,en+eg-enue,s(m-1)" ,hh6,exp(glm-hh6),ff1,enue,en+eg-enue,s(m-1)
  stop "z1=NaN"
  endif
  fh=fh1/(1.d0+fh)
  hfak=en*0.5d0
  hfakn=-enue*0.5d0
  z1p1=-fh*(hfak*epsp1+hfakn*enuep1+egp1)
  z1p=-fh*(hfak*epsp+hfakn*enuep+egp)
  z1t1=-fh*(hfak*epst1+hfakn*enuet1+egt1)
  z1t=-fh*(hfak*epst+hfakn*enuet+egt)
  d2=exp (p(j)-p(j1))
! f2: calcul de Pc selon Eq. (10.6) p. 69 de Kippenhahn & Weigert
  f2=exp (ccz2+(4.d0*rh1+2.d0*f1)/3.d0-p(j1))
  if (irot == 1 .or.omegi(j) > 1.d-20.or.omegi(j1) > 1.d-20) then
! Annexe K de "Rotation & Stellar Structure" private comm. Meynet
    z2rot=2.d0/3.d0*(0.6d0*xfpj+0.9d0)
    f2=f2*z2rot
  endif
  z2=d2+f2-1.d0
  z2p=d2
  z2p1=-d2+f2*(4.d0*rhp1/3.d0-1.d0)
  z2t1=4.d0*f2*rht1/3.d0
  z3=ccz3+f1-rh1-3.d0*r(j)
  z3p1=-rhp1
  z3t1=-rht1

!H-S      if (irot == 0.and.iledou == 0) then
  if (iledou == 0) then
    xnbrun=adim-radm
    xbruj1=adi1-rad1
    adgrad(j)=adi-rad
    adgrad(j1)=adi1-rad1
  endif
!H-S      if (irot == 0.and.iledou == 1) then
  if (iledou == 1) then
    admu=0.5d0*(adi +Nabla_mu(j)/(-rht)+adi1+Nabla_mu(j1)/(-rht1))
    admu1=adi1+Nabla_mu(j1)/(-rht1)
    xnbrun=admu-radm
    xbruj1=admu1-rad1
    adgrad(j)=adi +Nabla_mu(j)/(-rht)-rad
    adgrad(j1)=admu1-rad1
  endif

  if (xnbrun < 0.d0) then
    hfak=-0.5d0*(p(j)-p(j1))
    z4=t(j)-t(j1)+2.d0*adim*hfak
    z4p1=hfak*adip1+adim
    z4p=hfak*adip-adim
    z4t1=-1.d0+hfak*adit1
    z4t=1.d0+hfak*adit
    z4s=0.d0
  else
    f4=-radm*(p(j)-p(j1))
    z4=t(j)-t(j1)+f4
    z4p1=f4*(0.5d0*capp1+0.5d0)+radm
    z4p =f4*(0.5d0*capp+0.5d0)-radm
    z4t1=-1.d0+f4*(0.5d0*capt1-2.d0)
    z4t =1.d0+f4*(0.5d0*capt-2.d0)
    z4s = f4/exphi(-s(m-1))
    if (xbruj1 > 0.d0) then
      zensi(j)=-zensi(j)
    endif
  endif

  return

end subroutine zi
!***********************************************************************
subroutine henyey
! version avec ROTATION
!-----------------------------------------------------------------------
! Resolution des equations aux derivees partielles du premier ordre
!   pour l'interieur (0=< Mr =< MrF) par la methode de relaxation
!   de Henyey.

! Interieur stellaire divise en (m-2) couches.

! Systeme de (4m-2) equations a (4m-2) inconnues :
!         Bk(x1,s1,p1,teta1)=0 (k=1,2)
!           : Conditions limites pres de la surface.
!         Gi(xj,...,tetaj+1)=0 (i=1,2,3,4; j=1,...,m-2)
!           : Equations aux differences finies.
!         Zi(xm-1,...,sm)=0    (i=1,2,3,4)
!           : Conditions limites au centre.
! Les equations Gi et Zi sont non lineaires,
!   les fonctions Bk sont lineaires.

! Pour le modele initial approche :
!         Bk0#0, Gi0#0, Zi0#0.

! Determination des dBk, dGi et dZi tels que :
!         Bk0+dBk=0,
!         Gi0+dGi=0,
!         Zi0+dZi=0.

! Systeme non lineaire ---> Plusieurs iterations.

! La matrice A(k,l) est le BLOC1. k=1,2 ---> Bk
!                                 k=3,4,5,6 ---> Gk
! Pour l'interieur stellaire, on admet la quasi-adiabaticite
!   de la convection (et non plus la super-adiabaticite).

! HENYEY est appellee par MAIN,
!         et appelle UNDERS, OVER1, GI, GIRL, ZI, NETWKI,
!                    DICHTE, KAPPA, NABLA, ENERG.
! HENYEY appelle NUCAL et NUINT pour le calcul des neutrinos solaires.
! IPRC  =1: on imprime le modele
!       =0: on n'imprime pas le modele

! Derniere version : 2 decembre 2009
!-----------------------------------------------------------------------
  use evol, only: ldi
  use inputparam, only: modanf,alph,iout,imagn,isol,istati,iledou,idiff,idifcon,iover,iunder,gkorm,phase, &
    agdr,agds,agdp,agdt,ichem,idebug,plot,refresh,idebug,Add_Flux,nwseq
  use caramodele, only: rhoc,tc,hh6,PrintError,teff,gls
  use abundmod,only: epsn1,enuet,enuet1,enuep,enuep1,epsp,epsp1,epst,epst1
  use equadiffmod,only: gkor,iter,iprc,g1,g2,g3,g4,g1s,g1p,g1t,g1s1,g1p1,g1t1,g2r,g2p,g2r1,g2p1,g3r,g3p,g3t,g3r1,g3p1,g3t1,g4r, &
    g4s,g4p,g4t,g4r1,g4s1,g4p,g4p1,g4t1,z1,z2,z3,z4,z1p,z1t,z1p1,z1t1,z2p,z2p1,z2t1,z3p1,z3t1,z4p,z4s,z4t,z4p1,z4t1
  use EOS,only: dichte,rh,rh1,rhp,rhp1,rht,rht1,psi,num,invert_helm_pt,toni,rhe,gamma1_Timmes,gamma1_dichte
  use strucmod,only: m,j,j1,beta,beta1,vmy1,cap,cap1,capp,capp1,capt,capt1,rad,rad1,zrad,zrad1,adi,adi1,adip,adip1,adit,adit1, &
    xnabj,xnabj1,t,zensi,adgrad,xbruj1,Nabla_rad,Nabla_ad,delt,bet,opac,opact,epsit,rho,r,p,s,q,vr,vp,vt,rrp,rrt, &
    rrc,rlp,rlt,rlc,x_env
  use magmod,only: D_magx
  use omegamod,only: omenew,dlonew,omconv,omesta,vomcon
  use rotmod,only: dlelexsave,BTotal_EndAdvect,btotal_startmodel,Flux_remaining,vsuminenv,vvsuminenv
  use convection,only: over1,unders
  use diffadvmod,only: xnabyy,D_conv,D_shear,D_eff
  use PGPlotModule,only: Struc_Plotted,PlotStruc
  use SmallFunc,only: exphi,girl
  use advection,only: advect
  use opacity,only: kappa,ioutable,rout,tout
  use energy,only: energ,vmassen,rvect,t9n,pvect,epstot1,epsneut,dcoeff
  use chemicals,only: netnew,chemeps,chemold
  use diffusion,only: coedif,diffbr,diffom
  use nablas,only: nabla,nabgam,grapmui
  use timestep,only: alter,dzeit
  use PrintAll,only: Teff_save,Lum_save,mass_save,time_save,C12_save,C13_save,N14_save,O16_save
  use ionisation,only: ionized
! for ifort compiler, uncomment the next line:
!  use, INTRINSIC:: IEEE_ARITHMETIC, only: isnan => IEEE_IS_NAN

  implicit none

  integer:: ic,ii,jgg1,jgg2,jgg3,jgg4,j1v,jv,i,jgdr,jgds,jgdp,jgdt,iterlim1,iterlim2,flag_girl,iSE,jSE
  real(kindreal):: fred,vgdt,alph1,vmy,vrhoc,xm,egc,drhoc,zwi1,gg1,gg2,gg3,gg4,dp,dt,dp1,dt1,dr,ds,gdr,gds,gdp,gdt,t6
  real(kindreal), dimension(ldi):: ar,as,ap,at,br,bs,bp,bt,ccr,ccs,ccp,cct
  real(kindreal), dimension(6,9):: a
  real(kindreal), dimension(6,3):: u_hen
  real(kindreal), dimension(4,7),save:: ha
  real(kindreal), dimension(4,3):: hu

  logical::endIter
!-----------------------------------------------------------------------
  fred=1.d0
  vgdt=20.d0
  endIter=.false.

  alph1=0.125d0*alph
! alph : facteur de relaxation
!      : 1/4, a la premiere iteration
!      : 1/2, a la seconde iteration
!      : 1, a la troisieme iteration
! xj1 = xj0 + alph*delta(xj)
! ...
! tetaj1 = tetaj0 + alph*delta(tetaj)

  write(3,'(////,a,//)') ' henyey-methode'

  gkor = 0.d0
  iter=0
  ic=0

  do while (.not. endIter)
   iter=iter+1
   PrintError = .true.
! iter est le compteur d'iterations effectuees.
! iiter est le compteur d'iterations effectuees utilise pour la diffusion

! Les 3 lignes suivantes remettent a zero apres chaque iteration
! le compteur de depassement des limites des tables d'opacite
   ioutable = 0
   rout = 0.d0
   tout = 10.d0

   if (iter < 5) then
     alph1=2.d0*alph1
     alph1=min(alph1,alph)
   endif

   if (iter >= 4)alph1=fred*alph

! TRAITEMENT COUCHE PAR COUCHE
!--------------------------------------------------------------------
! Calcul du bloc 1. Couche 1 de l'interieur stellaire.
!------------------------------------------------------
   if (idebug > 0) then
     write(*,*) 'Henyey: layers treatment - going down'
   endif
   j=1
   j1=1
   if (iter > 1 .and. iunder /= 0) then
     if (idebug > 1) then
       write(*,*) 'call unders'
     endif
     call unders
   endif

   do
    if (idebug > 1) then
      write(*,*) 'Layer: ',j1
    endif
    if (j1 >= 2) then
      beta = beta1
      vmy = vmy1
      rh = rh1
      rhp = rhp1
      rht = rht1
      cap = cap1
      capp = capp1
      capt = capt1
      epsn=epsn1
      enuet=enuet1
      enuep=enuep1
      epsp = epsp1
      epst = epst1
      rad = rad1
      zrad = zrad1
      adi = adi1
      adip = adip1
      adit = adit1
! nabla calcule avec le gamma
     xnabj=xnabj1
    endif

    if (idebug > 1) then
      write(*,*) 'call dichte'
    endif

!!!!!!!!!!!!!!!!!  SWITCH EOS TIMMES/DICHTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!


if (EOS == 0) then
    call dichte

endif


if (EOS == 1) then
    ! On utilise dichte la plupart du temps, mais on switch
    ! sur Timmes EOS lorsque l'on atteint des régimes de hautes
    ! températures et/ou densités.
    if (j == 1) THEN

      call DICHTE
      continue

    else if ( (exp(rh1) .lt. 10**2.8d0) .or. (exp(t(j1)) .lt. 10**7.55d0) ) then !! Domaine du switch utilisé par MESA (Paper I 2011) !!

    call DICHTE

    ELSE

    call invert_helm_pt



    ENDIF
endif

open (177,file='EOS_check.dat',status='unknown',form='formatted')
if (henyey_last) then
if ( (exp(rh1) .lt. 10**2.8d0) .or. (exp(t(j1)) .lt. 10**7.55d0) ) then
  write(177,*),nwseq,j1,0,gamma1_dichte
ELSE
  write(177,*),nwseq,j1,1,gamma1_Timmes
endif
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calcul des opacites
    if (idebug > 1) then
      write(*,*) "Kappa model for advanced phases recent developements, if bugs contact adam.griffiths@uv.es"
      write(*,*) 'call kappa'
    endif

  if ( ( x(j1) .ne. 0.d0 )  .and. ( t(j1) .ge. log(3e8) ) )  then !Advanced phase opacity do not include protons in opacity computation.
      call kappa(rh1,t(j1),rhp1,rht1,0.d0,y(j1),cap1,capp1,capt1,j1)
  else
    call kappa(rh1,t(j1),rhp1,rht1,x(j1),y(j1),cap1,capp1,capt1,j1)
  endif

! Calcul des gradients adiabatiques et radiatifs et de leurs derivees
    if (idebug > 1) then
      write(*,*) 'call nabla'
    endif
    call nabla

! Calcul du gradient thermique du milieu
    if (irot == 1) then
      if (idebug > 1) then
        write(*,*) 'call nabgam'
      endif
      call nabgam
    endif


! Calcul des energies
    if (idebug > 1) then
      write(*,*) "Advance phase network under recent developements, if bugs contact adam.griffiths@uv.es"
      write(*,*) 'call energ'
    endif
    call energ

    adgrad(j1)=xbruj1
    Nabla_rad(j1)=rad1
    Nabla_ad(j1)=adi1
    xnabyy(j1)=xnabj1
    delt(j1)=-rht1
    bet(j1)=beta1
    opac(j1)=exp(cap1)
    opact(j1)=capt1
    epsit(j1)=epst1
    rho(j1)=rh1
    if (itminc == 1) then
      iout=0
      if (j > 1 .and. j <= m-3 .and. zensi(j1) > 0.d0 .and. zensi(j) < 0.d0) iout=j
      if (j1 >= m) then
        vrhoc=0.0d0
        if (rhoc > 0.d0) vrhoc=rhoc
        rhoc=rh1/um
        tc=t(j1)/um
        xm=x(j1)
        egc=eg
        if (vrhoc > 0.d0) drhoc=rhoc-vrhoc
      endif
    endif

!-----------------------------------------------------------------
    if (j1 == 1) then
!-----------------------------------------------------------------
! Initialisation et lien avec les conditions en bas de l'enveloppe
!-----------------------------------------------------------------
      j1 = 2
      cycle

!-------------------------------------------------------------
    else if (j1 == 2) then
!-------------------------------------------------------------
! Calcul du bloc 2 (n=2). Couches 2 + conditions au bord
!-------------------------------------------------------------
      a(1,1)=1.d0      ! dB1/dxj
      a(1,2)=0.d0      ! dB1/dsj
      a(1,3)=-rrp      ! dB1/dpj
      a(1,4)=-rrt      ! dB1/dtetaj
      a(1,5)=0.d0      ! dB1/dxj+1
      a(1,6)=0.d0      ! dB1/dsj+1
      a(1,7)=0.d0      ! dB1/dpj+1
      a(1,8)=0.d0      ! dB1/dtetaj+1
      a(1,9)=-r(1)+rrt*t(1)+rrp*p(1)+rrc      ! B1
      a(2,1)=0.d0      ! dB2/dxj
      a(2,2)=1.d0/(1.d0 - exp(-s(1)))      ! dB2/dsj
      a(2,3)=-rlp      ! dB2/dpj
      a(2,4)=-rlt      ! dB2/dtetaj
      a(2,5)=0.d0      ! dB2/dxj+1
      a(2,6)=0.d0      ! dB2/dsj+1
      a(2,7)=0.d0      ! dB2/dpj+1
      a(2,8)=0.d0      ! dB2/dtetaj+1
      a(2,9)=rlt*t(1)+rlp*p(1)+rlc- log(exp(s(1))-1.d0)-hh6   ! B2
      if (isnan(a(2,9))) then
        if (exp(s(1))-1.d0 <= 0.d0) then
          write(*,*) 'a(2,9)=NaN, log of negative number: exp(s(1))-1=',exp(s(1))-1.d0
        else
          write(*,*)'a(2,9)=NaN - rlt,t(1),rlp,p(1),rlc,exp(s(1))-1,hh6:',rlt,t(1),rlp,p(1),rlc,exp(s(1))-1.d0,hh6
        endif
      endif

! Calcul de Gi, dGi/dxj, dGi/dsj, dGi/dtetaj, dGi/dpj.
!------------------------------------------------------
! i=1,2,3,4, j=1,2.

      if (idebug > 1) then
        write(*,*) 'call gi'
      endif
      call gi

! [mod xfile]
      vmassen(j)=gms*(1.d0-exp(q(j)))
      rvect(j)= 10.d0**(r(j)/um)
      t9n(j)=10.d0**(t(j)/um)
      pvect(j)=10.d0**(p(j)/um)
      epstot1(j)=eps(j)+epsy(j)+epsc(j)
      epsneut(j)=-epsn

      if (imagn == 1) then
        dcoeff(j)=D_conv(j)+D_shear(j)+D_eff(j)+D_magx(j)
      else
        dcoeff(j)= D_conv(j)+D_shear(j)+D_eff(j)
      endif
! [/mod xfile]

      if (iprc  >  0) then

        write(3,'(/,"  j",4x,"vm",6x,"p",7x,"t",8x,"r",7x,"vl",6x,"x",5x,"y",6x,"xc12",6x,"xo16",8x,"eps",8x,"epsy",7x,&
                     &"epsc",7x,"radm",4x,"rh",5x,"zensi",/6x,"epsn",4x,"capp",4x,"capt",5x,"epsp",4x,"epst",4x,"alpha",&
                     &1x,"delta",3x,"psi",5x,"epsyy",5x,"epsyc",6x,"epsyo",6x,"eg",9x,"adim",4x,"cap",5x,"beta",/7x,"y3",&
                     &8x,"xc13",4x,"xn14",4x,"xn15",10x,"xo17",4x,"xo18",4x,"xne20",7x,"xne22",8x,"xmg24",8x,"xmg25",8x,&
                     &"xmg26",11x,"omega",/2x,"grad mu",3x,"gamma",2x,"Ri",3x,"Nabla",3x,"Dconv",3x,"Dshear",3x,"Deff"//)')

        if (ialflu == 1) then
          write(3,'(7x,"F19",7x,"NE21",4x,"NA23",4x,"AL26",10x,"AL27",4x,"SI28",4x,"C14",7x,"F18",8x,"NEU",8x,&
                    &"PRO",8x,"BID",//)')
        endif

        write(3,'("    i) Z(i) A(i):   ",77(i4,")",2i4)//)') (ii,nbzel(ii),nbael(ii),ii=1,nbelx)
        write(3,'(//)')

        zwi1=1.d0/(exp(s(1))-1.d0)

        call printhenyey(rh/um,cap/um,capp,capt,epsp,epst,rhp,-rht,beta,zwi1)
! save main data for full printing
        Teff_save = teff
        Lum_save = gls
        mass_save = gms
        time_save = alter
      else
        call Calcvmyhelio
      endif

      gg1=g1
      gg2=g2
      gg3=g3
      gg4=g4
      jgg1=1
      jgg2=1
      jgg3=1
      jgg4=1

      if (idebug == 2) then
        if (isnan(g1).or.isnan(g2).or.isnan(g3).or.isnan(g4)) then
          write(*,*)'iter,j,g1,g2,g3,g4',iter,j,g1,g2,g3,g4
          stop
        endif
      endif

! Equation avec G1
      a(3,1)=0.d0
      a(3,2)=g1s
      a(3,3)=g1p
      a(3,4)=g1t
      a(3,5)=0.d0
      a(3,6)=g1s1
      a(3,7)=-g1p1
      a(3,8)=-g1t1
      a(3,9)=-g1

! Equation avec G2
      a(4,1)=g2r
      a(4,2)=0.d0
      a(4,3)=g2p
      a(4,4)=0.d0
      a(4,5)=g2r1
      a(4,6)=0.d0
      a(4,7)=-g2p1
      a(4,8)=0.d0
      a(4,9)=-g2

! Equation avec G3
      a(5,1)=g3r
      a(5,2)=0.d0
      a(5,3)=g3p
      a(5,4)=g3t
      a(5,5)=g3r1
      a(5,6)=0.d0
      a(5,7)=-g3p1
      a(5,8)=-g3t1
      a(5,9)=-g3

! Equation avec G4
      a(6,1)=g4r
      a(6,2)=g4s
      a(6,3)=g4p
      a(6,4)=g4t
      a(6,5)=g4r1
      a(6,6)=g4s1
      a(6,7)=-g4p1
      a(6,8)=-g4t1
      a(6,9)=-g4

      if (idebug == 2) then
        do iSE=1,6
         do jSE=1,9
          write(3,'(a,2i3,a,i1,a,i1,a,d22.12)')'iter,j:',iter,j,', a(',iSE,',',jSE,') : ',a(iSE,jSE)
          if (isnan(a(iSE,jSE))) then
            write(*,'(a,2i3,a,i1,a,i1,a,d22.12)')'iter,j:',iter,j,', a(',iSE,',',jSE,') : ',a(iSE,jSE)
            stop
          endif
         enddo
        enddo
      endif

! La matrice A (bloc 1) est formee de 6 equations a 8 inconnues
! (dxi, dsi, dpi, dtetai, i=1,2).
! On exprime les variables par combinaisons lineaires de dp2
! et dteta2 :
!             dx1=U1*dp2+V1*dteta2+W1
!             ...
!             ds2=U6*dp2+V6*dteta2+W6
! Le tableau u(m,n) contient ces coefficients :
!             u(1,1)=U1,...,u(6,3)=W6
! girl calcule u(m,n), soit (Ui,Vi,Wi)(i=1,...,6) par inversion
! de matrice.

      if (idebug > 3) then
        write(*,*) 'call girl(a,u_hen)'
      endif
      call girl(a,u_hen,6,3,flag_girl)
      if (flag_girl /= 0) then
        if (idebug>0) then
          write(*,*) 'henyey - matrix a(6,9),flag:',flag_girl
          do iSE=1,6
           do jSE=1,9
            write(*,'("a(",i1,",",i1,") :",d22.12)') iSE,jSE,a(iSE,jSE)
          enddo
         enddo
        endif
        rewind(222)
        write(222,*) nwmd,':girl crash in henyey with matrix a(6,3)'
        stop
      endif


      if (idebug == 2) then
        do iSE=1,6
         do jSE=1,3
          write(3,'(a,2i3,a,i1,a,i1,a,d22.12)')'iter,j:',iter,j,', u(',iSE,',',jSE,') : ',u_hen(iSE,jSE)
          if (isnan(u_hen(iSE,jSE))) then
            write(*,'(a,2i3,a,i1,a,i1,a,d22.12)')'iter,j:',iter,j,', u(',iSE,',',jSE,') : ',u_hen(iSE,jSE)
            stop
          endif
         enddo
        enddo
      endif

      ar(1)=u_hen(1,1)
      br(1)=u_hen(1,2)
      ccr(1)=u_hen(1,3)
      as(1)=u_hen(2,1)
      bs(1)=u_hen(2,2)
      ccs(1)=u_hen(2,3)
      ap(1)=u_hen(3,1)
      bp(1)=u_hen(3,2)
      ccp(1)=u_hen(3,3)
      at(1)=u_hen(4,1)
      bt(1)=u_hen(4,2)
      cct(1)=u_hen(4,3)
      ar(2)=u_hen(5,1)
      br(2)=u_hen(5,2)
      ccr(2)=u_hen(5,3)
      as(2)=u_hen(6,1)
      bs(2)=u_hen(6,2)
      ccs(2)=u_hen(6,3)

      j=j+1
      j1=j+1

      cycle

!--------------------------------------------------------------
    else if (j1 < m) then
!--------------------------------------------------------------
! Calcul du bloc n (2<=n<=m-2). Couches 2 a m-2 de l'interieur.
!--------------------------------------------------------------
! Les variables sont exprimees sous forme de combinaisons lineaires
! de dpn et dtetan :
!                    dxn=U4n-3*dpn+V4n-3*dtetan+W4n-3
!                    dsn=U4n-2*dpn+V4n-2*dtetan+W4n-2
      if (idebug > 1) then
        write(*,*) 'call gi'
      endif
      call gi

      if (idebug == 2) then
        if (isnan(g1).or.isnan(g2).or.isnan(g3).or.isnan(g4)) then
          write(*,*)'iter,j,g1,g2,g3,g4',iter,j,g1,g2,g3,g4
          stop
        endif
      endif

! [mod xfile]
      vmassen(j)=gms*(1.d0-exp(q(j)))
      rvect(j)= 10.d0**(r(j)/um)
      t9n(j)=10.d0**(t(j)/um)
      pvect(j)=10.d0**(p(j)/um)
      epstot1(j)=eps(j)+epsy(j)+epsc(j)
      epsneut(j)=-epsn

      if (imagn == 1) then
        dcoeff(j)=D_conv(j)+D_shear(j)+D_eff(j)+D_magx(j)
      else
        dcoeff(j)= D_conv(j)+D_shear(j)+D_eff(j)
      endif
! [/mod xfile]

      if (iprc > 0) then
        call printhenyey(rh/um,cap/um,capp,capt,epsp,epst,rhp,-rht,beta,zwi1)
      else
        call Calcvmyhelio
      endif

      if (abs(gg1)-abs(g1)  <  0.d0) then
        gg1=g1
        jgg1=j
      endif
      if (abs(gg2)-abs(g2)  <  0.d0) then
        gg2=g2
        jgg2=j
      endif
      if (abs(gg3)-abs(g3)  <  0.d0) then
        gg3=g3
        jgg3=j
      endif
      if (abs(gg4)-abs(g4)  <  0.d0) then
        gg4=g4
        jgg4=j
      endif

! Remplissage du bloc n = HA(k,l) (n=2,...,m-2).
!                     k=4 : G1, G2, G3, G4.
!                     l=7
      ha(1,1)=g1p+as(j)*g1s
!            : dG1/dpj +U6*dg1/dsj
      ha(1,2)=g1t+bs(j)*g1s
      ha(1,3)=0.d0
      ha(1,4)=g1s1
      ha(1,5)=-g1p1
      ha(1,6)=-g1t1
      ha(1,7)=  -ccs(j)*g1s-g1
      ha(2,1)=g2p+ar(j)*g2r
      ha(2,2)=br(j)*g2r
      ha(2,3)=g2r1
      ha(2,4)=0.d0
      ha(2,5)=-g2p1
      ha(2,6)=0.d0
      ha(2,7)=  -ccr(j)*g2r-g2
!         : gama 2
      ha(3,1)=g3p+ar(j)*g3r
      ha(3,2)=g3t+br(j)*g3r
      ha(3,3)=g3r1
      ha(3,4)=0.d0
      ha(3,5)=-g3p1
      ha(3,6)=-g3t1
      ha(3,7)=  -ccr(j)*g3r-g3
!         : gama 3
      ha(4,1)=g4p+ar(j)*g4r+as(j)*g4s
      ha(4,2)=g4t+br(j)*g4r+bs(j)*g4s
      ha(4,3)=g4r1
      ha(4,4)=g4s1
      ha(4,5)=-g4p1
      ha(4,6)=-g4t1
      ha(4,7)=  -ccr(j)*g4r-ccs(j)*g4s-g4
!         : gama 4

      if (idebug == 2) then
        do iSE=1,4
         do jSE=1,7
          write(3,'(a,2(1x,i3),a,i1,a,i1,a,d22.12)')'iter,j:',iter,j,', ha(',iSE,',',jSE,') : ',ha(iSE,jSE)
          if (isnan(ha(iSE,jSE)) .or. abs(ha(iSE,jSE))>HUGE(ha(iSE,jSE))) then
            write(*,'(a,2(1x,i3),a,i1,a,i1,a,d22.12)')'iter,j:',iter,j,', ha(',iSE,',',jSE,') : ',ha(iSE,jSE)
            write(*,*) 'g1p,as(j),g1s:',g1p,as(j),g1s
            stop
          endif
         enddo
        enddo
      endif

! Calcul de (Ui,Vi,Wi) (i=4n-2,4n,4n+2)
!---------------------------------------
! Le tableau hu(m,n) contient ces coefficients :
!            hu(1,1) = U4n-1
!            ...
!            hu(4,3) = W4n+2

      if (idebug > 3) then
        write(*,*) 'call girl(ha,hu)'
      endif
      call girl(ha,hu,4,3,flag_girl)
      if (flag_girl /= 0) then
        if (idebug>0) then
          write(*,*) 'henyey - matrix ha(4,7),flag:',flag_girl
          do iSE=1,4
           do jSE=1,7
            write(*,'("ha(",i1,",",i1,") :",d22.12)') iSE,jSE,ha(iSE,jSE)
          enddo
         enddo
        endif
        rewind(222)
        write(222,*) nwmd,':girl crash in henyey with matrix ha(4,7)'
        stop
      endif

      if (idebug == 2) then
        do iSE=1,4
         do jSE=1,3
          write(3,'(a,2(1x,i3),a,i1,a,i1,a,d22.12)')'iter,j:',iter,j,', hu(',iSE,',',jSE,') : ',hu(iSE,jSE)
          if (isnan(hu(iSE,jSE))) then
            write(*,'(a,2(1x,i3),a,i1,a,i1,a,d22.12)')'iter,j:',iter,j,', hu(',iSE,',',jSE,') : ',hu(iSE,jSE)
            stop
          endif
         enddo
        enddo
      endif

      ap(j)=hu(1,1)
      bp(j)=hu(1,2)
      ccp(j)=hu(1,3)
      at(j)=hu(2,1)
      bt(j)=hu(2,2)
      cct(j)=hu(2,3)
      ar(j1)=hu(3,1)
      br(j1)=hu(3,2)
      ccr(j1)=hu(3,3)
      as(j1)=hu(4,1)
      bs(j1)=hu(4,2)
      ccs(j1)=hu(4,3)

      if (j < m-2) then
        j=j+1
        j1=j+1
      else
        j=j+1
        j1=m
      endif

      cycle
!-----------------------------------------------------
    else   ! CENTRE
!-----------------------------------------------------
! Calcul du bloc m (centre). Couches m de l'interieur.
!-----------------------------------------------------
      if (idebug > 0) then
        write(*,*) 'Henyey: layers treatment - centre'
      endif

      if (idebug > 1) then
        write(*,*) 'call zi'
      endif
      call zi

      if (idebug == 2) then
        if (isnan(z1).or.isnan(z2).or.isnan(z3).or.isnan(z4) &
            .or. abs(z1)>HUGE(z1).or. abs(z2)>HUGE(z2).or. &
             abs(z3)>HUGE(z3).or. abs(z4)>HUGE(z4)) then
          write(*,*)'iter,j,z1,z2,z3,z4',iter,j,z1,z2,z3,z4
          stop
        endif
      endif

! [mod xfile]
      vmassen(j)=gms*(1.d0-exp(q(j)))
      rvect(j)= 10.d0**(r(j)/um)
      t9n(j)=10.d0**(t(j)/um)
      pvect(j)=10.d0**(p(j)/um)
      epstot1(j)=eps(j)+epsy(j)+epsc(j)
      epsneut(j)=-epsn
      vmassen(j1)=gms*(1.d0-exp(q(j1)))
      rvect(j1)= 10.d0**(r(j1)/um)
      t9n(j1)=10.d0**(t(j1)/um)
      pvect(j1)=10.d0**(p(j1)/um)
      epstot1(j1)=eps(j1)+epsy(j1)+epsc(j1)
      epsneut(j1)=-epsn

      if (imagn == 1) then
        dcoeff(j)=D_conv(j)+D_shear(j)+D_eff(j)+D_magx(j)
        dcoeff(j1) = D_conv(j1)+D_shear(j1)+D_eff(j1)+D_magx(j1)
      else
        dcoeff(j)= D_conv(j)+D_shear(j)+D_eff(j)
        dcoeff(j1) = D_conv(j1)+D_shear(j1)+D_eff(j1)
      endif
! [/mod xfile]

      if (iprc > 0) then
        call printhenyey(rh/um,cap/um,capp,capt,epsp,epst,rhp,-rht,beta,zwi1)
        j1v=j1
        jv=j
        j=j1
        write(3,*)'j1v,jv,j,j1:',j1v,jv,j,j1
        call printhenyey(rh1/um,cap1/um,capp1,capt1,epsp1,epst1,rhp1,-rht1,beta1,zwi1)
        write (3,'(" at centre   psi=",f7.1,"  vmy=",f7.3//" COMPUTED WITH IONIZATION UNTIL J1=",i5///)') psi,vmy1/um,num-1
        j1=j1v
        j = jv
      else
        call Calcvmyhelio
        j1v=j1
        jv=j
        j=j1
        call Calcvmyhelio
        j1=j1v
        j = jv
      endif

      if (abs(gg1)-abs(z1)  <  0.d0) then
        gg1=z1
        jgg1=j
      endif
      if (abs(gg2)-abs(z2)  <  0.d0) then
        gg2=z2
        jgg2=j
      endif
      if (abs(gg3)-abs(z3)  <  0.d0) then
        gg3=z3
        jgg3=j
      endif
      if (abs(gg4)-abs(z4)  <  0.d0) then
        gg4=z4
        jgg4=j
      endif
      write(3,'(23x,"biggest gi",14x,4(i6,e12.5))') jgg1,gg1,jgg2,gg2,jgg3,gg3,jgg4,gg4

      ha(1,1)=as(j)+z1p
!          : U4m-6 + dZ1/dpm-1
      ha(1,2)=bs(j)+z1t
      ha(1,3)=z1p1
      ha(1,4)=z1t1
      ha(1,5)=  -ccs(j)-z1
      ha(2,1)=z2p
      ha(2,2)=0.d0
      ha(2,3)=z2p1
      ha(2,4)=z2t1
      ha(2,5)=-z2
      ha(3,1)=-ar(j)*3.d0
      ha(3,2)=-br(j)*3.d0
      ha(3,3)=z3p1
      ha(3,4)=z3t1
      ha(3,5)=  +ccr(j)*3.d0-z3
      ha(4,1)=z4p+as(j)*z4s
      ha(4,2)=z4t+bs(j)*z4s
      ha(4,3)=z4p1
      ha(4,4)=z4t1
      ha(4,5)=-z4-ccs(j)*z4s
      do i=1,4
       ha(i,6)=0.d0
       ha(i,7)=0.d0
      enddo

      if (idebug == 2) then
        do iSE=1,4
         do jSE=1,7
          write(3,'(a,2(1x,i3),a,i1,a,i1,a,d22.12)')'iter,j:',iter,j,', ha(',iSE,',',jSE,') : ',ha(iSE,jSE)
          if (isnan(ha(iSE,jSE))) then
            write(*,'(a,2(1x,i3),a,i1,a,i1,a,d22.12)')'iter,j:',iter,j,', ha(',iSE,',',jSE,') : ',ha(iSE,jSE)
            stop
          endif
         enddo
        enddo
      endif

! HU est une matrice (4,1).
! On ajoute deux colonnes contenant des 0.
! La premiere colonne donne les corrections :
!    dpm-1, dtetam-1, dpm, dtetam.

! girl calcule hu(m,n) par inversion de matrice.
      if (idebug > 3) then
        write(*,*) 'call girl(ha,hu)'
      endif
      call girl(ha,hu,4,3,flag_girl)
      if (flag_girl /= 0) then
        if (idebug>0) then
          write(*,*) 'henyey - matrix ha(4,7),flag:',flag_girl
          do iSE=1,4
           do jSE=1,7
            write(*,'("ha(",i1,",",i1,") :",d22.12)') iSE,jSE,ha(iSE,jSE)
          enddo
         enddo
        endif
        rewind(222)
        write(222,*) nwmd,':girl crash in henyey with matrix ha(4,3)'
        stop
      endif

      if (idebug == 2) then
        do iSE=1,4
         write(3,'(a,2(1x,i3),a,i1,a,d22.12)')'iter,j:',iter,j,', hu(',iSE,',1) : ',hu(iSE,1)
         if (isnan(hu(iSE,1))) then
           write(*,'(a,2(1x,i3),a,i1,a,d22.12)')'iter,j:',iter,j,', hu(',iSE,',1) : ',hu(iSE,1)
           stop
         endif
        enddo
      endif

      dp=hu(1,1)      ! dpm-1
      dt=hu(2,1)      ! dtetam-1
      dp1=hu(3,1)     ! dpm
      dt1=hu(4,1)     ! dtetam

! 0 < alph1 <= 1
      p(m)=p(m)+alph1*dp1
!           : p(couche m) modifie : pm + alpha1*dpm
      vp(m)=vp(m)+alph1*dp1
!           : teta(couche m) modifie
      t(m)=t(m)+alph1*dt1
!           : p(couche m-1) modifie : pm-1 - alpha1*dpm
      vt(m)=vt(m)+alph1*dt1
      p(j)=p(j)+alph1*dp
      vp(j)=vp(j)+alph1*dp
      t(j)=t(j)+alph1*dt
      vt(j)=vt(j)+alph1*dt
      dr=ar(j)*dp+br(j)*dt+ccr(j)
      ds=as(j)*dp+bs(j)*dt+ccs(j)
      r(j)=r(j)+alph1*dr
      s(j)=s(j)+alph1*ds
      gdr=dr
      gds=ds
      gdp=dp
      gdt=dt
      jgdr=m-1
      jgds=m-1
      jgdp=m-1
      jgdt=m-1

    endif

!=======================================================================
! Remontee de m-2 a 2.
!----------------------
! Pour 2 =< n =< m-2 :
!      dpn=U4n-1*dpn+1 + V4n-1*dtetan+1 + W4n-1
!      dtetan=U4n*dpn+1 + V4n*dtetan+1 + W4n
!      dxn+1=U4n+1*dpn+1 + V4n+1*dtetan+1 + W4n+1
!      dsn+1=U4n+2*dpn+1 + V4n+2*dtetan+1 + W4n+2
! Arrivee en j = 1.
!-------------------
! Les equations sont :
! dx1 = U1*dp2 +V1*dteta2 +W1
!        ...
! ds2 = U6*dp2 +V6*dteta2 +W6
      if (idebug > 0) then
        write(*,*) 'Henyey: layers treatment - back up'
      endif

      do while (j > 1)
       j = j-1
       dp1=dp
       dt1=dt
       dp=ap(j)*dp1+bp(j)*dt1+ccp(j)
       dt=at(j)*dp1+bt(j)*dt1+cct(j)
       if (j > 1) then
         dr=ar(j)*dp+br(j)*dt+ccr(j)
         ds=as(j)*dp+bs(j)*dt+ccs(j)
       else
         dr=ar(j)*dp1+br(j)*dt1+ccr(j)
         ds=as(j)*dp1+bs(j)*dt1+ccs(j)
       endif

! Stockage du numero de couche ou l'on a la plus importante correction
! en r, s, p, teta, et de la valeur de la plus grande correction.
       if (abs(dr)-abs(gdr) > 0.d0) then
         gdr=dr
         jgdr=j
       endif
       if (abs(ds)-abs(gds) > 0.d0) then
         gds= ds
         jgds=j
       endif
       if (abs(dp)-abs(gdp) > 0.d0) then
         gdp=dp
         jgdp=j
       endif
       if (abs(dt)-abs(gdt) > 0.d0) then
         gdt=dt
         jgdt=j
       endif
       p(j)=p(j)+alph1*dp
       vp(j)=vp(j)+alph1*dp
       t(j)=t(j)+alph1*dt
       vt(j)=vt(j)+alph1*dt
       r(j)=r(j)+alph1*dr
       s(j)=s(j)+alph1*ds
      enddo

!=======================================================================
! TRAITEMENT DE L'ENSEMBLE DES COUCHES
      if (idebug > 0) then
        write(*,*) 'Henyey: global'
      endif
      fred=0.5*(1.0+abs(gdt+vgdt)/(abs(gdt)+abs(vgdt)))
      vgdt=gdt
! : max(dteta1,...,dtetam)
      write(correction_message,'(a,4(i6,f10.5))') 'biggest correction p t r s:',&
                                  jgdp,gdp,jgdt,gdt,jgdr,gdr,jgds,gds
      if (itminc <= 1) then
! Ecriture du numero de couche ou l'on a la plus importante correction
! en r, s, p, teta, et de la valeur de cette plus grande correction.
        write(3,'(13x,a//)') trim(correction_message)

        if (iprc == 0) then
          if (iout > 0) then ! Si iout#0, on ajoute des couches au bord du noyau convectif
             write(3,'(1x,a,a,i4,a,f9.6)') 'RADIATIVE LEVEL ADJACENT TO CONVECTIVE CORE BOUNDARY AT', &
                                                                    ' IOUT=',iout,'  MCC/MR=',exphi(q(iout+1))
          endif
          write(3,'(1x,a,f6.3,2(2x,a,f7.4),a,1pe10.2)') 'CENTRAL VALUES RHOC=',rhoc,'TC=',tc,'X(M)=',xm,'  EGC=',egc
          write(3,'(1x,3f10.7,4e8.2/8e8.2)') x(1),y3(1),y(1),xc12(1),xc13(1),xn14(1),xn15(1),xo16(1),xo17(1), &
                                                        xo18(1),xne20(1),xne22(1),xmg24(1),xmg25(1),xmg26(1)
        endif
      else
! Ecriture du numero de couche ou l'on a la plus importante correction
! en r, s, p, teta, et de la valeur de cette plus grande correction.
        write(3,'(13x,a//)') trim(correction_message)
      endif
      gkor=max(gkor,abs(gdr),abs(gds),abs(gdp),abs(gdt))

      if (.not. henyey_last) then
        vsuminenv = vvsuminenv * (r(1)/vr(1))**2.0d0
        if (isnan(vsuminenv)) then
          write(*,*) 'vsuminenv=NaN'
          write(*,*) 'vvsuminenv,r(1),vr(1):',vvsuminenv,r(1),vr(1)
          rewind(222)
          write(222,*) nwmd,': vsuminenv=NaN'
          stop
        endif
        write(3,*) '--> adjustment of vsuminenv:',(r(1)/vr(1))**2.0d0
      endif

      if (iover /= 0) then
        call over1
      endif

      C12_save = xc12(1)
      C13_save = xc13(1)
      N14_save = xn14(1)
      O16_save = xo16(1)

!-------------------------------------------------------------------
! Calcul de la variation des abondances et de Omega
!--------------------------------------------------
! calcul des abondances de depart
! homogeneisee sur les nouvelles ZC
      if (idebug > 0) then
        write(*,*) 'call chemold'
      endif
      call chemold

! calcul des eps si ialflu=0 utilisation de netnew.
! Netnew considere les ZC comme une seule coquille pour les reactions nucleaires
      if (ialflu==0 .and. ichem==1 .and. idifcon==0) then
        if (idebug > 0) then
          write(*,*) 'call chemeps'
        endif
        call chemeps
      endif

! calcul du omega de depart, a partir des moment cinetiques du modele
! precedent et conservation locale du moment cinetique.
! Omega homogeneise dans les ZC
      if (irot == 1 .and. isol == 0) then
        if (idebug > 0) then
          write(*,*) 'call vomcon'
        endif
        call vomcon
        if (istati == 1) then
          if (idebug > 0) then
            write(*,*) 'call omesta'
          endif
          call omesta
        else
          if (idebug > 0) then
            write(*,*) 'call omenew'
          endif
          call omenew
        endif
        if (idebug > 0) then
          write(*,*) 'call dlonew'
        endif
        call dlonew
      endif

      if (ichem == 1 .and. idifcon == 0) then
        if (idebug > 0) then
          write(*,*) 'call chemeps'
        endif
        call chemeps
      endif
      if (idebug > 0) then
        write(*,*) 'call netnew'
      endif
      call netnew

      if (irot==1 .and. idiff/=0 .and. isol==0 .or. idifcon==1) then
        if (idebug > 0) then
          write(*,*) 'call coedif'
        endif
        call coedif
        if (idebug > 0) then
          write(*,*) 'call diffbr'
        endif
        call diffbr
      endif

      if (irot == 1 .and. isol == 0 .and. istati == 0) then
        if (idebug > 0) then
          write(*,*) 'call diffom'
        endif
        call diffom
        if (idebug > 0) then
          write(*,*) 'call advect'
        endif
        call advect
        if (idifcon == 0) then
          if (idebug > 0) then
            write(*,*) 'call omconv'
          endif
          call omconv
        endif
      endif

      !if (iledou == 1 .or. irot == 1 .or. idifcon == 1) then
      if (idebug > 0) then
        write(*,*) 'call grapmui'
      endif
      call grapmui
      !endif

      if (modanf > 1 .and. Add_Flux) then
        Flux_remaining = (BTotal_StartModel-dlelexsave-BTotal_EndAdvect)/dzeit
        if (idebug > 0) then
          write(*,*) 'Flux_remaining: ',Flux_remaining
          write(*,*) 'BTotal_StartModel: ', BTotal_StartModel
          write(*,*) 'dlelexsave: ', dlelexsave
          write(*,*) 'BTotal_EndAdvect: ', BTotal_EndAdvect
        endif
      else
        Flux_remaining = 0.d0
      endif

      if (plot .and. refresh .and. .not.Struc_Plotted) then
        Struc_Plotted = .true.
        call PlotStruc
      endif

! Si la valeur absolue de la plus grande correction est superieure
! a gkorm (valeur absolue de la plus grande correction toleree), il
! faut retourner dans main et repartir d'un autre modele.
      if (gkor >= gkorm .and. iter /= 1) then
         endIter=.true.
         exit
      endif
      if (gkor >= 1.25d0*gkorm .and. iter == 1) then
        endIter=.true.
        exit
      else if (gkor >= gkorm .and. iter == 1) then
        if (verbose) then
          write(*,*)'Premiere iteration: remise a zero de gkor'
        endif
        gkor = 0.d0
      endif

      if (iter < itminc) then
        exit
      endif

! On effectue au moins itminc iterations. C'est le nombre minimum
! d'iterations a effectuer dans henyey lorsque toutes les corrections
! calculees sont deja inferieures a gkorm, afin de prevenir une
! divergence.
! Si itminc = 1 (initialisation dans main), une seule iteration est
! necessaire : Calcul du modele definitif
!              ou approximation du modele suivant.

      if (phase < 3) then
        iterlim1=15
        iterlim2=15
      else
        iterlim1=20
        iterlim2=30
      endif
      if (max(abs(gg1),abs(gg2),abs(gg3),abs(gg4)) > 1.d-3 .and. itminc > 1 .and. iter <= iterlim1) then
        exit
      endif

      if (abs(gdr)>=agdr .or. abs(gds)>=agds .or. abs(gdp)>=agdp .or. abs(gdt)>=agdt) then
        if ((itminc <= 1.and.iter <= 1).or.(itminc > 1.and.iter <= iterlim1)) then
            exit
        endif
        if (max(abs(gg1),abs(gg2),abs(gg3),abs(gg4)) > 1.d2 .and. itminc > 1 .and. iter <= 50 ) then
          exit
        endif
      endif
      return

   enddo
  enddo
  return

  end subroutine henyey

!-----------------------------------------------------------------------
end module henyey_solver
