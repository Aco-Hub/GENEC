module energy

  use evol,only: ldi,kindreal
  use const,only: convMeVerg,cst_avo,cst_ecgs,pi,cst_k,cst_mh,cst_e
  use inputparam,only: phase,ialflu,ibasnet,ipop3,z,verbose,iapprox21
  use caramodele,only: gms,nwmd
  use abundmod,only: x,y3,y,xc12,xc13,xc14,xn14,xn15,xo16,xo17,xo18,xf18,xf19,xne20,xne21,xne22,xna23,xmg24,xmg25,xmg26, &
    xal26,xal27,xsi28,xprot,xneut,xbid,xbid1,eps,epsy,epsyy,epsyc,epsyo,epsc,b11,b33,b34,b112,b113,b114,b115a,b115g,b116, &
    b117a,b117g,b118a,b118g,epcne,eps20,c144,c184,c224,c134,c12ago,o16agn,epcna,b119a,b119g,b120,b121,b122,b123g,b123a, &
    b124,b125g,b125m,b1mg26,b1al26,b127g,b127a,e24ag,e17ag,e21ag,e18an,e21na,e25an,e20ng,e21ng,e22ng,e23ng,e24ng,e25ng, &
    e26ng,e27ng,e28ng,a26ga,a26gp,e14np,ec14pg,ec14ag,ef18na,e15ag,ef18np,e18pa,ec14ng,e19ap,e14be,e18be,e26be,e18ng,e17an, &
    c224g,e20ag,e12ng,e14ng,e19ng,enpp,encno,densityj,en13,eg,egp,egp1,egt,egt1,en,enue,enuet,enuet1,enuep,enuep1,epsn, &
    epsn1,epsp,epsp1,epst,epst1,tauxbe7pg,tauxbe7el,snu,ybe7,yb8,xnube7,xnub8,xnbrbe7,xnbrb8,b34neu,xfluxn,snube7,snub8, &
    eps_grav,eps_nu,nbelx,abelx,zabelx,nbael,nbzel
  use timestep,only: alter,dzeit
  use EOS,only: rh,rh1,rhp,rhp1,rht,rht1,rhe,rhpsi,rhpsip,rhpsit
  use strucmod,only: j,j1,m,t,p,vt,vp,zensi,beta,beta1,adi,adi1,adip,adip1,vmye,vmyo
  use nagmod,only: e02acf

  implicit none

  integer,parameter:: neutri=2,ippcno=0,itestx=0,maxz=105,maxel=100
! structure & time
  integer, save:: shellnb
  real(kindreal), save:: t9,rho,ane,pme,tmax

  real(kindreal),dimension(ldi),save:: taunucl
! elements network
  character (256), save:: namenet,namereac

  integer, save:: nbel,nbzmin,nbzmax,nbnmin,nbnmax,ineut,iprot,ialpha
  integer, dimension(maxel), save:: nbz,nba,nbn
  integer, dimension (0:maxz,0:maxz), save:: posel   !posel(Z,N)

  real(kindreal), dimension(maxel,maxel+1), save:: mata
  real(kindreal), dimension(maxel), save:: abunx,abuny
  real(kindreal), dimension(maxel), save:: vabunx,vabuny

  character (2), dimension (0:maxz), save:: SYM= (/ 'n ','H ','HE','LI','BE','B ','C ','N ','O ','F ','NE',&
      'NA','MG','AL','SI','P ','S ','CL','AR','K ','CA','SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',&
      'GA','GE','AS','SE','BR','KR','RB','SR','Y ','ZR','NB','MO','TC','RU','RH','PD','AG','CD','IN','SN',&
      'SB','TE','I ','XE','CS','BA','LA','CE','PR','ND','PM','SM','EU','GD','TB','DY','HO','ER','TM','YB',&
      'LU','HF','TA','W ','RE','OS','IR','PT','AU','HG','TL','PB','BI','PO','AT','RN','FR','RA','AC','TH',&
      'PA','U ','NP','PU','AM','CM','BK','CF','ES','FM','MD','NO','LR','RF','HA' /)
! reactions' network
  integer, parameter:: ngrid=70, nre=100
  integer, save::  ireac, kgrid
  integer, dimension (nre,4), save::  nsnb, elps

  real(kindreal), dimension (nre,ldi), save:: rrate
  real(kindreal), dimension (nre), save:: vrate,flag,qrad,qnu
  real(kindreal), dimension (ngrid), save:: tgrid
  real(kindreal), dimension (ngrid,nre,4), save:: vgrid
  real(kindreal), dimension (nre,4), save:: anb,znb

  character (6), dimension (nre,4), save:: zed
  character (37), dimension (nre), save:: reaction

  real(kindreal),dimension(ldi),save:: vmassen,rvect,t9n,pvect,epstot1,epsneut,dcoeff

! TAUX DE FOWLER ET AL. 1975 AVEC LES CORRECTIONS DE 1981 et 1988
! AT. DATA & NUCL. DATA TABLES 40, 283 (1988)
  real(kindreal),parameter:: cv2=0.9216d0,ca2=0.25d0,cvp2=1.6d-3,cap2=0.25d0,ecne=1.115d+18, &
    e23=9.980d+16,eo=1.142d+18,e20=2.246d+18,ec=1.728d+18,q11=6.666d0,q33=12.859d0,q34=1.587d0,&
    q17=9.547d0,qe7=17.347d0,q112=3.452d0,q113=7.551d0,q114=9.054d0,q115a=4.965d0,q115g=12.127d0,&
    q116=2.361d0,q117a=1.192d0,q117g=6.865d0,q118a=3.981d0,q118g=7.994d0,q119a=8.114d0,&
    q119g=12.845d0,q120=2.432d0,q121=6.740d0,q122=8.793d0,q123g=11.69d0,q123a=2.377d0,q124=2.271d0, &
!--- Q125G AND Q125M ARE GIVEN FOR THE GROUND STATE AND RESPECTIVELY
!       FOR THE METASTABLE STATE OF AL26
    q125g=6.306d0,q125m=6.078d0,q1mg26=8.271d0,q1al26=7.465d0,q127g=11.586d0,q127a=1.601d0,&
    xkmc2=1.68628d-10,q12g=5.494d0,q22n=3.269d0,q16=4.020d0,q17a=17.348d0,q17g=17.254d0
  real(kindreal),dimension(5),parameter:: dojl=(/0.d0,0.d0,0.d0,0.d0,0.d0/), &
    doji=(/-1.135d+08,1.256d+08,5.149d+07,3.436d+07,1.005d+07/),&
    dojh=(/4.724d+08,2.976d+08,2.242d+08,7.937d+07,4.859d+07/),&
    d1jl=(/-1.879d+10,-9.667d+09,-5.602d+09,-3.370d+09,-1.825d+09/),&
    d1ji=(/1.652d+09,-3.119d+09,-1.839d+09,-1.458d+09,-8.956d+08/),&
    d1jh=(/-7.094d+11,-3.697d+11,-2.189d+11,-1.273d+11,-5.705d+10/),&
    d2jl=(/-2.919d+10,-1.185d+10,-7.270d+09,-4.222d+09,-1.560d+09/),&
    d2ji=(/-1.549d+10,-9.338d+09,-5.899d+09,-3.035d+09,-1.598d+09/),&
    d2jh=(/-2.254d+10,-1.551d+10,-7.793d+09,-4.489d+09,-2.185d+09/),&
    zwy=(/12.d0,17.123d0,10.288d0,3.343d0,-1.46d0/),&
    zwc=(/12.d0,16.192d0,9.014d0,2.578d0,-0.889d0/),&
    zwo=(/16.d0,20.978d0,11.241d0,3.025d0,-0.946d0/),&
    zw12=(/36.d0,45.6635d0,23.275d0,5.668d0,-1.362d0/),&
    zw20=(/20.d0,25.616d0,13.308d0,3.409d0,-0.9875d0/),&
    zw14=(/14.d0,18.606d0,10.150d0,2.811d0,-.919d0/),&
    zwcp=(/6.d0,8.302d0,4.804d0,1.488d0,-0.643d0/),&
    zw18=(/8.d0,10.716d0,5.941d0,1.721d0,-0.673d0/),&
    zw19=(/18.d0,23.314d0,12.291d0,3.223d0,-0.968d0/),&
    zw24=(/24.d0,30.136d0,15.250d0,3.749d0,-1.020d0/),&
    zwn=(/7.d0,9.520d0,5.385d0,1.609d0,-0.659d0/)
! Neutrinos energy loss from Itoh & al. (1989). Same analytical formulae
! and same numerical factors in Itoh & al. (1996).
  real(kindreal),dimension(7),parameter:: cojl=(/1.008d+11,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),&
    coji=(/9.889d+10,-4.524d+08,-6.088d+06,4.269d+07,5.172d+07,4.910d+07,4.388d+07/),&
    cojh=(/9.581d+10,4.107d+08,2.305d+08,2.236d+08,1.580d+08,2.165d+08,1.721d+08/),&
    c1jl=(/8.156d+10,9.728d+08,-3.806d+09,-4.384d+09,-5.774d+09,-5.249d+09,-5.153d+09/),&
    c1ji=(/1.813d+11,-7.556d+09,-3.304d+09,-1.031d+09,-1.764d+09,-1.851d+09,-1.928d+09/),&
    c1jh=(/1.459d+12,1.314d+11,-1.169d+11,-1.765d+11,-1.867d+11,-1.983d+11,-1.896d+11/),&
    c2jl=(/1.067d+11,-9.782d+09,-7.193d+09,-6.936d+09,-6.893d+09,-7.041d+09,-7.193d+09/),&
    c2ji=(/9.750d+10,3.484d+10,5.199d+09,-1.695d+09,-2.865d+09,-3.395d+09,-3.418d+09/),&
    c2jh=(/2.424d+11,-3.669d+09,-8.691d+09,-7.967d+09,-7.932d+09,-7.987d+09,-8.333d+09/)

private
public :: energ,enint,netinit,netburning,nucal
public :: nbel,ireac,rrate
public :: vmassen,rvect,t9n,pvect,epstot1,epsneut,dcoeff
public :: namenet,namereac,nbz,nbn,ineut,iprot,ialpha,nbzmax,nbzmin,nbnmax,nbnmin,taunucl
public :: abuny,posel,tmax,abunx,nba,vabuny,vabunx,t9,pme,ane,rho,itestx,mata,shellnb,maxel
public :: kgrid,ngrid,vgrid,tgrid,nre,anb,znb,nsnb,elps,zed,vrate,qrad,qnu,flag,reaction,maxz,sym

contains
!======================================================================
subroutine energ
!-----------------------------------------------------------------------
  use SmallFunc,only: primus,drimus,secun,dsecun,tertiu,dterti

  implicit none

  integer:: k,nbouc=0,nj,njj,kp,ii,n1216,n1620,n2024
  real(kindreal):: convMeVAvo,cstnumdebye,conver,hfak1,hfak2,hfak3,hfak4,hfak5,zeta,sz2,t6ln,t8=0.d0,t8ln,t9,t9ln,t912=0.d0,&
    t913=0.d0,t923=0.d0,t932=0.d0,t943,t953=0.d0,t915,t954,t9a,dt9a,dt9a1,dt9a2,ft9a,dft9a,fpt9a,dfpt9a,gt9,dgt9,gpt9,dgpt9,&
    gam,gamlnt,qgam,glam,dglam,e2lam,ppp,rho,rhom,rhos9,rhom2,rhom3,rhotp,rhopt,ttt,xlttt,uxlam,uxlam2,uxlam3,uxlam4,xlam,xlam2,&
    xlam3,xlam4,xlam5,xlam6,xlam7,xlam8,xxi,xxitp,xxipt,uno,due,tre,quat,cinq,cinqa,cinqb,six,sixa,sixb,sept,septa,septb,huit,&
    neuf,cent,centa,duno,ddue,dtre,dquat,dcinq,dcinqa,dcinqb,dsix,dsixa,dsixb,dsept,dsepta,dseptb,dhuit,dneuf,dcent,tcent,u,&
    aa,al2,aln,ba,bb,bg,b1,b2,b3,b4,b5,b6,bet_en,bet2,b3f,cplus,cmoins,cna3,cc,cya,dd,d1,d2,d3,d4,d5,d6,dsnd,dsnc,dsnb,dsna,&
    dy,dy2,eb,ebr,ebp,ebt,ecp,ect,ee,ef,eot,eop,eyp=0.d0,eyt=0.d0,fb,fcp,fcpt,fcpp,ffa,ffat,ffap,fn,fnp,fnpt,fnpp,fop,fopt,fopp,&
    fp,fw,fy,fyt,fyp,fyc,fytc,fypc,fyo,fyto,fypo,fy12,fyt12,fyp12,fyn14,fyt14,fyp14,fy20,fyt20,fyp20,fyn22,fyt22,fyp22,fy24,fyt24,&
    fyp24,pf,f33,f112,f116,f1be7,be7pg,dbe7pg,be7el,dbe7el,ch4ab8,ch4abt,cbe8ag,dbe8ag,c12c12p,c12c12a,ect12,ecp12,ec14at,ec14ap,&
    ec14nt,ec14pt,ec14pp,e144,e144t,e144p,e14npt,e14npp,e15agt,e15agp,e17anp,e17ant,e184,e184t,e184p,ef18nt,f18nap,f18npt,f18npp,&
    e224,e224t,e224p,e224g,e22tg,e20agt,e20agp,e134,e13t,e13p,e18ngt,e18ngp,e18pat,e18pap,e18ant,e18anp,e19apt,e19app,e24agt,&
    e24agp,e17agt,e17agp,e20agw,e20t1,e20p1,e21agt,e21agp,e21nat,e21nap,e25ant,e25anp,e28ngh=0.d0,e28nht=0.d0,e28nhp=0.d0,e28ngt,&
    e28ngp,al26tn,dal26t,al26mn,dal26m,a26gat,a26gap,a26gpt,a26gpp,wpsyy=0.d0,wpsyc,wpsyo,w12ng=0.d0,wc14ag,wc14ng,wc14pg,&
    w14ng=0.d0,w14np,w15ag,w17ag,w17an,w18an,w18ng,w18pa,wf18na,wf18np,w19ap,w19ng=0.d0,w20ag,w20ng=0.d0,w21na,w21ng=0.d0,w24ag,&
    w21ag,w22ng=0.d0,w23ng=0.d0,w24ng=0.d0,w25an,w25ng=0.d0,w26ga,w26gp,w26ng=0.d0,w27ng=0.d0,w28ng,w14be=0.d0,w18be=0.d0,&
    w26be=0.d0,epsp2,epst2,ep23,e23p1,e23t1,etot,etott,etotp,eptrg,epg,eptra,epa,revra,drevra,revrat,fluneu,flunet,flunep,sna,snb,&
    snc,snd,sntot,tsntot,qpai1,qpai2,qpai3,qpai,qpabps,qpair,qpho3,qpho4,qpho5,qpho,qpho6,qphot,qpla1,qplas,qtot,dtpai1,dtpai2,&
    dtpai3,dtpai,dtpaib,dtpair,dtpho3,dtpho4,dtpho5,dtpho6,dtpho,dtphot,dtplas,dtqtot,dppai3,dppai,dppaib,dppair,dppho5,dppho,&
    dppho6,dpphot,dpplas,dpqtot,dqpt,dqtp,taup,xe1,xe2,xe3,xe4,xe5,xious1,dkiou1,xious2,dkiou2,xious3,dkiou3,xntp,xnpt,xdtp,xdpt

  real(kindreal),dimension(20):: yab=0.d0
  real(kindreal),dimension(ldi):: zensi2
  real(kindreal),dimension(23):: ep,eprt,eptr
  real(kindreal),dimension(3):: cna1,cna2,cnao,cnb1,cnb2,cnb3,cnc,wn,wd,wf,wftp,wfpt
!-----------------------------------------------------------------------
!    IN EPS THE NUMERICAL FACTOR IS 1.6022E-6*AVOGADRO=9.649E+17
  convMeVAvo = convMeVerg * cst_avo
! ippcno defini dans netb.inc
  if(ippcno == 1) then
    enpp(j1)=0.d0
    encno(j1)=0.d0
    densityj(j1)=exp(rh1)
  endif

  eg=0.d0
  egp1=0.d0
  egp=0.d0
  egt1=0.d0
  egt=0.d0
  eps_grav(j1) = 0.0d0
  eps_nu(j1) = 0.0d0

  if (j1 /= 1 .and. alter > 0.d0) then
    hfak1=0.25d0*exp(0.5d0*(p(j1)+p(j)-rh1-rh))/dzeit
    hfak2=vp(j)+vp(j1)
    hfak3=vt(j1)+vt(j)
    eg=hfak1*((-rht1-rht)*hfak2+(rht1/adi1+rht/adi)*hfak3)
    hfak4=-4.d0*(1.d0-beta1)/(beta1**2.d0)*(hfak2-hfak3/adi1)-rht1*adip1*hfak3/(adi1**2d0)
    hfak5=-4.d0*(1.d0-beta)/(beta**2.d0)*(hfak2-hfak3/adi)-rht*adip*hfak3/(adi**2.d0)
    egp1=0.5d0*eg*(1.d0-rhp1)+hfak1*(hfak4-rht1-rht)
    egp=0.5d0*eg*(1.d0-rhp)+hfak1*(hfak5-rht1-rht)
    egt1=-0.5d0*eg*rht1+hfak1*(-4.d0*hfak4+rht1/adi1+rht/adi)
    egt=-0.5d0*eg*rht+hfak1*(-4.d0*hfak5+rht1/adi1+rht/adi)
  endif
  eps_grav(j1) = eg

  ep(:)=0.d0
  eprt(:)=0.d0
  eptr(:)=0.d0

  b11(j1)=0.d0
  b33(j1)=0.d0
  b34(j1)=0.d0
  b34neu(j1)=0.d0
  tauxbe7pg(j1)=0.d0
  tauxbe7el(j1)=0.d0

  b112(j1)=0.d0
  b113(j1)=0.d0
  b114(j1)=0.d0
  b115a(j1)=0.d0
  b115g(j1)=0.d0
  b116(j1)=0.d0
  b117a(j1)=0.d0
  b117g(j1)=0.d0
  b118a(j1)=0.d0
  b118g(j1)=0.d0
  if(ialflu == 1) then
    xfluxn(j1)=0.d0
    b119a(J1)=0.d0
    b119g(J1)=0.d0
    b120(J1)=0.d0
    b121(J1)=0.d0
    b122(J1)=0.d0
    b123g(J1)=0.d0
    b123a(J1)=0.d0
    b124(J1)=0.d0
    b125g(J1)=0.d0
    b125m(J1)=0.d0
    b1mg26(J1)=0.d0
    b1aL26(J1)=0.d0
    b127g(J1)=0.d0
    b127a(J1)=0.d0
    e24ag(j1)=0.d0
    e17ag(j1)=0.d0
    e21ag(j1)=0.d0
    ec14ag(j1)=0.d0
    e15ag(j1)=0.d0
    e18an(j1)=0.d0
    e25an(j1)=0.d0
    e21na(j1)=0.d0
    ef18na(j1)=0.d0
    e12ng(j1)=0.d0
    e14ng(j1)=0.d0
    e19ng(j1)=0.d0
    e20ng(j1)=0.d0
    e21ng(j1)=0.d0
    e22ng(j1)=0.d0
    e23ng(j1)=0.d0
    e24ng(j1)=0.d0
    e25ng(j1)=0.d0
    e26ng(j1)=0.d0
    e27ng(j1)=0.d0
    e28ng(j1)=0.d0
    e18ng(j1)=0.d0
    ec14ng(j1)=0.d0
    e14np(j1)=0.d0
    ef18np(j1)=0.d0
    a26ga(j1)=0.d0
    a26gp(j1)=0.d0
    ec14pg(j1)=0.d0
    e18pa(j1)=0.d0
    e19ap(j1)=0.d0
    e14be(j1)=0.d0
    e18be(j1)=0.d0
    e26be(j1)=0.d0
  endif

  eps(j1)= 0.d0
  epsy(j1)= 0.d0
  epsyy(j1)=0.d0
  epsyc(j1)=0.d0
  epsyo(j1)=0.d0
  e17an(j1)=0.d0
  epsc(j1)=0.d0
  epcne(j1)=0.d0
  eps20(j1)=0.d0
  c144(j1)=0.d0
  c184(j1)=0.d0
  c224(j1)=0.d0
  c224g(j1)=0.d0
  e20ag(j1)=0.d0
  c134(j1)=0.d0

  epsp1= 0.d0
  epst1= 0.d0
  en= 0.d0
  epsn1 = 0.d0
  eps_nu(j1) = 0.d0
  enuet1 = 0.d0
  enuep1 = 0.d0
  enue = 0.d0
  eb=0.d0
  ebr=0.d0
  ebp=0.d0
  ebt=0.d0
  zensi(j1)=1.d0
! modif PopIII: rajout des variables epsp2, epst2 et zensi2(j1)
  epsp2=0.d0
  epst2=0.d0
  zensi2(j1)=1.d0

  yab(1)=x(j1)
  if (ipop3 == 1) then
    if (x(j1) <= 0.d0 .and. y(j1) <= 0.d0) then
      goto 27
    endif
  endif
  if (x(j1) <= 0.d0) then
    go to 23
  endif
  if ((t(j1)-log(4.d6)) < 0.d0) then
    return
  endif

  t6ln= t(j1)-log(1.d6)
  t9ln=t(j1)-log(1.d9)
  t9=exp(t9ln)
  t913=exp(t9ln/3.d0)
  t923=t913*t913
  t932=t9**1.5d0
  t943=t913*t9
  t953=t913*t943
  t912=sqrt(t9)
  t915=EXP(t9ln/5.d0)
  t954=t9**1.25d0

!--- PP-CHAINS AND CNO-TRICYCLE ASSUME H2(P,G)HE3 INSTANTANEOUS
!--- ASSUME WEAK ELECTRON SCREENING
! zeta = sum (X_i / A_i) * Z_i * (Z_i + 1) (cf Maeder 2009 page 205).
  zeta=2.d0*(x(j1)+y3(j1))+1.5d0*y(j1)+3.5d0*xc12(j1)+4.d0*(xn14(j1)+xo18(j1))+4.5d0*xo16(j1)+42.d0/13.d0*xc13(j1)+ &
         56.d0/15.d0*xn15(j1)+72.d0/17.d0*xo17(j1)
  if (ialflu == 1) then
    zeta=zeta+5.5d0*xne20(j1)+5.d0*xne22(j1)+6.5d0*xmg24(j1)+6.d0*xmg26(j1)+7.d0*xal26(j1)+7.5d0*xsi28(j1)+ &
          110.d0/21.d0*xne21(j1)+132.d0/23.d0*xna23(j1)+156.d0/25.d0*xmg25(j1)+182.d0/27.d0*xal27(j1)
  endif

! V_0 / (kT Z_i Z_j) : cf. Maeder 2009, page 206, eq. 9.49
! V_0/(kT) = 2 sqrt(pi) e^3 / (k^(3/2) sqrt(m_H) (10^6)^(3/2))
!           * Z_i Z_j sqrt(zeta rho) / T_6^(3/2)
  cstnumDebye = cst_ecgs**3.d0*2.d0*sqrt(pi)/(cst_k**(3.d0/2.d0)*sqrt(cst_mh)* 1.d9)
  fp=cstnumDebye*exp(0.5d0*rh1-1.5d0*t6ln)*sqrt(zeta)
  f33=4.d0*fp
  f112=6.d0*fp
  fn=7.d0*fp
  f116=8.d0*fp
  f1be7=4.d0*fp

!--- IN THE FOLLOWING GAM=NOMBRE PP III/NOMBRE PP II
!    IN EPS THE NUMERICAL FACTOR IS 1.6022E-6*AVOGADRO=9.649E+17
!--- IF ALL THE Q'S ARE IN MEV.
  yab(2)=y3(j1)/3.d0
  yab(3)=y(j1)/4.d0
  yab(4)=xc12(j1)/12.d0
  yab(5)=xc13(j1)/13.d0
  yab(6)=xn14(j1)/14.d0
  yab(7)=xn15(j1)/15.d0
  yab(8)=xo16(j1)/16.d0
  yab(9)=xo17(j1)/17.d0
  yab(10)=xo18(j1)/18.d0
  if (ialflu == 1) then
    yab(11)=xf19(J1)/19.d0
    yab(12)=xne20(J1)/20.d0
    yab(13)=xne21(J1)/21.d0
    yab(14)=xne22(J1)/22.d0
    yab(15)=xna23(J1)/23.d0
    yab(16)=xmg24(J1)/24.d0
    yab(17)=xmg25(J1)/25.d0
    yab(18)=xmg26(J1)/26.d0
    yab(19)=xal26(J1)/26.d0
    yab(20)=xal27(J1)/27.d0
  endif

! ----------------------------------------------------------------------
! H1(P,E+NU)H2 AT. DATA & NUCL. DATA TABLES 40, 283 (1988)
! pour ep(x),nacre pour taux de reaction
  aa=4.01d-15
  bb=3.380d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  uno=primus(aa,bb,xe1,xe2,t9)
  duno=drimus(bb,xe1,xe2,t9)
  aa=1.d0
  bb=0.123d0
  cc=1.09d0
  dd=0.938d0
  ee=0.d0
  fw=0.d0
  xe1=1.d0/3.d0
  xe2=2.d0/3.d0
  xe3=1.0d0
  xe4=0.d0
  xe5=0.d0
  due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  cent=uno*due
  dcent=uno*duno*due+uno*ddue

  call interpol(41,t9,u)   ! fichier='ppnu'
  tcent=u
  b11(j1)=exp(rh1+fp)*tcent
  ep(1)=0.5d0*yab(1)*yab(1)*q11*b11(j1)
  eprt(1)=1.d0+0.5d0*fp
  eptr(1)=-1.5d0*fp+dcent/cent

!  HE3(HE3,2P)HE4 REF IDEM H1(P,E+NU)H2
  aa=6.04d+10
  bb=12.276d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  uno=primus(aa,bb,xe1,xe2,t9)
  duno=drimus(bb,xe1,xe2,t9)
  aa=1.d0
  bb=0.034d0
  cc=-0.522d0
  dd=-0.124d0
  ee=0.353d0
  fw=0.213d0
  xe1=1.d0/3.d0
  xe2=2.d0/3.d0
  xe3=1.0d0
  xe4=4.d0/3.d0
  xe5=5.d0/3.d0
  due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  cent=uno*due
  dcent=uno*duno*due+uno*ddue

  call interpol(13,t9,u)   ! fichier='he3h'
  tcent=u
  b33(j1)=exp(rh1+f33)*tcent
  ep(2)=0.5d0*yab(2)*yab(2)*q33*b33(j1)
  eprt(2)=1.d0+0.5d0*f33
  eptr(2)=-1.5d0*f33+dcent/cent

! HE4(HE3,G)BE7 REF IDEM H1(P,E+NU)H2
  call interpol(12,t9,u)   !  fichier='he3a'
  tcent=u
  t9a=t9/(1.d0+4.95d-02*t9)
  cent=5.61d+06*t9a**(5.d0/6.d0)/t9**(3.d0/2.d0)*exp(-12.826d0/t9a**(1.d0/3.d0))
  dt9a=1.d0-4.95d-02*t9a
  dcent=5.d0/6.d0*dt9a-3.d0/2.d0+1.d0/3.d0*12.826d0/t9a**(1.d0/3.d0)*dt9a
  b34(j1)=exp(rh1+f33)*tcent

! BE7(P,G)B8 REF IDEM H1(P,E+NU)H2
! Caughlan and Fowler (1988) (cf article p.285)
! avec S17(0)=0.024keV-b (Filippone 1986)
  aa=3.11d+05
  bb=10.262d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  uno=primus(aa,bb,xe1,xe2,t9)
  duno=drimus(bb,xe1,xe2,t9)
! Caughlan and Fowler (1988) (cf article p.285)
! avec S17(0)=0.024keV-b (Filippone 1986)
  aa=2.53d+03
  bb=7.306d0
  xe1=3.d0/2.d0
  xe2=1.d0
  due=primus(aa,bb,xe1,xe2,t9)
  ddue=drimus(bb,xe1,xe2,t9)

  call interpol(5,t9,u)   ! fichier='be7g'
  be7pg=u
  tauxbe7pg(j1)=be7pg*exp(rh1+f1be7)
  b34neu(j1)=b34(j1)
  dbe7pg=uno*duno+due*ddue

! BE7(E-,NU+G)LI7 REF IDEM H1(P,E+NU)H2
  aa=0.0027d0
  bb=-2.515d-03
  xe1=1.d0
  xe2=1.d0
  uno=primus(aa,bb,xe1,xe2,t9)
  duno=drimus(bb,xe1,xe2,t9)
  be7el=1.34d-10/(t9**0.5d0)*(1.d0-0.537d0*(t9**(1.d0/3.d0))+3.86d0*(t9**(2.d0/3.d0))+uno)
  tauxbe7el(j1)=be7el*exp(rh1)
  dbe7el=-0.5d0*be7el+1.34d-10/(t9**0.5d0)*(-0.537d0*1.d0/3.d0*(t9**(1.d0/3.d0))+3.86d0*2.d0/3.d0*(t9**(2.d0/3.d0))+uno*duno)
  gam=x(j1)/(1.d0+x(j1))*2.d0*exp(f33)*be7pg/be7el
  if(be7pg==0) then
    gamlnt=-1.5d0*f33-dbe7el/be7el
  else
    gamlnt=-1.5d0*f33+dbe7pg/be7pg-dbe7el/be7el
  endif
  ep(3)=yab(2)*yab(3)*(q34+(gam*q17+qe7)/(1.d0+gam))*b34(j1)
  qgam=gam*(q17-qe7)/((1.d0+gam)*(q34+qe7+gam*(q34+q17)))
  eprt(3)=1.d0+0.5d0*f33*(1.d0+qgam)
  eptr(3)=-1.5d0*f33+dcent+qgam*gamlnt

! C12(P,G)N13 AT. DATA & NUCL. DATA TABLES 40, 283 (1988)/Nacre
  aa=2.04d+07
  bb=13.690d0
  cc=1.500d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  xe3=2.d0
  uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
  duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
  aa=1.d0
  bb=0.030d0
  cc=1.19d0
  dd=0.254d0
  ee=2.06d0
  fw=1.12d0
  xe1=1.d0/3.d0
  xe2=2.d0/3.d0
  xe3=1.0d0
  xe4=4.d0/3.d0
  xe5=5.d0/3.d0
  due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  aa=1.08d+05
  bb=4.925d0
  xe1=3.d0/2.d0
  xe2=1.d0
  tre=primus(aa,bb,xe1,xe2,t9)
  dtre=drimus(bb,xe1,xe2,t9)
  aa=2.15d+05
  bb=18.179d0
  xe1=3.d0/2.d0
  xe2=1.d0
  quat=primus(aa,bb,xe1,xe2,t9)
  dquat=drimus(bb,xe1,xe2,t9)
  cent=uno*due+tre+quat
  dcent=uno*duno*due+uno*ddue+tre*dtre+quat*dquat

  call interpol(7,t9,u)   ! fichier='c12g'
  tcent=u
  b112(j1)=exp(rh1+f112)*tcent
  ep(4)=yab(1)*yab(4)*q112*b112(j1)
  eprt(4)=1.d0+0.5d0*f112
  eptr(4)=-1.5d0*f112+dcent/cent

! C13(P,G)N14 AT. DATA & NUCL. DATA TABLES 40, 283 (1988)/Nacre
  aa=8.01d+07
  bb=13.717d0
  cc=2.d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  xe3=2.d0
  uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
  duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
  aa=1.d0
  bb=0.030d0
  cc=0.958d0
  dd=0.204d0
  ee=1.39d0
  fw=0.753d0
  xe1=1.d0/3.d0
  xe2=2.d0/3.d0
  xe3=1.0d0
  xe4=4.d0/3.d0
  xe5=5.d0/3.d0
  due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  aa=1.21d+06
  bb=5.701d0
  xe1=6.d0/5.d0
  xe2=1.d0
  tre=primus(aa,bb,xe1,xe2,t9)
  dtre=drimus(bb,xe1,xe2,t9)
  cent=uno*due+tre
  dcent=uno*duno*due+uno*ddue+tre*dtre

  call interpol(8,t9,u)   ! fichier='c13g'
  tcent=u
  b113(j1)=exp(rh1+f112)*tcent
  ep(5)=yab(1)*yab(5)*q113*b113(j1)
  eprt(5)=1.d0+0.5d0*f112
  eptr(5)=-1.5d0*f112+dcent/cent

! N14(P,G)O15 AT. DATA & NUCL. DATA TABLES 40, 283 (1988)/Nacre
  aa=4.90d+07
  bb=15.228d0
  cc=3.294d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  xe3=2.d0
  uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
  duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
  aa=1.d0
  bb=0.027d0
  cc=-0.778d0
  dd=-0.149d0
  ee=0.261d0
  fw=0.127d0
  xe1=1.d0/3.d0
  xe2=2.d0/3.d0
  xe3=1.0d0
  xe4=4.d0/3.d0
  xe5=5.d0/3.d0
  due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  aa=2.37d+03
  bb=3.011d0
  xe1=3.d0/2.d0
  xe2=1.d0
  tre=primus(aa,bb,xe1,xe2,t9)
  dtre=drimus(bb,xe1,xe2,t9)
  aa=2.19d+04
  bb=12.530d0
  xe1=0.d0
  xe2=1.d0
  quat=primus(aa,bb,xe1,xe2,t9)
  dquat=drimus(bb,xe1,xe2,t9)
  cent=uno*due+tre+quat
  dcent=uno*duno*due+uno*ddue+tre*dtre+quat*dquat

  call interpol(20,t9,u)   ! fichier='n14g'
  tcent=u
  b114(j1)=exp(rh1+fn)*tcent
  ep(6)=yab(1)*yab(6)*q114*b114(j1)
  eprt(6)=1.d0+0.5d0*fn
  eptr(6)=-1.5d0*fn+dcent/cent

! N15(P,G)O16 AT. DATA & NUCL. DATA TABLES 40, 283 (1988)/Nacre
  aa=9.78d+08
  bb=15.251d0
  cc=0.450d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  xe3=2.d0
  uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
  duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
  aa=1.d0
  bb=0.027d0
  cc=0.219d0
  dd=0.042d0
  ee=6.83d0
  fw=3.32d0
  xe1=1.d0/3.d0
  xe2=2.d0/3.d0
  xe3=1.0d0
  xe4=4.d0/3.d0
  xe5=5.d0/3.d0
  due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  aa=1.11d+04
  bb=3.328d0
  xe1=3.d0/2.d0
  xe2=1.d0
  tre=primus(aa,bb,xe1,xe2,t9)
  dtre=drimus(bb,xe1,xe2,t9)
  aa=1.49d+04
  bb=4.665d0
  xe1=3.d0/2.d0
  xe2=1.d0
  quat=primus(aa,bb,xe1,xe2,t9)
  dquat=drimus(bb,xe1,xe2,t9)
  aa=3.80d+06
  bb=11.048d0
  xe1=3.d0/2.d0
  xe2=1.d0
  cinq=primus(aa,bb,xe1,xe2,t9)
  dcinq=drimus(bb,xe1,xe2,t9)
  cent=uno*due+tre+quat+cinq
  dcent=uno*duno*due+uno*ddue+tre*dtre+quat*dquat+cinq*dcinq

  call interpol(22,t9,u)   ! fichier='n15g'
  tcent=u
  b115g(j1)=exp(rh1+fn)*tcent
  ep(7)=yab(1)*yab(7)*q115g*b115g(j1)
  eprt(7)=1.d0+0.5d0*fn
  eptr(7)=-1.5d0*fn+dcent/cent

! N15(P,A)C12 AT. DATA & NUCL. DATA TABLES 40, 283 (1988)/Nacre
  aa=1.08d+12
  bb=15.251d0
  cc=0.522d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  xe3=2.d0
  uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
  duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
  aa=1.d0
  bb=0.027d0
  cc=2.62d0
  dd=0.501d0
  ee=5.36d0
  fw=2.60d0
  xe1=1.d0/3.d0
  xe2=2.d0/3.d0
  xe3=1.0d0
  xe4=4.d0/3.d0
  xe5=5.d0/3.d0
  due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  aa=1.19d+08
  bb=3.676d0
  xe1=3.d0/2.d0
  xe2=1.d0
  tre=primus(aa,bb,xe1,xe2,t9)
  dtre=drimus(bb,xe1,xe2,t9)
  aa=5.41d+08
  bb=8.926d0
  xe1=1.d0/2.d0
  xe2=1.d0
  quat=primus(aa,bb,xe1,xe2,t9)
  dquat=drimus(bb,xe1,xe2,t9)
  aa=4.72d+08
  bb=7.721d0
  xe1=3.d0/2.d0
  xe2=1.d0
  cinq=primus(aa,bb,xe1,xe2,t9)
  dcinq=drimus(bb,xe1,xe2,t9)
  aa=2.20d+09
  bb=11.418d0
  xe1=3.d0/2.d0
  xe2=1.d0
  six=primus(aa,bb,xe1,xe2,t9)
  dsix=drimus(bb,xe1,xe2,t9)
  cent=uno*due+tre+quat+0.1d0*(cinq+six)
  dcent=uno*duno*due+uno*ddue+tre*dtre+quat*dquat+0.1d0*(cinq*dcinq+six*dsix)

  call interpol(21,t9,u)   ! fichier='n15a'
  tcent=u
  b115a(j1)=exp(rh1+fn)*tcent
  ep(8)=yab(1)*yab(7)*q115a*b115a(j1)
  eprt(8)=1.d0+0.5d0*fn
  eptr(8)=-1.5d0*fn+dcent/cent

! O16(P,G)F17
  aa=1.5d+08
  bb=16.692d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  uno=primus(aa,bb,xe1,xe2,t9)
  duno=drimus(bb,xe1,xe2,t9)
  due=1.d0+2.13d0*(1.d0-exp(-0.728d0*t9**(2.d0/3.d0)))
  cent=uno/due
  ddue=2.13d0*0.728d0*2.d0/3.d0*t9**(2.d0/3.d0)*exp(-0.728d0*t9**(2.d0/3.d0))/due
  dcent=duno-ddue

  call interpol(33,t9,u)   ! fichier='o16g'
  tcent=u
  b116(j1)=exp(rh1+f116)*tcent
  ep(9)=yab(1)*yab(8)*q116*b116(j1)
  eprt(9)=1.d0+0.5d0*f116
  eptr(9)=-1.5d0*f116+dcent

! O17(P,A)N14 AT. LANDRE ET AL. A&A 240 85 (1990)/nacre
  aa=1.53d+07
  bb=16.712d0
  cc=0.565d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  xe3=2.d0
  uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
  duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
  aa=1.d0
  bb=0.025d0
  cc=5.39d0
  dd=0.940d0
  ee=13.5d0
  fw=5.98d0
  xe1=1.d0/3.d0
  xe2=2.d0/3.d0
  xe3=1.0d0
  xe4=4.d0/3.d0
  xe5=5.d0/3.d0
  due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  aa=2.92d+06
  bb=4.247d0
  xe1=-1.d0
  xe2=1.d0
  tre=primus(aa,bb,xe1,xe2,t9)
  dtre=drimus(bb,xe1,xe2,t9)
  aa=1.78d+05
  bb=16.67d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  quat=primus(aa,bb,xe1,xe2,t9)
  dquat=drimus(bb,xe1,xe2,t9)
  quat=quat/(0.479d0*t9**(2.d0/3.d0)+3.12d-3)**2.d0
  dquat=dquat-2.d0*0.479d0*2.d0/3.d0*t9**(2.d0/3.d0)/(0.479d0*t9**(2.d0/3.d0)+3.12d-3)
  aa=2.8d+11
  bb=16.67d0
  cc=0.040d0
  xe1=-1.d0
  xe2=1.d0/3.d0
  xe3=2.d0
  cinq=secun(aa,bb,cc,xe1,xe2,xe3,t9)
  dcinq=dsecun(bb,cc,xe1,xe2,xe3,t9)
  aa=2.94d-03
  bb=0.767d0
  xe1=3.d0/2.d0
  xe2=1.d0
  six=primus(aa,bb,xe1,xe2,t9)
  dsix=drimus(bb,xe1,xe2,t9)
  aa=98.d0
  bb=2.077d0
  xe1=3.d0/2.d0
  xe2=1.d0
  sept=primus(aa,bb,xe1,xe2,t9)
  dsept=drimus(bb,xe1,xe2,t9)
  cent=uno*due+tre+quat+0.5d0*(cinq+six)+0.1d0*sept
  dcent=uno*duno*due+uno*ddue+tre*dtre+quat*dquat+0.5d0*(cinq*dcinq+six*dsix)+0.1d0*sept*dsept

  call interpol(34,t9,u)   ! fichier='o17a'
  tcent=u
  b117a(j1)=exp(rh1+f116)*tcent
  ep(10)=yab(1)*yab(9)*q117a*b117a(j1)
  eprt(10)=1.d0+0.5d0*f116
  eptr(10)=-1.5d0*f116+dcent/cent

! O17(P,G)F18 LANDRE ET AL A&A 240, 85 (1990) / nacre
  t9a=t9/(1.d0+2.69d0*t9)
  dt9a=1.d0-2.69d0*t9a
  uno=7.97d+07*t9a**(5.d0/6.d0)/t9**(3.d0/2.d0)*exp(-16.712d0/t9a**(1.d0/3.d0))
  duno=5.d0/6.d0*dt9a-3.d0/2.d0+1.d0/3.d0*16.712d0/t9a**(1.d0/3.d0)*dt9a
  aa=1.51d+08
  bb=16.712d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  due=primus(aa,bb,xe1,xe2,t9)
  ddue=drimus(bb,xe1,xe2,t9)
  aa=1.d0
  bb=0.025d0
  cc=-0.051d0
  dd=-8.82d-03
  ee=0.d0
  fw=0.d0
  xe1=1.d0/3.d0
  xe2=2.d0/3.d0
  xe3=1.d0
  xe4=0.d0
  xe5=0.d0
  tre=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  dtre=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  aa=1.56d+05
  bb=6.272d0
  xe1=1.d0
  xe2=1.d0
  quat=primus(aa,bb,xe1,xe2,t9)
  dquat=drimus(bb,xe1,xe2,t9)
  aa=3.16d-05
  bb=0.767d0
  xe1=3.d0/2.d0
  xe2=1.d0
  cinq=primus(aa,bb,xe1,xe2,t9)
  dcinq=drimus(bb,xe1,xe2,t9)
  aa=98.d0
  bb=2.077d0
  xe1=3.d0/2.d0
  xe2=1.d0
  six=primus(aa,bb,xe1,xe2,t9)
  dsix=drimus(bb,xe1,xe2,t9)
  cent=uno+due*tre+quat+0.5d0*cinq+0.1d0*six
  dcent=uno*duno+due*ddue*tre+due*dtre+0.5d0*cinq*dcinq+0.1d0*six*dsix

  call interpol(35,t9,u)   ! fichier='o17g'
  tcent=u
  b117g(j1)=exp(rh1+f116)*tcent
  ep(11)=yab(1)*yab(9)*q117g*b117g(j1)
  eprt(11)=1.d0+0.5d0*f116
  eptr(11)=-1.5d0*f116+cent/dcent

! O18(P,A)N15 AT. DATA & NUCL. DATA TABLES 40, 283 (1988)/nacre
  aa=3.63d+11
  bb=16.729d0
  cc=1.361d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  xe3=2.d0
  uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
  duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
  aa=1.d0
  bb=0.025d0
  cc=1.88d0
  dd=0.327d0
  ee=4.66d0
  fw=2.06d0
  xe1=1.d0/3.d0
  xe2=2.d0/3.d0
  xe3=1.0d0
  xe4=4.d0/3.d0
  xe5=5.d0/3.d0
  due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  aa=9.90d-14
  bb=0.231d0
  xe1=3.d0/2.d0
  xe2=1.d0
  tre=primus(aa,bb,xe1,xe2,t9)
  dtre=drimus(bb,xe1,xe2,t9)
  aa=2.66d+04
  bb=1.670d0
  xe1=3.d0/2.d0
  xe2=1.d0
  quat=primus(aa,bb,xe1,xe2,t9)
  dquat=drimus(bb,xe1,xe2,t9)
  aa=2.41d+09
  bb=7.638d0
  xe1=3.d0/2.d0
  xe2=1.d0
  cinq=primus(aa,bb,xe1,xe2,t9)
  dcinq=drimus(bb,xe1,xe2,t9)
  aa=1.46d+09
  bb=8.310d0
  xe1=1.d0
  xe2=1.d0
  six=primus(aa,bb,xe1,xe2,t9)
  dsix=drimus(bb,xe1,xe2,t9)
  cent=uno*due+tre+quat+cinq+six
  dcent=uno*duno*due+uno*ddue+tre*dtre+quat*dquat+cinq*dcinq+six*dsix

  call interpol(37,t9,u)   ! fichier='o18a'
  tcent=u
  b118a(j1)=exp(rh1+f116)*tcent
  ep(12)=yab(1)*yab(10)*q118a*b118a(j1)
  eprt(12)=1.d0+0.5d0*f116
  eptr(12)=-1.5d0*f116+dcent/cent

! O18(P,G)F19 AT. DATA & NUCL. DATA TABLES 40, 283 (1988)/nacre
  aa=3.45d+08
  bb=16.729d0
  cc=0.139d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  xe3=2.d0
  uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
  duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
  aa=1.d0
  bb=0.025d0
  cc=2.26d0
  dd=0.394d0
  ee=30.56d0
  fw=13.55d0
  xe1=1.d0/3.d0
  xe2=2.d0/3.d0
  xe3=1.0d0
  xe4=4.d0/3.d0
  xe5=5.d0/3.d0
  due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  aa=1.25d-15
  bb=0.231d0
  xe1=3.d0/2.d0
  xe2=1.d0
  tre=primus(aa,bb,xe1,xe2,t9)
  dtre=drimus(bb,xe1,xe2,t9)
  aa=1.64d+02
  bb=1.670d0
  xe1=3.d0/2.d0
  xe2=1.d0
  quat=primus(aa,bb,xe1,xe2,t9)
  dquat=drimus(bb,xe1,xe2,t9)
  aa=1.28d+04
  bb=5.098d0
  xe1=-1.d0/2.d0
  xe2=1.d0
  cinq=primus(aa,bb,xe1,xe2,t9)
  dcinq=drimus(bb,xe1,xe2,t9)
  cent=uno*due+tre+quat+cinq
  dcent=uno*duno*due+uno*ddue+tre*dtre+quat*dquat+cinq*dcinq

  call interpol(38,t9,u)   ! fichier='o18g'
  tcent=u
  b118g(j1)=exp(rh1+f116)*tcent
  ep(13)=yab(1)*yab(10)*q118g*b118g(j1)
! including Q-value of F19(p,a)O16 when ialflu=0
  if (ialflu == 0) then
    ep(13)=yab(1)*yab(10)*(q118g+8.114d0)*b118g(j1)
  endif
  eprt(13)=1.d0+0.5d0*f116
  eptr(13)=-1.5d0*f116+dcent/cent

! ----------------------------------------------------------------------
!---  RATES FOR NE-NA AND MG-AL CYCLES
  if (ialflu == 1) then

!---  F19(P,G)NE20 (FOWLER 75)/nacre
    b1=1.0d0+0.023d0*t913+2.06d0*t923+0.332d0*t9+3.16d0*t943+1.30d0*t953
    d1=-2.d0/3.d0+6.04d0/t913-11.56d0*t9**2.d0
    d1=d1+(0.0077d0*t913+1.37d0*t923+0.33d0*t9+4.21d0*t943+2.17d0*t953)/b1
    b1=1.d0/t923 * exp(17.92d0-18.11d0/t913-5.78d0*t9**2.)*b1
    if ((6.45d0-3.75d0/t9) < -500.0d0) then
      b2=0.0d0
    else
      b2=1.d0/t932 * exp(6.45d0-3.75d0/t9)
    endif
    if ((11.23d0-5.72d0/t9) < -500.0d0) then
      b3=0.0d0
    else
      b3=1.d0/t9**(2.d0/7.d0) * exp(11.23d0-5.72d0/t9)
    endif
    d2=-1.5d0+3.75d0/t9
    d3=-2.d0/7.d0+5.72d0/t9

    call interpol(11,t9,u)   ! fichier='f19g'
    bg=u
    if (bg /= 0.0d0) then
      eptrg=(b1*d1+b2*d2+b3*d3)/bg
    else
      eptrg=0.0d0
    endif
    bg=bg*exp(rh1+9.d0*fp)
    epg=yab(1)*yab(11)*bg*q119g
    b119g(j1)=bg

!---  F19(P,A)O16 KIOUS 90
    gt9=1.0d0
    if ((-2.090d0/t9) >= 500.0d0) then
      gt9=gt9+4.d0*exp(-2.090d0/t9)
    endif
    if ((-16.44d0/t9) >= 500.0d0) then
      gt9=gt9+7.d0*exp(-16.44d0/t9)
    endif
    if ((-2.090d0/t9) < -500.0d0) then
      dgt9=0.0d0
    else
      dgt9=4.d0*exp(-2.090d0/t9)*2.090d0/t9
    endif
    if ((-16.44d0/t9) >= 500.0d0) then
      dgt9=dgt9+7.d0*exp(-16.44d0/t9)*16.44d0/t9
    endif
    aa=0.776d+13
    bb=21.95d0
    cc=0.845d0
    xe1=2.d0
    xe2=1.d0/3.d0
    xe3=2.d0
    uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
    duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
    due=(1.d0-1.738d0*t923)/(1.d0-0.432d-01/t923)
    ddue=due/(1.d0-0.432d-01/t923)*(-1.738d0*2.d0/3.d0*t923/due-0.432d-01*(2.d0/3.d0)/t923)
    aa=0.136d+12
    bb=21.95d0
    xe1=2.d0/3.d0
    xe2=1.d0/3.d0
    tre=primus(aa,bb,xe1,xe2,t9)
    dtre=drimus(bb,xe1,xe2,t9)
    xious1=uno*due*due+tre
    dkiou1=uno*duno*due*due+uno*2*due*ddue+tre*dtre
    aa=0.2138d+12
    bb=19.40d0
    cc=0.845d0
    xe1=2.d0
    xe2=1.d0/3.d0
    xe3=2.d0
    uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
    duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
    due=(1.d0-4.198d0*t923)/(1.d0-0.3678d-01/t923)
    ddue=due/(1.d0-0.3678d-01/t923)*(-4.198d0*2.d0/3.d0*t923/due-0.3678d-01*(2.d0/3.d0)/t923)
    xious2=uno*due*due
    dkiou2=uno*duno*due*due+uno*2.d0*due*ddue
    aa=9.23d-24
    bb=0.1268d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    uno=primus(aa,bb,xe1,xe2,t9)
    duno=drimus(bb,xe1,xe2,t9)
    aa=888.7d0
    bb=2.47d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    due=primus(aa,bb,xe1,xe2,t9)
    ddue=drimus(bb,xe1,xe2,t9)
    aa=5.83d+06
    bb=3.748d0
    xe1=3.d0/2.d0
    xe2=1.d0
    tre=primus(aa,bb,xe1,xe2,t9)
    dtre=drimus(bb,xe1,xe2,t9)
    aa=3.07d+08
    bb=6.019d0
    xe1=0.d0
    xe2=1.d0
    quat=primus(aa,bb,xe1,xe2,t9)
    dquat=drimus(bb,xe1,xe2,t9)
    xious3=uno+due+tre+quat
    dkiou3=uno*duno+due*ddue+tre*dtre+quat*dquat
    cent=sqrt((xious3+xious1)*(xious3+xious2))/gt9
    dcent=0.5d0/(gt9*gt9*cent)*((dkiou3+dkiou1)*(xious3+xious2)+(xious3+xious1)*(dkiou3+dkiou2))-cent/gt9*dgt9
    eptra=dcent/cent

    call interpol(10,t9,u)   ! fichier='f19a'
    ba=u
    epa=yab(1)*yab(11)*ba*q119a
    ba=ba*exp(rh1+9.d0*fp)
    b119a(j1)=ba
    ep(14)=epa+epg
    if (ep(14) /= 0.d0) then
      eprt(14)=1.d0+4.5d0*fp
      eptr(14)=(epg*eptrg+epa*eptra)/ep(14)-13.5d0*fp
    endif

!---  NE20(P,G)NA21 (FOWLER 75)/nacre
    b1=1.0d0+0.0127d0/t923
    d1=-2.d0+6.482d0/t913-0.017d0/b1**3.d0*t923
    b1=1.d0/t9**2.d0 * exp(16.07d0-19.45d0/t913)/b1**2.d0
    b2=1.d0+2.67d0 *exp(-2.18d0*t912)
    d2=-2.d0/3.d0+6.482d0/t913-2.91d0/b2*exp(-2.18d0*t912)*t912
    b2=1.d0/t923 * exp(19.14d0-19.45d0/t913)*b2
    if ((2.89d0-4.24d0/t9) < -500.0d0) then
      b3=0.0d0
    else
      b3=1.d0/t932 * exp(2.89d0-4.24d0/t9)
    endif
    if ((2.32d0-4.61d0/t9) < -500.0d0) then
      b4=0.0d0
    else
      b4=1.d0/t932 * exp(2.32d0-4.61d0/t9)
    endif
    if ((10.49d0-11.25d0/t9) < -500.0d0) then
      b5=0.0d0
    else
      b5=1.d0/sqrt(t912)*exp(10.49d0-11.25d0/t9)
    endif
    d3=-1.5d0+4.24d0/t9
    d4=-1.5d0+4.61d0/t9
    d5=-0.25d0+11.25d0/t9
    bg=b1+b2+b3+b4+b5
    eptr(15)=(b1*d1+b2*d2+b3*d3+b4*d4+b5*d5)/bg
    eptr(15)=eptr(15)-15.d0*fp
    eprt(15)=1.d0+5.d0*fp

    call interpol(27,t9,u)   ! fichier='ne20'
    bg=u
    bg=bg*exp(rh1+10.d0*fp)
    ep(15)=yab(1)*yab(12)*bg*q120
    b120(j1)=bg

!---  NE21(P,G)NA22 CF88/nacre
    aa=4.37d+08
    bb=19.462d0
    xe1=2.d0/3.d0
    xe2=1.d0/3.d0
    uno=primus(aa,bb,xe1,xe2,t9)
    duno=drimus(bb,xe1,xe2,t9)
    aa=5.85d0
    bb=1.399d0
    xe1=3.d0/2.d0
    xe2=1.d0
    due=primus(aa,bb,xe1,xe2,t9)
    ddue=drimus(bb,xe1,xe2,t9)
    aa=1.29d+04
    bb=3.009d0
    xe1=3.d0/2.d0
    xe2=1.d0
    tre=primus(aa,bb,xe1,xe2,t9)
    dtre=drimus(bb,xe1,xe2,t9)
    aa=3.15d+05
    bb=5.763d0
    xe1=3.d0/5.d0
    xe2=1.d0
    quat=primus(aa,bb,xe1,xe2,t9)
    dquat=drimus(bb,xe1,xe2,t9)
    aa=2.95d+08
    bb=19.462d0
    cc=0.058d0
    xe1=2.d0/3.d0
    xe2=1.d0/3.d0
    xe3=2.d0
    cinq=secun(aa,bb,cc,xe1,xe2,xe3,t9)
    dcinq=dsecun(bb,cc,xe1,xe2,xe3,t9)
    aa=1.d0
    bb=0.021d0
    cc=13.29d0
    dd=1.99d0
    ee=124.1d0
    fw=47.29d0
    xe1=1.d0/3.d0
    xe2=2.d0/3.d0
    xe3=1.d0
    xe4=4.d0/3.d0
    xe5=5.d0/3.d0
    six=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
    dsix=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
    aa=7.80d-01
    bb=1.085d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    sept=primus(aa,bb,xe1,xe2,t9)
    dsept=drimus(bb,xe1,xe2,t9)
    cent=uno+due+tre+quat+0.1d0*(cinq*six+sept)
    dcent=duno*uno+due*ddue+tre*dtre+quat*dquat+0.1d0*(cinq*dcinq*six+cinq*dsix+sept*dsept)
    eptr(16)=dcent/cent
    eptr(16)=eptr(16)-15.d0*fp
    eprt(16)=1.d0+5.d0*fp

    call interpol(28,t9,u)   ! fichier='ne21'
    bg=u
    bg=bg*exp(rh1+10.d0*fp)
    ep(16)=yab(1)*yab(13)*bg*q121
    b121(j1)=bg

!---  NE22(P,G)NA23 (FOWLER 83)/nacre
    b1=1.d0/t923 * exp(20.86d0-19.47d0/t913)
    b2=1.d0/t932 * exp(-25.35d0-0.348d0/t9)
    if ((9.10d0-4.84d0/t9) < -500.0d0) then
      b3=0.0d0
    else
      b3=1.d0/t932 * exp(9.10d0-4.84d0/t9)
    endif
    if ((11.09d0-5.32d0/t9) < -500.0d0) then
      b4=0.0d0
    else
      b4=1.d0/t932 * exp(11.09d0-5.32d0/t9)
    endif
    if ((13.59d0-7.418d0/t9) < -500.0d0) then
      b5=0.0d0
    else
      b5=1.d0/t912 * exp(13.59d0-7.418d0/t9)
    endif
!---  facteur 0.1
    if ((-1.81d0-1.78d0/t9) < -500.0d0) then
      b6=0.0d0
    else
      b6=0.1d0/t932 * exp (-1.81d0-1.78d0/t9)
    endif
    d1=-2.d0/3.d0+6.49d0/t913
    d2=-1.5d0+0.348d0/t9
    d3=-1.5d0+4.84d0/t9
    d4=-1.5d0+5.32d0/t9
    d5=-0.5d0+7.42d0/t9
    d6=-1.5d0+1.78d0/t9
    bg=b1+b2+b3+b4+b5+b6
    eptr(17)=(b1*d1+b2*d2+b3*d3+b4*d4+b5*d5+b6*d6)/bg
    eptr(17)=eptr(17)-15.d0*fp
    eprt(17)=1.d0+5.d0*fp

    call interpol(29,t9,u)   ! fichier='ne2g'
    bg=u
    bg=bg*exp(rh1+10.d0*fp)
    ep(17)=yab(1)*yab(14)*bg*q122
    b122(j1)=bg

!---  NA23(P,A)NE20 (FOWLER 83)/nacre
    b1=1.d0+0.020d0*t913+8.21d0*t923+1.15d0*t9+44.36d0*t943+15.84d0*t953
    d1=-2.d0/3.d0+6.92d0/t913-116.5d0*t9**2.d0
    d1=d1+(0.0067d0*t913+5.47d0*t923+1.15d0*t9+59.15d0*t943+26.4d0*t953)/b1
    b1=1.d0/t923 * exp(22.87d0-20.76d0/t913-58.27d0*t9**2.d0)*b1
    if ((1.38d0-1.99d0/t9) < -500.0d0) then
      b2=0.0d0
    else
      b2=1.d0/t932 * exp(1.38d0-1.99d0/t9)
    endif
    if ((9.38d0-3.15d0/t9) < -500.0d0) then
      b3=0.0d0
    else
      b3=1.d0/t954 * exp(9.38d0-3.15d0/t9)
    endif
    if ((13.66d0-4.37d0/t9) < -500.0d0) then
      b4=0.0d0
    else
      b4=t943 * exp(13.66d0-4.37d0/t9)
    endif
!---  facteur 0.1
    if ((-26.51d0-0.447d0/t9) < -500.0d0) then
      b5=0.0d0
    else
      b5=0.1d0/t932 * exp (-26.51d0-0.447d0/t9)
    endif
    d2=-1.5d0+1.99d0/t9
    d3=-1.25d0+3.15d0/t9
    d4=4.d0/3.d0+4.37d0/t9
    d5=-1.5d0+0.447d0/t9
    ba=b1+b2+b3+b4+b5
    eptra=(b1*d1+b2*d2+b3*d3+b4*d4+b5*d5)/ba

    call interpol(24,t9,u)   ! fichier='na3a'
    ba=u
    ba=exp(rh1+11.d0*fp)*ba
    epa=yab(1)*yab(15)*ba*q123a
    b123a(j1)=ba

!---  NA23(P,G)MG24 (FOWLER 75)/nacre
    b1=1.d0+0.020d0*t913+1.61d0*t923+0.226d0*t9+4.94d0*t943+1.76d0*t953
    d1=-2.d0/3.d0+6.92d0/t913-22.67d0*t9**2.d0
    d1=(0.007d0*t913+1.07d0*t923+0.226d0*t9+6.59d0*t943+2.93d0*t953)/b1+d1
    b1=1.d0/t923 * exp(19.46d0-20.76d0/t913-11.34d0*t9**2.d0)*b1
    if ((4.54d0-2.79d0/t9) < -500.0d0) then
      b2=0.0d0
    else
      b2=1.d0/t932 * exp(4.54d0-2.79d0/t9)
    endif
    if ((9.85d0-3.43d0/t9) < -500.0d0) then
      b3=0.0d0
    else
      b3=1.d0/t932 * exp(9.85d0-3.43d0/t9)
    endif
    if ((10.84d0-5.51d0/t9) < -500.0d0) then
      b4=0.0d0
    else
      b4=t915*exp(10.84d0-5.51d0/t9)
    endif
    d2=-1.5d0+2.79d0/t9
    d3=-1.5d0+3.43d0/t9
    d4=0.2d0+5.51d0/t9
    bg=b1+b2+b3+b4
    eptrg=(b1*d1+b2*d2+b3*d3+b4*d4)/bg

    call interpol(25,t9,u)   ! fichier='na3g'
    bg=u
    bg=exp(rh1+11.d0*fp)*bg
    epg=yab(1)*yab(15)*bg*q123g
    b123g(j1)=bg
    ep(18)=epa+epg
    if (ep(18) /= 0.d0) then
      eprt(18)=1.d0+5.5d0*fp
      eptr(18)=(epg*eptrg+epa*eptra)/ep(18)-16.5d0*fp
    endif

!---  MG24(P,G)AL25 (FOWLER 83)/nacre
    b1=1.d0+0.019d0*t913-0.173d0*t923-0.023d0*t9
    d1=-2.d0/3.d0+7.34d0/t913
    d1=d1+(0.006d0*t913-0.115d0*t923-0.023d0*t9)/b1
    b1=1.d0/t923 * exp(20.14d0-22.02d0/t913)*b1
    if ((7.30d0-2.48d0/t9) < -500.0d0) then
      b2=0.0d0
    else
      b2=1.d0/t932 * exp(7.30d0-2.48d0/t9)
    endif
    if ((8.29d0-4.18d0/t9) < -500.0d0) then
      b3=0.0d0
    else
      b3= exp(8.29d0-4.18d0/t9)
    endif
    d2=-1.5d0+2.48d0/t9
    d3=4.18d0/t9
    bg=b1+b2+b3
    eptr(19)=(b1*d1+b2*d2+b3*d3)/bg-18.d0*fp
    eprt(19)=1.d0+6.d0*fp

    call interpol(14,t9,u)   ! fichier='mg4g'
    bg=u
    bg=exp(rh1+12.d0*fp)*bg
    ep(19)=yab(1)*yab(16)*bg*q124
    b124(j1)=bg

!---  MG25(P,G)AL26G ILLIADIS 90/nacre
    if (t9 <= 0.02d0) then
      cent=1254.7d0*t9-44.526d0
      dcent=1254.7d0*t9
    endif
    if (t9 > 0.02d0 .and. t9 <= 0.04d0) then
     cent=-11765.0d0*t9*t9+1048.75d0*t9-35.701d0
     dcent=cent+35.701d0
    endif
    if (t9 > 0.04d0 .and. t9 <= 0.06d0) then
      cent=-1975.0d0*t9*t9+313.05d0*t9-21.937d0
      dcent=cent+21.937d0
    endif
    if (t9 > 0.06d0 .and. t9 <= 0.08d0) then
      cent=+69.55d0*t9-14.437d0
      dcent=69.55d0*t9
    endif
    if (t9 > 0.08d0 .and. t9 <= 0.10d0) then
      cent=49.2d0*t9-12.809d0
      dcent=49.2d0*t9
    endif
    if (t9 > 0.10d0 .and. t9 <= 0.20d0) then
       cent=49.79d0*t9-12.868d0
       dcent=49.79d0*t9
    endif
    if (t9 > 0.20d0 .and. t9 <= 0.60d0) then
       cent=-28.8875d0*t9*t9+35.0475d0*t9-8.764d0
       dcent=cent+8.764d0
    endif
    if (t9 > 0.60d0 .and. t9 <= 0.80d0) then
      cent=3.235d0*t9-0.076d0
      dcent=3.235d0*t9
    endif
    if (t9 > 0.80d0 .and. t9 <= 1.00d0) then
       cent=2.005d0*t9+0.908d0
       dcent=2.005d0*t9
    endif
    if (t9 < 0.010d0 .or. t9 > 1.00d0) then
       cent=0.0d0
       dcent=0.0d0
    else
       cent=10.d0**cent
    endif
    bg=0.81d0*cent
    eptrg=2.3026d0*dcent

    call interpol(15,t9,u)   ! fichier='mg5g'
    bg=u
    bg=bg*exp(rh1+12.d0*fp)
    epg=yab(1)*yab(17)*bg*q125g
    b125g(j1)=bg

! ---  MG25(P,G)AL26M ILIADIS et al. 90/nacre
    ba=0.19d0*cent
    eptra=2.3026d0*dcent

    call interpol(16,t9,u)   ! fichier='mg5m'
    ba=u
    ba=ba*exp(rh1+12.d0*fp)
    epa=yab(1)*yab(17)*ba*q125m
    b125m(j1)=ba
    ep(20)=epa+epg
    if (ep(20) /= 0.d0) then
      eprt(20)=1.d0+6.d0*fp
      eptr(20)=(epg*eptrg+epa*eptra)/ep(20)-18.d0*fp
    endif

!---  AL26G(P,G)SI27 Champagne 1993/nacre
    aa=1.53d+09
    bb=23.19d0
    xe1=1.75d0
    xe2=1.d0/3.d0
    uno=primus(aa,bb,xe1,xe2,t9)
    duno=drimus(bb,xe1,xe2,t9)
    aa=8.97d0
    bb=2.191d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    due=primus(aa,bb,xe1,xe2,t9)
    ddue=drimus(bb,xe1,xe2,t9)
    aa=473.d0
    bb=3.220d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    tre=primus(aa,bb,xe1,xe2,t9)
    dtre=drimus(bb,xe1,xe2,t9)
    aa=7763.d0
    bb=3.944d0
    xe1=1.0d0
    xe2=1.0d0
    quat=primus(aa,bb,xe1,xe2,t9)
    dquat=drimus(bb,xe1,xe2,t9)
    aa=4.67d-10
    bb=0.789d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    cinqa=primus(aa,bb,xe1,xe2,t9)
    dcinqa=drimus(bb,xe1,xe2,t9)
    aa=3.71e-08
    bb=1.079
    xe1=3./2.
    xe2=1.0
    sixa=primus(aa,bb,xe1,xe2,t9)
    dsixa=drimus(bb,xe1,xe2,t9)
    aa=3.54d-08
    bb=0.789d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    cinqb=primus(aa,bb,xe1,xe2,t9)
    dcinqb=drimus(bb,xe1,xe2,t9)
    aa=3.87d-03
    bb=1.485d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    sixb=primus(aa,bb,xe1,xe2,t9)
    dsixb=drimus(bb,xe1,xe2,t9)
    septa=cinqa+sixa
    dsepta=cinqa*dcinqa+sixa*dsixa
    septb=cinqb+sixb
    dseptb=cinqb*dcinqb+sixb*dsixb
    cent=uno+due+tre+quat+(septa*septb)**0.5d0
    dcent=uno*duno+ddue*due+dtre*tre+dquat*quat+0.5d0*(septa*septb)**0.5d0*(dsepta+dseptb)
    bg=cent
    eptrg=dcent/cent

    call interpol(2,t9,u)   ! fichier='al6g'
    bg=u
    bg=bg*exp(rh1+13.d0*fp)
    epg=yab(1)*yab(19)*bg*q1al26
    b1al26(j1)=bg

! al26g(b+)mg26 (beta+ decay assumes instantaneous annihilation e,e+)
!     tau = 1/2.98e-14 s.  1<q<2.2 we take q=1.7 for mid value
    epa=yab(19)*3.05d-14*2.361d0
    ep(22)=epa+epg
    if (ep(22) /= 0.d0) then
      eprt(22)=epg/ep(22)*(1.0d0+6.5d0*fp)
      eptr(22)=epg/ep(22)*(eptrg-19.5d0*fp)
    endif

!---  MG26(P,G)AL27 ILIADIS ET AL 90/nacre
    if (t9 <= 0.02d0) then
       cent=1208.4d0*t9-45.812d0
       dcent=1208.4d0*t9
    endif
    if (t9 > 0.02d0 .and. t9 <= 0.04d0) then
       cent=-7930.0d0*t9*t9+820.4d0*t9-34.88d0
       dcent=cent+34.88d0
    endif
    if (t9 > 0.04d0 .and. t9 <= 0.06d0) then
       cent=-1505.0d0*t9*t9+321.05d0*t9-25.186d0
       dcent=cent+25.186d0
    endif
    if (t9 > 0.06d0 .and. t9 <= 0.10d0) then
       cent=-901.25d0*t9*t9+239.425d0*t9-22.462d0
       dcent=cent+22.462d0
    endif
    if (t9 > 0.10d0 .and. t9 <= 0.30d0) then
       cent=-137.935d0*t9*t9+92.96d0*t9-15.4487d0
       dcent=cent+15.4487d0
    endif
    if (t9 > 0.30d0 .and. t9 <= 0.70d0) then
       cent=-13.7587d0*t9*t9+20.6505d0*t9-4.9316d0
       dcent=cent+4.9316d0
    endif
    if(t9 > 0.70d0 .and. t9 <= 1.00d0) then
       cent=-2.3d0*t9*t9+5.95d0*t9-0.256d0
       dcent=cent+0.256d0
    endif
    if (t9 < 0.010d0 .or. t9 > 1.00d0) then
       cent=0.0d0
       dcent=0.0d0
    else
       cent=10.d0**cent
    endif
    bg=cent
    eptr(21)=2.3026d0*dcent-18.d0*fp

    call interpol(18,t9,u)   ! fichier='mg6g'
    bg=u
    bg=bg*exp(rh1+12.d0*fp)
    ep(21)=yab(1)*yab(18)*bg*q1mg26
    eprt(21)=1.d0+6.d0*fp
    b1mg26(j1)=bg

! AL27(p,a)Mg24 CHAMPAGNE ET AL. 88/nacre
    aa=1.82d-10
    bb=0.853d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    uno=primus(aa,bb,xe1,xe2,t9)
    duno=drimus(bb,xe1,xe2,t9)
    aa=6.6d-11
    bb=1.d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    due=primus(aa,bb,xe1,xe2,t9)
    ddue=drimus(bb,xe1,xe2,t9)
    aa=2.81d-03
    bb=2.27d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    tre=primus(aa,bb,xe1,xe2,t9)
    dtre=drimus(bb,xe1,xe2,t9)
    aa=2.64d-02
    bb=2.51d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    quat=primus(aa,bb,xe1,xe2,t9)
    dquat=drimus(bb,xe1,xe2,t9)
    aa=0.198d0
    bb=3.3d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    cinq=primus(aa,bb,xe1,xe2,t9)
    dcinq=drimus(bb,xe1,xe2,t9)
    aa=0.528d0
    bb=3.67d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    six=primus(aa,bb,xe1,xe2,t9)
    dsix=drimus(bb,xe1,xe2,t9)
    aa=2.81d0
    bb=4.55d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    sept=primus(aa,bb,xe1,xe2,t9)
    dsept=drimus(bb,xe1,xe2,t9)
    aa=495.d0
    bb=7.34d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    huit=primus(aa,bb,xe1,xe2,t9)
    dhuit=drimus(bb,xe1,xe2,t9)
    cent=uno+due+tre+quat+cinq+six+sept+huit
    dcent=uno*duno+due*ddue+tre*dtre+quat*dquat+cinq*dcinq+six*dsix+sept*dsept+huit*dhuit
    BA=cent
    EPTRA=dcent/cent

    call interpol(3,t9,u)   ! fichier='al7a'
    BA=u
    BA=BA*EXP(RH1+13.d0*FP)
    EPA=YAB(1)*YAB(20)*BA*Q127A
    B127A(J1)=BA

! AL27(p,g)SI28 CF88/nacre
    gt9=1.0d0
    if ((-9.792d0/t9) >= -500.0d0) then
      gt9=gt9+exp(-9.792d0/t9)/3.d0
    endif
    if ((-11.773d0/t9) >= -500.0d0) then
      gt9=gt9+2.d0*exp(-11.773d0/t9)/3.d0
    endif
    if ((-9.792d0/t9) < 500.0d0) then
      dgt9=0.0d0
    else
      dgt9=exp(-9.792d0/t9)/3.d0*9.792d0/t9
    endif
    if ((-11.773d0/t9) >= -500.0d0) then
      dgt9=dgt9+2.d0*exp(-11.773d0/t9)/3.d0*11.773d0/t9
    endif
    aa=1.67d+08
    bb=23.261d0
    cc=0.155d0
    xe1=2.d0/3.d0
    xe2=1.d0/3.d0
    xe3=2.0d0
    uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
    duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
    aa=1.d0
    bb=0.018d0
    cc=5.81d0
    dd=0.728d0
    ee=27.31d0
    fw=8.71d0
    xe1=1.d0/3.d0
    xe2=2.d0/3.d0
    xe3=1.0d0
    xe4=4.d0/3.d0
    xe5=5.d0/3.d0
    due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
    ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
    aa=2.2d0
    bb=2.269d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    tre=primus(aa,bb,xe1,xe2,t9)
    dtre=drimus(bb,xe1,xe2,t9)
    aa=12.2d0
    bb=2.491d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    quat=primus(aa,bb,xe1,xe2,t9)
    dquat=drimus(bb,xe1,xe2,t9)
    aa=1.5d+04
    bb=4.112d0
    xe1=-1.0d0
    xe2=1.0d0
    cinq=primus(aa,bb,xe1,xe2,t9)
    dcinq=drimus(bb,xe1,xe2,t9)
    aa=6.5d-10
    bb=0.853d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    six=primus(aa,bb,xe1,xe2,t9)
    dsix=drimus(bb,xe1,xe2,t9)
    aa=1.63d-10
    bb=1.001d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    sept=primus(aa,bb,xe1,xe2,t9)
    dsept=drimus(bb,xe1,xe2,t9)
    cent=uno*due+tre+quat+cinq+0.1d0*(six+sept)
    dcent=uno*duno*due+uno*ddue+tre*dtre+quat*dquat+cinq*dcinq+0.1d0*(six*dsix+sept*dsept)
    cent=cent/gt9
    dcent=dcent/gt9-cent/(gt9*gt9)*dgt9
    BG=cent
    EPTRG=dcent/cent

    call interpol(4,t9,u)   ! fichier='al7g'
    BG=u
    BG=BG*EXP(RH1+13.d0*FP)
    EPG=YAB(1)*YAB(20)*BG*Q127G
    B127G(J1)=BG
    EP(23)=EPA+EPG
    if (ep(23) /= 0.d0) then
      EPRT(23)=0.5d0+6.5d0*FP
      EPTR(23)=(EPG*EPTRG+EPA*EPTRA)/EP(23)-19.5d0*FP
    endif
  endif      !  IALFLU
! ----------------------------------------------------------------------
! mars 2003: pour fusion H+He, on retire l'ordre d'aller a la fin et on enchaine sur la fusion He (instr.23).
!            Les calculs qui suivent ne seront repris qu'a la fin.
  if (ipop3 == 0) then
    do k=1,13
     eps(j1)=eps(j1)+ep(k)
     epst1=epst1+ep(k)*(eptr(k)+rht1*eprt(k))
     epsp1=epsp1+ep(k)*rhp1*eprt(k)
    enddo
    if (ialflu == 1) then
      do k=14,23
       eps(j1)=eps(j1)+ep(k)
       epst1=epst1+ep(k)*(eptr(k)+rht1*eprt(k))
       epsp1=epsp1+ep(k)*rhp1*eprt(k)
      enddo
    endif
    if (ippcno == 1) then
      do k=1,13
       en13(k,j1)=convMeVAvo*ep(k)
      enddo
      enpp(j1)=convMeVAvo*(ep(1)+ep(2)+ep(3))
      encno(j1)=convMeVAvo*(ep(4)+ep(5)+ep(6)+ep(7)+ep(8)+ep(9)+ep(10)+ep(11)+ep(12)+ep(13))
    endif
    if (eps(j1) /= 0.d0) then
      epst1=epst1/eps(j1)
      epsp1=epsp1/eps(j1)
    else
      epst1=0.d0
      epsp1=0.d0
    endif
    eps(j1)=convMeVAvo*eps(j1)
    en=sqrt(abs(eps(j1)*eps(j)))
    dy=abs(yab(2)*b33(j1)*yab(2)+yab(1)*(-1.5d0*b11(j1)*yab(1)-b112(j1)*yab(4)-b113(j1)*yab(5)-b114(j1)*yab(6)-(b115a(j1)+ &
           b115g(j1))*yab(7)-b116(j1)*yab(8)-(b117a(j1)+b117g(j1))*yab(9)-(b118a(j1)+2.*b118g(j1))*yab(10)))
    if (ialflu == 1) then
      dy=abs(yab(2)*b33(J1)*yab(2)+yab(1)*(-1.5d0*b11(J1)*yab(1)-b112(J1)*yab(4)-b113(J1)*yab(5)-b114(J1)*yab(6)-(b115a(J1)+ &
            b115g(J1))*yab(7)-b116(J1)*yab(8)-(b117a(J1)+b117g(J1))*yab(9)-(b118a(J1)+b118g(J1))*yab(10)-(b119a(J1)+ &
            b119g(J1))*yab(11)-b120(J1)*yab(12)-b121(J1)*yab(13)-b122(J1)*yab(14)-(b123G(J1)+b123a(J1))*yab(15)- &
            b124(J1)*yab(16)-(b125G(J1)+b125m(J1))*yab(17)-b1al26(J1)*yab(19)-b1mg26(J1)*yab(18)-(b127G(J1)+b127a(J1))*yab(20)))
    endif
    if (dy /= 0.d0) then
      zensi(j1)=eps(j1)/(convMeVAvo*dy)
    else
      zensi(j1)=1.0d0
    endif
    return
  endif

23 t8ln=t(j1)-log(1.d8)
  t8=exp(t8ln)
  t9ln=t(j1)-log(1.d9)
  t9=exp(t9ln)
  t912=sqrt(t9)
  t913=exp(t9ln/3.d0)
  t923=t913*t913
  t932=exp(1.5d0*t9ln)
  t943=t923*t923
  t953=t943*t913

! modif PopIII: la temperature n'est testee que pour 22Ne
  if (ipop3 == 0) then
    if (t8ln+0.3567d0) 24,25,25
  endif
25 if (y(j1))27,27,26

! ==================== COMBUSTION HELIUM SEULEMENT ====================
26 continue

  if (t8ln-0.2303d0 > 0.d0 .and. y(j1) < 1.d-7 .or. t8ln-2.3d0 > 0.d0) then
    goto 27
  endif

  call screen (y(j1),xc12(j1),xo16(j1),xne20(j1),xmg24(j1),rh1,rhpsi,rhpsit,rhp1,rht1,vmyo,vmye,zwy(1),zwy(2),zwy(3), &
    zwy(4),zwy(5),t8ln,rhpsip,fy,fyt,fyp,xc13(j1),xn14(j1),xn15(j1),xo17(j1),xo18(j1),xne22(j1),xmg25(j1),xmg26(j1),z)
!--- HE4(2A,G)C12, CAUGHLAN & FOWLER AT. DAT. & NUCL. DAT. TA. 40,283,1988
! V0TO1=0.1
! modif PopIII: le test de temperature suivant est reporte aux calculs pour 8Be(a,g)12C
  if (ipop3 == 1) then
    goto 261
  endif
  if (t9-0.08d0) 261,261,267
!--- HE4(A)BE8
261 aa=7.40d+05
  bb=1.0663d0
  xe1=3.d0/2.d0
  xe2=1.d0
  uno=primus(aa,bb,xe1,xe2,t9)
  duno=drimus(bb,xe1,xe2,t9)
  aa=4.164d+09
  bb=13.490d0
  cc=0.098d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  xe3=2.d0
  due=secun(aa,bb,cc,xe1,xe2,xe3,t9)
  ddue=dsecun(bb,cc,xe1,xe2,xe3,t9)
  aa=1.d0
  bb=0.031d0
  cc=8.009d0
  dd=1.732d0
  ee=49.883d0
  fw=27.426d0
  xe1=1.d0/3.d0
  xe2=2.d0/3.d0
  xe3=1.0d0
  xe4=4.d0/3.d0
  xe5=5.d0/3.d0
  tre=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  dtre=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  ch4ab8=uno+due*tre
  ch4abt=uno*duno+due*ddue*tre+due*dtre
! BE8(A,G)C12
  aa=1.30d+02
  bb=3.3364d0
  xe1=3.d0/2.d0
  xe2=1.d0
  uno=primus(aa,bb,xe1,xe2,t9)
  duno=drimus(bb,xe1,xe2,t9)
  aa=2.510d+07
  bb=23.570d0
  cc=0.235d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  xe3=2.d0
  due=secun(aa,bb,cc,xe1,xe2,xe3,t9)
  ddue=dsecun(bb,cc,xe1,xe2,xe3,t9)
  aa=1.0d0
  bb=0.018d0
  cc=5.249d0
  dd=0.650d0
  ee=19.176d0
  fw=6.034d0
  xe1=1.d0/3.d0
  xe2=2.d0/3.d0
  xe3=1.0d0
  xe4=4.d0/3.d0
  xe5=5.d0/3.d0
  tre=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  dtre=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  cbe8ag=uno+due*tre
  dbe8ag=uno*duno+due*ddue*tre+due*dtre

! modif PopIII: on reintroduit la condition de temperature mise en commentaire precedemment
  if (ipop3 == 0) then
    goto 262
  endif
  if (t9-0.03d0) 262,262,267
262 due=1.d0+4.d0*exp(-(t9/0.025d0)**9.227d0)
  uno=0.01d0+0.2d0*(1.d0+4.d0*exp(-(0.025d0/t9)**3.263d0))/due
  aa=1.35d-07
  bb=24.811d0
  xe1=3.d0/2.d0
  xe2=1.d0
  tre=primus(aa,bb,xe1,xe2,t9)
  dtre=drimus(bb,xe1,xe2,t9)
  ddue=-4.d0*exp(-(t9/0.025d0)**9.227d0)*9.227d0*(t9/0.025d0)**9.227d0
  duno=0.8d0/due*exp(-(0.025d0/t9)**3.263d0)*3.263d0*(0.025d0/t9)**3.263d0-0.2d0*(1.d0+4.d0*exp(-(0.025d0/t9)**3.263d0))/ &
          (due*due)*ddue
  centa=2.90d-16*ch4ab8*cbe8ag*uno
  cent=centa+0.1d0*tre
  dcent=centa*(ch4abt/ch4ab8+dbe8ag/cbe8ag+duno/uno)+0.1d0*tre*dtre

  call interpol(1,t9,u)   ! fichier='aaag'
  tcent=u
  epsyy(j1)=exp(2.d0*rh1+fy)*tcent
  wpsyy=1.828d+16*(y(j1)**3.d0)*epsyy(j1)
  eyt=2.d0*rht1+fy*fyt+dcent/cent
  eyp=2.d0*rhp1+fy*fyp
  goto 266

267 aa=1.35d-07
  bb=24.811d0
  xe1=3.d0/2.d0
  xe2=1.d0
  tre=primus(aa,bb,xe1,xe2,t9)
  dtre=drimus(bb,xe1,xe2,t9)
  aa=2.79d-08
  bb=4.4027d0
  xe1=3.d0
  xe2=1.d0
  uno=primus(aa,bb,xe1,xe2,t9)
  duno=drimus(bb,xe1,xe2,t9)
  cent=uno+0.1d0*tre
  dcent=uno*duno+0.1d0*tre*dtre

  call interpol(1,t9,u)   ! fichier='aaag'
  tcent=u
  epsyy(j1)=exp(2.d0*rh1+fy)*tcent
  wpsyy=1.828d+16*(y(j1)**3.d0)*epsyy(j1)
  eyt=2.d0*rht1+fy*fyt+dcent/cent
  eyp=2.d0*rhp1+fy*fyp

266 n1216=1

!---  C12(a,g)O16, Caughlan et al. 1985/nacre
350 sna=(2.93d+08/(t9*t9))/((1.d0+0.0489d0/t923)**2.d0)*dexp((-32.120d0/t913)-(t9/3.496d0)**2.d0)
  snb=3.14d+08/(t9**2.d0)/((1.d0+0.2654d0/t923)**2.d0)*exp(-32.120d0/t913)
  snc=1.25d+03*dexp(-27.499d0/t9)/t932
  snd=1.43d-02*(t9**5.d0)*dexp(-15.541d0/t9)
  sntot=sna+snb+snc+snd
  call screen (y(j1),xc12(j1),xo16(j1),xne20(j1),xmg24(j1),rh1,rhpsi,rhpsit,rhp1,rht1,vmyo,vmye,zwc(1),zwc(2),zwc(3), &
    zwc(4),zwc(5),t8ln,rhpsip,fyc,fytc,fypc,xc13(j1),xn14(j1),xn15(j1),xo17(j1),xo18(j1),xne22(j1),xmg25(j1),xmg26(j1),z)

  call interpol(6,t9,u)   ! fichier='c12a'
  tsntot=u
  epsyc(j1)=exp(rh1+fyc)*tsntot
  wpsyc=1.440d+17*xc12(j1)*y(j1)*epsyc(j1)
  c12ago(j1)=wpsyc
  ecp=rhp1+fyc*fypc
  dsnd=5.d0+15.541d0/t9
  dsnc=-1.5d0+27.499d0/t9
  dsnb=-2.d0+2.d0/(1.d0+0.2654d0/t923)*2.d0/3.d0*0.2654d0/t923+1.d0/3.d0*32.120d0/t913
  dsna=-2.d0+2.d0/(1.d0+0.0489d0/t923)*2.d0/3.d0*0.0489d0/t923+1.d0/3.d0*32.120d0/t913-2.d0*(t9/3.496d0)**2.d0
  ect=rht1+(sna/sntot)*dsna+(snb/sntot)*dsnb+(snc/sntot)*dsnc+(snd/sntot)*dsnd+fyc*fytc

  goto (351,352),n1216

351 n1620=1

300 continue
  call screen (y(j1),xc12(j1),xo16(j1),xne20(j1),xmg24(j1),rh1,rhpsi,rhpsit,rhp1,rht1,vmyo,vmye,zwo(1),zwo(2),zwo(3), &
    zwo(4),zwo(5),t8ln,rhpsip,fyo,fyto,fypo,xc13(j1),xn14(j1),xn15(j1),xo17(j1),xo18(j1),xne22(j1),xmg25(j1),xmg26(j1),z)

!--- O16(A,G)NE20, CAUGHLAN & FOWLER AT. DAT. & NUCL. DAT. TA. 40,283,1988/nacre
  aa=9.37d+09
  bb=39.757d0
  cc=1.586d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  xe3=2.d0
  uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
  duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
  aa=6.21d+01
  bb=10.297d0
  xe1=3.d0/2.d0
  xe2=1.d0
  due=primus(aa,bb,xe1,xe2,t9)
  ddue=drimus(bb,xe1,xe2,t9)
  aa=5.38d+02
  bb=12.226d0
  xe1=3.d0/2.d0
  xe2=1.d0
  tre=primus(aa,bb,xe1,xe2,t9)
  dtre=drimus(bb,xe1,xe2,t9)
  aa=1.30d+01
  bb=20.093d0
  xe1=-2.d0
  xe2=1.d0
  quat=primus(aa,bb,xe1,xe2,t9)
  dquat=drimus(bb,xe1,xe2,t9)
  cent=uno+due+tre+quat
  dcent=uno*duno+due*ddue+tre*dtre+quat*dquat

  call interpol(32,t9,u)   ! fichier='o16a'
  tcent=u
  epsyo(j1)=exp(rh1+fyo)*tcent
  wpsyo=7.137d+16*xo16(j1)*y(j1)*epsyo(j1)
  eot=rht1+fyo*fyto+dcent/cent
  eop=rhp1+fyo*fypo
  o16agn(j1)=wpsyo

  goto (301,302),n1620

301 continue

!--- O17(A,N)NE20  MEME SCREENING QUE O16(A,G)NE20/nacre
  t9a=t9/(1.d0+0.0268d0*t9+0.0232d0*t923*t923*t913/(1.d0+0.0268d0*t9)**(2.d0/3.d0))
  gt9=1.d0+dexp(-10.106d0/t9)/3.d0
  cent=1.03d+18/gt9*(t9a**(5.d0/6.d0))/t932*dexp(-39.914d0/(t9a**(1.d0/3.d0)))
  dt9a=1.d0-t9a*(0.0268d0+0.0232d0*t923*(5.d0+3.d0*0.0268d0*t9)/(3.d0*(1.d0+0.0268d0*t9)**(5.d0/3.d0)))
  dgt9=10.106d0*dexp(-10.106d0/t9)/(3.d0*t9*gt9)
  dcent=-dgt9+5.d0*dt9a/6.d0-3.d0/2.d0+39.914d0*dt9a/(3.d0*(t9a**(1.d0/3.d0)))

  call interpol(36,t9,u)   ! fichier='o17n'
  tcent=u
  e17an(j1)=exp(rh1+fyo)*tcent
  w17an=8.372d+15*xo17(j1)*y(j1)*e17an(j1)
  e17ant=rht1+fyo*fyto+dcent
  e17anp=rhp1+fyo*fypo
  call screen(y(j1),xc12(j1),xo16(j1),xne20(j1),xmg24(j1),rh1,rhpsi,rhpsit,rhp1,rht1,vmyo,vmye,zw14(1),zw14(2),zw14(3), &
    zw14(4),zw14(5),t8ln,rhpsip,fyn14,fyt14,fyp14,xc13(j1),xn14(j1),xn15(j1),xo17(j1),xo18(j1),xne22(j1),xmg25(j1),xmg26(j1),z)

!--- N14(A,G)F18, CAUGHLAN & FOWLER AT. DAT. & NUCL. DAT. TA. 40,283,1988/nacre
  aa=7.78d+09
  bb=36.031d0
  cc=0.881d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  xe3=2.d0
  uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
  duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
  aa=1.d0
  bb=0.012d0
  cc=1.45d0
  dd=0.117d0
  ee=1.97d0
  fw=0.406d0
  xe1=1.d0/3.d0
  xe2=2.d0/3.d0
  xe3=1.0d0
  xe4=4.d0/3.d0
  xe5=5.d0/3.d0
  due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  aa=2.36d-10
  bb=2.798d0
  xe1=3.d0/2.d0
  xe2=1.d0
  tre=primus(aa,bb,xe1,xe2,t9)
  dtre=drimus(bb,xe1,xe2,t9)
  aa=2.03d0
  bb=5.054d0
  xe1=3.d0/2.d0
  xe2=1.d0
  quat=primus(aa,bb,xe1,xe2,t9)
  dquat=drimus(bb,xe1,xe2,t9)
  aa=1.15d+04
  bb=12.310d0
  xe1=2.d0/3.d0
  xe2=1.d0
  cinq=primus(aa,bb,xe1,xe2,t9)
  dcinq=drimus(bb,xe1,xe2,t9)
  cent=uno*due+tre+quat+cinq
  dcent=uno*duno*due+uno*ddue+tre*dtre+quat*dquat+cinq*dcinq

  call interpol(19,t9,u)   ! fichier='n14a'
  tcent=u
  c144(j1)=exp(rh1+fyn14)*tcent
  if(ialflu == 1) then
    e144=7.607d+16*xn14(j1)*y(j1)*c144(j1)
  else
    e144=9.774d+16*xn14(j1)*y(j1)*c144(j1)
  endif
  e144t=rht1+fyn14*fyt14+dcent/cent
  e144p=rhp1+fyn14*fyp14

!     MEME SCREENING: FYO
!--- O18(A,G)NE22, Giessen et al. 1994, Nucl. Phys. A567, 146/nacre
  if (t9 >= 0.06d0 .and. t9 <= 0.08d0) then
    cent=177.7285d0*t9-36.51758d0
    dcent=177.7285d0*t9
  endif
  if (t9 > 0.08d0 .and. t9 <= 0.12d0) then
    cent=142.033d0*t9-33.66194d0
    dcent=142.033d0
  endif
  if (t9 > 0.12d0 .and. t9 <= 0.16d0) then
    cent=-700.13750006d0*t9*t9+292.37625002d0*t9-41.62115d0
    dcent=2.d0*(-700.13750006d0)*t9*t9+292.37625002d0*t9
  endif
  if (t9 > 0.16d0 .and. t9 <= 0.20d0) then
    cent=-317.75d0*t9*t9+171.8525d0*t9-32.12647d0
    dcent=2.d0*(-317.75d0)*t9*t9+171.8525d0*t9
  endif
  if (t9 > 0.20d0 .and. t9 <= 0.40d0) then
    cent=-65.0205d0*t9*t9+65.48225d0*t9-20.9616d0
    dcent=2.d0*(-65.0205d0)*t9*t9+65.48225d0*t9
  endif
  if (t9 > 0.40d0 .and. t9 <= 1.00d0) then
    cent=-10.82333d0*t9*t9+21.527983d0*t9-12.05144d0
    dcent=2.d0*(-10.82333d0)*t9*t9+21.527983d0*t9
  endif
  if (t9 < 0.06d0 .or. t9 > 1.00d0) then
    cent=0.0d0
    dcent=0.0d0
    if (ipop3 /= 1 .and. verbose) then
      write(*,*) ' temp. hors des limites o18(a,g)'
    endif
  else
    cent=10.d0**cent
  endif

  call interpol(40,t9,u)   ! fichier='o8ag'
  tcent=u
  c184(j1)=exp(rh1+fyo)*tcent
  e184=1.296d+17*y(j1)*xo18(j1)*c184(j1)
  e184t=dcent*2.3026d0+rht1+fyo*fyto
  e184p=rhp1+fypo*fyo
  call screen(y(j1),xc12(j1),xo16(j1),xne20(j1),xmg24(j1),rh1,rhpsi,rhpsit,rhp1,rht1,vmyo,vmye,zw20(1),zw20(2),zw20(3), &
    zw20(4),zw20(5),t8ln,rhpsip,fyn22,fyt22,fyp22,xc13(j1),xn14(j1),xn15(j1),xo17(j1),xo18(j1),xne22(j1),xmg25(j1),xmg26(j1),z)

! modif PopIII: si la temperature n'est pas suffisante, on saute les reactions du 22Ne
  if (ipop3 == 1) then
    e224=0.d0
    e224t=0.d0
    e224p=0.d0
    e224g=0.d0
    e22tg=0.d0
  endif
  if ((t8ln+0.3567d0) > 0.d0 .or. ipop3 == 0) then
! NE22(A,N)MG25, CAUGHLAN, FOWLER, AT. DAT. & NUCL. DAT. TA. 40,283,88/nacre
    gt9=1.d0+5.0d0*exp(-14.791d0/t9)
    t9a=t9/(1.d0+0.0548d0*t9)
    ft9a=exp(-(0.197d0/t9a)**4.82d0)
    uno=4.16d+19*ft9a/gt9*(t9a**(5.d0/6.d0))/t932*exp(-47.004d0/(t9a**(1.d0/3.d0)))
    due=1.44d-04/gt9*exp(-5.577d0/t9)
    cent=uno+due
    dgt9=5.d0*exp(-14.791d0/t9)*14.791d0/t9
    dt9a=t9a*(1.d0-0.0548d0*t9a)
    dft9a=ft9a*4.82d0*(0.197d0/t9a)**4.82d0*dt9a/t9a
    duno=dft9a/ft9a-dgt9/gt9+5.d0/6.d0*dt9a/t9a-3.d0/2.d0+1.d0/3.d0*47.004d0/(t9a**(4.d0/3.d0))*dt9a
    ddue=-dgt9/gt9+5.577d0/t9
    dcent=uno*duno+due*ddue

    call interpol(30,t9,u)   ! fichier='ne2n'
    tcent=u
    c224(j1)=exp(rh1+fyn22)*tcent
    e224=-5.274d+15*xne22(j1)*y(j1)*c224(j1)
    e224t=rht1+fyt22*fyn22+dcent/cent
    e224p=rhp1+fyp22*fyn22

! NE22(A,G)MG26, CAUGHLAN, FOWLER, AT. DAT. & NUCL. DAT. TA. 40,283,88/nacre
    gt9=1.d0+5.0d0*exp(-14.791d0/t9)
    t9a=t9/(1.d0+0.0548d0*t9)
    ft9a=exp(-(0.197d0/t9a)**4.82d0)
    fpt9a=exp(-(t9a/0.249d0)**2.31d0)
    uno=4.16d+19*fpt9a/gt9*(t9a**(5.d0/6.d0))/t932*exp(-47.004d0/(t9a**(1.d0/3.d0)))
    due=2.08d+16*ft9a/gt9*(t9a**(5.d0/6.d0))/t932*exp(-47.004d0/(t9a**(1.d0/3.d0)))
    cent=uno+due
    dgt9=5.d0*exp(-14.791d0/t9)*14.791d0/t9
    dt9a=t9a*(1.d0-0.0548d0*t9a)
    dft9a=ft9a*4.82d0*(0.197d0/t9a)**4.82d0*dt9a/t9a
    dfpt9a=-fpt9a*2.31d0*((t9a/0.249)**1.31d0)*dt9a/0.249
    duno=dfpt9a/fpt9a-dgt9/gt9+5.d0/6.d0*dt9a/t9a-3.d0/2.d0+1.d0/3.d0*47.004d0/(t9a**(4.d0/3.d0))*dt9a
    ddue=dft9a/ft9a-dgt9/gt9+5.d0/6.d0*dt9a/t9a-3.d0/2.d0+1.d0/3.d0*47.004d0/(t9a**(4.d0/3.d0))*dt9a
    dcent=uno*duno+due*ddue

    call interpol(31,t9,u)   ! fichier='neag'
    tcent=u
    c224g(j1)=exp(rh1+fyn22)*tcent
    e224g=10.612d0*convMeVAvo/(22.d0*4.d0)*xne22(j1)*y(j1)*c224g(j1)
    e22tg=rht1+fyt22*fyn22+dcent/cent
  endif

!--- NE20(A,G)MG24,CAUGHLAN & FOWLER AT. DAT. & NUCL. DAT. TA. 40,283,1988
! V0TO1 0.1/nacre
  aa=4.11d+11
  bb=46.766d0
  cc=2.219d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  xe3=2.d0
  uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
  duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
  aa=1.d0
  bb=0.009d0
  cc=0.882d0
  dd=0.055d0
  ee=0.749d0
  fw=0.119d0
  xe1=1.d0/3.d0
  xe2=2.d0/3.d0
  xe3=1.0d0
  xe4=4.d0/3.d0
  xe5=5.d0/3.d0
  due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  aa=5.27d+03
  bb=15.869d0
  xe1=3.d0/2.d0
  xe2=1.0d0
  tre=primus(aa,bb,xe1,xe2,t9)
  dtre=drimus(bb,xe1,xe2,t9)
  aa=6.51d+03
  bb=16.223d0
  xe1=-1.d0/2.d0
  xe2=1.0d0
  quat=primus(aa,bb,xe1,xe2,t9)
  dquat=drimus(bb,xe1,xe2,t9)
  aa=4.21d+01
  bb=9.115d0
  xe1=3.d0/2.d0
  xe2=1.0d0
  cinq=primus(aa,bb,xe1,xe2,t9)
  dcinq=drimus(bb,xe1,xe2,t9)
  aa=3.2d+01
  bb=9.383d0
  xe1=2.d0/3.d0
  xe2=1.0d0
  six=primus(aa,bb,xe1,xe2,t9)
  dsix=drimus(bb,xe1,xe2,t9)
  gt9=1.d0+5.d0*exp(-18.960d0/t9)
  dgt9=5.d0*exp(-18.960d0/t9)*18.960d0/t9
  centa=uno*due+tre+quat+0.1d0*(cinq+six)
  cent=centa/gt9
  dcent=(uno*duno*due+uno*ddue+tre*dtre+quat*dquat+0.1d0*(cinq*dcinq+six*dsix))/gt9-centa*dgt9/(gt9*gt9)

  call interpol(26,t9,u)   ! fichier='ne0g'
  tcent=u
  e20ag(j1)=exp(rh1+fyn22)*tcent
  w20ag=9.312d0*convMeVAvo/(20.d0*4.d0)*xne20(j1)*y(j1)*e20ag(j1)
  e20agt=rht1+fyn22*fyt22+dcent/cent
  e20agp=rhp1+fyn22*fyp22

!--- C13(A,N)O16, CAUGHLAN & FOWLER AT. DAT. & NUCL. DAT. TA. 40,283,1988/nacre
  aa=6.77d+15
  bb=32.329d0
  cc=1.284d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  xe3=2.d0
  uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
  duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
  aa=1.d0
  bb=0.013d0
  cc=2.04d0
  dd=0.184d0
  ee=0.0d0
  fw=0.0d0
  xe1=1.d0/3.d0
  xe2=2.d0/3.d0
  xe3=1.0d0
  xe4=0.0d0
  xe5=0.0d0
  due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  aa=3.82d+05
  bb=9.373d0
  xe1=3.d0/2.d0
  xe2=1.0d0
  tre=primus(aa,bb,xe1,xe2,t9)
  dtre=drimus(bb,xe1,xe2,t9)
  aa=1.41d+06
  bb=11.873d0
  xe1=3.d0/2.d0
  xe2=1.d0
  quat=primus(aa,bb,xe1,xe2,t9)
  dquat=drimus(bb,xe1,xe2,t9)
  aa=2.00d+09
  bb=20.409d0
  xe1=3.d0/2.d0
  xe2=1.0d0
  cinq=primus(aa,bb,xe1,xe2,t9)
  dcinq=drimus(bb,xe1,xe2,t9)
  aa=2.92d+09
  bb=29.283d0
  xe1=3.d0/2.d0
  xe2=1.0d0
  six=primus(aa,bb,xe1,xe2,t9)
  dsix=drimus(bb,xe1,xe2,t9)
  cent=uno*due+tre+quat+cinq+six
  dcent=uno*duno*due+uno*ddue+tre*dtre+quat*dquat+cinq*dcinq+six*dsix

  call interpol(9,t9,u)   ! fichier='c13n'
  tcent=u
  c134(j1)=exp(rh1+fyc)*tcent
  e134=4.112d+16*xc13(j1)*y(j1)*c134(j1)
  e13t=rht1+fyc*fytc+dcent/cent
  e13p=rhp1+fyc*fypc

  if(ialflu == 1) then
! population III: si on est encore dans la fusion centrale H, on saute les reactions qui ont deja ete calculees
!                 ainsi que les reactions de capture de neutrons
    if (j1 == m .and. verbose) then
      write(3,'(a)')'ENERG: cycle de l''alu'
    endif
    if (ipop3 == 0 .or. x(j1) <= 1.0d-7) then
! O18(P,A)N15 AT. DATA & NUCL. DATA TABLES 40, 283 (1988)/nacre
      call screen(y(j1),xc12(j1),xo16(j1),xne20(j1),xmg24(j1),rh1,rhpsi,rhpsit,rhp1,rht1,vmyo,vmye,zw18(1),zw18(2),zw18(3), &
        zw18(4),zw18(5),t8ln,rhpsip,fop,fopt,fopp,xc13(j1),xn14(j1),xn15(j1),xo17(j1),xo18(j1),xne22(j1),xmg25(j1),xmg26(j1),z)
      aa=3.63d+11
      bb=16.729d0
      cc=1.361d0
      xe1=2.d0/3.d0
      xe2=1.d0/3.d0
      xe3=2.d0
      uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
      duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
      aa=1.d0
      bb=0.025d0
      cc=1.88d0
      dd=0.327d0
      ee=4.66d0
      fw=2.06d0
      xe1=1.d0/3.d0
      xe2=2.d0/3.d0
      xe3=1.0d0
      xe4=4.d0/3.d0
      xe5=5.d0/3.d0
      due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
      ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
      aa=9.90d-14
      bb=0.231d0
      xe1=3.d0/2.d0
      xe2=1.d0
      tre=primus(aa,bb,xe1,xe2,t9)
      dtre=drimus(bb,xe1,xe2,t9)
      aa=2.66d+04
      bb=1.670d0
      xe1=3.d0/2.d0
      xe2=1.d0
      quat=primus(aa,bb,xe1,xe2,t9)
      dquat=drimus(bb,xe1,xe2,t9)
      aa=2.41d+09
      bb=7.638d0
      xe1=3.d0/2.d0
      xe2=1.d0
      cinq=primus(aa,bb,xe1,xe2,t9)
      dcinq=drimus(bb,xe1,xe2,t9)
      aa=1.46d+09
      bb=8.310d0
      xe1=1.d0
      xe2=1.d0
      six=primus(aa,bb,xe1,xe2,t9)
      dsix=drimus(bb,xe1,xe2,t9)
      cent=uno*due+tre+quat+cinq+six
      dcent=uno*duno*due+uno*ddue+tre*dtre+quat*dquat+cinq*dcinq+six*dsix

      call interpol(37,t9,u)   ! fichier='o18a'
      tcent=u
      e18pa(j1)=exp(rh1+fop)*tcent
      w18pa=3.980d0*convMeVAvo/18.d0*xf18(j1)*xprot(j1)*e18pa(j1)
      e18pat=rht1+fop*fopt+dcent/cent
      e18pap=rhp1+fop*fopp
      if(j1 == m .and. verbose) then
        write(3,'(1x,a,4(1x,e24.18))') '18O(p,a):  ',e18pa(j1),w18pa,e18pat,e18pap
      endif
    endif

!--- C14(p,g)N15, WGT90
    call screen(y(j1),xc12(j1),xo16(j1),xne20(j1),xmg24(j1),rh1,rhpsi,rhpsit,rhp1,rht1,vmyo,vmye,zwcp(1),zwcp(2),zwcp(3), &
      zwcp(4),zwcp(5),t8ln,rhpsip,fcp,fcpt,fcpp,xc13(j1),xn14(j1),xn15(j1),xo17(j1),xo18(j1),xne22(j1),xmg25(j1),xmg26(j1),z)
    aa=1.09d+8
    bb=13.71d0
    cc=4.694d0
    xe1=2.d0/3.d0
    xe2=1.d0/3.d0
    xe3=2.d0
    uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
    duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
    aa=1.d0
    bb=0.0304d0
    cc=0.105d0
    dd=2.24d-02
    ee=0.109d0
    fw=5.94d-02
    xe1=1.d0/3.d0
    xe2=2.d0/3.d0
    xe3=1.d0
    xe4=4.d0/3.d0
    xe5=5.d0/3.d0
    due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
    ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
    aa=49.53d0
    bb=2.832d0
    xe1=3.d0/2.d0
    xe2=1.d0
    tre=primus(aa,bb,xe1,xe2,t9)
    dtre=drimus(bb,xe1,xe2,t9)
    aa=6.32d+3
    bb=3.795d0
    xe1=3.d0/2.d0
    xe2=1.d0
    quat=primus(aa,bb,xe1,xe2,t9)
    dquat=drimus(bb,xe1,xe2,t9)
    aa=529.5d0
    bb=5.64d0
    xe1=3.d0/2.d0
    xe2=1.d0
    cinq=primus(aa,bb,xe1,xe2,t9)
    dcinq=drimus(bb,xe1,xe2,t9)
    aa=1.435d+05
    bb=5.721d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    six=primus(aa,bb,xe1,xe2,t9)
    dsix=drimus(bb,xe1,xe2,t9)
    aa=4.612d+4
    bb=6.882d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    sept=primus(aa,bb,xe1,xe2,t9)
    dsept=drimus(bb,xe1,xe2,t9)
    cent=uno*due+tre+quat+cinq+six+sept
    dcent=duno*uno*due+uno*ddue+tre*dtre+quat*dquat+cinq*dcinq+six*dsix+sept*dsept

    call interpol(42,t9,u)   ! fichier='c14g'
    tcent=u
    ec14pg(j1)=exp(rh1+fcp)*tcent
    wc14pg=10.207d0*convMeVAvo/14.d0*xc14(j1)*xprot(j1)*ec14pg(j1)
    ec14pt=rht1+fcp*fcpt+dcent/cent
    ec14pp=rhp1+fcp*fcpp
    if (j1 == m .and. verbose) then
      write(3,'(1x,a,4(1x,e24.18))') '14C(p,g):  ',ec14pg(j1),wc14pg,ec14pt,ec14pp
    endif
    if (ipop3 == 0 .or. x(j1) <= 1.0d-7) then
! C12(P,G)N13 AT. DATA & NUCL. DATA TABLES 40, 283 (1988)*******
      aa=2.04d+07
      bb=13.690d0
      cc=1.500d0
      xe1=2.d0/3.d0
      xe2=1.d0/3.d0
      xe3=2.d0
      uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
      duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
      aa=1.d0
      bb=0.030d0
      cc=1.19d0
      dd=0.254d0
      ee=2.06d0
      fw=1.12d0
      xe1=1.d0/3.d0
      xe2=2.d0/3.d0
      xe3=1.0d0
      xe4=4.d0/3.d0
      xe5=5.d0/3.d0
      due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
      ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
      aa=1.08d+05
      bb=4.925d0
      xe1=3.d0/2.d0
      xe2=1.d0
      tre=primus(aa,bb,xe1,xe2,t9)
      dtre=drimus(bb,xe1,xe2,t9)
      aa=2.15d+05
      bb=18.179d0
      xe1=3.d0/2.d0
      xe2=1.d0
      quat=primus(aa,bb,xe1,xe2,t9)
      dquat=drimus(bb,xe1,xe2,t9)
      cent=uno*due+tre+quat
      dcent=uno*duno*due+uno*ddue+tre*dtre+quat*dquat

      call interpol(7,t9,u)   ! fichier='c12g'
      tcent=u
      b112(j1)=exp(rh1+fcp)*tcent
      ep(4)=q112*b112(j1)*convMeVAvo/12.d0*xc12(j1)*xprot(j1)
      eprt(4)=rhp1+fcp*fcpp
      eptr(4)=rht1+fcp*fcpt+dcent/cent

!--- C13(P,G)N14 AT. DATA & NUCL. DATA TABLES 40, 283 (1988)/nacre
      aa=8.01d+07
      bb=13.717d0
      cc=2.d0
      xe1=2.d0/3.d0
      xe2=1.d0/3.d0
      xe3=2.d0
      uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
      duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
      aa=1.d0
      bb=0.030d0
      cc=0.958d0
      dd=0.204d0
      ee=1.39d0
      fw=0.753d0
      xe1=1.d0/3.d0
      xe2=2.d0/3.d0
      xe3=1.0d0
      xe4=4.d0/3.d0
      xe5=5.d0/3.d0
      due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
      ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
      aa=1.21d+06
      bb=5.701d0
      xe1=6.d0/5.d0
      xe2=1.d0
      tre=primus(aa,bb,xe1,xe2,t9)
      dtre=drimus(bb,xe1,xe2,t9)
      cent=uno*due+tre
      dcent=uno*duno*due+uno*ddue+tre*dtre

      call interpol(8,t9,u)   ! fichier='c13g'
      tcent=u
      b113(j1)=exp(rh1+fcp)*tcent
      ep(5)=q113*b113(j1)*convMeVAvo/13.d0*xc13(j1)*xprot(j1)
      eprt(5)=rhp1+fcp*fcpp
      eptr(5)=rht1+fcp*fcpt+dcent/cent
! N14(P,G)O15 AT. DATA & NUCL. DATA TABLES 40, 283 (1988)/nacre
      call screen(y(j1),xc12(j1),xo16(j1),xne20(j1),xmg24(j1),rh1,rhpsi,rhpsit,rhp1,rht1,vmyo,vmye,zwn(1),zwn(2),zwn(3), &
        zwn(4),zwn(5),t8ln,rhpsip,fnp,fnpt,fnpp,xc13(j1),xn14(j1),xn15(j1),xo17(j1),xo18(j1),xne22(j1),xmg25(j1),xmg26(j1),z)
      aa=4.90d+07
      bb=15.228d0
      cc=3.294d0
      xe1=2.d0/3.d0
      xe2=1.d0/3.d0
      xe3=2.d0
      uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
      duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
      aa=1.d0
      bb=0.027d0
      cc=-0.778d0
      dd=-0.149d0
      ee=0.261d0
      fw=0.127d0
      xe1=1.d0/3.d0
      xe2=2.d0/3.d0
      xe3=1.0d0
      xe4=4.d0/3.d0
      xe5=5.d0/3.d0
      due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
      ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
      aa=2.37d+03
      bb=3.011d0
      xe1=3.d0/2.d0
      xe2=1.d0
      tre=primus(aa,bb,xe1,xe2,t9)
      dtre=drimus(bb,xe1,xe2,t9)
      aa=2.19d+04
      bb=12.530d0
      xe1=0.d0
      xe2=1.d0
      quat=primus(aa,bb,xe1,xe2,t9)
      dquat=drimus(bb,xe1,xe2,t9)
      cent=uno*due+tre+quat
      dcent=uno*duno*due+uno*ddue+tre*dtre+quat*dquat

      call interpol(20,t9,u)   ! fichier='n14g'
      tcent=u
      b114(j1)=exp(rh1+fnp)*tcent
      ep(6)=q114*b114(j1)*convMeVAvo/14.d0*xn14(j1)*xprot(j1)
      eprt(6)=rhp1+fnp*fnpp
      eptr(6)=dcent/cent+rht1+fnp*fnpt

!--- N15(P,G)O16 AT. DATA & NUCL. DATA TABLES 40, 283 (1988)/nacre
      aa=9.78d+08
      bb=15.251d0
      cc=0.450d0
      xe1=2.d0/3.d0
      xe2=1.d0/3.d0
      xe3=2.d0
      uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
      duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
      aa=1.d0
      bb=0.027d0
      cc=0.219d0
      dd=0.042d0
      ee=6.83d0
      fw=3.32d0
      xe1=1.d0/3.d0
      xe2=2.d0/3.d0
      xe3=1.0d0
      xe4=4.d0/3.d0
      xe5=5.d0/3.d0
      due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
      ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
      aa=1.11d+04
      bb=3.328d0
      xe1=3.d0/2.d0
      xe2=1.d0
      tre=primus(aa,bb,xe1,xe2,t9)
      dtre=drimus(bb,xe1,xe2,t9)
      aa=1.49d+04
      bb=4.665d0
      xe1=3.d0/2.d0
      xe2=1.d0
      quat=primus(aa,bb,xe1,xe2,t9)
      dquat=drimus(bb,xe1,xe2,t9)
      aa=3.80d+06
      bb=11.048d0
      xe1=3.d0/2.d0
      xe2=1.d0
      cinq=primus(aa,bb,xe1,xe2,t9)
      dcinq=drimus(bb,xe1,xe2,t9)
      cent=uno*due+tre+quat+cinq
      dcent=uno*duno*due+uno*ddue+tre*dtre+quat*dquat+cinq*dcinq

      call interpol(22,t9,u)   ! fichier='n15g'
      tcent=u
      b115g(j1)=exp(rh1+fnp)*tcent
      ep(7)=q115g*b115g(j1)*convMeVAvo/15.d0*xn15(j1)*xprot(j1)
      eprt(7)=rhp1+fnp*fnpp
      eptr(7)=dcent/cent+rht1+fnp*fnpt

!--- N15(P,A)C12 AT. DATA & NUCL. DATA TABLES 40, 283 (1988)/nacre
      aa=1.08d+12
      bb=15.251d0
      cc=0.522d0
      xe1=2.d0/3.d0
      xe2=1.d0/3.d0
      xe3=2.d0
      uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
      duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
      aa=1.d0
      bb=0.027d0
      cc=2.62d0
      dd=0.501d0
      ee=5.36d0
      fw=2.60d0
      xe1=1.d0/3.d0
      xe2=2.d0/3.d0
      xe3=1.0d0
      xe4=4.d0/3.d0
      xe5=5.d0/3.d0
      due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
      ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
      aa=1.19d+08
      bb=3.676d0
      xe1=3.d0/2.d0
      xe2=1.d0
      tre=primus(aa,bb,xe1,xe2,t9)
      dtre=drimus(bb,xe1,xe2,t9)
      aa=5.41d+08
      bb=8.926d0
      xe1=1.d0/2.d0
      xe2=1.d0
      quat=primus(aa,bb,xe1,xe2,t9)
      dquat=drimus(bb,xe1,xe2,t9)
      aa=4.72d+08
      bb=7.721d0
      xe1=3.d0/2.d0
      xe2=1.d0
      cinq=primus(aa,bb,xe1,xe2,t9)
      dcinq=drimus(bb,xe1,xe2,t9)
      aa=2.20d+09
      bb=11.418d0
      xe1=3.d0/2.d0
      xe2=1.d0
      six=primus(aa,bb,xe1,xe2,t9)
      dsix=drimus(bb,xe1,xe2,t9)
      cent=uno*due+tre+quat+0.1d0*(cinq+six)
      dcent=uno*duno*due+uno*ddue+tre*dtre+quat*dquat+0.1d0*(cinq*dcinq+six*dsix)

      call interpol(21,t9,u)   ! fichier='n15a'
      tcent=u
      b115a(j1)=exp(rh1+fnp)*tcent
      ep(8)=q115a*b115a(j1)*convMeVAvo/15.d0*xn15(j1)*xprot(j1)
      eprt(8)=rhp1+fnp*fnpp
      eptr(8)=dcent/cent+rht1+fnp*fnpt

!--- O16(P,G)F17 nacre
      aa=1.5d+08
      bb=16.692d0
      xe1=2.d0/3.d0
      xe2=1.d0/3.d0
      uno=primus(aa,bb,xe1,xe2,t9)
      duno=drimus(bb,xe1,xe2,t9)
      due=1.d0+2.13d0*(1.d0-exp(-0.728d0*t9**(2.d0/3.d0)))
      cent=uno/due
      ddue=2.13d0*0.728d0*2.d0/3.d0*t9**(2.d0/3.d0)*exp(-0.728d0*t9**(2.d0/3.d0))/due
      dcent=duno-ddue

      call interpol(33,t9,u)   ! fichier='o16g'
      tcent=u
      b116(j1)=exp(rh1+fop)*tcent
      ep(9)=q116*b116(j1)*convMeVAvo/16.d0*xo16(j1)*xprot(j1)
      eprt(9)=rhp1+fop*fopp
      eptr(9)=dcent+rht1+fop*fopt

!--- O17(P,A)N14 AT. LANDRE ET AL. A&A 240 85 (1990)/nacre
      aa=1.53d+07
      bb=16.712d0
      cc=0.565d0
      xe1=2.d0/3.d0
      xe2=1.d0/3.d0
      xe3=2.d0
      uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
      duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
      aa=1.d0
      bb=0.025d0
      cc=5.39d0
      dd=0.940d0
      ee=13.5d0
      fw=5.98d0
      xe1=1.d0/3.d0
      xe2=2.d0/3.d0
      xe3=1.0d0
      xe4=4.d0/3.d0
      xe5=5.d0/3.d0
      due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
      ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
      aa=2.92d+06
      bb=4.247d0
      xe1=-1.d0
      xe2=1.d0
      tre=primus(aa,bb,xe1,xe2,t9)
      dtre=drimus(bb,xe1,xe2,t9)
      aa=1.78d+05
      bb=16.67d0
      xe1=2.d0/3.d0
      xe2=1.d0/3.d0
      quat=primus(aa,bb,xe1,xe2,t9)
      dquat=drimus(bb,xe1,xe2,t9)
      quat=quat/(0.479d0*t9**(2.d0/3.d0)+0.00312d0)**2.d0
      dquat=dquat-2.d0*0.479d0*2.d0/3.d0*t9**(2.d0/3.d0)/(0.479d0*t9**(2.d0/3.d0)+0.00312d0)
      aa=2.8d+11
      bb=16.67d0
      cc=0.040d0
      xe1=-1.d0
      xe2=1.d0/3.d0
      xe3=2.d0
      cinq=secun(aa,bb,cc,xe1,xe2,xe3,t9)
      dcinq=dsecun(bb,cc,xe1,xe2,xe3,t9)
      aa=2.94d-03
      bb=0.767d0
      xe1=3.d0/2.d0
      xe2=1.d0
      six=primus(aa,bb,xe1,xe2,t9)
      dsix=drimus(bb,xe1,xe2,t9)
      aa=98.d0
      bb=2.077d0
      xe1=3.d0/2.d0
      xe2=1.d0
      sept=primus(aa,bb,xe1,xe2,t9)
      dsept=drimus(bb,xe1,xe2,t9)
      cent=uno*due+tre+quat+0.5d0*(cinq+six)+0.1d0*sept
      dcent=uno*duno*due+uno*ddue+tre*dtre+quat*dquat+0.5d0*(cinq*dcinq+six*dsix)+0.1d0*sept*dsept

      call interpol(34,t9,u)   ! fichier='o17a'
      tcent=u
      b117a(j1)=exp(rh1+fop)*tcent
      ep(10)=q117a*b117a(j1)*convMeVAvo/17.d0*xo17(j1)*xprot(j1)
      eprt(10)=rhp1+fop*fopp
      eptr(10)=dcent/cent+rht1+fop*fopt

!--- O17(P,G)F18 LANDRE ET AL A&A 240, 85 (1990)/nacre
      t9a=t9/(1.d0+2.69d0*t9)
      dt9a=1.d0-2.69d0*t9a
      uno=7.97d+07*t9a**(5.d0/6.d0)/t9**(3.d0/2.d0)*exp(-16.712d0/t9a**(1.d0/3.d0))
      duno=5.d0/6.d0*dt9a-3.d0/2.d0+1.d0/3.d0*16.712d0/t9a**(1.d0/3.d0)*dt9a
      aa=1.51d+08
      bb=16.712d0
      xe1=2.d0/3.d0
      xe2=1.d0/3.d0
      due=primus(aa,bb,xe1,xe2,t9)
      ddue=drimus(bb,xe1,xe2,t9)
      aa=1.d0
      bb=0.025d0
      cc=-0.051d0
      dd=-8.82d-03
      ee=0.d0
      fw=0.d0
      xe1=1.d0/3.d0
      xe2=2.d0/3.d0
      xe3=1.d0
      xe4=0.d0
      xe5=0.d0
      tre=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
      dtre=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
      aa=1.56d+05
      bb=6.272d0
      xe1=1.d0
      xe2=1.d0
      quat=primus(aa,bb,xe1,xe2,t9)
      dquat=drimus(bb,xe1,xe2,t9)
      aa=3.16d-05
      bb=0.767d0
      xe1=3.d0/2.d0
      xe2=1.d0
      cinq=primus(aa,bb,xe1,xe2,t9)
      dcinq=drimus(bb,xe1,xe2,t9)
      aa=98.d0
      bb=2.077d0
      xe1=3.d0/2.d0
      xe2=1.d0
      six=primus(aa,bb,xe1,xe2,t9)
      dsix=drimus(bb,xe1,xe2,t9)
      cent=uno+due*tre+quat+0.5d0*cinq+0.1d0*six
      dcent=uno*duno+due*ddue*tre+due*dtre+0.5d0*cinq*dcinq+0.1d0*six*dsix

      call interpol(35,t9,u)   ! fichier='o17g'
      tcent=u
      b117g(j1)=exp(rh1+fop)*tcent
      ep(11)=q117g*b117g(j1)*convMeVAvo/17.d0*xo17(j1)*xprot(j1)
      eprt(11)=rhp1+fop*fopp
      eptr(11)=cent/dcent+rht1+fop*fopt

!--- O18(P,G)F19 AT. DATA & NUCL. DATA TABLES 40, 283 (1988)/nacre
      aa=3.45d+08
      bb=16.729d0
      cc=0.139d0
      xe1=2.d0/3.d0
      xe2=1.d0/3.d0
      xe3=2.d0
      uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
      duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
      aa=1.d0
      bb=0.025d0
      cc=2.26d0
      dd=0.394d0
      ee=30.56d0
      fw=13.55d0
      xe1=1.d0/3.d0
      xe2=2.d0/3.d0
      xe3=1.0d0
      xe4=4.d0/3.d0
      xe5=5.d0/3.d0
      due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
      ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
      aa=1.25d-15
      bb=0.231d0
      xe1=3.d0/2.d0
      xe2=1.d0
      tre=primus(aa,bb,xe1,xe2,t9)
      dtre=drimus(bb,xe1,xe2,t9)
      aa=1.64d+02
      bb=1.670d0
      xe1=3.d0/2.d0
      xe2=1.d0
      quat=primus(aa,bb,xe1,xe2,t9)
      dquat=drimus(bb,xe1,xe2,t9)
      aa=1.28d+04
      bb=5.098d0
      xe1=-1.d0/2.d0
      xe2=1.d0
      cinq=primus(aa,bb,xe1,xe2,t9)
      dcinq=drimus(bb,xe1,xe2,t9)
      cent=uno*due+tre+quat+cinq
      dcent=uno*duno*due+uno*ddue+tre*dtre+quat*dquat+cinq*dcinq

      call interpol(38,t9,u)   ! fichier='o18g'
      tcent=u
      b118g(j1)=exp(rh1+fop)*tcent
      ep(13)=q118g*b118g(j1)*convMeVAvo/18.d0*xo18(j1)*xprot(j1)
      eprt(13)=rhp1+fop*fopp
      eptr(13)=dcent/cent+rht1+fop*fopt
    endif

!--- F19(a,p)NE22
    call screen(y(j1),xc12(j1),xo16(j1),xne20(j1),xmg24(j1),rh1,rhpsi,rhpsit,rhp1,rht1,vmyo,vmye,zw19(1),zw19(2),zw19(3), &
      zw19(4),zw19(5),t8ln,rhpsip,ffa,ffat,ffap,xc13(j1),xn14(j1),xn15(j1),xo17(j1),xo18(j1),xne22(j1),xmg25(j1),xmg26(j1),z)
    aa=4.50d+18
    bb=43.467d0
    cc=0.637d0
    xe1=2.d0/3.d0
    xe2=1.d0/3.d0
    xe3=2.d0
    uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
    duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
    aa=7.98d+04
    bb=12.760d0
    xe1=-3.d0/2.d0
    xe2=1.d0
    due=primus(aa,bb,xe1,xe2,t9)
    ddue=drimus(bb,xe1,xe2,t9)
    cent=uno+due
    dcent=uno*duno+due*ddue

    call interpol(43,t9,u)   ! fichier='f9ap'
    tcent=u
    e19ap(j1)=exp(rh1+ffa)*tcent
    w19ap=1.675d0*convMeVAvo/(19.d0*4.d0)*xf19(j1)*y(j1)*e19ap(j1)
    e19apt=rht1+ffat*ffa+dcent/cent
    e19app=rhp1+ffa*ffap
    if (j1 == m .and. verbose) then
      write(3,'(1x,a,4(1x,e24.18))') '19F(a,p):  ',e19ap(j1),w19ap,e19apt,e19app
    endif

!--- MG24(a,g)SI28, CF88
    call screen(y(j1),xc12(j1),xo16(j1),xne20(j1),xmg24(j1),rh1,rhpsi,rhpsit,rhp1,rht1,vmyo,vmye,zw24(1),zw24(2),zw24(3), &
      zw24(4),zw24(5),t8ln,rhpsip,fy24,fyt24,fyp24,xc13(j1),xn14(j1),xn15(j1),xo17(j1),xo18(j1),xne22(j1),xmg25(j1),xmg26(j1),z)
    aa=4.78d+01
    bb=13.506d0
    xe1=1.5d0
    xe2=1.d0
    uno=primus(aa,bb,xe1,xe2,t9)
    duno=drimus(bb,xe1,xe2,t9)
    aa=2.38d+03
    bb=15.218d0
    xe1=1.5d0
    xe2=1.d0
    due=primus(aa,bb,xe1,xe2,t9)
    ddue=drimus(bb,xe1,xe2,t9)
    aa=2.47d+02
    bb=15.147d0
    xe1=-1.5d0
    xe2=1.d0
    tre=primus(aa,bb,xe1,xe2,t9)
    dtre=drimus(bb,xe1,xe2,t9)
    aa=1.72d-09
    bb=5.028d0
    xe1=1.5d0
    xe2=1.d0
    quat=primus(aa,bb,xe1,xe2,t9)
    dquat=drimus(bb,xe1,xe2,t9)
    aa=1.25d-03
    bb=7.929d0
    xe1=1.5d0
    xe2=1.d0
    cinq=primus(aa,bb,xe1,xe2,t9)
    dcinq=drimus(bb,xe1,xe2,t9)
    aa=2.43d+01
    bb=11.523d0
    xe1=1.d0
    xe2=1.d0
    six=primus(aa,bb,xe1,xe2,t9)
    dsix=drimus(bb,xe1,xe2,t9)
    gt9=1.d0+5.d0*exp(-15.882d0/t9)
    dgt9=+15.882d0/t9*5.d0*exp(-15.882d0/t9)
    cent=uno+due+tre+0.1*(quat+cinq+six)
    dcent=uno*duno+due*ddue+tre*dtre+0.1d0*(quat*dquat+cinq*dcinq+six*dsix)
    cent=cent/gt9
    dcent=dcent/gt9-cent/(gt9*gt9)*dgt9

    call interpol(44,t9,u)   ! fichier='mg4a'
    tcent=u
    e24ag(j1)=exp(rh1+fy24)*tcent
    w24ag=9.984d0*convMeVAvo/(24.d0*4.d0)*xmg24(j1)*y(j1)*e24ag(j1)
    e24agt=rht1+fy24*fyt24+dcent/cent
    e24agp=rhp1+fy24*fyp24
    if (j1 == m .and. verbose) then
      write(3,'(1x,a,4(1x,e24.18))') '24Mg(a,g): ',e24ag(j1),w24ag,e24agt,e24agp
    endif

!--- O17(a,g)NE21, CF88
    gt9=1.d0+exp(-10.106d0/t9)/3.d0
    t9a=t9/(1.d0+0.1646d0*t9)
    ft9a=exp(-(0.786d0/t9a)**3.51d0)
    fpt9a=exp(-(t9a/1.084d0)**1.69d0)
    uno=1.73d+17*fpt9a/gt9*t9a**(5.d0/6.d0)/t932*exp(-39.914d0/t9a**(1.d0/3.d0))
    due=3.50d+15*ft9a/gt9*t9a**(5.d0/6.d0)/t932*exp(-39.914d0/t9a**(1.d0/3.d0))
    cent=uno+due
    dgt9=exp(-10.106d0/t9)/3.d0*10.106d0/t9
    dt9a=t9a*(1.d0-0.1646d0*t9a)
    dft9a=ft9a*3.51d0*(0.786d0/t9a)**3.51d0*dt9a/t9a
    dfpt9a=-fpt9a*1.69d0*((t9a/1.084d0)**0.69d0)*dt9a/1.084d0
    duno=dfpt9a/fpt9a-dgt9/gt9+5.d0/6.d0*dt9a/t9a-3.d0/2.d0+1.d0/3.d0*39.914d0/(t9a**(4.d0/3.d0))*dt9a
    if (ft9a == 0.d0) then
      ddue=-dgt9/gt9+5.d0/6.d0*dt9a/t9a-3.d0/2.d0+1.d0/3.d0*39.914d0/(t9a**(4.d0/3.d0))*dt9a
    else
      ddue=dft9a/ft9a-dgt9/gt9+5.d0/6.d0*dt9a/t9a-3.d0/2.d0+1.d0/3.d0*39.914d0/(t9a**(4.d0/3.d0))*dt9a
    endif
    dcent=duno*uno+ddue*due

    call interpol(47,t9,u)   ! fichier='o7ag'
    tcent=u
    e17ag(j1)=exp(rh1+fyo)*tcent
    w17ag=7.351d0*convMeVAvo/(17.d0*4.d0)*xo17(j1)*y(j1)*e17ag(j1)
    e17agt=rht1+fyo*fyto+dcent/cent
    e17agp=rhp1+fypo*fyo
    if (j1 == m .and. verbose) then
      write(3,'(1x,a,4(1x,e24.18))') '17O(a,g):  ',e17ag(j1),w17ag,e17agt,e17agp
    endif

!--- NE21(a,g)MG25, CF88
    aa=2.66d+07
    bb=22.049d0
    xe1=1.5d0
    xe2=1.d0
    uno=primus(aa,bb,xe1,xe2,t9)
    duno=drimus(bb,xe1,xe2,t9)
    t9a=t9/(1.d0+0.0537d0*t9)
    dt9a=t9a*(1.d0-0.0537d0*t9a)
    due=4.94d+19*t9a**(5.d0/6.d0)/t932*exp(-46.890d0/t9a**(1.d0/3.d0))
    tre=8.72d-03*t9-6.87d-04*t9*t9+2.15d-05*t9*t9*t9
    quat=1.52d-04*exp(-46.90d0/t913*tre)
    cinq=1.d0+1.5d0*exp(-4.068d0/t9)+2.0d0*exp(-20.258d0/t9)
    cent=(due+uno)*quat/cinq
    ddue=5.d0/6.d0*dt9a/t9a-3.d0/2.d0+1.d0/3.d0*46.890d0/(t9a**(4.d0/3.d0))*dt9a
    dtre=8.72d-03*t9-2.d0*6.87d-04*t9*t9+3.d0*2.15d-05*t9*t9*t9
    dquat=46.90d0*tre*1.d0/3.d0/t913-46.90d0/t913*dtre
    dcinq=1.5d0*exp(-4.068d0/t9)*4.068d0/t9+2.0d0*exp(-20.258d0/t9)*20.258d0/t9
    dcent=(due*ddue+uno*duno)*quat/cinq+(due+uno)*(quat*dquat/cinq-quat/(cinq*cinq)*dcinq)

    call interpol(45,t9,u)   ! fichier='ne1a'
    tcent=u
    e21ag(j1)=exp(rh1+fyn22)*tcent
    w21ag=9.882d0*convMeVAvo/(21.d0*4.d0)*xne21(j1)*y(j1)*e21ag(j1)
    e21agt=rht1+fyn22*fyt22+dcent/cent
    e21agp=rhp1+fyn22*fyp22
    if (j1 == m .and. verbose) then
     write(3,'(1x,a,4(1x,e24.18))') '21Ne(a,g): ',e21ag(j1),w21ag,e21agt,e21agp
    endif

!--- C14(a,g)O18 FL89
    aa=1.274d+4
    bb=10.29d0
    xe1=3.d0/2.d0
    xe2=1.d0
    uno=primus(aa,bb,xe1,xe2,t9)
    duno=drimus(bb,xe1,xe2,t9)
    aa=3.449d+4
    bb=16.15d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    due=primus(aa,bb,xe1,xe2,t9)
    ddue=drimus(bb,xe1,xe2,t9)
    aa=1.319d+4
    bb=18.94d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    tre=primus(aa,bb,xe1,xe2,t9)
    dtre=drimus(bb,xe1,xe2,t9)
    aa=9.035d+4
    bb=21.03d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    quat=primus(aa,bb,xe1,xe2,t9)
    dquat=drimus(bb,xe1,xe2,t9)
    aa=8.185d+4
    bb=22.03d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    cinq=primus(aa,bb,xe1,xe2,t9)
    dcinq=drimus(bb,xe1,xe2,t9)
    aa=5.727d+4
    bb=23.05d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    six=primus(aa,bb,xe1,xe2,t9)
    dsix=drimus(bb,xe1,xe2,t9)
    aa=8.555d+04
    bb=23.86d0
    xe1=3.d0/2.d0
    xe2=1.0d0
    sept=primus(aa,bb,xe1,xe2,t9)
    dsept=drimus(bb,xe1,xe2,t9)
    aa=2.839d+11
    bb=32.51d0
    xe1=2.d0/3.d0
    xe2=1.d0/3.d0
    huit=primus(aa,bb,xe1,xe2,t9)
    dhuit=drimus(bb,xe1,xe2,t9)
    aa=1.d0
    bb=2.564d-3
    cc=-0.3182d0
    dd=-2.856d-2
    ee=-6.808d-2
    fw=-1.554d-2
    xe1=1.d0/3.d0
    xe2=2.d0/3.d0
    xe3=1.d0
    xe4=4.d0/3.d0
    xe5=5.d0/3.d0
    neuf=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
    dneuf=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
    cent=uno+due+tre+quat+cinq+six+sept+huit*neuf
    dcent=uno*duno+due*ddue+tre*dtre+quat*dquat+cinq*dcinq+six*dsix+sept*dsept+huit*dhuit*neuf+huit*dneuf
    ec14ag(j1)=exp(rh1+fyc)*cent
    wc14ag=6.227d0*convMeVAvo/(14.d0*4.d0)*xc14(j1)*y(j1)*ec14ag(j1)
    ec14at=rht1+fyc*fytc+dcent/cent
    ec14ap=rhp1+fyc*fypc
    if (j1 == m .and. verbose) then
      write(3,'(1x,a,4(1x,e24.18))') '14C(a,g):  ',ec14ag(j1),wc14ag,ec14at,ec14ap
    endif

!--- N15(a,g)F19, CF88 nacre
    aa=2.54d+10
    bb=36.211d0
    cc=0.616d0
    xe1=2.d0/3.d0
    xe2=1.d0/3.d0
    xe3=2.d0
    uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
    duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
    aa=1.d0
    bb=0.012d0
    cc=1.69d0
    dd=0.136d0
    ee=1.91d0
    fw=0.391d0
    xe1=1.d0/3.d0
    xe2=2.d0/3.d0
    xe3=1.d0
    xe4=4.d0/3.d0
    xe5=5.d0/3.d0
    due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
    ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
    aa=9.83d-03
    bb=4.232d0
    xe1=3.d0/2.d0
    xe2=1.d0
    tre=primus(aa,bb,xe1,xe2,t9)
    dtre=drimus(bb,xe1,xe2,t9)
    aa=1.52d+03
    bb=9.747d0
    xe1=-1.d0
    xe2=1.d0
    quat=primus(aa,bb,xe1,xe2,t9)
    dquat=drimus(bb,xe1,xe2,t9)
    cent=uno*due+tre+quat
    dcent=uno*duno*due+uno*ddue+tre*dtre+quat*dquat

    call interpol(23,t9,u)   ! fichier='n5ag'
    tcent=u
    e15ag(j1)=exp(rh1+fyn14)*tcent
    w15ag=4.014d0*convMeVAvo/(15.d0*4.d0)*xn15(j1)*y(j1)*e15ag(j1)
    e15agt=rht1+fyn14*fyt14+dcent/cent
    e15agp=rhp1+fyn14*fyp14
    if (j1 == m .and. verbose) then
      write(3,'(1x,a,4(1x,e24.18))') '15N(a,g):  ',e15ag(j1),w15ag,e15agt,e15agp
    endif

!--- O18(a,n)NE21, CF 88/nacre
    t9a=t9/(1.d0+0.0483d0*t9+0.00569d0*t953/(1.d0+0.0483d0*t9)**(2.d0/3.d0))
    gt9=1.d0+5.d0*exp(-23.002d0/t9)
    ft9a=exp(-(0.431d0/t9a)**3.89d0)
    uno=7.22d+17*ft9a/gt9*t9a**(5.d0/6.d0)/t932*exp(-40.056d0/t9a**(1.d0/3.d0))
    due=150.31d0/gt9*exp(-8.045d0/t9)
    cent=uno+due
    dt9a1=1.d0/(1.d0+0.0483d0*t9)**(4.d0/3.d0)*(0.00569d0*5.d0/3.d0*t9**(5.d0/3.d0)*(1.d0+0.0483d0*t9)**(2.d0/3.d0)- &
          0.00569d0*t9**(5.d0/3.d0)*2.d0/3.d0*(1.d0+0.0483d0*t9)**(-1.d0/3.d0)*0.0483d0*t9)
    dt9a2=0.0483d0*t9+dt9a1
    dt9a=t9a*(1.d0-t9a/t9*dt9a2)
    dgt9=5.d0*exp(-23.002d0/t9)*23.002d0/t9
    dft9a=ft9a*3.89d0*(0.431d0/t9a)**3.89d0*dt9a/t9a
    if (ft9a == 0.0d0) then
      duno=-dgt9/gt9+5.d0/6.d0*dt9a/t9a-3.d0/2.d0+40.056d0*1.d0/3.d0/(t9a**(4.d0/3.d0))*dt9a
    else
      duno=dft9a/ft9a-dgt9/gt9+5.d0/6.d0*dt9a/t9a-3.d0/2.d0+40.056d0*1.d0/3.d0/(t9a**(4.d0/3.d0))*dt9a
    endif
    ddue=-dgt9/gt9+8.045/t9
    dcent=uno*duno+due*ddue

    call interpol(39,t9,u)   ! fichier='o18n'
    tcent=u
    e18an(j1)=exp(rh1+fyo)*tcent
    w18an=-0.693d0*convMeVAvo/(18.d0*4.d0)*xo18(j1)*y(j1)*e18an(j1)
    e18ant=rht1+fyo*fyto+dcent/cent
    e18anp=rhp1+fyo*fypo
    if (j1 == m .and. verbose) then
      write(3,'(1x,a,4(1x,e24.18))') '18O(a,n):  ',e18an(j1),w18an,e18ant,e18anp
    endif

! population III: si on est encore dans la fusion H, on saute les reactions de capture de neutrons
    if (ipop3 == 0 .or. x(j1) <= 1.0d-7) then

!--- NE21(n,a)O18, CF88
      aa=7.84d-01
      bb=-8.045d0
      xe1=0.d0
      xe2=1.d0
      revra=primus(aa,bb,xe1,xe2,t9)
      drevra=drimus(bb,xe1,xe2,t9)

      call interpol(46,t9,u)   ! fichier='ne1n'
      tcent=u
      e21na(j1)=exp(rh1)*tcent
      w21na=0.693d0*convMeVAvo/21.d0*xne21(j1)*xneut(j1)*e21na(j1)
      e21nat=rht1+dcent/cent+drevra
      e21nap=rhp1
      if (j1 == m .and. verbose) then
        write(3,'(1x,a,4(1x,e24.18))') '21Ne(n,a): ',e21na(j1),w21na,e21nat,e21nap
      endif
    endif

!--- MG25(a,n)SI28, CF88/nacre
    t9a=t9/(1.d0+0.063d0*t9)
    gt9=1.d0+10.d0*exp(-13.180d0/t9)/3.d0
    cent=3.59d+20/gt9*t9a**(5.d0/6.d0)/t932*exp(-53.410d0/t9a**(1.d0/3.d0))
    dt9a=t9a*(1.d0-0.063d0*t9a)
    dgt9=10.d0*exp(-13.180d0/t9)/3.d0*13.180d0/t9
    dcent=-dgt9/gt9+5.d0/6.d0*dt9a/t9a-3.d0/2.d0+1.d0/3.d0*53.410d0/(t9a**(4.d0/3.d0))*dt9a

    call interpol(17,t9,u)   ! fichier='mg5n'
    tcent=u
    e25an(j1)=exp(rh1+fy24)*tcent
    w25an=2.653d0*convMeVAvo/(25.d0*4.d0)*xmg25(j1)*y(j1)*e25an(j1)
    e25ant=rht1+fy24*fyt24+dcent
    e25anp=rhp1+fy24*fyp24
    if (j1 == m .and. verbose) then
      write(3,'(1x,a,4(1x,e24.18))') '25Mg(a,n): ',e25an(j1),w25an,e25ant,e25anp
    endif

! population III: si on est encore dans la fusion H, on saute les reactions de capture de neutrons
    if (ipop3 == 0 .or. (x(j1) <= 1.0d-7 .or. phase > 1)) then

!--- F18(n,a)N15, CF 88
      uno=3.14d+08*(1.d0-0.641d0*t912+0.108d0*t9)
      revrat=2.00d+00
      cent=uno*revrat
      duno=3.14d+08*(-1.d0/2.d0*0.641d0*t912+0.108d0*t9)
      dcent=duno*revrat
      ef18na(j1)=exp(rh1)*cent
      wf18na=6.418d0*convMeVAvo/18.d0*xf18(j1)*xneut(j1)*ef18na(j1)
      ef18nt=rht1+dcent/cent
      f18nap=rhp1
      if (j1 == m .and. verbose) then
        write(3,'(1x,a,4(1x,e24.18))') '18F(n,a):  ',ef18na(j1),wf18na,ef18nt,f18nap
      endif
! AL26G(N,A)NA23, CF88
      al26tn=3.38d+06*exp(0.388d0*t9+9.08d-03*t9*t9-2.07d-03*t9*t9*t9)
      dal26t=0.388d0*t9+2.d0*9.08d-03*t9*t9-3.d0*2.07d-03*t9*t9*t9
      gt9=1.d0+exp(-4.612d0/t9-5.623d-04+7.460d-02*t9)
      dgt9=(gt9-1.d0)*(4.612d0/t9+7.460d-02*t9)
      gpt9=1.d0+exp(-3.573d0/t9-1.008d0+0.1357d0*t9)
      dgpt9=(gpt9-1.d0)*(3.573d0/t9+0.1357d0*t9)
      aa=5.43d+07
      bb=0.9653d0
      xe1=3.d0/2.d0
      xe2=1.d0
      uno=primus(aa,bb,xe1,xe2,t9)
      duno=drimus(bb,xe1,xe2,t9)
      aa=6.97d+07
      bb=1.494d0
      xe1=-2.d0/7.d0
      xe2=1.d0
      due=primus(aa,bb,xe1,xe2,t9)
      ddue=drimus(bb,xe1,xe2,t9)
      al26mn=4.18d+06+uno+due
      dal26m=uno*duno+due*ddue
      aa=9.09d-02
      bb=2.651d0
      xe1=0.0d0
      xe2=1.0d0
      tre=primus(aa,bb,xe1,xe2,t9)
      dtre=drimus(bb,xe1,xe2,t9)
      cent=al26tn*gpt9/gt9-al26mn*tre
      dcent=al26tn*dal26t*gpt9/gt9+al26tn*(dgpt9/gt9-gpt9/(gt9*gt9)*dgt9)-dal26m*tre-al26mn*tre*dtre
      a26ga(j1)=exp(rh1)*cent
      w26ga=2.968d0*convMeVAvo/26.d0*xal26(j1)*xneut(j1)*a26ga(j1)
      a26gat=rht1+dcent/cent
      a26gap=rhp1
      if (j1 == m .and. verbose) then
       write(3,'(1x,a,4(1x,e24.18))') '26Al(n,a): ',a26ga(j1),w26ga,a26gat,a26gap
      endif

!--- AL26G(n,p)MG26, CF88
      gt9=1.d0+exp(-19.30d0/t9+0.6642d0+0.1386d0*t9)
      gpt9=1.d0+exp(-3.573d0/t9-1.008d0+0.1357d0*t9)
      al26tn=3.09d+07*exp(0.0731d0*t9+0.0381d0*t9*t9-3.22d-03*t9*t9*t9)
      aa=1.84d+05
      bb=0.043d0
      xe1=1.5d0
      xe2=1.d0
      uno=primus(aa,bb,xe1,xe2,t9)
      duno=drimus(bb,xe1,xe2,t9)
      aa=2.28d+07
      bb=0.348d0
      xe1=1.5d0
      xe2=1.d0
      due=primus(aa,bb,xe1,xe2,t9)
      ddue=drimus(bb,xe1,xe2,t9)
      aa=6.54d+08
      bb=0.826d0
      xe1=3.d0/8.d0
      xe2=1.d0
      tre=primus(aa,bb,xe1,xe2,t9)
      dtre=drimus(bb,xe1,xe2,t9)
      quat=1.60d+04/t932
      al26mn=uno+due+tre+0.1d0*quat
      cent=al26tn*gpt9/gt9-al26mn*9.09d-02*exp(-2.651d0/t9)
      dgt9=exp(-19.30d0/t9+0.6642d0+0.1386d0*t9)*(19.30d0/t9+0.1386d0*t9)
      dgpt9=exp(-3.573d0/t9-1.008d0+0.1357d0*t9)*(3.573d0/t9+0.1357d0*t9)
      dal26t=0.0731d0*t9+2.d0*0.0381d0*t9*t9-3.d0*3.22d-03*t9*t9*t9
      dquat=-3.d0/2.d0
      dal26m=uno*duno+due*ddue+tre*dtre+0.1d0*quat*dquat
      dcent=al26tn*dal26t*gpt9/gt9+al26tn*dgpt9/gt9-al26tn*gpt9/(gt9*gt9)*dgt9-dal26m*9.09d-02*exp(-2.651d0/t9)- &
            al26mn*9.09d-02*exp(-2.651d0/t9)*2.651d0/t9
      a26gp(j1)=exp(rh1)*cent
      w26gp=4.786d0*convMeVAvo/26.d0*xal26(j1)*xneut(j1)*a26gp(j1)
      a26gpt=rht1+dcent/cent
      a26gpp=rhp1
      if (j1 == m .and. verbose) then
        write(3,'(1x,a,4(1x,e24.18))') '26Al(n,p): ',a26gp(j1),w26gp, a26gpt,a26gpp
      endif
! On a <sigma v>/v_th [mbarn] pour v_th= 30 keV, on le convertit en cgs:
! facteur de conversion: sqrt(2*30*1000*e*1e7/m_H)*1e-27/m_H
      conver=sqrt(2.d0*30.d0*1000.d0*cst_e*1.d7/cst_mh)*1.d-27/cst_mh

!--- N14(n,p)C14 cf SMOKER
      cent=0.74d0*conver
      e14np(j1)=exp(rh1)*cent
      w14np=0.626d0*convMeVAvo/14.d0*xn14(j1)*xneut(j1)*e14np(j1)
      e14npt=rht1
      e14npp=rhp1
      if (j1 == m .and. verbose) then
        write(3,'(1x,a,4(1x,e24.18))') '14N(n,p):  ',e14np(j1),w14np,e14npt,e14npp
      endif

!--- NE20(n,g)NE21, BEER & VOSS 91
      e20ng(j1)=exp(rh1)*1.190d-01*conver
      w20ng=6.761d0*convMeVAvo/20.d0*xne20(j1)*xneut(j1)*e20ng(j1)

!--- C12(n,g)C13 NIT91
      e12ng(j1)=exp(rh1)*1.68d-02*conver
      w12ng=4.946d0*convMeVAvo/12.d0*xc12(j1)*xneut(j1)*e12ng(j1)

!--- N14(n,g)N15
      e14ng(j1)=exp(rh1)*4.1d-02*conver
      w14ng=10.833d0*convMeVAvo/14.d0*xn14(j1)*xneut(j1)*e14ng(j1)

!--- F19(n,g)F20(B-)NE20
      e19ng(j1)=exp(rh1)*5.7d0*conver
      w19ng=9.943d0*convMeVAvo/19.d0*xf19(j1)*xneut(j1)*e19ng(j1)

!--- NE21(n,g)NE22 BAO & KAEPPELER 86
      e21ng(j1)=exp(rh1)*1.5d0*conver
      w21ng=10.370d0*convMeVAvo/21.d0*xne21(j1)*xneut(j1)*e21ng(j1)

!---- NE22(n,g)NA23 BEER & VOSS 91
      e22ng(j1)=exp(rh1)*6.0d-02*conver
      w22ng=9.576d0*convMeVAvo/22.d0*xne22(j1)*xneut(j1)*e22ng(j1)

!--- NA23(n,g)MG24 BEER & VOSS 91
      e23ng(j1)=exp(rh1)*2.1d0*conver
      w23ng=9.548d0*convMeVAvo/23.d0*xna23(j1)*xneut(j1)*e23ng(j1)

!--- MG24(n,g)MG25 BEER & VOSS 91
      e24ng(j1)=exp(rh1)*4.2d0*conver
      w24ng=7.332d0*convMeVAvo/24.d0*xmg24(j1)*xneut(j1)*e24ng(j1)

!--- MG25(n,g)MG26 BEER & VOSS 91
      e25ng(j1)=exp(rh1)*6.5d0*conver
      w25ng=11.094d0*convMeVAvo/25.d0*xmg25(j1)*xneut(j1)*e25ng(j1)

!--- MG26(n,g)AL27 BEER & VOSS 91
      e26ng(j1)=exp(rh1)*6.6d-02*conver
      w26ng=7.601d0*convMeVAvo/26.d0*xmg26(j1)*xneut(j1)*e26ng(j1)

!--- AL27(n,g)SI28 BEER & VOSS 91
      e27ng(j1)=exp(rh1)*3.8d0*conver
      w27ng=5.568d0*convMeVAvo/27.d0*xal27(j1)*xneut(j1)*e27ng(j1)

!--- SI28(n,g)SI29 BIDONS
      fluneu=xneut(j1)*exp(rh1)*cst_avo*dzeit*sqrt(163d+06*exp(t(j1)))/1.0d+27
      xfluxn(j1)=fluneu
      flunet=1.d0/2.d0+rht1
      flunep=rhp1
      if (fluneu <= 1.72d0) then
         e28ngh=43.8372d0*fluneu+10.4d0
         e28nht=43.8372d0*fluneu*flunet
         e28nhp=43.8372d0*fluneu*flunep
      endif
      if (fluneu > 1.72d0 .and. fluneu < 2.92d0) then
         e28ngh=-69.8417d0*fluneu+205.93d0
         e28nht=-69.8417d0*fluneu*flunet
         e28nhp=-69.8417d0*fluneu*flunep
      endif
      if (fluneu >= 2.92d0) then
         e28ngh=1.99d0
         e28nht=0.d0
         e28nhp=0.d0
      endif
      e28ng(j1)=(5.11404d0+e28ngh)*conver*exp(rh1)
      w28ng=8.221d0*convMeVAvo/40.d0*xbid(j1)*xneut(j1)*e28ng(j1)
      e28ngt=rht1+e28nht/(5.11404d0+e28ngh)
      e28ngp=rhp1+e28nhp/(5.11404d0+e28ngh)
      if (j1 == m .and. verbose) then
       write(3,'(1x,a,4(1x,e24.18))') '28Si(n,g): ',e28ng(j1),w28ng,e28ngt,e28ngp
      endif

!--- C14(n,g)C15 THIELEMANN ET AL 91, Nuclei in the COSMOS
      uno=3.24d+03*t9
      aa=2.05d+03
      bb=21.14d0
      xe1=3.d0/2.d0
      xe2=1.d0
      due=primus(aa,bb,xe1,xe2,t9)
      ddue=drimus(bb,xe1,xe2,t9)
      cent=uno+due
      dcent=uno+due*ddue
      ec14ng(j1)=exp(rh1)*cent
      wc14ng=10.990d0*convMeVAvo/14.d0*xc14(j1)*xneut(j1)*ec14ng(j1)
      ec14nt=rht1+dcent/cent
      if (j1 == m .and. verbose) then
        write(3,'(1x,a,4(1x,e24.18))') '14C(n,g):  ',ec14ng(j1),wc14ng,ec14nt,ec14nt
      endif

!--- O18(n,g)O19 THIELEMANN ET AL.
      uno=2.12d+01
      aa=2.55d+03
      bb=1.769d0
      xe1=3.d0/2.d0
      xe2=1.d0
      due=primus(aa,bb,xe1,xe2,t9)
      ddue=drimus(bb,xe1,xe2,t9)
      cent=uno+due
      dcent=due*ddue
      e18ng(j1)=exp(rh1)*cent
      w18ng=6.199d0*convMeVAvo/18.d0*xo18(j1)*xneut(j1)*e18ng(j1)
      e18ngt=rht1+dcent/cent

! e18ngp jamais calcule. Mis a 0 pour l'impression.
      e18ngp = 0.d0
      if (j1 == m .and. verbose) then
        write(3,'(1x,a,4(1x,e24.18))') '18O(n,g):  ',e18ng(j1),w18ng,e18ngt,e18ngp
      endif

!--- F18(n,p)O18 SMOKER
      ef18np(j1)=0.36d+08*exp(rh1)
      wf18np=2.438d0*convMeVAvo/18.d0*xf18(j1)*xneut(j1)*ef18np(j1)
      f18npt=rht1
      f18npp=rhp1
      if (j1 == m .and. verbose) then
       write(3,'(1x,a,4(1x,e24.18))') '18F(n,p):  ',ef18np(j1),wf18np,f18npt,f18npp
      endif
    endif

!--- AL26G(B+)MG26 (BETA+ DECAY ASSUME INSTANTANEOUS ANIHILATION E,E+)
    e26be(j1)=3.05d-14
    w26be=xal26(j1)*2.361d0*convMeVAvo/26.d0*e26be(j1)

!--- F18(B+)O18
    e18be(j1)=1.05d-04
    w18be=1.258d0*convMeVAvo/18.d0*xf18(j1)*e18be(j1)

!--- C14(B-)N14
    e14be(j1)=3.83d-12
    w14be=0.054d0*convMeVAvo/14.d0*xc14(j1)*e14be(j1)
  endif
!-----------------------------------------------------------------------
! modif PopIII: on teste si on est encore dans la fusion H.
!               Si oui, on ne prend pas en compte epspt2
  if (ipop3 == 0) then
    goto 263
  endif
  if (x(j1)-1.0d-7) 263,263,264

! ---   Fusion He seul   ---
263 epsy(j1)=wpsyy+wpsyc+wpsyo+e144+e184+e224+e224g+e134+w17an+w20ag
  epsp2=(eyp*wpsyy+ecp*wpsyc+eop*wpsyo+e144p*e144+w20ag*e20agp+e184p*e184+e224p*(e224+e224g)+ecp*e134+w17an*e17anp)
  epst2=(eyt*wpsyy+ect*wpsyc+eot*wpsyo+e144t*e144+e184t*e184+w20ag*e20agt+e224t*e224+e22tg*e224g+e13t*e134+w17an*e17ant)

  if (ialflu == 1 .and. t(j1) > 18.42d0) then
    epsy(j1)=epsy(j1)+w24ag+w17ag+w21ag+w18an+w21na+w25an+w20ng+w21ng+w22ng+w23ng+w24ng+w25ng+w26ng+w27ng+w28ng+w26ga+ &
             w26gp+w14np+wc14pg+wc14ag+wf18na+w15ag+wf18np+w18pa+wc14ng+w19ap+w14be+w18be+w26be+w18ng+w12ng+w14ng+w19ng+ &
             ep(4)+ep(5)+ep(6)+ep(7)+ep(8)+ep(9)+ep(10)+ep(11)+ep(13)
  epsp2=epsp2+e24agp*w24ag+e17agp*w17ag+e21agp*w21ag+e18anp*w18an+e21nap*w21na+e25anp*w25an+rhp1*(w20ng+w12ng+w14ng+w19ng+ &
        w21ng+w22ng+w23ng+w24ng+w25ng+w26ng+w27ng+wc14ng+w18ng)+e28ngp*w28ng+a26gap*w26ga+a26gpp*w26gp+e14npp*w14np+ &
        ec14pp*wc14pg+ec14ap*wc14ag+f18nap*wf18na+ep(4)*eprt(4)+ep(5)*eprt(5)+ep(6)*eprt(6)+ep(7)*eprt(7)+ep(8)*eprt(8)+ &
        ep(9)*eprt(9)+ep(10)*eprt(10)+ep(11)*eprt(11)+ep(13)*eprt(13)+e15agp*w15ag+f18npp*wf18np+e18pap*w18pa+e19app*w19ap
  epst2=epst2+e24agt*w24ag+e17agt*w17ag+e21agt*w21ag+e18ant*w18an+e21nat*w21na+e25ant*w25an+rht1*(w20ng+w12ng+w14ng+w19ng+ &
        w21ng+w22ng+w23ng+w24ng+w25ng+w26ng+w27ng)+e28ngt*w28ng+ec14nt*wc14ng+e18ngt*w18ng+a26gat*w26ga+a26gpt*w26gp+ &
        e14npt*w14np+ec14pt*wc14pg+ec14at*wc14ag+ef18nt*wf18na+ep(4)*eptr(4)+ep(5)*eptr(5)+ep(6)*eptr(6)+ep(7)*eptr(7)+ &
        ep(8)*eptr(8)+ep(9)*eptr(9)+ep(10)*eptr(10)+ep(11)*eptr(11)+ep(13)*eptr(13)+e15agt*w15ag+f18npt*wf18np+e18pat*w18pa+ &
        e19apt*w19ap
  endif

  if (epsy(j1) /= 0.d0) then
    epsp2=epsp2/epsy(j1)
    epst2=epst2/epsy(j1)
  else
    epsp2 = 0.d0
    epst2 = 0.d0
  endif
  epsp1=epsp2
  epst1=epst2
  en=sqrt(abs(epsy(j1)*epsy(j)))

! Calcul de l'energie produite, en Mev, par particule alpha brulee
  dy=abs(y(j1)/4.d0*(0.5d0*y(j1)/4.d0*y(j1)/4.d0*epsyy(j1)+xc12(j1)/12.d0*epsyc(j1)+xn14(j1)/14.d0*c144(j1)+ &
         xo16(j1)/16.d0*epsyo(j1)+xo18(j1)/18.d0*c184(j1)+xne20(j1)/20.d0*e20ag(j1)+xne22(j1)/22.d0*c224g(j1)+ &
         xc13(j1)/13.d0*c134(j1)+xo17(j1)/17.d0*e17an(j1)+xne22(j1)/22.d0*c224(j1)))
  if (ialflu == 1 .and. t(j1) > 18.42d0) then
    dy=abs(y(j1)/4.d0*(0.5d0*y(j1)/4.d0*y(j1)/4.d0*epsyy(j1)+xc12(j1)/12.d0*epsyc(j1)+xc13(j1)/13.d0*c134(j1)+ &
           xc14(j1)/14.d0*ec14ag(j1)+xn14(j1)/14.d0*c144(j1)+xn15(j1)/15.d0*e15ag(j1)+xo16(j1)/16.d0*epsyo(j1)+ &
           xo17(j1)/17.d0*(e17ag(j1)+e17an(j1))+xo18(j1)/18.d0*(c184(j1)+e18an(j1))+xf19(j1)/19.d0*e19ap(j1)+ &
           xne20(j1)/20.d0*e20ag(j1)+xne21(j1)/21.d0*e21ag(j1)+xne22(j1)/22.d0*(c224(j1)+c224g(j1))+ &
           xmg24(j1)/24.d0*e24ag(j1)+xmg25(j1)/25.d0*e25an(j1)))
  endif
  if (dy /= 0.d0) then
    zensi(j1)=epsy(j1)/(convMeVAvo*dy)
  else
    zensi(j1)=1.d0
  endif
  goto 265

! ---   Fusion H   ---
264 do k=1,13
   eps(j1)=eps(j1)+ep(k)
   epst1=epst1+ep(k)*(eptr(k)+rht1*eprt(k))
   epsp1=epsp1+ep(k)*rhp1*eprt(k)
  enddo
  if (ialflu == 1) then
    do k=14,23
     eps(j1)=eps(j1)+ep(k)
     epst1=epst1+ep(k)*(eptr(k)+rht1*eprt(k))
     epsp1=epsp1+ep(k)*rhp1*eprt(k)
    enddo
  endif

  if (ippcno == 1) then
    do k=1,13
     en13(k,j1)=convMeVAvo*ep(k)
    enddo
    enpp(j1)=convMeVAvo*(ep(1)+ep(2)+ep(3))
    encno(j1)=convMeVAvo*(ep(4)+ep(5)+ep(6)+ep(7)+ep(8)+ep(9)+ep(10)+ep(11)+ep(12)+ep(13))
  endif

  epsy(j1)=wpsyy+wpsyc+wpsyo+e144+e184+e224+e224g+e134+w17an+w20ag
  epsp2=(eyp*wpsyy+ecp*wpsyc+eop*wpsyo+e144p*e144+w20ag*e20agp+e184p*e184+e224p*(e224+e224g)+ecp*e134+w17an*e17anp)
  epst2=(eyt*wpsyy+ect*wpsyc+eot*wpsyo+e144t*e144+e184t*e184+w20ag*e20agt+e224t*e224+e22tg*e224g+e13t*e134+w17an*e17ant)

  if (ialflu == 1 .and. t(j1) > 18.42d0) then
    epsy(j1)=epsy(j1)+w24ag+w17ag+w21ag+w18an+w25an+wc14pg+wc14ag+w15ag+w19ap+w14be+w18be+w26be
    epsp2=epsp2+e24agp*w24ag+e17agp*w17ag+e21agp*w21ag+e18anp*w18an+e25anp*w25an+ec14pp*wc14pg+ec14ap*wc14ag+ &
                e15agp*w15ag+e19app*w19ap
    epst2=epst2+e24agt*w24ag+e17agt*w17ag+e21agt*w21ag+e18ant*w18an+e25ant*w25an+ec14pt*wc14pg+ec14at*wc14ag+ &
                e15agt*w15ag+e19apt*w19ap
  endif

  if (eps(j1) /= 0.d0 .and. epsy(j1) /= 0.d0) then
    epst1=epst1/(eps(j1)+(epsy(j1)/convMeVAvo))
    epsp1=epsp1/(eps(j1)+(epsy(j1)/convMeVAvo))
    eps(j1)=convMeVAvo*eps(j1)
    epsp2=epsp2/(epsy(j1)+eps(j1))
    epst2=epst2/(epsy(j1)+eps(j1))
  else
    if(epsy(j1) /= 0.d0) then
      epsp2=epsp2/epsy(j1)
      epst2=epst2/epsy(j1)
      epst1=0.d0
      epsp1=0.d0
    else
      if(eps(j1) /= 0.d0) then
        epsp1=epsp1/eps(j1)
        epst1=epst1/eps(j1)
        eps(j1)=convMeVAvo*eps(j1)
        epsp2=0.d0
        epst2=0.d0
      else
        epsp1=0.d0
        epst1=0.d0
        epsp2=0.d0
        epst2=0.d0
      endif
    endif
  endif
  epsp1=epsp1+epsp2
  epst1=epst1+epst2
  dy=abs(yab(2)*b33(j1)*yab(2)+yab(1)*(-1.5d0*b11(j1)*yab(1)-b112(j1)*yab(4)-b113(j1)*yab(5)-b114(j1)*yab(6)-(b115a(j1)+ &
         b115g(j1))*yab(7)-b116(j1)*yab(8)-(b117a(j1)+b117g(j1))*yab(9)-(b118a(j1)+2.d0*b118g(j1))*yab(10)))

  if(ialflu == 1) then
    dy=abs(yab(2)*b33(j1)*yab(2)+yab(1)*(-1.5d0*b11(j1)*yab(1)-b112(j1)*yab(4)-b113(j1)*yab(5)-b114(J1)*yab(6)-(b115a(j1)+ &
           b115g(j1))*yab(7)-b116(j1)*yab(8)-(b117a(j1)+b117g(j1))*yab(9)-(b118a(j1)+b118g(j1))*yab(10)-(b119a(J1)+ &
           b119g(j1))*yab(11)-b120(j1)*yab(12)-b121(j1)*yab(13)-b122(j1)*yab(14)-(b123g(j1)+b123a(j1))*yab(15)- &
           b124(j1)*yab(16)-(b125g(j1)+b125m(j1))*yab(17)-b1al26(j1)*yab(19)-b1mg26(j1)*yab(18)-(b127g(j1)+ &
           b127a(j1))*yab(20))+e19ap(j1)*yab(11))
  endif
  if (isnan(eps(j1)) .or. isnan(epsy(j1)) .or. isnan(dy)) then
    write(*,*) 'j1:eps,epsy,dy',j1,eps(j1),epsy(j1),dy
    stop
  endif

  if (dy /= 0.d0) then
    zensi(j1)=eps(j1)/(convMeVAvo*dy)
  else
    zensi(j1)=1.d0
  endif

! Calcul de l'energie produite, en Mev, par particule alpha brulee
  dy2=abs(y(j1)/4.d0*(0.5d0*y(j1)/4.d0*y(j1)/4.d0*epsyy(j1)+xc12(j1)/12.d0*epsyc(j1)+xn14(j1)/14.d0*c144(j1)+ &
          xo16(j1)/16.d0*epsyo(j1)+xo18(j1)/18.d0*c184(j1)+xne20(j1)/20.d0*e20ag(j1)+xne22(j1)/22.d0*(c224(j1)+ &
          c224g(j1))+xc13(j1)/13.d0*c134(j1)+xo17(j1)/17.d0*e17an(j1)))
  if (ialflu == 1 .and. t(j1) > 18.42d0) then
    dy2=abs(y(j1)/4.d0*(0.5d0*y(j1)/4.d0*y(j1)/4.d0*epsyy(j1)+xc12(j1)/12.d0*epsyc(j1)+xc13(j1)/13.d0*c134(j1)+ &
            xc14(j1)/14.d0*ec14ag(j1)+xn14(j1)/14.d0*c144(j1)+xn15(j1)/15.d0*e15ag(j1)+xo16(j1)/16.d0*epsyo(j1)+ &
            xo17(j1)/17.d0*(e17ag(j1)+e17an(j1))+xo18(j1)/18.d0*(c184(j1)+e18an(j1))+xf19(j1)/19.d0*e19ap(j1)+ &
            xne20(j1)/20.d0*e20ag(j1)+xne21(j1)/21.d0*e21ag(j1)+xne22(j1)/22.d0*(c224(j1)+c224g(j1))+ &
            xmg24(j1)/24.d0*e24ag(j1)+xmg25(j1)/25.d0*e25an(j1)))
  endif
  if (dy2 /= 0.d0) then
    zensi2(j1)=epsy(j1)/(convMeVAvo*dy2)
  else
    zensi2(j1)=1.d0
  endif
  en=sqrt(abs(eps(j1)+epsy(j1))*abs(eps(j)+epsy(j)))
  zensi(j1)= max(zensi(j1),zensi2(j1))

265 continue
  go to 24

! ==================== COMBUSTION CARBONE SEULEMENT ===================
27 continue
  if (t8ln-0.2303d0) 24,24,28
28 continue
  if (xc12(j1)) 24,36,36
36 continue
  yab(1)= y(j1)/4.d0
! C12(C12,A)NE20, CF88
  call screen(0.d0,xc12(j1),xo16(j1),xne20(j1),xmg24(j1),rh1,rhpsi,rhpsit,rhp1,rht1,vmyo,vmye,zw12(1),zw12(2),zw12(3),&
    zw12(4),zw12(5),t8ln,rhpsip,fy12,fyt12,fyp12,xc13(j1),xn14(j1),xn15(j1),xo17(j1),xo18(j1),xne22(j1),xmg25(j1),xmg26(j1),z)
  t9a=t9/(1.d0+0.0396d0*t9)
  dt9a=1.d0-0.0396d0*t9a
  cent=4.27d+26*t9a**(5.d0/6.d0)/t932*exp(-84.165d0/t9a**(1.d0/3.d0)-2.12d-03*t9**3.d0)
  dcent=5.d0/6.d0*dt9a-3.d0/2.d0+1.d0/3.d0*84.165d0/t9a**(1.d0/3.d0)*dt9a-3.d0*2.12d-03*t9**3.d0
  if (t9 < 1.75d0) then
     c12c12p=0.44d0
     c12c12a=0.56d0
  else
     c12c12p=0.45d0
     c12c12a=0.50d0
  endif
  if (t9 > 3.3d0) then
     c12c12p=0.40d0
     c12c12a=0.53d0
  endif
  epcne(j1)=1.548d+16*c12c12a*exp(rh1+fy12)*xc12(j1)*xc12(j1)*cent
! C12(C12,P)NA23, CF88
  epcna(j1)=0.751d+16*c12c12p*exp(rh1+fy12)*xc12(j1)*xc12(j1)*cent
  ect12=rht1+fy12*fyt12+dcent
  ecp12=rhp1+fy12*fyp12
! YAC EST L ABONDANCE EN PARTICULES ALPHA
  nbouc=1
  cya=1.d-10
901 n1620=2
  y(j1)=cya
! O16(A,G)NE20, CF88
  goto 300

302 continue
  n1216=2
  y(j1)=cya
! C12(A,G)O16, Fowler 75 * 3
  goto 350

352 continue

! NA23(P,A)NE20
! L'abondance des protons est supposee a l'equilibre
! en consequence on a que epsilon(NA23(P,A)NE20)=
! Q(NA23(P,A)NE20)/Q(C12(C12,P)NA23)*epsilon(C12(C12,P)NA23)
  ep23=1.0611d0*epcna(j1)
  e23p1=ecp12
  e23t1=ect12
  n2024=1
  call screen(0.d0,xc12(j1),xo16(j1),xne20(j1),xmg24(j1),rh1,rhpsi,rhpsit,rhp1,rht1,vmyo,vmye,zw20(1),zw20(2),zw20(3),&
    zw20(4),zw20(5),t8ln,rhpsip,fy20,fyt20,fyp20,xc13(j1),xn14(j1),xn15(j1),xo17(j1),xo18(j1),xne22(j1),xmg25(j1),xmg26(j1),z)
! NE20(A,G)MG24,CAUGHLAN & FOWLER AT. DAT. & NUCL. DAT. TA. 40,283,1988
! V0TO1 0.1
  aa=4.11d+11
  bb=46.766d0
  cc=2.219d0
  xe1=2.d0/3.d0
  xe2=1.d0/3.d0
  xe3=2.d0
  uno=secun(aa,bb,cc,xe1,xe2,xe3,t9)
  duno=dsecun(bb,cc,xe1,xe2,xe3,t9)
  aa=1.d0
  bb=0.009d0
  cc=0.882d0
  dd=0.055d0
  ee=0.749d0
  fw=0.119d0
  xe1=1.d0/3.d0
  xe2=2.d0/3.d0
  xe3=1.0d0
  xe4=4.d0/3.d0
  xe5=5.d0/3.d0
  due=tertiu(aa,bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  ddue=dterti(bb,cc,dd,ee,fw,xe1,xe2,xe3,xe4,xe5,t9)
  aa=5.27d+03
  bb=15.869d0
  xe1=3.d0/2.d0
  xe2=1.0d0
  tre=primus(aa,bb,xe1,xe2,t9)
  dtre=drimus(bb,xe1,xe2,t9)
  aa=6.51d+03
  bb=16.223d0
  xe1=-1.d0/2.d0
  xe2=1.0d0
  quat=primus(aa,bb,xe1,xe2,t9)
  dquat=drimus(bb,xe1,xe2,t9)
  aa=4.21d+01
  bb=9.115d0
  xe1=3.d0/2.d0
  xe2=1.0d0
  cinq=primus(aa,bb,xe1,xe2,t9)
  dcinq=drimus(bb,xe1,xe2,t9)
  aa=3.2d+01
  bb=9.383d0
  xe1=2.d0/3.d0
  xe2=1.0d0
  six=primus(aa,bb,xe1,xe2,t9)
  dsix=drimus(bb,xe1,xe2,t9)
  gt9=1.d0+5.d0*exp(-18.960d0/t9)
  dgt9=5.d0*exp(-18.960d0/t9)*18.960d0/t9
  centa=uno*due+tre+quat+0.1d0*(cinq+six)
  cent=centa/gt9
  dcent=(uno*duno*due+uno*ddue+tre*dtre+quat*dquat+0.1d0*(cinq*dcinq+six*dsix))/gt9-centa*dgt9/(gt9*gt9)
  e20agw=exp(rh1+fy20)*cent
  eps20(j1)=9.312d0*convMeVAvo/(20.d0*4.d0)*xne20(j1)*cya*e20agw
  e20t1=rht1+fy20*fyt20+dcent/cent
  e20p1=rhp1+fy20*fyp20

  if (nbouc /= 2) then
! L'abondance des particules alpha est supposee a l'equilibre
    cya=cya*(epcne(j1)/ecne+0.1739d0*ep23/e23)/(wpsyo/eo+eps20(j1)/e20+wpsyc/ec)
    nbouc=2
    goto 901
  endif

  epsc(j1)=epcne(j1)+epcna(j1)+wpsyo+ep23+eps20(j1)+wpsyc
  if (epsc(j1) /= 0.d0) then
    epsp1=(ecp12*(epcne(j1)+epcna(j1))+eop*wpsyo+e23p1*ep23+ecp*wpsyc+e20p1*eps20(j1))/epsc (j1)
    epst1=(ect12*(epcne(j1)+epcna(j1))+eot*wpsyo+e23t1*ep23+ect*wpsyc+e20t1*eps20(j1))/epsc (j1)
  else
    epsp1 = 0.d0
    epst1 = 0.d0
  endif
  if (epsc(j1) > 0.d0 .and. epsc(j) > 0.d0) then
    en=sqrt(epsc(j1)*epsc(j))
  else
    en = 0.d0
  endif

  if (j1 >= m) then
    write(3,*) 'energy C-burning'
    write(3,'(1x,i4,2x,9(1x,1pe11.3))') j1,cya,p(j1),t9,epsc(j1),epst1,epsp1
    write(3,*)  epcna(j1)+ep23,epcne(j1),wpsyc,wpsyo,eps20(j1)
    write(3,*) 'dE/dT*(T/E) C-burning'
    write(3,'(1x,i4,2x,9(1x,1pe11.3))') j1,e23t1,ect12,ect,eot,e20t1
    write(3,*) 'dE/dp*(T/E) C-burning'
    write(3,'(1x,i4,2x,9(1x,1pe11.3))')j1,e23p1,ecp12,ecp,eop,e20p1
    write(3,*) 'screening C-burning: ', fy12,fyo,fy20,fyc
  endif

! use 4.d0*yab(1) instead of y(j1)
  y(j1) = 4.d0* yab(1)
  call calcrates(j1,m,t9,exp(rh1),x(j1),y3(j1),y(j1),xc12(j1),xo16(j1),xne20(j1),xmg24(j1),rh1,rhpsi,rhpsit,rhp1,rht1,&
          vmyo,vmye,rhpsip,xc13(j1),xn14(j1),xn15(j1),xo17(j1),xo18(j1),xne22(j1),xmg25(j1),xmg26(j1),etot,etott,etotp)
  epsc(j1)= etot
  epsp1= etotp
  epst1= etott
  en=sqrt(abs(epsc(j1)*epsc(j)))

! NEUTRINOS,ITOH ET AL. (1989) APJ,339,354
!     same values as ITOH ET AL. (1996) APJS,102,411
24 continue
  ppp=exp(p(j1))
  ttt=exp(t(j1))
  xlttt=log10(ttt)
  xlam=xkmc2*ttt
  uxlam=1.d0/xlam
  uxlam2=uxlam*uxlam
  uxlam3=uxlam2*uxlam
  uxlam4=uxlam3*uxlam
  rho=exp(rh1)
  rhom=rho/vmye
  rhos9=1.0d-09*rhom
  xxi=uxlam*(rhos9**(1.d0/3.d0))

! DETERMINATION DE AO,A1,A2,A3,C,B1,B2,B3
! PAIR
  if (xlttt < 10.d0) then
   cnao(1)=6.002d+19
   cna1(1)=2.084d+20
   cna2(1)=1.872d+21
   cnc(1) =5.5924d0
   cnb1(1)=9.383d-01
   cnb2(1)=-4.141d-01
   cnb3(1)=5.829d-02
  else
   cnao(1)=6.002d+19
   cna1(1)=2.084d+20
   cna2(1)=1.872d+21
   cnc(1) =4.9924d0
   cnb1(1)=1.2383d0
   cnb2(1)=-0.8141d0
   cnb3(1)=0.0d0
  endif
! PHOTO
  if (xlttt >= 7.d0 .and. xlttt < 8.d0) then
    taup=xlttt-7.d0
    cnao(2)=0.5d0*cojl(1)+0.5d0*cojl(7)*cos(10.d0*pi*taup)
    do nj=1,5
     njj=nj+1
     cnao(2)=cnao(2)+cojl(njj)*cos(5.d0/3.d0*pi*nj*taup)+dojl(nj)*sin(5.d0/3.d0*pi*nj*taup)
    enddo
    cna1(2)=0.5d0*c1jl(1)+0.5d0*c1jl(7)*cos(10.d0*pi*taup)
    do nj=1,5
     njj=nj+1
     cna1(2)=cna1(2)+c1jl(njj)*cos(5.d0/3.d0*pi*nj*taup)+d1jl(nj)*sin(5.d0/3.d0*pi*nj*taup)
    enddo
    cna2(2)=0.5d0*c2jl(1)+0.5d0*c2jl(7)*cos(10.d0*pi*taup)
    do nj=1,5
     njj=nj+1
     cna2(2)=cna2(2)+c2jl(njj)*cos(5.d0/3.d0*pi*nj*taup)+d2jl(nj)*sin(5.d0/3.d0*pi*nj*taup)
    enddo
    cnc(2)=0.5654d0+taup
  elseif (xlttt >= 8.d0 .and. xlttt < 9.d0) then
    taup=xlttt-8.d0
    cnao(2)=0.5d0*coji(1)+0.5d0*coji(7)*cos(10.d0*pi*taup)
    do nj=1,5
     njj=nj+1
     cnao(2)=cnao(2)+coji(njj)*cos(5.d0/3.d0*pi*nj*taup)+doji(nj)*sin(5.d0/3.d0*pi*nj*taup)
    enddo
    cna1(2)=0.5d0*c1ji(1)+0.5d0*c1ji(7)*cos(10.d0*pi*taup)
    do nj=1,5
     njj=nj+1
     cna1(2)=cna1(2)+c1ji(njj)*cos(5.d0/3.d0*pi*nj*taup)+d1ji(nj)*sin(5.d0/3.d0*pi*nj*taup)
    enddo
    cna2(2)=0.5d0*c2ji(1)+0.5d0*c2ji(7)*cos(10.d0*pi*taup)
    do nj=1,5
     njj=nj+1
     cna2(2)=cna2(2)+c2ji(njj)*cos(5.d0/3.d0*pi*nj*taup)+d2ji(nj)*sin(5.d0/3.d0*pi*nj*taup)
    enddo
    cnc(2)=1.5654d0
  elseif (xlttt >= 9.d0) then
    taup=xlttt-9.d0
    cnao(2)=0.5d0*cojh(1)+0.5d0*cojh(7)*cos(10.d0*pi*taup)
    do nj=1,5
     njj=nj+1
     cnao(2)=cnao(2)+cojh(njj)*cos(5.d0/3.d0*pi*nj*taup)+dojh(nj)*sin(5.d0/3.d0*pi*nj*taup)
    enddo
    cna1(2)=0.5d0*c1jh(1)+0.5d0*c1jh(7)*cos(10.d0*pi*taup)
    do nj=1,5
     njj=nj+1
     cna1(2)=cna1(2)+c1jh(njj)*cos(5.d0/3.d0*pi*nj*taup)+d1jh(nj)*sin(5.d0/3.d0*pi*nj*taup)
    enddo
    cna2(2)=0.5d0*c2jh(1)+0.5d0*c2jh(7)*cos(10.d0*pi*taup)
    do nj=1,5
     njj=nj+1
     cna2(2)=cna2(2)+c2jh(njj)*cos(5.d0/3.d0*pi*nj*taup)+d2jh(nj)*sin(5.d0/3.d0*pi*nj*taup)
    enddo
    cnc(2)=1.5654d0
  endif
  cnb1(2)=6.290d-03
  cnb2(2)=7.483d-03
  cnb3(2)=3.061d-04
! PLASMA
  if (xlttt < 7.8d0) then
    cnao(3)=2.787d-07*xlttt*xlttt-3.936d-06*xlttt+1.408d-05
    cna1(3)=5.408d-07*xlttt*xlttt-7.715d-06*xlttt+2.751d-05
    cna2(3)=-1.812d-07*xlttt*xlttt+2.731d-06*xlttt-1.024d-05
    cna3=1.639d-08*xlttt*xlttt-2.445d-07*xlttt+9.137d-07
    cnc(3)=-7.683d-01+1.798d-01*xlttt
    cnb1(3)=1.148d-02
    cnb2(3)=1.209d-01
    cnb3(3)=2.432d-04
  else
    cnao(3)=2.320d-07
    cna1(3)=8.449d-08
    cna2(3)=1.787d-08
    cna3=0.d0
    cnc(3)=0.56457d0
    cnb1(3)=2.581d-02
    cnb2(3)=1.734d-02
    cnb3(3)=6.990d-04
  endif
! PAIR  PHOTO  PLASMA
  do kp=1,3
   wn(kp)=(cnao(kp)+cna1(kp)*xxi+cna2(kp)*xxi*xxi)*exp(-cnc(kp)*xxi)
   if(kp == 3) then
     wn(kp)=wn(kp)+cna3*xxi*xxi*xxi*exp(-cnc(kp)*xxi)
   endif
   wd(kp)=(xxi*xxi*xxi+cnb1(kp)*uxlam+cnb2(kp)*uxlam2+cnb3(kp)*uxlam3)
   wf(kp)=wn(kp)/wd(kp)
  enddo

! CALCUL DU TAUX DE PERTE D ENERGIE PAR UNITE DE VOLUME ET DE TEMPS
! QPAIR, QPHOT, QPLAS

! QPAIR
  xlam2=xlam*xlam
  xlam3=xlam2*xlam
  xlam4=xlam2*xlam2
  xlam5=xlam4*xlam
  xlam6=xlam4*xlam2
  xlam7=xlam6*xlam
  xlam8=xlam4*xlam4
  glam=1.d0-13.04d0*xlam2+133.5d0*xlam4+1534.d0*xlam6+918.6d0*xlam8
  dglam=-2.d0*13.04d0*xlam+4.d0*133.5d0*xlam3+6.d0*1534.d0*xlam5+8.d0*918.6d0*xlam7
  rhom2=rhom*rhom
  rhom3=rhom2*rhom
  e2lam=exp(-2.d0/xlam)
  qpai1=10.748d0*xlam2+0.3967d0*sqrt(xlam)+1.005d0
  qpai2=7.692d+07*xlam3+9.715d+06*sqrt(xlam)
  qpai3=1.d0+rhom/qpai2
  qpai=1.d0/(qpai1*(qpai3**0.3d0))
  qpabps=glam*e2lam*wf(1)
  cplus=(cv2+ca2)+neutri*(cvp2+cap2)
  cmoins=(cv2-ca2)+neutri*(cvp2-cap2)
  qpair=0.5d0*cplus*(1.d0+(cmoins/cplus)*qpai)*qpabps
! QPHOT
  qpho3=1.d0+2.045d0*xlam
  qpho4=1.875d+08*xlam+1.653d+08*xlam2+8.499d+08*xlam3-1.604d+08*xlam4
  qpho5=1.d0+rhom/qpho4
  qpho=0.666d0/((qpho3**2.066d0)*qpho5)
  qpho6=wf(2)*rhom*xlam5
  qphot=0.5d0*cplus*(1.d0-(cmoins/cplus)*qpho)*qpho6
! QPLAS
  qpla1=cv2+neutri*cvp2
  qplas=qpla1*rhom3*wf(3)
  qtot=qpair+qphot+qplas
  epsn1=qtot/rho
  if (epsn1 > HUGE(epsn1)) then
    write(*,*) 'j1,epsn1,qtot,rho:',j1,epsn1,qtot,rho
  endif
! DERIVEES
  rhotp=rho*rht1/ttt
  xxitp=(1.d0/3.d0)*(rhos9**(-(2.d0/3.d0)))*(1.d0/(vmye*1.d+09))*rhotp*uxlam-(rhos9**(1.d0/3.d0))*uxlam2*xkmc2
  do kp=1,3
   xntp=(cna1(kp)*xxitp+2.d0*cna2(kp)*xxi*xxitp)*exp(-cnc(kp)*xxi)
   if (kp == 3) then
     xntp=xntp+3.d0*cna3*xxi*xxi*xxitp*exp(-cnc(kp)*xxi)
   endif
   xntp=xntp-wn(kp)*cnc(kp)*xxitp
   xdtp=3.d0*xxi*xxi*xxitp-cnb1(kp)*uxlam2*xkmc2
   xdtp=xdtp-2.d0*cnb2(kp)*uxlam3*xkmc2-3.d0*cnb3(kp)*uxlam4*xkmc2
   wftp(kp)=(xntp*wd(kp)-wn(kp)*xdtp)/(wd(kp)*wd(kp))
  enddo
  dtpai1=(2.d0*10.748d0*xlam+0.5d0*0.3967d0/sqrt(xlam))*xkmc2
  dtpai2=(3.d0*7.692d+07*xlam2+0.5d0*9.715d+06/sqrt(xlam))*xkmc2
  dtpai3=rhotp/(vmye*qpai2)-rhom*dtpai2/(qpai2*qpai2)
  dtpai=-(qpai*qpai)*(dtpai1*(qpai3**0.3d0)+qpai1*0.3d0*dtpai3/(qpai3**0.7d0))
  dtpaib=dglam*xkmc2*e2lam*wf(1)+2.d0*glam*e2lam*uxlam2*xkmc2*wf(1)+glam*e2lam*wftp(1)
  dtpair=0.5d0*cplus*(dtpaib*(1.d0+(cmoins/cplus)*qpai)+qpabps*(cmoins/cplus)*dtpai)
  dtpho3=2.045d0*xkmc2
  dtpho4=(1.875d+08+2.d0*1.653d+08*xlam+3.d0*8.499d+08*xlam2-4.d0*1.604d+08*xlam3)*xkmc2
  dtpho5=rhotp/(vmye*qpho4)-rhom*dtpho4/(qpho4*qpho4)
  dtpho6=wftp(2)*rhom*xlam5+wf(2)*rhotp/vmye*xlam5+5.d0*wf(2)*rhom*xlam4*xkmc2
  dtpho=-(qpho*qpho/0.666d0)*(qpho3**1.066d0)*(2.066d0*dtpho3*qpho5+qpho3*dtpho5)
  dtphot=0.5d0*cplus*(dtpho6*(1.d0-(cmoins/cplus)*qpho)-qpho6*(cmoins/cplus)*dtpho)
  dtplas=qpla1*(3.d0*rhom2*rhotp/vmye*wf(3)+rhom3*wftp(3))
  dtqtot=dtpair+dtphot+dtplas
  dqtp=ttt*dtqtot/qtot
  enuet1=dqtp-rht1
  rhopt=rho*rhp1/ppp
  xxipt=(1.d0/3.d0)*1.d-09*(rhos9**(-(2.d0/3.d0)))*(1.d0/vmye)*rho*rhp1*(1.d0/ppp)*uxlam
  do kp=1,3
   xnpt=(cna1(kp)*xxipt+2.d0*cna2(kp)*xxi*xxipt)*exp(-cnc(kp)*xxi)-wn(kp)*cnc(kp)*xxipt
   if (kp == 3) then
     xnpt=xnpt+3.d0*cna3*xxi*xxi*xxipt*exp(-cnc(kp)*xxi)
   endif
   xdpt=3.d0*xxi*xxi*xxipt
   wfpt(kp)=(xnpt*wd(kp)-wn(kp)*xdpt)/(wd(kp)*wd(kp))
  enddo
  dppai3=rhopt/(vmye*qpai2)
  dppai=-(qpai*qpai)*qpai1*0.3d0*dppai3/(qpai3**0.7d0)
  dppaib=glam*e2lam*wfpt(1)
  dppair=0.5d0*cplus*(dppaib*(1.d0+(cmoins/cplus)*qpai)+qpabps*(cmoins/cplus)*dppai)
! DPPHO1=XIANT*XXIPT
  dppho5=rhopt/(vmye*qpho4)
  dppho=-(qpho*qpho/0.666d0)*(qpho3**2.066d0)*dppho5
  dppho6=wfpt(2)*rhom*xlam5+wf(2)*rhopt/vmye*xlam5
  dpphot=0.5d0*cplus*(dppho6*(1.d0-(cmoins/cplus)*qpho)-qpho6*(cmoins/cplus)*dppho)
  dpplas=qpla1*(3.d0*rhom2*rhopt/vmye*wf(3)+rhom3*wfpt(3))
  dpqtot=dppair+dpphot+dpplas
  dqpt=ppp*dpqtot/qtot
  enuep1=dqpt-rhp1
!---  NEUTRINO-PAIR BREMSSTRAHLUNG FROM DEGENERATE ELECTRON GAS
!     (CF. FESTA AND RUDERMAN 1969) ALLOWING FOR CORRECTION FACTOR OF
!     DICUS ET AL (1976) TO TAKE INTO ACCOUNT NEUTRAL CURRENT EFFECTS
!     PF AND EF ARE IN UNITS OF MC AND MC**2 RESP.
!---  EF INCLUDES THE ELECTRON REST-MASS ENERGY
  if ((p(j1)-2.5d0*t(j1)+3.680d0) > 0.d0) then
    pf=1.01d-2*rhe**(1.d0/3.d0)
    ef=sqrt(1.d0+pf*pf)
    if (ef-1.d0 > xlam) then
      bet_en=ef/pf
      al2=bet_en/215.d0
      bet2=bet_en*bet_en
      b2=log((2.d0+al2)/al2)
      b3=b2-2.d0/(2.d0+al2)
      b2=-4.d0+2.d0*b2*(1.d0+al2)
      aln=log((bet_en+1.d0)/(bet_en-1.d0))
      b1=(1.d0/3.d0)+(1.d0-bet2)*(0.5d0*bet_en*aln-1.d0)
      b3f=(2.d0/3.d0)+0.5d0*(bet2-0.5d0*bet_en*(bet2+1.d0)*aln)
      fb=0.14d0*(b1*b2-(bet2-1.d0)*b3*b3f)
!---  SZ2=2.*SUM OF ZI**2*XI/AI
!---  MAKE SURE TO ADD ABUNDANCES XNE20(J1),X23(J1),X24(J1),X28(J1) WHEN A
!---  MORE COMPLETE NETWORK IS INCLUDED FOR HE-,C-,NE-BURNINGS ETC.
!---  IN THE DEFINITION OF SZ2.
      sz2=2.d0*(x(j1)+y(j1)+3.d0*xc12(j1)+4.d0*xo16(j1)+3.5d0*xn14(j1)+6d0**2.d0*xc13(j1)/13.d0+7.d0**2.d0*xn15(j1)/15.d0+ &
          8.d0**2.d0*(xo17(j1)/17.d0+xo18(j1)/18.d0)+5.d0*xne20(j1)+10.d0**2.d0*xne22(j1)/22.d0+6.d0*xmg24(j1)+ &
          12.d0**2.d0*(xmg25(j1)/25.d0+xmg26(j1)/26.d0)+10.d0*z)

!  z elements taken ~ as Ni56: A=56 & Z=28
      sz2=sz2-2.d0*10.d0*z+2.d0*zabelx*56.d0**2.d0/28.d0
      do ii=1,nbelx
       sz2=sz2+2.d0*abelx(ii,j1)*nbzel(ii)**2.d0/nbael(ii)
      enddo

! Z EST CONSIDERE COMME  CA40 EN MOYENNE
      eb=(sz2*fb*t8**6.d0)*1.3d0
      ebr=-(1.d0/3.d0)*(bet2-1.d0)/(bet_en*fb)*(2.d0*b1*(b3-2.d0/(al2*(2.d0+al2)))/215.d0+ &
           b2*(3.d0*bet_en-0.5d0*(3.d0*bet2-1.d0)*aln)-2.d0*b3f*(bet_en*b3-2.d0*(bet2-1.d0)/(bet_en*(2.d0+al2)**2.d0))- &
           0.25d0*b3*(2.d0*bet_en*(3.d0*bet2-1.d0)-(3.d0*bet2+1.d0)*(bet2-1.d0)*aln))
      ebt=6.d0+rht1*ebr
      ebp=rhp1*ebr
    endif
  endif
  epsn1=-(epsn1+eb)
  enuet1=enuet1+ebt
  enuep1=enuep1+ebp
  if (j1 /= 1) then
    enue=sqrt(abs(epsn1*epsn))
    if (enue>HUGE(enue)) then
      write(*,*) 'enue,epsn,epsn1:',enue,epsn,epsn1
      stop
    endif
  endif
  eps_nu(j1) = epsn1

  return

end subroutine energ
!=======================================================================
subroutine screen(y,xc,xo,x20,x24,rh1,rhpsi,rhpsit,rhp1,rht1,vmyo,vmye,zw1,zw2,zw3,zw4,zw5,t8ln,rhpsip,fy,fyt,fyp, &
                  x13,x14,x15,x17,x18,x22,x25,x26,z)
!-----------------------------------------------------------------------
  implicit none

  real(kindreal),intent(in):: y,xc,xo,x20,x24,rh1,rhpsi,rhpsit,rhp1,rht1,vmyo,vmye,x13,x14,x15,x17,x18,x22,x25,x26, &
                              z,zw1,zw2,t8ln,rhpsip,zw3,zw4,zw5
  real(kindreal),intent(out):: fy,fyp,fyt

  real(kindreal):: zeta,ztild,xlamo,zbar,xl12,z3b1,eta,z13,xl23,zeb,sfy,zeb1,sfyp,sfyt
!-----------------------------------------------------------------------
! FACTEURS D ECRAN D APRES GRABOSKE APJ 181,457  (1973)
! LES TERMES ZW1-ZW5 SONT LES DIFFERENTS FACTEURS ZETA B
  zeta=sqrt(y+3.d0*xc+4.d0*xo+5.d0*x20+6.d0*x24+36.d0*x13/13.d0+3.5d0*x14+49.d0*x15/15.d0+64.d0*x17/17.d0+64.d0*x18/18.d0+ &
            100.d0*x22/22.d0+144.d0*x25/25.d0+144.d0*x26/26.d0+10.d0*z+0.5d0*rhpsi)
! ztild: rms charge average (Eq. (4) Dewitt et al. 1973 ApJ 181, 439)
  ztild=zeta*sqrt(vmyo)
! xlamo: Lambda_0 Eq. (19) Graboske et al. 1973 ApJ 181, 457
  xlamo=1.88d-04*exp(0.5d0*rh1-1.5d0*t8ln)/sqrt(vmyo)
! zbar: average ionic charge (Eq. (4) Dewitt et al. 1973 ApJ 181, 439)
  zbar=vmyo/vmye
  xl12=zw1*ztild*xlamo

  if (xl12 <= 0.1d0) then
! WEAK SCREENING
! fy: H_12(0) Eq. (19) Graboske et al. 1973 ApJ 181, 457
    fy=zw1*ztild*xlamo
    fyp=0.5d0*rhp1+0.25d0*rhpsip/(zeta**2.d0)
    fyt=0.5d0*rht1-1.5d0+0.25d0*rhpsit/(zeta**2.d0)
    return
  endif
  if (xl12 < 5.d0) then
! INTERMEDIATE SCREENING
! TOUS LES ELEMENTS LOURDS SONT CONSIDERES COMME X40
! z3b1: Z^(3b-1) selon Eq. (65) Dewitt et al. 1973 ApJ 181, 439
! avec b= 0.86 pour intermediate screening (Graboske et al. 1973 ApJ 181, 457)
    z3b1=(2.d0**1.58d0*y/4.d0+6.d0**1.58d0*(xc/12.d0+x13/13.d0)+7.d0**1.58d0*(x14/14.d0+x15/15.d0)+ &
          8.d0**1.58d0*(xo/16.d0+x17/17.d0+x18/18.d0)+10.d0**1.58d0*(x20/20.d0+x22/22.d0)+ &
          12.d0**1.58d0*(x24/24.d0+x25/25.d0+x26/26.d0)+20.d0**1.58d0*z/40.d0)*vmyo
! Ci-dessous: 0.58=3b-2 et 0.28=2-2b (cas b+0.86) Eq. (65) Dewitt et al. 1973 ApJ 181, 439
    eta=z3b1/((ztild**0.58d0)*(zbar**0.28d0))
    fy=0.38d0*eta*zw2*(xlamo**0.86d0)
    fyp=0.43d0*rhp1-0.145d0*rhpsip/(zeta*zeta)
    fyt=-1.290d0+0.43d0*rht1-0.145d0*rhpsit/(zeta*zeta)
    if (xl12 <= 2.d0) then
      return
    endif
  endif
! STRONG AND INTERMEDIATE STRONG SCREENING
  z13=zbar**(1.d0/3.d0)
  xl23=xlamo**(2.d0/3.d0)
  zeb=zw3+0.316d0*zw4*z13+0.737d0*zw5/(zbar*xl23)
  sfy=0.624d0*z13*zeb*xl23
  zeb1=zeb*zbar*xl23
  sfyp=-0.2457d0*zw5*rhp1/zeb1+(1.d0/3.d0)*rhp1
  sfyt=-0.4913d0*zw5*(0.5d0*rht1-1.5d0)/zeb1+(1.d0/3.d0)*rht1-1.d0
  if (xl12 >= 5.d0) then
! STRONG ONLY
    fy=sfy
    fyp=sfyp
    fyt=sfyt
    return
  endif
  if (fy <= sfy) then
    return
  endif
! INTERMEDIATE STRONG
  fy=sfy
  fyp=sfyp
  fyt=sfyt

  return

end subroutine screen
!=======================================================================
subroutine calcrates(j1,m,temp9,rh,xx,xy3,xy,xc,xo,x20,x24,rh1,rhpsi,rhpsit,rhp1,rht1,vmyo,vmye,rhpsip, &
                     xc13,xn14,xn15,xo17,xo18,x22,x25,x26,etot,etott,etotp)
!-----------------------------------------------------------------------
  use abundmod,only: nbelx,nbzel,nbael,abelx,zabelx,eps_c_adv,eps_ne_adv,eps_o_adv,eps_si_adv

  implicit none

  integer, parameter:: nn=5
  integer:: iqse,j1,m,i,ii,j,klo,khi,k,jlo,jhi,mm,ks

  real(kindreal):: temp9,rh,t8,t8ln,zeta
  real(kindreal):: bb,h,anegridhigh,deltagrid
  real(kindreal):: v,vbetalow,vbetahigh
  real(kindreal):: rh1,rhpsi,rhpsit,rhp1,rht1,vmyo,vmye,rhpsip
  real(kindreal):: xx,xy3,xy,xc,xo,x20,x24
  real(kindreal):: xc13,xn14,xn15,xo17,xo18,x22,x25,x26
  real(kindreal):: ztild,xlamo,zbar,xl12
  real(kindreal):: zw1,zw2,zw3,zw4,zw5,b,z3b1,eta,z13,xl23,zeb,zeb1
  real(kindreal):: fy,fyp,fyt, sfy,sfyp,sfyt,dedt,ref
  real(kindreal):: eprod,eprodt,eprodp,etot,etott,etotp,e2e
  real(kindreal), dimension(nn):: logt,logrr,coef
  real(kindreal), dimension(0:4):: f
!-----------------------------------------------------------------------
! energy production [MeV/mH] --> [erg/g]
  e2e = cst_avo*convMeVerg
  etot=0.d0
  eprod=0.d0
  etott=0.d0
  eprodt=0.d0
  etotp=0.d0
  eprodp=0.d0
  eps_c_adv(j1)=0.d0
  eps_ne_adv(j1)=0.d0
  eps_o_adv(j1)=0.d0
  eps_si_adv(j1)=0.d0

!  initialisation
! temperature in 10^8 [K], density in [g/cm^3]
  t8= temp9*10.d0
  rho= rh
! abundances
  do i=1,nbel
   abuny(i)=0.0d0
  enddo

  if (posel( 1, 1- 1) > 0) abuny(posel( 1, 1- 1)) = xx
  if (posel( 2, 3- 2) > 0) abuny(posel( 2, 3- 2)) = xy3 / 3.d0
  if (posel( 2, 4- 2) > 0) abuny(posel( 2, 4- 2)) = xy  / 4.d0
  if (posel( 6,12- 6) > 0) abuny(posel( 6,12- 6)) = xc  /12.d0
  if (posel( 6,13- 6) > 0) abuny(posel( 6,13- 6)) = xc13/13.d0
  if (posel( 7,14- 7) > 0) abuny(posel( 7,14- 7)) = xn14/14.d0
  if (posel( 7,15- 7) > 0) abuny(posel( 7,15- 7)) = xn15/15.d0
  if (posel( 8,16- 8) > 0) abuny(posel( 8,16- 8)) = xo  /16.d0
  if (posel( 8,17- 8) > 0) abuny(posel( 8,17- 8)) = xo17/17.d0
  if (posel( 8,18- 8) > 0) abuny(posel( 8,18- 8)) = xo18/18.d0
  if (posel(10,20-10) > 0) abuny(posel(10,20-10)) = x20 /20.d0
  if (posel(10,22-10) > 0) abuny(posel(10,22-10)) = x22 /22.d0
  if (posel(12,24-12) > 0) abuny(posel(12,24-12)) = x24 /24.d0
  if (posel(12,25-12) > 0) abuny(posel(12,25-12)) = x25 /25.d0
  if (posel(12,26-12) > 0) abuny(posel(12,26-12)) = x26 /26.d0

  do ii=1,nbelx
   if (posel(nbzel(ii),nbael(ii)-nbzel(ii)) > 0) abuny(posel(nbzel(ii),nbael(ii)-nbzel(ii))) = abelx(ii,j1)/nbael(ii)
  enddo

! Z, PME
!  calculation of PME & ANE
!     ane is the electron density in units of 1E26 cm-3
!     PME is the mean molecular weight per electron
!         =  [ SUM (XZ/A) ]**(-1)
  pme=xx+xy3*2.d0/3.d0+xy*2.d0/4.d0+xc*6.d0/12.d0+ xc13*6.d0/13.d0+xn14*7.d0/14.d0+xn15*7.d0/15.d0+xo*8.d0/16.d0+ &
    xo17*8.d0/17.d0+xo18*8.d0/18.d0+x20*10.d0/20.d0+x22*10.d0/22.d0+x24*12.d0/24.d0+x25*12.d0/25.d0+x26*12.d0/26.d0+0.5d0*zabelx
  do ii=1,nbelx
   pme = pme + nbzel(ii)*abelx(ii,j1)/nbael(ii)
  enddo
  pme =1.d0/pme
  ane   = rho * cst_avo / (1.d26 * pme)

! electron screening - NB: all heavy elements considered as Ni56
  t8ln=log(t8)
  zeta=xx+4.d0/3.d0*xy3+xy+xc*36.d0/12.d0+xc13*36.d0/13.d0+xn14*49.d0/14.d0+xn15*49.d0/15.d0+xo*64.d0/16.d0+xo17*64.d0/17.d0+ &
       xo18*64.d0/18.d0+x20*100.d0/20.d0+x22*100.d0/22.d0+x24*144.d0/24.d0+x25*144.d0/25.d0+x26*144.d0/26.d0+784.d0/56.d0*zabelx

  do ii=1,nbelx
   zeta = zeta + nbzel(ii)**2.d0*abelx(ii,j1)/nbael(ii)
  enddo
! rhpsi for degeneracy contribution
  zeta= zeta + 0.5d0*rhpsi
  zeta=sqrt(zeta)
  ztild =zeta*sqrt(vmyo)
  xlamo=1.88d-04*exp(0.5d0*rh1-1.5d0*t8ln)/sqrt(vmyo)
  zbar=vmyo/vmye

! factorials
  f(0) = 1.d0
  f(1) = 1.d0
  f(2) = 2.d0
  f(3) = 6.d0
  f(4) = 24.d0

! interpolate reaction rate

! locate the position of current T within the grid
! and compute step, a, b, c and d parameters for cubic spline interpolation
  klo = 1
  khi = kgrid
100 if (khi-klo > 1) then
    k = (khi+klo)/2
    if (tgrid(k) > T8) then
      khi=k
    else
      klo=k
    endif
    goto 100
  endif
  h=tgrid(khi)-tgrid(klo)
  bb=(T8 - tgrid(klo))/h

! locate the position of current rho within the density grid (1.,3.,10.,30.)
  if (ane <= 3.d0) then
    anegridhigh=3.d0
    deltagrid=2.d0
    jlo=1
    jhi=2
  else if (ane <= 10.d0) then
    anegridhigh=10.d0
    deltagrid=7.d0
    jlo=2
    jhi=3
  else
    anegridhigh=30.d0
    deltagrid=20.d0
    jlo=3
    jhi=4
  endif

! perform linear interpolation of log(v)
  do i=1,ireac
   eprod=0.d0
   eprodt=0.d0
   eprodp=0.d0

   if (flag(i) <= 0.d0) then
! reaction is not electron-density-dependent beta-decay
     v= log10(vgrid(klo,i,1))+ bb * (log10(vgrid(khi,i,1)) - log10(vgrid(klo,i,1)))

! dlnE/dlnT calculation
     if (bb > 0.5d0) then
       ks=klo-(nn/2-1)
     else
       ks=klo-nn/2
     endif
     if (ks < 1) then
       ks=1
     endif
     if (ks > kgrid-nn+1) then
       ks=kgrid-nn+1
     endif
     do j=1,nn
      logt(j)=log10(tgrid(ks+j-1)/10.d0)
      logrr(j)=log10(vgrid(ks+j-1,i,1))
      coef(j)=0.0d0
     enddo
     mm=3
     ref=0.0d0
     call e02acf(logt,logrr,nn,coef,mm,ref)
     dedt=coef(2)+2.d0*coef(3)*log10(temp9)

! test the reaction kind
     if (flag(i) == -14.d0) then
!    electron capture on Be7
       rrate(i,j1) = 10.d0**v * RHO / PME
       if (rrate(i,j1) > 1.51d-7.and.T8 < 0.01d0) then
         rrate(i,j1) = 1.51d-7
       endif

! en. prod. = e2e*Y1*[1e]*Qreac
       eprod= e2e*abuny(elps(i,1))*rrate(i,j1)*qrad(i)
     else if (flag(i) == -13.d0) then
!     electron capture
       rrate(i,j1) = 10.d0**v *RHO / PME

! en. prod. = e2e*Y1*[1e]*Qreac
       eprod= e2e*abuny(elps(i,1))*rrate(i,j1)*qrad(i)
       eprodt = rht1 + dedt
       eprodp = rhp1
     else if (flag(i) == -11.d0) then
! photodisintegration or beta-decay
       rrate(i,j1) = 10.d0**v

! en. prod. = e2e*Y1*[1]*Qreac
       eprod= e2e*abuny(elps(i,1))*rrate(i,j1)*qrad(i)
       eprodt =  dedt
     else if (flag(i) == -10.d0) then
! two-particle reaction   !if identical particles: factorials!
       rrate(i,j1) = 10.d0**v *RHO /f(nsnb(i,1))/f(nsnb(i,2))

! screening:    cf 1973PaJ...181..457G by Graboske, DeWitt, ... p.465
       zw1= 2.d0*znb(i,1)**nsnb(i,1)*znb(i,2)**nsnb(i,2)
       xl12=0.5d0*zw1*ztild*xlamo
! WEAK SCREENING     b=1,kb=0.5,
       if (xl12 <= 0.1d0) then
         fy=0.5d0*zw1*ztild*xlamo
         fyp=0.5d0*rhp1+0.25d0*rhpsip/(zeta**2.d0)
         fyt=0.5d0*rht1-1.5d0+0.25d0*rhpsit/(zeta**2.d0)
! INTERMEDIATE SCREENING  b=0.860,kb=0.38,
       else if (xl12 <= 5.d0) then
         b=0.860d0
         zw2= (znb(i,1)*nsnb(i,1)+ znb(i,2)*nsnb(i,2))**(1.d0+b)-nsnb(i,1)*znb(i,1)**(1.d0+b)-nsnb(i,2)*znb(i,2)**(1.d0+b)
         z3b1=xx+2.d0**1.58d0*xy/4.d0+6.d0**1.58d0*(xc/12.d0+xc13/13.d0)+7.d0**1.58d0*(xn14/14.d0+xn15/15.d0)+ &
              8.d0**1.58d0*(xo/16.d0+xo17/17.d0+xo18/18.d0)+10.d0**1.58d0*(x20/20.d0+x22/22.d0)+12.d0**1.58d0*(x24/24.d0+ &
              x25/25.d0+x26/26.d0)+28.d0**1.58d0*zabelx/56.d0
! all heavy elements considered as Ni56 Ai:56 Zi:28 2.842*z-->3.454*zabelx
         do ii=1,nbelx
          z3b1= z3b1+ nbzel(ii)**(3.d0*b-1.d0)*abelx(ii,j1)/nbael(ii)
         enddo
         z3b1=z3b1*vmyo
         eta=z3b1/((ztild**(3.d0*b-2.d0))*(zbar**(2.d0-2.d0*b)))
         fy=0.38d0*eta*zw2*(xlamo**b)
         fyp=0.43d0*rhp1-0.145d0*rhpsip/(zeta*zeta)
         fyt=-1.290d0+0.43d0*rht1-0.145d0*rhpsit/(zeta*zeta)
! INTERMEDIATE STRONG SCREENING    kb=0.624
         if (xl12 >= 2.d0) then
           b=2.d0/3.d0
           z13=zbar**(1.d0/3.d0)
           xl23=xlamo**b
           zw3= (znb(i,1)*nsnb(i,1)+ znb(i,2)*nsnb(i,2))**(5.d0/3.d0)-nsnb(i,1)*znb(i,1)**(5.d0/3.d0)- &
                 nsnb(i,2)*znb(i,2)**(5.d0/3.d0)
           zw4= (znb(i,1)*nsnb(i,1)+ znb(i,2)*nsnb(i,2))**(4.d0/3.d0)-nsnb(i,1)*znb(i,1)**(4.d0/3.d0)- &
                 nsnb(i,2)*znb(i,2)**(4.d0/3.d0)
           zw5= (znb(i,1)*nsnb(i,1)+ znb(i,2)*nsnb(i,2))**(2.d0/3.d0)-nsnb(i,1)*znb(i,1)**(2.d0/3.d0)- &
                 nsnb(i,2)*znb(i,2)**(2.d0/3.d0)

           zeb=zw3+0.316d0*zw4*z13+0.737d0*zw5/(zbar*xl23)
           sfy=0.624d0*z13*zeb*xl23
           zeb1=zeb*zbar*xl23
           sfyp=-0.2457d0*zw5*rhp1/zeb1+0.3333d0*rhp1
           sfyt=-0.4913d0*zw5*(0.5d0*rht1-1.5d0)/zeb1+(1.d0/3.d0)*rht1-1.d0
! INTERMEDIATE STRONG
           if (fy > sfy) then
             fy=sfy
             fyp=sfyp
             fyt=sfyt
           endif
         endif
! STRONG ONLY  (xl12 >= 5.)     kb=0.624
       else
         b=2.d0/3.d0
         z13=zbar**(1.d0/3.d0)
         xl23=xlamo**b
         zw3= (znb(i,1)*nsnb(i,1)+ znb(i,2)*nsnb(i,2))**(5.d0/3.d0)-nsnb(i,1)*znb(i,1)**(5.d0/3.d0)- &
               nsnb(i,2)*znb(i,2)**(5.d0/3.d0)
         zw4= (znb(i,1)*nsnb(i,1)+ znb(i,2)*nsnb(i,2))**(4.d0/3.d0)-nsnb(i,1)*znb(i,1)**(4.d0/3.d0)- &
               nsnb(i,2)*znb(i,2)**(4.d0/3.d0)
         zw5= (znb(i,1)*nsnb(i,1)+ znb(i,2)*nsnb(i,2))**(2.d0/3.d0)-nsnb(i,1)*znb(i,1)**(2.d0/3.d0)- &
               nsnb(i,2)*znb(i,2)**(2.d0/3.d0)

         zeb=zw3+0.316d0*zw4*z13+0.737d0*zw5/(zbar*xl23)
         fy=0.624d0*z13*zeb*xl23
         zeb1=zeb*zbar*xl23
         fyp=-0.2457d0*zw5*rhp1/zeb1+0.3333d0*rhp1
         fyt=-0.4913d0*zw5*(0.5d0*rht1-1.5d0)/zeb1+0.3333d0*rht1-1.d0

       endif
! end screening;

       rrate(i,j1) = rrate(i,j1)*exp(fy)
! en. prod. = e2e*(1/(n1!*n2!)*[12])*Y1*Y2*Qreac

!        2 a --> ... : nsnb(i,1)=2, nsnb(i,2)=0
       if (elps(i,2) == 0.or.nsnb(i,2) == 0) then
         eprod = e2e*abuny(elps(i,1))**2.d0*rrate(i,j1)*qrad(i)
       else
!      a + b --> ... : nsnb(i,1)=1, nsnb(i,2)=1
         eprod = e2e*abuny(elps(i,1))*abuny(elps(i,2))*rrate(i,j1)*qrad(i)
       endif
       eprodt = rht1+ fyt*fy+ dedt
       eprodp = rhp1+ fyp*fy
     else if (flag(i) == -100.d0) then
! three-particle reactions   !if identical particles: factorials!
       rrate(i,j1) = 10.d0**v *RHO**2.d0 /f(nsnb(i,1))/f(nsnb(i,2))
       eprod = e2e*abuny(elps(i,1))**3.d0*rrate(i,j1)*qrad(i)
     else if (flag(i) == -200.d0) then
! four-particle reactions   !if identical particles: factorials!
       rrate(i,j1) = 10.d0**v *RHO**3.d0 /f(nsnb(i,1))/f(nsnb(i,2))
     endif

   else

! beta-decay rate has to be interpolated in density

! linear interpolation in temperature for two density grid points
!                (note: log of beta decay rate is handled)

     vbetalow =log10(vgrid(klo,i,jlo))+ bb*(log10(vgrid(khi,i,jlo))-log10(vgrid(klo,i,jlo)))
     vbetahigh=log10(vgrid(klo,i,jhi))+ bb*(log10(vgrid(khi,i,jhi))-log10(vgrid(klo,i,jhi)))

     vbetalow    = 10.d0**vbetalow
     vbetahigh   = 10.d0**vbetahigh

! linear interpolation in density
     rrate(i,j1)=(vbetahigh-vbetalow)/deltagrid*(ane-anegridhigh)+vbetahigh

! extrapolation often leads to negative rates
! check if so, and then put it to zero
! (supplementary table values would be needed for those cases)
     if (rrate(i,j1) < 0.d0) then
       rrate(i,j1) = 0.d0
     endif
   endif

   if (elps(i,1) == posel(10,20-10).and.elps(i,4) == posel(8,8)) then
     taunucl(j1)=1.d0/(abs(rrate(i,j1))+1.d-50)
   endif

   if (phase >= 6.and.xo < 0.1d0) then
     iqse=1
   else
     iqse=0
   endif
   if (iqse == 1) then
     if (elps(i,1)==posel(22,22) .and. elps(i,4)==posel(24,24) .or. elps(i,1)==posel(24,24) .and. elps(i,4)==posel(22,22)) then
! 49.383d0=Q(Si->Ni) & 7.692d0=Q(Ti->Cr)
       etot=etot + eprod*49.383d0/7.692d0
       etott=etott + eprodt*eprod*49.383d0/7.692d0
       etotp=etotp + eprodp*eprod*49.383d0/7.692d0
       eps_si_adv(j1)=eps_si_adv(j1)+abs(eprod)*49.383d0/7.692d0
! count energy from elements lighter than Si
     elseif (nbz(elps(i,1)) <= 14 .and. nbz(elps(i,4)) <= 14) then
       etot=etot + eprod
       etott=etott + eprodt*eprod
       etotp=etotp + eprodp*eprod
     endif
   else
     etot=etot + eprod
     etott=etott + eprodt*eprod
     etotp=etotp + eprodp*eprod
     if (reaction(i)(1:1) == '2' .and. reaction(i)(6:7) == '12') then
       eps_c_adv(j1)=eps_c_adv(j1)+eprod
     else if (reaction(i)(1:1) == '2' .and. reaction(i)(6:7) == '16') then
       eps_o_adv(j1)=eps_o_adv(j1)+eprod
     else if (reaction(i)(6:7) == '20' .and. reaction(i)(11:11) == '0') then
       eps_ne_adv(j1)=eps_ne_adv(j1)+abs(eprod)
     endif
   endif

   if (j1 >= m) then
     if (rrate(i,j1) <= 1.d-30) then
       write(3,'("energy prod.i,e,t,p,f,r: X",i3,5(1p,e12.5),a40)')i,eprod,eprodt,eprodp,fy,rrate(i,j1),reaction(i)
     else
       write(3,'("energy prod.i,e,t,p,f,r: ",i4,5(1p,e12.5),a40)')i,eprod,eprodt,eprodp,fy,rrate(i,j1),reaction(i)
     endif
   endif
  enddo
  etott=etott/etot
  etotp=etotp/etot

  if (j1 >= m) then
    write(3,'("energy prod. tot: ",i4,4(1p,e12.5))')j1,t8/10.d0,etot,etott,etotp
  endif

  return

end subroutine calcrates
!======================================================================
subroutine interpol(it,x,v)
!----------------------------------------------------------------------
  use evol,only: input_dir
  use inputparam,only: var_rates

  implicit none

  integer,intent(in)::it
  real(kindreal),intent(in)::x
  real(kindreal),intent(out)::v

  integer, parameter:: dim=64,nbtables=47
  integer, save:: first
  integer, dimension(dim), save:: nblignes
  real(kindreal), dimension(2,nbtables,dim), save:: fichier
  character(256), dimension(nbtables), save:: tables3

  integer:: j,jligne=0,nn,n,i,ns,ierror
  real(kindreal):: hh,h,dy,y
  real(kindreal), dimension(dim):: xa,ya
!----------------------------------------------------------------------
  xa(:) = 0.d0
  ya(:) = 0.d0

  if (first  ==  0) then
    if (var_rates) then
      open(33,file='liste_tables3',status='old')
    else
      open(33,file=trim(input_dir)//'inputs/liste_tables3',status='old')
    endif
    ierror = 0
    do i=1,nbtables
     read(33,'(a64)', iostat=ierror) tables3(i)
     if (ierror /= 0) then
       exit
     endif
    enddo
    close(33)

    do i=1,nbtables
     if (var_Rates) then
       open(34,file=tables3(i),status='old',form='formatted')
     else
       open(34,file=trim(input_dir)//'taux/'//tables3(i),status='old',form='formatted')
     endif
     ierror = 0
     do j=1,dim
      read(34,'(f7.3,d13.2)',iostat=ierror) fichier(1,i,j),fichier(2,i,j)
      if (ierror /= 0) then
        exit
      endif
      jligne=j
     enddo
     nblignes(i) = jligne
    enddo
    close(34)
  endif

  first = 1

  do i=1,nblignes(it)
   xa(i)=fichier(1,it,i)
   ya(i)=log(fichier(2,it,i))  ! to read rate
  enddo
  nn=nblignes(it)
  n=nblignes(it)
  ns=1
  hh=abs(x-xa(1))
  do i=1,n
   h=abs(x-xa(i))
   if (h == 0d0) then
     v=exp(ya(i))
     dy=0.0d0
     return
   else if (h < hh) then
     ns=i
     hh=h
   endif
  enddo
  y=ya(ns)
  if (ns == nn) then
    dy=((ya(ns)-ya(ns-1))/(xa(ns)-xa(ns-1)))*(xa(ns)-x)
    y=y-dy
  else
    if (ns == 1) then
      if (x < xa(1)) then
        v=0.d0
        return
      else
        dy=((ya(ns+1)-ya(ns))/(xa(ns+1)-xa(ns)))*(x-xa(ns))
        y=y+dy
      endif
    else
      if (x > xa(ns)) then
        dy=((ya(ns+1)-ya(ns))/(xa(ns+1)-xa(ns)))*(x-xa(ns))
        y=y+dy
      else
        dy=((ya(ns)-ya(ns-1))/(xa(ns)-xa(ns-1)))*(xa(ns)-x)
        y=y-dy
      endif
    endif
  endif

  v=exp(y)

  return

end subroutine interpol
!======================================================================
subroutine netburning(l,temp9,ddeit,vxab,onetwo)
!----------------------------------------------------------------------
  use abundmod,only: fnucdif,vvabelx
  use inputparam,only: idifcon

  implicit none

  integer:: i,ii,j,l,onetwo
  real(kindreal):: temp9,ddeit,sumdxdt
  real(kindreal), dimension(15):: vxab
!----------------------------------------------------------------------
  do i=1,nbel
   abuny(i)=0.0d0
  enddo

  if (posel( 1, 1- 1) > 0) abuny(posel( 1, 1- 1)) = vxab(1)

  if (posel( 2, 3- 2) > 0) abuny(posel( 2, 3- 2)) = vxab(2) / 3.d0
  if (posel( 2, 4- 2) > 0) abuny(posel( 2, 4- 2)) = vxab(3) / 4.d0

  if (posel( 6,12- 6) > 0) abuny(posel( 6,12- 6)) = vxab(4) /12.d0
  if (posel( 6,13- 6) > 0) abuny(posel( 6,13- 6)) = vxab(5) /13.d0

  if (posel( 7,14- 7) > 0) abuny(posel( 7,14- 7)) = vxab(6) /14.d0
  if (posel( 7,15- 7) > 0) abuny(posel( 7,15- 7)) = vxab(7) /15.d0

  if (posel( 8,16- 8) > 0) abuny(posel( 8,16- 8)) = vxab(8) /16.d0
  if (posel( 8,17- 8) > 0) abuny(posel( 8,17- 8)) = vxab(9) /17.d0
  if (posel( 8,18- 8) > 0) abuny(posel( 8,18- 8)) = vxab(10)/18.d0

  if (posel(10,20-10) > 0) abuny(posel(10,20-10)) = vxab(11)/20.d0
  if (posel(10,22-10) > 0) abuny(posel(10,22-10)) = vxab(12)/22.d0

  if (posel(12,24-12) > 0) abuny(posel(12,24-12)) = vxab(13)/24.d0
  if (posel(12,25-12) > 0) abuny(posel(12,25-12)) = vxab(14)/25.d0
  if (posel(12,26-12) > 0) abuny(posel(12,26-12)) = vxab(15)/26.d0


  fnucdif =0.00d0
  if (phase >= 5 .and. idifcon == 1) then
    fnucdif=0.5d0
  endif
  if (onetwo == 1) then
    tmax=ddeit*(1.d0-fnucdif)
    do ii=1,nbelx
     if (posel(nbzel(ii),nbael(ii)-nbzel(ii)) > 0) abuny(posel(nbzel(ii),nbael(ii)-nbzel(ii))) = vvabelx(ii,l)/nbael(ii)
    enddo
  else if (onetwo == 2) then
    tmax=ddeit*fnucdif
    do ii=1,nbelx
     if (posel(nbzel(ii),nbael(ii)-nbzel(ii)) > 0) abuny(posel(nbzel(ii),nbael(ii)-nbzel(ii))) = abelx(ii,l)/nbael(ii)
    enddo
  else
    print*, 'stop in netburning: onetwo'
    stop
  endif
  if (tmax == 0.d0) return

! X & previous abundances
  do i=1,nbel
   abunx(i)=abuny(i)*nba(i)
   vabuny(i)=abuny(i)
   vabunx(i)=abunx(i)
  enddo

! temperature, density, time step
  t9= temp9

!  calculation of PME & ANE
!     ane is the electron density in units of 1E26 cm-3
!     PME is the mean molecular weight per electron
!         =  [ SUM (XZ/A) ]**(-1)
  pme=vxab(1)+ vxab(2)*2.d0/3.d0+vxab(3)*2.d0/4.d0+vxab(4)*6.d0/12.d0+vxab(5)*6.d0/13.d0+vxab(6)*7.d0/14.d0+ &
      vxab(7)*7.d0/15.d0+vxab(8)*8.d0/16.d0+vxab(9)*8.d0/17.d0+vxab(10)*8.d0/18.d0+vxab(11)*10.d0/20.d0+ &
      vxab(12)*10.d0/22.d0+vxab(13)*12.d0/24.d0+vxab(14)*12.d0/25.d0+vxab(15)*12.d0/26.d0+0.5d0*zabelx
  do ii=1,nbelx
   pme = pme + nbzel(ii)*abelx(ii,l)/nbael(ii)
  enddo
  pme =1.d0/pme
  ane   = rho * cst_avo/(1.d26*pme)

  if (itestx == 1 .and. verbose) then
    write(*,*) ' t9         rho        tmax       pme        ane  '
    write(*,'(5(e10.3,1p))') t9,rho,tmax,pme,ane
  endif

! start of calculations for matrix A, size: NBEL x (NBEL + 1)
  do i=1,nbel
   do j=1,nbel
    if (j /= i) then
      mata(i,j) =0.d0
    else
      mata(i,j)=-1.d0/tmax
    endif
   enddo
   mata(i,nbel+1)=-vabuny(i)/tmax
  enddo

  shellnb=l
  call contribreac

  sumdxdt=0.d0
  if (itestx == 0) then
    do i=1,nbel
     do j=1,nbel
      sumdxdt=sumdxdt+mata(i,j)*nba(i)*vabuny(j)
     enddo
     sumdxdt=sumdxdt-mata(i,nbel+1)*nba(i)
    enddo
    if (phase < 6) then
      if (abs(sumdxdt) > 1.d-15 .and. verbose) then
        write(*,*) l,': sumdxdt = ',sumdxdt
      endif
    else
      if (abs(sumdxdt) > 1.d8 .and. verbose) then
        write(*,*) l,': sumdxdt = ',sumdxdt
      endif
    endif
  endif

! inverse matrix A
  call inversemat(nbel,mata,abuny,maxel)
  do i=1,nbel
   abunx(i)=abuny(i)*nba(i)
  enddo

! output results
  if (itestx == 1 .and. verbose) then
   write(*,*) '  abuny  ,  vabuny'
   do i=1,nbel
    write(*,'(2(1p,1x,e12.5))')abuny(i) , vabuny(i)
   enddo
  endif

  if (posel( 1, 1- 1) > 0) vxab(1)  = abuny(posel( 1, 1- 1))

  if (posel( 2, 3- 2) > 0) vxab(2)  = abuny(posel( 2, 3- 2))* 3.d0
  if (posel( 2, 4- 2) > 0) vxab(3)  = abuny(posel( 2, 4- 2))* 4.d0

  if (posel( 6,12- 6) > 0) vxab(4)  = abuny(posel( 6,12- 6))*12.d0
  if (posel( 6,13- 6) > 0) vxab(5)  = abuny(posel( 6,13- 6))*13.d0

  if (posel( 7,14- 7) > 0) vxab(6)  = abuny(posel( 7,14- 7))*14.d0
  if (posel( 7,15- 7) > 0) vxab(7)  = abuny(posel( 7,15- 7))*15.d0

  if (posel( 8,16- 8) > 0) vxab(8)  = abuny(posel( 8,16- 8))*16.d0
  if (posel( 8,17- 8) > 0) vxab(9)  = abuny(posel( 8,17- 8))*17.d0
  if (posel( 8,18- 8) > 0) vxab(10) = abuny(posel( 8,18- 8))*18.d0

  if (posel(10,20-10) > 0) vxab(11) = abuny(posel(10,20-10))*20.d0
  if (posel(10,22-10) > 0) vxab(12) = abuny(posel(10,22-10))*22.d0

  if (posel(12,24-12) > 0) vxab(13) = abuny(posel(12,24-12))*24.d0
  if (posel(12,25-12) > 0) vxab(14) = abuny(posel(12,25-12))*25.d0
  if (posel(12,26-12) > 0) vxab(15) = abuny(posel(12,26-12))*26.d0

  do ii=1,nbelx
   if (posel(nbzel(ii),nbael(ii)-nbzel(ii)) > 0) abelx(ii,l) = nbael(ii)*abuny(posel(nbzel(ii),nbael(ii)-nbzel(ii)))
  enddo

  return

end subroutine netburning
!=======================================================================
subroutine netinit(z)
!----------------------------------------------------------------------
  use evol,only: input_dir
  use inputparam,only: idebug, iapprox21
  use abundmod,only: mbelx,abels,xlostneu

  implicit none

  integer:: i,ii,ierror
  real(kindreal),intent(in):: z
  character(256):: vit_fileCNE, vit_fileCNEO, netinit_fileCNE, netinit_fileCNEO
!----------------------------------------------------------------------
! Reading network information (elements, ...)
! first add elements to the program
  ierror = 0
  i = 1
  open (unit=76,file='netdef.in')
  read (76,*)
  read (76,*) xlostneu
  read (76,*)
  read (76,*)
  do while (ierror == 0)
   read (76,'(3x,i3,1x,i3,1x,1p,d23.15)',iostat=ierror)nbzel(i),nbael(i),abels(i)
   if (ierror /= 0) then
     close(76)
     exit
   endif
   if (verbose) then
     write(*,*) "THE VARS I NEED",nbzel(i),nbael(i),abels(i)
   endif
   i = i+1
  enddo

  nbelx=i-1


  zabelx=z
  do ii=1,nbelx
   zabelx=zabelx-abels(ii)
  enddo

  if (nbelx > mbelx) then
    write(*,*) 'nbelx= ',nbelx,' > mbelx= ',mbelx
    stop 'stop in netinit/netrates.f'
  endif

! then decide which element are followed in netnewr.f
! Use input name
!Approx21 flag to either use new rates or basic ones. Note that I enfore that all isotopes are included even in non approx21 for simplicity in rest of code and outputs.
  if ( iapprox21 == 1 ) then
    netinit_fileCNE = 'netinit_approx21.inCNE'
    netinit_fileCNEO = 'netinit_approx21.inCNEO'
    vit_fileCNE = 'vit_approx21.datCNE'
    vit_fileCNEO = 'vit_approx21.datCNEO'
  elseif (iapprox21 == 2 ) then
      netinit_fileCNE = 'netinit_approx21.inCNE'
      netinit_fileCNEO = 'netinit_approx21.inCNEO'
      vit_fileCNE = 'vit_approx21.datCNE'
      vit_fileCNEO = 'vit_approx21_noenerg.datCNEO'

  else
    netinit_fileCNE = 'netinit_approx21.inCNE'
    netinit_fileCNEO = 'netinit_approx21.inCNEO'
    vit_fileCNE = 'vit.datCNE'
    vit_fileCNEO = 'vit.datCNEO'
  endif


  if (phase < 4) then
    namenet=trim(input_dir)//'inputs/'//netinit_fileCNE
    namereac=trim(input_dir)//'inputs/'//vit_fileCNE   
  else
    namenet=trim(input_dir)//'inputs/'//netinit_fileCNEO
    namereac=trim(input_dir)//'inputs/'//vit_fileCNEO
  endif

  if (idebug > 0) then
    write(*,*) 'call readnetZA'
  endif
  call  readnetZA
  if (idebug > 0) then
    write(*,*) 'call readreac'
  endif
  call readreac

  return

end subroutine netinit
!=======================================================================
subroutine contribreac
!----------------------------------------------------------------------
  implicit none

  integer:: i,klo,khi,k,jlo,jhi

  real(kindreal):: t8,b,h,anegridhigh,deltagrid,vbetalow,vbetahigh
  real(kindreal), dimension(0:4):: f
!----------------------------------------------------------------------
  t8 = t9*10.d0

! factorials
  f(0) = 1.d0
  f(1) = 1.d0
  f(2) = 2.d0
  f(3) = 6.d0
  f(4) = 24.d0

! interpolate reaction rate
! locate the position of current T within the grid
! and compute step, a, b, c and d parameters for cubic spline interpolation
  klo = 1
  khi = kgrid
  do while (khi-klo > 1)
   k = (khi+klo)/2
   if (tgrid(k) > t8) then
     khi=k
   else
     klo=k
   endif
  enddo
  h=tgrid(khi)-tgrid(klo)
  b=(t8 - tgrid(klo))/h

! locate the position of current rho within the density grid (1.,3.,10.,30.)
  if (ane <= 3.d0) then
    anegridhigh=3.d0
    deltagrid=2.d0
    jlo=1
    jhi=2
  else if (ane <= 10.d0) then
    anegridhigh=10.d0
    deltagrid=7.d0
    jlo=2
    jhi=3
  else
    anegridhigh=30.d0
    deltagrid=20.d0
    jlo=3
    jhi=4
  endif

! perform linear interpolation of log(v)
  do i=1,ireac
   vrate(i)=rrate(i,shellnb)
   if (vrate(i) > 1.d-30) then
     if (flag(i) <= 0.d0) then
! reaction is not electron-density-dependent beta-decay
! test the reaction kind
       if (flag(i) == -14.d0) then
! electron capture on Be7
         if (itestx == 1 .and. verbose) then
           write(*,*)'special treatment for Be -> He4 in H-burning:'
           write(*,*) reaction(i)
         endif
         mata(elps(i,1),elps(i,1))=mata(elps(i,1),elps(i,1))-vrate(i)
! - instead of + -> H1
         mata(elps(i,3),elps(i,1))=mata(elps(i,3),elps(i,1))-vrate(i)*nsnb(i,3)
         mata(elps(i,4),elps(i,1))=mata(elps(i,4),elps(i,1))+vrate(i)*nsnb(i,4)
       else if (flag(i) == -13.d0) then
! electron capture
         mata(elps(i,1),elps(i,1))=mata(elps(i,1),elps(i,1))-vrate(i)
         mata(elps(i,3),elps(i,1))=mata(elps(i,3),elps(i,1))+vrate(i)*nsnb(i,3)
         mata(elps(i,4),elps(i,1))=mata(elps(i,4),elps(i,1))+vrate(i)*nsnb(i,4)
       else if (flag(i) == -11.d0) then
! photodisintegration or beta-decay
         mata(elps(i,1),elps(i,1))=mata(elps(i,1),elps(i,1))-vrate(i)
         if (elps(i,3) > 0.and.nsnb(i,3) > 0) then
           mata(elps(i,3),elps(i,1))=mata(elps(i,3),elps(i,1))+vrate(i)*nsnb(i,3)
         endif
         mata(elps(i,4),elps(i,1))=mata(elps(i,4),elps(i,1))+vrate(i)*nsnb(i,4)
       else if (flag(i) == -10.d0) then
! two-particle reaction   !if identical particles: factorials!
! linearisation used: Yi(n+1)*Yj(n+1)
!                     = Yi(n+1)*Yj(n)+Yi(n)*Yj(n+1)-Yi(n)*Yj(n)

! special treatment for deuterium in H-burning:
         if (reaction(i) == '2 H   1 ( 0 OOOOO, 0 OOOOO)  1 HE  3 ') then
           if (itestx == 1 .and. verbose) then
             write(*,*) 'special treatment for deuterium in H-burning:'
             write(*,*) reaction(i)
           endif

! 1st element variation: H  rate-->rate*3./2.
! dependence on 1st element only
           mata(elps(i,1),elps(i,1))=mata(elps(i,1),elps(i,1))-nsnb(i,1)*vrate(i)*2.d0*vabuny(elps(i,1)) *3.d0/2.d0
! right hand side term (RHS)
           mata(elps(i,1),nbel+1)=mata(elps(i,1),nbel+1)-nsnb(i,1)*vrate(i)*vabuny(elps(i,1))**2.d0 *3.d0/2.d0
! 4th element variation: He3
! dependence on 1st element only
           mata(elps(i,4),elps(i,1))=mata(elps(i,4),elps(i,1))+nsnb(i,4)*vrate(i)*2.d0*vabuny(elps(i,1))
! right hand side term (RHS)
           mata(elps(i,4),nbel+1)=mata(elps(i,4),nbel+1)+nsnb(i,4)*vrate(i)*vabuny(elps(i,1))**2.d0

! special treatment for Li, Be & B in H-burning:
         else if (reaction(i) == '1 HE  4 ( 1 HE  3, 0 OOOOO)  1 H   1 ') then
           if (itestx == 1 .and. verbose) then
             write(*,*) 'special treatment for Li, Be & B in H-burning:'
             write(*,*) reaction(i)
           endif
! 1st el. variation: He4 + instead of -
! dep. on 1st el.
           mata(elps(i,1),elps(i,1))=mata(elps(i,1),elps(i,1))+vrate(i)*vabuny(elps(i,2))
! dep. on 2nd el.
           mata(elps(i,1),elps(i,2))=mata(elps(i,1),elps(i,2))+vrate(i)*vabuny(elps(i,1))
! RHS
           mata(elps(i,1),nbel+1)=mata(elps(i,1),nbel+1)+vrate(i)*vabuny(elps(i,1))*vabuny(elps(i,2))
! 2nd el. variation: He3 same
! dep. on 1st el.
           mata(elps(i,2),elps(i,1))=mata(elps(i,2),elps(i,1))-vrate(i)*vabuny(elps(i,2))
! dep. on 2nd el.
           mata(elps(i,2),elps(i,2))=mata(elps(i,2),elps(i,2))-vrate(i)*vabuny(elps(i,1))
! RHS
           mata(elps(i,2),nbel+1)=mata(elps(i,2),nbel+1)-vrate(i)*vabuny(elps(i,1))*vabuny(elps(i,2))
! 4th element variation: proton - instead of +
! dep. on 1st el.
           mata(elps(i,4),elps(i,1))=mata(elps(i,4),elps(i,1))-nsnb(i,4)*vrate(i)*vabuny(elps(i,2))
! dep. on 2nd el.
           mata(elps(i,4),elps(i,2))=mata(elps(i,4),elps(i,2))-nsnb(i,4)*vrate(i)*vabuny(elps(i,1))
! RHS
           mata(elps(i,4),nbel+1)=mata(elps(i,4),nbel+1)-nsnb(i,4)*vrate(i)*vabuny(elps(i,1))*vabuny(elps(i,2))

! 2 a --> ... : nsnb(i,1)=2, nsnb(i,2)=0
         else if (elps(i,2) == 0.or.nsnb(i,2) == 0) then
! 1st element variation
! dependence on 1st element only
           mata(elps(i,1),elps(i,1))=mata(elps(i,1),elps(i,1))-nsnb(i,1)*vrate(i)*2.d0*vabuny(elps(i,1))
! right hand side term (RHS)
           mata(elps(i,1),nbel+1)=mata(elps(i,1),nbel+1)-nsnb(i,1)*vrate(i)*vabuny(elps(i,1))**2.d0
           if (elps(i,3) > 0.and.nsnb(i,3) > 0) then
! 3rd element variation
! dependence on 1st element only
             mata(elps(i,3),elps(i,1))=mata(elps(i,3),elps(i,1))+nsnb(i,3)*vrate(i)*2.d0*vabuny(elps(i,1))
! right hand side term (RHS)
             mata(elps(i,3),nbel+1)=mata(elps(i,3),nbel+1)+nsnb(i,3)*vrate(i)*vabuny(elps(i,1))**2.d0
           endif
! 4th element variation
! dependence on 1st element only
           mata(elps(i,4),elps(i,1))=mata(elps(i,4),elps(i,1))+nsnb(i,4)*vrate(i)*2.d0*vabuny(elps(i,1))
! right hand side term (RHS)
           mata(elps(i,4),nbel+1)=mata(elps(i,4),nbel+1)+nsnb(i,4)*vrate(i)*vabuny(elps(i,1))**2.d0

! a + b  --> ... : nsnb(i,1)=nsnb(i,2)=1
         else if (elps(i,2) > 0.and.nsnb(i,2) > 0) then
! 1st el. variation
! dep. on 1st el.
           mata(elps(i,1),elps(i,1))=mata(elps(i,1),elps(i,1))-vrate(i)*vabuny(elps(i,2))
! dep. on 2nd el.
           mata(elps(i,1),elps(i,2))=mata(elps(i,1),elps(i,2))-vrate(i)*vabuny(elps(i,1))
! RHS
           mata(elps(i,1),nbel+1)=mata(elps(i,1),nbel+1)-vrate(i)*vabuny(elps(i,1))*vabuny(elps(i,2))
! 2nd el. variation
! special treatment for F19 in H-burning:
           if (reaction(i) == '1 O  18 ( 1 PROT , 1 HE  4)  1 O  16 ') then
             if (itestx == 1) then
               write(*,*) 'special treatment for F19 in H-burning:'
               write(*,*) reaction(i)
             endif
! dep. on 1st el.
             mata(elps(i,2),elps(i,1))=mata(elps(i,2),elps(i,1))-2.d0* vrate(i)*vabuny(elps(i,2))
! dep. on 2nd el.
             mata(elps(i,2),elps(i,2))=mata(elps(i,2),elps(i,2))-2.d0* vrate(i)*vabuny(elps(i,1))
! RHS
             mata(elps(i,2),nbel+1)=mata(elps(i,2),nbel+1)-2.d0* vrate(i)*vabuny(elps(i,1))*vabuny(elps(i,2))
           else
! dep. on 1st el.
             mata(elps(i,2),elps(i,1))=mata(elps(i,2),elps(i,1))-vrate(i)*vabuny(elps(i,2))
! dep. on 2nd el.
             mata(elps(i,2),elps(i,2))=mata(elps(i,2),elps(i,2))-vrate(i)*vabuny(elps(i,1))
! RHS
             mata(elps(i,2),nbel+1)=mata(elps(i,2),nbel+1)-vrate(i)*vabuny(elps(i,1))*vabuny(elps(i,2))
           endif

           if (elps(i,3) > 0.and.nsnb(i,3) > 0) then
! 3rd element variation
! dep. on 1st el.
             mata(elps(i,3),elps(i,1))=mata(elps(i,3),elps(i,1))+nsnb(i,3)*vrate(i)*vabuny(elps(i,2))
! dep. on 2nd el.
             mata(elps(i,3),elps(i,2))=mata(elps(i,3),elps(i,2))+nsnb(i,3)*vrate(i)*vabuny(elps(i,1))
! RHS
             mata(elps(i,3),nbel+1)=mata(elps(i,3),nbel+1)+nsnb(i,3)*vrate(i)*vabuny(elps(i,1))*vabuny(elps(i,2))
           endif
! 4th element variation
! dep. on 1st el.
           mata(elps(i,4),elps(i,1))=mata(elps(i,4),elps(i,1))+nsnb(i,4)*vrate(i)*vabuny(elps(i,2))
! dep. on 2nd el.
           mata(elps(i,4),elps(i,2))=mata(elps(i,4),elps(i,2))+nsnb(i,4)*vrate(i)*vabuny(elps(i,1))
! RHS
           mata(elps(i,4),nbel+1)=mata(elps(i,4),nbel+1)+nsnb(i,4)*vrate(i)*vabuny(elps(i,1))*vabuny(elps(i,2))
         endif


       else if (flag(i) == -100.d0) then
! three-particle reactions   !if identical particles: factorials!
! 3 a --> ... : nsnb(i,1)=3, nsnb(i,2)=0
         if (elps(i,2) == 0.or.nsnb(i,2) == 0) then
! 1st element variation
! dependence on 1st element only
           mata(elps(i,1),elps(i,1))=mata(elps(i,1),elps(i,1))-3.d0*nsnb(i,1)*vrate(i)*vabuny(elps(i,1))**2.d0
! right hand side term (RHS)
           mata(elps(i,1),nbel+1)=mata(elps(i,1),nbel+1)-2.d0*nsnb(i,1)*vrate(i)*vabuny(elps(i,1))**3.d0
           if (elps(i,3) > 0.and.nsnb(i,3) > 0) then
! 3rd element variation
! dependence on 1st element only
             mata(elps(i,3),elps(i,1))=mata(elps(i,3),elps(i,1))+3.d0*nsnb(i,3)*vrate(i)*vabuny(elps(i,1))**2.d0
! right hand side term (RHS)
             mata(elps(i,3),nbel+1)=mata(elps(i,3),nbel+1)+2.d0*nsnb(i,3)*vrate(i)*vabuny(elps(i,1))**3.d0
           endif
! 4th element variation
! dependence on 1st element only
           mata(elps(i,4),elps(i,1))=mata(elps(i,4),elps(i,1))+3.d0*nsnb(i,4)*vrate(i)*vabuny(elps(i,1))**2.d0
! right hand side term (RHS)
           mata(elps(i,4),nbel+1)=mata(elps(i,4),nbel+1)+2.d0*nsnb(i,4)*vrate(i)*vabuny(elps(i,1))**3.d0
         endif
       endif
     else
! beta-decay rate has to be interpolated in density
! linear interpolation in temperature for two density grid points
!                (note: log of beta decay rate is handled)
       vbetalow =log10(vgrid(klo,i,jlo))+ b * (log10(vgrid(khi,i,jlo))-log10(vgrid(klo,i,jlo)))
       vbetahigh=log10(vgrid(klo,i,jhi))+ b * (log10(vgrid(khi,i,jhi))-log10(vgrid(klo,i,jhi)))
       vbetalow = 10.d0**vbetalow
       vbetahigh= 10.d0**vbetahigh

! linear interpolation in density
!  extrapolation often leads to negative rates
! check if so, and then put it to zero
!       (supplementary table values would be needed for those cases)
       if (vrate(i) < 0.d0) then
         vrate(i) = 0.d0
       endif
     endif
     if (itestx == 1) then
       if (ireac == 1 .and. verbose) then
         write(*,*)'T8: ',t8,' rho: ',rho,' ane: ',ane
         write(*,*) reaction(i),', vitr= ',rrate(i,shellnb),', vitv= ',vrate(i) , ' & flag= ',flag(i)
       endif
     endif
   endif
  enddo

  return

end subroutine contribreac
!=======================================================================
subroutine readnetZA
!----------------------------------------------------------------------
  implicit none

  integer:: i,ii,j,iel,ierror
!----------------------------------------------------------------------
  nbel=0
  iel=1
  do i=1,maxel
    nbz(i)=-222
    nba(i)=-222
    nbn(i)=-222
  enddo

! read network
  ierror = 0
  open(unit=77,file=namenet)
  read(77,*)
  do while(ierror == 0)
   read(77,'(2x,2i4)',iostat=ierror) nbz(iel),nba(iel)
   if (ierror /= 0) then
     close(77)
     exit
   endif
   nbn(iel)=nba(iel)-nbz(iel)
   iel=iel+1
  enddo
  nbel=iel-1
  if (nbel > maxel) then
   write(*,'(2(a,i5),/,a)') 'Nbr of nuclei in network= ',nbel,' > maxel =', maxel,' maxel has to be increased ---> STOP'
   stop
  endif

! calculation of ineut, iprot & ialpha
  ineut=-2
  iprot=-2
  ialpha=-2
  do i=1,nbel
    if (nbz(i) == 0 .and. nba(i) == 1) ineut=i
    if (nbz(i) == 1 .and. nba(i) == 1) iprot=i
    if (nbz(i) == 2 .and. nba(i) == 4) ialpha=i
  enddo
  if (ineut < 1 .or. iprot < 1 .or. ialpha < 1 .and. verbose) then
    write(*,'(3(a,i3))') 'ineut=',ineut,' iprot=',iprot,' ialpha=',ialpha
    write(*,*) '  ---> STOP in readnet'
  endif

! calculation of nzmax, nnmax, nzmin & nnmin
  nbzmax=0
  nbnmax=0
  nbzmin=100
  nbnmin=100
  do i=0,maxz
    do j=0,maxz
      posel(i,j)=-2
    enddo
  enddo
  posel(0,0)=0
  do i=1,nbel
    if (nbz(i) > nbzmax) nbzmax=nbz(i)
    if (nbz(i) < nbzmin) nbzmin=nbz(i)
    if (nbn(i) > nbnmax) nbnmax=nbn(i)
    if (nbn(i) < nbnmin) nbnmin=nbn(i)
    posel(nbz(i),nbn(i))=i
  enddo

  checkel: do i=1,nbel
    if (nbz(i) == 0  .and. nba(i) == 1)  cycle  !Cycle for neutrons
   if (nbz(i) == 1  .and. nba(i) == 1)  cycle
   if (nbz(i) == 2  .and. nba(i) == 3)  cycle
   if (nbz(i) == 2  .and. nba(i) == 4)  cycle
   if (nbz(i) == 6  .and. nba(i) == 12) cycle
   if (nbz(i) == 6  .and. nba(i) == 13) cycle
   if (nbz(i) == 7  .and. nba(i) == 14) cycle
   if (nbz(i) == 7  .and. nba(i) == 15) cycle
   if (nbz(i) == 8  .and. nba(i) == 16) cycle
   if (nbz(i) == 8  .and. nba(i) == 17) cycle
   if (nbz(i) == 8  .and. nba(i) == 18) cycle
   if (nbz(i) == 10 .and. nba(i) == 20) cycle
   if (nbz(i) == 10 .and. nba(i) == 22) cycle
   if (nbz(i) == 12 .and. nba(i) == 24) cycle
   if (nbz(i) == 12 .and. nba(i) == 25) cycle
   if (nbz(i) == 12 .and. nba(i) == 26) cycle

   do ii=1,nbelx
    if (nbz(i)==nbzel(ii) .and. nba(i)==nbael(ii)) cycle checkel
   enddo
   print*,'element ',i,nbz(i),nba(i),' not followed in the prog.'
   stop
  enddo checkel

  if (nbzmax > maxz .and. verbose) then
    write(*,'(2i5)') 'Attention: nbzmax= ',nbzmax,' > maxz= ',maxz
  endif

  if (itestx == 2 .and. verbose) then
    write(*,*)"table of elements:"
    do i=nbzmax,nbzmin,-1
     write(*,'(a2,1p,99(i3))') sym(i),(posel(i,j),j=nbnmin,nbnmax)
    enddo
    do i=1,nbel
     write(*,'(a2,1p,i3,a,i3)') sym(nbz(i)),nba(i),', nn= ',nbn(i)
    enddo
  endif

  return

end subroutine readnetZA
!=======================================================================
subroutine readreac
!----------------------------------------------------------------------
! This routine provides the reaction rates v(i) (i=1,nreac)
! for temperature T and density RHO
! from a linear interpolation of the log of the grid-point rates read
! from the file 'vit.dat' containing the output from Netgen.

! Note that the factorials accounting for identical particles have
! already been included in v(i).
! To obtain the evolution dYj/dt of species j, simply multiply v(i) by
! the molar fractions of the reacting species j1 + j2: Yj1^n1 Yj2^n2,
! where n1, n2 are the stoechiometric factors (stored in the vectors
! with the same name).

! The last lines of subroutine vit print all variables to standard output,
! to check that it worked properly

! The first call to vit should be done with
! iread = 0: to read the vit.dat file, and then store them in array vgrid

! Subsequent calls may be done with
! iread = 1: to compute rates at other temperatures, after vgrid
!            has been initialised by the first call with iread = 0

! The parameters
!      NGRID = maximum number of temperature grid points
!      NRE   = maximum number of reactions
! may be lowered to better match your network.

! The actual number of grid points and reaction rates are computed
! automatically, and stored in variables IREAC and KGRID

! PME = mean molecular weight per electron
! T = temperature in K
! RHO = density in g/cm3

! The rates on grid points are stored in the array
!     vgrid(ngrid,nre,4)
! where the last index refers to the density grid for the beta-decays
! dependent upon electron number density:
! vgrid(ngrid,nre,1) corresponds to the rate at n_e = 1 x 10^26 cm-3
! vgrid(ngrid,nre,2) corresponds to the rate at n_e = 3 x 10^26 cm-3
! vgrid(ngrid,nre,3) corresponds to the rate at n_e = 10 x 10^26 cm-3
! vgrid(ngrid,nre,4) corresponds to the rate at n_e = 30 x 10^26 cm-3

! For the other reactions
! (with flags different from 1., 3., 10. or 30.), only
! vgrid(ngrid,nre,1) is relevant

! Coding of reactions:
!   n1 to n4: stoechiometric factors
!   z1 to z4: chemical symbol
!   a1 to a4: atomic mass

! The following table allows to convert proton number into
! chemical symbol, if required:
!----------------------------------------------------------------------
  implicit none

  integer:: iskip,i,j,ll,l,leng,k,i1,m,zel,icancel,mgrid
  integer,dimension(10,4):: nn

  real(kindreal),dimension(10):: qqrad,qqnu
  real(kindreal),dimension(0:4):: f
  real(kindreal),dimension(ngrid,10):: vdum

  character(2) elsym
  character(121),dimension(7):: longline
  character(8),dimension(10):: aflag
  character(6),dimension(10,4):: zz
!----------------------------------------------------------------------
  ireac=0
  kgrid=0
  do ll=1,nre
   do m=1,4
    do k=1,ngrid
     vgrid(k,ll,m)=-1.d+77
    enddo
    anb(ll,m)=-222.d0
    znb(ll,m)=-222.d0
    nsnb(ll,m)=-222
    elps(ll,m)=-222
    zed(ll,m)='      '
   enddo
   vrate(ll)=-1.d22
   qrad(ll)=-1.d22
   qnu(ll)=-1.d21
   flag(ll)=12345.d0
   reaction(ll)='a + b --> c + d'
  enddo

! factorials
  f(0) = 1.d0
  f(1) = 1.d0
  f(2) = 2.d0
  f(3) = 6.d0
  f(4) = 24.d0

! T8 = T * 1.E-08
! ane is the electron density in units of 1E26 cm-3
!     = rho * 6.02E-03 / PME
! PME is the mean molecular weight per electron
!     =  [ SUM (XZ/A) ]**(-1)

! if iread = 0, start by reading the data file, then compute the rates
!      if (iread == 0) then
  open(unit=78,file=namereac)

! check the number of header lines (iskip)
! check the number of grid points (kgrid)
! read the temperature grid Tgrid
  mgrid = ngrid + 30

  readloop: do i=1,mgrid
   read(78,'(a80)') longline(1)(1:80)
   if (longline(1)(1:1) /= '#') then
     iskip = i
     backspace(78)
     do j=1,ngrid
      read(78,'(a80)') longline(1)(1:80)
      if (longline(1)(1:1) /= '#') then
         read(longline(1),'(3x,f8.4)') Tgrid(j)
      else
         kgrid=j-1
         rewind(78)
         exit readloop
      endif
     enddo
   endif
  enddo readloop

! read the data table
  ireac = 0
  icancel = 0

  do ll=1,10000

! read the header
   read(78,'(4(a121,//),3(a121,/))',end=9999) (longline(j),j=1,7)

! l = number of data records on line
   leng = index(longline(5),'        ')
   if (leng == 0) leng = 121
   l = (leng - 11)/11
   read(longline(7),'(14x,10(a8,3x))') (aflag(j),j=1,l)
   do m=1,4
    read(longline(m),'(14x,10(i1,1x,a6,3x))')(nn(j,m),zz(j,m),j=1,l)
   enddo
   read(longline(5),'(14x,10(f8.3,3x))',end=11) (QQrad(j),j=1,l)
11 read(longline(6),'(14x,10(f8.3,3x))',end=12) (QQnu(j), j=1,l)
12 continue

! read the data
   read(78,*)
   read(78,*)

   do i=1,kgrid
    read(78,'(12x,9(e10.4,1x),e10.4)',end=13) (vdum(i,j),j=1,l)
   enddo
   read(78,*)
13 continue

! transfer to vgrid
   do j=1,l
    if (aflag(j)(8:8) /= '3' .and. aflag(j)(7:8) /= '10' .and. aflag(j)(7:8) /= '30') then
      if (icancel == 0) then
        ireac = ireac + 1
      else if (icancel == 1) then
        do m=1,4
         do k=1,ngrid
          vgrid(k,ireac,m)=-1.d+77
         enddo
         anb(ireac,m)=-222.d0
         znb(ireac,m)=-222.d0
         nsnb(ireac,m)=-222
         elps(ireac,m)=-222
         zed(ireac,m)='      '
        enddo
        vrate(ireac)=-1.d22
        qrad(ireac)=-1.d22
        qnu(ireac)=-1.d21
        flag(ireac)=12345.d0
        reaction(ireac)='a + b --> c + d'
        icancel = 0
      endif

      do k=1,kgrid
       vgrid(k,ireac,1) = vdum(k,j)
      enddo

      do m=1,4
       nsnb(ireac,m) = nn(j,m)
       zed(ireac,m) = zz(j,m)
      enddo
      Qrad(ireac)=QQrad(j)
      Qnu(ireac) =QQnu(j)

      if (aflag(j)(5:8) == '----') then
        flag(ireac) = -200.d0
      else if (aflag(j)(6:8) == '---') then
        flag(ireac) = -100.d0
      else if (aflag(j)(7:8) == '--') then
        flag(ireac) = -10.d0
      else if (aflag(j)(6:8) == '+++') then
        flag(ireac) = -14.d0
      else if (aflag(j)(7:8) == '++') then
        flag(ireac) = -13.d0
      else if (aflag(j)(8:8) == '+') then
        flag(ireac) = -11.d0
      endif

      do m=1,4
       if (zz(j,m)(1:4) == 'NEUT') then
         anb(ireac,m) = 1.d0
         znb(ireac,m) = 0.d0
       else if (zz(j,m)(1:4) == 'PROT') then
         anb(ireac,m) = 1.d0
         znb(ireac,m) = 1.d0
       else if (zz(j,m)(1:5) == 'OOOOO') then
         anb(ireac,m) = 0.d0
         znb(ireac,m) = 0.d0
       else if (zz(j,m)(1:4) == 'DEUT') then
         anb(ireac,m) = 2.d0
         znb(ireac,m) = 1.d0
       else if (zz(j,m)(1:4) == 'TRIT') then
         anb(ireac,m) = 3.d0
         znb(ireac,m) = 1.d0
       else if (zz(j,m)(1:5) == 'HE  4') then
         anb(ireac,m) = 4.d0
         znb(ireac,m) = 2.d0
       else if (zz(j,m)(1:5) == 'HE  3') then
         anb(ireac,m) = 3.d0
         znb(ireac,m) = 2.d0
       else if (m /= 3) then
         read(zz(j,m)(3:5),'(i3)') i1
         anb(ireac,m) = real(i1)
         read(zz(j,m)(1:2),'(a2)') elsym
         do zel=0,maxz
          if (sym(zel) == elsym) znb(ireac,m) = zel*1.d0
         enddo
       endif
       elps(ireac,m) = posel(int(znb(ireac,m)),int(anb(ireac,m)-znb(ireac,m)))
       if (elps(ireac,m) < 0) then
         icancel=1
       endif
      enddo

      write(reaction(ireac),'(i1,1x,a6,"(",i1,1x,a5,",",1x,i1,1x,a5,")",2x,i1,1x,a6)') nn(j,1),zz(j,1)(1:6),nn(j,2),&
                                                                 zz(j,2)(1:5),nn(j,3),zz(j,3)(1:5),nn(j,4),zz(j,4)(1:6)

    else if (aflag(j)(8:8) == '3') then
      flag(ireac) = 1.d0
      do k=1,kgrid
       vgrid(k,ireac,2) = vdum(k,j)
      enddo

    else if (aflag(j)(7:8) == '10') then
      do k=1,kgrid
       vgrid(k,ireac,3) = vdum(k,j)
      enddo

    else if (aflag(j)(7:8) == '30') then
      do k=1,kgrid
       vgrid(k,ireac,4) = vdum(k,j)
      enddo

    endif
   enddo

  enddo

9999 continue

  close(78)
  if (icancel == 1) then
    do m=1,4
     do k=1,ngrid
      vgrid(k,ireac,m)=-1.d22
     enddo
     anb(ireac,m)=-222.d0
     znb(ireac,m)=-222.d0
     nsnb(ireac,m)=-222
     elps(ireac,m)=-222
     zed(ireac,m)='      '
    enddo
    vrate(ireac)=-1.d22
    qrad(ireac)=-1.d22
    qnu(ireac)=-1.d21
    flag(ireac)=12345.d0
    reaction(ireac)='a + b --> c + d'
    ireac=ireac-1
    icancel = 0
  endif

  if (itestx == 1) then
    do i = 1,ireac+1
     write(6,'(1x,a37,1x,a,f5.0)') reaction(i),'flag= ',flag(i)
     write(6,'(1x,"i1 = ",i1,1x,"a1 = ",f4.0,1x,"z1 = ",a6,f4.0,i3,/, &
             & 1x,"i2 = ",i1,1x,"a2 = ",f4.0,1x,"z2 = ",a6,f4.0,i3,/, &
             & 1x,"i3 = ",i1,1x,"a3 = ",f4.0,1x,"z3 = ",a6,f4.0,i3,/, &
             & 1x,"i4 = ",i1,1x,"a4 = ",f4.0,1x,"z4 = ",a6,f4.0,i3)') &
                   nsnb(i,1),anb(i,1),zed(i,1),znb(i,1),elps(i,1), &
                   nsnb(i,2),anb(i,2),zed(i,2),znb(i,2),elps(i,2), &
                   nsnb(i,3),anb(i,3),zed(i,3),znb(i,3),elps(i,3), &
                   nsnb(i,4),anb(i,4),zed(i,4),znb(i,4),elps(i,4)
     write(6,'(2(1x,a,1x,f7.3))') 'Qrad (MeV) =',Qrad(i),'Qnu (MeV) =',Qnu(i)
     do j=1,10
      write(6,'(10(e12.6,1x))')(log10(vgrid(k,i,1)),k=(j-1)*6+1,j*6)
     enddo
    enddo
  endif

  if (kgrid > ngrid .and. verbose) then
    write(*,'(2i5)') 'Attention: kgrid= ',kgrid,' > ngrid= ',ngrid
  endif
  if (ireac > nre .and. verbose) then
    write(*,'(2i5)') 'Attention: ireac= ',ireac,' > nre= ',nre
  endif

  return

end subroutine readreac
!=======================================================================
subroutine inversemat(nbel,mata,abuny,maxel2)
!----------------------------------------------------------------------
  use inputparam,only: idebug
  use caramodele,only: nwmd
  use SmallFunc,only: girl

  implicit none

  integer,intent(in):: nbel,maxel2
  real(kindreal),dimension(maxel,maxel+1),intent(in):: mata
  real(kindreal),dimension(maxel),intent(out):: abuny

  integer:: i,j,flag_girl,iSE
  real(kindreal),dimension(maxel):: b
  real(kindreal),dimension(maxel*(maxel+1)):: aa
!----------------------------------------------------------------------
  if (maxel /= maxel2) then
    print*,'maxel= ',maxel,'# maxel2= ',maxel2
    stop
  endif
  do i=1,nbel
   do j=1,nbel+1
    aa(i+(j-1)*(nbel))=mata(i,j)
   enddo
  enddo

  call girl(aa,b,nbel,1,flag_girl)
  if (flag_girl /= 0) then
    if (idebug>0) then
      write(*,*) 'energy, inversemat - matrix aa,flag:',flag_girl
      do iSE=1,maxel*(maxel+1)
       write(*,'("aa(",i1,") :",d22.12)') aa(iSE)
      enddo
    endif
    rewind(222)
    write(222,*) nwmd,':girl crash in inversemat with matrix aa'
    stop
  endif


  do i=1,nbel
   abuny(i)=b(i)
  enddo

  return

end subroutine inversemat
!======================================================================
subroutine nucal
!-----------------------------------------------------------------
!   CALCUL DU FLUX DES NEUTRINOS SOLAIRES

! Auteur : Y. Lebreton
! Modifications : C.Charbonnel
! Derniere version : 18 novembre 1992
!-----------------------------------------------------------------
  use const,only: Msol,cst_avo,pi,uastr
  use abundmod ,only: b11,x,tauxbe7el,ybe7,tauxbe7pg,b112,xc12,b114,xn14,b116,xo16
  use strucmod,only: m,q,t,rho

  implicit none

  integer:: ncouch,kf,j
! *** LES SIGMA 37CL ET 71GA SONT EN 10E-46 CM2  ***
  real(kindreal),dimension(7), parameter:: sigcl=(/11.8d0,215.d0,73.2d0,2.43d4,61.8d0,116.d0,117.d0/),&
                                           sigga=(/0.0d0,16.d0,2.38d0,1.06d+4,1.66d0,6.6d0,6.67d0/)

  real(kindreal):: a,b,c,d,coef,t6ln,t63,t623,t6,coef1,coef2,conv,tot
  real(kindreal), dimension(ldi):: qq,aflu
  real(kindreal), dimension(7):: flux,flix

  ncouch=M

  a=q(1)
  b=exp(a)
  qq(1)=a
  do j=1,ncouch-1
   c=q(j+1)
   d=exp(c)
   qq(j+1)=(b-d)/(c-a)
   a=c
   b=d
  enddo

! *** FLUX 1H(P,E+ NU) 2D   ***
  kf=1
  do j=1,ncouch
   aflu(j)=b11(j)*x(j)*x(j)*msol/2.d0
  enddo

  call nuint(qq,aflu,ncouch,flux(kf))

! *** FLUX 1H(PE+,NUE) 2D   ***
  kf=2
  coef=-6.d0*log(10.d0)
  do j=1,ncouch
   t6ln=t(j)+coef
   t63=exp(t6ln/3.d0)
   t623=t63**2.d0
   t6=t623*t63
   coef1=(1.d0-.0729d0*t63+.0982d0*t623)*exp(rho(j))*(1.d0+x(j))
   coef2=1.d0+.94d-3*t6+.0109d0*t623+.0123d0*t63
   aflu(j)=aflu(j)*coef1*5.51d-5/(coef2*sqrt(t6))
  enddo

  call nuint(qq,aflu,ncouch,flux(kf))

! *** FLUX 7BE(E-  NU) 7LI  ***
  kf=3
  do j=1,ncouch
   aflu(j)=tauxbe7el(j)*7.d0*ybe7(j)*(1.d0+x(j))*msol/14.d0
  enddo

  call nuint(qq,aflu,ncouch,flux(kf))

! *** FLUX 7BE(P,GAM) 8B (E+ NU) 8BE*   ***
  kf=4
  do j=1,ncouch
   aflu(j)=tauxbe7pg(j)*x(j)*7.d0*ybe7(j)*msol/7.d0
  enddo

  call nuint(qq,aflu,ncouch,flux(kf))

! *** FLUX 12C(P,GAM)N13(E+ NU) 13C   ***
  kf=5
  do j=1,ncouch
   aflu(j)=b112(j)*xc12(j)*x(j)*msol/12.d0
  enddo

  call nuint(qq,aflu,ncouch,flux(kf))

! *** FLUX 14N(P,GAM)O15(E+ NU) 15N   ***
  kf=6
  do j=1,ncouch
   aflu(j)=b114(j)*xn14(j)*x(j)*msol/14.d0
  enddo

  call nuint(qq,aflu,ncouch,flux(kf))

! *** FLUX 16O(P,GAM)F17(E+ NU) 17O   ***
  kf=7
  do j=1,ncouch
   aflu(j)=b116(j)*xo16(j)*x(j)*msol/16.d0
  enddo

  call nuint(qq,aflu,ncouch,flux(kf))

! ***  FLUX SUR TERRE ***
  conv=((cst_avo/1.d23)/(4.d0*pi*uastr*uastr/1.d26))*1.d-13
  write(3,'(1x/10x,a/25x,a,8x,a,7x,a,7x,a,8x,a,7x,a,7x,a,7x,a/)') 'FLUX DES NEUTRINOS SOLAIRES','PP','PEP','7BE',&
                                                                  '8B','13N','15O','17F','TOTAL SNU'
  do kf=1,7
   flux(kf)=conv*flux(kf)
  enddo

  write(3,'(5x,a,5x,7(1pe8.1,2x)/8x,a/)') 'PHI CM-2 S-1',(flux(kf), kf=1,7),'10**10'
  tot=0.d0

  do kf=1,7
   flix(kf)=flux(kf)
   flux(kf)=sigcl(kf)*flux(kf)
   tot=tot+flux(kf)
  enddo
  write(3,'(15x,a,7(2x,f8.4),6x,f8.4//)') '37CL',(flux(kf), kf=1,7),tot

  tot=0.d0
  do kf=1,7
   flux(kf)=sigga(kf)*flix(kf)
   tot=tot+flux(kf)
  enddo
  write(3,'(15x,a,7(2x,f8.4),6x,f8.4//)') '71GA',(flux(kf),kf=1,7),tot

  return

end subroutine nucal
!======================================================================
subroutine nuint(qq,aflu,ncouch,flux)
!-----------------------------------------------------------------
!  CALCUL DE L'INTEGRALE (METHODE DES TRAPEZES)
!   POUR LE FLUX DES NEUTRINOS SOLAIRES
!-----------------------------------------------------------------
  implicit none

  integer,intent(in):: ncouch
  real(kindreal),intent(in), dimension(ldi):: qq,aflu
  real(kindreal),intent(out):: flux

  integer:: j
  real(kindreal):: a,b=0.d0
!-----------------------------------------------------------------
  a=aflu(1)
  flux=-a*qq(1)
  do j=1,ncouch-1
   b=aflu(j+1)
   flux=flux+(b-a)*qq(j+1)
   a=b
  enddo
  flux=flux+b

  return

end subroutine nuint
!======================================================================
subroutine enint
! sous-routine de calcule de l energie potentielle gravitationnelle, de
! l'energie rotationnelle, de l'energie cinetique du gaz et de l'energie
! de rayonnement
  use const,only: Msol,um,cst_G,cst_a
  use abundmod,only: epote,ekrote,ekine,erade
  use strucmod,only: q,r,vmyhelio
  use rotmod,only: bmomin,vomegi

  implicit none

  integer::imb,nmb,jm
  real(kindreal):: xtmu,xtmum1
  real(kindreal), dimension(ldi):: xq,xmasr,dmasr,ray,xtemp,xxmu
  real(kindreal), dimension(ldi):: epot,ekrot,ekin,erad,dpot,dkrot,dekin,drad

  do imb=1,m
   xq(imb)=1.d0-exp(q(imb))
   xmasr(imb)=xq(imb)*gms*Msol ! mass of sphere in grams
   ray(imb)=10.d0**(r(imb)/um) ! radius of sphere in cms
   xtemp(imb)=exp(t(imb))      ! temperature
   if (vmyhelio(imb) /= 0.d0) then
     xxmu(imb)=vmyhelio(imb)
   else
     if (imb /= 1) then
       xxmu(imb) = vmyhelio(imb-1)
     else
       xxmu(imb) = vmyhelio(imb+1)
     endif
   endif
  enddo

  dmasr(1)=gms*(1.d0-xq(1))*Msol ! mass of shell in grams
  do imb=2,m
   dmasr(imb)=(xq(imb-1)-xq(imb))*gms*Msol
  enddo

! ENERGIE POTENTIELLE
  epote=-cst_G*xmasr(m-1)*xmasr(m-1)/ray(m-1)*3.d0/5.d0
  dpot(m)=epote
  epot(m)=epote
  do imb=2,m-1
   nmb=m-imb+1
   dpot(nmb)=-cst_G*dmasr(nmb)/2.d0*(xmasr(nmb)/ray(nmb)+xmasr(nmb-1)/ray(nmb-1))
   epote=epote+dpot(nmb)
   epot(nmb)=epote
  enddo

! ENERGIE ROTATIONNELLE
  ekrote=0.d0
  do imb=1,m-1
   nmb=m-imb+1
   dkrot(nmb)=0.5d0*bmomin(nmb)*vomegi(nmb)*vomegi(nmb)
   ekrote=ekrote+dkrot(nmb)
   ekrot(nmb)=ekrote
  enddo

! ENERGIE CINETIQUE DU GAZ
  ekine=0.d0
  do jm=1,m-1
   nmb=m-jm+1
! numerical factor = 3/2 k/m_h)
   if (xxmu(nmb) /= 0.d0) then
     xtmu = xtemp(nmb)/xxmu(nmb)
   else
     xtmu = 0.d0
   endif
   if (xxmu(nmb-1) /= 0.d0) then
     xtmum1 = xtemp(nmb-1)/xxmu(nmb-1)
   else
     xtmum1 = 0.d0
   endif
   dekin(nmb)=(1.50d0*cst_k/cst_mh)*dmasr(nmb)/2.d0*(xtmu+xtmum1)
   ekine=ekine+dekin(nmb)
   ekin(nmb)=ekine
  enddo

! ENERGIE DU RAYONNEMENT
  erade=0.d0
  do jm=1,m-1
   nmb=m-jm+1
   drad(nmb)=(4.d0*pi*cst_a/3.d0)*(ray(nmb-1)**3.d0-ray(nmb)**3.d0)*((xtemp(nmb)+xtemp(nmb-1))/2.d0)**4.d0
   erade=erade+drad(nmb)
   erad(nmb)=erade
  enddo

  return

end subroutine enint
!======================================================================
end module energy
