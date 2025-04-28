! module contenant les abondances, taux de reaction, opacites
module strucmod

  use evol, only: kindreal,ldi

  implicit none

! COUCHES
  integer, save:: m,j,j1,NPcoucheEff
  real(kindreal), save:: shl,dmsh,vmx

! STRUCTURE
  real(kindreal), save:: beta,beta1,vmy1,vmyo,vmye,rad,rad1,radm,adi,adi1,adim,adip,adip1,adit,adit1,ccrad1,&
    xnabj,xnabj1,xnabm,zrad,zrad1,zradm,xbruj1,k2_AMC

  real(kindreal), dimension(ldi), save:: bet,adgrad,Nabla_rad,Nabla_ad,delt,Nabla_mu,H_P,gravi,q,e,p,t,r,s,vp,vt,vr,vs,&
    rho,zensi,rprov,pb,tb,rb,sb,xb,amu,amub,vmyhelio,epsit,gradp,gradt,xmufit,xomegafit

! OPACITES:
  real(kindreal), save:: cap,cap1,capp,capp1,capt,capt1
  real(kindreal), dimension(ldi), save:: opac,opact

! VARIABLES DE L'ENVELOPPE
  integer, save:: izsa,it2,ih,ihv,nr,konv,konvv,kzone,id1,id2,neudr
  integer, save:: ProfIon,nenvel
  integer, dimension(3), save:: kk

  real(kindreal), save:: Eddmax
  real(kindreal), save:: chem,ychem,uvlp,g,vmol,h,vlttau,vlka,vkap,vkat,beta_env,cp, &
    vna,vmion,vmionp,vmiont,vlro,vlte,vll,vnr,vlt,vlr,vlm,u,xpsi_atm,rap1_atm,xft_atm, &
    rap2_atm,xgmoym_atm,patmos,tatmos,taulim,pradlim,vlmm,vlrr,vlll,uvlpm,vltm,vlrm,tau_conv
  real(kindreal), save:: rlp,rlt,rlc,rrp,rrt,rrc,rtp,rtt,rtc,dk
  real(kindreal), save:: fitmIon
  real(kindreal), dimension(3), save:: vny,y4,x_env,y5int
  real(kindreal), dimension(3), save:: drp,drt,drr,drl,drte
  real(kindreal), dimension(5,3), save:: f
  real(kindreal), dimension(1000,50), save:: envel

! VARIABLES POUR LES ROUTINES THERMO
  real(kindreal), save:: vlro_thermo,rhp_thermo,rht_thermo,xmol_thermo,psi_thermo
  real(kindreal), save:: vlmix,dwdo2,vlhp,vlgr

end module strucmod
