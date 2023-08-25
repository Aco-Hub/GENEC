module rotmod

  use evol, only: kindreal,ldi,npondcouche

  implicit none

! TAUX DE ROTATION
  real(kindreal), save:: rapcri,rapom2,rapvco,rapvc2,xobla,xoblaj
  real(kindreal), save:: vpsi,fMdot_rot,rrro,ygmoye
  logical,save:: ivcalc

! OMEGA
  real(kindreal), save:: do1dr
  real(kindreal), dimension(ldi), save:: omegi,vomegi,vvomeg,omegp,omegd,dlodlr,wom
  real(kindreal), dimension(npondcouche), save:: CorrOmega

! MOMENT CINETIQUE
  real(kindreal), save:: xlexcs,xldoex,dlelex,dlelexprev,bdotis,bmomit,btot,&
                         btotatm,Flux_remaining,BTotal_EndAdvect, &
                         BTotal_StartModel,dlelexsave,timestep_control
  real(kindreal), save:: thext1,thext2,condbe,suminenv,vsuminenv,vvsuminenv
  real(kindreal), dimension(ldi), save:: vxmoci,xinert,bmomin,btotq

! ADVECTION
  real(kindreal), dimension(ldi), save:: xmeg,xlo,theta,ht,xmu,gtilde,ebem,aux,ur,ur1,deladv,vcirc

end module rotmod
