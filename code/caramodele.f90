module caramodele

  use evol,only: kindreal,ldi

  implicit none

  logical,save:: firstmods,PrintError
  integer,save:: nwmd,inum,modell,nwseqini

! GLOBAL CHARACTERISTICS
  real(kindreal),save:: xmini,gms,gls,glm,teff,glsv,teffv,zams_radius,radius,xini
  real(kindreal),save:: hh1,hh6,rayequat
  real(kindreal),save:: rhoc,tc
  real(kindreal),save:: xteffprev,xtefflast,xlprev,xllast,xrhoprev,xrholast,&
                        xcprev,xclast,xtcprev,xtclast

! MASS LOSS
  integer, save:: iwr
  real(kindreal),save:: dm_lost,xmdot,eddesc,xLstarbefHen,xltotbeg,Mdot_NotCorrected

! BINARIES
  real(kindreal),save:: ab

  public

end module caramodele
