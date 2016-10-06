module caramodele

  use evol, only: kindreal,ldi

  implicit none

  logical, save:: firstmods,PrintError
  integer, save:: nwmd,inum,iprezams

! CARACTERISTIQUES GLOBALES
  real(kindreal), save:: xmini,gms,gls,glm,teff,glsv,teffv
  real(kindreal), save:: hh1,hh6,rayequat
  real(kindreal), save:: rhoc,tc

! PERTE DE MASSE
  integer, save:: iwr
  real(kindreal), save:: dm_lost,xmdot,eddesc,xLstarbefHen,xltotbeg

! BINAIRES
  real(kindreal), save:: ab

end module caramodele
