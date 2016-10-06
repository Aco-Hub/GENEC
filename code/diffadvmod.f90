module diffadvmod

  use evol, only: ldi,kindreal

  implicit none

! CONDITIONS AU BORD DE L'ADVECTION
  integer, save:: mtu,npasr,ibext,ibint,inoadv
  real(kindreal), save:: xbint,xbext,xcint,xcext,xomint,xomext

! CIRCULATION MERIDIENNE
  real(kindreal), save:: xueff
  real(kindreal), dimension(ldi), save:: ucicoe,vcicoe,ursmooth

! COEFFICIENTS DE DIFFUSION (surface=1, centre=m)
  real(kindreal), dimension(ldi), save:: D_shear,D_conv,D_eff,Richardson,xnabyy,K_ther
! COEFFICIENTS DE DIFFUSION (centre=1, surface=m)
  real(kindreal), dimension(ldi), save:: dchim,D_h

! TEMPS DE DIFFUSION
  integer, save:: jdiff
  real(kindreal), save:: tdiff

end module diffadvmod
