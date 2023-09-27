! module contenant les abondances, taux de reaction, opacites
module equadiffmod

  use evol, only: kindreal

  implicit none

  integer, save::  iter,iprc,jterma,iadnok
  real(kindreal), save:: gkor,gkorv
  real(kindreal), save:: ccg1,ccg2,ccg3,ccz2,ccz3

! ELEMENTS D'EQUATIONS DANS GI ET ZI
  real(kindreal), save:: g1,g1s,g1p,g1t,g1s1,g1p1,g1t1,g2,g2r,g2p,g2r1,g2p1, &
    g3,g3r,g3p,g3t,g3r1,g3p1,g3t1,g4,g4r,g4s,g4p,g4t,g4r1,g4s1,g4p1,g4t1
  real(kindreal), save:: z1,z1p,z1t,z1p1,z1t1,z2,z2p,z2p1,z2t1,z3,z3p1,z3t1, &
    z4,z4s,z4p,z4t,z4p1,z4t1

! ELEMENTS D'EQUATIONS DANS BADV, GIADV ET ZADV
  real(kindreal), save:: ag1,ag1t,ag1t1,ag1a,ag1a1,ag1u,ag1u1,ag1o,ag1o1,ag2,ag2t, &
    ag2t1,ag2a,ag2a1,ag2u,ag2u1,ag2o,ag2o1,ag3,ag3t,ag3t1,ag3a,ag3a1,ag3u,ag3u1,ag3o, &
    ag3o1,ag4,ag4t,ag4t1,ag4a,ag4a1,ag4u,ag4u1,ag4o,ag4o1
  real(kindreal), save:: b1,b1o,b1o1,b1a,b1a1,b1u,b1u1,b1t,b1t1
  real(kindreal), save:: az1,az1o,az1o1,az1a,az1a1,az1u,az1u1,az1t,az1t1

end module equadiffmod
