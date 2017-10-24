module bintidemod
! module computing the tidal forces if the star is a binary.
! the formulas are explained in Song+ 2013 (A&A 556, A100)
! Here the period and orbital angular momentum are kept constant, but the orbital parameters
! are computed in subroutine orbitalevol, so it is possible to compute the 'real' situation
! modulo some light changes

  use evol,only: kindreal
  use const,only: pi,cst_G,cst_sigma,Lsol,Rsol,Msol,year,day
  use inputparam,only: binm2,periodini,verbose

  implicit none

  real(kindreal):: ab,romorb,romini,rstar,period


private
public:: period,dLtidcalc

contains
!======================================================================
subroutine dLtidcalc(dLtid)
!----------------------------------------------------------------------
! formalism taken from Zahn (1977), A&A 57, 383
! with the tidal coefficient E_2 as in Yoon (2010), ApJ 725, 940
!----------------------------------------------------------------------
  use inputparam,only: const_per,isol
  use caramodele,only: gms,gls,teff,nwmd
  use rotmod,only: omegi,bmomit,vsuminenv
  use convection,only: r_core
  use strucmod,only: vnr,vna,r
  use timestep,only: alter,dzeit

  implicit none

  real(kindreal),intent(out):: dLtid
  real(kindreal):: spiper,regra,qr1,qr2,rltid,fact1,fact2,fact3,fact4,e2,tsyn,roint,&
       rroche,e21,tsyn1,fact5,s22,dltido,fas1,fas2,&
       difsav,trot
!----------------------------------------------------------------------
  ab=(cst_G*(binm2+gms)*Msol*period**2.d0/(4.d0*pi**2.d0))**(1.d0/3.d0)
  romorb=2.d0*pi/period
  rstar=sqrt(gls*Lsol/(4.d0*pi*cst_sigma))/(teff**2.d0)
  qr1=binm2/gms
  qr2=1.0d0/qr1
  regra=gms*Msol*rstar**2.d0/bmomit
  fact1=5.0d0*2.d0**(5.d0/3.d0)*regra
  fact2=sqrt(cst_G*gms*Msol/rstar**3.d0)

! Roche radius according to Eggleton 1983
  rroche=ab*0.49d0*qr2**(2.d0/3.d0)/(0.6d0*qr2**(2.d0/3.d0)+log(1.d0+qr2**(1.d0/3.d0)))

  fact3=qr1**2.d0*(1.d0+qr1)**(5.d0/6.d0)
  fact4=(rstar/ab)**8.50d0
  fact5=abs(1.d0-romorb/omegi(1))

! e2: tidal coefficient, depends on the structure of the star,
!     especially the radius of the convective core (see Yoon et al. 2010)
  e2=10.d0**(-1.37d0)*(r_core/rstar)**8.d0
  e21=1.592d-9*gms**2.84d0
  tsyn=1.0d0/(fact1*fact2*fact3*fact4*e2)
  tsyn1=1.0d0/(fact1*fact2*fact3*fact4*e21)
  roint=gms*Msol*rstar**2.d0
  spiper=2.d0*pi/(omegi(1)*day)

  s22=2.d0*abs(omegi(1)-romorb)/fact2
  if (abs(period/day-spiper)>2.0d-4)then
    if (vnr-vna <0.0d0) then
      if (verbose) then
        write(*,*) 'radiative envelope',vnr,vna
      endif
      difsav=romorb-omegi(1)
      fas1=cst_G*(gms*Msol)**2.d0/rstar
      fas2=qr1**2.d0*(rstar/ab)**6.d0
      dltid=1.5d0*fas1*e2*sign(fas2,difsav)*s22**(8.d0/3.d0)*dzeit
      rltid=dltid/dzeit
      trot=1.0d0/(3.d0*fact2*regra*e2*fas2*s22**(5.d0/3.d0))
      dltido=vsuminenv*(romorb-omegi(1))/tsyn*dzeit
      if (verbose) then
        write(*,*)'previous tidtor',dltido,'current tidtor',dltid,'trot',trot/year
      endif
    else
      if (verbose) then
        write(*,*) 'convective envelope',vnr,vna
      endif
      tsyn=((1.d0/qr1**2.d0)*(ab/rstar)**6.d0)*year
      dltido=vsuminenv*(romorb-omegi(1))/tsyn*dzeit
      rltid=vsuminenv*(romorb-omegi(1))/tsyn
    endif
  else
    if (vnr-vna <0.0d0) then
      if (verbose) then
        write(*,*) 'radiative envelope',vnr,vna
      endif
      dltid=0.0d0
      rltid=0.0d0
    else
      if (verbose) then
        write(*,*) 'convective envelope',vnr,vna
      endif
      dltid=0.0d0
      rltid=0.0d0
    endif
  endif

  if (.not. const_per .and. isol==0) then
    call orbitalevol(dLtid)
  endif
  if (rstar/rroche > 1.01d0)then
    rewind(222)
    write(222,*) nwmd,'the star overfills the roche lobe'
      stop 'the star overfills the roche lobe'
  endif
  if (verbose) then
    write(*,*)'Model',nwmd,'alter=',alter,'period=',period/day,'delta t=',dzeit/year,'spin period=',&
         spiper,'synch timescale1=',tsyn/year,'synch timescale2=',tsyn1/year,'rstar/rroche=',rstar/rroche,&
         'rcon/rstar=',r_core/rstar,'rcon1/Rsol=',r_core/Rsol,'1.0/regra=',1.0d0/regra,'gms=',gms,&
         'gyr/rstar=', exp(r(1))/rstar,'e2=',e2,'e21=',e21
  endif
  write(3,*)'Model',nwmd,'alter=',alter,'period=',period/day,'delta t=',dzeit/year,'spin period=',&
         spiper,'synch timescale1=',tsyn/year,'synch timescale2=',tsyn1/year,'rstar/rroche=',rstar/rroche,&
         'rcon/rstar=',r_core/rstar,'rcon1/Rsol=',r_core/Rsol,'1.0/regra=',1.0d0/regra,'gms=',gms,&
         'gyr/rstar=', exp(r(1))/rstar,'e2=',e2,'e21=',e21

  return

end subroutine dLtidcalc
!======================================================================
subroutine orbitalevol(dLtid)
!----------------------------------------------------------------------
  use caramodele, only: gms,dm_lost,nwmd
  use rotmod,only: omegi
  use timestep,only: alter,dzeit

  implicit none

  real(kindreal),intent(in)::  dLtid
  real(kindreal):: dab,qr1,qr2,orbang,orstwi,rorstw,fact1,fact2,fact3,fact4,rltid
!----------------------------------------------------------------------
  qr1=binm2/(binm2+gms)
  qr2=1.0d0/qr1
  orbang=(gms*Msol*(qr1*ab)**2.d0+binm2*Msol*((1.d0-qr1)*ab)**2.d0)*romorb
  orstwi=dm_lost*Msol*(qr1*ab)**2.d0*romorb
  rorstw=orstwi/dzeit
  fact1=-2.d0*dLtid/orbang
  rltid=dLtid/dzeit
  fact2=-2.d0*dm_lost/gms
  fact3=dm_lost/(binm2+gms)
  fact4=2.d0*orstwi/orbang
  dab=ab*(fact1+fact2+fact3+fact4)
  ab=ab+dab
  romorb=(cst_G*(binm2+gms)*Msol/ab**3.0d0)**0.5d0
  period=2.0*pi/romorb
  if (verbose) then
    write(*,*) 'alter=',alter,'dltid=',dltid,'orstwi=',orstwi,'romini=',romini,'rmw=',dm_lost/(dzeit/year),&
       'fact1=',fact1,'fact2=',fact2,'fact3=',fact3,'dab=', dab/ab,'binm2=',binm2,'rstar/rsun=',rstar/Rsol,&
       'ab*qr1/rsun=',ab*qr1/Rsol
  endif
  write(3,*) 'alter=',alter,'dltid=',dltid,'orstwi=',orstwi,'romini=',romini,'rmw=',dm_lost/(dzeit/year),&
       'fact1=',fact1,'fact2=',fact2,'fact3=',fact3,'dab=', dab/ab,'binm2=',binm2,'rstar/rsun=',rstar/Rsol,&
       'ab*qr1/rsun=',ab*qr1/Rsol
  write(3,'(a11,e12.6)') 'new period=',period/day
  write(81,'(i7,1x,e22.15,1x,f9.3,2(1x,es13.6))'),nwmd,alter,rstar/Rsol,omegi(1)/day,period/day

  return

end subroutine orbitalevol
!======================================================================
end module bintidemod
