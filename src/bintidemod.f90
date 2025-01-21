module bintidemod
! module computing the tidal forces if the star is a binary.
! the formulas are explained in Song+ 2013 (A&A 556, A100) and in Sciarini+ 2024
! Module updated by Luca Sciarini in 2024 (luca.sciarini@unige.ch) to include eccentric orbits


! Here the period and orbital angular momentum are kept constant, but the orbital parameters
! are computed in subroutine orbitalevol, so it is possible to compute the 'real' situation
! modulo some light changes

  use io_definitions
  use evol,only: kindreal
  use const,only: pi,cst_G,cst_sigma,Lsol,Rsol,Msol,year,day
  use inputparam,only: binm2,periodini,verbose,eccentricity_ini,ie2_prescription

  implicit none

  real(kindreal):: ab,romorb,romini,rstar,period,eccentricity


private
public:: period,dLtidcalc,eccentricity

contains
!======================================================================
subroutine dLtidcalc(dLtid)
!----------------------------------------------------------------------
! formalism taken from Zahn (1977), A&A 57, 383 and Sciarini et al. (2024), A&A 681, L1 for eccentric orbits
! The dynamical tides model consists of an expansion valid for low eccentricities, not recommended for very
! eccentric orbits, of say e >~ 0.5. 
! The set of equations in the case of eccentric orbits is the Eq. (9) of the letter.
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
       rroche,fact5,s10,s12,s22,s32,dltido,fas1,fas2,&
       difsav,difsav10,difsav12,difsav32,trot
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

  select case(ie2_prescription)
    case(0)
      !Qin 2018 (recommended) prescription, used in Sciarini+24
      e2=10.d0**(-0.42d0)*(r_core/rstar)**(15.d0/2.d0)
    case(1)
      !Yoon et al. 2010 prescription, used in Song+13
      e2=10.d0**(-1.37d0)*(r_core/rstar)**8.d0
    case(2)
      !Hurley 2002 prescription, only mass dependent (not recommended).
      e2=1.592d-9*gms**2.84d0
  end select
  tsyn=1.0d0/(fact1*fact2*fact3*fact4*e2)
  roint=gms*Msol*rstar**2.d0
  spiper=2.d0*pi/(omegi(1)*day)

  s22=2.d0*abs(omegi(1)-romorb)/fact2
  s12=abs(romorb-2.d0*omegi(1))/fact2
  s32=abs(3.d0*romorb-2.d0*omegi(1))/fact2
  s10=abs(romorb)/fact2
  if (abs(period/day-spiper)>2.0d-4)then
    if (vnr-vna <0.0d0) then
      if (verbose) then
        write(*,*) 'radiative envelope',vnr,vna
      endif
      difsav=romorb-omegi(1)
      difsav12=romorb-2.d0*omegi(1)
      difsav32=3.d0*romorb-2.d0*omegi(1)
      difsav10=romorb
      fas1=cst_G*(gms*Msol)**2.d0/rstar
      fas2=qr1**2.d0*(rstar/ab)**6.d0
      
      ! Correct expression for the tidal torque in circular orbits. Formula used eg. in Song+13,16,18
      ! dltid=1.5d0*fas1*e2*sign(fas2,difsav)*s22**(8.d0/3.d0)*dzeit
      
      ! For eccentric orbits (general case), we use equations 9 in Sciarini+24
      ! dltid obtained with first equation, eccentricity evolution given by the second one (implemented in subroutine orbitalevol)
      
      dltid=1.5d0*fas1*e2*fas2*(sign(s22**(8.d0/3.d0),difsav)+eccentricity**2.d0*((1.d0/4.d0)*sign(s12**(8.d0/3.d0),difsav12)-&
      5.d0*sign(s22**(8.d0/3.d0),difsav)+(49.d0/4.d0)*sign(s32**(8.d0/3.d0),difsav32)))*dzeit
      
      rltid=dltid/dzeit
      trot=1.0d0/(3.d0*fact2*regra*e2*fas2*s22**(5.d0/3.d0))
      
      ! Incorrect expression for the dynamical tide, which has however been used a lot in the literature
      ! (Binstar, MESA, Posydon, Tres, ...). We don't recommend this expression, and it is not used in this code
      ! dltido is anyway not used further in the code.
      ! It is probably still in GENEC because in the first GENEC paper with tides (Song+13), the mistake was originally
      ! made, then corrected by Zahn. But this expression is not correct. See Sciarini+24 for more details.
      ! We might as well remove it at this point
      dltido=vsuminenv*(romorb-omegi(1))/tsyn*dzeit
      
      if (verbose) then
        write(*,*) 'previous tidtor',dltido,'current tidtor',dltid,'trot',trot/year
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
    rewind(io_runfile)
    write(io_runfile,*) nwmd,'the star overfills the roche lobe'
      stop 'the star overfills the roche lobe'
  endif
  if (verbose) then
    write(*,*) 'Model',nwmd,'alter=',alter,'period=',period/day,'delta t=',dzeit/year,'spin period=',&
         spiper,'synch timescale=',tsyn/year,'rstar/rroche=',rstar/rroche,&
         'rcon/rstar=',r_core/rstar,'rcon1/Rsol=',r_core/Rsol,'1.0/regra=',1.0d0/regra,'gms=',gms,&
         'gyr/rstar=', exp(r(1))/rstar,'e2=',e2
  endif
  write(io_logs,*) 'Model',nwmd,'alter=',alter,'period=',period/day,'delta t=',dzeit/year,'spin period=',&
         spiper,'synch timescale=',tsyn/year,'rstar/rroche=',rstar/rroche,&
         'rcon/rstar=',r_core/rstar,'rcon1/Rsol=',r_core/Rsol,'1.0/regra=',1.0d0/regra,'gms=',gms,&
         'gyr/rstar=', exp(r(1))/rstar,'e2=',e2

  return

end subroutine dLtidcalc
!======================================================================
subroutine orbitalevol(dLtid)
!----------------------------------------------------------------------
  !Module updated to correctly evolve the orbits in eccentric systems (Sciarini+24)
  use caramodele, only: gms,gls,dm_lost,nwmd,teff
  use rotmod,only: omegi
  use timestep,only: alter,dzeit
  
  use convection,only: r_core

  implicit none

  real(kindreal),intent(in)::  dLtid
  real(kindreal):: dab,qr1,qr2,qtot,qtot_inv,orbang,orstwi,rorstw,term1,term2,term3,term4,term_eccentricity,deccentricity,rltid,&
  s10,s12,s22,s32,difsav,difsav10,difsav12,difsav32,e2,fact2
!----------------------------------------------------------------------
  rstar=sqrt(gls*Lsol/(4.d0*pi*cst_sigma))/(teff**2.d0)
  fact2=sqrt(cst_G*gms*Msol/rstar**3.d0)
  
  s22=2.d0*abs(omegi(1)-romorb)/fact2
  s12=abs(romorb-2.d0*omegi(1))/fact2
  s32=abs(3.d0*romorb-2.d0*omegi(1))/fact2
  s10=abs(romorb)/fact2
  
  difsav=romorb-omegi(1)
  difsav12=romorb-2.d0*omegi(1)
  difsav32=3.d0*romorb-2.d0*omegi(1)
  difsav10=romorb
  
! e2: tidal coefficient, depends on the structure of the star,
!     especially the radius of the convective core (see Yoon et al. 2010)

  select case(ie2_prescription)
    case(0)
      !Qin 2018 (recommended) prescription, used in Sciarini+24
      e2=10.d0**(-0.42d0)*(r_core/rstar)**(15.d0/2.d0)
    case(1)
      !Yoon et al. 2010 prescription, used in Song+13
      e2=10.d0**(-1.37d0)*(r_core/rstar)**8.d0
    case(2)
      !Hurley 2002 prescription, only mass dependent (not recommended).
      e2=1.592d-9*gms**2.84d0
  end select
  
  qtot=binm2/(binm2+gms)
  qtot_inv=1.0d0/qtot
  qr1 = binm2/gms
  qr2 = 1.0d0/qr1
  orbang=(gms*Msol*(qtot*ab)**2.d0+binm2*Msol*((1.d0-qtot)*ab)**2.d0)*romorb
  orstwi=dm_lost*Msol*(qtot*ab)**2.d0*romorb
  rorstw=orstwi/dzeit
  term1=-2.d0*dLtid/orbang
  rltid=dLtid/dzeit
  term2=-2.d0*dm_lost/gms
  term3=dm_lost/(binm2+gms)
  term4=2.d0*orstwi/orbang
  
  ! deccentricity computed using equation 9 in Sciarini+24
  deccentricity=-(3.d0/4.d0)*eccentricity*fact2*qr1*sqrt(1+qr1)*e2*(rstar/ab)**(13.d0/2.d0)*((3.d0/2.d0)*&
  sign(s10**(8.d0/3.d0),difsav10)-(1.d0/4.d0)*sign(s12**(8.d0/3.d0),difsav12)-sign(s22**(8.d0/3.d0),difsav)+(49.d0/4.d0)*&
  sign(s32**(8.d0/3.d0),difsav32))*dzeit
  term_eccentricity=(2.d0*eccentricity*deccentricity)/(1-eccentricity**2.d0)
  eccentricity=eccentricity+deccentricity
  
  ! Variation of separation in the circular case
  !dab=ab*(term1+term2+term3+term4)
  
  
  ! Variation of separation in the eccentric case (general case)
  dab=ab*(term1+term2+term3+term4+term_eccentricity)
  ab=ab+dab
  romorb=(cst_G*(binm2+gms)*Msol/ab**3.0d0)**0.5d0
  period=2.0*pi/romorb
  write(*,*)'semi major axis: a = ',ab
  write(*,*)'eccentricity= ',eccentricity
  if (verbose) then
    write(*,*) 'alter=',alter,'dltid=',dltid,'orstwi=',orstwi,'romini=',romini,'rmw=',dm_lost/(dzeit/year),&
       'term1=',term1,'term2=',term2,'term3=',term3,'dab=', dab/ab,'binm2=',binm2,'rstar/rsun=',rstar/Rsol,&
       'ab*qtot/rsun=',ab*qtot/Rsol
  endif
  write(io_logs,*) 'alter=',alter,'dltid=',dltid,'orstwi=',orstwi,'romini=',romini,'rmw=',dm_lost/(dzeit/year),&
       'term1=',term1,'term2=',term2,'term3=',term3,'dab=', dab/ab,'binm2=',binm2,'rstar/rsun=',rstar/Rsol,&
       'ab*qtot/rsun=',ab*qtot/Rsol
  write(io_logs,'(a11,e12.6)') 'new period=',period/day
  write(io_period_evol,'(i7,1x,e22.15,1x,f9.3,5(1x,es13.6))') nwmd,alter,rstar/Rsol,omegi(1),romorb,period/day,ab,eccentricity

  return

end subroutine orbitalevol
!======================================================================
end module bintidemod

