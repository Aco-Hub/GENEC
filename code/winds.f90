!======================================================================
!> Module computing all aspects of mass loss
!!
!!  @version 211
!!  @date 14.2.2013
!======================================================================
module winds

  use evol,only: kindreal
  use inputparam,only: verbose

  implicit none

private
public :: aniso,corrwind,xloss,xldote

contains
!======================================================================
!> Computation of the radiative mass loss
!! @param[in] WRNoJump (skips the bistability jump in case of WR winds)
!! @param[out] checkVink (checks whether IMLOSS 6 > IMLOSS 7 or 8)
!!
!! @brief Computes the radiative mass loss according to the following recipes given the value of IMLOSS
!!        1. de Jager et al. (1988) and Sylvester (1998), van Loon 1999 for RSG (cf Crowther 2000)
!!        2. mass loss in Msol/yr given by FMLOS
!!        3. Reimers formula with etaR given by FMLOS
!!        4. WR mass loss : as in papier V
!!        5. for log(Teff) > 3.95 and Mini > 15 Msol: Kudritzki et Puls (2000);
!!           otherwise switches to IMLOSS=1
!!        6. for log(Teff) > 3.9 and Mini > 15 Msol: Vink et al (2001);
!!           otherwise switches to IMLOSS=1
!!        7. WR mass loss : Nugis et Lamers 2000 with Z-dependance as in Eldridge & Vink (2006)
!!        8. WR mass loss : Graefener and Hamann (2008);
!!           if not in the validity domain, switches to IMLOSS=7
!!        9. for log(Teff) > 4.6, very low Z and Mini >= 100 Msol : Kudritzki (2002)
!!        10. for RSG and AGB : van Loon & al. (2005)
!!        11. for log(Teff) > 3.9 and Mini > 15 Msol: Vink et al (2001) modified by
!!           Markova & Puls (2008) + priv. comm. Puls (nov. 2010)
!!           otherwise switches to IMLOSS=1
!----------------------------------------------------------------------
subroutine xloss(checkVink,WRNoJump)
!----------------------------------------------------------------------
  use const, only: lgLsol,lgpi,cstlg_sigma,lgRsol,cst_thomson,cst_avo,xlsomo,qapicg,cst_G,Msol,Rsol, &
                   Lsol,cst_sigma,pi,year
  use inputparam, only: ipop3,zsol,imloss,zinit,fmlos,irot,B_initial,frein
  use caramodele, only: teff,gls,iwr,xmini,eddesc,gms,xmdot,teffv,nwmd,zams_radius,Mdot_NotCorrected
  use strucmod, only: m
  use abundmod, only: x,y,y3,xc12,xo16,xn14
  use rotmod, only: alpro6,omegi,vomegi

  implicit none
  logical,intent(in):: WRNoJump
  logical,intent(out):: checkVink

  integer:: imlosscalc
  real(kindreal):: xteff,ygls,zheavy,zlim,xlgfz,xlogz,zeta,gbetaz,ggam0,xxx,yyy,t2x,t2y,t3x,t3y,t4x,t4y,t5x,dotm, &
    xrsol,xmdotn,xmdvir,xqhe,xepsi,xsigme,xgame,xmasef,xlgrrs,xrrs,xvescp,xvinfi,charrho,teffjump1,teffjump2,ratio, &
    xlmdot,gtest,gledd,xteffcond,azs,als,als2,azmin,aqmin,aq0,aq1,xxtt,xxll,teffjump,rstar,Bsurf,v_inf,v_esc,r_K,Correction_factor

  real(kindreal),parameter:: gram0=-3.763d0
! constants used for Jager et al 1988
  real(kindreal),parameter:: a00=6.34916d0,a01=-5.04240d0,a02=-0.83426d0,a03=-1.13925d0,a04=-0.12202d0,a10=3.41678d0, &
    a11=0.15629d0,a12=2.96244d0,a13=0.33659d0,a14=0.57576d0,a20=-1.08683d0,a21=0.41952d0,a22=-1.37272d0,a23=-1.07493d0, &
    a30=0.13095d0,a31=-0.09825d0,a32=0.13025d0,a40=0.22427d0,a41=0.46591d0,a50=0.11968d0

  rstar = 0.d0
  Bsurf = 0.d0
  v_inf = 0.d0
  v_esc = 0.d0

!----------------------------------------------------------------------
  xteff=log10(teff)
  ygls=log10(gls)
  write (3,*) 'xteff= ',xteff
  ! calcul de la dependance en metallicite
  ! le Log facteur=xlgfz
  ! ce facteur est utilise partout dans le DHR
  ! sauf si l'etoile est WR
  zheavy=1.d0-x(1)-y(1)-y3(1)
  if (ipop3 == 1) then
    zlim=1.d-04*zsol
    if (zheavy < zlim) zheavy=zlim
  endif
!1b calcul de Log f(Z)
  if (zheavy > (zsol-1.0d-3).and.zheavy < (zsol+1.0d-3)) then
    xlgfz=0.0d0
  else
    xlgfz=0.5d0*log10(zheavy/zsol)
  endif

  iwr=0
  select case (imloss)
    case (4,7,8)
      iwr=1
  end select

  if (ipop3 == 1) then
    xlogz=log10(zheavy/zsol)
  else
    if (iwr  ==  1) then
      xlogz=log10(zinit/zsol)
    else
      xlogz=log10(zheavy/zsol)
    endif
  endif
  zeta=(xc12(1)+xo16(1))/y(1)

  select case (imloss)
  case (1)
    imlosscalc = 1
  case (2)
    imlosscalc = 2
    if (x(1) < 0.3d0.and.xteff > 4.d0) iwr=1
  case (3)
    imlosscalc = 3
  case (4)
    imlosscalc = 4
  case (5)
    if (xmini > 15.d0 .and. xteff >= 3.95d0) then
      imlosscalc = 5
    else
      imlosscalc = 1
    endif
  case (6)
    if (xmini > 15.d0 .and. xteff >= 3.90d0) then
      imlosscalc = 6
      if (x(1) < 0.3d0.and.xteff > 4.d0) iwr=1
    else
      imlosscalc = 1
    endif
  case (7)
    imlosscalc = 7
  case (8)
    if ((xteff < 4.477d0.or.xteff > 4.845d0) .or.(xlogz < -3.d0.or.xlogz > 0.30d0)) then
      write(*,*) 'GRAEF transferred to Nugis (Teff,Z)'
      write(3,*) 'GRAEF transferred to Nugis (Teff,Z)'
      imlosscalc = 7
    else
      gbetaz=1.727d0+0.25d0*xlogz
      ggam0=0.326d0-0.301d0*xlogz-0.045d0*xlogz*xlogz
      gledd=log10(eddesc-ggam0)
      if (eddesc <= ggam0) then
        write(*,*) 'GRAEF transferred to Nugis (Teff,Z)'
        write(3,*) 'GRAEF transferred to Nugis (Teff,Z)'
        imlosscalc = 7
      else
        imlosscalc = 8
      endif
    endif
  case (9)
    imlosscalc = 9
  case (10)
    imlosscalc = 10
  case (11)
    if (xmini > 15.d0 .and. xteff >= 3.90d0) then
      imlosscalc = 11
    else
      imlosscalc = 1
    endif
  case default
    stop 'Bad IMLOSS value, must be between 1 - 10'
  end select

!=======================================================================
! Calcul de la perte de masse
  select case (imlosscalc)
!-----------------------------------------------------------------------
  case (1)
!***de Jager et al 88 est pris pour log Teff plus grand que 3.7
    if (xteff > 3.7d0) then
      xxx = (xteff-4.05d0)/0.75d0
      yyy = (ygls-4.6d0)/2.1d0
      t2x = cos(2.d0*acos(xxx))
      t2y = cos(2.d0*acos(yyy))
      t3x = cos(3.d0*acos(xxx))
      t3y = cos(3.d0*acos(yyy))
      t4x = cos(4.d0*acos(xxx))
      t4y = cos(4.d0*acos(yyy))
      t5x = cos(5.d0*acos(xxx))
      dotm = a00+a01*yyy+a10*xxx+a02*t2y+a11*xxx*yyy+a20*t2x+a03*t3y+a12*xxx*t2y+a21*t2x*yyy+a30*t3x+a04*t4y+ &
             a13*xxx*t3y+a22*t2x*t2y+a31*t3x*yyy+a40*t4x+a14*xxx*t4y+a23*t2x*t3y+a32*t3x*t2y+a41*t4x*yyy+a50*t5x
      xmdot = 10.d0**(-dotm)*10.d0**xlgfz
    else   ! xteff <= 3.7
!*** taux propose par Maeder sur la base des figures dans le
! papier de Crowther (2000), observations de
! Sylvester et al 1998 et van Loon et al. (LMC) 1999
! NB: pas de dependance en Z quand Teff <= 3.7
      dotm = -(1.7d0*ygls-13.83d0)
      xmdot = 10.d0**(-dotm)
    endif   ! xteff
!-----------------------------------------------------------------------
  case (2)
!*** perte de masse donnee en Msol/yr par FMLOS
    xmdot=1.d0
!-----------------------------------------------------------------------
  case (3)
!*** formule de Reimers, etaR donne par fmlos
! rayon en unite de rayon solaire:
    xrsol = 0.5d0*(log10(gls)-4.d0*log10(teff)+lgLsol-log10(4.d0)-lgpi-cstlg_sigma)-lgRsol
    xrsol = 10.d0**xrsol
    xmdot = 4.0d-13*gls*xrsol/gms
!-----------------------------------------------------------------------
  case (4)
    if (x(1) > 1.0d-3) then
! CAS DES ETOILES WNL cf Nugis et al. 1998 AA,333,956
      xmdot = 3.0d-05
      xmdotn = 2.4d-08*(gms**2.5d0)
      xmdot = min(xmdot,xmdotn)
    else
      if (xc12(1) <= xn14(1)) then
! CAS DES ETOILES WNE, cf Schmutz 1997
        xmdot = 2.4d-08*(gms**2.5d0)
      else
! CAS DES ETOILES WC, cf Schmutz 1997
        xmdot = 2.4d-08*(gms**2.5d0)
      endif
    endif
!-----------------------------------------------------------------------
  case (5)
!*** Kudritzki et Puls (2000)
    if (x(m) > 0.d0) then
      xmdvir = 19.87d0+1.57d0*ygls+xlgfz-30.799531d0
    else
      if (xteff > 4.5d0) then
        xmdvir = 20.69d0+1.51d0*ygls+xlgfz-30.799531d0
      else if (xteff <= 4.5d0.and.xteff > 4.4d0) then
        xmdvir = 21.24d0+1.34d0*ygls+xlgfz-30.799531d0
      else if (xteff <= 4.4d0.and.xteff > 4.275d0) then
        xmdvir = 17.07d0+1.95d0*ygls+xlgfz-30.799531d0
      else
        xmdvir = 14.22d0+2.64d0*ygls+xlgfz-30.799531d0
      endif
    endif
!2 calcul de la masse effective
!2a calcul de qHe
    if (teff >= 35000.d0) then
      xqhe = 2.0d0
    else if (teff >= 30000.d0.and.teff < 35000.d0) then
      xqhe = 1.5d0
    else if (teff >= 25000.d0.and.teff < 30000.d0) then
      xqhe = 1.0d0
    else
      xqhe = 0.0d0
    endif
!2b calcul de epsilon
    xepsi = y(1)/(4.d0*x(1)+y(1))
!2c calcul de Sigmae
    xsigme = cst_thomson*cst_avo*(1.d0+(xqhe-1.d0)*xepsi)/(1.d0+3.d0*xepsi)
!2d calcul de Gammae
    xgame = xlsomo*xsigme*(gls/gms)/qapicg
!2e calcul de la masse effective
    xmasef = gms*(1.d0-xgame)
!3 calcul de Vinf
!3a calcul de Rstar
    xlgrrs = 0.5d0*(ygls-4.d0*xteff+lgLsol-log10(4.d0)-lgpi-cstlg_sigma-2.d0*lgRsol)
    xrrs = 10.d0**xlgrrs
!3b calcul de Vescp en km/s.
    xvescp = sqrt(2.d0*cst_G*Msol/Rsol)*sqrt(xmasef/xrrs)/1.d5
!3c calcul de Vinf
    if (teff >= 21000.d0) then
      xvinfi = 2.65d0*xvescp
    else if (teff <= 10000.d0) then
      xvinfi = xvescp
    else
      xvinfi = 1.4d0*xvescp
    endif
    xmdot = 10.d0**(xmdvir)/(xvinfi*sqrt(xrrs))
!-----------------------------------------------------------------------
  case (6)
!*** Vink et al (2001)
    charrho = -14.94d0+3.1857d0*eddesc+0.85d0*xlogz
    teffjump1 = 61.2d0+2.59d0*charrho
    teffjump1 = teffjump1*1000.d0
    teffjump2 = 100.d0+6.d0*charrho
    teffjump2 = teffjump2*1000.d0

    if (teffjump1 <= teffjump2) then
      write(*,*) ' These stellar parameters are unrealistic'
      checkVink = .false.
    endif

    if (.not. WRNoJump) then
      if (teff < teffjump1) then
        if (teffv >= teffjump1) then
         write(3,'(a,f8.5)')'XLOSS - Teff_jump,1 reached',log10(teffjump1)
         write(*,'(a,f8.5)')'XLOSS - Teff_jump,1 reached',log10(teffjump1)
         write(997,'(i7.7,a,f8.5)')nwmd,': Teff_jump,1 reached',log10(teffjump1)
        endif
        if (teff < teffjump2) then
          if (teffv >= teffjump2) then
            write(3,'(a,f8.5)')'XLOSS - Teff_jump,2 reached',log10(teffjump2)
            write(*,'(a,f8.5)')'XLOSS - Teff_jump,2 reached',log10(teffjump2)
            write(997,'(i7.7,a,f8.5)')nwmd,': Teff_jump,2 reached',log10(teffjump2)
          endif
          ratio = 0.7d0
          xlmdot = -5.99d0+2.210d0*log10(gls/1.0d+5)-1.339d0*log10(gms/30.d0)-1.601d0*log10(ratio/2.d0)+ &
                   1.07d0*log10(teff/20000.d0)+0.85d0*xlogz
        else if (teff > teffjump2) then
          if (teffv <= teffjump2) then
            write(3,'(a,f8.5)')'XLOSS - Teff_jump,2 reached',log10(teffjump2)
            write(*,'(a,f8.5)')'XLOSS - Teff_jump,2 reached',log10(teffjump2)
            write(997,'(i7.7,a,f8.5)')nwmd,': Teff_jump,2 reached',log10(teffjump2)
          endif
          ratio = 1.3d0
          xlmdot = -6.688d0+2.210d0*log10(gls/1.0d+5)-1.339d0*log10(gms/30.d0)-1.601d0*log10(ratio/2.d0)+ &
                   1.07d0*log10(teff/20000.d0)+0.85d0*xlogz
        else
          stop ' STAR at the second jump'
        endif
      else if (teff > teffjump1) then
        if (teffv <= teffjump1) then
         write(3,'(a,f8.5)')'XLOSS - Teff_jump,1 reached',log10(teffjump1)
         write(*,'(a,f8.5)')'XLOSS - Teff_jump,1 reached',log10(teffjump1)
         write(997,'(i7.7,a,f8.5)')nwmd,': Teff_jump,1 reached',log10(teffjump1)
        endif
        ratio = 2.6d0
        xlmdot = -6.697d0+2.194d0*log10(gls/1.0d+5)-1.313d0*log10(gms/30.d0)-1.226d0*log10(ratio/2.d0)+ &
                 0.933d0*log10(teff/40000.d0)-10.92d0*(log10(teff/40000.d0))*(log10(teff/40000.d0))+0.85d0*xlogz
      else
        stop ' STAR at the first jump'
      endif
    else
      ratio = 2.6d0
      xlmdot = -6.697d0+2.194d0*log10(gls/1.0d+5)-1.313d0*log10(gms/30.d0)-1.226d0*log10(ratio/2.d0)+ &
               0.933d0*log10(teff/40000.d0)-10.92d0*(log10(teff/40000.d0))*(log10(teff/40000.d0))+0.85d0*xlogz

    endif
    xmdot = 10.d0**xlmdot
!-----------------------------------------------------------------------
  case (7)
!*** Nugis & Lamers 2000
! pour WN:
    if (x(1) > 0.d0.or.zeta <= 0.03d0) then
      xlmdot=-13.60d0+1.63d0*ygls+2.22d0*log10(y(1))+0.85d0*xlogz
! pour WC + WO:
    else
      if (zinit > zsol) then
        xlmdot=-8.30d0+0.84d0*ygls+2.04d0*log10(y(1))+1.04d0*log10(zheavy)+0.40d0*xlogz
      else if (zinit  >  0.002d0) then
        xlmdot=-8.30d0+0.84d0*ygls+2.04d0*log10(y(1))+1.04d0*log10(zheavy)+0.66d0*xlogz
      else
        if (ipop3 == 1) then
          xlmdot = -8.30d0+0.84d0*ygls+2.04d0*log10(y(1))+1.04d0*log10(zheavy)+0.66d0*log10(0.002d0/zsol)+ &
                   0.35d0*log10(zheavy/0.002d0)
        else
          xlmdot = -8.30d0+0.84d0*ygls+2.04d0*log10(y(1))+1.04d0*log10(zheavy)+0.66d0*log10(0.002d0/zsol)+ &
                   0.35d0*log10(zinit/0.002d0)
        endif
      endif
    endif
    xmdot=10.d0**xlmdot
!-----------------------------------------------------------------------
  case (8)
! formulae 3, 5, and 6 of Graefener & Hamann 2008 A&A 482,945
! http://ukads.nottingham.ac.uk/abs/2008A%26A...482..945G
    if (x(1) >= 0.05d0) then
      gtest=xteff-4.65d0
      xlmdot=gram0+gbetaz*gledd-3.5d0*gtest+0.42d0*(ygls-6.3d0)-0.45d0*(x(1)-0.4d0)
    else if (x(1) > 0.d0) then
! pour WNE:
      xlmdot=-13.60d0+1.63d0*ygls+2.22d0*log10(y(1))+0.85d0*xlogz
    else
      if (zeta <= 0.03d0) then
        xlmdot=-13.60d0+1.63d0*ygls+2.22d0*log10(y(1))+0.85d0*xlogz
! pour WC + WO:
      else
        if (zinit > zsol) then
          xlmdot=-8.30d0+0.84d0*ygls+2.04d0*log10(y(1))+1.04d0*log10(zheavy)+0.40d0*xlogz
        else if (zinit > 0.002d0) then
          xlmdot=-8.30d0+0.84d0*ygls+2.04d0*log10(y(1))+1.04d0*log10(zheavy)+0.66d0*xlogz
        else
          if (ipop3 == 1) then
            xlmdot = -8.30d0+0.84d0*ygls+2.04d0*log10(y(1))+1.04d0*log10(zheavy)+0.66d0*log10(0.002d0/zsol)+ &
                     0.35d0*log10(zheavy/0.002d0)
          else
            xlmdot = -8.30d0+0.84d0*ygls+2.04d0*log10(y(1))+1.04d0*log10(zheavy)+0.66d0*log10(0.002d0/zsol)+ &
                     0.35d0*log10(zinit/0.002d0)
          endif
        endif
      endif
    endif
    xmdot=10.d0**xlmdot
!-----------------------------------------------------------------------
  case (9)
!*** Kudritzki 2002
! on contraint logTeff au domaine de validite de la formule de Kudritzki
    if (xteff > 4.778d0) then
      xteffcond=4.778d0
    else
      xteffcond=xteff
    endif
! perte de masse selon Kudritzki 2002
    azs=log10(zheavy/zsol)
    als=ygls-6.d0
    als2=als*als
    if (xteffcond > 4.699d0.and.xteffcond <= 4.778d0) then
      azmin=-3.4d0-0.4d0*als-0.65d0*als2
      aqmin=-8.d0-1.2d0*als+2.15d0*als2
      aq0=-5.99d0+als+1.5d0*als2
    else if (xteffcond > 4.602d0) then
      azmin=-3.85d0-0.05d0*als-0.60d0*als2
      aqmin=-10.35d0+3.25d0*als
      aq0=-4.85d0+0.5d0*als+als2
    else
      azmin=-4.45d0+0.35d0*als-0.80d0*als2
      aqmin=-11.75d0+3.65d0*als
      aq0=-5.20d0+0.93d0*als+0.85d0*als2
    endif
    aq1=(aq0-aqmin)*((-azmin)**(-0.5d0))
    if ((azs-azmin) > 0.0d0) then
      xlmdot=aq1*((azs-azmin)**0.5d0)+aqmin
      xmdot=10.d0**xlmdot
    else
      xmdot = 0.d0
      write(*,*) 'IMLOSS 9: xmdot set to 0.'
      write(3,*) 'IMLOSS 9: azs-azmin<0 --> xmdot set to 0.'
    endif
!-----------------------------------------------------------------------
  case (10)
!*** van Loon & al. (2005) for RSG and AGB
    xxtt=log10(teff/3500.d0)
    xxll=log10(gls/10000.d0)
    if (ygls > 4.9d0) then
      xlmdot=-5.3d0+0.82d0*xxll-10.8d0*xxtt
    else
      xlmdot=-5.6d0+1.1d0*xxll-5.2d0*xxtt
    endif
    xmdot=10.d0**xlmdot
!-----------------------------------------------------------------------
  case (11)
!*** Vink et al (2001, IMLOSS 6) modified by Markova & Puls (2008) + priv. comm. Puls (nov. 2010)
    teffjump = 10.d0**(4.3d0)
    if (teff <= teffjump) then
      if (teffv > teffjump) then
       write(3,'(a,f8.5)')'XLOSS - Teff_jump reached, B -> R',log10(teffjump)
        write(*,'(a,f8.5)')'XLOSS - Teff_jump reached, B -> R',log10(teffjump)
        write(997,'(i7.7,a,f8.5)')nwmd,': Teff_jump reached, B -> R',log10(teffjump)
      endif
      ratio = 1.4d0
    else
      if (teffv <= teffjump) then
        write(3,'(a,f8.5)')'XLOSS - Teff_jump reached, R -> B',log10(teffjump)
        write(*,'(a,f8.5)')'XLOSS - Teff_jump reached, R -> B',log10(teffjump)
        write(997,'(i7.7,a,f8.5)')nwmd,': Teff_jump reached, R -> B',log10(teffjump)
      endif
      ratio=3.0d0
    endif

    xlmdot = -6.697d0+2.194d0*log10(gls/1.0d+5)-1.313d0*log10(gms/30.d0)-1.226d0*log10(ratio/2.d0)+ &
             0.933d0*log10(teff/40000.d0)-10.92d0*(log10(teff/40000.d0))*(log10(teff/40000.d0))+0.85d0*xlogz
    xmdot=10.d0**xlmdot
!-----------------------------------------------------------------------
  case default
    stop 'Bad IMLOSSCALC value. Problem with IMLOSS ??'
  end select
!=======================================================================
  xmdot=fmlos*xmdot
  write(3,*) 'fmlos= ',fmlos,'  xmdot= ',xmdot
  if(irot == 1) then
    if (alpro6 /= 0.d0) xmdot = alpro6*xmdot
    write(3,'(2x,a,f14.7,a,f11.7)') 'facteur du a la rotation=',alpro6,' Gamma el. sc.=',eddesc
  endif

  if (B_initial > 1.d-5 .and. zams_radius > 0.d0) then
! Note here that we neglect the deformation of the stellar surface. It is easy to account for it
! in case it would be needed by taking its value in aniso.
    rstar = sqrt(gls*Lsol/(4.d0*pi*cst_sigma))/teff**2.d0
    Bsurf = B_initial*(zams_radius/rstar)**2.d0
! For consistency, the magnetic field intensity accounted for in the computation of the
! angular momentum loss is adapted here as well.
    frein = Bsurf
    v_esc = sqrt(2.d0*cst_G*gms*Msol/rstar)
    if (teff >= 21000.d0) then
      v_inf = 2.65d0*v_esc
    else if (teff <= 10000.d0) then
      v_inf = v_esc
    else
      v_inf = 1.4d0*v_esc
    endif
    r_K = (cst_G*gms*Msol/omegi(1)**2.d0)**(1.d0/3.d0)
    Correction_factor = Compute_MagCorrection(rstar,xmdot*Msol/year,Bsurf,v_inf,r_K)
  else
    Correction_factor = 1.d0
  endif
  if (xmdot > Mdot_NotCorrected) then
     write(3,*) 'Mdot (no corrections) = ', xmdot
     Mdot_NotCorrected = xmdot
  endif
  xmdot = xmdot*Correction_factor
  write(3,*) 'fmlos= ',fmlos,'  xmdot= ',xmdot, 'Magnetic correction: ', Correction_factor

  return

end subroutine xloss

!=======================================================================
!> Computation of the angular-momentum loss linked to the radiative mass loss
!! @param[in]  dmdot (dm_lost in main)
!! @param[out] dmneed, dlelex, xldoex
!----------------------------------------------------------------------
subroutine xldote(dmdot,dmneed)
!----------------------------------------------------------------------
  use evol, only: npondcouche
  use const, only: Msol
  use inputparam, only: ianiso,modanf,rapcrilim,bintide,xcn,diff_only
  use caramodele, only: gms,rayequat,xltotbeg,firstmods,nwmd,inum
  use strucmod, only: q,r
  use rotmod, only: xlexcs,rapom2,omegi,dlelex,xldoex,bdotis,vsuminenv
  use timestep,only: dzeit
  use bintidemod,only: dLtidcalc

  implicit none

  real(kindreal), intent(out)::dmneed
  real(kindreal), intent(in):: dmdot

  integer:: jo
  real(kindreal):: DeltaMCG,omcrit,dLmag,dL_Kawaler,dLtid

  real(kindreal):: omlim, newomega, dLisotrop,dLmeca, xLe, qcorr, qnum, qdenom, dmneednum, dmneeddenom
  real(kindreal), dimension(npondcouche)::Li
!> facteur multiplicatif: sqrt(r_out/R_star) cf. Owocki (perte L par disque)
  real(kindreal), parameter:: factordisk = 1.d0, omega_min = 1.d-20
!----------------------------------------------------------------------

  dmneed= 0.0d0
  dLmeca = 0.0d0
  dLtid = 0.0d0

! Si l'anisotropie n'est pas prise en compte, xlexcs est nul.
  if (ianiso == 0) xlexcs= 0.d0

! masse de l'enveloppe= 1.-FITM
  DeltaMCG = gms*Msol*exp(q(1))

! En premier lieu: calcul de la perte de moment cinetique isotrope.
! Plutot que de considerer la surface comme spherique, on reprend ici le moment
! cinetique de la perte de masse calculee dans anisotrop et qui prend en compte
! la deformation de l'etoile.
  dLisotrop = bdotis

  if (rapom2 == 0.d0) then
! Dans le cas ou le rapport omega/omega crit est strictement nul,
! on attribue des valeurs tres superieures a omlim et omcrit afin
! de ne pas poser de problemes. Cela ne devrait arriver que lors des premiers
! modeles.
    omlim = 1.d0
    omcrit=2.d0*omegi(1)
    if (omegi(1) >= 3.d-20) then
      if (modanf > 2 .and. verbose) then
        write(*,*) ' !!!!! WARNING !!!!!'
         write(*,*) 'rapom2 = 0.0 in xldote. If this model is',' not the 1st, this is an error!'
      endif
      write(3,*) ' !!!!! WARNING !!!!!'
      write(3,*) 'rapom2 = 0.0 in xldote. If this model is',' not the 1st, this is an error!'
    endif
  else
! Si rapcrilim est nul, la correction n'est pas appliquee et son calcul
! est inutile et source de bug.
    if (rapcrilim > 1.d-5) then
      omlim= rapcrilim*omegi(1)/rapom2
      omcrit= omegi(1)/rapom2
    else
      omlim = 1.d0
      omcrit=2.d0*omegi(1)
    endif
    do jo=1,1!NPcoucheEff
     if (omegi(jo)  <=  0.d0) then
       omegi(jo) = omega_min
     endif
    enddo

! On calcule ici DANS LA SITUATION ACUTELLE DE STRUCTURE la masse perdue
! par l'equateur si on est au-dela de la vitesse critique (on SUPPOSE que
! cette masse sera approximativement correcte lorsque la structure aura
! un peu change).
! On commence par calculer le "nouvel omega" (omega que le systeme enveloppe
! + couche 1-"NPcoucheEff" aurait par conservation du moment cinetique (sans tenir compte de
! la diffusion) avec la masse perdue par les vents).
! dLisotrop = 2.d0/3.d0 *omegi(1) *dmdot*Msol *rayequat**2.d0

! Application de la correction anisotrope (elle est nulle si ianiso = 0) (ici, on prend
! un signe - afin que dlelex soit positif)
    dlelex = -dLisotrop * (1.d0 + xlexcs)

! Calcul des moments cinetiques des couches 1 et de l'enveloppe.
    xLe = vsuminenv*omegi(1)
    Li(1) = 2.d0/3.d0 * (exp(q(2)) - exp(q(1)))/2.d0 * gms * Msol * exp(2.d0*r(1)) * omegi(1)

    qnum = -dmdot / (gms*exp(q(1))) * xLe - dlelex
    qdenom = Li(1) + (gms*exp(q(1)) + dmdot)/(gms*exp(q(1)))*xLe

    qcorr = qnum/qdenom

    newomega = omegi(1)*(1.d0+qcorr)

! Si cette nouvelle vitesse de surface est plus grande que la vitesse critique,
! alors on calcul ici la masse qu'il faudra perdre a l'equateur pour la ramener
! a sa valeur critique.
    write(3,*)
    write(3,*) 'In xldote: newomega, omega limit: ',newomega,omlim
    if (newomega >= omlim) then
! Le detail de cette formule se trouve dans la documentation.
      dmneednum = Li(1) * (omlim/omegi(1)-1.d0)
      dmneednum = dmneednum + xLe*(omlim/omegi(1)*(1.d0 + dmdot/(gms*exp(q(1)))) - 1.d0) + dlelex
      dmneeddenom = -factordisk*omegi(1)*rayequat**2.d0 + xLe*omlim / (gms*exp(q(1))*Msol*omegi(1))
      dmneed = dmneednum / dmneeddenom
      if (omegi(1)  <  1.d-25) then
        dmneed = 0.d0
      endif
    endif

! dmneed doit etre positif. Si ce n'est pas le cas, message d'erreur et arret de l'execution.
    if (dmneed < 0.d0) then
      rewind(222)
      write (222,*) nwmd,': dmneed negative in xldote'
      stop 'dmneed negative in xldote. Aborting...'
    endif

! Si la vitesse de surface est superieure de 0.25% a la valeur maximale toleree
! (ou 1 le cas echeant), on multiplie par 1.5 la perte de masse equatoriale.
    write(3,*) 'dmneed = ', dmneed

    if (dmneed > 0.d0) then
      if (rapom2  <=  0.995d0 .and. rapcrilim  >  0.d0) then
        if (rapom2  >  (rapcrilim + 0.0025d0)) then
          dmneed = 2.0d0*dmneed
          write(*,*) '!!! WARNING: equatorial mass loss increased by a factor 2.0!'
          write(3,*) 'dmneed multiplied by 2.0. New dmneed = ',dmneed
        else if (rapom2  >  (rapcrilim + 0.005d0)) then
          dmneed = 4.d0*dmneed
          write(*,*) '!!! WARNING: equatorial mass loss increased by a factor 4.0!'
          write(3,*) 'dmneed multiplied by 4. New dmneed = ',dmneed
        endif
      else
        if (rapom2  >  1.05d0 .and. rapcrilim  >  0.d0) then
          dmneed = 10.d0*dmneed
          write(*,*) '!!! WARNING: equatorial mass loss increased by a factor 10!'
          write(3,*) 'dmneed multiplied by 10. New dmneed = ',dmneed
        else if (rapom2  >  1.d0 .and. rapcrilim  >  0.d0) then
          dmneed = 1.5d0*dmneed
          write(*,*) '!!! WARNING: equatorial mass loss increased by a factor 1.5!'
          write(3,*) 'dmneed multiplied by 1.5. New dmneed = ',dmneed
        endif
      endif
    endif

! Pour appliquer la correction du moment cinetique de la surface dans la routine
! tridog, il faut connaitre le moment cinetique emporte par les vents et la perte
! mecanique equatoriale.
    if (dmneed > 0.d0) then
      dLmeca = factordisk*omegi(1)*dmneed*rayequat**2.d0
    endif
  endif

! dlelex est le moment cinetique total que doit perdre l'etoile sous l'influence de la
! perte de masse (signe - ici afin que dlelex soit positif dans le cas d'une perte de masse).
  call dLmagcalc(dLisotrop,dLmag)
! The total spherical angular momentum loss is anyway computed in dLmagcalc.
! It returns simply -dLisotrop if magnetic braking is not accounted for.
! The anisotropic correction is for now only applied to dLisotrop. The compatibility
! bewteen wind anisotropy and magnetic braking is NOT VERIFIED on theoretical basis.

! Computation of the Kawaler loss of angular momentum (for low-mass stars). The value is negative
! in the subroutine, so we take a minus sign to have a positive value here.
  dL_Kawaler = -dLmagLM()

  if (bintide .and. nwmd/=1) then
    call dLtidcalc(dLtid)
    if (diff_only) then
! Possibility to compute the torque only when diffusion is applied.
! In that case, dLtid=0 when advection is computed (odd number timestep)
! and dLtid=dLtid*2 when diffusion is computed (even number timestep)
      if (mod(nwmd,2)==1) then
        dLtid=0.0d0
      else
        if (inum==0) then
          dLtid=dLtid*(1.d0+xcn)
        else
          dLtid=dLtid*2.d0
        endif
      endif
    endif
  endif

  dlelex = dLmag - dLisotrop*xlexcs + dLmeca + dL_Kawaler - dLtid
  if (.not.firstmods) then
    if (abs(dlelex) > 0.05d0*xltotbeg) then
      write(*,*) 'MORE THAN 5% of total angular momentum removed !'
      dlelex = dlelex*omega_min/omegi(1)
      omegi(1) = omega_min
    endif
  endif

  if (verbose) then
    write(*,*) 'dLiso, xlexcs, dLmeca,  dLaniso, dLmag, dLtide: '
    write(*,*) dLisotrop,xlexcs,dLmeca,dlelex,dLmag+dLisotrop,dLtid
  endif
  write(3,*) 'dLiso, xlexcs, dLmeca,  dLaniso, dL_Kawaler, dLtide: '
  write(3,*) dLisotrop,xlexcs,dLmeca,dlelex,dL_Kawaler,dLtid

! dmneed doit etre retournee en unites solaires, et negatif dans le cas d'une perte de masse:
  if (dmneed /= 0.d0) then
     dmneed = -dmneed / Msol
  endif

! calcul de Ldot_excess utilise comme condition au bord dans l'advection. Attention aux
! signes de dLisotrop et dlelex, qui sont inverses.
  xldoex = -(dlelex + dLisotrop)/dzeit

  return

end subroutine xldote
!======================================================================
!> Computation of the surface magnetic braking
!!
!! @brief Computes the surface magnetic braking according to ud-Doula & Owocki 2002 and ud-Doula et al. 2008:
!! \f$ \frac{dJ}{dt}=\frac{2}{3} \dot{M} \Omega R_\star^2 \left[ 0.29 + (\eta_\star+0.25)^{1/4} \right]^2 \f$
!! where \f$ \eta_\star=B_{eq}^2 R_\star^2/(\dot{M} V_\infty) \f$ is the magnetic confinement parameter
!!
!! @param[in]  dL_isotrop
!! @param[out] dLmag
!----------------------------------------------------------------------
subroutine dLmagcalc(dL_isotrop,dLmag)
!----------------------------------------------------------------------
  use const,only: pi,cst_G,cst_sigma,Lsol,Rsol,Msol,year
  use inputparam,only: frein
  use caramodele,only: gms,gls,teff,xmdot,Mdot_NotCorrected

  implicit none

  real(kindreal), intent(in):: dL_isotrop
  real(kindreal), intent(out):: dLmag

  real(kindreal):: rstar,mdot,vinf,vesc,etastar,Constant
!----------------------------------------------------------------------
  if (frein > 1.d-5) then
    rstar = sqrt(gls*Lsol/(4.d0*pi*cst_sigma))/(teff**2.d0)
    write(3,*) 'Mass loss used in dLmagcalc = ', Mdot_NotCorrected
    mdot = Mdot_NotCorrected*Msol/year
    vesc = sqrt(2.d0*cst_G*gms*Msol/rstar)
    if (teff >= 21000.d0) then
      vinf = 2.65d0*vesc
    else if (teff <= 10000.d0) then
      vinf = vesc
    else
      vinf = 1.4d0*vesc
    endif

    etastar = frein**2.d0*rstar**2.d0/(mdot*vinf)
    write(*,*) 'ETA XLDOTE: ', etastar
    if (verbose) then
      write(*,*)'FREIN: eta=',etastar
    endif

    Constant = 1.d0 - 0.25d0**(1.d0/4.d0)
! In case the mass-loss rate correction is used, according to ud-Doula, the loss of angular
! momentum should be computed with the non-corrected mass-loss rate.
! However, for consistency, dL_isotrop should not be modified (as it concerns the mass
! that is really removed from the star, and is thus correct if computed with the
! reduced mass-loss rate). In order to recover the angular momentum removed by the
! non-corrected mass-loss rate, we have to multiply dL_isotrop(reduced mass-loss) by
! Mdot_NotCorrected/xmdot.
    write(*,*) 'Mupltiplying factor: ', Mdot_NotCorrected/xmdot
    dLmag = -dL_isotrop*Mdot_NotCorrected/xmdot*(Constant+(etastar+0.25d0)**(1.d0/4.d0))**2.d0
  else
    dLmag = -dL_isotrop
  endif

  return

end subroutine dLmagcalc
!======================================================================
!> Computation of the surface magnetic braking for low-mass stars
!!
!! @brief Computes the surface magnetic braking for low-mass stars according to
!! Kawaler (1988).
!!
!! @param[out] dLdtmagLM
!----------------------------------------------------------------------
double precision function dLmagLM()
!----------------------------------------------------------------------
  use const,only: pi,cst_sigma,Lsol,Rsol,Omega_sol
  use inputparam,only: K_Kawaler,Omega_saturation
  use caramodele,only: gms,gls,teff
  use timestep,only: dzeit
  use rotmod, only: omegi

  implicit none

  real(kindreal):: rstar_Rsol,common_factor
!----------------------------------------------------------------------
   if (K_Kawaler /= 0.d0) then
     rstar_Rsol = sqrt(gls*Lsol/(4.d0*pi*cst_sigma))/(teff**2.d0)/Rsol
     common_factor = -omegi(1)*K_Kawaler*sqrt(rstar_Rsol/gms)
     write(*,*) 'rstar = ',rstar_Rsol
     write(*,*) 'mass = ', gms
     write(*,*) 'omega_surf = ', omegi(1)
     write(*,*) 'K, omega_sat = ', K_Kawaler, Omega_sol*Omega_saturation
     if (omegi(1) <= Omega_sol*Omega_saturation) then
       dLmagLM = omegi(1)**2.d0*common_factor
     else
       dLmagLM = (Omega_sol*Omega_saturation)**2.d0*common_factor
     endif
   else
     dLmagLM = 0.d0
   endif
   dLmagLM = dLmagLM*dzeit
  return

end function dLmagLM
!=======================================================================
!> Computation of the anisotropy of the winds
subroutine aniso(f,yyygmo,rrro)
! f=oblateness Req/Rpo
!----------------------------------------------------------------------
  use const,only: cst_G,Msol,lgLsol,cstlg_sigma,Lsol,pi,lgpi
  use inputparam,only: ianiso,isol
  use caramodele,only: teff,gls,gms,dm_lost,eddesc,rayequat
  use rotmod,only: omegi,bdotis,xlexcs
  use timestep,only: dzeit

  implicit none

  real(kindreal),intent(in)::rrro,f,yyygmo

  real(kindreal), dimension(1001):: rayeq,sint,cost,theta_an
  real(kindreal), dimension(1000):: dtheta_an,thetam,dmai,dsurf,disax
  real(kindreal):: xlogtc,xl,xo,xuni,gguni,sguni,rguni,vpsi,dmis,surf,xmoman,xmdani,xmstar,xlpgmp,edding,xequa,xpole,pota,dx, &
    sint2,rmean,themea,themea1,themea2,ger,gthe,gnor,rapg,tfmean,xogtef,allam,xkkk,xa,xb,xc,cosep,dl,aaaa,bbbb,cccc,xccmdo, &
    bdotan,xrad

  integer:: i
!----------------------------------------------------------------------
! constantes
  xlogtc=log10(teff)
  xl=log10(gls)
  xo=omegi(1)

  if (f > 1.00000001d0 .and. isol == 0) then
    xuni=cst_G*gms*Msol/(xo*xo)
    xuni=xuni**(1.d0/3.d0)
    gguni=(cst_G*gms*Msol*xo*xo*xo*xo)**(1.d0/3.d0)
    sguni=(cst_G*gms*Msol/(xo*xo))**(2.d0/3.d0)
    rguni=(cst_G*gms*Msol/(xo*xo))**(1.d0/3.d0)
    vpsi=10.d0**(xl+lgLsol-4.d0*xlogtc-cstlg_sigma)
    dmis=dm_lost/vpsi
    surf=0.d0
    xmoman=0.d0
    xmdani=0.d0
    xmstar=gms*Msol*(1.d0-rrro)
    xlpgmp=gls*Lsol/(4.d0*pi*cst_G*xmstar)
    edding=eddesc/(1.d0-rrro)
!******************************************************
! (1) calcul de la forme de la surface r theta
!******************************************************
    xequa=(2.d0*(f-1.d0))**(1.d0/3.d0)
    xpole=xequa/f
    pota=1.d0/xpole
! pour calculer la forme d'une equipotentielle, on decoupe
! l'intervalle [xpole;xequa] en 1000 points
    dx=(xequa-xpole)/1000.d0
! on initialise les rayons et les angles sur l'equipotentielle.
! Par angle on entend ici sin theta et cos theta
!***********************************
    rayeq(1)=xpole
    sint(1)=0.0d0
    cost(1)=1.0d0
    theta_an(1)=0.0d0
!***********************************
    do i=2,1000
!***********************************
     rayeq(i)=rayeq(i-1)+dx
     sint2 =(pota-1.d0/rayeq(i))*2.d0/(rayeq(i)*rayeq(i))
     sint(i)=sqrt(sint2)
     cost(i)=sqrt(1.d0-sint2)
     theta_an(i)=asin(sint(i))
     dtheta_an(i-1)=theta_an(i)-theta_an(i-1)
!***********************************
     if(sint2<0.d0.or.sint2>1.d0) then
       stop 'problem with sint2 in aniso'
     endif
    enddo
!***********************************
    rayeq(1001)=xequa
    sint(1001)=1.0d0
    cost(1001)=0.0d0
    theta_an(1001)=pi/2.d0
    dtheta_an(1000)=theta_an(1001)-theta_an(1000)
!***************************************************************
! (2) calcul de la gravite effective en fonction de theta
!***************************************************************
    do i=1,1000
     rmean=(rayeq(i)+rayeq(i+1))/2.d0
     themea=(pota-1.d0/rmean)*2.d0/(rmean*rmean)
     themea1=sqrt(themea)
     themea2=sqrt(1.d0-themea)
     ger=-1.d0/(rmean*rmean)+themea*rmean
     gthe=rmean*themea1*themea2
     gnor=sqrt(ger*ger+gthe*gthe)
     thetam(i)=asin(themea1)
!***************************************************************
! (3) calcul de la Teff effective en fonction de theta
!***************************************************************
     rapg=(gnor/yyygmo)**0.25d0
     tfmean=rapg*teff
!***************************************************************
! (4) calcul du alpha(Teff) en fonction de theta
!***************************************************************
     xogtef=log10(tfmean)

! choix du alpha et du k
! Multiplicateurs de force pour la perte de masse due ‡ la rotation
! allam: Communication privÈe de Lamers (2004).
     if (xogtef >= 4.3d0) then
       allam = 0.6d0
     else if (xogtef >= 4.05d0 .and. xogtef < 4.3d0) then
       allam = 0.43d0
     else
       allam = 0.33d0
     endif
! xkkk: Paper IV (1999).
     if (xogtef >= 4.6d0) then
       xkkk = 0.124d0
     else if (xogtef >= 4.48d0 .and. xogtef < 4.6d0) then
       xkkk = 0.17d0
     else
       xkkk = 0.32d0
     endif
!***************************************************************
! (5) calcul de la surface Spsi
!***************************************************************
     xa=1.d0/(rmean*rmean)-rmean*themea
     xb=xa*xa
     xc=(rmean*themea1*themea2)*(rmean*themea1*themea2)
     cosep=xa/sqrt(xb+xc)
     dl=rmean*dtheta_an(i)/cosep
     dsurf(i)=2.d0*pi*rmean*themea1*dl*sguni
     disax(i)=rmean*rmean*themea1*themea1
     disax(i)=rguni*rguni*disax(i)
     surf=surf+2.d0*pi*rmean*themea1*dl*sguni
!***************************************************************
! (6) calcul de la
!     perte de masse anisotrope par unite de surface dmai
!     en masses solaires par annee
!***************************************************************
     aaaa=(xkkk*allam)**(1.d0/allam)
     aaaa=aaaa*((1.d0-allam)/allam)**((1.d0-allam)/allam)
     bbbb=xlpgmp**(1.d0/allam - 1.d0/8.d0)
     cccc=(1.d0-edding)**((1.d0-allam)/allam)
     dmai(i)=aaaa*bbbb*(gnor*gguni)**(7.d0/8.d0)/cccc*dzeit/Msol
! calcul de la constante multipliant dmai
     xmdani=xmdani+dmai(i)*dsurf(i)
    enddo
    dmis=dm_lost/(2.d0*surf)
    xmdani=2.d0*xmdani
    xccmdo=dm_lost/xmdani
    bdotan=0.d0
    bdotis=0.d0
    do i=1,1000
     dmai(i)=dmai(i)*xccmdo
     bdotan=bdotan+dmai(i)*disax(i)*dsurf(i)*Msol
     bdotis=bdotis+dmis   *disax(i)*dsurf(i)*Msol
    enddo

    if (bdotis /= 0.0d0) then
      xlexcs=(bdotan/bdotis-1.d0)
    else
      xlexcs= 0.0d0
    endif

    bdotis = 2.d0*bdotis*xo
    if (ianiso == 0 .or. dm_lost == 0.0d0) then
       xlexcs = 0.d0
       bdotan = bdotis
    endif

! calcul du vrai rayon equatorial
    rayequat = xequa*xuni
    write(3,'(a,es14.7)') 'Equatorial radius after aniso:',rayequat

  else if ((f >= 1.d0 .and. f <= 1.00000001d0) .or. isol >= 1) then
     xlexcs = 0.d0
     xrad = (xl + lgLsol - log10(4.d0) - lgpi - cstlg_sigma -4.d0*xlogtc)/2.d0
     xrad = 10.d0**xrad
     bdotis = 2.d0*xo*dm_lost*Msol*xrad**2.d0/3.d0
     rayequat=xrad
  else
     stop 'f smaller than 1, aborting...'
  endif

  return

end subroutine aniso
!=======================================================================
subroutine corrwind(teffpr,teffel,xmdot,teff,raysl)
!----------------------------------------------------------------------
! Sous-routine calculant la correction a la temperature effective due
! au vent stellaire tenant compte de la diffusion par electrons libres
! et l'influence des raies. (Theorie CAK amelioree)
! (D. Schaerer, Fin juillet 1990)
!-----------------------------------------------------------------------
!  Parametres comme dans MAIN.
!  Parametre change: TEFFPR.
  use const,only: Msol,year,Rsol,cst_thomson,cst_mh,cst_k,pi
  use strucmod,only: patmos,tatmos,vmion,vmol,g
  use ionisation,only: ionpart

  implicit none

  real(kindreal),intent(in):: xmdot,teff,raysl
  real(kindreal),intent(out):: teffpr,teffel

  integer::i
  real(kindreal):: xMsolyear,vinf,bet,ck,alp,del,r0,xmax,step,reff,xmd,se,tintlines,tintelec,tautot,tmcak,xold, &
                   tauelecold,tauold,ex,vth,xmin,xlen,x_corr,gr,qu,xh,xtwo,cf,rh,xeff
!-----------------------------------------------------------------------
! Les variables PATMOS et TATMOS representent le log de la pression
! resp. temperature a tau=2/3 de l'atmosphere du programme d'evolution.
! On s'en sert pour calculer la concentration d'electrons. De cette
! maniere on peut tenir compte du changement pour les etoiles WR.
!----  Constantes:  ----------------------------------------------------
! xMSOLYEAR: Perte de masse en unites 'Masse solaire par an'
  xMsolyear=Msol/year
!***** Debut de la sous-routine: ***************************************
! Parametres du vent et des raies:
  vinf = 2000.d0
  bet =  2.0d0
  ck =   0.124d0
  alp =  0.64d0
  del =  0.07d0
  r0 =   1.0d0
  xmax = 1000.d0
  step = 600.d0

!-----Conversions en bonnes unites:----------------------------------
  reff=raysl*Rsol
  xmd=(10.d0**(xmdot))*xMsolyear
  vinf=vinf*1.d+5
  r0=r0*reff
!----------Initialisations:-------------------------------------------
! Langer posait SE=0.22 dans la correction appliquee jusqu' a present.
! Ici cependant on recalcule SE:
  call ionpart(patmos,tatmos)
! Lamers & Cassinelli, page 223, eq. 8.93 + def poids moleculaire moyen electrons.
  se=cst_thomson/cst_mh*(1.d0/vmion-1.d0/vmol)

! On met TEFFPR,TEFFEL = 0 au debut. Si la routine donne ces valeurs
! ceci signalerait que tau<2/3 jusqu'a X=1.0001 * REFF.
! Dans ce cas les erreurs deviendrait trop grandes avec ce champ de
! vitesse adopte.
  teffpr=0.d0
  teffel=0.d0

  tintlines=0.d0
  tintelec=0.d0
  tautot=0.d0
  tmcak=0.d0
  xold=0.d0
  tauelecold=0.d0
  tauold=0.d0

  ex=(2.d0*alp*bet-alp-bet+1.d0)
! cf. Lamers & Cassinelli, page 217, eq. 8.83
  vth=dsqrt(2.d0*cst_k * teff/cst_mh)
!------------------Methode MCAK avec integration:-------------------
!  Avec le choix de XLEN=XMAX+4. on s'approchera au maximum jusqu'a
!  1.0001 * REFF.
  xmin=-4.d0
  xmax=log10(xmax)
  xlen=xmax+abs(xmin)

  do i=nint(abs(step)),0,-1
   x_corr=reff*(1.d0+10.d0**( i*xlen/abs(step) + xmin))

!  Densite, Gradient dans le vent
!      (differentes versions pour BETA >1, <1)
   if (bet >= 1.d0) then
     gr=vinf*bet*((1.d0-r0/x_corr)**(bet-1.d0))*r0/(x_corr*x_corr)
   else
     gr=(vinf*bet*r0)/((1.d0-r0/x_corr)**(1.d0-bet)*x_corr*x_corr)
   endif
   tmcak=se*vth*xmd/(4.d0*pi*x_corr*x_corr*vinf*(1.d0-r0/x_corr)**bet*gr)

! Partie dependante de delta:
   if (bet >= 2.d0) then
     qu=22.5d0**del-1.d0
   elseif (bet >= 1.d0) then
     qu=7.5d0**del-1.d0
   elseif (bet >= 0.7d0) then
     qu=4.d0**del-1.d0
   elseif (bet >= 0.5d0) then
     qu=2.5d0**del-1.d0
   else
     qu=1.18d0**del-1.d0
   endif
   g=((xmd/(4.d0*pi*reff*reff*vinf))**del)*((1.02436d+13)**del)*2.d0**del*(qu*(reff/x_corr)**2.d0+1.d0)
! Partie du Correction factor:
   xh=((x_corr/reff)-1.d0)/bet
   xtwo=(x_corr/reff)*(x_corr/reff)
   cf=(xtwo/((alp+1.d0)*(1.d0-xh)))*(1.d0-(1.d0+(xh-1.d0)/xtwo)**(alp+1.d0))

   tmcak=(ck/(tmcak**alp)) * g * cf
   rh=xmd/(4.d0*pi*x_corr*x_corr*vinf*(1.d0-r0/x_corr)**bet)

   if (xold /= 0.d0) then
     tintlines=tintlines+(xold-x_corr)*rh*se*tmcak
! Si l'on utilise un champ de vitesse different il faudra aussi integrer
! la prof.opt. des electrons numeriquement (prochaine ligne).
! Ici on utilise la version analytique...
!          TINTELEC=TINTELEC+(XOLD-X)*RH*SE

!  Definition des 2 profondeurs optiques (electrons, raies)
!        (differentes versions pour BETA >1, <1)
     if (bet >= 1.d0) then
       tintelec=se*xmd/(4.d0*pi*vinf*r0*(1.d0-bet))*(1.d0-1.d0/(1.d0-r0/x_corr)**(bet-1.d0))
     else
       tintelec=se*xmd/(4.d0*pi*vinf*r0*(1.d0-bet))*(1.d0-(1.d0-r0/x_corr)**(1.d0-bet))
     endif

     tautot=tintlines+tintelec
   endif

! Quand l'integration a depassee tau=2/3 on trouve le XEFF par interpola
! lineaire.
   if (tautot >= (2.d0/3.d0).and.teffpr == 0.d0) then
     xeff=xold-(xold-x_corr)*((2.d0/3.d0)-tauold)/(tautot-tauold)
     teffpr=log10(teff) + 0.5d0*log10(reff/xeff)
   endif

! Comme comparaison on fait la meme chose pour la profondeur optique des
! electrons (sans raies) et puis on peut TERMINER l'integration.
   if (tintelec >= (2.d0/3.d0)) then
     xeff=xold-(xold-x_corr)*((2.d0/3.d0)-tauelecold)/(tintelec-tauelecold)
     teffel=log10(teff) + 0.5d0*log10(reff/xeff)

     return
   endif

   tauelecold=tintelec
   tauold=tautot
   xold=x_corr

  enddo
!----------------fin MCAK-------------------------------------------
  return

end subroutine corrwind
!=======================================================================
real(8) function Compute_MagCorrection(rstar,mdot,Bsurf,v_inf,r_K)
!----------------------------------------------------------------------
! The purpose of this function is to compute the reduction factor to be applied to
! the mass-loss rate due to external magnetic fields, according to Petit et al.
! MNRAS (arXiv 1611.08964).
!----------------------------------------------------------------------

implicit none

real(kind=kindreal), intent(in)::rstar,mdot,Bsurf,v_inf,r_K
real(kind=kindreal):: r_a,r_c,mdot_factor,Constant

Constant = 1.d0 - 0.25d0**(1.d0/4.d0)
r_a = (Constant + ((Bsurf*rstar)**2.d0/(mdot*v_inf) + 0.25d0)**0.25d0)*rstar
write(*,*) 'ETA CMC: ', (Bsurf*rstar)**2.d0/(mdot*v_inf)
r_c = rstar + 0.7d0*(r_a-rstar)
! Here we compare with the Keplerian co-rotation radius. In case it is smaller than r_c,
! it is accounted for instead in the relation.
if (r_K < r_c) then
  r_c = r_K
  write(*,*) 'KEPLERIAN CORROTATION RADIUS ACCOUNTED FOR.'
endif

mdot_factor = 1.d0 - sqrt(1.d0-rstar/r_c)

Compute_MagCorrection = mdot_factor

return

end function Compute_MagCorrection
!=======================================================================

end module winds
