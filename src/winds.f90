!======================================================================
!> Module computing all aspects of mass loss
!!
!!  @version 211
!!  @date 14.2.2013
!======================================================================
module winds

  use io_definitions
  use evol,only: kindreal, input_dir
  use inputparam,only: verbose,Xs_WR,RSG_Mdot,WR_Mdot,OB_Mdot,Fallback_Mdot,zsol,&
                       zinit,fmlos,hardJump,imloss,force_prescription
  use caramodele,only: nwmd,xmini,gms,teff,teffv,gls,glsv,eddesc,is_MS,is_OB,is_RSG,is_WR

  implicit none

  integer,save:: imloss_fallback,imloss_ob,imloss_wr,imloss_rsg
  real(kindreal),save:: zheavy,xlogz,alpha_winds
  logical:: WRNoJump,checkVink

  integer,save:: lenf
  integer,dimension(:),allocatable,save:: imdot
  character(256),dimension(:),allocatable,save:: bmdot,smdot


  private
  public :: aniso,corrwind,xloss,xldote
  public:: read_Mdot_prescriptions,print_Mdot_prescription

contains
!======================================================================
subroutine xloss
!> Computation of the radiative mass loss
! ------------------------------
! Possible values for OB_MDOT
! ------------------------------
!  0: none
!  1: de Jager+ 1988
!  2: mass loss in Msol/yr given by FMLOS
!  3: Vink+ 2001
!  4: Vink+ 2001 modified by Markova & Puls 2008 + priv. comm. Puls (nov. 2010)
!  5: Kudritzki & Puls 2000
!  6: Kudritzki 2002
!  7: Bestenlehner+ 2020
!  8: Bjorklund+ 2022
!  9: Gormaz-Matamala+ 2022
! ------------------------------
! Possible values for RSG_MDOT
! ------------------------------
!  0: none
!  1: Sylvester (1998), van Loon 1999 (cf Crowther 2001)
!  2: de Jager+ 1988
!  3: Beasor & Davies 2020
!  4: Kee+ 2021
!  5: van Loon+ 2005
!  6: Reimers formula with etaR given by FMLOS
! ------------------------------
! Possible values for WR_MDOT
! ------------------------------
!  0: none
!  1: Graefener & Hammann (2008)
!  2: Nugis & Lamers (2000)
!  3: Schmutz (1997) except for WNL = Nugis+ (1998)
! ------------------------------
! Possible values for FALLBACK_MDOT
! ------------------------------
!  0: none
!  1: de Jager+ 1988
!  2: mass loss in Msol/yr given by FMLOS
!----------------------------------------------------------------------
  use const, only: lgLsol,lgpi,cstlg_sigma,lgRsol,cst_thomson,cst_avo,xlsomo,qapicg,cst_G,Msol,Rsol, &
                   Lsol,cst_sigma,pi,year,cst_k,cst_c,cst_mh,Teffsol
  use inputparam, only: phase,ipop3,irot,B_initial,frein
  use caramodele, only: xmdot,zams_radius,Mdot_NotCorrected
  use abundmod, only: x,y,y3
  use rotmod, only: fMdot_rot,omegi

  implicit none

  integer:: imlossprev
  real(kindreal):: xteff,ygls,zlim,rstar, &
                   Bsurf,v_inf,v_esc,r_K,Correction_factor, &
                   xmdotwr,xmdotrsg,xmdotob,xmdotfallback
  character(256):: mdotpresc
!----------------------------------------------------------------------
  rstar = 0.d0
  Bsurf = 0.d0
  v_inf = 0.d0
  v_esc = 0.d0
  imlossprev = imloss
  WRNoJump = .false.
  xmdotrsg = 0.d0
  xmdotwr = 0.d0
  xmdotob = 0.d0
  xmdotfallback = 0.d0
  xmdot = 0.d0
  imloss_ob = 0
  imloss_rsg = 0
  imloss_wr = 0
  imloss_fallback = 0

  xteff=log10(teff)
  ygls=log10(gls)
  write(io_logs,*) 'xteff= ',xteff

  call Star_type

  ! calcul de la dependance en metallicite
  ! le Log facteur=xlgfz
  ! ce facteur est utilise partout dans le DHR
  ! sauf si l'etoile est WR
  zheavy=1.d0-x(1)-y(1)-y3(1)
  zlim=1.d-04*zsol
  if (zinit <= zlim) then
    zheavy = min(zheavy,zlim)
  endif

  if (ipop3 == 1) then
    xlogz=log10(zheavy/zsol)
  else
    if (is_WR  > epsilon(is_WR)) then
      xlogz=log10(zinit/zsol)
    else
      xlogz=log10(zheavy/zsol)
    endif
  endif

  fMdot_rot = fMdot_rot_calc()
  write(io_logs,*) 'fMdot_rot = ',fMdot_rot,'eddesc = ',eddesc

!-----------------------------------------------------------------------------
! computation of the various Mdot
  xmdotfallback = Fallback_Mdot_calc()
  xmdotrsg = RSG_Mdot_calc()
  xmdotwr = WR_Mdot_calc()
  xmdotob = OB_Mdot_calc(xmdotfallback,imloss_fallback)
  write(io_logs,*) 'imloss_fallback:',imloss_fallback,', Mdot_fallback:',xmdotfallback
  write(io_logs,*) 'imloss_ob:',imloss_ob,', Mdot_OB:',xmdotob
  write(io_logs,*) 'imloss_rsg:',imloss_rsg,', Mdot_RSG:',xmdotrsg
  write(io_logs,*) 'imloss_wr:',imloss_wr,'Mdot_WR:',xmdotwr
!-----------------------------------------------------------------------------
  if (hardJump) then
    if (is_RSG > epsilon(is_RSG)) then
      xmdot = xmdotrsg
      imloss = imloss_rsg
      if (RSG_Mdot == 6) then
        ! fmlos is no more the eta_Reimers, which is now included in the recipe,
        ! so we have to make sure that fmlos=1.
        fmlos = 1.0d0
      endif
    elseif (is_WR > epsilon(is_WR)) then
      if (xmdotob > xmdotwr) then
        write(io_logs,*) 'Mdot OB larger than Mdot WR'
        xmdot = xmdotob
        imloss = imloss_ob
      else
        xmdot = xmdotwr
        imloss = imloss_wr
        if (WR_Mdot == 1 .and. imloss_wr == 202) then
          write(*,*) 'GRAEF transferred to Nugis'
          write(io_logs,*) 'GRAEF transferred to Nugis'
        endif
      endif
    elseif (is_OB > epsilon(is_OB)) then
      xmdot = xmdotob
      imloss = imloss_ob
    else
      ! Vink's and KP prescriptions valid a bit beyond OB type
      if ( (OB_Mdot == 3 .or. OB_Mdot == 4) .and. xteff >= 3.90d0 ) then
        xmdot = xmdotob
        imloss = imloss_ob
      elseif (OB_Mdot == 5 .and. xteff >= 3.95d0) then
        xmdot = xmdotob
        imloss = imloss_ob
      else
        xmdot = xmdotfallback
        imloss = imloss_fallback
      endif
    endif
  else
    write(*,*) 'Soft jumps between mass-loss rates not yet implemented'
    stop
  endif
  write(io_logs,*) 'imloss applied: ',imloss
  if (imloss /= imlossprev) then
    write(io_input_changes,'(i0,a,i0,a,i0)') nwmd,': imloss changed from ',imlossprev,' to ',imloss
  endif
!=======================================================================
  mdotpresc = print_Mdot_prescription(imloss)
  write(*,'(a,i0,a3,a)') 'imloss: ',imloss,' = ',trim(mdotpresc)
  xmdot=fmlos*xmdot
  write(io_logs,*) 'fmlos= ',fmlos,'  xmdot*fmlos= ',xmdot
  xmdot = fMdot_rot*xmdot

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
     write(io_logs,*) 'Mdot (no corrections) = ', xmdot
     Mdot_NotCorrected = xmdot
  endif
  xmdot = xmdot*Correction_factor
  write(io_logs,*) 'fmlos= ',fmlos,'  xmdot= ',xmdot, 'Magnetic correction: ', Correction_factor

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
  use inputparam, only: ianiso,modanf,rapcrilim,bintide,xcn,diff_only,Be_mdotfrac,start_mdot
  use caramodele, only: gms,rayequat,xltotbeg,firstmods,inum
  use strucmod, only: q,r
  use rotmod, only: xlexcs,rapom2,omegi,dlelex,xldoex,bdotis,vsuminenv
  use timestep,only: dzeit
  use bintidemod,only: dLtidcalc

  implicit none

  real(kindreal), intent(out)::dmneed
  real(kindreal), intent(in):: dmdot

  integer:: jo
  real(kindreal):: DeltaMCG,omcrit,dLmag,dL_Kawaler,dLtid

  real(kindreal):: omlim, newomega, dLisotrop,dLmeca, xLe, qcorr, qnum, qdenom, &
                   dmneednum, dmneeddenom,rapcrilim_calc,Be_mdot_factor
  real(kindreal), dimension(npondcouche)::Li
!> facteur multiplicatif: sqrt(r_out/R_star) cf. Owocki (perte L par disque)
  real(kindreal), parameter:: factordisk = 1.d0, omega_min = 1.d-20
!----------------------------------------------------------------------

  dmneed= 0.0d0
  dLmeca = 0.0d0
  dLtid = 0.0d0
  rapcrilim_calc = rapcrilim
  Be_mdot_factor = 1.0d0

  if (Be_mdotfrac > 0.d0 .and. rapom2 >= 0.99*rapcrilim) then
    Be_mdotfrac = 0.d0
  endif

  ! Be_mdotfrac allows for a progressive mechanical mass loss from O/Oc=start_mdot to rapcrilim
  ! At O/Oc=start_mdot, only Be_mdotfrac of dmneed is applied, at rapcrilim the full correction
  ! is applied, and in between, a linear progression is used
  if (Be_mdotfrac > 0.d0 .and. rapom2 >= start_mdot) then
    rapcrilim_calc = rapom2 * 0.99d0
    Be_mdot_factor = Be_mdotfrac + ((1.d0 - Be_mdotfrac)*((rapom2 - start_mdot)/(rapcrilim - start_mdot))**64.d0)
    write(*,*) 'Be_Mdot_factor: ',Be_mdot_factor
    write(io_logs,*) 'Be_Mdot_factor: ',Be_mdot_factor
  endif

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
      write(io_logs,*) ' !!!!! WARNING !!!!!'
      write(io_logs,*) 'rapom2 = 0.0 in xldote. If this model is',' not the 1st, this is an error!'
    endif
  else
! Si rapcrilim est nul, la correction n'est pas appliquee et son calcul
! est inutile et source de bug.
    if (rapcrilim_calc > 1.d-5) then
      omlim= rapcrilim_calc*omegi(1)/rapom2
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
    write(io_logs,*)
    write(io_logs,*) 'In xldote: newomega, omega limit: ',newomega,omlim
    if (newomega >= omlim) then
! Le detail de cette formule se trouve dans la documentation.
      dmneednum = Li(1) * (omlim/omegi(1)-1.d0)
      dmneednum = dmneednum + xLe*(omlim/omegi(1)*(1.d0 + dmdot/(gms*exp(q(1)))) - 1.d0) + dlelex
      dmneeddenom = -factordisk*omegi(1)*rayequat**2.d0 + xLe*omlim / (gms*exp(q(1))*Msol*omegi(1))
      dmneed = Be_mdot_factor * dmneednum / dmneeddenom
      if (omegi(1)  <  1.d-25) then
        dmneed = 0.d0
      endif
    endif

! dmneed doit etre positif. Si ce n'est pas le cas, message d'erreur et arret de l'execution.
    if (dmneed < 0.d0) then
      rewind(io_runfile)
      write(io_runfile,*) nwmd,': dmneed negative in xldote'
      stop 'dmneed negative in xldote. Aborting...'
    endif

! Si la vitesse de surface est superieure de 0.25% a la valeur maximale toleree
! (ou 1 le cas echeant), on multiplie par 1.5 la perte de masse equatoriale.
    write(io_logs,*) 'dmneed = ', dmneed

    if (dmneed > 0.d0) then
      if (rapom2  <=  0.995d0 .and. rapcrilim_calc  >  0.d0) then
        if (rapom2  >  (rapcrilim + 0.0025d0)) then
          dmneed = 2.0d0*dmneed
          write(*,*) '!!! WARNING: equatorial mass loss increased by a factor 2.0!'
          write(io_logs,*) 'dmneed multiplied by 2.0. New dmneed = ',dmneed
        else if (rapom2  >  (rapcrilim + 0.005d0)) then
          dmneed = 4.d0*dmneed
          write(*,*) '!!! WARNING: equatorial mass loss increased by a factor 4.0!'
          write(io_logs,*) 'dmneed multiplied by 4. New dmneed = ',dmneed
        endif
      else
        if (rapom2  >  1.05d0 .and. rapcrilim_calc  >  0.d0) then
          dmneed = 10.d0*dmneed
          write(*,*) '!!! WARNING: equatorial mass loss increased by a factor 10!'
          write(io_logs,*) 'dmneed multiplied by 10. New dmneed = ',dmneed
        else if (rapom2  >  1.d0 .and. rapcrilim_calc  >  0.d0) then
          dmneed = 1.5d0*dmneed
          write(*,*) '!!! WARNING: equatorial mass loss increased by a factor 1.5!'
          write(io_logs,*) 'dmneed multiplied by 1.5. New dmneed = ',dmneed
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

! dlelex contains all changes in angular momentum on a time-step
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
  write(io_logs,*) 'dLiso, xlexcs, dLmeca,  dLaniso, dL_Kawaler, dLtide: '
  write(io_logs,*) dLisotrop,xlexcs,dLmeca,dlelex,dL_Kawaler,dLtid

! dmneed set in solar units, and >0 as a mass loss:
  if (dmneed /= 0.d0) then
     dmneed = -dmneed / Msol
  endif

! computation of Ldot_excess, used as a border condition in advection (badv).
! Beware the signs of dLisotrop and dlelex, that are inverted.
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
    write(io_logs,*) 'Mass loss used in dLmagcalc = ', Mdot_NotCorrected
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
! f = oblateness = Req/Rpo = 1/xobla
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
    sint2,rmean,themea,themea1,themea2,ger,gthe,gnor,rapg,tfmean,xogtef,xkkk,xa,xb,xc,cosep,dl,aaaa,bbbb,cccc,xccmdo, &
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
! Force multipliers used for rotational enhancement of the mass loss
! alpha_winds: private communication from Lamers (2004), already computed in fMdot_rot_calc called in xloss.
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
     aaaa=(xkkk*alpha_winds)**(1.d0/alpha_winds)
     aaaa=aaaa*((1.d0-alpha_winds)/alpha_winds)**((1.d0-alpha_winds)/alpha_winds)
     bbbb=xlpgmp**(1.d0/alpha_winds - 1.d0/8.d0)
     cccc=(1.d0-edding)**((1.d0-alpha_winds)/alpha_winds)
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
    write(io_logs,'(a,es14.7)') 'Equatorial radius after aniso:',rayequat

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
double precision function Compute_MagCorrection(rstar,mdot,Bsurf,v_inf,r_K)
!----------------------------------------------------------------------
! The purpose of this function is to compute the reduction factor to be applied to
! the mass-loss rate due to external magnetic fields, according to Petit et al.
! MNRAS (arXiv 1611.08964).
!----------------------------------------------------------------------

  implicit none

  real(kind=kindreal), intent(in)::rstar,mdot,Bsurf,v_inf,r_K
  real(kind=kindreal):: r_a,r_c,mdot_factor,Constant
!----------------------------------------------------------------------
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
subroutine Star_type
! is_MS stores the percentage of the MS
! is_RSG is 0 if Teff>3.7
!       and 1 if Teff < 3.66
!       with a linear evolution in-between
! is_WR is 0 if Teff<4.0 or Xsurf>=Xs_WR+delta_Xs_WR
!      and 1 if Teff>4.0 and Xsurf<X=s_WR
!      with a linear evolution wrt Xs in-between (only if Teff>4.0)
! is_OB is 1-is_WR
  use caramodele,only: xini
  use inputparam,only: phase
  use strucmod, only: m
  use abundmod, only: x

  implicit none

  real(kindreal),parameter:: delta_Xs_WR = 0.1d0
!----------------------------------------------------------------------
  is_MS = 0.d0
  is_OB = 0.d0
  is_WR = 0.d0
  is_RSG = 0.d0

  is_MS = x(m)/xini

  if (xmini >= 8.5d0) then
    if (log10(teff) <= 3.7d0) then
      is_RSG = (3.7d0-log10(teff))/(3.7d0-3.66d0)
    endif
    if (x(1) < (Xs_WR+delta_Xs_WR) .and. log10(teff) > 4.d0) then
      if (x(1) > Xs_WR) then
        is_WR = 1.d0-(x(1)-Xs_WR)/delta_Xs_WR
      else
        is_WR = 1.d0
      endif
    endif
    if (log10(teff) > 3.9d0) then
      is_OB = 1.d0 - is_WR
    endif
  else
    if (phase /= 1 .and. log10(teff) < 3.8d0 .and. gls > glsv) then
      is_RSG = 1.d0
    endif
    if (log10(teff) > 3.9d0) then
      is_OB = 1.d0 - is_RSG
    endif
  endif
  write(io_logs,*) 'STAR TYPE:'
  write(io_logs,*) '---------------'
  write(io_logs,*) 'is_MS =',is_MS
  write(io_logs,*) 'is_OB =',is_OB
  write(io_logs,*) 'is_RSG=',is_RSG
  write(io_logs,*) 'is_WR =',is_WR
  write(io_logs,*) '---------------'

end subroutine Star_type
!=======================================================================
double precision function RSG_Mdot_calc()
  implicit none

  real(kindreal):: mdot
!----------------------------------------------------------------------
  select case (RSG_Mdot)
    case (0)
      mdot = 0d0
      imloss_rsg = 0
    case (1)
      mdot = deJager88()
      imloss_rsg = 1
    case (2)
      mdot = Crowther01()
      imloss_rsg = 302
    case (3)
      mdot = Crowther01()
      imloss_rsg = 303
    case (4)
      mdot = Beasor20()
      imloss_rsg = 304
    case (5)
      mdot = Kee21()
      imloss_rsg = 305
    case (6)
      mdot = Reimers75()
      imloss_rsg = 306
    case (7)
      mdot = vanLoon05()
      imloss_rsg = 307
    case (8)
      mdot = deJager88_lin()
      imloss_rsg = 308
    case (9)
      mdot = Nieuwenhuijzen90()
      imloss_rsg = 309
    case (10)
      mdot = Vanbeveren98()  ! /!\ Only valid for Teff < 10kK - [MM]
      imloss_rsg = 310
    case (11)
      mdot = Salasnich99()
      imloss_rsg = 311
    case (12)
      mdot = Decin21()
      imloss_rsg = 312
    case (13)
      mdot = Decin24()
      imloss_rsg = 313
    case (14)
      mdot = Yang23()
      imloss_rsg = 314
    case (15)
      mdot = Wachter02()
      imloss_rsg = 315
    case (16)
      mdot = Schroder05()
      imloss_rsg = 316
    case (17)
      mdot = Vink23()
      imloss_rsg = 317
    case default
      write(*,*) 'Bad RSG_Mdot value, should be:'
      write(*,*) '    0 none'
      write(*,*) '    1 (de Jager+ 1988)'
      write(*,*) '    2 (mass loss in Msol/yr given by FMLOS)'
      write(*,*) '    3 (Crowther 2001 (standard GENEC))'
      write(*,*) '    4 (Beasor & Davies 2020)'
      write(*,*) '    5 (Kee+ 2021)'
      write(*,*) '    6 (Reimers 1975)'
      write(*,*) '    7 (van Loon+ 2005)'
      write(*,*) '    8 (de Jager 1988 (linear))'
      write(*,*) '    9 (Nieuwanhuijzen 1990)'
      write(*,*) '   10 (Vanbeveren 1998)'
      write(*,*) '   11 (Salasnich 1999)'
      write(*,*) '   12 (Decin 2021)'
      write(*,*) '   13 (Decin 2024)'
      write(*,*) '   14 (Yang 2023)'
      write(*,*) '   15 (Wachter 2002)'
      write(*,*) '   16 (Schroder 2005)'
      write(*,*) '   17 (Vink+ 2023)'
      stop
  end select
  RSG_Mdot_calc = mdot

end function RSG_Mdot_calc
!=======================================================================
double precision function WR_Mdot_calc()
  use inputparam,only: ipop3
  use abundmod, only: x,y,y3,xc12,xo16,xn14

  implicit none

  real(kindreal):: xteff,mdot,ggam0
!----------------------------------------------------------------------
  xteff = log10(teff)

  if( x(1)>0.d0 .and. (OB_Mdot==7 .or. OB_Mdot==12 .or. OB_Mdot==13) ) then
    select case (OB_Mdot)
      case (7)
        mdot = Bestenlehner20()
        imloss_wr = 107
      case (12)
        mdot = Sabhahit22()
        imloss_wr = 112
      case (13)
        mdot = Grafener21(10.d0)
        imloss_wr = 113
      case default
        write(*,*) "Wrong prescription number in WR_Mdot_calc()"
        stop
    end select
      write(io_logs,*) "---> OB mass loss prescription has been used in WR prescription."
      WR_Mdot_calc = mdot
      return
  endif

  select case (WR_Mdot)
    case (0)
      mdot = 0.d0
      imloss_wr = 0
    case (1)
      mdot = deJager88()
      imloss_wr = 1
    case (2)
      mdot = 1.d0
      imloss_wr = 2
    case (3)
! formulae 3, 5, and 6 of Graefener & Hamann (2008A&A...482..945G)
      if (.not. force_prescription) then
        if ((xteff < 4.477d0.or.xteff > 4.845d0) .or.(xlogz < -3.d0.or.xlogz > 0.30d0)) then
          mdot = Nugis00(x(1),y(1),xc12(1),xo16(1))
          imloss_wr = 204
        else
          ggam0=0.326d0-0.301d0*xlogz-0.045d0*xlogz*xlogz
          if (eddesc <= ggam0) then
            mdot = Nugis00(x(1),y(1),xc12(1),xo16(1))
            imloss_wr = 204
          else
            mdot = Graefener08(x(1),y(1),xc12(1),xo16(1))
            imloss_wr = 203
          endif
        endif
      else
        mdot = Graefener08(x(1),y(1),xc12(1),xo16(1))
        imloss_wr = 203
      endif
    case (4)
      mdot = Nugis00(x(1),y(1),xc12(1),xo16(1))
      imloss_wr = 204
    case (5)
      mdot = Schmutz97(x(1),xc12(1),xn14(1))
      imloss_wr = 205
    case (6)
      mdot = Hainich15(y(1), y3(1))
      imloss_wr = 206
    case (7)
      mdot = Langer89(y(1)+y3(1), xc12(1), xo16(1))
      imloss_wr = 207
    case (8)
      mdot = Yoon06(x(1))
      imloss_wr = 208
    case (9)
      mdot = Nugis00_bis(y(1), y3(1))
      imloss_wr = 209
    case (10)
      mdot = Sander20()
      imloss_wr = 210
    case (11)
      mdot = Vink17()
      imloss_wr = 211
    case (12)
      mdot = Shenar19(y(1), y3(1))
      imloss_wr = 212
    case (13)
      mdot = Tramper16(y(1), y3(1))
      imloss_wr = 213
    case (14)
      print*, '!!! Sander23() not implemented yet. mdot = 0.0 !!!'
      stop
      !mdot = Sander23()
      imloss_wr = 214
    case default
      write(*,*) 'Bad WR_Mdot value, should be:'
      write(*,*) '    0 (none)'
      write(*,*) '    1 (de Jager+ 1988)'
      write(*,*) '    2 (mass loss in Msol/yr given by FMLOS)'
      write(*,*) '    3 (Graefener & Hammann 2008)'
      write(*,*) '    4 (Nugis & Lamers 2000)'
      write(*,*) '    5 (Schmutz 1997, except for WNL = Nugis+ 1998)'
      write(*,*) '    6 (Hainich 2015)'
      write(*,*) '    7 (Langer 1989)'
      write(*,*) '    8 (Yoon+ 2006)'
      write(*,*) '    9 (Nugis & Lamers 2000, combined eq. for WN and WC)'
      write(*,*) '   10 (Sander 2020)'
      write(*,*) '   11 (Vink 2017)'
      write(*,*) '   12 (Shenar 2019)'
      write(*,*) '   13 (Tramper 2016)'
  end select
  WR_Mdot_calc = mdot

end function WR_Mdot_calc
!=======================================================================
double precision function OB_Mdot_calc(mdotfallback,imloss_fallback)
  use const,only: cst_G,Msol,Rsol,Teffsol
  use strucmod, only: m
  use abundmod, only: x,y,y3

  implicit none

  integer,intent(in):: imloss_fallback
  real(kindreal),intent(in):: mdotfallback
  real(kindreal):: logg,mdot
!-----------------------------------------------------------------------
  logg = 0.d0

  select case(OB_Mdot)
  case (0)
    mdot = 0.d0
    imloss_ob = 0
  case (1)
    mdot = deJager88()
    imloss_ob = 1
  case (2)
    mdot = 1.d0
    imloss_ob = 2
  case (3)
    if (.not. force_prescription) then
      if (xmini > 15.d0 .and. log10(teff) >= 3.90d0) then
        mdot = Vink01()
        imloss_ob = 103
      else
        mdot = mdotfallback
        imloss_ob = imloss_fallback
      endif
    else
      mdot = Vink01()
      imloss_ob = 103
    endif
  case (4)
    if (.not. force_prescription) then
      if (xmini > 15.d0 .and. log10(teff) >= 3.90d0) then
        mdot = V01MP08()
        imloss_ob = 104
      else
        mdot = mdotfallback
        imloss_ob = imloss_fallback
      endif
    else
        mdot = V01MP08()
        imloss_ob = 104
    endif
  case (5)
    if (.not. force_prescription) then
      if (xmini > 15.d0 .and. log10(teff) >= 3.95d0) then
        mdot = KP2000(x(m),x(1),y(1))
        imloss_ob = 105
      else
        mdot = mdotfallback
        imloss_ob = imloss_fallback
      endif
    else
      mdot = KP2000(x(m),x(1),y(1))
      imloss_ob = 105
    endif
  case (6)
    mdot = Kudritzki02()
    imloss_ob = 106
  case (7)
    if (.not. force_prescription) then
      if (teff>=3.0d4) then
        mdot = Bestenlehner20()
        imloss_ob = 107
      else
        mdot = mdotfallback
        imloss_ob = imloss_fallback
      endif
    else
      mdot = Bestenlehner20()
      imloss_ob = 107
    endif
  case (8)
    if (.not. force_prescription) then
      if (log10(gls)>=4.5d0 .and. log10(gls)<=6.0d0 .and. xmini>=15.0d0 .and. xmini<=80.0d0 &
        .and. teff>=1.5d4 .and. teff<=5.0d4 .and. zinit/zsol>=0.2 .and. zinit/zsol<=1.0) then
        mdot = Bjorklund23()
        imloss_ob = 108
      else
        mdot = mdotfallback
        imloss_ob = imloss_fallback
      endif
    else
      mdot = Bjorklund23()
      imloss_ob = 108
    endif
  case (9)
    if (.not. force_prescription) then
      logg=log10(cst_G)+log10(gms*Msol)-2.d0*log10((sqrt(gls)*(Teffsol/teff)**2.d0)*Rsol)
      if (logg > 3.20d0 .and. log10(teff)>=4.48d0) then
        mdot = Gormaz22(x(1),y(1),y3(1))
        fmlos = 0.850d0
        imloss_ob = 109
      else
        WRNoJump = .true.
        mdot = Vink01()
        fmlos = 0.850d0
        imloss_ob = 103
      endif
    else
      mdot = Gormaz22(x(1),y(1),y3(1))
      fmlos = 0.850d0
      imloss_ob = 109
    endif
  case (10)
      mdot = Krticka21()
      imloss_ob = 110
  case (11)
      mdot = Sabhahit22()
      imloss_ob = 111
  case (12)
      mdot = Grafener21(10.d0)
      imloss_ob = 112
  case default
      write(*,*) 'Bad OB_Mdot value, should be:'
      write(*,*) '    0 (none)'
      write(*,*) '    1 (de Jager+ 1988)'
      write(*,*) '    2 (mass loss in Msol/yr given by FMLOS)'
      write(*,*) '    3 (Vink+ 2001)'
      write(*,*) '    4 {Vink+ 2001 modified by Markova & Puls 2008}'
      write(*,*) '    5 (Kudritzki & Puls 2000)'
      write(*,*) '    6 (Kudritzki 2002)'
      write(*,*) '    7 (Bestenlehner+ 2020)'
      write(*,*) '    8 (Bjorklund+ 2022)'
      write(*,*) '    9 (Gormaz-Matamala+ 2022)'
      write(*,*) '   10 (Krticka+ 2021)'
      write(*,*) '   11 (Sabhahit+ 2022)'
      write(*,*) '   12 (Grafener 2021)'

  end select
  OB_Mdot_calc = mdot

end function OB_Mdot_calc
!=======================================================================
double precision function Fallback_Mdot_calc()
  implicit none

  real(kindreal):: mdot
!-----------------------------------------------------------------------
  select case (Fallback_Mdot)
    case (0)
      mdot = 0.d0
      imloss_fallback = 0
    case (1)
      mdot = deJager88()
      imloss_fallback = 1
    case (2)
      mdot = 1.d0
      imloss_fallback = 2
    case (3)
      mdot = deJager88_lin()
      imloss_fallback = 308
    case default
      write(*,*) 'Bad Fallback_Mdot value, should be:'
      write(*,*) '    0 (none)'
      write(*,*) '    1 (de Jager+ 1988)'
      write(*,*) '    2 (mass loss in Msol/yr given by FMLOS)'
      write(*,*) '    3 (de Jager+ 1988 linear)'
  end select
  Fallback_Mdot_calc = mdot

end function Fallback_Mdot_calc
!=======================================================================
double precision function fMdot_rot_calc()
! Computation of the rotation factor correcting the mass loss fMdot_rot
  use rotmod,only: ivcalc,xobla,vpsi,rapom2,rrro,ygmoye
  use geomod,only: geomat,geomeang,rpsi_min
  implicit none

  real(kindreal) :: xpsi,rap1,xft,rap2,xgmoym,xrequa,ygequa,rapg, &
                    teffeq,xogtef,fffff
!-----------------------------------------------------------------------
  alpha_winds = 0.d0
  fffff = 1.d0/xobla
  if (ivcalc .or. abs(1.d0-xobla)>=1.0d-10) then
! calcul de O^2/(2 pi G rhom)
    call geomat(vpsi,xpsi,rap1,xft,rap2,xgmoym)
    rrro=2.d0/3.d0*xpsi*xpsi*xpsi

! Calcul de la Teff a l equateur
    if (xpsi >= rpsi_min) then
       call geomeang(xpsi,ygmoye)
       xrequa=(2.d0*(fffff-1.d0))**(1.d0/3.d0)
       ygequa=1.d0/(xrequa**2.d0)-xrequa
       rapg  =(ygequa/ygmoye)**0.25d0
       teffeq=teff*rapg
    else
       teffeq=teff
    endif

! Multiplicateurs de force pour la perte de masse due a la rotation
! Communication privee de Lamers (2004).
! choix du alpha
    xogtef=log10(teffeq)
    if (xogtef < 4.05d0) then
      alpha_winds = 0.33d0
    else if (xogtef >= 4.05d0.and.xogtef < 4.3d0) then
      alpha_winds = 0.43d0
    else if (xogtef >= 4.3d0) then
      alpha_winds = 0.6d0
    endif

! FACTEUR CORRECTIF DE LA PERTE DE MASSE
    if (imloss == 2) then
       fMdot_rot_calc=1.0d0
    else
      if (rapom2 < 1.d0) then
        fMdot_rot_calc=((1.d0-eddesc)/(1.d0-rrro-eddesc))**(1.d0/alpha_winds-1.d0)
      else
! Cas d'un modele surcritique. Dans ce cas, on augmente fortement la perte de masse (sensee diverger).
        write(*,'(a)') 'Warning: star overcritical. Mass loss increased by a factor of 100'
        write(io_logs,'(a)') 'Warning: star overcritical. Mass loss increased by a factor of 100'
        fMdot_rot_calc = 100.d0
      endif
    endif
    write(io_logs,*) 'rrro = ',rrro
  else   !< not ivcalc
! Si la rotation n'est pas traitee, on initialise tout de meme les variables utilisees ci-dessus.
! Certaines etant imprimee, le resultat est plus propre.
    rrro=0.d0
    teffeq=teff
    alpha_winds=1.d0
    fMdot_rot_calc=1.d0
  endif   !   ivcalc
  return
end function fMdot_rot_calc
!=======================================================================
double precision function Beasor20()
!*** mass-loss rates proposed by Maeder on the basis of figures in Crowther (2001),
! observations of Sylvester et al 1998 and van Loon et al. (LMC) 1999
  implicit none

  real(kindreal):: dotm
!----------------------------------------------------------------------
  write(io_logs,*) 'Beasor20 Mdot'
  dotm = -26.4 - 0.23*xmini + 4.8*log10(gls)
  Beasor20 = 10.d0**dotm

end function Beasor20
!=======================================================================
real(8) function Bestenlehner20()
!*** Bestenlehner (2020) prescription for hot stars
!*** with fitting parameters from Brands et al. (2022)
  implicit none

  real(kindreal):: dotm
!-----------------------------------------------------------------------
    dotm = -5.19d0 + 2.69d0 * log10(eddesc) - 3.19d0 * log10(1-eddesc)
    Bestenlehner20 = 10.0d0**dotm

end function Bestenlehner20

!=======================================================================
double precision function Bjorklund23()
!*** Bjorklund et al. (2023) prescription for hot stars (Eq. 7)
  use inputparam,only: zsol
  implicit none

  real(kindreal):: dotm

!----------------------------------------------------------------------
  dotm = -5.52d0 + 2.39d0*(log10(gls)-6.0d0) - 1.48d0*log10(gms*(1.0d0-eddesc)/45.0d0) &
      + 2.12d0*(log10(teff)-log10(4.5d4)) &
      + (0.75d0-1.87d0*(log10(teff)-log10(4.5d4))) * log10(zheavy/zsol)
  Bjorklund23 = 10.0d0**dotm

end function Bjorklund23
!=======================================================================
double precision function Crowther01()
!*** mass-loss rates proposed by Maeder on the basis of figures in Crowther (2001),
! observations of Sylvester et al 1998 and van Loon et al. (LMC) 1999
  implicit none

  real(kindreal):: dotm
!----------------------------------------------------------------------
  write(io_logs,*) 'Crowther01 Mdot'
  dotm = -(1.7d0*log10(gls)-13.83d0)
  Crowther01 = 10.d0**(-dotm)

end function Crowther01
!=======================================================================
double precision function Decin24() ! - [MM]
! 2024A&A...681A..17D

  implicit none

  real(kindreal) :: dotm
!----------------------------------------------------------------------

  dotm = -20.63d0 - 0.16d0 * xmini + 3.47d0 * log10(gls)
  Decin23 = 10.d0**dotm

end function Decin24
!=======================================================================
double precision function deJager88()
!***de Jager et al 88 mass loss
  implicit none

  real(kindreal),parameter:: a00=6.34916d0,a01=-5.04240d0,a02=-0.83426d0,a03=-1.13925d0, &
    a04=-0.12202d0,a10=3.41678d0,a11=0.15629d0,a12=2.96244d0,a13=0.33659d0,a14=0.57576d0, &
    a20=-1.08683d0,a21=0.41952d0,a22=-1.37272d0,a23=-1.07493d0,a30=0.13095d0,a31=-0.09825d0, &
    a32=0.13025d0,a40=0.22427d0,a41=0.46591d0,a50=0.11968d0
  real(kindreal):: xlgfz,xxx,yyy,t2x,t2y,t3x,t3y,t4x,t4y,t5x,dotm
!----------------------------------------------------------------------
  if (zheavy > (zsol-1.0d-3).and.zheavy < (zsol+1.0d-3)) then
    xlgfz=0.0d0
  else
    xlgfz=0.5d0*log10(zheavy/zsol)
  endif

  xxx = (log10(teff)-4.05d0)/0.75d0
  yyy = min(((log10(gls)-4.6d0)/2.1d0),1.d0)
  t2x = cos(2.d0*acos(xxx))
  t2y = cos(2.d0*acos(yyy))
  t3x = cos(3.d0*acos(xxx))
  t3y = cos(3.d0*acos(yyy))
  t4x = cos(4.d0*acos(xxx))
  t4y = cos(4.d0*acos(yyy))
  t5x = cos(5.d0*acos(xxx))
  dotm = a00+a01*yyy+a10*xxx+a02*t2y+a11*xxx*yyy+a20*t2x+a03*t3y+a12*xxx*t2y+a21*t2x*yyy+a30*t3x+a04*t4y+ &
         a13*xxx*t3y+a22*t2x*t2y+a31*t3x*yyy+a40*t4x+a14*xxx*t4y+a23*t2x*t3y+a32*t3x*t2y+a41*t4x*yyy+a50*t5x
  deJager88 = 10.d0**(-dotm)*10.d0**xlgfz

end function deJager88

!=======================================================================
double precision function deJager88_lin() ! - [MM]
!***de Jager et al 1988 mass loss, linear approximation
  implicit none

  real(kindreal):: dotm
!----------------------------------------------------------------------
  dotm = -8.158d0 + 1.769d0 * log10(gls) - 1.676d0 * log10(teff)

  deJager88_lin = 10.d0**dotm

end function deJager88_lin

!=======================================================================
double precision function Decin21() ! - [MM]
!*** Mass loss according to Decin (2021)
  implicit none

  real(kindreal) :: dotm
!----------------------------------------------------------------------

  dotm = -1.28d0 + 1.62d0 * log10(gls) - 2.91d0 * log10(teff) - 0.675d0 * log10(gms)

  Decin21 = 10.d0**dotm

end function Decin21


!=======================================================================
double precision function Gormaz22(xsurf,ysurf,y3surf)
!*** Rates from Gormaz-Matamala 2022A&A...665A.133G
  use const,only: cst_G,Msol,Rsol
  implicit none

  real(kindreal),intent(in):: xsurf,ysurf,y3surf
  real(kindreal):: gmrstar,gmlogg,lteff,hehratio,xlmdot
!----------------------------------------------------------------------
  gmrstar=sqrt(gls)*(5777.d0/teff)**2 ! Rstar/Rsun
  gmlogg=log10(cst_G)+log10(gms*Msol)-2*log10(gmrstar*Rsol) !cstlogG + Log10[u2*Msol] - 2 xrad - 2 Log10[Rsol]
  lteff=log10(teff/1000)
  xlmdot=-40.314+15.438*lteff+45.838/gmlogg-8.284*lteff/gmlogg+1.0564*gmrstar
  xlmdot=xlmdot-lteff*gmrstar/2.36-1.1967*gmrstar/gmlogg+11.6*xlogz
  xlmdot=xlmdot-4.223*lteff*xlogz-16.377*xlogz/gmlogg+(gmrstar*xlogz)/81.735
  hehratio=0.25*(ysurf+y3surf)/xsurf
  xlmdot=xlmdot+0.0475-0.559*hehratio
  if (verbose) then
    write(*,*)'T_eff:',teff
    write(*,*)'log g:',gmlogg
    write(*,*)'xlmdot:',xlmdot
    write(*,*)'He/H:',hehratio
  endif
  Gormaz22 = 10.d0**xlmdot

end function Gormaz22
!======================================================================
double precision function Graefener08(xsurf,ysurf,c12surf,o16surf)
! formulae 3, 5, and 6 of Graefener & Hamann (2008A&A...482..945G)
  use inputparam,only: ipop3,zinit,zsol,Z_dep

  real(kindreal),intent(in):: xsurf,ysurf,c12surf,o16surf
  real(kindreal),parameter:: gram0=-3.763d0
  real(kindreal):: zeta,gbetaz,gtest,ggam0,gledd,xlmdot
!----------------------------------------------------------------------
  zeta=(c12surf+o16surf)/ysurf
  gbetaz=1.727d0+0.25d0*xlogz
  ggam0=0.326d0-0.301d0*xlogz-0.045d0*xlogz*xlogz
  gledd=log10(eddesc-ggam0)

  if (xsurf >= 0.05d0) then
    gtest=log10(teff)-4.65d0
    xlmdot=gram0+gbetaz*gledd-3.5d0*gtest+0.42d0*(log10(gls)-6.3d0)-0.45d0*(xsurf-0.4d0)
  else if (xsurf > 0.d0) then
! pour WNE:
    xlmdot=-13.60d0+1.63d0*log10(gls)+2.22d0*log10(ysurf)+Z_dep*xlogz
  else
    if (zeta <= 0.03d0) then
      xlmdot=-13.60d0+1.63d0*log10(gls)+2.22d0*log10(ysurf)+Z_dep*xlogz
! pour WC + WO:
    else
      if (zinit > zsol) then
        xlmdot=-8.30d0+0.84d0*log10(gls)+2.04d0*log10(ysurf)+1.04d0*log10(zheavy)+0.40d0*xlogz
      else if (zinit > 0.002d0) then
        xlmdot=-8.30d0+0.84d0*log10(gls)+2.04d0*log10(ysurf)+1.04d0*log10(zheavy)+0.66d0*xlogz
      else
        if (ipop3 == 1) then
          xlmdot = -8.30d0+0.84d0*log10(gls)+2.04d0*log10(ysurf)+1.04d0*log10(zheavy) &
                   +0.66d0*log10(0.002d0/zsol)+0.35d0*log10(zheavy/0.002d0)
        else
          xlmdot = -8.30d0+0.84d0*log10(gls)+2.04d0*log10(ysurf)+1.04d0*log10(zheavy) &
                   +0.66d0*log10(0.002d0/zsol)+0.35d0*log10(zinit/0.002d0)
        endif
      endif
    endif
  endif
  Graefener08=10.d0**xlmdot

end function Graefener08
!=======================================================================
double precision function Grafener21(D) ! - [MM]
!*** Mass loss according to Grafener (2021)

  use rotmod, only: omegi
  use const, only: cst_G, Msol, Rsol

  implicit none

  real(kindreal), intent(in) :: D ! Clumping factor
  real(kindreal) :: dotm, gmrstar

!----------------------------------------------------------------------
  gmrstar = sqrt(gls)*(5777.d0/teff)**2 ! Rstar/Rsun
  dotm = 5.22d0 * (eddesc + (0.5d0 * omegi(1)**2.d0 * (gmrstar * Rsol)**3.d0) / (cst_G * gms * Msol)) - 0.5d0 * log10(D) - 2.6d0

  Grafener21 = 10.d0**dotm

end function Grafener21
!=======================================================================
double precision function Hainich15(ysurf, ysurf3) ! - [MM]

  implicit none

  real(kindreal), intent(in) :: ysurf, ysurf3
  real(kindreal):: dotm
!----------------------------------------------------------------------

  dotm = -5.13d0 + 0.63d0 * log10(gls) - 0.23d0 * log10(teff) + 1.3d0 * log10(ysurf + ysurf3) + 1.02d0 * xlogz

  Hainich15 = 10.d0**dotm

end function Hainich15


!=======================================================================
double precision function Kee21()
!*** mass-loss rates proposed by Kee+ 2021 (2021A&A...646A.180K)
  use const,only: pi,Msol,Lsol,cst_G,year,cst_k,cst_mh,cst_c,lgLsol,cstlg_sigma,lgpi

  implicit none

  real(kindreal):: vturb,cs,lgrstar3,rstar3,mstar3,cseff,kappa,rpmod,a,b,b1,c,rho3,dotm
!----------------------------------------------------------------------
  write(io_logs,*) 'Kee21 Mdot'
  vturb=18.2d5 ! recommended value for turbulence velocity in Kee+ 2021
  cs = sqrt(cst_k*teff/cst_mh) ! sound velocity as in Kee+ 2021 p.4 (in text)
  lgrstar3 = 0.5d0*(log10(gls)-4.d0*log10(teff)+lgLsol-log10(4.d0)-lgpi-cstlg_sigma)
  rstar3 = 10.d0**lgrstar3
  mstar3 =gms*Msol
  cseff = sqrt(cs**2.d0+vturb**2.d0) ! effective sound velocity as in Kee+ 2021 p.2 (in text)
  kappa = (eddesc*4.d0*pi*cst_G*mstar3*cst_c)/(gls*Lsol) ! mean opacity from Edd factor def by Kee p.2
  rpmod = (cst_G*mstar3*(1.d0-eddesc))/(2.d0*cseff**2.d0) ! modified Parker radius by Kee Eq.5 p.2

  a = (rpmod/(kappa*rstar3**2.d0))  ! step to compute Eq.14
  b = -((2.d0*rpmod)/rstar3)+(3.d0/2.d0)  ! step to compute Eq.14
  b1 = exp(b)  ! step to compute Eq.14
  c = 1.d0-exp(-(2.d0*rpmod)/rstar3)  ! step to compute Eq.14
  rho3 = (4.d0/3.d0)*a*b1/c  ! density in terms of rpmod by Kee Eq.14 p.3

  dotm=4.d0*pi*rho3*cseff*rpmod**2.d0 ! mass loss in cgs by Kee Eq.13 p.3
  Kee21 = (dotm*year)/Msol ! mass loss in Msol/yr

end function Kee21

!=======================================================================
double precision function Krticka21() ! - [MM]
  !*** Mass loss according to Krticka & al. (2021)
  implicit none

  real(kindreal) :: dotm, TeffkK
!----------------------------------------------------------------------

  TeffkK = Teff / 1000.d0 ! Effective temperature in kilo Kelvin

  dotm = - 24.228d0 + 1.5d0 * (log10(gls) - 6.d0) &
         + 24.228d0 * log10( exp(-((TeffkK-14.1d0)/4.88d0)**2.d0) + 5.82d0 * exp(-((TeffkK-37.3d0)/58.8d0)**2.d0))

  Krticka21 = 10.d0**dotm

end function Krticka21


!=======================================================================
double precision function Kudritzki02()
!*** Kudritzki 2002 rates
  use inputparam,only: zsol
  implicit none

  real(kindreal):: xteffcond,xlmdot,azs,als,als2,azmin,aqmin,aq0,aq1
!----------------------------------------------------------------------
! on contraint logTeff au domaine de validite de la formule de Kudritzki

  if (log10(teff) > 4.778d0) then
    xteffcond=4.778d0
  else
    xteffcond=log10(teff)
  endif

  azs=log10(zheavy/zsol)
  als=log10(gls)-6.d0
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
    Kudritzki02=10.d0**xlmdot
  else
    Kudritzki02 = 0.d0
    write(*,*) 'OB_MDOT 6: xmdot set to 0.'
    write(io_logs,*) 'OB_MDOT 6: azs-azmin<0 --> xmdot set to 0.'
  endif

end function Kudritzki02
!======================================================================
double precision function KP2000(xcen,xsurf,ysurf)
!*** Kudritzki et Puls (2000ARA&A..38..613K)
  use const,only: lgpi,qapicg,cst_G,Msol,Rsol,lgRsol,lgLsol,xlsomo,cstlg_sigma,cst_thomson,cst_avo
  implicit none

  real(kindreal),intent(in):: xcen,xsurf,ysurf
  real(kindreal):: xlgfz,xmdvir,xqhe,xepsi,xsigme,xgame,xmasef,xlgrrs,xrrs,xvescp,xvinfi
!----------------------------------------------------------------------
  if (zheavy > (zsol-1.0d-3).and.zheavy < (zsol+1.0d-3)) then
    xlgfz=0.0d0
  else
    xlgfz=0.5d0*log10(zheavy/zsol)
  endif

! Computation of D_mom according to Eq. 10 + Table 2
  if (xcen > 0.d0) then
    xmdvir = 19.87d0+1.57d0*log10(gls)
  else
    if (log10(teff) > 4.5d0) then
      xmdvir = 20.69d0+1.51d0*log10(gls)
    else if (log10(teff) <= 4.5d0 .and. log10(teff) > 4.4d0) then
      xmdvir = 21.24d0+1.34d0*log10(gls)
    else if (log10(teff) <= 4.4d0 .and. log10(teff) > 4.275d0) then
      xmdvir = 17.07d0+1.95d0*log10(gls)
    else
      xmdvir = 14.22d0+2.64d0*log10(gls)
    endif
  endif
  xmdvir = xmdvir + xlgfz - 30.799531d0
! Computation of V_inf:
! a - computation of the Eddington factor
  if (teff >= 35000.d0) then
    xqhe = 2.0d0
  else if (teff >= 30000.d0.and.teff < 35000.d0) then
    xqhe = 1.5d0
  else if (teff >= 25000.d0.and.teff < 30000.d0) then
    xqhe = 1.0d0
  else
    xqhe = 0.0d0
  endif
  xepsi = ysurf/(4.d0*xsurf+ysurf)
  xsigme = cst_thomson*cst_avo*(1.d0+(xqhe-1.d0)*xepsi)/(1.d0+3.d0*xepsi)
  xgame = xlsomo*xsigme*(gls/gms)/qapicg
! b - computation of the effective mass
  xmasef = gms*(1.d0-xgame)
! c - computation of Rstar
  xlgrrs = 0.5d0*(log10(gls)-4.d0*log10(teff)+lgLsol-log10(4.d0)-lgpi-cstlg_sigma-2.d0*lgRsol)
  xrrs = 10.d0**xlgrrs
! d - computation of V_esc in km/s.
  xvescp = sqrt(2.d0*cst_G*Msol/Rsol)*sqrt(xmasef/xrrs)/1.d5
! e - Computation of V_inf according to Eq. 9
  if (teff >= 21000.d0) then
    xvinfi = 2.65d0*xvescp
  else if (teff <= 10000.d0) then
    xvinfi = xvescp
  else
    xvinfi = 1.4d0*xvescp
  endif
  KP2000 = 10.d0**(xmdvir)/(xvinfi*sqrt(xrrs))

end function KP2000


!======================================================================
double precision function Langer89(ysurf, c12surf, o16surf) ! - [MM]
!*** Mass loss according to Langer & al. (1989a) /!\ Only for H free WR

implicit none

real(kindreal), intent(in) :: ysurf, c12surf, o16surf
real(kindreal) :: dotm, zeta
!----------------------------------------------------------------------

  zeta = (c12surf + o16surf) / ysurf

  if( zeta <= 0.03d0 ) then ! for WN
    dotm = -7.2d0 + 2.5d0 * log10(gms)
  else ! for WC
    dotm = -7.d0 + 2.5d0 * log10(gms)
  endif

  Langer89 = 10.d0**dotm

end function Langer89
!======================================================================
double precision function Vink17() ! - [MM]
!*** Mass loss according to Vink (2017)
  implicit none

  real(kindreal) :: dotm
!----------------------------------------------------------------------

  dotm = -13.3d0 + 1.36d0 * log10(gls) + 0.61d0 * log10(zheavy/zsol)
  Vink17 = 10.d0**dotm

end function Vink17
!======================================================================
double precision function Nieuwenhuijzen90() ! - [MM]
!*** Mass loss according to Nieuwenhuijzen & al. (1990)
  implicit none

  real(kindreal) :: dotm
!----------------------------------------------------------------------

  dotm = 1.64d0 * log10(gls) + 0.16d0 * log10(gms) - 1.61d0 * log10(teff) - 7.93d0

  Nieuwenhuijzen90 = 10.d0**dotm

end function Nieuwenhuijzen90
!======================================================================
double precision function Nugis00(xsurf,ysurf,c12surf,o16surf)
!*** Nugis & Lamers 2000
  use inputparam,only: ipop3,zinit,zsol,Z_dep
  implicit none

  real(kindreal),intent(in):: xsurf,ysurf,c12surf,o16surf
  real(kindreal):: ygls,zeta,xlmdot
!----------------------------------------------------------------------
  ygls = log10(gls)
  zeta=(c12surf+o16surf)/ysurf
! pour WN:
  if (xsurf > 0.d0.or.zeta <= 0.03d0) then
    xlmdot=-13.60d0+1.63d0*ygls+2.22d0*log10(ysurf)+Z_dep*xlogz
! pour WC + WO:
  else
    if (zinit > zsol) then
      xlmdot=-8.30d0+0.84d0*ygls+2.04d0*log10(ysurf)+1.04d0*log10(zheavy)+0.40d0*xlogz
    else if (zinit  >  0.002d0) then
      xlmdot=-8.30d0+0.84d0*ygls+2.04d0*log10(ysurf)+1.04d0*log10(zheavy)+0.66d0*xlogz
    else
      if (ipop3 == 1) then
        xlmdot = -8.30d0+0.84d0*ygls+2.04d0*log10(ysurf)+1.04d0*log10(zheavy)+0.66d0*log10(0.002d0/zsol)+ &
                 0.35d0*log10(zheavy/0.002d0)
      else
        xlmdot = -8.30d0+0.84d0*ygls+2.04d0*log10(ysurf)+1.04d0*log10(zheavy)+0.66d0*log10(0.002d0/zsol)+ &
                 0.35d0*log10(zinit/0.002d0)
      endif
    endif
  endif
  Nugis00=10.d0**xlmdot

end function Nugis00
!=======================================================================
double precision function Nugis00_bis(ysurf, ysurf3) ! - [MM]
!*** Unified version of the mass loss according to Nugis & al. (2000)

  implicit none

  real(kindreal), intent(in) :: ysurf, ysurf3
  real(kindreal) :: dotm
!----------------------------------------------------------------------

  dotm = -11.d0 + 1.29d0 * log10(gls) + 1.73d0 * log10(ysurf + ysurf3) + 0.47d0 * xlogz

  Nugis00_bis = 10.d0**dotm

end function Nugis00_bis
!=======================================================================
double precision function Reimers75()
!*** formule de Reimers, etaR donne par fmlos
! rayon en unite de rayon solaire
  use const,only: lgLsol,lgpi,cstlg_sigma,lgRsol

  implicit none

  real(kindreal):: xrsol,eta_Reimers
!----------------------------------------------------------------------
  if (xmini < 5.5d0) then
    eta_Reimers = 0.5d0
  elseif (xmini >= 5.5d0) then
    eta_Reimers = 0.6d0
  endif

  xrsol = 0.5d0*(log10(gls)-4.d0*log10(teff)+lgLsol-log10(4.d0)-lgpi-cstlg_sigma)-lgRsol
  xrsol = 10.d0**xrsol
  Reimers75 = eta_Reimers*4.0d-13*gls*xrsol/gms
end function Reimers75
!======================================================================
double precision function Sabhahit22() ! - [MM]
!*** Mass loss according to Sabhahit et al. (2022)

  use const, only : cst_G, cst_sigma, Lsol, Msol, pi

  implicit none

  real(kindreal) :: dotm, rstar, vesc, vesc_eff, vinf
!----------------------------------------------------------------------

  rstar    = sqrt(gls*Lsol/(4.d0*pi*cst_sigma))/(teff**2.d0)
  vesc     = sqrt(2.d0*cst_G*gms*Msol/rstar)
  vesc_eff = vesc * sqrt(1 - eddesc)
  vinf     = 2.25d0 * alpha_winds / (1 - alpha_winds) * vesc_eff
  dotm     = -26.032d0 + 4.77d0 * log10(gls) - 3.99d0 * log(gms) - 1.226d0 * vinf / vesc_eff + 0.761d0 * log10(zheavy/zsol)

  Sabhahit22 = 10.d0**dotm

end function Sabhahit22

!======================================================================
double precision function Salasnich99() ! - [MM]
!*** Mass loss according to Salashich & al. (1999)

  implicit none

  real(kindreal) :: dotm
  !----------------------------------------------------------------------

  dotm = -14.5d0 + 2.1d0 * log10(gls)
  Salasnich99 = 10.d0**dotm

end function Salasnich99


!======================================================================
double precision function Sander20() ! - [MM]
!*** Mass loss according to Sander & Vink (2020)
  implicit none

  real(kindreal) :: dotm, cbd, gammae
  !----------------------------------------------------------------------

  cbd = 9.15d0 - 0.44d0 * xlogz
  gammae = 10.d0**(-4.813) * gls / gms

  dotm = 2.932d0 * log10(-log10(1.d0-gammae)) - log10(2.d0) *  ((0.244d0 - 0.324d0 * xlogz) / gammae)**cbd &
         + 0.23d0 * xlogz - 2.61d0

  Sander20 = 10.d0**dotm

end function Sander20


!======================================================================
double precision function Schmutz97(xsurf,c12surf,n14surf)
!*** WR mass loss according to Schmutz (1997A&A...321..268S) except for WNL = Nugis+ (1998A&A...333..956N)
  implicit none

  real(kindreal),intent(in):: xsurf,c12surf,n14surf
  real(kindreal):: xmdot,xmdotn
!----------------------------------------------------------------------
  if (xsurf > 1.0d-3) then
    xmdot = 3.0d-05
    xmdotn = 2.4d-08*(gms**2.5d0)
    xmdot = min(xmdot,xmdotn)
  else
    if (c12surf <= n14surf) then
      xmdot = 2.4d-08*(gms**2.5d0)
    else
      xmdot = 2.4d-08*(gms**2.5d0)
    endif
  endif
  Schmutz97 = xmdot
end function Schmutz97

!======================================================================
double precision function Schroder05() ! - [MM]
!*** Mass loss according to Schröder & al. (2005)
!*** The factor -13.1d0 corresponds to $\log(\eta_{sc})$

  use const,      only : Msol, Rsol, cst_G, Teffsol ! /!\ cgs units

  implicit none

  real(kindreal) :: dotm, gsol, g, gmrstar, logetaSC
!----------------------------------------------------------------------

  logetaSC = -13.1d0
  gmrstar  = sqrt(gls)*(Teffsol/teff)**2 ! Rstar/Rsun
  gsol     = cst_G * Msol / Rsol**2.d0
  g        = cst_G * gms * Msol / (gmrstar * Rsol)**2.d0
  dotm     = logetaSC + log10(gls) + log10(gmrstar) - log10(gms) &
             +  3.5d0 * log10(teff/4000.d0) + log10(1 + gsol / (4300.d0 * g))

  Schroder05 = 10.d0**dotm

end function  Schroder05


!======================================================================
double precision function Shenar19(ysurf, ysurf3) ! - [MM]
  !*** Mass loss according to Shenar & al. (2019)
  implicit none

  real(kindreal), intent(in) :: ysurf, ysurf3
  real(kindreal) :: dotm
!----------------------------------------------------------------------

  dotm = -6.22d0 + 0.74d0 * log10(gls) - 0.21d0 * log10(Teff) + 1.42d0 * log10(ysurf + ysurf3) + 0.83d0 * xlogz

  Shenar19 = 10.d0*dotm

end function Shenar19

!======================================================================
double precision function Tramper16(ysurf, ysurf3) ! - [MM]
  !*** Mass loss according to Tramper & al. (2016)
  implicit none

  real(kindreal), intent(in) :: ysurf, ysurf3
  real(kindreal) :: dotm
!----------------------------------------------------------------------

  dotm = -9.2d0 + 0.85d0 * log10(gls) + 0.44d0 * log10(ysurf + ysurf3) + 0.25d0 * xlogz

  Tramper16 = 10.d0**dotm

end function Tramper16


!======================================================================
double precision function Vanbeveren98() ! - [MM]
!*** Vanbeveren & al. (1998). Is applied if Teff < 10kK

  implicit none

  real(kindreal):: dotm
!----------------------------------------------------------------------

  dotm = -8.3d0 + 0.8d0 * log10(gls) + 0.5d0 * log10(zheavy/zsol)
  Vanbeveren98 = 10.d0**dotm


end function Vanbeveren98

!======================================================================
double precision function vanLoon05()
!*** van Loon & al. (2005) for RSG and AGB
  implicit none

  real(kindreal):: xxtt,xxll,xlmdot
!----------------------------------------------------------------------
  write(io_logs,*) 'vanLoon05 Mdot'
  xxtt=log10(teff/3500.d0)
  xxll=log10(gls/10000.d0)
  if (log10(gls) > 4.9d0) then
    xlmdot=-5.3d0+0.82d0*xxll-10.8d0*xxtt
  else
    xlmdot=-5.6d0+1.1d0*xxll-5.2d0*xxtt
  endif
  vanLoon05 = 10.d0**xlmdot

end function vanLoon05
!======================================================================
double precision function Vink01()
!*** Vink et al (2001) mass loss
  use caramodele,only: nwmd
  use inputparam,only: Z_dep

  implicit none

real(kindreal):: charrho,teffjump1,teffjump2,ratio,xlmdot
!----------------------------------------------------------------------
  write(io_logs,*) 'Vink01 Mdot'
  charrho = -14.94d0+3.1857d0*eddesc+Z_dep*xlogz
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
       write(io_logs,'(a,f8.5)')'XLOSS - Teff_jump,1 reached',log10(teffjump1)
       write(*,'(a,f8.5)')'XLOSS - Teff_jump,1 reached',log10(teffjump1)
       write(io_input_changes,'(i7.7,a,f8.5)')nwmd,': Teff_jump,1 reached',log10(teffjump1)
      endif
      if (teff < teffjump2) then
        if (teffv >= teffjump2) then
          write(io_logs,'(a,f8.5)')'XLOSS - Teff_jump,2 reached',log10(teffjump2)
          write(*,'(a,f8.5)')'XLOSS - Teff_jump,2 reached',log10(teffjump2)
          write(io_input_changes,'(i7.7,a,f8.5)')nwmd,': Teff_jump,2 reached',log10(teffjump2)
        endif
        ratio = 0.7d0
        xlmdot = -5.99d0+2.210d0*log10(gls/1.0d+5)-1.339d0*log10(gms/30.d0)-1.601d0*log10(ratio/2.d0)+ &
                 1.07d0*log10(teff/20000.d0)+Z_dep*xlogz
      else if (teff > teffjump2) then
        if (teffv <= teffjump2) then
          write(io_logs,'(a,f8.5)')'XLOSS - Teff_jump,2 reached',log10(teffjump2)
          write(*,'(a,f8.5)')'XLOSS - Teff_jump,2 reached',log10(teffjump2)
          write(io_input_changes,'(i7.7,a,f8.5)')nwmd,': Teff_jump,2 reached',log10(teffjump2)
        endif
        ratio = 1.3d0
        xlmdot = -6.688d0+2.210d0*log10(gls/1.0d+5)-1.339d0*log10(gms/30.d0)-1.601d0*log10(ratio/2.d0)+ &
                 1.07d0*log10(teff/20000.d0)+Z_dep*xlogz
      else
        stop ' STAR at the second jump'
      endif
    else if (teff > teffjump1) then
      if (teffv <= teffjump1) then
       write(io_logs,'(a,f8.5)')'XLOSS - Teff_jump,1 reached',log10(teffjump1)
       write(*,'(a,f8.5)')'XLOSS - Teff_jump,1 reached',log10(teffjump1)
       write(io_input_changes,'(i7.7,a,f8.5)')nwmd,': Teff_jump,1 reached',log10(teffjump1)
      endif
      ratio = 2.6d0
      xlmdot = -6.697d0+2.194d0*log10(gls/1.0d+5)-1.313d0*log10(gms/30.d0)-1.226d0*log10(ratio/2.d0)+ &
               0.933d0*log10(teff/40000.d0)-10.92d0*(log10(teff/40000.d0))*(log10(teff/40000.d0))+0.42d0*xlogz ! 0.42d0 = Z dependency from Vink & Sander (2021) - [MM]
    else
      stop ' STAR at the first jump'
    endif
  else
    ratio = 2.6d0
    xlmdot = -6.697d0+2.194d0*log10(gls/1.0d+5)-1.313d0*log10(gms/30.d0)-1.226d0*log10(ratio/2.d0)+ &
             0.933d0*log10(teff/40000.d0)-10.92d0*(log10(teff/40000.d0))*(log10(teff/40000.d0))+Z_dep*xlogz

  endif
  Vink01 = 10.d0**xlmdot

end function Vink01
!=======================================================================
double precision function V01MP08()
!*** Vink et al (2001, IMLOSS 6) modified by Markova & Puls (2008) + priv. comm. Puls (nov. 2010)
  use inputparam,only: Z_dep
  implicit none

  real(kindreal):: teffjump,ratio,xlmdot
!----------------------------------------------------------------------
    teffjump = 10.d0**(4.3d0)
    if (teff <= teffjump) then
      if (teffv > teffjump) then
       write(io_logs,'(a,f8.5)')'XLOSS - Teff_jump reached, B -> R',log10(teffjump)
        write(*,'(a,f8.5)')'XLOSS - Teff_jump reached, B -> R',log10(teffjump)
        write(io_input_changes,'(i7.7,a,f8.5)')nwmd,': Teff_jump reached, B -> R',log10(teffjump)
      endif
      ratio = 1.4d0
    else
      if (teffv <= teffjump) then
        write(io_logs,'(a,f8.5)')'XLOSS - Teff_jump reached, R -> B',log10(teffjump)
        write(*,'(a,f8.5)')'XLOSS - Teff_jump reached, R -> B',log10(teffjump)
        write(io_input_changes,'(i7.7,a,f8.5)')nwmd,': Teff_jump reached, R -> B',log10(teffjump)
      endif
      ratio=3.0d0
    endif

    xlmdot = -6.697d0+2.194d0*log10(gls/1.0d+5)-1.313d0*log10(gms/30.d0)-1.226d0*log10(ratio/2.d0)+ &
             0.933d0*log10(teff/40000.d0)-10.92d0*(log10(teff/40000.d0))*(log10(teff/40000.d0))+Z_dep*xlogz
    V01MP08 = 10.d0**xlmdot

end function V01MP08
!=======================================================================
double precision function Vink23()
!*** Vink & Sabhahit 2023 (2023A&A...678L...3V)
  implicit none
  real(kindreal):: xlmdot
!----------------------------------------------------------------------
  if (gls <= 4.77) then
    xlmdot = -8.d0 + 0.7d0*log10(gls) - 0.7d0*log10(gms)
  else
    xlmdot = -24.d0 + 4.77d0*log10(gls) - 3.99d0*log10(gms)
  endif

  Vink23 = 10.d0**xlmdot

end function Vink23
!=======================================================================
double precision function Wachter02() ! - [MM]
  !*** Mass loss according to Wachter & al. (2002)
  implicit none

  real(kindreal):: dotm
!----------------------------------------------------------------------

  dotm = 8.86d0 + 2.47d0 * log10(gls) - 1.95d0 * log10(gms) - 6.81d0 * log10(teff)

  Wachter02 = 10.d0**dotm

end function Wachter02



!=======================================================================
double precision function Yang23() ! - [MM]
!*** Mass loss according to Yang & al. (2023)
  implicit none

  real(kindreal):: dotm, loggls
!----------------------------------------------------------------------

  loggls = log10(gls)

  dotm = 0.45d0 * loggls**3.d0 - 5.26d0 * loggls**2.d0 + 20.93d0 * loggls - 34.56d0

  Yang23 = 10.d0**dotm

end function Yang23

!=======================================================================
double precision function Yoon06(xsurf) ! - [MM]
!*** Mass loss according to Yoon & al. (2006)

  implicit none

  real(kindreal), intent(in) :: xsurf
  real(kindreal) :: dotm, loggls
!----------------------------------------------------------------------

  loggls = log10(gls)

  if ( loggls > 4.5 ) then   ! Condition according to Eq. (1) from Yoon & al. (2006)
    dotm = -12.95d0 + 1.5d0 * loggls - 2.85d0 * xsurf + 0.85d0 * xlogz
  else
    dotm = -36.8d0 + 6.8d0 * loggls - 2.85d0 * xsurf + 0.85d0 * xlogz
  endif

  Yoon06 = 10.d0**dotm

end function Yoon06
!=======================================================================
subroutine read_Mdot_prescriptions
  implicit none

  integer:: io_error,i
  character(256):: line
!----------------------------------------------------------------------
  open(file=trim(input_dir)//'inputs/Mdot_recipes.dat',unit=33)
  io_error = 0
  lenf = 1
  do
    read(33,'(a)',iostat=io_error) line
    if (io_error/=0) then
      exit
    else
      lenf = lenf+1
    endif
  enddo
  rewind(33)
  allocate(imdot(lenf-1))
  allocate(bmdot(lenf-1))
  allocate(smdot(lenf-1))

  do i = 1,lenf-1
    read(33,'(i3,2x,a19,2x,a)') imdot(i),bmdot(i),smdot(i)
  enddo

  return
end subroutine read_Mdot_prescriptions
!=======================================================================
function print_Mdot_prescription(iloss) result(MdotRecipe)
  implicit none

  character(:), allocatable:: MdotRecipe
  integer,intent(in):: iloss
  integer:: i
!----------------------------------------------------------------------
  do i = 1,lenf-1
    if (iloss == imdot(i)) then
!      write(io_sfile,*) nwmd,': Mdot prescription from ',trim(smdot(i))
!      write(io_logs,*) 'Mdot prescription from ',trim(smdot(i))
!      write(*,*) '*** Mdot prescription from ',trim(smdot(i))
      MdotRecipe = trim(smdot(i))
      exit
    endif
  enddo
  return

end function print_Mdot_prescription
!=======================================================================
end module winds
