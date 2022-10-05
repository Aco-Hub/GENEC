module timestep

use io_definitions
use evol,only: ldi,kindreal
use const,only: cst_avo,convMeVerg,um,Q_H,Q_He,Q_C
use inputparam,only: verbose,fitm,phase,idifcon,islow,xcn,iadvec
use caramodele,only: gms,hh6,hh1,nwmd
use abundmod,only: x,y,eps,epsy,epsc
use strucmod,only: m,q,s,zensi

implicit none

! AGE, PAS DE TEMPS
real(kindreal),save:: dzeit,dzeitj,dzeitv,alter,xcnwant

private
public :: zeit,TimestepControle
public :: dzeit,dzeitj,dzeitv,alter,xcnwant

contains
!======================================================================
subroutine zeit
!----------------------------------------------------------------------
! modifiee pour que les pas temporels suivent l'evolution de la couche
! la plus active
  use caramodele,only: inum
  use equadiffmod,only: gkorv
  use SmallFunc,only: exphi

  implicit none

  integer:: i, kl, klmax
  real(kindreal):: dzeitvzz,convfactor,cartime,tchem,tchem1,tchem2,fitmoldz,epsxcn,epsxcnk,ratxcn
!----------------------------------------------------------------------
! ConvFactor = 1/(Na*(Mev --> erg conversion))
  ConvFactor=1.d0/(cst_avo*convMeVerg)

  dzeitvzz=dzeit
  if (gkorv < 0.15d0 .and. inum > 0) dzeit=dzeit*2.d0
  if (gkorv > 0.75d0) dzeit=dzeit/2.d0

  i=m-1

  do while (zensi(i+1) > 0.d0 .and. (i+1) > 1)
   i=i-1
  enddo

  if (i < m-3 .and. x(i+1) > 1.d-5) then
    tchem=-1.2d16*gms*exphi(q(i+1))*exp(36.d0*um-hh6)/exphi(s(i+1))
    tchem1=tchem*0.4d0*x(i+1)*0.5d0
    tchem2=tchem*0.1d0*sqrt(x(i+1))*0.5d0
    dzeit=min(tchem1,tchem2,dzeit)
  endif


  do i=m-3,2,-1
! Caracteristic H-burning time.
   if (eps(i) /= 0.d0) then
     Cartime=Q_H/(4.d0*ConvFactor)*x(i)/eps(i)*0.25d0
     if (eps(i) >= eps(i+1) .and. eps(i) > eps(i-1)) dzeit=min(Cartime,dzeit)
   endif
  enddo
  if (hh1 >= 0.02d0) dzeit=dzeit/2.d0
  dzeit=dzeitvzz

! dXi/dt=epsij/Eij, Eij=Qij/n/mH
!   Qij [erg] = 1.602e-6 Qij [MeV]
!   n= number of mass unit of the element
!      considered transformed by the reaction
!   mH = NA^-1 = 1.67e-24 [g]
! dXi(want)=dt*epsij*mH/(Qij/n)

! H-b.: (Qi/n=6.6 MeV)
! dXi(want)=dt*epsij*1.67e-24/(6.6*1.602e-6)=dt*epsij*1.58e-19
! => dt(want) = dXi(want)/epsij*6.33e18

! He-b.: (Qi/n=0.61 MeV)
! dXi(want)=dt*epsij*1.67e-24/(0.61*1.602e-6)=dt*epsij*1.71e-18
! => dt(want) = dXi(want)/epsij*5.85e17

! C-b.: (Qi/n=0.19 MeV)
! dXi(want)=dt*epsij*1.67e-24/(0.19*1.602e-6)=dt*epsij*5.49e-18
! => dt(want) = dXi(want)/epsij*1.82e17

! ratxcn =dt(want)/dt
! ratxcn=dXi(want)/dt/epsij*facxcn,
! facxcn=(6.33e18,5.85e17,1.82e17) for (H,He,C)-b.
! dXi(want)=0.002 (-->0.005)

  epsxcn =-10.d0
  klmax=-10
  ratxcn=12.0d0
  do kl=1,m
   epsxcnk=abs(4.d0*eps(kl)/Q_H+12.d0*epsy(kl)/Q_He+24.d0*epsc(kl)/Q_C)
   epsxcn=max(epsxcn,epsxcnk)
   if (epsxcn == epsxcnk) klmax=kl
  enddo
  epsxcn=abs(eps(klmax)+epsy(klmax)+epsc(klmax))
  if (epsxcn < 2.d0*abs(eps(klmax))) then
    write(io_logs,*) 'max H-b.'
    ratxcn=Q_H/(4.d0*ConvFactor)
  endif
  if (epsxcn < 2.d0*abs(epsy(klmax))) then
    write(io_logs,*) 'max He-b.'
    ratxcn=Q_He/(12.d0*ConvFactor)
  endif
  if (epsxcn < 2.d0*abs(epsc(klmax))) then
    write(io_logs,*) 'max C-b.'
    ratxcn=Q_C/(24.d0*ConvFactor)
  endif
  write(io_logs,*) klmax, 'en./s max= ', epsxcn,' dt ',dzeit
  write(io_logs,*) m,'en./s ', eps(m)+epsy(m)+epsc(m)
  write(io_logs,*) 'en. X,He,C:',eps(klmax),epsy(klmax),epsc(klmax)

! dXi(want)= 0.02/0.005/0.002; 0.002 =value in radiative zone
  ratxcn=0.010d0/dzeit/epsxcn*ratxcn
  ratxcn=0.1d0*ratxcn
  if (phase >= 5 .and. idifcon == 1) ratxcn=3.d0*ratxcn
  if (phase >= 6) ratxcn=0.5d0*ratxcn
  if (epsxcn > 2.d0*abs(epsc(klmax))) ratxcn=2.d0*ratxcn
  if (phase == 10) ratxcn=20.d0*ratxcn
  if (islow == 1) then !dt/3
    ratxcn=0.33d0*ratxcn
  elseif (islow == 2) then !dt/10
    ratxcn=0.1d0*ratxcn
  elseif (islow == 3) then !dt/100
    ratxcn=0.01d0*ratxcn
  elseif (islow == 4) then !dt/25
    ratxcn=0.04d0*ratxcn
  elseif (islow == 5) then
    ratxcn=0.005d0*ratxcn !dt/200
  elseif (islow == 6) then
    ratxcn=0.001d0*ratxcn
  endif
  write(io_logs,*) 'ratio dtwant/dt= ',ratxcn,'eps= ',epsxcn
  write(io_logs,*) 'en. prod= ', (eps(m)+epsy(m)+epsc(m))*dzeit
  write(io_logs,*) 'dtwant= ', ratxcn*dzeit
  if (dzeit /= dzeitvzz) then
    write(io_logs,*)'zeit.fa:dt=',dzeit,' dtv=',dzeitvzz,' xcn=',xcn,' new=',ratxcn
  endif

  dzeit=dzeitvzz

!  if(mod(nwmd,10).eq.1.or.mod(nwmd,10).eq.6) then
  if (mod(nwmd,10)==1 .or. (iadvec==0 .and. mod(nwmd,10)==6 .and. xcn>1.2d0) .or. (iadvec==0 .and. ratxcn<0.5d0)) then
    if (ratxcn < xcn) then
      write(*,*) 'zeit.fa: xcn=',xcn,' m=',m
      write(io_logs,*) 'zeit.fa: xcn=',xcn,' m=',m
      xcn=ratxcn
      xcn=max(nint(10.d0*xcn)/10.d0,0.1d0)
      write(*,*) 'zeit.fa: xcn=',xcn,' couche: ',klmax
      write(io_logs,*) 'zeit.fa: xcn=',xcn,' couche: ',klmax
    endif
    if (xcn < 0.1d0) xcn=0.15d0
    if (xcn > 1.4d0) xcn =1.3d0
    dzeit=dzeit*xcn
  endif
  if (ratxcn < 0.5d0) then
    write(io_input_changes,*) nwmd,': XCN(<0.5)=',xcn,'dzeit:',dzeit
    write(io_input_changes,*) 'New dzeit=',dzeit*xcn,'(ratxcn=',ratxcn
  endif

  if (ratxcn < 0.5d0) then
    if(verbose) then
      write(*,*) 'zeit.f test: ratxcn=',ratxcn
    endif
    write(io_logs,*) 'zeit.f test: ratxcn=',ratxcn
    write(io_sfile,*) 'zeit.f test: ratxcn=',ratxcn
  endif

  if (dzeit /= dzeitvzz) then
    write(io_logs,*) 'zeit.fb: dt=',dzeit,' dtv=',dzeitvzz,' xcn=',xcn
  endif

  fitmoldz=1.d0-exp(q(1))
  if ((fitmoldz > 0.990d0) .and. (abs(fitm-fitmoldz) > min((1.d0-fitmoldz),1.d-3))) then
    dzeit = 0.5d0*dzeit
    xcn = 0.50d0
    write(io_logs,*)'dzeit reduced due to fitm change too big:',dzeit
    write(*,*)'dzeit reduced due to fitm change too big:',dzeit
  elseif (fitmoldz > 0.980d0 .and. abs(fitm-fitmoldz) > 5.d-3) then
    dzeit = 0.5d0*dzeit
    xcn = 0.50d0
    write(io_logs,*)'dzeit reduced due to fitm change too big:',dzeit
    write(*,*)'dzeit reduced due to fitm change too big:',dzeit
  endif

  write(io_logs,*)'dzeit fin:',dzeit

  return

end subroutine zeit
!======================================================================
subroutine TimestepControle(nzmodini)
!-----------------------------------------------------------------------
  use caramodele,only: xmini,xteffprev,xtefflast,xlprev,xllast,xrhoprev,xrholast,&
                       xcprev,xclast
  use inputparam,only: isol,irot,rapcrilim,imloss,icncst,tauH_fit,iprezams
  use rotmod,only: vomegi,rapcri,rapom2,CorrOmega,timestep_control

  implicit none

  integer,intent(in):: nzmodini

  real(kindreal):: ratio_max = 0.05d0
  real(kindreal):: varprev,varlast,newxcnwant,xcnteff,xcnlum,xcnrhoc,RapCorr,stepCritmax,xcnNearCrit,xTolerance,xcnMloss
!-----------------------------------------------------------------------
  stepCritmax = 1.d6

  if (icncst == 0) then
    xcnwant=1.4d0
  elseif (icncst == 1) then
    xcnwant=1.0d0
    write(*,*) '***** ICNCST=1 --> XCN=1.00 *****'
    return
  endif
  varprev=xcprev
  varlast=xclast
  newxcnwant=0.0010d0
  if (xclast >= 0.0d0) then
    if (xclast > 0.3d0) then
      newxcnwant=1.d-3
    else if (xclast < 2.d-4) then
      newxcnwant=1.d0
    else if (xclast < 4.d-4) then
      newxcnwant=1.d-5
    else if (xclast < 8.d-4) then
      newxcnwant=2.d-5
    else if (xclast < 8.d-3) then
      newxcnwant=7.d-5
    else if (xclast < 0.01d0) then
      newxcnwant=2.d-4
    else if (xclast < 0.02d0) then
      newxcnwant=5.d-4
    else if (xclast < 0.04d0) then
      newxcnwant=6.d-4
    else if (xclast < 0.08d0) then
      newxcnwant=8.d-4
    endif

    if (phase >= 4 .and. xclast > 0.02d0) newxcnwant=min(newxcnwant,6.d-4)

    if (varprev >= varlast) then
      if (xclast < 2.d-4) then
        xcnwant=1.4d0
      else
        xcnwant=sqrt(1.d0/((abs(varprev-varlast)+1.d-15)/newxcnwant))
      endif
      if (xcnwant <= 0.d0) then
        write(io_logs,*)'main:xcn<=0',xcnwant
        xcnwant=1.d0
      endif
    endif

    xcnwant=nint(10.d0*xcnwant)/10.d0
    write(*,*)' Critere sur Xc'
    write(*,'(2(a,f20.16),a,f15.5)')' xcprev=',varprev,' xclast=',varlast,' xcn=',xcnwant
  endif

! Assurer limite sur Delta log Teff:
  if (rapom2 < 0.95d0) then
    varprev=xteffprev
    varlast=xtefflast
    newxcnwant=4.d-3
    if (abs(xtefflast-xteffprev) > 5.d-3) then      ! trop grand pas
      newxcnwant=4.d-3
    else if (abs(xtefflast-xteffprev) < 2.d-3) then ! pas trop petit
      newxcnwant=3.d-3
    else                                    ! laisser pas de temps...
      varprev=2.d0
      varlast=1.d0
      newxcnwant=1.d0
    endif
  endif
  xcnteff=sqrt(1.d0/((abs(varprev-varlast)+1.d-15)/newxcnwant))
  xcnteff=nint(10.d0*xcnteff)/10.d0
  if (imloss /= 7 .and. imloss /= 8 .and. phase < 3 .and. iprezams /= 1 .and. rapom2 < 0.90d0) then
    xcnwant=min(xcnwant,xcnteff)
  endif
  write(*,'(a,f10.5)')' Criterion on Teff ',xcnteff
  write(*,'(2(a,f20.16),a,f15.5)')' Teffprev=',xteffprev,' Tefflast=',xtefflast,' xcn=',xcnwant

! Assurer que log l ne change pas de plus que 0.05:
  xcnlum=sqrt(1.d0/((abs(xlprev-xllast)+1.d-15)/0.04d0))
  xcnlum=nint(10.d0*xcnlum)/10.d0
  xcnwant=min(xcnwant,xcnlum)
  write(*,'(a,f10.5)')' Criterion on L ',xcnlum
  write(*,'(2(a,f20.16),a,f15.5)')' xlprev=',xlprev,' xllast=',xllast,' xcn=',xcnwant

! Assurer que log rho_c ne change pas de plus que 0.02:
  xcnrhoc=sqrt(1.d0/((abs(xrhoprev-xrholast)+1.d-15)/0.02d0))
  xcnrhoc=nint(10.d0*xcnrhoc)/10.d0
  xcnwant=min(xcnwant,xcnrhoc)
  write(*,'(a,f10.5)')' Criterion on rho_c ',xcnrhoc
  write(*,'(2(a,f20.16),a,f15.5)')' rhop=',xrhoprev,' rhol=',xrholast,' xcn=',xcnwant

! Limiter XCN a < 1.7:
  if (phase == 2) then
    xcnwant=min(xcnwant,1.7d0)
  else
    if (icncst == 0) then
      xcnwant=min(xcnwant,1.4d0)
    else
      if (icncst == 1) then
        xcnwant=1.0d0
      endif
    endif
  endif

  ! We should not need these lines anymore qith the new way of dealing with the timeseries.
  ! if (nzmodini < 2) then
  !   xcnwant = 1.d0
  ! endif

! [Modif CG / Update SM 2021]
! Close to the chosen maximal velocity, the timestep needs to be decreased so we don't
! crash into the maximal velocity. The maximal timestep is computed as a fraction of the MS lifetime.
! The MS lifetime is found with the following relation:
! if tauH_fit set to 1 (by default):
!   log(tau_H) = A * log(M) + B,
!             with A = -2.632 et B = 9.827 pour M <= 10 Msol
!             and  A = -0.715 et B = 7.819 pour M  > 10 Msol
!   The chosen fraction is dt_crit/tau_H = 2.372e-5 for M <= 10 Msol
!                      and dt_crit/tau_H = 1.902e-5 for M  > 10 Msol
!
! if tauH_fit set to 2:
!   log(tau_H) = A * log(M)**3 + B * log(M)**2 + C * log(M) + D,
!             with A = -0.28, B = 1.96, C = -4.75 et D = 10.4
!   The chosen fraction is dt_crit/tau_H = 2.5e-5

  xcnNearCrit = 10.d0
  if (tauH_fit == 1) then
    if (xmini <= 10.d0) then
      stepCritmax = 10.0d0**(-2.632d0*log10(xmini)+5.202d0)
    else
      stepCritmax = 10.0d0**(-0.715d0*log10(xmini)+3.098d0)
    endif
  elseif (tauH_fit == 2) then
      stepCritmax = 2.d-5*10.0d0**(-0.28*(log10(xmini)**3)+1.96*(log10(xmini)**2)-4.75*(log10(xmini))+10.4)
  endif
  if (rapcrilim > 1.d-5 .and. rapom2 >= 0.98d0*rapcrilim .and. isol == 0) then
    if (dzeitj > stepCritmax) then
      xcnNearCrit = 0.9d0
    else if (dzeitj > 0.75d0*stepCritmax) then
      xcnNearCrit = 1.d0
    endif
  else if (rapcrilim > 1.d-5 .and. rapom2 >= 0.95d0*rapcrilim .and. isol == 0) then
    if (dzeitj > 3.d0*stepCritmax) then
      xcnNearCrit = 0.9d0
    else if (dzeitj > 0.75d0*stepCritmax) then
      xcnNearCrit = 1.d0
    endif
  else if (rapcrilim > 1.d-5 .and. rapom2 >= 0.90d0*rapcrilim .and. isol == 0) then
    if (dzeitj > 6.d0*stepCritmax) then
      xcnNearCrit = 0.9d0
    else if (dzeitj > 0.75d0*stepCritmax) then
      xcnNearCrit = 1.d0
    endif
  endif
  xcnwant = min(xcnwant,xcnNearCrit)

! Avec le nouveau traitement de la perte de masse, il faut faire attention que la correction a appliquer sur
! la premiere couche ne soit pas trop grande. On la limite a  1% de la vitesse de surface.
  if (irot == 1 .and. isol < 1) then
    xcnMloss = 10.d0
    RapCorr = abs(CorrOmega(1)/vomegi(1))
    if (rapom2 < 0.05d0*rapcrilim) then
      xTolerance = 1.d-2
    else if (rapom2 < 0.2d0*rapcrilim) then
      xTolerance = 1.d-2
    else if (rapom2 < 0.5d0*rapcrilim) then
      xTolerance = 5.d-3
    else if (rapom2 < 0.8d0*rapcrilim) then
      xTolerance = 3.d-3
    else if (rapom2 < 0.9d0*rapcrilim) then
      xTolerance = 4.d-3
    else
      xTolerance = 5.d-3
    endif
    if (phase >= 2 .and. rapom2 < 0.7d0*rapcrilim) then
      if (rapom2 >= 0.05d0*rapcrilim) then
        xTolerance = 5.d-3
      else
        xTolerance = 0.2d0
      endif
    endif
    if (RapCorr > xTolerance .and. dzeitj >= stepCritmax/10.d0) then
      if (xtefflast <= 4.0d0 .and. (y(m) > 0.9d0 .or. rapcri < 0.1d0)) then
        xcnMloss = 1.0d0
        if (rapcri < 0.05d0) then
          xcnMloss = 10.d0
        endif
      else
        if (rapcri > 0.005d0) then
          xcnMloss = 0.8d0
          write(io_input_changes,'(i7.7,a)')nwmd+1,': XCN=0.80 (RapCorr)'
        else
          xcnMloss = 1.0d0
        endif
      endif
    endif
    if (RapCorr > 0.1d0 .and. dzeitj >= stepCritmax/10.d0 .and. rapcri >= 5.d-2) then
      if (xtefflast <= 4.0d0 .and. (y(m) > 0.9d0 .or. rapcri < 0.1d0)) then
        xcnMloss = 1.0d0
      else
        xcnMloss = 0.8d0
        write(io_input_changes,'(i7.7,a)')nwmd+1,': XCN=0.80 (RapCorr)'
      endif
    endif
    write(*,*) 'Criterion on CorrOmega:'
    write(*,'(a,d14.8,a,f15.5)') 'CorrOmega(1)/omega(1): ',RapCorr,'  XCN: ',xcnMloss
    if (xcnMloss < 1.0d0 .and. rapom2 > 0.90d0) then
      xcnMloss = 1.0d0
      write(*,*) 'close to critical, xcn floored to 1.'
    endif
    xcnwant = min(xcnwant,xcnMloss)
  endif

  if (timestep_control > ratio_max) then
    if (xcnwant >= 1.0d0 .and. dzeitj > 1.d0) then ! prevents from cutting too much the timestep
      xcnwant = min(xcnwant,0.3d0)
    endif
    write(*,*) 'Criterion on the envelope flux: ', xcnwant
  endif

  write(*,'(a,f4.2,a)') '***** NEW XCN: ',xcnwant,' *****'

  return

end subroutine TimestepControle
!=======================================================================
end module timestep
