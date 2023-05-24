!Added Adam's MRI subroutine into TS subrorutine of GENEC
module magmod

  use evol, only: ldi,kindreal

  implicit none

  real(kindreal), dimension(ldi):: D_mago,D_magx,etask,Nmag,bphi,alven,qmin,D_circh

private
public:: D_mago,D_magx,etask,Nmag,bphi,alven,qmin,D_circh
public:: Mag_diff

contains
!=======================================================================
!> Computes the magnetorotational instability written September 2020
!! Reference: Stellar evolution with rotation and magnetic fields
!!            Paper 1: A&A (2003) 411, 543
!!            Paper 2: A&A (2004) 422, 225
!!            Paper 3: A&A (2005) 440, 1041
!! Reference: The role of the magnetorotational instability in massive stars Wheeler
!!            Paper 4: Wheeler(2014)
subroutine Mag_diff(k,zensi,H_P,gravi,Nabla_mu,delt,Nabla_rad,Nabla_ad,rb,omegi,tb,dlodlr,rho,K_ther)
!-----------------------------------------------------------------------
  use const,only: pi
! new MRI switch and coefficient to toggle chemical gradient efficency
  use inputparam,only: mri !mri=0 -> No magnetic instabilities, mri=1 -> Just MRI included, mri=2-> Just TS included, mri=3 -> Both MRI and TS
  use inputparam,only: fmu !In Paper 4 set at 0.05 for mri. 

  implicit none

  integer,intent(in):: k
  real(kindreal),dimension(ldi),intent(in):: zensi,H_P,gravi,Nabla_mu,delt,Nabla_rad,Nabla_ad,rb,omegi,tb,dlodlr,rho,K_ther

  integer:: n
  real(kindreal):: bnmu,bnte,bmos,bq2,bomu,bote,bkr,xhs,xbvmag



  real(kindreal),dimension(ldi):: dmago_TS,dmagx_TS,dmago_mri,dmagx_rest,dmagx_mri,etask_cond,Nvais_cond,bphi_cond,& 
  alven_cond,qmin_cond_TS,lambdab,qmin_cond_mri,dlodlr_avg


  logical,parameter:: scale=.true.
  logical:: mag_instab,mag_instab_mri

! Facteur additionnel pour le champ magnetique, cf travaux avec
! Patrick. Enlever ou mettre a 1 pour que ce soit "normal".
  real(kind=kindreal), parameter::f_factor = 1.d0 !0.04d0

!-----------------------------------------------------------------------
  D_mago(:)=0.0d0
  D_magx(:)=0.0d0
  etask(:)=0.0d0
  Nmag(:)=0.0d0
  bphi(:)=0.0d0
  alven(:)=0.0d0
  qmin(:)=0.0d0
  dmago_TS(:)=0.0d0
  dmagx_TS(:)=0.0d0
  etask_cond(:)=0.0d0
  Nvais_cond(:)=0.0d0
  bphi_cond(:)=0.0d0
  alven_cond(:)=0.0d0
  qmin_cond_TS(:)=0.0d0
  lambdab(:)=0.0d0
  dmago_mri(:)=0.0d0
  dmagx_mri(:)=0.0d0
  dmagx_rest(:)=0.0d0
  qmin_cond_mri(:)=0.0d0

  !The shear is averaged to give a smoother profile that works better with the qmin condition and avoids too much on/off behaviour between shells
  do n=1,k
   mag_instab=.false.
   mag_instab_mri=.false.
   if (n>2 .and. n<k-1) then
      dlodlr_avg(n)=(dlodlr(n)+dlodlr(n-1)+dlodlr(n-2)+dlodlr(n+1)+dlodlr(n+2))/5.d0
   elseif (n==1) then
      dlodlr_avg(n)=(dlodlr(n)+dlodlr(n+1)+dlodlr(n+2)+dlodlr(n+3)+dlodlr(n+4))/5.d0
   elseif (n==2) then
      dlodlr_avg(n)=(dlodlr(n)+dlodlr(n-1)+dlodlr(n+1)+dlodlr(n+2)+dlodlr(n+3))/5.d0
   elseif (n==k-1) then
      dlodlr_avg(n)=(dlodlr(n)+dlodlr(n+1)+dlodlr(n-1)+dlodlr(n-2)+dlodlr(n-3))/5.d0
   else 
      dlodlr_avg(n)=(dlodlr(n)+dlodlr(n-1)+dlodlr(n-2)+dlodlr(n-3)+dlodlr(n-4))/5.d0
    endif 

   if (zensi(n) > 0.0d0) cycle
   if (H_P(n) /= 0.0d0) then
! bnmu: N_mu^2 (Paper 1, Eq. 1)
     bnmu=gravi(n)*Nabla_mu(n)/H_P(n)
! bnte: N_T^2 (Paper 1, Eq. 2)
     bnte=gravi(n)*delt(n)/H_P(n)*abs(Nabla_rad(n)-Nabla_ad(n))
   else
     bnmu=0.0d0
     bnte=0.0d0
   endif
! bmos: r^2 Omega
   bmos=exp(rb(n))*exp(rb(n))*omegi(n)
! bq2: (dlnOmega/dlnr)^2 = q^2 (Paper 1, Eq. 10)
   bq2=dlodlr_avg(n)*dlodlr_avg(n)
! bomu: Omega/N_mu (Paper 1, Eq. 10)
   if (bnmu /= 0.0d0) then
     bomu=omegi(n)/(bnmu**(0.5d0))
   else
     bomu=0.0d0
   endif
! bote: Omega/N_T
! bkr: K/(r^2 N_T) (Paper 1, Eq. 11)
   if (bnte /= 0) then
     bote=omegi(n)/(bnte**(0.5d0))
     bkr  = K_ther(n)/(exp(rb(n))*exp(rb(n))*bnte**(0.5d0))
   else
     bote=0.0d0
     bkr=0.0d0
   endif

!! dmagx: magnetic diffusivity according to Spitzer
    lambdab(n)=-12.7d0+tb(n)-0.5d0*rho(n)
    dmagx_rest(n)=5.2d0*(10.d0**11.d0)*lambdab(n)*exp(-1.5*tb(n))  !Diffusivity at rest(Paper 4 Eq. 5) 
!! etask: eta/K
    etask_cond(n)=dmagx_rest(n)/K_ther(n)
!! Nvais: N^2 in case eta/K << 1 (Paper 2, Eq. 14)
   Nvais_cond(n)=etask_cond(n)*bnte+bnmu
!! alven: omega_A in case eta/K << 1 (Paper 2, Eq. 18)
! alven_cond(n)=sqrt((bq2*omegi(n)*omegi(n)*omegi(n)*omegi(n))/Nvais_cond(n))
    alven_cond(n)=(dmagx_rest(n)*dlodlr_avg(n)*dlodlr_avg(n)/bmos)**(1./6.)*omegi(n)
    xhs=alven_cond(n)/omegi(n)
!! bphi: B_phi (Paper 2, Eq. 40)
    bphi_cond(n)=sqrt(4.d0*pi*exp(rho(n)))*exp(rb(n))*alven_cond(n)
    xbvmag=sqrt(Nvais_cond(n))

!Implemnation of magnetic instabilities
    if (mri==2 .or. mri==3) then !TS is added into model 
       dmago_TS(n)=bmos*bq2*(omegi(n)/xbvmag)**(4.d0)
       dmagx_TS(n)=dmagx_rest(n)

       qmin_cond_TS(n)=((xbvmag/omegi(n))**(7.d0/4.d0)*(dmagx_rest(n)/(exp(rb(n))*exp(rb(n))*xbvmag))**(0.25d0)) !(Paper 4 Eq. 7)
       
       if (abs(dlodlr_avg(n))>qmin_cond_TS(n)) then !TS is active
           mag_instab=.true.
       endif
    else
       dmago_TS(n)=0.0d0
       dmagx_TS(n)=0.0d0
    endif
    if (mri==1 .or. mri==3) then !MRI is added into model 
       !! dmago: magnetic viscosity nu (Paper 4, Eq. 13 )
       dmago_mri(n)=0.02d0*abs(dlodlr_avg(n))*omegi(n)*exp(rb(n))*exp(rb(n))
       dmagx_mri(n)=dmago_mri(n)
        !qmin: condition in (Paper 4 Eq. 9) 
       qmin_cond_mri(n)=abs(-(etask_cond(n)*bnte+fmu*bnmu)/(2d0*omegi(n)*omegi(n))) !fmu parameter added here in minimim condition
      if ((abs(dlodlr_avg(n))>qmin_cond_mri(n)) .and. (abs(dlodlr_avg(n))<4) ) then !MRI is active
           mag_instab_mri=.true.
           qmin(n) = 1.
      else
         qmin(n) = -1.
      endif
    else
       dmago_mri(n)=0.0d0
       dmagx_mri(n)=0.0d0
     endif


     if (mag_instab_mri .and. mag_instab) then
          D_magx(n)=dmagx_mri(n)+dmagx_TS(n)
          D_mago(n)=dmago_mri(n)+dmago_TS(n)
          if (D_mago(n)>1d+12) then
              D_mago(n)=1d+12
          endif
          if (D_magx(n)>1d+12) then
             D_magx(n)=1d+12
         endif
          etask(n)=etask_cond(n)
          Nmag(n)=Nvais_cond(n)
          alven(n)=alven_cond(n)
          bphi(n)=bphi_cond(n)
      else if ((mag_instab_mri .eqv. .false.) .and. mag_instab) then !MRI not active
          D_magx(n)=dmagx_TS(n)
          D_mago(n)=dmago_TS(n)
          if (D_mago(n)>1d+12) then
              D_mago(n)=1d+12
          endif
          if (D_magx(n)>1d+12) then
             D_magx(n)=1d+12
         endif
          etask(n)=etask_cond(n)
          Nmag(n)=Nvais_cond(n)
          alven(n)=alven_cond(n)
          bphi(n)=bphi_cond(n)
      else if (mag_instab_mri .and. (mag_instab .eqv. .false.)) then !TS not active
          D_magx(n)=dmagx_mri(n)
          D_mago(n)=dmago_mri(n)
          if (D_mago(n)>1d+12) then
              D_mago(n)=1d+12
          endif
          if (D_magx(n)>1d+12) then
             D_magx(n)=1d+12
         endif
          etask(n)=etask_cond(n)
          Nmag(n)=Nvais_cond(n)
          alven(n)=alven_cond(n)
          bphi(n)=bphi_cond(n)
      else !Neither active
          D_magx(n)=0.0d0
          D_mago(n)=0.0d0
          etask(n)=etask_cond(n)
          Nmag(n)=Nvais_cond(n)
          alven(n)=alven_cond(n)
          bphi(n)=bphi_cond(n)
      endif
  enddo

  return

end subroutine Mag_diff
!=======================================================================
end module magmod

  
