! Routine that computes the MRI, equations based on Wheeler(2014) paper and routine based on George's for ST
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
! new MRI switch
  use inputparam,only: mri

  implicit none

  integer,intent(in):: k
  real(kindreal),dimension(ldi),intent(in):: zensi,H_P,gravi,Nabla_mu,delt,Nabla_rad,Nabla_ad,rb,omegi,tb,dlodlr,rho,K_ther

  integer:: n
! MRI - 2 new variables: bontot,limite
  real(kindreal):: bnmu,bnte,bmos,bq2,bomu,bote,bkr,xhs,xbvmag,bontot,limite
! MRI - new array lambdab
  real(kindreal),dimension(ldi):: dmago_cond,dmagx_cond,etask_cond,Nvais_cond,bphi_cond,alven_cond,qmin_cond,lambdab

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
  dmago_cond(:)=0.0d0
  dmagx_cond(:)=0.0d0
  etask_cond(:)=0.0d0
  Nvais_cond(:)=0.0d0
  bphi_cond(:)=0.0d0
  alven_cond(:)=0.0d0
  qmin_cond(:)=0.0d0
  lambdab(:)=0.0d0

  do n=1,k
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
   bq2=dlodlr(n)*dlodlr(n)
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
!dmagx_cond(n)=1.0d+13*exp(-1.5*tb(n))

!Adam Change
!lambdab : Ln(Lambda)=-12.7+ln(T)-0.5ln(rho) as in Wheeler eq(5)  
    lambdab(n)=-12.7d0+tb(n)-0.5d0*rho(n)
    dmagx_cond(n)=5.2d0*(10.d0**11.d0)*lambdab(n)*exp(-1.5*tb(n))

!! etask: eta/K
    etask_cond(n)=dmagx_cond(n)/K_ther(n)
!! Nvais: N^2 in case eta/K << 1 (Paper 2, Eq. 14)
    Nvais_cond(n)=etask_cond(n)*bnte+bnmu
!! alven: omega_A in case eta/K << 1 (Paper 2, Eq. 18)
! alven_cond(n)=sqrt((bq2*omegi(n)*omegi(n)*omegi(n)*omegi(n))/Nvais_cond(n))
    alven_cond(n)=(dmagx_cond(n)*dlodlr(n)*dlodlr(n)/bmos)**(1./6.)*omegi(n)
    xhs=alven_cond(n)/omegi(n)
!! bphi: B_phi (Paper 2, Eq. 40)
    bphi_cond(n)=sqrt(4.d0*pi*exp(rho(n)))*exp(rb(n))*alven_cond(n)
    xbvmag=sqrt(Nvais_cond(n))
!!  !Adam Change
!bontot: Omega/N_tot
    if (Nvais_cond(n) /= 0) then
      bontot=omegi(n)/(xbvmag)
    else
      bontot=0.0d0
    endif



!Adam Change
!MRI switch
    if (mri) then
!! dmago: magnetic viscosity nu (Paper 4, Eq. 13 )
      dmago_cond(n)=0.02d0*abs(dlodlr(n))*omegi(n)*exp(rb(n))*exp(rb(n))
!qmin: condition in Paper 4 Eq. 9 
      qmin_cond(n)=Nvais_cond(n)+2d0*dlodlr(n)*omegi(n)*omegi(n)
      limite=0.0d0
    else 
!!dmago: r^2*Omega*(Omega/N)^4*q^2
      dmago_cond(n)=exp(rb(n))*exp(rb(n))*omegi(n)*bontot*bontot*bontot*bontot*dlodlr(n)*dlodlr(n)
!qmin: condition in Paper 2, above Eq. 17 with omega_A/Omega as above Eq. 18
      qmin_cond(n)=(xbvmag/omegi(n))**(7.d0/4.d0)*(dmagx_cond(n)/(exp(rb(n))*exp(rb(n))*xbvmag))**(0.25d0)
      limite=abs(dlodlr(n))
    endif

  !! Instability condition from Paper 4 Eq. 9
     

    if (limite > qmin_cond(n) .and. omegi(n) > alven_cond(n)) then
      D_magx(n)=dmagx_cond(n)
      D_mago(n)=dmago_cond(n)
      etask(n)=etask_cond(n)
      Nmag(n)=Nvais_cond(n)
      alven(n)=alven_cond(n)
      bphi(n)=bphi_cond(n)
      qmin(n)=qmin_cond(n)
    endif
  enddo

  return

end subroutine Mag_diff
!=======================================================================
end module magmod
