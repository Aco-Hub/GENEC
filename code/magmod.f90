module magmod

  use evol, only: ldi,kindreal

  implicit none

  real(kindreal), dimension(ldi):: D_mago,D_magx,etask,Nmag,bphi,alven,qmin,D_circh

private
public:: D_mago,D_magx,etask,Nmag,bphi,alven,qmin,D_circh
public:: Mag_diff

contains
!=======================================================================
!> Computes the magnetic diffusivity
!! Reference: Stellar evolution with rotation and magnetic fields
!!            Paper 1: A&A (2003) 411, 543
!!            Paper 2: A&A (2004) 422, 225
!!            Paper 3: A&A (2005) 440, 1041
subroutine Mag_diff(k,zensi,H_P,gravi,Nabla_mu,delt,Nabla_rad,Nabla_ad,rb,omegi,dlodlr,rho,K_ther)
!-----------------------------------------------------------------------
  use const,only: pi
  use caramodele,only: nwmd
  use nagmod,only: c02agf
! Modif B_param
!  use strucmod,only:q_mass => q

  implicit none

  integer,intent(in):: k
  real(kindreal),dimension(ldi),intent(in):: zensi,H_P,gravi,Nabla_mu,delt,Nabla_rad,Nabla_ad,rb,omegi,dlodlr,rho,K_ther

  integer:: n,ifail,jpos,jpo,ndegre,nroot
  real(kindreal):: bnmu,bnte,bmos,bq2,bomu,bote,bkr,xhs,q0,xbvmag
  real(kindreal),dimension(0:4):: apol4
  real(kindreal),dimension(5):: xsolur
  real(kindreal),dimension(10):: www4
  real(kindreal),dimension(ldi):: dmago_fast,dmagx_fast,etask_fast,Nvais_fast,bphi_fast,alven_fast,qmin_fast, &
    dmago_slow,dmagx_slow,etask_slow,Nvais_slow,bphi_slow,alven_slow,qmin_slow
  real(kindreal), dimension(2,4):: zero4

  logical,parameter:: scale=.true.
  logical:: fast_rot,slow_rot,mag_instab

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
  dmago_fast(:)=0.0d0
  dmago_slow(:)=0.0d0
  dmagx_fast(:)=0.0d0
  dmagx_slow(:)=0.0d0
  etask_fast(:)=0.0d0
  etask_slow(:)=0.0d0
  Nvais_fast(:)=0.0d0
  Nvais_slow(:)=0.0d0
  bphi_fast(:)=0.0d0
  bphi_slow(:)=0.0d0
  alven_fast(:)=0.0d0
  alven_slow(:)=0.0d0
  qmin_fast(:)=0.0d0
  qmin_slow(:)=0.0d0

! Modif B_param
!  open(unit=47, file= "Bdata.dat")

  do n=1,k
   fast_rot=.false.
   slow_rot=.false.
   mag_instab=.false.

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

   if (bq2 <= 0.0d0) cycle
   ifail=0
   jpos=0
   do jpo=1,5
    xsolur(jpo)=0.d0
   enddo
   ndegre=4
! Polynomial for x=(omega_A/Omega)^2 (Paper 2, Eq. 25)
   apol4(0)=f_factor**5.d0*bmos*(bnmu+bnte)/(K_ther(n)*bq2)
   apol4(1)=-f_factor**3.d0*bmos*omegi(n)*omegi(n)/K_ther(n)
   apol4(2)=0.d0
   apol4(3)=2.d0*f_factor**2.d0*bnmu
   apol4(4)=-2.d0*omegi(n)*omegi(n)*bq2
   call c02agf(apol4,ndegre,scale,zero4,www4,ifail)
   nroot=1
! xsolur contains the real and positive roots of the above polynomial.
   do while (nroot <= ndegre)
    if (zero4(2,nroot) == 0.d0 .and. zero4(1,nroot) > 0.d0) then
      jpos=jpos+1
      xsolur(jpos)=zero4(1,nroot)
    endif
    nroot=nroot+1
   enddo

!******************************
! Case fast rotation (Paper 2)
!******************************
! Diffusion coefficients obtained from the root of the 4th degree polynomial
! cf. Paper 2
   if (jpos > 1) then
     write(*,*) " WARNING ! MORE THAN 1 ROOT IN MAG_DIFF "
     do jpo=1,jpos
      write(*,*) " root ",jpo," = ",xsolur(jpo)
     enddo
     rewind(222)
     write(222,*) nwmd," WARNING ! MORE THAN 1 ROOT IN MAG_DIFF"
     stop " WARNING ! MORE THAN 1 ROOT IN MAG_DIFF "
   else
     if (jpos == 1) then
! xhs: omega_A/Omega
       xhs=sqrt(xsolur(1))
! Modif B_param
!       write(47,'(14(3x,e16.8))') 1.d0-exp(q_mass),exp(rb(n))*exp(rb(n)),omegi(n),bmos,bnmu,bnte,K_ther(n),bq2, &
!                      apol4(0),apol4(1),apol4(2),apol4(3),apol4(4),xhs
! dmagx: magnetic diffusivity eta=(r^2 Omega)/q^2 * (omega_A/Omega)^6 (Paper 2, Eq. 19)
       dmagx_fast(n)=bmos/bq2*xhs*xhs*xhs*xhs*xhs*xhs
! etask: eta/K
       etask_fast(n)=dmagx_fast(n)/K_ther(n)
! Nvais: N^2 in case eta/K << 1 (Paper 2, Eq. 14)
       Nvais_fast(n)=etask_fast(n)/2.d0*bnte+bnmu
! alven: omega_A in case eta/K << 1 (Paper 2, Eq. 18)
       alven_fast(n)=sqrt((bq2*omegi(n)*omegi(n)*omegi(n)*omegi(n))/Nvais_fast(n))
! bphi: B_phi (Paper 2, Eq. 40)
       bphi_fast(n)=sqrt(4.d0*pi*exp(rho(n)))*exp(rb(n))*alven_fast(n)
       xbvmag=sqrt(Nvais_fast(n))
! dmago: magnetic viscosity nu (Paper 2, Eq. 45)
       dmago_fast(n)=bmos/abs(dlodlr(n))*xhs*xhs*xhs*omegi(n)/xbvmag
! qmin: condition in Paper 2, above Eq. 17 with omega_A/Omega as above Eq. 18
       qmin_fast(n)=(xbvmag/omegi(n))**(7.d0/4.d0)*(dmagx_fast(n)/(exp(rb(n))*exp(rb(n))*xbvmag))**(0.25d0)
     endif
   endif
!*************************
! Case slow rotation
! cf. Paper 3, appendix A
!*************************
! cf. qmin: Paper 3, Eq. A.6
   q0=bq2*omegi(n)*omegi(n)-bnmu
   qmin_slow(n)=sqrt(abs(bnmu))/omegi(n)
   if (q0 <= 0.d0) then
     dmagx_slow(n)=0.d0
     dmago_slow(n)=0.d0
   else
! cf. eta: Paper 3, Eq. A.8 and A.12, case eta/K << 1
     dmagx_slow(n)=2.d0*K_ther(n)*q0/bnte
     dmago_slow(n)=2.d0*K_ther(n)*q0/bnte
     etask_slow(n)=dmagx_slow(n)/K_ther(n)
! Nvais: N^2 in case eta/K << 1 (Paper 2, Eq. 14)
     Nvais_slow(n)=etask_slow(n)/2.d0*bnte+bnmu
! omega_A: cf Paper 3, Eq. A.9
     alven_slow(n)=(Nvais_slow(n)*dmagx_slow(n)/(exp(rb(n))*exp(rb(n))))**(1.0d0/3.0d0)
     bphi_slow(n)=sqrt(4.d0*pi*exp(rho(n)))*exp(rb(n))*alven_slow(n)
   endif

! Modif B_param
!   qmin_fast(n) = -1.d0
!   qmin_slow(n) = 1.d30

   if (dmagx_fast(n) /= 0.d0 .and. dmagx_slow(n) == 0.d0) then
     fast_rot=.true.
     if (abs(dlodlr(n)) > qmin_fast(n) .and. omegi(n) > alven_fast(n)) then
       mag_instab=.true.
     endif
   else if (dmagx_fast(n) == 0.d0 .and. dmagx_slow(n) /= 0.0d0) then
     slow_rot=.true.
     if (abs(dlodlr(n)) > qmin_slow(n) .and. omegi(n) < alven_slow(n)) then
       mag_instab=.true.
     endif
   else if (dmagx_fast(n) /= 0.d0 .and. dmagx_slow(n) /= 0.d0) then
     if (abs(dlodlr(n)) > qmin_slow(n) .and. omegi(n) < alven_slow(n) .and. &
         abs(dlodlr(n)) > qmin_fast(n) .and. omegi(n) > alven_fast(n)) then
! Both slow and fast rotation conditions: fast rot values applied
       write(3,*) " Mag_diff: conditions for slow and fast rot, layer ",n
       fast_rot=.true.
       mag_instab=.true.
     else if (abs(dlodlr(n)) > qmin_fast(n) .and. omegi(n) > alven_fast(n)) then
       fast_rot=.true.
       mag_instab=.true.
     else if (abs(dlodlr(n)) > qmin_slow(n) .and. omegi(n) < alven_slow(n)) then
       slow_rot=.true.
       mag_instab=.true.
     endif
   endif

   if (mag_instab) then
     if (fast_rot) then
       D_magx(n)=dmagx_fast(n)
       D_mago(n)=dmago_fast(n)
       etask(n)=etask_fast(n)
       Nmag(n)=Nvais_fast(n)
       alven(n)=alven_fast(n)
       bphi(n)=bphi_fast(n)
       qmin(n)=qmin_fast(n)
     else if (slow_rot) then
       D_magx(n)=dmagx_slow(n)
       D_mago(n)=dmago_slow(n)
       etask(n)=etask_slow(n)
       Nmag(n)=Nvais_slow(n)
       alven(n)=alven_slow(n)
       bphi(n)=bphi_slow(n)
       qmin(n)=qmin_slow(n)
     endif
   else
     D_magx(n)=0.0d0
     D_mago(n)=0.0d0
     etask(n)=0.0d0
     Nmag(n)=bnmu
     alven(n)=0.0d0
     bphi(n)=0.0d0
     if (slow_rot) then
       qmin(n)=qmin_slow(n)
     else
       qmin(n)=qmin_fast(n)
     endif
   endif
  enddo
! Modif B_param
!  close(47)

  return

end subroutine Mag_diff

!=======================================================================
!> Computes the magnetic diffusivity in a general case, according to the values of n_mag (n_mag=1 TS, n_mag=3 Fuller+2019) and alpha_F.
!> n_mag and alpha_F are input parameters (n=1 , alpha_F=1 by default)
!> n_mag is integer, alpha_F is real
!! Reference: Stellar evolution with rotation and magnetic fields
!!            Paper 1: A&A (2003) 411, 543
!!            Paper 2: A&A (2004) 422, 225
!!            Paper 3: A&A (2005) 440, 1041
!!            Fuller+2019: 2019MNRAS.485.3661F
subroutine Mag_diff_general(k,zensi,H_P,gravi,Nabla_mu,delt,Nabla_rad,Nabla_ad,rb,omegi,dlodlr,rho,K_ther,alpha_F,n_mag)
!-----------------------------------------------------------------------
  use const,only: pi
  use caramodele,only: nwmd
  use nagmod,only: c02agf
! Modif B_param
!  use strucmod,only:q_mass => q

  implicit none

  integer,intent(in):: k,n_mag
  real(kindreal),intent(in):: alpha_F
  real(kindreal),dimension(ldi),intent(in):: zensi,H_P,gravi,Nabla_mu,delt,Nabla_rad,Nabla_ad,rb,omegi,dlodlr,rho,K_ther

  integer:: n,ifail,jpos,jpo,ndegre,nroot
  real(kindreal):: bnmu,bnte,bmos,bq2,bomu,bote,bkr,xhs,q0,xbvmag,c
  real(kindreal),dimension(0:4):: apol4
  real(kindreal),dimension(5):: xsolur
  real(kindreal),dimension(10):: www4
  real(kindreal),dimension(ldi):: dmago_fast,dmagx_fast,etask_fast,Nvais_fast,bphi_fast,alven_fast,qmin_fast, &
    dmago_slow,dmagx_slow,etask_slow,Nvais_slow,bphi_slow,alven_slow,qmin_slow
  real(kindreal), dimension(2,4):: zero4

  logical,parameter:: scale=.true.
  logical:: fast_rot,slow_rot,mag_instab

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
  dmago_fast(:)=0.0d0
  dmago_slow(:)=0.0d0
  dmagx_fast(:)=0.0d0
  dmagx_slow(:)=0.0d0
  etask_fast(:)=0.0d0
  etask_slow(:)=0.0d0
  Nvais_fast(:)=0.0d0
  Nvais_slow(:)=0.0d0
  bphi_fast(:)=0.0d0
  bphi_slow(:)=0.0d0
  alven_fast(:)=0.0d0
  alven_slow(:)=0.0d0
  qmin_fast(:)=0.0d0
  qmin_slow(:)=0.0d0

! Modif B_param
!  open(unit=47, file= "Bdata.dat")

  do n=1,k
   fast_rot=.false.
   slow_rot=.false.
   mag_instab=.false.

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

   if (bq2 <= 0.0d0) cycle
   ifail=0
   jpos=0
   do jpo=1,5
    xsolur(jpo)=0.d0
   enddo
   ndegre=4
! Polynomial for x=(omega_A/Omega)^2 (Paper 2, Eq. 25)
   apol4(0)=f_factor**5.d0*bmos*(bnmu+bnte)/(K_ther(n)*bq2)
   apol4(1)=-f_factor**3.d0*bmos*omegi(n)*omegi(n)/K_ther(n)
   apol4(2)=0.d0
   apol4(3)=2.d0*f_factor**2.d0*bnmu
   apol4(4)=-2.d0*omegi(n)*omegi(n)*bq2
   call c02agf(apol4,ndegre,scale,zero4,www4,ifail)
   nroot=1
! xsolur contains the real and positive roots of the above polynomial.
   do while (nroot <= ndegre)
    if (zero4(2,nroot) == 0.d0 .and. zero4(1,nroot) > 0.d0) then
      jpos=jpos+1
      xsolur(jpos)=zero4(1,nroot)
    endif
    nroot=nroot+1
   enddo

!******************************
! Case fast rotation (Paper 2)
!******************************
! Diffusion coefficients obtained from the root of the 4th degree polynomial
! cf. Paper 2
   if (jpos > 1) then
     write(*,*) " WARNING ! MORE THAN 1 ROOT IN MAG_DIFF "
     do jpo=1,jpos
      write(*,*) " root ",jpo," = ",xsolur(jpo)
     enddo
     rewind(222)
     write(222,*) nwmd," WARNING ! MORE THAN 1 ROOT IN MAG_DIFF"
     stop " WARNING ! MORE THAN 1 ROOT IN MAG_DIFF "
   else
     if (jpos == 1) then
! xhs: omega_A/Omega
       xhs=sqrt(xsolur(1))
! Modif B_param
!       write(47,'(14(3x,e16.8))') 1.d0-exp(q_mass),exp(rb(n))*exp(rb(n)),omegi(n),bmos,bnmu,bnte,K_ther(n),bq2, &
!                      apol4(0),apol4(1),apol4(2),apol4(3),apol4(4),xhs
! dmagx: magnetic diffusivity eta=(r^2 Omega)/q^2 * (omega_A/Omega)^6 (Paper 2, Eq. 19)
       dmagx_fast(n)=bmos/bq2*xhs*xhs*xhs*xhs*xhs*xhs
! etask: eta/K
       etask_fast(n)=dmagx_fast(n)/K_ther(n)
! Nvais: N^2 in case eta/K << 1 (Paper 2, Eq. 14)
       Nvais_fast(n)=etask_fast(n)/2.d0*bnte+bnmu
! alven: omega_A in case eta/K << 1 (Paper 2, Eq. 18)
       alven_fast(n)=sqrt((bq2*omegi(n)*omegi(n)*omegi(n)*omegi(n))/Nvais_fast(n))
! bphi: B_phi (Paper 2, Eq. 40)
       bphi_fast(n)=sqrt(4.d0*pi*exp(rho(n)))*exp(rb(n))*alven_fast(n)
       xbvmag=sqrt(Nvais_fast(n))
! dmago: magnetic viscosity nu (Paper 2, Eq. 45)
       dmago_fast(n)=bmos/abs(dlodlr(n))*xhs*xhs*xhs*omegi(n)/xbvmag
! qmin: condition in Paper 2, above Eq. 17 with omega_A/Omega as above Eq. 18
       qmin_fast(n)=(xbvmag/omegi(n))**(7.d0/4.d0)*(dmagx_fast(n)/(exp(rb(n))*exp(rb(n))*xbvmag))**(0.25d0)
     endif
   endif
!*************************
! Case slow rotation
! cf. Paper 3, appendix A
!*************************
! cf. qmin: Paper 3, Eq. A.6
   q0=bq2*omegi(n)*omegi(n)-bnmu
   qmin_slow(n)=sqrt(abs(bnmu))/omegi(n)
   if (q0 <= 0.d0) then
     dmagx_slow(n)=0.d0
     dmago_slow(n)=0.d0
   else
! cf. eta: Paper 3, Eq. A.8 and A.12, case eta/K << 1
     dmagx_slow(n)=2.d0*K_ther(n)*q0/bnte
     dmago_slow(n)=2.d0*K_ther(n)*q0/bnte
     etask_slow(n)=dmagx_slow(n)/K_ther(n)
! Nvais: N^2 in case eta/K << 1 (Paper 2, Eq. 14)
     Nvais_slow(n)=etask_slow(n)/2.d0*bnte+bnmu
! omega_A: cf Paper 3, Eq. A.9
     alven_slow(n)=(Nvais_slow(n)*dmagx_slow(n)/(exp(rb(n))*exp(rb(n))))**(1.0d0/3.0d0)
     bphi_slow(n)=sqrt(4.d0*pi*exp(rho(n)))*exp(rb(n))*alven_slow(n)
   endif

! Modif B_param
!   qmin_fast(n) = -1.d0
!   qmin_slow(n) = 1.d30

   if (dmagx_fast(n) /= 0.d0 .and. dmagx_slow(n) == 0.d0) then
     fast_rot=.true.
     if (abs(dlodlr(n)) > qmin_fast(n) .and. omegi(n) > alven_fast(n)) then
       mag_instab=.true.
     endif
   else if (dmagx_fast(n) == 0.d0 .and. dmagx_slow(n) /= 0.0d0) then
     slow_rot=.true.
     if (abs(dlodlr(n)) > qmin_slow(n) .and. omegi(n) < alven_slow(n)) then
       mag_instab=.true.
     endif
   else if (dmagx_fast(n) /= 0.d0 .and. dmagx_slow(n) /= 0.d0) then
     if (abs(dlodlr(n)) > qmin_slow(n) .and. omegi(n) < alven_slow(n) .and. &
         abs(dlodlr(n)) > qmin_fast(n) .and. omegi(n) > alven_fast(n)) then
! Both slow and fast rotation conditions: fast rot values applied
       write(3,*) " Mag_diff: conditions for slow and fast rot, layer ",n
       fast_rot=.true.
       mag_instab=.true.
     else if (abs(dlodlr(n)) > qmin_fast(n) .and. omegi(n) > alven_fast(n)) then
       fast_rot=.true.
       mag_instab=.true.
     else if (abs(dlodlr(n)) > qmin_slow(n) .and. omegi(n) < alven_slow(n)) then
       slow_rot=.true.
       mag_instab=.true.
     endif
   endif

   if (mag_instab) then
     if (fast_rot) then
       D_magx(n)=dmagx_fast(n)
       D_mago(n)=dmago_fast(n)
       etask(n)=etask_fast(n)
       Nmag(n)=Nvais_fast(n)
       alven(n)=alven_fast(n)
       bphi(n)=bphi_fast(n)
       qmin(n)=qmin_fast(n)
     else if (slow_rot) then
       D_magx(n)=dmagx_slow(n)
       D_mago(n)=dmago_slow(n)
       etask(n)=etask_slow(n)
       Nmag(n)=Nvais_slow(n)
       alven(n)=alven_slow(n)
       bphi(n)=bphi_slow(n)
       qmin(n)=qmin_slow(n)
     endif
   else
     D_magx(n)=0.0d0
     D_mago(n)=0.0d0
     etask(n)=0.0d0
     Nmag(n)=bnmu
     alven(n)=0.0d0
     bphi(n)=0.0d0
     if (slow_rot) then
       qmin(n)=qmin_slow(n)
     else
       qmin(n)=qmin_fast(n)
     endif
   endif
  enddo
! Modif B_param
!  close(47)

  return

end subroutine Mag_diff_general
!=======================================================================

end module magmod
