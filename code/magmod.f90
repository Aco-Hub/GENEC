module magmod

  use evol, only: ldi,kindreal

  implicit none

  real(kindreal), dimension(ldi):: D_mago,D_magx,etask,Nmag,bphi,alven,qmin,D_circh

private
public:: D_mago,D_magx,etask,Nmag,bphi,alven,qmin,D_circh
public:: Mag_diff,mag_diff_general

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
!> I only consider the case of fast rotation as explained by Patrick
!! Reference: Stellar evolution with rotation and magnetic fields
!!            Paper 1: A&A (2003) 411, 543
!!            Paper 2: A&A (2004) 422, 225
!!            Paper 3: A&A (2005) 440, 1041
!!            Fuller+2019: 2019MNRAS.485.3661F
subroutine Mag_diff_general(k,zensi,H_P,gravi,Nabla_mu,delt,Nabla_rad,Nabla_ad,rb,omegi,dlodlr,rho,K_ther,alpha_F,n_mag,tb,&
     nsmooth,qminsmooth)
  !-----------------------------------------------------------------------
  use const,only: pi
  use caramodele,only: nwmd
  use nagmod,only: c02agf
  ! Modif B_param
  use strucmod,only: q
  use SmallFunc,only: SmoothProfile,CheckProfile,weighed_smoothing,threshold_smoothing

  implicit none

  integer,intent(in):: k,n_mag,nsmooth
  real(kindreal),intent(in):: alpha_F
  logical,intent(in):: qminsmooth
  real(kindreal),dimension(ldi),intent(in):: zensi,H_P,gravi,Nabla_mu,delt,Nabla_rad,Nabla_ad,rb,omegi,rho,K_ther,tb,dlodlr

  integer:: n,j,l,ifail,jpos,jpo,nroot,nterms,ndegre,mini,mupper,nsmootham,nsmooth_mu
  real(kindreal):: bnmu,bnte,bmos,bq2,bote,bkr,xhs,xbvmag,c_F,coulog,alven_crit,Nabla_mu_thresh,factor_smooth
  real(kindreal),dimension(0:2+2*n_mag):: apol4
  real(kindreal),dimension(3+2*n_mag):: xsolur
  real(kindreal),dimension(2*(2+2*(n_mag+1))):: www4
  real(kindreal),dimension(ldi):: dmago_fast,dmagx_fast,etask_fast,Nvais_fast,bphi_fast,alven_fast,qmin_fast,N2eff,dlodlr_avg &
       ,nabla_mu_avg,Nabla_mu_old,bnmu_avg,D_mago_old
  real(kindreal), dimension(2,2+2*n_mag):: zero4

  logical,parameter:: scale=.true.
  logical:: fast_rot,mag_instab,preserve_sign

  ! Facteur additionnel pour le champ magnetique, cf travaux avec
  ! Patrick. Enlever ou mettre a 1 pour que ce soit "normal".
  real(kind=kindreal), parameter::f_factor = 1.d0 !0.04d0
  save D_mago_old !we save this variable to take an average over time
  !-----------------------------------------------------------------------
  N2eff(:)=0.0d0
  apol4(:)=0.0d0
  D_mago(:)=0.0d0
  D_magx(:)=0.0d0
  etask(:)=0.0d0
  Nmag(:)=0.0d0
  bphi(:)=0.0d0
  alven(:)=0.0d0
  qmin(:)=0.0d0
  dmago_fast(:)=0.0d0
  dmagx_fast(:)=0.0d0
  etask_fast(:)=0.0d0
  Nvais_fast(:)=0.0d0
  bphi_fast(:)=0.0d0
  alven_fast(:)=0.0d0
  qmin_fast(:)=0.0d0
  ! Set up of smoothing variables
  if (nsmooth > 1) then
     mupper=k-(nsmooth+1)
     mini=nsmooth+1
  else
     mini=1
     mupper=k
  endif
  !BLOCK TO SMOOTH NABLA_MU BY USING THE SUBROUTINES 'weighed_smoothing' and 'threshold_smoothing' FROM MESA
  preserve_sign=.TRUE.
  Nabla_mu_avg=Nabla_mu
  Nabla_mu_old=Nabla_mu
  Nabla_mu_thresh=0.0d0 !threshold of Nabla_mu to start smoothing
  nsmooth_mu=5
  call threshold_smoothing(Nabla_mu_avg,Nabla_mu_thresh,k,nsmooth_mu,preserve_sign,Nabla_mu_old)

  !  do n=1,k
  do n=mini,mupper
     ! Smooth shear 
     if (nsmooth > 1) then
        dlodlr_avg(n)=abs(dlodlr(n)) / (2.d0*nsmooth+1.d0)
        do j=1,nsmooth
           dlodlr_avg(n)=dlodlr_avg(n) + ( abs(dlodlr(n-j)) + abs(dlodlr(n+j)) ) / (2.d0*nsmooth+1.d0)
        enddo
     endif

     fast_rot=.false.
     mag_instab=.false.

     if (zensi(n) > 0.0d0) cycle
     if (H_P(n) /= 0.0d0) then
        ! bnmu: N_mu^2 (Paper 1, Eq. 1)
        bnmu=gravi(n)*Nabla_mu_avg(n)/H_P(n)
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
     ! bote: Omega/N_T
     ! bkr: K/(r^2 N_T) (Paper 1, Eq. 11)
     if (bnte /= 0) then
        bote=omegi(n)/(bnte**(0.5d0))
        bkr  = K_ther(n)/(exp(rb(n))*exp(rb(n))*bnte**(0.5d0))
     else
        bote=0.0d0
        bkr=0.0d0
     endif
     
     ! c_F=alpha^3, where alpha is the dimensionless parameter introduced in Fuller+2019
     c_F=alpha_F**3.d0
     
     ! coulog: Coulomb logarithm ln(lambda) (see eq. 5, Wheeler+2015)
     !     if (tb(n) < 11.608235645d0) then ! if T < 1.1x10^5 K
     !        coulog= -17.4d0 + 1.5d0*tb(n) - 0.5d0*rho(n)
     !     else
     !        coulog= -12.7d0 + tb(n) - 0.5d0 * rho(n)
     !     endif
     
     ! Expression of the Coulomb logarithm from Wendell+1987 (1987ApJ...313..284W)
     coulog=log(12.d0* sqrt(4.2d5/exp(tb(n))) * 3.78d-9 * exp(1.5d0 * tb(n)) * exp(-0.5d0*rho(n)) )
     ! dmagx: magnetic diffusivity eta, calculated using Spitzer's formulae (e.g. eq. (5) Wheeler+2015)
     dmagx_fast(n)= 5.2d+11 * coulog * exp(-1.5d0 * tb(n))
     ! etask: eta/K
     etask_fast(n)=dmagx_fast(n)/K_ther(n)
     ! N2eff: effective Brunt-Vaisala frequency Neff^2= eta/k*N_T^2+N_mu^2
     N2eff(n)= etask_fast(n)*bnte + bnmu
     xbvmag=sqrt(N2eff(n))
     if (dmagx_fast(n) < 0.d0) then
        print*, "eta in magmod is negative"
        write(*,*) 'eta=',dmagx_fast(n)
        stop
     endif
     ! qmin: general condition for min q =  1/c * Neff/Omega * (Neff/Omega)^(n/2) * (eta/(r^2*Omega))^(n/4)
     qmin_fast(n)=1.d0/c_F * xbvmag/omegi(n) * (xbvmag/omegi(n))**(real(n_mag)*0.5d0) &
          * (dmagx_fast(n)/bmos) ** (real(n_mag)*0.25d0)
     ! q > qmin ?
     ! If we smooth the qmin condition, the qmin condition is taken into account in the computation of dmago
     if (qminsmooth .eqv. .True. ) then
        mag_instab= .true.
     else
        if (abs(dlodlr_avg(n)) > qmin_fast(n)) then
           mag_instab= .true.
        endif
     endif

     if (mag_instab) then
        if (bq2 <= 0.0d0) cycle
        ifail=0
        jpos=0
        nterms=3+2*n_mag !General number of terms
        ! Not needed with the new way to calculate eta   
        do jpo=1,nterms
           xsolur(jpo)=0.d0    
        enddo
        ndegre=2+2*n_mag !General degree of polynomial
        !  Non zero coefficients of polynomial for x=(omega_A/Omega)^2 (similar to Paper 2, Eq. 25)
        apol4(0)= bmos*bnte/K_ther(n) ! main term a_n
        apol4(2+n_mag)= c_F**2*bnmu*bq2 ! nth power term
        apol4(2+2*n_mag)= - c_F**4*omegi(n)*omegi(n)*bq2*bq2 ! independent term
        call c02agf(apol4,ndegre,scale,zero4,www4,ifail)
        nroot=1
        ! xsolur contains the real and positive roots of the above polynomial.  
        do while (nroot <= ndegre)
           if (zero4(2,nroot) == 0.d0 .and. zero4(1,nroot) > 0.d0) then
              if( zero4(1,nroot) < 1.e+30 .and. zero4(1,nroot) > 1.e-30 ) then ! to avoid over--under flow
                 jpos=jpos+1
                 xsolur(jpos)=zero4(1,nroot)
              endif
           endif
           nroot=nroot+1
        enddo

        ! I ask for the roots to be different by more than 2%
        if (jpos == 2) then
           print*,"TWO ROOTS"
           print*,xsolur(1),xsolur(2)
           if(abs(xsolur(1)-xsolur(2))/abs(xsolur(1)) < 0.02) then
              jpo=1
           endif
           jpos=jpo
        endif
        !  print*,"VALUE OF xhs=",xsolur(jpos)
        !******************************
        ! Case fast rotation (the only one considered, for Omegi >> omega_A)
        !******************************
        !   jpos=1
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
              ! alven: omega_A Alfven frequency for any value of n_mag with unsaturated value of the mag. diffusivity
              alven_fast(n)=omegi(n) * (c_F*abs(dlodlr_avg(n))*omegi(n)/xbvmag)**(1.d0/real(n_mag))
              ! alven_crit: Critical value of Alven frequency, with unsaturated value of mag. diffusivity
              alven_crit= sqrt(xbvmag/exp(rb(n))) * (omegi(n)*dmagx_fast(n))**0.25d0
              ! if Alfven frequency is over the critical one or we don't smooth the qmin condition we
              ! employ the saturated value of the mag. diffusivity and recompute the following quantities
              if (alven_fast(n) > alven_crit .or. qminsmooth .eqv. .False.) then
                 ! dmagx_fast: magnetic diffusivity eta, we use the saturated value if Alfven frequency is over the critical one
                 dmagx_fast(n)=bmos/(c_F**2 * bq2) * xhs**(4+2*n_mag)
                 ! etask: eta/K
                 etask_fast(n)=dmagx_fast(n)/K_ther(n)
                 ! N2eff: effective Brunt-Vaisala frequency Neff^2= eta/k*N_T^2+N_mu^2   
                 N2eff(n)= etask_fast(n)*bnte + bnmu
              endif
              ! bphi: B_phi (Paper 2, Eq. 40)
              bphi_fast(n)=sqrt(4.d0*pi*exp(rho(n)))*exp(rb(n))*alven_fast(n)
              ! dmago_fast: magnetic viscosity, nu(n=n_mag)= Omega*r^2/q * (c*q*Omega/Neff)^(3/n) * (Omega/Neff)
              if(n_mag==1) then  ! n=1 (TS)
                 dmago_fast(n)= min(1.d12, c_F ** 3 * bmos * bq2 * omegi(n)**4/N2eff(n)**2)
              else if (n_mag==2) then ! n=2 --> nu=Omega*r^2 * sqrt(q) * c^(3/2) * (Omega/Neff)^(5/2)
                 dmago_fast(n)= min(1.d12, dmago_fast(n) * bmos * sqrt(abs(dlodlr_avg(n))) * c_F**1.5d0 * &
                      (omegi(n)/sqrt(N2eff(n)))**2.5d0)
              else if (n_mag==3) then ! n=3 (Fuller) --> nu=c*r^2*omega*(omega/Neff)^2 simplified expression
                 dmago_fast(n)= min(1.d12, c_F * bmos * omegi(n)*omegi(n)/N2eff(n))
              endif
              ! Smoothing of dmago profile in non-active regions
              if (qminsmooth .eqv. .True.) then
                 factor_smooth=max(1.d-20, 0.5d0+0.5d0 * tanh( 5.d0*log(alven_fast(n)/alven_crit) ))
                 ! Error message in case of Nan
                 if (isnan( abs(dlodlr_avg(n))/ qmin_fast(n) )) then
                    write(*,*) "stop",abs(dlodlr_avg(n)),qmin_fast(n),n,1.d0-exp(q(n))
                 endif
                 dmago_fast(n)= factor_smooth*dmago_fast(n)
              endif
           endif
        endif
        D_magx(n)=dmagx_fast(n)
        D_mago(n)=dmago_fast(n)
        etask(n)=etask_fast(n)
        Nmag(n)=N2eff(n)
        alven(n)=alven_fast(n)
        bphi(n)=bphi_fast(n)
        qmin(n)=qmin_fast(n)
     else
        D_magx(n)=0.0d0
        D_mago(n)=0.0d0
        etask(n)=0.0d0
        Nmag(n)=bnmu
        alven(n)=0.0d0
        bphi(n)=0.0d0
        qmin(n)=qmin_fast(n)
     endif
  enddo

  do n=mini,mupper
     if(abs(D_mago_old(n)) < 1.d12 .and. D_mago_old(n) /= 0.d0 ) then
        D_mago(n)=0.5d0* (D_mago_old(n) + D_mago(n))
     endif
  enddo
     
  ! Smoothing of magnetic viscosity (nu), once it is calculated
  ! Number of layers used on one side to smooth the magnetic viscosity 
  nsmootham=nsmooth
  if (nsmootham > 1) then
     do n=nsmooth+1, k-nsmooth-1
        ! If the layer is convective or the dynamo is not active we skip that layer
        if ( zensi(n) > 0.d0 .or. D_mago(n)==0.d0) cycle
        D_mago(n)= log10(D_mago(n))
        ! Just to make sure we are not taking log(0)   
        if (isnan(D_mago(n))) stop 'D_mago=0 or Nan'
        j=1
        ! If we are in a radiative zone we multiply the values and add one to the counter j, to account properly for the power of the geometric mean
        do while (j < nsmootham)
           ! Same condition as before but for the layer n+j
           if (zensi(n+j) < 0.d0 .and. D_mago(n+j) .ne. 0.d0) then 
              D_mago(n)= D_mago(n) + log10(D_mago(n+j))
              j = j+1
           else
              exit
           endif
        enddo
        j=j-1 !to come back to the original number of layers taken for the geometric mean
        l=1
        ! If we are in a radiative zone we multiply the values and add one to the the counter l   
        do while (l < nsmootham)
           ! Same condition as before but for the layer n+j      
           if (zensi(n-l) < 0.d0 .and. D_mago(n-l) .ne. 0.d0) then
              D_mago(n)= D_mago(n) + log10(D_mago(n-l))
              l = l+1
           else
              exit
           endif
        enddo
        l=l-1
        ! We divide by the actual number of layers used in the mean, NOT by the maximum number of layers possible
        D_mago(n)=D_mago(n)/real(j+l+1.)
        D_mago(n)= 10.d0 ** D_mago(n)
     enddo
     ! Values near boundaries
     ! inner layers
     do n=k-nsmooth,k
        D_mago(n) = D_mago(n-1)
     enddo

     ! outer layers
     do n=nsmooth,1,-1
        D_mago(n) = D_mago(n+1)
     enddo
  endif

  ! We set eta=0 to avoid mixing of chemical elements
  D_magx(:)= 0.d0 ! we set it equal to 0 --> only consider AMT
  D_mago_old=D_mago ! save for next run
  return
end subroutine Mag_diff_general

end module magmod
