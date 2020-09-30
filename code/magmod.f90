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
subroutine Mag_diff_general(k,zensi,H_P,gravi,Nabla_mu,delt,Nabla_rad,Nabla_ad,rb,omegi,dlodlr,rho,K_ther,alpha_F,n_mag,tb)
!-----------------------------------------------------------------------
  use const,only: pi
  use caramodele,only: nwmd
  use nagmod,only: c02agf
! Modif B_param
!  use strucmod,only:q_mass => q

  implicit none

  integer,intent(in):: k,n_mag
  real(kindreal),intent(in):: alpha_F
  real(kindreal),dimension(ldi),intent(in):: zensi,H_P,gravi,Nabla_mu,delt,Nabla_rad,Nabla_ad,rb,omegi,dlodlr,rho,K_ther,tb

  integer:: n,ifail,jpos,jpo,ndegre,nroot,nterms
  real(kindreal):: bnmu,bnte,bmos,bq2,bomu,bote,bkr,xhs,q0,xbvmag,c_F,coulog,aux
  real(kindreal),dimension(0:2+2*n_mag):: apol4
  real(kindreal),dimension(5):: xsolur
  real(kindreal),dimension(10):: www4
  real(kindreal),dimension(ldi):: dmago_fast,dmagx_fast,etask_fast,Nvais_fast,bphi_fast,alven_fast,qmin_fast, &
    dmago_slow,dmagx_slow,etask_slow,Nvais_slow,bphi_slow,alven_slow,qmin_slow,neff
  real(kindreal), dimension(2,4):: zero4

  logical,parameter:: scale=.true.
  logical:: fast_rot,slow_rot,mag_instab

! Facteur additionnel pour le champ magnetique, cf travaux avec
! Patrick. Enlever ou mettre a 1 pour que ce soit "normal".
  real(kind=kindreal), parameter::f_factor = 1.d0 !0.04d0

  !-----------------------------------------------------------------------
  Neff(:)=0.0d0
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
!  open(40,file='dmago.dat')!,status='new')
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
! c_F=alpha^3, where alpha is the dimensionless parameter introducen in Fuller+2019
  c_F=alpha_F**3
! coulog: Coulomb logarithm ln(lambda) (see eq. 5, Wheeler+2015)
  if (tb(n) < 11.608235645d0) then ! if ln(T) < 1.1x10^5 K
     coulog= -17.4d0 + 1.5d0*tb(n) - 0.5d0*rho(n)
       print*,"low T ",coulog
  else
     coulog= -12.7d0 + tb(n) - 0.5d0 * rho(n)
!     print*,"high T ",coulog
  endif

! bmos: r^2 Omega
   bmos=exp(rb(n))*exp(rb(n))*omegi(n)
! bq2: (dlnOmega/dlnr)^2 = q^2 (Paper 1, Eq. 10)
   bq2=dlodlr(n)*dlodlr(n)
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
   nterms=3+2*n_mag !General number of terms
! Not needed with the new way to calculate eta   
!  do jpo=1,nterms
!     xsolur(jpo)=0.d0    
!  enddo
!  ndegre=2+2*n_mag !General degree of polynomial
!!  Non zero coefficients of polynomial for x=(omega_A/Omega)^2 (similar to Paper 2, Eq. 25)
!  apol4(0)= bmos*bnte/K_ther(n) ! main term a_n
!  apol4(2+n_mag)= c_F**2*bnmu*bq2 ! nth power term
!  apol4(2+2*n_mag)= - c_F**4*omegi(n)*omegi(n)*bq2*bq2 ! independent term
!  call c02agf(apol4,ndegre,scale,zero4,www4,ifail)
!  nroot=1
!  ! xsolur contains the real and positive roots of the above polynomial.  
!  do while (nroot <= ndegre)
!     if (zero4(2,nroot) == 0.d0 .and. zero4(1,nroot) > 0.d0) then
!        if( zero4(1,nroot) < 1.e+30 .and. zero4(1,nroot) > 1.e-30 ) then ! to avoid over--under flow
!           jpos=jpos+1
!           xsolur(jpos)=zero4(1,nroot)
!        endif
!     endif
!     nroot=nroot+1
!  enddo
!
!! I ask for the roots to be different by more than 2%
!  if (jpos == 2) then
!     print*,"TWO ROOTS"
!     print*,xsolur(1),xsolur(2)
!        if(abs(xsolur(1)-xsolur(2))/abs(xsolur(1)) < 0.02) then
!           jpo=1
!        endif
!        jpos=jpo
!  endif
!  
!  print*,"VALUE OF xhs=",xsolur(jpos)
!******************************
! Case fast rotation (the only one considered, for Omegi >> omega_A)
!******************************
   jpos=1
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
! dmagx: magnetic diffusivity eta, calculated using via Spitzer's formulae (e.g. eq. (5) Wheeler+2015)
         dmagx_fast(n)= 5.2d+11 * coulog * exp(-1.5d0 * tb(n))
! etask: eta/K
         etask_fast(n)=dmagx_fast(n)/K_ther(n)
! Neff: effective Brunt-Vaisala frequency Neff^2= eta/k*N_T^2+N_mu^2   
         Neff(n)= etask_fast(n)*bnte + bnmu
         xbvmag=sqrt(Neff(n))
! xhs: omega_A/Omega
         xhs=( c_F * abs(dlodlr(n)) * omegi(n)/xbvmag )** (1.d0/real(n_mag))
! alven: omega_A Alfven frequency
         alven_fast(n)=xhs*omegi(n)
!         alven_fast(n)=(c_F*abs(dlodlr(n))*omegi(n)/sqrt(Neff(n))) ** (1.d0/real(n)) * omegi(n) 
! bphi: B_phi (Paper 2, Eq. 40)
         bphi_fast(n)=sqrt(4.d0*pi*exp(rho(n)))*exp(rb(n))*alven_fast(n)
! dmago: magnetic viscosity nu= Omega*r^2/q * (c*q*Omega/Neff)^(3/n) * (Omega/Neff) (Paper 2, Eq. 45)
         if(n_mag==3) then ! n=3 --> nu=c*r^2*omega*(omega/Neff)^2 simplified expression
            dmago_fast(n)=c_F * bmos * omegi(n)*omegi(n)/Neff(n)
         else if (n_mag==1) then
            dmago_fast(n)= c_F ** 3 * bmos * bq2 * omegi(n)**4/Neff(n)**2! to avoid divide by q
         endif
!         write(40,*) exp(rb(n))/7.d10, log10(dmago_fast(n)),log10(dmagx_fast(n)),omegi(n)
!        dmago_fast(n)=bmos/abs(dlodlr(n)) * (c_F*abs(dlodlr(n))*omegi(n)/xbvmag)**(3.d0/real(n_mag)) * (omegi(n)/xbvmag)         
! bound for Dmago
         if (dmago_fast(n) > 1.d+12) then     !set to upper value
!            print*, "Neff^2=",Neff(n)
!            print*, "WARNING: dmago > 10^16, dmago=",dmago_fast(n),"layer=",n
            dmago_fast(n)=1.d+12
         endif
! qmin: general condition for min q =  1/c * Neff/Omega * (Neff/Omega)^(n/2) * (eta/(r^2*Omega))^(n/4)
         qmin_fast(n)=1.d0/c_F * xbvmag/omegi(n) * (xbvmag/omegi(n))**(real(n_mag)*0.5d0) &
              * (dmagx_fast(n)/bmos) ** (real(n_mag)*0.25d0)
      endif
   endif
!##############################################################
   ! CHOICE OF COEFFICIENTS ACCORDING TO EACH CASE
!##############################################################
   
   if (abs(dlodlr(n)) > qmin_fast(n)) then
      mag_instab=.true.
      if( omegi(n) > alven_fast(n)) then
         fast_rot=.true.
      else if ( omegi(n) < alven_fast(n) ) then
!         print*,"SLOW ROTATION",xhs
!         stop
      endif
   endif

   if (mag_instab) then
!      if (fast_rot) then ! at the moment if q > qmin we assign the values which correspond to fast rotation
         D_magx(n)=dmagx_fast(n)
         D_mago(n)=dmago_fast(n)
         etask(n)=etask_fast(n)
         !        Nmag(n)=Nvais_fast(n)
         Nmag(n)=Neff(n)
         alven(n)=alven_fast(n)
         bphi(n)=bphi_fast(n)
         qmin(n)=qmin_fast(n)
!      endif
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
! We set eta=0 to avoid mixing of chemical elements
D_magx(:)= 0.d0 ! we set it equal to 0 --> only consider AMT
!close(40)
  return
end subroutine Mag_diff_general
!=======================================================================

end module magmod
