module PrintAll
! Module used to collect, organise and print in a file the data used by Saio to compute pulsationnal
! properties of stars. Could be modified for a more general purpose.

use evol, only: kindreal,ldi
use const, only: pi,cst_c,cst_a,cst_G,cst_k,cst_u,Msol,Lsol,cst_sigma

implicit none

character(256), save:: DataAll_FileName

integer, parameter:: File_Unit = 59, Size_Env = 1000, Size_Atmos = 200

integer, save:: Atm_Layers,Env_Layers,Int_Layers

real(kind=kindreal), save:: Teff_save,Lum_save,mass_save,time_save,C12_save,C13_save,N14_save,O16_save
real(kind=kindreal), dimension(ldi), save:: logr_int,m_int,logT_int,logrho_int,logP_int,alpha_int,delta_int,nablaad_int, &
  Gamma1_int,cp_int,cv_int,dlnPdlnrho_int,dlnPdlnT_int,nablae_int,logkappa_int,Lrad_int,L_int,dlnkappadlnrhoT_int, &
  dlnkappadlnTrho_int,epsilon_int,dlnepsilondlnrhoT_int,dlnepsilondlnTrho_int,XH1_int,XHe4_int,omega_int,mu_int,mu0_int, &
  HI_int,HeI_int,HeII_int
real(kind=kindreal), dimension(Size_Env), save:: logr_env,m_env,logT_env,logrho_env,logPturb_env,logP_env,alpha_env, &
  delta_env,nablaad_env,Gamma1_env,cp_env,cv_env,dlnPdlnrho_env,dlnPdlnT_env,nablae_env,logkappa_env,Lrad_env,L_env, &
  dlnkappadlnrhoT_env,dlnkappadlnTrho_env,epsilon_env,dlnepsilondlnrhoT_env,dlnepsilondlnTrho_env,XH1_env,XHe4_env, &
  omega_env,dlnMdlnVar_env,dlnrdlnVar_env,mu_env,mu0_env,v_MLT_env,time_TO_env,HI_env,HeI_env,HeII_env
real(kind=kindreal), dimension(Size_Atmos), save:: logr_atm,m_atm,logT_atm,logrho_atm,logP_atm,alpha_atm,delta_atm,nablaad_atm, &
  Gamma1_atm,cp_atm,cv_atm,dlnPdlnrho_atm,dlnPdlnT_atm,nablae_atm,logkappa_atm,Lrad_atm,L_atm,dlnkappadlnrhoT_atm, &
  dlnkappadlnTrho_atm,epsilon_atm,dlnepsilondlnrhoT_atm,dlnepsilondlnTrho_atm,XH1_atm,XHe4_atm,omega_atm,tau_atm,mu_atm,mu0_atm, &
  HI_atm,HeI_atm,HeII_atm

character(16), dimension(Size_Env), save:: Var_env

public:: StoreStructure_int
public:: StoreStructure_env
public:: StoreStructure_atm
public:: PrintCompleteStructure

contains

!******************************************************************************
subroutine StoreStructure_int(layer,logr,Mr,logT,logrho,logP,alpha,delta,nablaad,nablarad,logkappa,L, &
                              dlnkappadlnPT,dlnkappadlnTP,epsi,dlnepsilondlnPT,dlnepsilondlnTP, &
                              XH1,XHe4,omega,mu,mu0)
!-----------------------------------------------------------------------
   implicit none

   integer, intent(in):: layer

   real(kind=kindreal), intent(in):: logr,Mr,logT,logrho,logP,alpha,delta,nablaad,nablarad,logkappa,L, &
                                     dlnkappadlnPT,dlnkappadlnTP,epsi,dlnepsilondlnPT,dlnepsilondlnTP, &
                                     XH1,XHe4,omega,mu,mu0
!-----------------------------------------------------------------------
   logr_int(layer) = logr
   m_int(layer) = Mr*Msol
   logT_int(layer) = logT
   logrho_int(layer) = logrho
   logP_int(layer) = logP
   alpha_int(layer) = alpha
   delta_int(layer) = delta
   nablaad_int(layer) = nablaad
   Gamma1_int(layer) = 1.d0/(alpha-delta*nablaad)
   mu_int(layer) = mu
   mu0_int(layer) = mu0
   cp_int(layer) = 10.d0**(logP-logrho-logT)*delta/nablaad*mu0*cst_u/cst_k
   cv_int(layer) = cp_int(layer)/(alpha*Gamma1_int(layer))
   dlnPdlnrho_int(layer) = 1.d0/alpha
   dlnPdlnT_int(layer) = delta/alpha
   if (nablarad >= nablaad) then
     nablae_int(layer) = nablaad
   else
     nablae_int(layer) = nablarad
   endif
   logkappa_int(layer) = logkappa
   Lrad_int(layer) = 16.d0*pi*cst_a*cst_c*cst_G/3.d0 * Mr*Msol*10.d0**(4.d0*logT-logkappa-logP) * &
                     nablae_int(layer)
   L_int(layer) = L*Lsol
   dlnkappadlnrhoT_int(layer) = dlnkappadlnPT/alpha
   dlnkappadlnTrho_int(layer) = delta*dlnkappadlnPT/alpha + dlnkappadlnTP
   epsilon_int(layer) = epsi
   dlnepsilondlnrhoT_int(layer) = dlnepsilondlnPT/alpha
   dlnepsilondlnTrho_int(layer) = delta*dlnepsilondlnPT/alpha + dlnepsilondlnTP
   XH1_int(layer) = XH1
   XHe4_int(layer) = XHe4
   omega_int(layer) = omega
   HI_int(layer) = 1.0d0
   HeI_int(layer) = 0.0d0
   HeII_int(layer) = 1.0d0
   Int_Layers = layer

   return

end subroutine StoreStructure_int
!******************************************************************************

!******************************************************************************
subroutine StoreStructure_env(layer,logr,Mr,logT,logrho,logPturb,logP,alpha,delta,nablaad,gamma1,cp,nablae, &
                              logkappa,dlnkappadlnPT,dlnkappadlnTP,dlnMdlnVar,dlnrdlnVar,mu,mu0, &
                              v_MLT,time_TO,xion,Var)
!-----------------------------------------------------------------------
   implicit none

   integer, intent(in):: layer

   real(kind=kindreal), intent(in):: logr,Mr,logT,logrho,logPturb,logP,alpha,delta,nablaad,gamma1,cp,nablae, &
                                     logkappa,dlnkappadlnPT,dlnkappadlnTP,dlnMdlnVar,dlnrdlnVar,mu,mu0,v_MLT,time_TO
   real(kind=kindreal),dimension(3),intent(in):: xion
   character(16), intent(in):: Var
!-----------------------------------------------------------------------
   logr_env(layer) = logr
   m_env(layer) = 10.d0**Mr
   logT_env(layer) = logT
   logrho_env(layer) = logrho
   logP_env(layer) = logP
   logPturb_env(layer) = logPturb
   alpha_env(layer) = alpha
   delta_env(layer) = delta
   nablaad_env(layer) = nablaad
   Gamma1_env(layer) = gamma1
   cp_env(layer) = cp
   cv_env(layer) = cp/(alpha*gamma1)
   dlnPdlnrho_env(layer) = 1.d0/alpha
   dlnPdlnT_env(layer) = delta/alpha
   nablae_env(layer) = nablae
   logkappa_env(layer) = logkappa
   Lrad_env(layer) = 16.d0*pi*cst_a*cst_c*cst_G/3.d0 *10.d0**(Mr+4.d0*logT-logkappa-logP) * &
                     nablae_env(layer)
   L_env(layer) = -1.d-20
   dlnkappadlnrhoT_env(layer) = dlnkappadlnPT/alpha
   dlnkappadlnTrho_env(layer) = delta*dlnkappadlnPT/alpha + dlnkappadlnTP
   epsilon_env(layer) = 1.d-32
   dlnepsilondlnrhoT_env(layer) = 0.d0
   dlnepsilondlnTrho_env(layer) = 0.d0
   XH1_env(layer) = -1.d20
   XHe4_env(layer) = -1.d20
   omega_env(layer) = -1.d20
   dlnMdlnVar_env(layer) = dlnMdlnVar
   dlnrdlnVar_env(layer) = dlnrdlnVar
   mu_env(layer) = mu
   mu0_env(layer) = mu0
   v_MLT_env(layer) = v_MLT
   time_TO_env(layer) = time_TO
   HI_env(layer) = xion(1)
   HeI_env(layer) = xion(2)
   HeII_env(layer) = xion(3)
   Var_env(layer) = Var
   Env_Layers = layer
   return

end subroutine StoreStructure_env
!******************************************************************************

!******************************************************************************
subroutine StoreStructure_atm(layer,logT,logrho,logP,alpha,delta,nablaad,cp,logkappa,dlnkappadlnPT, &
                              dlnkappadlnTP,tau,mu,mu0,xion)
!-----------------------------------------------------------------------
   implicit none

   integer, intent(in):: layer

   real(kind=kindreal), intent(in)::logT,logrho,logP,alpha,delta,nablaad,cp,logkappa,dlnkappadlnPT, &
                                    dlnkappadlnTP,tau,mu,mu0
   real(kind=kindreal),dimension(3),intent(in):: xion
!-----------------------------------------------------------------------
   logr_atm(layer) = -1.d-20
   m_atm(layer) = -1.d-20
   logT_atm(layer) = logT
   logrho_atm(layer) = logrho
   logP_atm(layer) = logP
   alpha_atm(layer) = alpha
   delta_atm(layer) = delta
   nablaad_atm(layer) = nablaad
   Gamma1_atm(layer) = 1.d0/(alpha-delta*nablaad)
   cp_atm(layer) = cp
   cv_atm(layer) = cp/(alpha*Gamma1_atm(layer))
   dlnPdlnrho_atm(layer) = 1.d0/alpha
   dlnPdlnT_atm(layer) = delta/alpha
   nablae_atm(layer) = -1.d-20
   logkappa_atm(layer) = logkappa
   Lrad_atm(layer) = -1.d-20
   L_atm(layer) = -1.d-20
   dlnkappadlnrhoT_atm(layer) = dlnkappadlnPT/alpha
   dlnkappadlnTrho_atm(layer) = delta*dlnkappadlnPT/alpha + dlnkappadlnTP
   epsilon_atm(layer) = 1.d-32
   dlnepsilondlnrhoT_atm(layer) = 0.d0
   dlnepsilondlnTrho_atm(layer) = 0.d0
   XH1_atm(layer) = -1.d-20
   XHe4_atm(layer) = -1.d-20
   omega_atm(layer) = -1.d-20
   tau_atm(layer) = tau
   mu_atm(layer) = mu
   mu0_atm(layer) = mu0
   HI_atm(layer) = xion(1)
   HeI_atm(layer) = xion(2)
   HeII_atm(layer) = xion(3)
   Atm_Layers = layer

   return

end subroutine StoreStructure_atm
!******************************************************************************

!******************************************************************************
subroutine PrintCompleteStructure
!-----------------------------------------------------------------------
   use caramodele, only: nwmd

   implicit none

   integer:: i

   real(kind=kindreal):: DeltaVar,kappa_mean,rho_mean,r_mean,nablarad,R_phot

   character(*), parameter:: &
   Title_Line = "   n            log(r)             M_int            log(T)          log(rho)          "// &
                "  log(P)                Cv (dln P/dln rho)_T (dln P/dln T)_rho           nabla_e      "// &
                "    nabla_ad             L_rad             L_tot        log(kappa) (dln k/dln rho)_T "// &
                "(dln k/dln T)_rho           epsilon (dln E/dln rho)_T (dln E/dln T)_rho              "// &
                "X_H1             X_He4                mu               mu0             Omega            P_turb"// &
                "             V_MLT     time_TurnOver       HII      HeII     HeIII"
!-----------------------------------------------------------------------
   ! Remove the useless line in the atmosphere:
   logr_atm(Atm_Layers-1) = logr_atm(Atm_Layers)
   m_atm(Atm_Layers-1) = m_atm(Atm_Layers)
   logT_atm(Atm_Layers-1) = logT_atm(Atm_Layers)
   logrho_atm(Atm_Layers-1) = logrho_atm(Atm_Layers)
   logP_atm(Atm_Layers-1) = logP_atm(Atm_Layers)
   alpha_atm(Atm_Layers-1) = alpha_atm(Atm_Layers)
   delta_atm(Atm_Layers-1) = delta_atm(Atm_Layers)
   nablaad_atm(Atm_Layers-1) = nablaad_atm(Atm_Layers)
   Gamma1_atm(Atm_Layers-1) = Gamma1_atm(Atm_Layers)
   cp_atm(Atm_Layers-1) = cp_atm(Atm_Layers)
   cv_atm(Atm_Layers-1) = cv_atm(Atm_Layers)
   dlnPdlnrho_atm(Atm_Layers-1) = dlnPdlnrho_atm(Atm_Layers)
   dlnPdlnT_atm(Atm_Layers-1) = dlnPdlnT_atm(Atm_Layers)
   nablae_atm(Atm_Layers-1) = nablae_atm(Atm_Layers)
   logkappa_atm(Atm_Layers-1) = logkappa_atm(Atm_Layers)
   Lrad_atm(Atm_Layers-1) = Lrad_atm(Atm_Layers)
   L_atm(Atm_Layers-1) = L_atm(Atm_Layers)
   dlnkappadlnrhoT_atm(Atm_Layers-1) = dlnkappadlnrhoT_atm(Atm_Layers)
   dlnkappadlnTrho_atm(Atm_Layers-1) = dlnkappadlnTrho_atm(Atm_Layers)
   epsilon_atm(Atm_Layers-1) = epsilon_atm(Atm_Layers)
   dlnepsilondlnrhoT_atm(Atm_Layers-1) = dlnepsilondlnrhoT_atm(Atm_Layers)
   dlnepsilondlnTrho_atm(Atm_Layers-1) = dlnepsilondlnTrho_atm(Atm_Layers)
   XH1_atm(Atm_Layers-1) = XH1_atm(Atm_Layers)
   XHe4_atm(Atm_Layers-1) = XHe4_atm(Atm_Layers)
   omega_atm(Atm_Layers-1) = omega_atm(Atm_Layers)
   tau_atm(Atm_Layers-1) = tau_atm(Atm_Layers)
   mu_atm(Atm_Layers-1) = mu_atm(Atm_Layers)
   mu0_atm(Atm_Layers-1) = mu0_atm(Atm_Layers)
   HI_atm(Atm_Layers-1) = HI_atm(Atm_Layers)
   HeI_atm(Atm_Layers-1) = HeI_atm(Atm_Layers)
   HeII_atm(Atm_Layers-1) = HeII_atm(Atm_Layers)
   Atm_Layers = Atm_Layers - 1

   ! remove the useless line in the envelope:

   do i=1,Env_Layers-1
      logr_env(i) = logr_env(i+1)
      m_env(i) = m_env(i+1)
      logT_env(i) = logT_env(i+1)
      logrho_env(i) = logrho_env(i+1)
      logP_env(i) = logP_env(i+1)
      logPturb_env(i) = logPturb_env(i+1)
      alpha_env(i) = alpha_env(i+1)
      delta_env(i) = delta_env(i+1)
      nablaad_env(i) = nablaad_env(i+1)
      Gamma1_env(i) = Gamma1_env(i+1)
      cp_env(i) = cp_env(i+1)
      cv_env(i) = cv_env(i+1)
      dlnPdlnrho_env(i) = dlnPdlnrho_env(i+1)
      dlnPdlnT_env(i) = dlnPdlnT_env(i+1)
      nablae_env(i) = nablae_env(i+1)
      logkappa_env(i) = logkappa_env(i+1)
      Lrad_env(i) = Lrad_env(i+1)
      L_env(i) = L_env(i+1)
      dlnkappadlnrhoT_env(i) = dlnkappadlnrhoT_env(i+1)
      dlnkappadlnTrho_env(i) = dlnkappadlnTrho_env(i+1)
      epsilon_env(i) = epsilon_env(i+1)
      dlnepsilondlnrhoT_env(i) = dlnepsilondlnrhoT_env(i+1)
      dlnepsilondlnTrho_env(i) = dlnepsilondlnTrho_env(i+1)
      XH1_env(i) = XH1_env(i+1)
      XHe4_env(i) = XHe4_env(i+1)
      omega_env(i) = omega_env(i+1)
      dlnMdlnVar_env(i) = dlnMdlnVar_env(i+1)
      dlnrdlnVar_env(i) = dlnrdlnVar_env(i+1)
      mu_env(i) = mu_env(i+1)
      mu0_env(i) = mu0_env(i+1)
      v_MLT_env(i) = v_MLT_env(i+1)
      time_TO_env(i) = time_TO_env(i+1)
      HI_env(i) = HI_env(i+1)
      HeI_env(i) = HeI_env(i+1)
      HeII_env(i) = HeII_env(i+1)
      Var_env(i) = Var_env(i+1)
   enddo
   Env_Layers = Env_Layers-1
   ! Removal of the envelope layer with log(r) smaller than the first internal layer.
   do while (logr_env(Env_Layers) < logr_int(1))
     Env_Layers = Env_Layers-1
   enddo

   ! in the envelope and atmosphere, store the data similar to the first layer of the interior:
   do i=1,Env_Layers
      L_env(i) = L_int(1)
      XH1_env(i) = XH1_int(1)
      XHe4_env(i) = XHe4_int(1)
      omega_env(i) = omega_int(1)
   enddo
   do i=1,Atm_Layers
      L_atm(i) = L_int(1)
      XH1_atm(i) = XH1_int(1)
      XHe4_atm(i) = XHe4_int(1)
      omega_atm(i) = omega_int(1)
   enddo

   ! Compute the current mass and radius in the second layer of the envelope :
   if (Var_env(1) == "logrho") then
      DeltaVar = logrho_env(2) - logrho_env(1)
   else if (Var_env(1) == "logP") then
      DeltaVar = logP_env(2) - logP_env(1)
   else
      write(*,*) 'Variable problem in PrintCompleteStructure, in module PrintAll.'
      stop
   endif
   logr_env(2) = logr_env(1) + dlnrdlnVar_env(1)*DeltaVar
   m_env(2) = 10.d0**(log10(m_env(1)) + dlnMdlnVar_env(1)*DeltaVar)
   Lrad_env(2) = 16.d0*pi*cst_a*cst_c*cst_G/3.d0 * m_env(2)*10.d0**(4.d0*logT_env(2)-logkappa_env(2) &
                     -logP_env(2)) * nablae_env(2)

   ! Compute the missing data in atmosphere :
   R_phot = sqrt(Lum_save*Lsol/(4.d0*pi*cst_sigma))/Teff_save**2.d0
   logr_atm(Atm_Layers) = log10(R_phot)
   m_atm(Atm_Layers) = m_env(1)
   do i=Atm_Layers-1,1,-1
      kappa_mean = (10.d0**(logkappa_atm(i))+10.d0**(logkappa_atm(i+1)))/2.d0
      rho_mean = (10.d0**(logrho_atm(i))+10.d0**(logrho_atm(i+1)))/2.d0
      logr_atm(i) = log10(10.d0**logr_atm(i+1) + (tau_atm(i+1) - tau_atm(i))/(kappa_mean*rho_mean))
      r_mean = (10.d0**(logr_atm(i))+10.d0**(logr_atm(i+1)))/2.d0
      m_atm(i) = m_atm(i+1)-4.d0*pi*r_mean**2.d0*rho_mean*(10.d0**logr_atm(i+1)-10.d0**logr_atm(i))
   enddo

   ! Computation of the other atmospheric data :
   do i = 1,Atm_Layers
      nablarad = 3.d0*L_atm(i)*10.d0**(logkappa_atm(i)+logP_atm(i)-4.d0*logT_atm(i))/ &
                (16.d0*pi*cst_a*cst_c*cst_G*m_atm(i))
      if (nablarad >= nablaad_atm(i)) then
        nablae_atm(i) = nablaad_atm(i)
      else
        nablae_atm(i) = nablarad
      endif
      Lrad_atm(i) = 16.d0*pi*cst_a*cst_c*cst_G/3.d0 * m_atm(i)*10.d0**(4.d0*logT_atm(i)-logkappa_atm(i)- &
                    logP_atm(i)) * nablae_atm(i)
   enddo

   ! Printing of the data
   write(File_Unit,'("Model num   : ",i6)'),nwmd
   write(File_Unit,'("Time [yr]   : ",es16.9)'),time_save
   write(File_Unit,'("Mass [Msun] : ",es16.9)'),mass_save
   write(File_Unit,'("Radius [cm] : ",es16.9)'),R_phot
   write(File_Unit,'("log(L/Lsun) : ",es16.9)'),log10(Lum_save)
   write(File_Unit,'("log(Teff/K) : ",es16.9)'),log10(Teff_save)
   write(File_Unit,'("C12_surf    : ",es16.9)'),C12_save
   write(File_Unit,'("C13_surf    : ",es16.9)'),C13_save
   write(File_Unit,'("N14_surf    : ",es16.9)'),N14_save
   write(File_Unit,'("O16_surf    : ",es16.9)'),O16_save
   write(File_Unit,*)
   write(File_Unit,'(a)') trim(Title_Line)
   do i=Int_Layers,1,-1
      write(File_Unit,'(i4,26(2x,es16.9),3(2x,f8.6))')Int_Layers-i+1,logr_int(i),m_int(i),logT_int(i), &
                                           logrho_int(i),logP_int(i),cv_int(i),dlnPdlnrho_int(i), &
                                           dlnPdlnT_int(i),nablae_int(i),nablaad_int(i),Lrad_int(i), &
                                           L_int(i),logkappa_int(i),dlnkappadlnrhoT_int(i), &
                                           dlnkappadlnTrho_int(i),epsilon_int(i),dlnepsilondlnrhoT_int(i), &
                                           dlnepsilondlnTrho_int(i),XH1_int(i),XHe4_int(i),mu_int(i), &
                                           mu0_int(i),omega_int(i),0.0d0,0.0d0,0.0d0,HI_int(i),HeI_int(i),&
                                           HeII_int(i)
   enddo
   do i=Env_Layers,1,-1
      write(File_Unit,'(i4,26(2x,es16.9),3(2x,f8.6))')Env_Layers-i+1+Int_Layers,logr_env(i),m_env(i),logT_env(i), &
                                           logrho_env(i),logP_env(i),cv_env(i),dlnPdlnrho_env(i), &
                                           dlnPdlnT_env(i),nablae_env(i),nablaad_env(i),Lrad_env(i), &
                                           L_env(i),logkappa_env(i),dlnkappadlnrhoT_env(i), &
                                           dlnkappadlnTrho_env(i),epsilon_env(i),dlnepsilondlnrhoT_env(i), &
                                           dlnepsilondlnTrho_env(i),XH1_env(i),XHe4_env(i),mu_env(i), &
                                           mu0_env(i),omega_env(i),logPturb_env(i),v_MLT_env(i),time_TO_env(i), &
                                           HI_env(i),HeI_env(i),HeII_env(i)
   enddo
   do i=Atm_Layers-1,1,-1
      write(File_Unit,'(i4,26(2x,es16.9),3(2x,f8.6))')Atm_Layers-i+Int_Layers+Env_Layers,logr_atm(i),m_atm(i), &
                                           logT_atm(i),logrho_atm(i),logP_atm(i),cv_atm(i), &
                                           dlnPdlnrho_atm(i),dlnPdlnT_atm(i),nablae_atm(i), &
                                           nablaad_atm(i),Lrad_atm(i),L_atm(i),logkappa_atm(i), &
                                           dlnkappadlnrhoT_atm(i),dlnkappadlnTrho_atm(i),epsilon_atm(i), &
                                           dlnepsilondlnrhoT_atm(i),dlnepsilondlnTrho_atm(i),XH1_atm(i), &
                                           XHe4_atm(i),mu_atm(i),mu0_atm(i),omega_atm(i),0.0d0,0.0d0,0.0d0, &
                                           HI_atm(i),HeI_atm(i),HeII_atm(i)
   enddo
   write(File_Unit,*)

   return

end subroutine PrintCompleteStructure
!***********************************************************************
end module PrintAll
