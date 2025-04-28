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
  use inputparam,only: binm2,periodini,verbose,eccentricity_ini,ie2_prescription,posyd_prescription,&
  twin_system

  implicit none

  real(kindreal):: ab,romorb,romini,rstar,period,eccentricity


private
public:: period,dLtidcalc,eccentricity,compute_k2_from_structure

contains
!======================================================================
subroutine remove_duplicate_elements(table1, table2, table3, subtable1, subtable2, subtable3)
!----------------------------------------------------------------------
  ! This subroutine removes:
  ! 1. Duplicate or unordered values in table1 (with corresponding values in table2 and table3)
  ! 2. Duplicate or unordered values in table3 (with corresponding values in table1 and table2)

  implicit none

  ! Input arrays
  real(kind=kindreal), intent(in) :: table1(:), table2(:), table3(:)
  ! Output arrays
  real(kind=kindreal), allocatable, intent(out) :: subtable1(:), subtable2(:), subtable3(:)

  integer :: i, j, n
  logical, allocatable :: is_valid(:)
  real(kind=kindreal) :: last_valid_value1, last_valid_value3

!----------------------------------------------------------------------
  ! Get the size of the input arrays
  n = size(table1)

  ! Allocate a logical array to track valid elements
  allocate(is_valid(n))
  is_valid = .true.  ! Initialize all elements as valid

  ! Track the last valid value
  last_valid_value1 = table1(1)
  last_valid_value3 = table3(1)

  ! Check for duplicates and ordering issues
  do i = 2, n
    if (table1(i) <= last_valid_value1 .or. table3(i) <= last_valid_value3) then
      is_valid(i) = .false.  ! Mark duplicates and out-of-order elements as invalid
    else
      last_valid_value1 = table1(i)  ! Update last valid value for table1
      last_valid_value3 = table3(i)  ! Update last valid value for table3
    end if
  end do

  ! Count the number of valid elements
  j = count(is_valid)

  ! Allocate the output arrays based on the count of valid elements
  allocate(subtable1(j), subtable2(j), subtable3(j))

  ! Populate the output arrays with valid elements 
  j = 0
  do i = 1, n
    if (is_valid(i)) then
      j = j + 1
      subtable1(j) = table1(i)
      subtable2(j) = table2(i)
      subtable3(j) = table3(i)
    end if
  end do

  ! Deallocate temporary arrays
  deallocate(is_valid)

end subroutine remove_duplicate_elements
!======================================================================
subroutine cubic_Lagrangian_polynomial(x_interpol, y_interpol, x_i1, x_i2, x_i3, x_i4, y_i1, y_i2, y_i3, y_i4)
!----------------------------------------------------------------------
  ! Subroutine that performs a cubic interpolation based on four points:
  ! (x_i1, y_i1), (x_i2, y_i2), (x_i3, y_i3), (x_i4, y_i4), with x_interpol as input and
  ! y_interpol as output.
  
  implicit none
  
  real(kindreal),intent(in)::  x_interpol,x_i1,x_i2,x_i3,x_i4,y_i1,y_i2,y_i3,y_i4
  real(kindreal),intent(out):: y_interpol
  
  real(kindreal):: lag_1, lag_2, lag_3, lag_4
!----------------------------------------------------------------------

  lag_1 = (x_interpol-x_i2)*(x_interpol-x_i3)*(x_interpol-x_i4)/((x_i1-x_i2)*(x_i1-x_i3)*(x_i1-x_i4))
  lag_2 = (x_interpol-x_i1)*(x_interpol-x_i3)*(x_interpol-x_i4)/((x_i2-x_i1)*(x_i2-x_i3)*(x_i2-x_i4))
  lag_3 = (x_interpol-x_i1)*(x_interpol-x_i2)*(x_interpol-x_i4)/((x_i3-x_i1)*(x_i3-x_i2)*(x_i3-x_i4))
  lag_4 = (x_interpol-x_i1)*(x_interpol-x_i2)*(x_interpol-x_i3)/((x_i4-x_i1)*(x_i4-x_i2)*(x_i4-x_i3))
  
  y_interpol = y_i1*lag_1 + y_i2*lag_2 + y_i3*lag_3 + y_i4*lag_4
  
end subroutine cubic_Lagrangian_polynomial
!======================================================================
subroutine function_eta_radius_rho_ratio(f_eta_r,eta,radius,rho_ratio)
!----------------------------------------------------------------------
  ! Subroutine that retrieves f(eta,r,rho_ratio(r)) of the Clairaut-Radau equation,
  ! satisfying the equation deta/dr = f(eta,r,rho_ratio(r)).
  
  implicit none
  

  real(kindreal),intent(in)::  eta,radius,rho_ratio
  real(kindreal),intent(out)::  f_eta_r
!----------------------------------------------------------------------

  f_eta_r = (6.d0 * (1.d0 - rho_ratio*(eta + 1.d0)) + eta * (1.d0 - eta))/radius
  
end subroutine function_eta_radius_rho_ratio
!======================================================================
subroutine RK4_solver(y_array,x_array,g_x_array)
!----------------------------------------------------------------------
! Subroutine that solves an equation of the form dy/dx = f(y, x, g(x)) using a fourth order
! Runge-Kutta method.
! 
! In our case, g(x) corresponds to the ratio of densities rho/rho_bar. We don't have an explicit
! form for g(x), but only discrete values. 
! In the Runge-Kutta approach, it is necessary to evalue f at some intermediate values x_i-1/2.
! Since we don't have an explicit form for g(x) we need to extrapolate the points. A simple
! linear interpolation bewtween x_i-1 and x_i would break down the 4th order convergence of the
! integrator. For this reason, we use Lagrangian cubic polynomials, relying on 4 points, which
! preserve the 4th order of convergence (actually this is true for equally spaced interpolation
! nodes, but this is the best we can do) : https://en.wikipedia.org/wiki/Polynomial_interpolation.

  implicit none
  
  real(kindreal),intent(in)::  x_array(:),g_x_array(:)
  real(kindreal),intent(out):: y_array(:)
  
  integer:: i, n_tot
  
  real(kindreal):: x_half,g_x_half,f_y_x,y_temp1_RK,y_temp2_RK,y_temp3_RK,k1_RK,k2_RK,k3_RK,k4_RK
!----------------------------------------------------------------------

  n_tot = size(x_array)

  do i = 1, n_tot
    if (i == 1) then
      ! center: boundary condition: y(0) = 0
      y_array(1) = 0.d0
    else
      ! interior -> solve equation
      
      ! start by building the intermediate points
      ! x_half and g_x_half, necessary for the 4th order Runge Kutta scheme
      
      ! middle point of interval x_i-1, x_i
      x_half = (x_array(i-1) + x_array(i))/2.d0
    
      ! cubic interpolation to find the corresponding value of g(x_half)
      if (i == 2) then
        ! second layer, take layers 1, 2, 3, 4 for cubic interpolation

        call cubic_Lagrangian_polynomial(x_half, g_x_half, x_array(1), x_array(2), x_array(3), x_array(4),&
        g_x_array(1), g_x_array(2),g_x_array(3),g_x_array(4))

      else if (i == n_tot) then
        ! last layer, take layers n_tot - 3, n_tot - 2, n_tot - 1, n_tot

        call cubic_Lagrangian_polynomial(x_half, g_x_half, x_array(n_tot-3), x_array(n_tot-2), x_array(n_tot-1),&
        x_array(n_tot), g_x_array(n_tot-3), g_x_array(n_tot-2),g_x_array(n_tot-1),g_x_array(n_tot))
      else
        ! intermediate layers, take layers i - 2, i - 1, i, i + 1
        ! ok since smallest i in this case is i = 3 -> take 1, 2, 3, 4 and largest i is i = n_tot - 1
        ! -> n_tot - 3, n_tot - 2, n_tot - 1, n_tot
        
        call cubic_Lagrangian_polynomial(x_half, g_x_half, x_array(i-2), x_array(i-1), x_array(i), x_array(i+1),&
        g_x_array(i-2), g_x_array(i-1),g_x_array(i),g_x_array(i+1))
      end if
      
      ! construction of k1, k2, k3 and k4 of the Runge Kutta method
      ! https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
      
      call function_eta_radius_rho_ratio(f_y_x,y_array(i-1),x_array(i-1),g_x_array(i-1))
      k1_RK = (x_array(i) - x_array(i-1)) * f_y_x
      
      y_temp1_RK = y_array(i-1) + 0.5d0*k1_RK
      
      call function_eta_radius_rho_ratio(f_y_x,y_temp1_RK,x_half,g_x_half)
      k2_RK = (x_array(i) - x_array(i-1)) * f_y_x
      
      y_temp2_RK = y_array(i-1) + 0.5d0*k2_RK
      
      call function_eta_radius_rho_ratio(f_y_x,y_temp2_RK,x_half,g_x_half)
      k3_RK = (x_array(i) - x_array(i-1)) * f_y_x
      
      y_temp3_RK = y_array(i-1) + k3_RK
      
      call function_eta_radius_rho_ratio(f_y_x,y_temp3_RK,x_array(i),g_x_array(i))
      k4_RK = (x_array(i) - x_array(i-1)) * f_y_x
      
      y_array(i) = y_array(i-1) + (k1_RK + 2.d0*k2_RK + 2.d0*k3_RK + k4_RK)/6.d0
      
    end if
  end do
  
end subroutine RK4_solver
!======================================================================
subroutine compute_k2_from_profiles(k2_AMC,radius_profile,rho_profile,mass_profile,rstar)
!----------------------------------------------------------------------
  ! Subroutine that retrieves the k2, when the radius (r/R), density (in g.cm-3)
  ! and mass (in g) profiles are given as input by solving the Clairaut-Radau equation.

  implicit none
  
  integer :: i, i_rmax

  real(kindreal),intent(in)::  radius_profile(:),rho_profile(:),mass_profile(:),rstar
  real(kindreal),intent(out)::  k2_AMC
  
  ! rho_ratio corresponds to rho/rho_bar, where rho_bar is the mean density of the sphere of radius r.

  real(kindreal):: rho_ratio_profile(size(radius_profile)), eta_profile(size(radius_profile))
  real(kindreal):: rho_bar, eta_surface
!----------------------------------------------------------------------

  i_rmax = 1 ! just to avoid segmentation fault in the first model
  do i = 1, size(rho_ratio_profile)
    if (mass_profile(i) < epsilon(0.d0)) then
      ! central layer (sometimes more than 1): mass = 0, so rho_bar set equal to rho -> ratio = 1
      rho_ratio_profile(i) = 1.d0
    else
      rho_bar = 3.d0 * mass_profile(i) / (4.d0 * pi * (radius_profile(i)*rstar)**3.d0)
      rho_ratio_profile(i) = rho_profile(i) / rho_bar
    end if
  end do

  ! The integration over the star is performed over the whole interior and envelope, i.e. up to
  ! tau(2/3) (base of the atmosphere).
  
  call RK4_solver(eta_profile,radius_profile,rho_ratio_profile)
  
  i_rmax = size(eta_profile)
  eta_surface = eta_profile(i_rmax)
  k2_AMC = (3.d0 - eta_surface)/(4.d0 + 2.d0 * eta_surface)
  
end subroutine compute_k2_from_profiles

!======================================================================
subroutine compute_k2_from_structure(k2_AMC)
!----------------------------------------------------------------------
  ! Subroutine to compute the apsidal motion constant, relevant for precession by distortion and
  ! for the equilibrium tides. In this subroutine the required profiles (radius and density) are
  ! built joining the interior and the envelope. Then the subroutine compute_k2_from_profiles is
  ! called with the built profiles.
  
  use caramodele, only: gls,glm,teff
  use strucmod, only: r,q,rho,envel

  implicit none
  
  integer :: i, i_bot_env, i_center, nsize_k2

  real(kindreal),intent(out)::  k2_AMC
  real(kindreal), allocatable :: radius_for_k2(:), rho_for_k2(:), mass_for_k2(:), rad_unique(:), rho_unique(:), mass_unique(:)
!----------------------------------------------------------------------
  rstar=sqrt(gls*Lsol/(4.d0*pi*cst_sigma))/(teff**2.d0)
  
  ! We start by building the radius and density profiles joining the interior and envelope   
  do i = 1, size(envel, 1)
    ! Envel(i,3) contains the different radii of the envelope.
    if (10.d0 ** envel(i,3) < exp(r(1))) then
      ! Retrieves the index of the bottom of the envelope. Going from surface to center, the first radius to satisfy the if 
      ! statement is considered to be the bottom of the envelope. The envelope may contain a few layer with radii
      ! smaller than the extend of the interior. These layers are ignored for the computation of k2, consistently with 
      ! these layers being removed when creating the strucdata files.
      i_bot_env = i - 1
      exit
    end if
    if (i > 1) then
      ! In some instances, it can happen that the deepest layers of the envelope are not correctly ordered
      ! in this case, we define as the bottom of the envelope the first layer for which the next layer has
      ! a higher radius. Typically, this only occurs very near the bottom, and we can smoothly transition to
      ! the interior there.
      if (10.d0 ** envel(i,3) > 10.d0 ** envel(i-1,3)) then
         i_bot_env = i - 1
         exit
      end if
    end if
  end do
  
  ! index of the center of the star
  do i = 1, size(r)
    if (r(i) < epsilon(0.0d0)) then
      i_center = i - 1
      exit
    end if
  end do
  
  ! Size of the tables that contain radius and other quantities from center to surface
  ! Number of layers of envelope + number of layers of interior
  nsize_k2 = i_bot_env + i_center

  ! Tables that contain the relevant quantities
  allocate(radius_for_k2(nsize_k2))
  allocate(rho_for_k2(nsize_k2))
  allocate(mass_for_k2(nsize_k2))
  
  do i = 1, nsize_k2
    if (i <= i_center) then
      ! Stores the quantities inside the interior
      radius_for_k2(i_center - i + 1) = exp(r(i))/rstar
      rho_for_k2(i_center - i + 1) = exp(rho(i))
      mass_for_k2(i_center - i + 1) = exp(glm)*(1.d0 - exp(q(i)))
    else
      ! Stores the quantities inside the envelope
      radius_for_k2(nsize_k2 - i + i_center + 1) = (10.d0**envel(i - i_center,3))/rstar
      rho_for_k2(nsize_k2 - i + i_center + 1) = 10.d0**envel(i - i_center,4)
      mass_for_k2(nsize_k2 - i + i_center + 1) = 10.d0**envel(i - i_center,5)
    end if
  end do
  
  ! Remove duplicate or unordered elements
  call remove_duplicate_elements(radius_for_k2, rho_for_k2, mass_for_k2, rad_unique, rho_unique, mass_unique)
  
  ! Call the routine to compute the k2 from the built profiles
  call compute_k2_from_profiles(k2_AMC,rad_unique, rho_unique, mass_unique,rstar)
 
 
end subroutine compute_k2_from_structure
!======================================================================
subroutine dLtidcalc(dLtid)
!----------------------------------------------------------------------
! Formalism taken from Zahn (1977), A&A 57, 383 and Sciarini et al. (2024), A&A 681, L1 for eccentric orbits
! The dynamical tides model consists of an expansion valid for low eccentricities, not recommended for very
! eccentric orbits, of say e >~ 0.5. 
! The set of equations in the case of eccentric orbits is the Eq. (9) of the letter.
!----------------------------------------------------------------------
  use inputparam,only: const_per,isol,include_dyn_tides,include_eq_tides
  use caramodele,only: gms,gls,teff,nwmd
  use rotmod,only: omegi,bmomit,vsuminenv
  use convection,only: r_core
  use strucmod,only: vnr,vna,r,envel,k2_AMC
  use timestep,only: alter,dzeit

  implicit none

  real(kindreal),intent(out):: dLtid
  real(kindreal):: spiper,regra,qr1,qr2,rltid,fact1,fact2,fact3,fact4,e2,tsyn,roint,&
       rroche,fact5,s10,s12,s22,s32,dltido,fas1,fas2,&
       difsav,difsav10,difsav12,difsav32,trot,dLtid_eq,f_conv,&
       f2_e,f5_e,P_tid,omega_eccen_term, min_tau_conv, mean_tau_conv, Mconv_fconv_ov_tau, Iconv_fconv_ov_tau_I,&
       temp_mass, temp_radius_squared
  integer:: i, j, top_conv_zone, count_convective_zones
  real(kindreal), allocatable :: min_tau_conv_zones(:), mean_tau_conv_zones(:), M_conv_zones(:), I_conv_zones(:),&
      tau_conv_zones_posyd(:)

!----------------------------------------------------------------------

  ! In this updated version we account for the equilibrium tides acting on small convective zones near the surface of the
  ! stars during the main sequence, in an approach very similar to that in Fragos et al. 2023ApJS..264...45F, except that 
  ! we obtain the convective turnover timescale within the MLT framework. 
  
  ! We identify small convective regions in the envelope and obtain the convective turnover timescale either taking the minimum
  ! values accross the layers of the region, or taking the geometric average. Both approaches provide a convective turnover timescale
  ! of the same order of magnitude, typically of the order of 10^4 seconds for a 15 Msun star.
   
  if (twin_system) then
    ! To mimic a twin system, we update the value of binm2 with the actual value of gms
    binm2 = gms
  endif
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
  
  ! Quantities relevant for the equilibrium tides
  
  count_convective_zones = 0
  do i=1,size(envel, 1)-1
    ! envel(i,2) is konv, i.e. tells if zone is convective
    if (envel(i,2) < epsilon(0.d0) .and. envel(i+1,2) > epsilon(0.d0)) then
      count_convective_zones = count_convective_zones + 1
    end if
  end do
  
  allocate(min_tau_conv_zones(count_convective_zones))
  allocate(mean_tau_conv_zones(count_convective_zones))
  allocate(M_conv_zones(count_convective_zones))
  allocate(I_conv_zones(count_convective_zones))
  allocate(tau_conv_zones_posyd(count_convective_zones))
  
  count_convective_zones = 0

  do i=1,size(envel, 1)-1
    if (envel(i,2) < epsilon(0.d0) .and. envel(i+1,2) > epsilon(0.d0)) then
      top_conv_zone = i + 1
      count_convective_zones = count_convective_zones + 1
    end if
    if (envel(i,2) > epsilon(0.d0) .and. envel(i+1,2) < epsilon(0.d0)) then
      do j = top_conv_zone, i
        if (j == top_conv_zone) then
          ! initialize tau_conv and I_conv_zone
          min_tau_conv = envel(j,20)
          mean_tau_conv = log10(envel(j,20))
          I_conv_zones(count_convective_zones) = 0.d0
        else
          ! minimum tau_conv of the region
          min_tau_conv = min(min_tau_conv,envel(j,20))
          ! geometric average of tau_conv of the region
          mean_tau_conv = mean_tau_conv + log10(envel(j,20))
        end if
        temp_mass = (10.d0**envel(j-1,5) - 10.d0**envel(j+1,5))/2
        temp_radius_squared = 10.d0**(2*envel(j,3))
        
        ! Obtain the moment of inertia of the convective shell by adding the contribution of each layer
        I_conv_zones(count_convective_zones) = I_conv_zones(count_convective_zones) + 2.d0*temp_mass*temp_radius_squared/3.d0
      end do
      
      ! Since the convective turnover varies logarithmically accros layers, we use the geometric 
      ! mean to obtain an average accros layers. Other option is to take the minimum (default)
      mean_tau_conv = mean_tau_conv / (i - top_conv_zone + 1)
      mean_tau_conv = 10.d0**(mean_tau_conv)
      
      ! In the default version, we actually take the minimum value of tau_conv accros the convective region
      min_tau_conv_zones(count_convective_zones) = min_tau_conv
      mean_tau_conv_zones(count_convective_zones) = mean_tau_conv

      ! Total mass of the convective region
      M_conv_zones(count_convective_zones) = 10.d0**envel(top_conv_zone,5) - 10.d0**envel(i,5)
      
      ! Posydon (Fragos 2023ApJS..264...45F) prescription: tau_conv calculated according to Eq. (7) 
      tau_conv_zones_posyd(count_convective_zones) = (M_conv_zones(count_convective_zones)*(10.d0**envel(top_conv_zone,3) -&
      10.d0**envel(i,3))*(10.d0**envel(top_conv_zone,3) + 10.d0**envel(i,3))/(6.d0*(gls*Lsol)))**(1.d0/3.d0)
      
    end if
  end do
  
  ! Tidal pumping period 1/P_tid = abs(1/P_orb - 1/P_spin)
  P_tid = 1.d0/(abs(1.d0/period - omegi(1)/(2.d0*pi)))
     
      
  Mconv_fconv_ov_tau = 0.d0
  Iconv_fconv_ov_tau_I = 0.d0

  do i=1,count_convective_zones
    ! f_conv to decrease the strength of the tides when the turnover timescale is greater than
    ! the pumping timescale (Golreich & Nicholson 1977Icar...30..301G).

    if (posyd_prescription) then
      ! f_conv consistent with the Posydon prescription
      f_conv = min(1.d0, (P_tid/(2.d0*tau_conv_zones_posyd(i)))**2.d0)
    else
      ! f_conv using tau_conv obtained from MLT (default)
      f_conv = min(1.d0, (P_tid/(2.d0*min_tau_conv_zones(i)))**2.d0)
    endif
    ! Add the contribution of all the intermediate convective zones. This can be seen as adding different contributions
    ! to the total torque. Each layer has its own tau_conv, f_conv and mass, and in this formalism all other terms in the torque
    ! are the same no matter the shells. So we add the contributions of all zones by adding I_conv_reg/I*f_conv/tau_conv, which 
    ! is equivalent to adding the torques of these regions. In practice, most of the time the total torque is largely dominated by the
    ! contribution of one of the convective zones.
    
    ! To be applied if the Posydon prescription is used: ratio of the mass M_conv.reg/M with f_conv_posyd consistent with
    ! tau_conv_posyd
    Mconv_fconv_ov_tau = Mconv_fconv_ov_tau + M_conv_zones(i)*f_conv/tau_conv_zones_posyd(i)
    
    ! Sciarini et al. 2025 (in prep.) prescription: ratio of the moment of inertia I_conv.reg/I with f_conv consistent with
    ! tau_conv obtained from the MLT. By default, we use the minimum tau_conv in the convective region.
    Iconv_fconv_ov_tau_I = Iconv_fconv_ov_tau_I + I_conv_zones(i)*f_conv/(min_tau_conv_zones(i)*(bmomit+vsuminenv))
  end do
  
  ! Deallocate temporary arrays
  deallocate(min_tau_conv_zones)
  deallocate(mean_tau_conv_zones)
  deallocate(M_conv_zones)
  deallocate(I_conv_zones)
  
  ! Functions fi(e^2), i = {2,3,4,5} in Hut 1981A&A....99..126H (Eq. 11)
  f2_e = 1.d0 + 15.d0/2.d0 * eccentricity**2.d0 + 45.d0/8.d0 * eccentricity**4.d0 + 5.d0/16.d0 * eccentricity**6.d0
  f5_e = 1.d0 + 3.d0 * eccentricity**2.d0 + 3.d0/8.d0 * eccentricity**4.d0
  
  ! Omega and eccentricity dependent term (different term for eccentricity derivative)
  omega_eccen_term = f2_e - f5_e*(omegi(1)/romorb)*(1.d0 - eccentricity**2.d0)**(1.5d0)
  
  ! Quantities relevant for the dynamical tides

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

  ! slm coefficients for the dynamical tides
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
      
      ! Dynamical tides with radiative damping
      if (include_dyn_tides) then
        dltid=1.5d0*fas1*e2*fas2*(sign(s22**(8.d0/3.d0),difsav)+eccentricity**2.d0*((1.d0/4.d0)*sign(s12**(8.d0/3.d0),difsav12)-&
        5.d0*sign(s22**(8.d0/3.d0),difsav)+(49.d0/4.d0)*sign(s32**(8.d0/3.d0),difsav32)))*dzeit
      else
        dltid = 0.d0
      end if
      
      ! Equilibrium tides with viscous friction, small convective shell, following Hut 81, Sciarini et al. 2025

      if (posyd_prescription) then
        ! Posydon prescription for equilibrium tides on small convective regions. Independent on k2. Ratio of mass, tau_conv
        ! according to Eq. (7)
        dLtid_eq = 2.d0/21.d0*3.d0*Mconv_fconv_ov_tau*rstar**2.d0*fas2*romorb*omega_eccen_term/&
        ((1.d0-eccentricity**2.d0)**6.d0)*dzeit
      else
        ! Default: Sciarini et al. 2025 (in prep.) prescription. Depends on k2, ratio of moment of inertia instead of ratio of
        ! masses, tau_conv self-consistent (obtained from MLT)
        dLtid_eq = k2_AMC* 3.d0*Iconv_fconv_ov_tau_I*(gms*Msol)*rstar**2.d0*fas2*romorb*omega_eccen_term/&
        ((1.d0-eccentricity**2.d0)**6.d0)*dzeit
      endif

      ! Add up the two contributions
      if (include_eq_tides) then
        dltid = dltid + dLtid_eq
      end if
      
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
    call orbitalevol(dLtid,dLtid_eq,Mconv_fconv_ov_tau,Iconv_fconv_ov_tau_I)
  endif
  if (rstar/rroche > 1.01d0 .and. isol == 0)then
    rewind(io_runfile)
    write(io_runfile,*) nwmd,'the star overfills the roche lobe'
      stop 'the star overfills the roche lobe'
  endif
  write(*,*)'########################################'
  write(*,*)'Roche lobe filling factor R/R_L = ',rstar/rroche
  write(*,*)'########################################'
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
subroutine orbitalevol(dLtid,dLtid_eq,Mconv_fconv_ov_tau,Iconv_fconv_ov_tau_I)
!----------------------------------------------------------------------
  !Module updated to correctly evolve the orbits in eccentric systems (Sciarini+24)
  use caramodele, only: gms,gls,dm_lost,nwmd,teff
  use rotmod,only: omegi
  use timestep,only: alter,dzeit
  use inputparam,only: include_dyn_tides,include_eq_tides
  use strucmod,only: k2_AMC
  
  use convection,only: r_core

  implicit none

  real(kindreal),intent(in)::  dLtid, dLtid_eq, Mconv_fconv_ov_tau,Iconv_fconv_ov_tau_I
  real(kindreal):: dab,qr1,qr2,qtot,qtot_inv,orbang,orstwi,rorstw,term1,term2,term3,term4,term_eccentricity,deccentricity,rltid,&
  s10,s12,s22,s32,difsav,difsav10,difsav12,difsav32,e2,fact2,deccentricity_eq,&
       f3_e,f4_e,omega_eccen_term
!----------------------------------------------------------------------
  rstar=sqrt(gls*Lsol/(4.d0*pi*cst_sigma))/(teff**2.d0)
  fact2=sqrt(cst_G*gms*Msol/rstar**3.d0)
  
  ! Quantities relevant for the equilibrium tides
  
  ! Functions fi(e^2), i = {2,3,4,5} in Hut 1981A&A....99..126H (Eq. 11)
  f3_e = 1.d0 + 15.d0/4.d0 * eccentricity**2.d0 + 15.d0/8.d0 * eccentricity**4.d0 + 5.d0/64.d0 * eccentricity**6.d0
  f4_e = 1.d0 + 3.d0/2.d0 * eccentricity**2.d0 + 1.d0/8.d0 * eccentricity**4.d0
  
  ! Omega and eccentricity dependent term (different term for AM derivative)
  omega_eccen_term = f3_e - (11.d0/18.d0)*f4_e*(omegi(1)/romorb)*(1.d0 - eccentricity**2.d0)**(1.5d0)
  
  ! slm coefficients for the dynamical tides
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
  ! Updated calculation of orbital angular momentum, relevant for eccentric orbits
  orbang = gms*Msol*binM2*Msol*sqrt(cst_G*ab*(1.d0 - eccentricity**2.d0)/((gms + binM2) * Msol))
  orstwi=dm_lost*Msol*(qtot*ab)**2.d0*romorb*sqrt(1.d0 - eccentricity**2.d0)
  rorstw=orstwi/dzeit
  term1=-2.d0*dLtid/orbang
  rltid=dLtid/dzeit
  term2=-2.d0*dm_lost/gms
  term3=dm_lost/(binm2+gms)
  term4=2.d0*orstwi/orbang
  
  if (include_dyn_tides) then
    ! deccentricity dynamical tides computed using equation 9 in Sciarini+24
    deccentricity=-(3.d0/4.d0)*eccentricity*fact2*qr1*sqrt(1+qr1)*e2*(rstar/ab)**(13.d0/2.d0)*((3.d0/2.d0)*&
    sign(s10**(8.d0/3.d0),difsav10)-(1.d0/4.d0)*sign(s12**(8.d0/3.d0),difsav12)-sign(s22**(8.d0/3.d0),difsav)+(49.d0/4.d0)*&
    sign(s32**(8.d0/3.d0),difsav32))*dzeit
  else
    deccentricity=0.d0
  end if

  ! deccentricity equilibrium tides following Hut 81, Sciarini et al. 2025

  if (posyd_prescription) then
    ! Posydon prescription for equilibrium tides on small convective regions. Independent on k2. Ratio of mass, tau_conv
    ! according to Eq. (7)
    deccentricity_eq = -2.d0/21.d0*27.d0*eccentricity*Mconv_fconv_ov_tau*qr1*(1+qr1)*(rstar/ab)**(8.d0)*omega_eccen_term/&
    (gms*Msol*(1-eccentricity**2.d0)**6.5d0)*dzeit
  else
  ! Default: Sciarini et al. 2025 (in prep.) prescription. Depends on k2, ratio of moment of inertia instead of ratio of
  ! masses, tau_conv self-consistent (obtained from MLT)  
    deccentricity_eq = -k2_AMC*27.d0*eccentricity*Iconv_fconv_ov_tau_I*qr1*(1+qr1)*(rstar/ab)**(8.d0)*omega_eccen_term/&
    ((1-eccentricity**2.d0)**6.5d0)*dzeit
  endif


  ! Add up the two contributions
  if (include_eq_tides) then
    deccentricity = deccentricity + deccentricity_eq
  end if
  
  term_eccentricity=(2.d0*eccentricity*deccentricity)/(1-eccentricity**2.d0)
    
  ! Variation of separation in the eccentric case
  dab=ab*(term1+term2+term3+term4+term_eccentricity)
  
  ! In case of a twin system (i.e. system of equal mass), we consider that the companion is identical 
  ! and brings the same contribution to the orbital evolution (same mass loss, same torques, ...) so 
  ! Delta a and Delta e can just be multiplied by 2.
  if (twin_system) then
    dab = dab*2.d0
    deccentricity = deccentricity*2.d0
  endif
  eccentricity=eccentricity+deccentricity
  ab=ab+dab
  romorb=(cst_G*(binm2+gms)*Msol/ab**3.0d0)**0.5d0
  period=2.0*pi/romorb
  write(*,*)'semi major axis: a = ',ab
  write(*,*)'eccentricity= ',eccentricity
  write(*,*)'period [days]= ',period/day
  if (verbose) then
    write(*,*) 'alter=',alter,'dltid=',dltid,'orstwi=',orstwi,'romini=',romini,'rmw=',dm_lost/(dzeit/year),&
       'term1=',term1,'term2=',term2,'term3=',term3,'dab=', dab/ab,'binm2=',binm2,'rstar/rsun=',rstar/Rsol,&
       'ab*qtot/rsun=',ab*qtot/Rsol
  endif
  write(io_logs,*) 'alter=',alter,'dltid=',dltid,'orstwi=',orstwi,'romini=',romini,'rmw=',dm_lost/(dzeit/year),&
       'term1=',term1,'term2=',term2,'term3=',term3,'dab=', dab/ab,'binm2=',binm2,'rstar/rsun=',rstar/Rsol,&
       'ab*qtot/rsun=',ab*qtot/Rsol
  write(io_logs,'(a11,e12.6)') 'new period=',period/day
  write(io_period_evol,'(i7,1x,e22.15,1x,f9.3,9(1x,es13.6))') nwmd,alter,rstar/Rsol,omegi(1),romorb,period/day,ab,eccentricity,&
  (dltid-dltid_eq)/dzeit,dltid_eq/dzeit,(deccentricity - deccentricity_eq)/dzeit,deccentricity_eq/dzeit

  return

end subroutine orbitalevol
!======================================================================
end module bintidemod

