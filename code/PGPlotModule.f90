!********************************************************
module PGPlotModule
! Interface to use PGplot in parallel with computation
! In order to be compatible with pgplot, all the value
! should be real(4) !!!
!********************************************************
use evol,only: kindreal
use inputparam,only: plot,refresh

implicit none

integer,private,parameter:: save_unit = 101,Data_Max = 10000, Data_Size = 31, CZ_resolution = 200
integer,private,parameter:: i_time = 1, i_mass = 2, i_L = 3, i_Teff = 4, i_Tc = 5, i_rhoc = 6, i_MaxE_H = 7, &
                              i_Mean_H_up = 8, i_Mean_H_down = 9, i_MaxE_He = 10, i_Mean_He_up = 11, &
                              i_Mean_He_down = 12, i_MaxE_C = 13, i_Mean_C_up = 14, i_Mean_C_down = 15, &
                              i_MaxE_Ne = 16, i_Mean_Ne_up = 17, i_Mean_Ne_down = 18, i_MaxE_O = 19, &
                              i_Mean_O_up = 20, i_Mean_O_down = 21, i_MaxE_Si = 22, i_Mean_Si_up = 23, &
                              i_Mean_Si_down = 24, i_H = 25, i_He = 26, i_C = 27, i_N = 28, i_O = 29, i_Ne = 30, &
                              i_Si = 31
integer,public,parameter:: Chem_Species_Number = 7
integer,public,save:: restart
integer,private,save:: Data_Number,Device_Number_1,Device_Number_2
integer,dimension(Data_Max),save:: DataLines
integer,dimension(Data_Max,CZ_resolution),private,save:: CZData

real(kindreal),dimension(CZ_resolution),private,save:: CZ_yaxis
real(kindreal),dimension(Data_Size),private,save:: PreviousData
real(kindreal),dimension(Data_Max,Data_Size),private,save:: SavedData

logical,private,save:: Continue_Writing
logical,public,save:: Struc_Plotted

character(256),public,save:: HRD_FileName

public:: InitPGplot
public:: SavePlotData
public:: EndPGplot
public:: PlotStruc
private:: PlotEvol
private:: FindCZ
private:: Mass_Vector
private:: Estimate_Lifetime
private:: Find_Max_Energy

contains
!======================================================================
subroutine InitPGplot
! Module initialisation
!----------------------------------------------------------------------
  implicit none

  integer:: ierror,Number,pgopen
  integer,dimension(CZ_resolution):: CZRead

  real(kindreal),dimension(Data_Size):: ReadData
!----------------------------------------------------------------------
  open(unit=save_unit,file=HRD_FileName,status='old',form='unformatted',iostat=ierror)
  if (ierror /= 0) then
    restart = 0
    open(unit=save_unit,file=HRD_FileName,status='new',form='unformatted')
  endif

  Continue_Writing = .true.
  Struc_Plotted = .false.

  if (restart > 0) then
    Data_Number = 1
    Number = -1
    do while (Number < restart)
      read(save_unit,iostat = ierror) Number,ReadData(1:Data_Size),CZRead(1:CZ_resolution)
      if (ierror /= 0) then
        if (ierror < 0) then
          PreviousData(:) = SavedData(Data_Number-1,:)
          backspace(save_unit)
          exit
        else
          write(*,*) 'Problem reading plotting data file. Aborting...'
          stop
        endif
      endif

      if (Number >= restart) then
! In case the file contains to much data, we need to remove the end of the file. This is done only if something
! is written in the file.
        Data_Number = Data_Number-1
        PreviousData(:) = SavedData(Data_Number,:)
        backspace(save_unit)
        backspace(save_unit)
        write(save_unit) DataLines(Data_Number),SavedData(Data_Number,:),CZdata(Data_Number,:)
        exit
      endif

      DataLines(Data_Number) = Number
      SavedData(Data_Number,:) = ReadData(:)
      CZdata(Data_Number,:) = CZRead(:)

      Data_Number = Data_Number + 1
    enddo
  else
! If Data_Number is set to 1, the plot is beginning at the current HRD position
    Data_Number = 1
    rewind(save_unit)
    CZData(:,:) = 0
  endif

  if (refresh) then
    Device_Number_1=pgopen('/XSERVE')
    Device_Number_2=pgopen('/XSERVE')
    if (Device_Number_1 < 0 .or. Device_Number_2 < 0) then
      write(*,*) 'problems opening plotting device, aborting...'
      stop
    endif
  endif

  if (refresh) then
    call pgask(.false.)
  endif
  call Mass_Vector
  if (refresh) then
    call PlotEvol
  endif

  return

end subroutine InitPGplot
!======================================================================
subroutine SavePlotData(mass,L,Teff,Number,time,Tc,rhoc,Species)
! Routine saving the HRD data and printing them in a .hrd file
!----------------------------------------------------------------------
  use inputparam,only: idebug

  implicit none

  integer, intent(in):: Number

  real(kindreal),dimension(Data_Size), save:: Tol = (/0.d0,0.d0,5.d-3,5.d-3,1.d-2,4.d-2,0.d0,0.d0,0.d0,0.d0,0.d0, &
                                                    0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
                                                    1.d-2,1.d-2,1.d-2,1.d-2,1.d-2,1.d-2,1.d-2/)

  real(kindreal), intent(in):: mass, L, Teff, time, Tc, rhoc
  real(kindreal), dimension(Chem_Species_Number), intent(in):: Species

  real(kindreal):: Diff_Species, Tau_H
!----------------------------------------------------------------------
! if no more space is left, abort.
  if (.not.Continue_Writing) then
    return
  endif

! Estimated lifetime:
  if (Data_Number == 1 .or. idebug > 0) then
    Tau_H = Estimate_Lifetime(mass)
  else
    Tau_H = Estimate_Lifetime(SavedData(1,i_mass))
  endif
  Tol(i_time) = Tau_H/100.d0

! Set the previous L and Teff
  if (Data_Number == 1) then
    PreviousData(:) = -99.d0
  endif

! Convert L and Teff into log (real4) values if the current values are suffiently different from the
! previous ones.
  Diff_Species = maxval(abs(Species(:)-PreviousData(i_H:i_H+Chem_Species_Number-1)))
  if (abs(log10(L) - PreviousData(i_L)) > Tol(i_L) .or. &
      abs(log10(Teff) - PreviousData(i_Teff)) > Tol(i_Teff) .or. &
      abs(time - PreviousData(i_time)*1.d6) > Tol(i_time) .or. &
      abs(Tc - PreviousData(i_Tc)) > Tol(i_Tc) .or. &
      abs(rhoc - PreviousData(i_rhoc)) > Tol(i_rhoc) .or. &
      Diff_Species > Tol(i_H)) then
    SavedData(Data_Number,i_mass) = mass
    SavedData(Data_Number,i_L) = log10(L)
    SavedData(Data_Number,i_Teff) = log10(Teff)
    SavedData(Data_Number,i_time) = time/1.d6
    SavedData(Data_Number,i_Tc) = Tc
    SavedData(Data_Number,i_rhoc) = rhoc
    SavedData(Data_Number,i_H:i_H+Chem_Species_Number-1) = Species(:)
    DataLines(Data_Number) = Number

    if (idebug > 0) then
      write(*,*) 'PGPlotModule: call FindCZ'
    endif
    call FindCZ
    if (idebug > 0) then
      write(*,*) 'PGPlotModule: call Find_Max_Energy'
    endif
    call Find_Max_Energy
    write(save_unit) DataLines(Data_Number),SavedData(Data_Number,1:Data_Size), &
                     CZdata(Data_Number,:)

! Save previous data
    PreviousData(:) = SavedData(Data_Number,:)

! Increase the data number
    Data_Number = Data_Number + 1

    if (Data_Number > Data_Max) then
      write(*,*) 'No more space in data vectors in PGplot module'
      Continue_Writing = .false.
    endif
    restart = Number
    if (refresh) then
      if (idebug > 0) then
        write(*,*) 'PGPlotModule: call PlotEvol'
      endif
      call PlotEvol
      if (idebug > 0) then
        write(*,*) 'PGPlotModule: call PlotStruc'
      endif
      call PlotStruc
    endif
  endif
  if (idebug > 0) then
    write(*,*) 'PGPlotModule: end'
  endif

  return

end subroutine SavePlotData
!======================================================================
subroutine FindCZ
! Data for Kippenhahn diagram
!----------------------------------------------------------------------
  use strucmod, only:q,zensi,shell_number => m

  implicit none

  integer:: i,j,Exit_Value
  integer, dimension(shell_number)::zensi_downup
  integer, dimension(CZ_resolution):: zensi_reduced

  real(kindreal):: mass_change
  real(kindreal), dimension(shell_number):: StellarModel_Mass
!----------------------------------------------------------------------
  call Mass_Vector

  zensi_reduced(:) = 0

  do i=shell_number,1,-1
    StellarModel_Mass(shell_number-i+1) = (1.d0 - exp(q(i)))*SavedData(Data_Number,i_mass)
    if (zensi(i) > 0.d0) then
      zensi_downup(shell_number-i+1) = 1
    else
      zensi_downup(shell_number-i+1) = 0
    endif
  enddo

  i=1
  Exit_Value = 0
  LargeLoop: do j=2,shell_number-1
    if (zensi_downup(j) /=zensi_downup(j+1)) then
      mass_change = StellarModel_Mass(j)
! If the size of an intermediate RZ is very small, we miss the change. In that case, force change.
      if (CZ_yaxis(i) > mass_change .and. CZ_yaxis(i ) <SavedData(Data_Number,i_mass)) then
        zensi_reduced(i) = zensi_downup(j)
        if (i == CZ_resolution) then
          Exit_Value = CZ_resolution
          exit LargeLoop
        endif
        Exit_Value = i
        i = i+1
      endif
      do while (CZ_yaxis(i) < mass_change .and. CZ_yaxis(i) < SavedData(Data_Number,i_mass))
        zensi_reduced(i) = zensi_downup(j)
        if (i == CZ_resolution) then
          Exit_Value = CZ_resolution
          exit LargeLoop
        endif
        Exit_Value = i
        i = i+1
      enddo
    endif
  enddo LargeLoop

! case of a fully radiative star:
  if (Exit_Value == 0) then
    zensi_reduced(:) = 0
  else if (Exit_Value /= CZ_resolution) then
    do j=Exit_Value+1,CZ_resolution
      if (CZ_yaxis(j) < SavedData(Data_Number,i_mass)) then
        if (zensi_reduced(Exit_Value) == 0) then
          zensi_reduced(j) = 1
        else
          zensi_reduced(j) = 0
        endif
      else
        zensi_reduced(j) = 0
      endif
    enddo
  endif

  CZData(Data_Number,:) = zensi_reduced(:)

  return

end subroutine FindCZ
!======================================================================
subroutine Mass_Vector
! Determination of the mass vector
!----------------------------------------------------------------------
  implicit none

  integer:: i

  real(kindreal):: Delta_m
!----------------------------------------------------------------------
! Define the y-axis for Kippenhahn diagram
  Delta_m = 1.d0/real(CZ_resolution)
  do i=1,CZ_resolution
    CZ_yaxis(i) = Delta_m/2.d0 + real(i-1)*Delta_m
  enddo
  CZ_yaxis(:) = CZ_yaxis(:)*SavedData(1,i_mass)

  return

end subroutine Mass_Vector
!======================================================================
real(kindreal) function Estimate_Lifetime(mass)
! Estimation of the lifetime (Ekstrom et al. 2012)
!----------------------------------------------------------------------
  implicit none

  real(kindreal),dimension(4),parameter::A_fit = (/-2.776d0,-2.444d0,-1.763d0,-0.775d0/), &
                                         B_fit = (/ 9.938d0, 9.761d0, 9.186d0, 8.004d0/)
  real(kindreal),intent(in):: mass
!----------------------------------------------------------------------
  if (mass < 3.d0) then
    Estimate_Lifetime = 10.d0**(A_fit(1)*log10(mass)+B_fit(1))
  else if (mass < 7.d0) then
    Estimate_Lifetime = 10.d0**(A_fit(2)*log10(mass)+B_fit(2))
  else if (mass < 15.d0) then
    Estimate_Lifetime = 10.d0**(A_fit(3)*log10(mass)+B_fit(3))
  else
    Estimate_Lifetime = 10.d0**(A_fit(4)*log10(mass)+B_fit(4))
  endif

  return

end function Estimate_Lifetime
!======================================================================
subroutine Find_Max_Energy
! Find maximum energy production by H, He burning
!----------------------------------------------------------------------
  use abundmod,only: eps,epsy,eps_c_adv,eps_ne_adv,eps_o_adv,eps_si_adv
  use strucmod,only: q,shell_number => m

  implicit none

  integer:: i_shell_H,i_shell_He,i_shell_C,i_shell_Ne,i_shell_O,i_shell_Si,i

  real(kindreal):: Max_H,Max_He,Max_C,Max_Ne,Max_O,Max_Si
  real(kindreal),parameter:: threshold = 0.1d0,min_H_Burning = 1.d2, min_He_Burning = 1.d3, &
                min_C_Burning = 1.d4, min_Ne_Burning = 1.d0, min_O_Burning = 1.d6, &
                min_Si_Burning = 1.d6
!----------------------------------------------------------------------
! Localisation of the maximal energy production (H-b, He-b)
  i_shell_H = maxloc(eps(1:shell_number),1)
  i_shell_He = maxloc(epsy(1:shell_number),1)
  i_shell_C = maxloc(eps_c_adv(1:shell_number),1)
  i_shell_Ne = maxloc(eps_ne_adv(1:shell_number),1)
  i_shell_O = maxloc(eps_o_adv(1:shell_number),1)
  i_shell_Si = maxloc(eps_si_adv(1:shell_number),1)

  Max_H = maxval(eps(1:shell_number))
  Max_He = maxval(epsy(1:shell_number))
  Max_C = maxval(eps_c_adv(1:shell_number))
  Max_Ne = maxval(eps_ne_adv(1:shell_number))
  Max_O = maxval(eps_o_adv(1:shell_number))
  Max_Si = maxval(eps_si_adv(1:shell_number))

! Don't plot whether energy production is too small
  if (Max_H >= min_H_Burning) then
    SavedData(Data_Number,i_MaxE_H) = (1.d0 - exp(q(i_shell_H)))*SavedData(Data_Number,i_mass)
  else
    SavedData(Data_Number,i_MaxE_H) = -0.05d0
  endif
  if (Max_He >= min_He_Burning) then
    SavedData(Data_Number,i_MaxE_He) = (1.d0 - exp(q(i_shell_He)))*SavedData(Data_Number,i_mass)
  else
    SavedData(Data_Number,i_MaxE_He) = -0.05d0
  endif
  if (Max_C >= min_C_Burning) then
    SavedData(Data_Number,i_MaxE_C) = (1.d0 - exp(q(i_shell_C)))*SavedData(Data_Number,i_mass)
  else
    SavedData(Data_Number,i_MaxE_C) = -0.05d0
  endif
  if (Max_Ne >= min_Ne_Burning) then
    SavedData(Data_Number,i_MaxE_Ne) = (1.d0 - exp(q(i_shell_Ne)))*SavedData(Data_Number,i_mass)
  else
    SavedData(Data_Number,i_MaxE_Ne) = -0.05d0
  endif
  if (Max_O >= min_O_Burning) then
    SavedData(Data_Number,i_MaxE_O) = (1.d0 - exp(q(i_shell_O)))*SavedData(Data_Number,i_mass)
  else
    SavedData(Data_Number,i_MaxE_O) = -0.05d0
  endif
  if (Max_Si >= min_Si_Burning) then
    SavedData(Data_Number,i_MaxE_Si) = (1.d0 - exp(q(i_shell_Si)))*SavedData(Data_Number,i_mass)
  else
    SavedData(Data_Number,i_MaxE_Si) = -0.05d0
  endif

! Locate the 10% production level
  if (SavedData(Data_Number,i_MaxE_H) >= 0.d0) then
    SavedData(Data_Number,i_Mean_H_down) = 0.d0
    do i=i_shell_H,shell_number
      if (eps(i) <= threshold*Max_H) then
        SavedData(Data_Number,i_Mean_H_down) = (1.d0 - exp(q(i)))*SavedData(Data_Number,i_mass)
        exit
      endif
    enddo
    SavedData(Data_Number,i_Mean_H_up) = (1.d0 - exp(q(1)))*SavedData(Data_Number,i_mass)
    do i=i_shell_H,1,-1
      if (eps(i) <= threshold*Max_H) then
        SavedData(Data_Number,i_Mean_H_up) = (1.d0 - exp(q(i)))*SavedData(Data_Number,i_mass)
        exit
      endif
    enddo
  else
    SavedData(Data_Number,i_Mean_H_down) = -0.05d0
    SavedData(Data_Number,i_Mean_H_up) = -0.05d0
  endif

  if (SavedData(Data_Number,i_MaxE_He) >= 0.d0) then
    SavedData(Data_Number,i_Mean_He_down) = 0.d0
    do i=i_shell_He,shell_number
      if (epsy(i) <= threshold*Max_He) then
        SavedData(Data_Number,i_Mean_He_down) = (1.d0 - exp(q(i)))*SavedData(Data_Number,i_mass)
        exit
      endif
    enddo
    SavedData(Data_Number,i_Mean_He_up) = (1.d0 - exp(q(1)))*SavedData(Data_Number,i_mass)
    do i=i_shell_He,1,-1
      if (epsy(i) <= threshold*Max_He) then
        SavedData(Data_Number,i_Mean_He_up) = (1.d0 - exp(q(i)))*SavedData(Data_Number,i_mass)
        exit
      endif
    enddo
  else
    SavedData(Data_Number,i_Mean_He_down) = -0.05d0
    SavedData(Data_Number,i_Mean_He_up) = -0.05d0
  endif

  if (SavedData(Data_Number,i_MaxE_C) >= 0.d0) then
    SavedData(Data_Number,i_Mean_C_down) = 0.d0
    do i=i_shell_C,shell_number
      if (eps_c_adv(i) <= threshold*Max_C) then
        SavedData(Data_Number,i_Mean_C_down) = (1.d0 - exp(q(i)))*SavedData(Data_Number,i_mass)
        exit
      endif
    enddo
    SavedData(Data_Number,i_Mean_C_up) = (1.d0 - exp(q(1)))*SavedData(Data_Number,i_mass)
    do i=i_shell_C,1,-1
      if (eps_c_adv(i) <= threshold*Max_C) then
        SavedData(Data_Number,i_Mean_C_up) = (1.d0 - exp(q(i)))*SavedData(Data_Number,i_mass)
        exit
      endif
    enddo
  else
    SavedData(Data_Number,i_Mean_C_down) = -0.05d0
    SavedData(Data_Number,i_Mean_C_up) = -0.05d0
  endif

  if (SavedData(Data_Number,i_MaxE_Ne) >= 0.d0) then
    SavedData(Data_Number,i_Mean_Ne_down) = 0.d0
    do i=i_shell_Ne,shell_number
      if (eps_ne_adv(i) <= threshold*Max_Ne) then
        SavedData(Data_Number,i_Mean_Ne_down) = (1.d0 - exp(q(i)))*SavedData(Data_Number,i_mass)
        exit
      endif
    enddo
    SavedData(Data_Number,i_Mean_Ne_up) = (1.d0 - exp(q(1)))*SavedData(Data_Number,i_mass)
    do i=i_shell_Ne,1,-1
      if (eps_ne_adv(i) <= threshold*Max_Ne) then
        SavedData(Data_Number,i_Mean_Ne_up) = (1.d0 - exp(q(i)))*SavedData(Data_Number,i_mass)
        exit
      endif
    enddo
  else
    SavedData(Data_Number,i_Mean_Ne_down) = -0.05d0
    SavedData(Data_Number,i_Mean_Ne_up) = -0.05d0
  endif

  if (SavedData(Data_Number,i_MaxE_O) >= 0.d0) then
    SavedData(Data_Number,i_Mean_O_down) = 0.d0
    do i=i_shell_O,shell_number
      if (eps_o_adv(i) <= threshold*Max_O) then
        SavedData(Data_Number,i_Mean_O_down) = (1.d0 - exp(q(i)))*SavedData(Data_Number,i_mass)
        exit
      endif
    enddo
    SavedData(Data_Number,i_Mean_O_up) = (1.d0 - exp(q(1)))*SavedData(Data_Number,i_mass)
    do i=i_shell_O,1,-1
      if (eps_o_adv(i) <= threshold*Max_O) then
        SavedData(Data_Number,i_Mean_O_up) = (1.d0 - exp(q(i)))*SavedData(Data_Number,i_mass)
        exit
      endif
    enddo
  else
    SavedData(Data_Number,i_Mean_O_down) = -0.05d0
    SavedData(Data_Number,i_Mean_O_up) = -0.05d0
  endif

  if (SavedData(Data_Number,i_MaxE_Si) >= 0.d0) then
    SavedData(Data_Number,i_Mean_Si_down) = 0.d0
    do i=i_shell_Si,shell_number
      if (eps_si_adv(i) <= threshold*Max_Si) then
        SavedData(Data_Number,i_Mean_Si_down) = (1.d0 - exp(q(i)))*SavedData(Data_Number,i_mass)
        exit
      endif
    enddo
    SavedData(Data_Number,i_Mean_Si_up) = (1.d0 - exp(q(1)))*SavedData(Data_Number,i_mass)
    do i=i_shell_Si,1,-1
      if (eps_si_adv(i) <= threshold*Max_Si) then
        SavedData(Data_Number,i_Mean_Si_up) = (1.d0 - exp(q(i)))*SavedData(Data_Number,i_mass)
        exit
      endif
    enddo
  else
    SavedData(Data_Number,i_Mean_Si_down) = -0.05d0
    SavedData(Data_Number,i_Mean_Si_up) = -0.05d0
  endif

  return

end subroutine Find_Max_Energy
!======================================================================
subroutine PlotEvol
! Plotting routine for evolution data
!----------------------------------------------------------------------
  use inputparam,only: phase

  implicit none

  integer::i,j,CZ_min,CZ_max

  real(kindreal):: L_min,L_max,T_min,T_max,Tc_min,Tc_max,rhoc_min,rhoc_max,time_min,time_max,kippen_min,kippen_max, &
                 delta
  real(kindreal),dimension(2):: Kippen_X_vector,Kippen_Y_vector
  real(kindreal),allocatable:: time_adv(:)

  logical::Is_Convective
  allocate(time_adv(10000))
!----------------------------------------------------------------------
  time_adv(:) = 0.d0

  if (Data_Number > 2) then
    L_max = maxval(SavedData(1:Data_Number-1,i_L))
    L_min = minval(SavedData(1:Data_Number-1,i_L))
    delta = (L_max - L_min)/20.d0
    L_max = L_max + delta
    L_min = L_min - delta
    T_max = maxval(SavedData(1:Data_Number-1,i_Teff))
    T_min = minval(SavedData(1:Data_Number-1,i_Teff))
    delta = (T_max - T_min)/20.d0
    T_max = T_max + delta
    T_min = T_min - delta
    Tc_max = maxval(SavedData(1:Data_Number-1,i_Tc))
    Tc_min = minval(SavedData(1:Data_Number-1,i_Tc))
    delta = (Tc_max - Tc_min)/20.d0
    Tc_max = Tc_max + delta
    Tc_min = Tc_min - delta
    rhoc_max = maxval(SavedData(1:Data_Number-1,i_rhoc))
    rhoc_min = minval(SavedData(1:Data_Number-1,i_rhoc))
    delta = (rhoc_max - rhoc_min)/20.d0
    rhoc_max = rhoc_max + delta
    rhoc_min = rhoc_min - delta
    if (phase > 1) then
      time_adv(:) = log10(2.d0*SavedData(Data_Number-1,i_time)-SavedData(Data_Number-2,i_time)-SavedData(:,i_time))
    else
      time_adv(:) = SavedData(:,i_time)
    endif
    time_max = maxval(time_adv(1:Data_Number-1))
    time_min = minval(time_adv(1:Data_Number-1))
    delta = (time_max - time_min)/20.d0
    time_max = time_max + delta
    time_min = time_min - delta
    kippen_min = 0.d0
    kippen_max = SavedData(1,i_mass)

    call pgask(.false.)
    call pgslct(Device_Number_1)
    call pgsubp(2,2)
    call pgpage
    call pgvstd
    call pgsls(1)
    call pgsci(1)
    call pgswin(real(T_max),real(T_min),real(L_min),real(L_max))
    call pgbox('BCNTS',0.0,5,'BCNTS',0.0,5)
    call pglab('log(Teff)','log(L)','HRD')
    call pgline(Data_Number-1,real(SavedData(1:Data_Number-1,i_Teff)),real(SavedData(1:Data_Number-1,i_L)))
    call pgsch(2.5)
    call pgsci(7)
    call pgpt1(real(SavedData(Data_Number-1,i_Teff)),real(SavedData(Data_Number-1,i_L)),12)
    call pgsci(1)
    call pgsch(1.)
    call pgpage
    call pgvstd
    call pgswin(real(rhoc_min),real(rhoc_max),real(Tc_min),real(Tc_max))
    call pgbox('BCNTS',0.0,5,'BCNTS',0.0,5)
    call pglab('log(rhoc)','log(Tc)','Tc-rhoc')
    call pgsls(1)
    call pgsci(1)
    call pgline(Data_Number-1,real(SavedData(1:Data_Number-1,i_rhoc)),real(SavedData(1:Data_Number-1,i_Tc)))
    call pgsch(2.5)
    call pgsci(7)
    call pgpt1(real(SavedData(Data_Number-1,i_rhoc)),real(SavedData(Data_Number-1,i_tc)),12)
    call pgsci(1)
    call pgsch(1.)
    call pgpage
    call pgvstd
    if (phase <= 1) then
      call pgswin(real(time_min),real(time_max),-0.05,1.05)
    else
      call pgswin(real(time_max),real(time_min),-0.05,1.05)
    endif
    call pgbox('BCNTS',0.,5,'BCNTS',0.2,4)
    call pglab('time','mass fraction','central abundances')
    call pgsls(1)
    call pgsci(1)
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_H)))
    call pgsci(3)
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_He)))
    call pgsci(2)
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_C)))
    call pgsci(4)
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_O)))
    call pgsci(5)
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_Ne)))
    call pgsci(8)
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_Si)))
    call pgpage
    call pgvstd
    call pgsci(1)
    call pgsls(1)
    call pgslw(10)
    call pgsci(14)
    if (phase <= 1) then
      call pgswin(real(time_min),real(time_max),real(kippen_min),real(kippen_max))
    else if (phase < 3) then
      call pgswin(real(time_max),real(time_min),real(kippen_min),real(kippen_max))
    else
      call pgswin(real(time_max),real(time_min),real(kippen_min),real(kippen_max)/3.)
    endif
    do i=1,Data_Number-1
      Is_Convective = .false.
      do j=1,CZ_resolution
        if (CZData(i,j) == 1 .and. .not.Is_Convective) then
          Is_Convective = .true.
          CZ_min = j
        endif
        if ((CZData(i,j) == 0 .or. j==CZ_resolution) .and. Is_Convective) then
          Is_Convective = .false.
          CZ_max = j
          Kippen_X_vector(:) = time_adv(i)
          Kippen_Y_vector(1) = CZ_yaxis(CZ_min)
          Kippen_Y_vector(2) = CZ_yaxis(CZ_max)
          call pgline(2,real(Kippen_X_vector),real(Kippen_Y_vector))
        endif
      enddo
    enddo
    call pgsci(1)
    call pgslw(1)
    call pgbox('BCNTS',0.,5,'BCNTS',0.,5)
    call pglab('time','mass','Kippenhahn diagram')
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_mass)))
    call pgslw(2)
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_MaxE_H)))
    call pgsls(2)
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_Mean_H_up)))
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_Mean_H_down)))
    call pgsls(1)
    call pgsci(3)
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_MaxE_He)))
    call pgsls(2)
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_Mean_He_up)))
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_Mean_He_down)))
    call pgsls(1)
    call pgsci(2)
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_MaxE_C)))
    call pgsls(2)
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_Mean_C_up)))
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_Mean_C_down)))
    call pgsls(1)
    call pgsci(5)
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_MaxE_Ne)))
    call pgsls(2)
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_Mean_Ne_up)))
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_Mean_Ne_down)))
    call pgsls(1)
    call pgsci(4)
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_MaxE_O)))
    call pgsls(2)
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_Mean_O_up)))
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_Mean_O_down)))
    call pgsls(1)
    call pgsci(8)
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_MaxE_Si)))
    call pgsls(2)
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_Mean_Si_up)))
    call pgline(Data_Number-1,real(time_adv(1:Data_Number-1)),real(SavedData(1:Data_Number-1,i_Mean_Si_down)))
    call pgsls(1)
    call pgsci(1)
    call pgslw(1)
  endif

  return

end subroutine PlotEvol
!======================================================================
subroutine PlotStruc
! Plotting routine for evolution data
!----------------------------------------------------------------------

  use const,only: um
  use strucmod,only: q,s,p,t,Nabla_ad,Nabla_rad,Nabla_mu,shell_number => m
  use abundmod,only: x,y,xc12,xn14,xo16,xne20,eps,epsy,eps_c_adv,eps_ne_adv,eps_o_adv,eps_si_adv,eps_grav,eps_nu,abelx

  implicit none

  integer:: i

  real(kind=8),dimension(shell_number):: StellarModel_Mass,H_Profile,He_Profile,C_Profile,N_Profile,O_Profile,Ne_Profile, &
    Si_Profile,Ni_Profile,Fe_Profile,Nabla_ad_plot,Nabla_rad_plot,Nabla_mu_plot,Zero,Luminosity,eps_H_log,eps_He_log, &
    eps_C_log,eps_Ne_log,eps_O_log,eps_Si_log,temperature,pressure,eps_grav_moins,eps_grav_plus,eps_nu_log
  real(kindreal):: Min_Nabla,Max_Nabla,NormLum,Max_Lum,Min_Lum,Max_eps,Min_eps,Max_T,Min_T,Max_P,Min_P
!----------------------------------------------------------------------
  if (Data_Number > 2) then
    do i=shell_number,1,-1
      StellarModel_Mass(shell_number-i+1) = (1.d0 - exp(q(i)))
      if (x(i) > 0.d0) then
        H_Profile(shell_number-i+1) = log10(x(i))
      else
        H_Profile(shell_number-i+1) = -32.d0
      endif
      if (y(i) > 0.d0) then
        He_Profile(shell_number-i+1) = log10(y(i))
      else
        He_Profile(shell_number-i+1) = -32.d0
      endif
      if (xc12(i) > 0.d0) then
        C_Profile(shell_number-i+1) = log10(xc12(i))
      else
        C_Profile(shell_number-i+1) = -32.d0
      endif
      if (xn14(i) > 0.d0) then
        N_Profile(shell_number-i+1) = log10(xn14(i))
      else
        N_Profile(shell_number-i+1) = -32.d0
      endif
      if (xo16(i) > 0.d0) then
        O_Profile(shell_number-i+1) = log10(xo16(i))
      else
        O_Profile(shell_number-i+1) = -32.d0
      endif
      if (xne20(i) > 0.d0) then
        Ne_Profile(shell_number-i+1) = log10(xne20(i))
      else
        Ne_Profile(shell_number-i+1) = -32.d0
      endif
      if (abelx(1,i) > 0.d0) then
        Si_Profile(shell_number-i+1) = log10(abelx(1,i))
      else
        Si_Profile(shell_number-i+1) = -32.d0
      endif
      if (abelx(7,i) > 0.d0) then
        Ni_Profile(shell_number-i+1) = log10(abelx(7,i))
      else
        Ni_Profile(shell_number-i+1) = -32.d0
      endif
      if (abelx(8,i) > 0.d0) then
        Fe_Profile(shell_number-i+1) = log10(abelx(8,i))
      else
        Fe_Profile(shell_number-i+1) = -32.d0
      endif
      if (eps(i) > 0.d0) then
        eps_H_log(shell_number-i+1) = log10(eps(i))
      else
        eps_H_log(shell_number-i+1) = -32.d0
      endif
      if (epsy(i) > 0.d0) then
        eps_He_log(shell_number-i+1) = log10(epsy(i))
      else
        eps_He_log(shell_number-i+1) = -32.d0
      endif
      if (eps_c_adv(i) > 0.d0) then
        eps_C_log(shell_number-i+1) = log10(eps_c_adv(i))
      else
        eps_C_log(shell_number-i+1) = -32.d0
      endif
      if (eps_ne_adv(i) > 0.d0) then
        eps_Ne_log(shell_number-i+1) = log10(eps_ne_adv(i))
      else
        eps_Ne_log(shell_number-i+1) = -32.d0
      endif
      if (eps_o_adv(i) > 0.d0) then
        eps_O_log(shell_number-i+1) = log10(eps_o_adv(i))
      else
        eps_O_log(shell_number-i+1) = -32.d0
      endif
      if (eps_si_adv(i) > 0.d0) then
        eps_Si_log(shell_number-i+1) = log10(eps_si_adv(i))
      else
        eps_Si_log(shell_number-i+1) = -32.d0
      endif
      if (eps_nu(i) > 0.d0) then
        eps_nu_log(i) = log10(eps_nu(i))
      else
        eps_nu_log(i) = -32.d0
      endif
      if (eps_grav(i) > 0.d0) then
        eps_grav_plus(i) = log10(eps_grav(i))
        eps_grav_moins(i) = -32.d0
      elseif (eps_grav(i) < 0.d0) then
        eps_grav_plus(i) = -32.d0
        eps_grav_moins(i) = log10(-eps_grav(i))
      else
        eps_grav_plus(i) = -32.d0
        eps_grav_moins(i) = -32.d0
      endif


      Nabla_ad_plot(shell_number-i+1) = Nabla_ad(i)
      Nabla_rad_plot(shell_number-i+1) = Nabla_rad(i)
      Nabla_mu_plot(shell_number-i+1) = Nabla_mu(i)
      NormLum = 1.d0/(exp(s(1))-1.d0)
      Luminosity(shell_number-i+1) = (exp(s(i))-1.d0)*NormLum
      temperature(shell_number-i+1) = t(i)/um
      pressure(shell_number-i+1) = p(i)/um
    enddo

    Min_Nabla = min(minval(Nabla_ad_plot),minval(Nabla_rad_plot))
    Max_Nabla = max(maxval(Nabla_ad_plot),maxval(Nabla_rad_plot))
    Zero(:) = 0.d0
    Max_Lum = min(maxval(Luminosity),10.d0)
    Min_Lum = max(minval(Luminosity),-0.5d0)
    Max_eps = max(maxval(eps_H_log),maxval(eps_He_log),maxval(eps_C_log),maxval(eps_Ne_log),maxval(eps_O_log),maxval(eps_Si_log))
    Min_eps = max(minval(eps_H_log),2.d0)
    Max_T = maxval(temperature)
    Min_T = minval(temperature)
    Max_P = maxval(pressure)
    Min_P = minval(pressure)

    call pgask(.false.)
    call pgslct(Device_Number_2)
    call pgsubp(2,2)
    call pgpage
    call pgvstd
    call pgsls(1)
    call pgsci(1)   ! blanc
    call pgswin(0.,1.,-4.05,0.05)
    call pgbox('BCNTS',0.0,5,'BCNTS',0.0,5)
    call pglab('Mr/Mtot','mass frac.','Abundance profile')
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(H_Profile(1:shell_number)))
    call pgsci(3)   ! vert
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(He_Profile(1:shell_number)))
    call pgsci(2)   ! rouge
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(C_Profile(1:shell_number)))
    call pgsci(6)   ! magenta
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(N_Profile(1:shell_number)))
    call pgsci(4)   ! bleu
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(O_Profile(1:shell_number)))
    call pgsci(5)   ! cyan
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(Ne_Profile(1:shell_number)))
    call pgsci(8)   ! orange
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(Si_Profile(1:shell_number)))
    call pgsci(15)  ! gris clair
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(Ni_Profile(1:shell_number)))
    call pgsci(14)  ! gris fonce
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(Fe_Profile(1:shell_number)))
    if (Min_Nabla /= Max_Nabla) then
      Min_Nabla = min(Min_Nabla,-0.2d0)
      Min_Nabla = max(Min_Nabla,-10.d0)
      Max_Nabla = min(Max_Nabla,10.d0)
      call pgpage
      call pgvstd
      call pgsls(1)
      call pgsci(1)   ! blanc
      call pgswin(0.,1.,real(Min_Nabla),real(Max_Nabla))
      call pgbox('BCNTS',0.0,5,'BCNTS',0.0,5)
      call pglab('Mr/Mtot','Nabla_(rad-ad-mu).','gradient')
      call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(Nabla_ad_plot(1:shell_number)))
      call pgsci(2)   ! rouge
      call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(Nabla_rad_plot(1:shell_number)))
      call pgsci(4)   ! bleu
      call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(Nabla_mu_plot(1:shell_number)))
      call pgsci(1)
      call pgsls(2)
      call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(Zero(1:shell_number)))
      call pgsls(1)
    endif
    call pgpage
    call pgvstd
    call pgsls(1)
    call pgsci(1)
    call pgswin(0.,1.,real(Min_Lum),real(Max_Lum))
    call pgbox('BCNTS',0.0,5,'BNTS',0.0,5)
    call pglab('Mr/Mtot','L/L_tot','Luminosity and energy production')
    call pgsci(15)
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(Luminosity(1:shell_number)))
    call pgsci(1)
    call pgswin(0.,1.,real(Min_eps),real(Max_eps))
    call pgbox('',0.0,5,'CMTS',0.0,5)
    call pgmtxt('R',3.,0.5,0.5,'log(epsilon)')
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(eps_H_log(1:shell_number)))
    call pgsci(3)
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(eps_He_log(1:shell_number)))
    call pgsci(2)
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(eps_C_log(1:shell_number)))
    call pgsci(5)
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(eps_Ne_log(1:shell_number)))
    call pgsci(4)
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(eps_O_log(1:shell_number)))
    call pgsci(8)
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(eps_Si_log(1:shell_number)))
    call pgsci(5)
    call pgsls(2)
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(eps_nu_log(1:shell_number)))
    call pgsci(4)
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(eps_grav_plus(1:shell_number)))
    call pgsci(2)
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(eps_grav_moins(1:shell_number)))
    call pgpage
    call pgvstd
    call pgsls(1)
    call pgsci(1)
    call pgswin(0.,1.,real(Min_T),real(Max_T))
    call pgbox('BCNTS',0.0,5,'BNTS',0.0,5)
    call pglab('Mr/Mtot','log(T)','Temperature and pressure')
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(temperature(1:shell_number)))
    call pgswin(0.,1.,real(Min_P),real(Max_P))
    call pgbox('',0.0,5,'CMTS',0.0,5)
    call pgmtxt('R',3.,0.5,0.5,'log(P)')
    call pgsci(2)
    call pgline(shell_number,real(StellarModel_Mass(1:shell_number)),real(pressure(1:shell_number)))
    call pgsci(1)
  endif

  return

end subroutine PlotStruc
!======================================================================
subroutine EndPGplot
! Closing file and ending module
!----------------------------------------------------------------------
  implicit none
!----------------------------------------------------------------------
  close(save_unit)
  call pgclos()

  return

end subroutine EndPGplot
!======================================================================

end module PGPlotModule
!********************************************************
