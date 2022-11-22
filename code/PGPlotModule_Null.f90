!********************************************************
module PGPlotModule
! Interface to use PGplot in parallel with computation
! In order to be compatible with pgplot, all the value
! should be real(4) !!!
!********************************************************
use evol,only: kindreal
use inputparam,only: display_plot

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
integer,private,save:: Data_Number
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
public:: PlotEvol
private:: FindCZ
public:: Mass_Vector
private:: Estimate_Lifetime
private:: Find_Max_Energy

contains
!======================================================================
subroutine InitPGplot
! Module initialisation
!----------------------------------------------------------------------
  implicit none

  integer:: ierror,Number
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

  call Mass_Vector
  if (display_plot) then
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
    if (display_plot) then
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
                min_Si_Burning = 1.d10
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
! Null version of the routine in case PGPlot is not active
!----------------------------------------------------------------------

  implicit none

  return

end subroutine PlotEvol
!======================================================================
subroutine PlotStruc
!  Null version of the routine in case PGPlot is not active
!----------------------------------------------------------------------

  implicit none

  return

end subroutine PlotStruc
!======================================================================
subroutine EndPGplot
! Null version of the routine in case PGPlot is not active
!----------------------------------------------------------------------
  implicit none
!----------------------------------------------------------------------

  return

end subroutine EndPGplot
!======================================================================

end module PGPlotModule
!********************************************************
