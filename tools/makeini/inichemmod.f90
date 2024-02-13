!> \file inichemmod.f95
!> \brief Initial chemical composition module
!> 
!> Determines the initial chemical composition according to the metallicity
!! and the mixture requested. Writes the netdef.in and netalu.dat files.
module inichemmod

  use inputparam, only: zsol,iopac

  implicit none

  integer:: idefaut
  integer, dimension(20), parameter:: &
     elemA=(/1,3,4,12,13,14,15,16,17,18,20,22,24,25,26,19,21,23,27,28/), &
     elemZ=(/1,2,2,6,6,7,7,8,8,8,10,10,12,12,12,9,10,11,13,14/)

  real(kind=8), save:: zini,znew
  real(kind=8), dimension(20), save:: xx

  character(len=6), dimension(20), parameter::  &
     mainnam=(/'     x','    y3','     y','  xc12','  xc13','  xn14','  xn15','  xo16','  xo17','  xo18',&
               ' xne20',' xne22',' xmg24',' xmg25',' xmg26','  xf19',' xne21',' xna23',&
               ' xal27',' xsi28'/)

  integer:: n2                     !< number of isotopes read from the isotopic percentage input file
  integer,parameter:: nalpha=9,&   !< number of alpha enhanced isotopes
                      nely=83,&    !< number of elements with measured abundances eps(i)/
                      niso=286     !< number of isotopes/
! Z and A of alpha-enhanced isotopes
!        (charge number and mass number of alpha enhanced isotopes)
  integer, dimension(nalpha), parameter:: Zalpha=(/ 6, 8,10,12,14,16,18,20,22/),&
                                          Aalpha=(/12,16,20,24,28,32,36,40,48/)
! parameters A and B of [alpha/Fe]=A*[Fe/H]+B
  integer,dimension(nely),parameter:: &
     elz=(/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,&
          26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,44,45,46,47,48,49,50,&
          51,52,53,54,55,56,57,58,59,60,62,63,64,65,66,67,68,69,70,71,72,73,74,75,&
          76,77,78,79,80,81,82,83,90,92/)
  integer,dimension(niso),parameter:: &
     isoz=(/1, 1, 2, 2, 3, 3, 4, 5, 5, 6, 6, 7, 7, 8, 8, 8, 9,10,10,10,11,12,12,12,&
            13,14,14,14,15,16,16,16,16,17,17,18,18,18,19,19,19,20,20,20,20,20,20,21,&
            22,22,22,22,22,23,23,24,24,24,24,25,26,26,26,26,27,28,28,28,28,28,29,29,&
            30,30,30,30,30,31,31,32,32,32,32,32,33,34,34,34,34,34,34,35,35,36,36,36,&
            36,36,36,37,37,38,38,38,38,39,40,40,40,40,40,41,42,42,42,42,42,42,42,44,&
            44,44,44,44,44,44,45,46,46,46,46,46,46,47,47,48,48,48,48,48,48,48,48,49,&
            49,50,50,50,50,50,50,50,50,50,50,51,51,52,52,52,52,52,52,52,52,53,54,54,&
            54,54,54,54,54,54,54,55,56,56,56,56,56,56,56,57,57,58,58,58,58,59,60,60,&
            60,60,60,60,60,62,62,62,62,62,62,62,63,63,64,64,64,64,64,64,64,65,66,66,&
            66,66,66,66,66,67,68,68,68,68,68,68,69,70,70,70,70,70,70,70,71,71,72,72,&
            72,72,72,72,73,73,74,74,74,74,74,75,75,76,76,76,76,76,76,76,77,77,78,78,&
            78,78,78,78,79,80,80,80,80,80,80,80,81,81,82,82,82,82,83,90,92,92/),&
     isoa=(/1,2,3,4,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,&
           29,30,31,32,33,34,36,35,37,36,38,40,39,40,41,40,42,43,44,46,48,45,46,47,&
           48,49,50,50,51,50,52,53,54,55,54,56,57,58,59,58,60,61,62,64,63,65,64,66,&
           67,68,70,69,71,70,72,73,74,76,75,74,76,77,78,80,82,79,81,78,80,82,83,84,&
           86,85,87,84,86,87,88,89,90,91,92,94,96,93,92,94,95,96,97,98,100,96,98,99,&
          100,101,102,104,103,102,104,105,106,108,110,107,109,106,108,110,111,112,&
          113,114,116,113,115,112,114,115,116,117,118,119,120,122,124,121,123,120,&
          122,123,124,125,126,128,130,127,124,126,128,129,130,131,132,134,136,133,&
          130,132,134,135,136,137,138,138,139,136,138,140,142,141,142,143,144,145,&
          146,148,150,144,147,148,149,150,152,154,151,153,152,154,155,156,157,158,&
          160,159,156,158,160,161,162,163,164,165,162,164,166,167,168,170,169,168,&
          170,171,172,173,174,176,175,176,174,176,177,178,179,180,180,181,180,182,&
          183,184,186,185,187,184,186,187,188,189,190,192,191,193,190,192,194,195,&
          196,198,197,196,198,199,200,201,202,204,203,205,204,206,207,208,209,232,&
          235,238/)

  real(kind=8):: xh                   !< hydrogen mass fraction
  real(kind=8):: xhsol                ! solar hydrogen mass fraction
  real(kind=8):: xfesol               ! solar iron mass fraction
  real(kind=8):: met                  ! target metallicity Z
  real(kind=8), dimension(nalpha):: A,B
  real(kind=8), dimension(nely),parameter:: &! elemental number abundance from input in astronomical log scale
     elab_AnGr89=(/12.00d0,10.99d0,3.27d0,1.42d0,2.88d0,8.56d0,8.05d0,8.93d0,4.48d0,&
                    8.09d0,6.32d0,7.58d0,6.48d0,7.55d0,5.57d0,7.27d0,5.27d0,6.56d0,&
                    5.13d0,6.35d0,3.09d0,4.94d0,4.02d0,5.68d0,5.53d0,7.51d0,4.91d0,&
                    6.25d0,4.27d0,4.65d0,3.13d0,3.63d0,2.37d0,3.35d0,2.63d0,3.23d0,&
                    2.41d0,2.93d0,2.23d0,2.61d0,1.40d0,1.96d0,1.83d0,1.10d0,1.70d0,&
                    1.24d0,1.76d0,0.82d0,2.14d0,1.05d0,2.24d0,1.15d0,2.23d0,1.12d0,&
                    2.21d0,1.20d0,1.61d0,0.78d0,1.47d0,0.97d0,0.54d0,1.07d0,0.33d0,&
                    1.15d0,0.50d0,0.95d0,0.13d0,0.95d0,0.12d0,0.73d0,0.13d0,0.68d0,&
                    0.27d0,1.38d0,1.37d0,1.68d0,0.83d0,1.09d0,0.82d0,2.05d0,0.71d0,&
                    0.08d0,-0.49d0/), &
     elab_GrNo93=(/12.00d0,10.99d0,1.16d0,1.15d0,2.60d0,8.55d0,7.97d0,8.87d0,4.56d0,&
                    8.08d0,6.33d0,7.58d0,6.47d0,7.55d0,5.45d0,7.21d0,5.50d0,6.52d0,&
                    5.12d0,6.36d0,3.17d0,5.02d0,4.00d0,5.67d0,5.39d0,7.50d0,4.92d0,&
                    6.25d0,4.27d0,4.65d0,3.13d0,3.63d0,2.37d0,3.35d0,2.63d0,3.23d0,&
                    2.41d0,2.93d0,2.23d0,2.61d0,1.40d0,1.96d0,1.83d0,1.10d0,1.70d0,&
                    1.24d0,1.77d0,0.82d0,2.14d0,1.05d0,2.24d0,1.15d0,2.23d0,1.12d0,&
                    2.21d0,1.20d0,1.61d0,0.78d0,1.50d0,1.01d0,0.54d0,1.07d0,0.33d0,&
                    1.15d0,0.50d0,0.95d0,0.13d0,0.95d0,0.12d0,0.73d0,0.13d0,0.68d0,&
                    0.27d0,1.38d0,1.37d0,1.68d0,0.83d0,1.09d0,0.82d0,1.95d0,0.71d0,&
                    0.27d0,-0.49d0/), &
     elab_asplund05=(/12.00d0,10.93d0,1.05d0,1.38d0,2.70d0,8.39d0,7.78d0,8.66d0,&
                       4.56d0,7.84d0,6.17d0,7.53d0,6.37d0,7.51d0,5.36d0,7.14d0,&
                       5.50d0,6.18d0,5.08d0,6.31d0,3.05d0,4.90d0,4.00d0,5.64d0,&
                       5.39d0,7.45d0,4.92d0,6.23d0,4.21d0,4.60d0,2.88d0,3.58d0,&
                       2.29d0,3.33d0,2.56d0,3.28d0,2.60d0,2.92d0,2.21d0,2.59d0,&
                       1.42d0,1.92d0,1.84d0,1.12d0,1.69d0,0.94d0,1.77d0,1.60d0,&
                       2.00d0,1.00d0,2.19d0,1.51d0,2.27d0,1.07d0,2.17d0,1.13d0,&
                       1.58d0,0.71d0,1.45d0,1.01d0,0.52d0,1.12d0,0.28d0,1.14d0,&
                       0.51d0,0.93d0,0.00d0,1.08d0,0.06d0,0.88d0,-0.17d0,1.11d0,&
                       0.23d0,1.45d0,1.38d0,1.64d0,1.01d0,1.13d0,0.90d0,2.00d0,&
                       0.65d0,0.06d0,-0.52d0/), &
     elab_Asp05Cun06=(/12.00d0,10.93d0,1.05d0,1.38d0,2.70d0,8.39d0,7.78d0,8.66d0,&
                        4.56d0,8.11d0,6.17d0,7.53d0,6.37d0,7.51d0,5.36d0,7.14d0,&
                        5.50d0,6.18d0,5.08d0,6.31d0,3.05d0,4.90d0,4.00d0,5.64d0,&
                        5.39d0,7.45d0,4.92d0,6.23d0,4.21d0,4.60d0,2.88d0,3.58d0,&
                        2.29d0,3.33d0,2.56d0,3.28d0,2.60d0,2.92d0,2.21d0,2.59d0,&
                        1.42d0,1.92d0,1.84d0,1.12d0,1.69d0,0.94d0,1.77d0,1.60d0,&
                        2.00d0,1.00d0,2.19d0,1.51d0,2.27d0,1.07d0,2.17d0,1.13d0,&
                        1.58d0,0.71d0,1.45d0,1.01d0,0.52d0,1.12d0,0.28d0,1.14d0,&
                        0.51d0,0.93d0,0.00d0,1.08d0,0.06d0,0.88d0,-0.17d0,1.11d0,&
                        0.23d0,1.45d0,1.38d0,1.64d0,1.01d0,1.13d0,0.90d0,2.00d0,&
                        0.65d0,0.06d0,-0.52d0/), &
     elab_AGSS09=(/12.00d0,10.93d0,1.05d0,1.38d0,2.70d0,8.43d0,7.83d0,8.69d0,&
                    4.56d0,7.93d0,6.24d0,7.60d0,6.45d0,7.51d0,5.41d0,7.12d0,&
                    5.50d0,6.40d0,5.03d0,6.34d0,3.15d0,4.95d0,3.93d0,5.64d0,&
                    5.43d0,7.50d0,4.99d0,6.22d0,4.19d0,4.56d0,3.04d0,3.65d0,&
                    2.30d0,3.34d0,2.54d0,3.25d0,2.52d0,2.87d0,2.21d0,2.58d0,&
                    1.46d0,1.88d0,1.75d0,0.91d0,1.57d0,0.94d0,1.71d0,0.80d0,&
                    2.04d0,1.01d0,2.18d0,1.55d0,2.24d0,1.08d0,2.18d0,1.10d0,&
                    1.58d0,0.72d0,1.42d0,0.96d0,0.52d0,1.07d0,0.30d0,1.10d0,&
                    0.48d0,0.92d0,0.10d0,0.84d0,0.10d0,0.85d0,-0.12d0,0.85d0,&
                    0.26d0,1.40d0,1.38d0,1.62d0,0.92d0,1.17d0,0.90d0,1.75d0,&
                    0.65d0,0.02d0,-0.54d0/)
  double precision, dimension(niso),parameter:: &
     isoperc100=(/99.99806d0,0.00194d0,0.016597d0,99.983403d0,7.589d0,92.411d0,&
                 100.d0,19.82d0,80.18d0,98.8922d0,1.1078d0,99.6337d0,0.3663d0,&
                 99.7628d0,0.0372d0,0.20004d0,100.d0,92.9431d0,0.2228d0,6.8341d0,&
                 100.d0,78.992d0,10.003d0,11.005d0,100.d0,92.22968d0,4.68316d0,&
                 3.08716d0,100.d0,95.018d0,0.75d0,4.215d0,0.017d0,75.771d0,24.229d0,&
                 84.5946d0,15.3808d0,0.0246d0,93.25811d0,0.01167d0,6.73022d0,&
                 96.941d0,0.647d0,0.135d0,2.086d0,0.004d0,0.187d0,100.d0,8.249d0,&
                 7.437d0,73.72d0,5.409d0,5.185d0,0.2497d0,99.7503d0,4.3452d0,&
                 83.7895d0,9.5006d0,2.3647d0,100.d0,5.845d0,91.754d0,2.119d0,0.282d0,&
                 100.d0,68.0769d0,26.2231d0,1.1399d0,3.6345d0,0.9256d0,69.174d0,&
                 30.826d0,48.63d0,27.9d0,4.10d0,18.75d0,0.62d0,60.1079d0,39.8921d0,&
                 21.234d0,27.662d0,7.717d0,35.943d0,7.444d0,100.d0,0.889d0,9.366d0,&
                 7.635d0,23.772d0,49.607d0,8.731d0,50.686d0,49.314d0,0.362d0,2.33d0,&
                 11.65d0,11.55d0,56.90d0,17.21d0,72.1654d0,27.8346d0,0.5551d0,&
                 9.8168d0,7.3771d0,82.2510d0,100.d0,51.452d0,11.223d0,17.146d0,&
                 17.38d0,2.799d0,100.d0,14.8362d0,9.2466d0,15.9201d0,16.6756d0,&
                 9.5551d0,24.1329d0,9.6335d0,5.542d0,1.8688d0,12.7579d0,12.5985d0,&
                 17.06d0,31.5519d0,18.621d0,100.d0,1.02d0,11.14d0,22.33d0,27.33d0,&
                 26.46d0,11.72d0,51.8392d0,48.1608d0,1.25d0,0.89d0,12.49d0,12.8d0,&
                 24.13d0,12.22d0,28.73d0,7.49d0,4.288d0,95.712d0,0.971d0,0.659d0,&
                 0.339d0,14.536d0,7.676d0,24.223d0,8.585d0,32.593d0,4.629d0,5.789d0,&
                 57.213d0,42.787d0,0.096d0,2.603d0,0.908d0,4.816d0,7.139d0,18.952d0,&
                 31.687d0,33.799d0,100.d0,0.129d0,0.112d0,2.23d0,27.46d0,4.38d0,&
                 21.80d0,26.36d0,9.66d0,7.87d0,100.d0,0.1058d0,0.1012d0,2.417d0,&
                 6.592d0,7.853d0,11.232d0,71.699d0,0.09017d0,99.90983d0,0.186d0,&
                 0.251d0,88.449d0,11.114d0,100.d0,27.16d0,12.19d0,23.83d0,8.30d0,&
                 17.17d0,5.74d0,5.62d0,3.0734d0,14.9934d0,11.2406d0,13.8189d0,&
                 7.3796d0,26.7421d0,22.752d0,47.81d0,52.19d0,0.2029d0,2.1809d0,&
                 14.7998d0,20.4664d0,15.6518d0,24.8347d0,21.8635d0,100.d0,0.056d0,&
                 0.096d0,2.34d0,18.91d0,25.51d0,24.9d0,28.19d0,100.d0,0.137d0,&
                 1.609d0,33.61d0,22.93d0,26.79d0,14.93d0,100d0,0.13d0,3.04d0,14.28d0,&
                 21.83d0,16.13d0,31.83d0,12.76d0,97.416d0,2.584d0,0.1620d0,5.2502d0,&
                 18.5973d0,27.2840d0,13.6225d0,35.0840d0,0.0123d0,99.9877d0,0.1198d0,&
                 26.4985d0,14.3136d0,30.6422d0,28.4259d0,37.398d0,62.602d0,0.0198d0,&
                 1.5922d0,1.644d0,13.2865d0,16.1992d0,26.3438d0,40.9142d0,37.272d0,&
                 62.728d0,0.013634d0,0.78266d0,32.967d0,33.83156d0,25.24166d0,&
                 7.16349d0,100.d0,0.15344d0,9.968d0,16.873d0,23.096d0,13.181d0,&
                 29.863d0,6.865d0,29.524d0,70.476d0,1.9820d0,18.7351d0,20.5900d0,&
                 58.6929d0,100.d0,100.d0,0.72d0,99.2745d0/)
  real(kind=8), dimension(niso):: isoabsol !< solar isotopic number abundance

  private
  public :: inichem
  public :: idefaut,mainnam,xx,zini,znew
  public :: elemZ,elemA

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine inichem
! last modifications 13.12.2008 Urs Frischknecht - No, it was no Friday!

    implicit none

    integer:: i,j,&
              formatx, &   !< output formatting mode
              readstatus,& !< fortran format read status
              n1, &        !< number of elements read from elemental input file
              alpha, &     !< alpha-enhaced mode (1) or solar scaled (0)
              source, &    !< input file choice mode
              check        !< used to distinguish between alpha-enhanced/not-enhaced isotopes in a loop
    integer, dimension(5)::   selectzalu,selectaalu
    integer, dimension(26)::   selectz,selecta
    integer, dimension(15)::  mainz,maina
    integer, dimension(19)::  opalz

    real(kind=8), dimension(3)::    psol     !< contains X,Y,Z for protosolar case
    real(kind=8), dimension(19)::   elab,avmass
    real(kind=8), dimension(nely):: ely      !< elemental number abundance from input in astronomical log scale
    real(kind=8), dimension(niso):: isoab, & !< isotopic abundance
                                    isoperc  !< isotopic abundance
    real(kind=8):: normX, &                  !< X which is need due to astronomical log scaled input
                   sumX,sumX2, &             !< sum of all mass fractions before and after renormalisation
                   zold, &                   !< metallicity Z of the input abundance pattern
                   scalf, &                  !< scale factor Z_new/Z_start
                   xfe, &                    !< new iron mass fraction for alpha-enhanced case
                   afe, &                    !< "mean" mass number <A(Fe)> = sum(isoperc(i)*A(i)) over all i iron isotopes
                   zlimit, &                 !< beyond this metallicity is the alpha-enhanced plateau
                   FeHlimit, &               !< specifies where the plateau starts
                   xfelimit, &               !< iron mass fraction at FeHlimit
                   sumelab, &
                   dydz                      !< dY/dZ helium mass fraction change with metallicity - dY/dZ out of yprim and protosolar abundancess

    character(len=2)::    elname(95),elnam    !< elemental name vector / elemental name for writing output
    character(len=4), dimension(3):: fnend    !< output file endings
    character(len=4), dimension(4):: sourceid !< source identifier used in output file name

    character(*),parameter:: &
        headeralu='## abondances initiales du reseau alu pour Z= '
    character(len=5)::  elnam1      !< elemental name for writing output
    character(len=7)::  metname     !< metallicity for filename creation
    character(len=20):: outputfile  !< output file name


! primordial He mass fraction - derived from WMAP (Cyburt et al. 2003)
    real(kind=8),parameter:: yprim= 0.2484d0

! default settings
! zini is read in makeini
    met = zini

! source identifier
    sourceid=(/'AG89','GN93','As05','As09'/)
!-0.370 Mishenina et al. 2000 for oxygen
    A=(/-0.562d0,-0.886d0,-0.50d0,-0.411d0,-0.307d0,-0.435d0,-0.30d0,-0.222d0,-0.251d0/)
    B=0.d0   ! 0.2368

    source=3
    alpha=0
    formatx=2
    write(*,*)'Default settings are:'
    write(*,*)'         - Asplund-Cunha abundances'
    write(*,*)'         - scaled solar'
    write(*,*)'         - Geneva format'
    write(*,*)'         - structure from pre-calculated model'
    write(*,*)'Is it ok (yes:1, no:0)'
    read(5,*) idefaut
    if(idefaut/=0 .and. idefaut/=1) then
      stop 'Answer should be 0 or 1: aborting...'
    else if(idefaut ==0) then
! choose input file
      write(*,*) 'Choose source'
      write(*,*)'Anders & Grevesse 1989: 1 / Grevesse and Noels 1993: 2 / Asplund 2005: 3 / Asplund 2009: 4'
      read(5,*) source
      write(*,*) sourceid(int(source))
! choose alpha enhanced or solar scaled composition
      write(*,*)'scaled solar abundances (0) or alpha-enhanced (1)?'
      read(5,*) alpha
      write(*,*) 'alpha-enhanced:',alpha
      if(alpha/=0 .and. alpha/=1) then
        stop 'wrong choice - program stopped!'
      endif
! choose output format
      write(*,*)'choose format! [1=basnet/2=GENEC/3=PPN input format]'
      read(5,*) formatx
      if(formatx/=1.and.formatx/=2.and.formatx/=3) then
        stop 'wrong format choice'
      endif
    endif  ! idefaut

    select case (source)
! psol contains the proto-solar abundances
! normX is the actual solar X(H)
      case (1)
        ely=elab_AnGr89
        psol=(/0.680d0,0.30d0,0.020d0/)
        normX=0.70683d0 ! absolute value of H  - X(H) for Anders and Grevesse 1989
        zsol=0.020d0
        iopac=1
      case (2)
        ely=elab_GrNo93
        psol=(/0.680d0,0.30d0,0.020d0/)
        normX=0.70332d0 ! absolute value of H  - X(H) for Grevesse and Noels 1993
        zsol=0.020d0
        iopac=1
      case (3)
        ely=elab_Asp05Cun06
        psol=(/0.720d0,0.266d0,0.014d0/)
        normX=0.73990d0 ! absolute value of H  - X(H) for Asplund 2005 + Ne Cunha 2006
        zsol=0.0140d0
        iopac=3
      case (4)
        ely=elab_AGSS09
        psol=(/0.7198d0,0.2655d0,0.01470d0/)
        normX=0.73810d0 ! absolute value of H  - X(H) for Asplund+ 2009
        zsol=0.01420d0
        iopac=5
      case default
        stop 'wrong choice - stop program'
    end select
    ! ATTENTION: on devrait utiliser He4 et non He(tot) pour determiner dY/dZ
    ! protosolar abundances from Lodders 2003
    !   psol=(/0.7110,0.2741,0.0149/)
    dydz = (psol(2)-yprim)/psol(3)
    write(*,*) 'dY/dZ: ',dydz
!---------------------------------------------------------------------------------
    n1=size(elz)
    n2=size(isoa)

!calculate isotopic number abundances
    isoperc=isoperc100/100.d0

! protosolar abundances from Geneva Group calculated out of Lodders 2003
    sumX=0.d0
    xfe=0.d0
    do i=1,n2
     do j=1,n1
      if (isoz(i)==elz(j)) then
        isoab(i)=10.d0**(ely(j)-12.d0)*isoperc(i)*normX
        sumX=sumX+isoab(i)*isoa(i)
      endif
     enddo
    enddo
    write(*,*) 'sum of mass fractions:',sumX
    write(*,*) 'normalise the sum of the mass fraction to one'

! normalise the total mass fraction to one
    sumX2=0.d0
    xfesol=0.d0
    afe=0.d0
    do i=1,n2
     isoab(i)=isoab(i)/sumX
     isoabsol(i)=isoab(i) ! abundance for solar mixture
     sumX2=sumX2+isoab(i)*isoa(i)
     if (isoz(i)>2) then
         zold=zold+isoab(i)*isoa(i)
     endif
     if (isoz(i)== 26) then
       xfesol=xfesol+isoab(i)*isoa(i) !solar Fe mass fraction
       afe=afe + isoperc(i)*isoa(i)
     endif
    enddo
    write(*,*) 'normalisation check:',sumX2

!-----------------------------------------------------------------------
! calculate the new abundances for different than sol metallicity
!-----------------------------------------------------------------------
    scalf=1.d0
    znew=0.d0
    if(alpha/=1) then
      scalf= met/zold
    endif
!scale abundances (for Z>2) to the new metallicity
    do i=1,n2
     if (isoz(i)>2) then
       isoab(i)=isoab(i)*scalf
       znew=znew+isoab(i)*isoa(i)
     endif
    enddo

! calculate new X and Y
!----------------------------------------------------------------------------
! solar mixture
    if(alpha== 0) then
! ATTENTION: si on utilise He4 et non He(tot) pour dY/dZ, il faut enlever
! la correction du percentile isotopique
      isoab(3)=(yprim+dydz*met)*isoperc(3)/isoa(3)
      isoab(4)=(yprim+dydz*met)*isoperc(4)/isoa(4)
      isoab(1)=isoperc(1)*(1.d0-isoab(3)*isoa(3)-isoab(4)*isoa(4)-znew)
      isoab(2)=isoperc(2)*(1.d0-isoab(3)*isoa(3)-isoab(4)*isoa(4)-znew)/isoa(2)
!----------------------------------------------------------------------------
! alpha enhanced abundances
    else
      iopac=iopac+1
      xhsol= isoabsol(1)+isoabsol(2)*isoa(2)
      isoab(3)=(yprim+dydz*met)*isoperc(3)/isoa(3)
      isoab(4)=(yprim+dydz*met)*isoperc(4)/isoa(4)
      isoab(1)=isoperc(1)*(1.d0-isoab(3)*isoa(3)-isoab(4)*isoa(4)-met)
      isoab(2)=isoperc(2)*(1.d0-isoab(3)*isoa(3)-isoab(4)*isoa(4)-met)/isoa(2)
      xh= isoab(1)+isoab(2)*isoa(2)

! calculate the limit of Z for which the plateau should start
      FeHlimit=-1.d0
      xfelimit = 10.d0**(FeHlimit)*xfesol*xh/xhsol
      zlimit=func(xfelimit)+met

! check if metallicity in region of plateau
      if(met<zlimit)then
        do j=1,nalpha
         B(j)=B(j)+FeHlimit*A(j)
         A(j)=0.d0
        enddo
      endif

!calculate Xfe with the secant method
      xfe=xfesol ! initial guess
      xfe=rtsec(xfe,1.d-25,1.d-25)
      write(*,'(a4,1pe13.6)') 'xfe=',xfe

! calculate new abundances out of Xfe
      do i=1,n2
       check=0
       if(isoz(i)>2 .and. isoz(i)/=26) then
         do j=1,nalpha
!calculate abundances of alpha-enhanced isotopes
          if(isoa(i)== Aalpha(j) .and. isoz(i)==Zalpha(j)) then
            isoab(i)=isoabsol(i)*xfe/xfesol*(10.d0**B(j))*((xfe*xhsol/xh/xfesol)**A(j))
            check=1
          endif
         enddo
! calculate abundances for Z>2 except those of alpha-enhanced isotopes and Fe
         if (check==0) isoab(i)=isoabsol(i)*xfe/xfesol
! calculate abundances of Fe-isotopes
       elseif(isoz(i)==26) then
         isoab(i)=xfe/afe*isoperc(i)
       endif
      enddo
    endif

    write(*,'(a42,3(1pe13.6))')' final sum of mass fractions, metallicity:',&
                               sumX2, znew
    if(alpha==1) then
      print*,'[Fe/H]=',log10(xfe/xfesol*xhsol/xh),'[O/Fe]=',&
              log10(isoab(14)/isoabsol(14)*xfesol/xfe)
    endif
!------------------------------------------------------------------------------------------
!write number abundances/mass fractions in the output file

    sumX2=0.d0
    znew=0.d0
    do i=1,n2
     sumX2=sumX2+isoab(i)*isoa(i)
     if(isoz(i)>2) znew=znew+isoab(i)*isoa(i)
    enddo

! write to output file
    write(metname,'(1pe7.1)') znew
    fnend =(/'.bas','.gva','.ppn'/)
    elname=(/'n ','h ','he','li','be','b ','c ','n ','o ','f ','ne','na','mg','al','si',&
             'p ','s ','cl','ar','k ','ca','sc','ti','v ','cr','mn','fe','co','ni',&
             'cu','zn','ga','ge','as','se','br','kr','rb','sr','y ','zr','nb','mo',&
             'tc','ru','rh','pd','ag','cd','in','sn','sb','te','i ','xe','cs','ba',&
             'la','ce','pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb',&
             'lu','hf','ta','w ','re','os','ir','pt','au','hg','tl','pb','bi','po',&
             'at','rn','fr','ra','ac','th','pa','u ','np','pu'/)

    outputfile='iniab'//metname//sourceid(source)//fnend(formatx)
    write(*,*) 'output file: ',outputfile

    select case (formatx)
! write basnet output
      case (1)
        open(10,file=outputfile)
        do i=1,n2
         write(10,'(2i6,2(1pe14.6))') isoz(i),isoa(i),isoab(i),isoab(i)*isoa(i)
        enddo
        close(10)

! write Genec output
      case(2)
        open(11,file='netdef.in')
        write(11,'(2a7)') '### Z= ',metname
        write(11,'(e21.15)') isoab(26)*isoa(26)
        write(11,'(a83)') '### above: initial Mg25 abundance (to calculate the neutrons lost in Ne22(a,n)Mg25)'
        write(11,'(a37)') '#name,Z, A,   xabun of the elements'
        selectz=(/ 0,14,14,15,16,16,17,18,18,19,20,20,22,22,24,24,24,26,26,26,26,26,27,27,27,28 /)
        selecta=(/ 1,28,30,31,32,34,35,36,38,39,40,42,44,46,48,50,56,52,53,54,55,56,55,56,57,56 /)

        ! selectz=(/ 0,14,16,18,20,22,24,24,26,26,26,26,26,27,27,27,28 /)
        ! selecta=(/ 1,28,32,36,40,44,48,56,52,53,54,55,56,55,56,57,56 /)
        netd: do i=1,26
               do j=1,n2
                if( selectz(i)==isoz(j) .and. selecta(i)==isoa(j) ) then
                  write(*,*), "ADAM",elname(selectz(i)+1),selectz(i),selecta(i),isoab(j),isoa(j)
                  write(11,'(a2,2i4,1pe23.15)') elname(isoz(j)+1),isoz(j),isoa(j),isoab(j)*isoa(j)
                  cycle netd
                endif
               enddo
               write(11,'(a2,2i4,1pe23.15)') elname(selectz(i)+1),selectz(i),selecta(i),0.d0
              enddo netd
        close(11)
        write(6,*) 'netdef.in is done!'

        mainz=(/1,2,2, 6, 6, 7, 7, 8, 8, 8,10,10,12,12,12/)
        maina=(/1,3,4,12,13,14,15,16,17,18,20,22,24,25,26/)
        znew=1.d0
        main: do i=1,15
               do j=1,n2
                if( mainz(i)==isoz(j) .and. maina(i)==isoa(j) ) then
                  xx(i)=isoab(j)*isoa(j)
                  if(i==1) xx(i)=xx(i)+isoab(j+1)*isoa(j+1)
                  write(*,'(a6,1x,a4,1pe21.15,a1)') mainnam(i),'=50*',xx(i),','
                  znew=znew-xx(i)
                  cycle main
                endif
               enddo
              enddo main
        write(*,*) 'znew in inichem : ', znew

        open(12,file='netalu.dat')
        write(12,'(a,a7)') headeralu,metname
        selectzalu=(/  9,10,11,13,14 /)
        selectaalu=(/ 19,21,23,27,28 /)
        netalu: do i=1,5
                 do j=1,n2
                  if(selectzalu(i)==isoz(j) .and. selectaalu(i)==isoa(j)) then
                    write(12,'(a6,1pe23.15)') mainnam(i+15),isoab(j)*isoa(j)
                    xx(i+15)=isoab(j)*isoa(j)
                    cycle netalu
                  endif
                 enddo
                enddo netalu
        write(6,*) 'netalu.dat is done!'
        close(12)

! write PPN format
      case (3)
        open(10,file=outputfile)
        do i=1,n2
         if(isoz(i)==0) then
           elnam1='NEUT '
           write(10,'(i3,x,a5,x,1pe24.10)') isoz(i),elnam1,isoab(i)*isoa(i)
         elseif(isoz(i)==1.and.isoa(i)==1) then
           elnam1='PROT '
           write(10,'(i3,x,a5,x,1pe24.10)') isoz(i),elnam1,isoab(i)*isoa(i)
         else
           elnam=elname(isoz(i))
           write(10,'(i3,x,a2,i3,x,1pe24.10)') isoz(i),elnam,isoa(i),isoab(i)*isoa(i)
         endif
        enddo
        close(10)
      case default
        stop 'Bad format choice'
    end select

    return

  end subroutine inichem

!***********************************************************************************************
!> find root of func with the secant method
  real(kind=8) function rtsec(x1,x2,xacc)
    implicit none
    integer, parameter:: maxit=200 !maximum number of iterations.
    integer:: j
    real(kind=8):: x1,x2,xacc
    real(kind=8):: dx,f,fl,swap,xl
!using the secant method, find the root of a function func thought to lie between x1 and x2.
!The root, returned as rtsec, is refined until its accuracy is ±xacc.
    fl=func(x1)
    f=func(x2)
    if(abs(fl)<abs(f))then !pick the bound with the smaller function value as the most
      rtsec=x1 !recent guess.
      xl=x2
      swap=fl
      fl=f
      f=swap
    else
      xl=x1
      rtsec=x2
    endif
    do j=1,maxit !secant loop.
     dx=(xl-rtsec)*f/(f-fl) !increment with respect to latest value.
     xl=rtsec
     fl=f
     rtsec=rtsec+dx
     !print*,j,rtsec,func(rtsec)
     if(rtsec<0.d0) rtsec=abs(rtsec)
     f=func(rtsec)
     if(abs(dx)<xacc.or.f==0.d0)return !convergence.
     if(f-fl==0.d0) return
    enddo
    print*, rtsec
    stop 'rtsec exceed maximum iterations'
  end function rtsec

!***********************************************************************************************
!function of which the root is searched
  real(kind=8) function func(x)
    implicit none
    real(kind=8),intent(in) :: 	x   ! iron mass fraction
    integer:: i, j   ! runtime variables

    func=0.d0
    main: do i=1,n2
           if (isoz(i)>2) then
             do j=1,nalpha
              if(isoa(i)== Aalpha(j) .and. isoz(i)==Zalpha(j)) then
                func=func + isoabsol(i)*isoa(i)*x/xfesol*((x*xhsol/xh/xfesol)**A(j))&
                      *(10.d0**B(j))
                cycle main
              endif
             enddo
             func = func + isoabsol(i)*isoa(i)*x/xfesol
           endif
          enddo main

    func=func-met

  end function func

end module inichemmod
