module inputparam

  implicit none

! Definition de kindreal, pour eviter de devoir ajouter le module evol:
  integer, parameter:: kindreal = 8
  integer,parameter:: imagn_default=0,ianiso_default=0,ipop3_default=0,ibasnet_default=0,iopac_default=3,&
    ikappa_default=5,istati_default=0,igamma_default=0,nndr_default=1,iledou_default=0,idifcon_default=0,&
    iover_default=1,iunder_default=0,nbchx_default=200,nrband_default=1,icncst_default=0,iprn_default=10,&
    iout_default=0,itmin_default=5,idebug_default=0,itests_default=0,tauH_fit_default=1,RSG_Mdot_default=0,&
    n_mag_default=1,nsmooth_default=1,end_at_phase_default=4,end_at_model_default=0,iprezams_default=1,&
    n_snap_default=10
  real(kindreal),parameter:: fenerg_default=1.0d0,richac_default=1.0d0,zsol_default=1.40d-2,frein_default=0.0d0,&
    K_Kawaler_default=0.d0,Omega_saturation_default=14.d0,vwant_default=0.0d0,xfom_default=1.0d0, &
    dunder_default=0.0d0,dgro_default=0.010d0,dgr20_default=0.010d0,binm2_default=0.d0,periodini_default=0.d0,&
    B_initial_default=0.d0,add_diff_default=0.0d0,Be_mdotfrac_default=0.0d0,start_mdot_default=0.80d0,&
    alpha_F_default=1.d0
  logical,parameter:: xyfiles_default=.false.,bintide_default=.false.,const_per_default=.true.,&
    var_rates_default=.false.,verbose_default=.false.,Add_Flux_default=.true.,&
    diff_only_default=.false.,stop_deg_default=.true.,SupraEddMdot_default=.true.,qminsmooth_default=.false., &
    superv_default=.false.

! NAMELISTS VARIABLES
! **** Model characteristics
  integer,save:: nwseq,modanf,nzmod,end_at_phase=end_at_phase_default,end_at_model=end_at_model_default
  character(256),save:: starname
!-----------------------------------------------------------------------
  namelist /CharacteristicsParams/starname,nwseq,modanf,nzmod,end_at_phase,end_at_model
!-----------------------------------------------------------------------

! **** Physical inputs
  integer,save:: irot,isol,imagn=imagn_default,ialflu,ianiso=ianiso_default,ipop3=ipop3_default,&
      ibasnet=ibasnet_default,phase,iprezams=iprezams_default
  real(kindreal),save:: binm2=binm2_default,periodini=periodini_default
  logical,save:: var_rates=var_rates_default,bintide=bintide_default,const_per=const_per_default
!-----------------------------------------------------------------------
  namelist /PhysicsParams/irot,isol,imagn,ialflu,ianiso,ipop3,ibasnet,phase,var_rates,&
                          bintide,binm2,periodini,const_per,iprezams
!-----------------------------------------------------------------------

! **** Chemical composition
  integer,save:: iopac=iopac_default,ikappa=ikappa_default
  real(kindreal),save:: zinit,zsol=zsol_default,z
!-----------------------------------------------------------------------
  namelist /CompositionParams/zinit,zsol,z,iopac,ikappa
!-----------------------------------------------------------------------

! **** Rotation-linked parameters
  integer,save:: idiff,iadvec,istati=istati_default,icoeff,igamma=igamma_default,idialo,idialu, &
                 n_mag=n_mag_default,nsmooth=nsmooth_default
  real(kindreal),save:: fenerg=fenerg_default,richac=richac_default,frein=frein_default,K_Kawaler=K_Kawaler_default, &
                          Omega_saturation=Omega_saturation_default,rapcrilim,vwant=vwant_default, &
                          xfom=xfom_default,omega,xdial,B_initial=B_initial_default,add_diff=add_diff_default, &
                          alpha_F=alpha_F_default
  logical,save:: Add_Flux=Add_Flux_default,diff_only=diff_only_default,qminsmooth=qminsmooth_default
!-----------------------------------------------------------------------
  namelist /RotationParams/idiff,iadvec,istati,icoeff,fenerg,richac,igamma,frein,K_Kawaler,Omega_saturation,rapcrilim, &
                           vwant,xfom,omega,xdial,idialo,idialu,Add_Flux,diff_only,B_initial,add_diff, &
                           n_mag,alpha_F,nsmooth,qminsmooth
!-----------------------------------------------------------------------

! **** Surface parameters
  integer,save:: imloss,ifitm,nndr=nndr_default,RSG_Mdot=RSG_Mdot_default
  real(kindreal),save:: fmlos,fitm,fitmi,deltal,deltat,fitmi_default,Be_mdotfrac=Be_mdotfrac_default,start_mdot=start_mdot_default
  logical,save:: SupraEddMdot=SupraEddMdot_default
!-----------------------------------------------------------------------
  namelist /SurfaceParams/imloss,fmlos,RSG_Mdot,SupraEddMdot,ifitm,fitm,fitmi,deltal,deltat,nndr,Be_mdotfrac,start_mdot
!-----------------------------------------------------------------------

! **** Convection-linked parameters
  integer,save:: iledou=iledou_default,idifcon=idifcon_default,my,iover=iover_default,iunder=iunder_default
  real(kindreal),save:: elph,dovhp,dunder=dunder_default
!-----------------------------------------------------------------------
  namelist /ConvectionParams/iledou,idifcon,iover,elph,my,dovhp,iunder,dunder
!-----------------------------------------------------------------------

! **** Convergence-linked parameters
  integer,save:: nbchx=nbchx_default,nrband=nrband_default
  real(kindreal),save:: gkorm,alph,agdr,faktor,dgrp,dgrl,dgry,dgrc,dgro=dgro_default,dgr20=dgr20_default
!-----------------------------------------------------------------------
  namelist /ConvergenceParams/gkorm,alph,agdr,faktor,dgrp,dgrl,dgry,dgrc,dgro,dgr20,nbchx,nrband
!-----------------------------------------------------------------------

! **** Timestep controle
  integer,save:: islow,icncst=icncst_default,tauH_fit=tauH_fit_default
  real(kindreal),save:: xcn
!-----------------------------------------------------------------------
  namelist /TimeControle/xcn,islow,icncst,tauH_fit
!-----------------------------------------------------------------------

! **** Other controles
  integer,save:: iauto,iprn=iprn_default,iout=iout_default,itmin=itmin_default,&
      idebug=idebug_default,itests=itests_default,n_snap=n_snap_default
  logical,save:: display_plot,xyfiles=xyfiles_default,verbose=verbose_default,&
      stop_deg=stop_deg_default,superv=superv_default
!-----------------------------------------------------------------------
  namelist /VariousSettings/display_plot,iauto,iprn,iout,itmin,xyfiles,idebug,&
      itests,verbose,stop_deg,n_snap,superv
!-----------------------------------------------------------------------

  public
!  private:: Write_param_int,Write_param_real,Write_param_logical
  private:: imagn_default,ianiso_default,ipop3_default,ibasnet_default,iopac_default,ikappa_default,istati_default,&
    igamma_default,nndr_default,iledou_default,iover_default,iunder_default,nbchx_default,nrband_default, &
    icncst_default,iprn_default,iout_default,itmin_default,fenerg_default,richac_default,zsol_default, &
    frein_default,K_Kawaler_default,Omega_saturation_default,vwant_default,xfom_default,dunder_default,dgr20_default, &
    xyfiles_default,idebug_default,bintide_default,binm2_default,periodini_default,const_per_default, &
    var_rates_default,verbose_default,stop_deg_default,tauH_fit_default,SupraEddMdot_default,RSG_Mdot_default,&
    Be_mdotfrac_default,start_mdot_default,n_snap_default,superv_default

contains
!=======================================================================
subroutine Write_namelist(Unit,nwseqnew,modanfnew,nzmodnew,xcnwant)
!-----------------------------------------------------------------------
  implicit none

  integer,intent(in):: Unit,nwseqnew,modanfnew,nzmodnew
  real(kindreal),intent(in):: xcnwant
!-----------------------------------------------------------------------
  write(Unit,'(a)') "&CharacteristicsParams"
  write(Unit,'(1x,a,a)') "starname=","'"//trim(starname)//"'"
  write(Unit,'(1x,a,i0)') "nwseq=",nwseqnew
  write(Unit,'(1x,a,i0)') "modanf=",modanfnew
  write(Unit,'(1x,a,i0)') "nzmod=",nzmodnew
  write(Unit,'(1x,a,i0)') "end_at_phase=",end_at_phase
  write(Unit,'(1x,a,i0)') "end_at_model=",end_at_model
  write(Unit,'("&END"/)')

  write(Unit,'(a)') "&PhysicsParams"
  write(Unit,'(1x,2(a,i0))') "irot=",irot,", isol=",isol
  write(Unit,'(1x,a,i0)') "imagn=",imagn
  write(Unit,'(1x,a,i0)') "ialflu=",ialflu
  write(Unit,'(1x,a,i0)') "ianiso=",ianiso
  if (abs(zinit) < epsilon(0.d0)) then
      write(Unit,'(1x,a,i0)') "ipop3=1"
  else
      write(Unit,'(1x,a,i0)') "ipop3=0"
  endif
  write(Unit,'(1x,a,i0)') "ibasnet=",ibasnet
  write(Unit,'(1x,a,i0)') "phase=",phase
  write(Unit,'(1x,a,i0)') "iprezams=",iprezams
  write(Unit,'(1x,a,l2)') "var_rates=",var_rates
  write(Unit,'(1x,a,l2)') "bintide=",bintide
  write(Unit,'(1x,a,es9.2)') "binM2=",binm2
  write(Unit,'(1x,a,es9.2)') "periodini=",periodini
  write(Unit,'(1x,a,l2)') "const_per=",const_per
  write(Unit,'("&END"/)')

  write(Unit,'(a)') "&CompositionParams"
  write(Unit,'(1x,a,es9.2)') "zinit=",zinit
  write(Unit,'(1x,a,d10.3)') "zsol=",zsol
  write(Unit,'(1x,a,d21.15)') "z=",z
  write(Unit,'(1x,a,i0)') "iopac=",iopac
  write(Unit,'(1x,a,i0)') "ikappa=",ikappa
  write(Unit,'("&END"/)')

  write(Unit,'(a)') "&RotationParams"
  write(Unit,'(1x,2(a,i0))') "idiff=",idiff,", iadvec=",iadvec
  write(Unit,'(1x,a,i0)') "istati=",istati
  write(Unit,'(1x,a,i0)') "icoeff=",icoeff
  write(Unit,'(1x,a,d10.3)') "fenerg=",fenerg
  write(Unit,'(1x,a,d10.3)') "richac=",richac
  write(Unit,'(1x,a,i0)') "igamma=",igamma
  write(Unit,'(1x,a,d10.3)') "frein=",frein
  write(Unit,'(1x,a,d10.3)') "K_Kawaler=",K_Kawaler
  write(Unit,'(1x,a,d10.3)') "Omega_saturation=",Omega_saturation
  write(Unit,'(1x,a,f7.5)') "rapcrilim=",rapcrilim
  write(Unit,'(1x,a,d10.3)') "vwant=",vwant
  write(Unit,'(1x,a,d10.3)') "xfom=",xfom
  write(Unit,'(1x,a,es21.15)') "omega=",omega
  write(Unit,'(1x,a,f5.3,2(a,i0))') "xdial=",xdial,", idialo=",idialo,", idialu=",idialu
  write(Unit,'(1x,a,l1)') "Add_Flux=",Add_Flux
  write(Unit,'(1x,a,l1)') "diff_only=",diff_only
  write(Unit,'(1x,a,d10.3)') "B_initial=",B_initial
  write(Unit,'(1x,a,d10.3)') "add_diff=",add_diff
  write(Unit,'(1x,a,i0)') "n_mag=",n_mag
  write(Unit,'(1x,a,d10.3)') "alpha_F=",alpha_F
  write(Unit,'(1x,a,i0)') "nsmooth=",nsmooth
  write(Unit,'(1x,a,l1)') "qminsmooth=",qminsmooth
  write(Unit,'("&END"/)')

  if (irot > 0) then
    fitmi_default = 0.9990d0
  else
    fitmi_default = 0.980d0
  endif
  write(Unit,'(a)') "&SurfaceParams"
  write(Unit,'(1x,a,i0,a,d10.3)') "imloss=",imloss,", fmlos=",fmlos
  write(Unit,'(1x,a,i0,1x,a,l2)') "RSG_Mdot=",RSG_Mdot,"SupraEddMdot=",SupraEddMdot
  write(Unit,'(1x,a,f4.2)') "Be_mdotfrac=",Be_mdotfrac
  write(Unit,'(1x,a,f4.2)')"start_mdot=",start_mdot
  write(Unit,'(1x,a,i0,a,f11.9,a,f11.9)') "ifitm=",ifitm,", fitmi=",fitmi_default,", fitm=",fitmi_default
  write(Unit,'(1x,2(a,f7.5))') "deltal=",deltal,", deltat=",deltat
  write(Unit,'(1x,a,i0)') "nndr=",nndr
  write(Unit,'("&END"/)')

  write(Unit,'(a)') "&ConvectionParams"
  write(Unit,'(1x,a,i0)') "iledou=",iledou
  write(Unit,'(1x,a,i0)') "idifcon=",idifcon
  write(Unit,'(1x,a,f0.3,a,i0)') "elph=",elph,", my=",my
  write(Unit,'(1x,a,i0,a,f5.3)') "iover=",iover,", dovhp=",dovhp
  write(Unit,'(1x,a,i0)') "iunder=",iunder
  write(Unit,'(1x,a,f5.3)') "dunder=",dunder
  write(Unit,'("&END"/)')

  write(Unit,'(a)') "&ConvergenceParams"
  write(Unit,'(1x,a,f0.3,a,f5.3)') "gkorm=",gkorm,", alph=",alph
  write(Unit,'(1x,a,d9.2,a,es9.2)') "agdr=",agdr,", faktor=",faktor
  write(Unit,'(1x,2(a,f6.4))') "dgrp=",dgrp,", dgrl=",dgrl
  write(Unit,'(1x,3(a,f7.5))') "dgry=",dgry,", dgrc=",dgrc,", dgro=",dgro
  write(Unit,'(1x,a,d10.3)') "dgr20=",dgr20
  write(Unit,'(1x,a,i0)') "nbchx=",nbchx
  write(Unit,'(1x,a,i0)') "nrband=",nrband
  write(Unit,'("&END"/)')

  write(Unit,'(a)') "&TimeControle"
  write(Unit,'(1x,a,i0,a,f0.3)') "islow=",islow,", xcn=",xcnwant
  write(Unit,'(1x,a,i0)') "icncst=",icncst
  write(Unit,'(1x,a,i0)') "tauH_fit=",tauH_fit
  write(Unit,'("&END"/)')

  write(Unit,'(a)') "&VariousSettings"
  write(Unit,'(1x,2(a,l2))') "display_plot=",display_plot
  write(Unit,'(1x,a,i2)') "iauto=",iauto
  write(Unit,'(1x,a,i0)') "n_snap=",n_snap
  write(Unit,'(1x,a,i0)') "iprn=",iprn
  write(Unit,'(1x,a,i0)') "iout=",iout
  write(Unit,'(1x,a,i0)') "itmin=",itmin
  write(Unit,'(1x,a,l1)') "superv=",superv
  write(Unit,'(1x,a,l1)') "xyfiles=",xyfiles
  write(Unit,'(1x,a,i0)') "idebug=",idebug
  write(Unit,'(1x,a,i0)') "itests=",itests
  write(Unit,'(1x,a,l1)') "verbose=",verbose
  write(Unit,'(1x,a,l1)') "stop_deg=",stop_deg
  write(Unit,'("&END")')

  return

end subroutine Write_namelist
!=======================================================================
subroutine Ask_changes
!-----------------------------------------------------------------------
  implicit none

  integer:: Change_params,Category_change,Temp_Var_Int
  real(kindreal):: Temp_Var_real
  character:: answer,Temp_Var_char
!-----------------------------------------------------------------------
  Category_change = 99
  Change_params = 99
  answer = ''

  do while (answer /= 'y' .and. answer /= 'n')
    write(*,*) 'Do you want to change some input parameters ? (y)es (n)o '
    read(5,*) answer
    if (answer /= 'y' .and. answer /= 'n') then
      write(*,*) 'Please type y or n...'
    endif
  enddo

  if (answer == 'y') then
    Category_change = 99
    do while (Category_change /= 0)
      write(*,*) '------------------------------'
      write(*,*) '*** CATEGORIES ***'
      write(*,*) '  1: CHARACTERISTICS inputs'
      write(*,*) '  2: PHYSICS inputs'
      write(*,*) '  3: ROTATION inputs'
      write(*,*) '  4: SURFACE inputs'
      write(*,*) '  5: CONVECTION inputs'
      write(*,*) '  6: CONVERGENCE inputs'
      write(*,*) '  7: TIME CONTROL inputs'
      write(*,*) '  8: VARIOUS SETTINGS inputs'
      write(*,*) '------------------------------'
      write(*,*) 'Enter the category number (0 to skip or exit):'
      read(5,*) Category_change
      select case(Category_change)
      case (0)
        write(*,*) 'No more changes...'
      case (1) ! *** change of CHARACTERISTICS inputs
        Change_params = 99
        do while (Change_params /= 0)
          write(*,*) '------------------------------'
          write(*,*) '*** CHARACTERISTICS inputs ***'
          write(*,'(a,i5)') ' 1: nzmod        :',nzmod
          write(*,'(a,i2)') ' 2: end_at_phase :',end_at_phase
          write(*,'(a,i5)') ' 3: end_at_model :',end_at_model
          write(*,*) '------------------------------'
          write(*,*) 'Parameters to change (0 to skip or exit):'
          read(5,*) Change_params
          select case (Change_params)
          case (0)
            write(*,*) 'No more changes of CHARACTERISTICS parameters'
          case (1)
            Temp_Var_Int = 99999
            write(*,*) 'Enter the desired value for nzmod (default 1000):'
            read(5,*) Temp_Var_Int
            nzmod = Temp_Var_Int
          case (2)
            Temp_Var_Int = 99
            do while (Temp_Var_Int<2 .or. Temp_Var_Int>7)
              write(*,*) 'Possible values for END_AT_PHASE'
              write(*,*) '------------------------------'
              write(*,*) ' 2: stop at the end of H-b'
              write(*,*) ' 3: stop at the end of He-b'
              write(*,*) ' 4: stop at the end of C-b'
              write(*,*) ' 5: stop at the end of Ne-b'
              write(*,*) ' 6: stop at the end of O-b'
              write(*,*) '------------------------------'
              write(*,*)'Enter the desired value (default 4):'
              read(5,*) Temp_Var_Int
            enddo
            end_at_phase = Temp_Var_Int
          case (3)
            Temp_Var_Int = 0
            write(*,*) 'Enter the desired value for end_at_model (default 0)'
            read(5,*) Temp_Var_Int
            end_at_model = Temp_Var_Int
          case default
            write(*,*) 'Wrong number, should be 0,1,2, or 3'
          end select ! end CHARACTERISTICS inputs selection
        enddo
      case (2) ! *** change of PHYSICS inputs
        Change_params = 99
        do while (Change_params /= 0)
          write(*,*) '------------------------------'
          write(*,*) '*** PHYSICS inputs ***'
          write(*,'(a,i2)') ' 1: isol     :',isol
          write(*,'(a,i2)') ' 2: imagn    :',imagn
          write(*,'(a,i2)') ' 3: ialflu   :',ialflu
          write(*,'(a,i2)') ' 4: ianiso   :',ianiso
          write(*,'(a,l2)') ' 5: var_rates:',var_rates
          write(*,'(a,l2)') ' 6: bintide  :',bintide
          write(*,'(a,l2)') ' 7: const_per:',const_per
          write(*,'(a,f7.3)') ' 8: binM2    :',binM2
          write(*,'(a,f7.3)') ' 9: periodini:',periodini
          write(*,*) '------------------------------'
          write(*,*) 'Parameters to change (0 to skip or exit):'
          read(5,*) Change_params
          select case (Change_params)
          case (0)
            write(*,*) 'No more changes of PHYSICS parameters'
          case (1)
            if (irot == 1) then
              Temp_Var_Int = 99
              do while (Temp_Var_Int/=0 .and. Temp_Var_Int/=1 .and. Temp_Var_Int/=2)
                write(*,*) 'Enter the desired value for isol (0,1,2):'
                read(5,*) Temp_Var_Int
              enddo
              isol = Temp_Var_Int
            else
              write(*,*) 'You defined a non-rotating star, you should not touch isol.'
            endif
          case (2)
            Temp_Var_Int=99
            do while (Temp_Var_Int/=0 .and. Temp_Var_Int/=1)
              write(*,*) 'Enter the desired value for imagn (0,1):'
              read(5,*) Temp_Var_Int
            enddo
            imagn = Temp_Var_Int
          case (3)
            Temp_Var_Int = 99
            do while (Temp_Var_Int/=0 .and. Temp_Var_Int/=1)
              write(*,*)'Enter the desired value for ialflu (0,1):'
              read(5,*) Temp_Var_Int
            enddo
            ialflu = Temp_Var_Int
          case (4)
            Temp_Var_Int = 99
            do while (Temp_Var_Int/=0 .and. Temp_Var_Int/=1)
              write(*,*)'Enter the desired value for ianiso (0,1):'
              read(5,*) Temp_Var_Int
            enddo
            ianiso = Temp_Var_Int
          case (5)
            Temp_Var_char = ''
            do while (Temp_Var_char/='t' .and. Temp_Var_char/='f' &
                 .and. Temp_Var_char/='T' .and. Temp_Var_char/= 'F')
              write(*,*)'Enter the desired value for var_rates (T/F):'
              read(5,*) Temp_Var_char
            enddo
            if (Temp_Var_char=='t' .or. Temp_Var_char=='T') then
              var_rates = .true.
            elseif (Temp_Var_char=='f' .or. Temp_Var_char=='F') then
              var_rates = .false.
            endif
          case(6)
            Temp_Var_char = ''
            do while (Temp_Var_char/='t' .and. Temp_Var_char/='f' &
                 .and. Temp_Var_char/='T' .and. Temp_Var_char/= 'F')
              write(*,*)'Enter the desired value for bintide (T/F):'
              read(5,*) Temp_Var_char
            enddo
            if (Temp_Var_char=='t' .or. Temp_Var_char=='T') then
              bintide = .true.
            elseif (Temp_Var_char=='f' .or. Temp_Var_char=='F') then
              bintide = .false.
            endif
          case (7)
            if (bintide) then
              Temp_Var_char = ''
              do while (Temp_Var_char/='t' .and. Temp_Var_char/='f' &
                   .and. Temp_Var_char/='T' .and. Temp_Var_char/= 'F')
                write(*,*)'Enter the desired value for const_per (T/F):'
                read(5,*) Temp_Var_char
              enddo
              if (Temp_Var_char=='t' .or. Temp_Var_char=='T') then
                const_per = .true.
              elseif (Temp_Var_char=='f' .or. Temp_Var_char=='F') then
                const_per = .false.
              endif
            else
              write(*,*) 'bintide is set to F, you should not touch const_per'
            endif
          case (8)
            if (bintide) then
              Temp_Var_real = -2.d0
              do while (Temp_Var_real < 0.d0)
                write(*,*)'Enter the desired value for binM2 (in Msol):'
                read(5,*) Temp_Var_real
              enddo
              binM2 = Temp_Var_real
            else
              write(*,*) 'bintide is set to F, you should not touch binM2'
            endif
          case (9)
            if (bintide) then
              Temp_Var_real = -2.d0
              do while (Temp_Var_real < 0.d0)
                write(*,*)'Enter the desired value for periodini (in days):'
                read(5,*) Temp_Var_real
              enddo
              periodini = Temp_Var_real
            else
              write(*,*) 'bintide is set to F, you should not touch periodini'
            endif
          case default
            write(*,*) 'Wrong number, should be an integer between 0 and 9'
          end select ! end PHYSICS inputs selection
        enddo
      case (3) ! *** change of ROTATION inputs
        Change_params = 99
        do while (Change_params /= 0)
          write(*,*) '*** ROTATION inputs ***'
          write(*,*) '------------------------------'
          write(*,'(a,i2)') ' 1: icoeff   :',icoeff
          write(*,'(a,i2)') ' 2: istati   :',istati
          write(*,'(a,d11.5)') ' 3: frein    :',frein
          write(*,'(a,d11.5)') ' 4: K_Kawaler:',K_Kawaler
          write(*,'(a,l2)') ' 5: Add_Flux :',Add_Flux
          write(*,'(a,d11.5)') ' 6: B_initial:',B_initial
          write(*,'(a,d11.5)') ' 7: add_diff :',add_diff
          write(*,'(a,i2)') ' 8: n_mag    :',n_mag
          write(*,'(a,f7.3)') ' 9: alpha_F  :',alpha_F
          write(*,'(a,i2)') '10: nsmooth  :',nsmooth
          write(*,'(a,l2)') '11: qminsmooth:',qminsmooth
          write(*,*) '------------------------------'
          write(*,*) 'Parameters to change (0 to skip or exit):'
          read(5,*) Change_params
          select case (Change_params)
          case (0)
            write(*,*) 'No more changes of ROTATION parameters'
          case (1)
            Temp_Var_Int=99
            do while (Temp_Var_Int/=11 .and. Temp_Var_Int/=12 .and. Temp_Var_Int/=13 &
                .and. Temp_Var_Int/=21 .and. Temp_Var_Int/=22 .and. Temp_Var_Int/=23 &
                .and. Temp_Var_Int/=31 .and. Temp_Var_Int/=32 .and. Temp_Var_Int/=33)
              write(*,*) 'Possible values for ICOEFF'
              write(*,*) ' ___________________________________________________'
              write(*,*) '|     Dshear      |               Dh                |'
              write(*,*) '|                 | Zahn 92 | Maeder 03 | Mathis 04 |'
              write(*,*) '|---------------------------------------------------|'
              write(*,*) '| Maeder 97       |    11   |     12    |    13     |'
              write(*,*) '| Talon & Zahn 97 |    21   |     22    |    23     |'
              write(*,*) '| Dtot Maeder 13  |    31   |     32    |    33     |'
              write(*,*) '|___________________________________________________|'
              write(*,*) 'Enter the desired value (default 11):'
              read(5,*) Temp_Var_Int
            enddo
            icoeff = Temp_Var_Int
          case (2)
            Temp_Var_Int = 99
            do while (Temp_Var_Int/=0 .and. Temp_Var_Int/=1)
              write(*,*)'Enter the desired value for istati (0,1):'
              read(5,*) Temp_Var_Int
            enddo
            istati = Temp_Var_Int
          case (3)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*)'Enter the desired value for frein (in Gauss):'
              read(5,*) Temp_Var_real
            enddo
            frein = Temp_Var_real
          case (4)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*)'Enter the desired value for K_Kawaler:'
              read(5,*) Temp_Var_real
            enddo
            K_Kawaler = Temp_Var_real
          case (5)
            Temp_Var_char = ''
            do while (Temp_Var_char/='t' .and. Temp_Var_char/='f' &
                 .and. Temp_Var_char/='T' .and. Temp_Var_char/= 'F')
              write(*,*)'Enter the desired value for Add_Flux (T/F):'
              read(5,*) Temp_Var_char
            enddo
            if (Temp_Var_char=='t' .or. Temp_Var_char=='f') then
              Add_Flux = .true.
            elseif (Temp_Var_char=='0' .or. Temp_Var_char=='F') then
              Add_Flux = .false.
            endif
          case (6)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*)'Enter the desired value for B_initial (in Gauss):'
              read(5,*) Temp_Var_real
            enddo
            B_initial = Temp_Var_real
          case(7)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*)'Enter the desired value for add_diff:'
              read(5,*) Temp_Var_real
            enddo
            add_diff = Temp_Var_real
          case(8)
            Temp_Var_Int = 99
            do while (Temp_Var_Int/= 1 .and. Temp_Var_Int/=2 .and. Temp_Var_Int/=3)
              write(*,*) 'Possible values for N_MAG'
              write(*,*) '------------------------------'
              write(*,*) ' 1: pure Taylor-Spruit (2002A&A...381..923S)'
              write(*,*) ' 2: modified TS (Geneva group development)'
              write(*,*) ' 3: Fuller+ 2019 modified TS (2019MNRAS.485.3661F)'
              write(*,*) '------------------------------'
              write(*,*)'Enter the desired value (default 1):'
              read(5,*) Temp_Var_Int
            enddo
            n_mag = Temp_Var_Int
            if (n_mag == 3) then
              write(*,*) 'With this settings we advise you to change nsmooth=5'
            endif
          case(9)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*)'Enter the desired value for alpha_F (default 1.0):'
              read(5,*) Temp_Var_real
            enddo
            alpha_F = Temp_Var_real
          case(10)
            Temp_Var_Int = 99
            do while (Temp_Var_Int > 20)
              write(*,*)'Recommended values for NSMOOTH:'
              write(*,*) '------------------------------'
              write(*,*) ' 1: default value'
              write(*,*) ' 5: Fuller+ 2019 implementation (n_mag=3)'
              write(*,*) '------------------------------'
              write(*,*)'Enter the desired value:'
              read(5,*) Temp_Var_Int
            enddo
            nsmooth = Temp_Var_Int
          case (11)
            Temp_Var_char = ''
            do while (Temp_Var_char/='t' .and. Temp_Var_char/='f' &
                 .and. Temp_Var_char/='T' .and. Temp_Var_char/= 'F')
              write(*,*)'Enter the desired value for qminsmooth (T/F):'
              read(5,*) Temp_Var_char
            enddo
            if (Temp_Var_char=='t' .or. Temp_Var_char=='f') then
              qminsmooth = .true.
            elseif (Temp_Var_char=='0' .or. Temp_Var_char=='F') then
              qminsmooth = .false.
            endif
          case default
            write(*,*) 'Wrong number, should be an integer between 0 and 11'
          end select ! end ROTATION inputs selection
        enddo
      case(4) ! *** change of SURFACE inputs
        Change_params = 99
        do while (Change_params /= 0)
          write(*,*) '------------------------------'
          write(*,*) '*** SURFACE inputs ***'
          write(*,'(a,i2)') ' 1: imloss        :',imloss
          write(*,'(a,d12.5)') ' 2: fmlos         :',fmlos
          write(*,'(a,i2)') ' 3: RSG_Mdot      :',RSG_Mdot
          write(*,'(a,l2)') ' 4: SupraEddMdot  :',SupraEddMdot
          write(*,'(a,f6.2)') ' 5: Be_Mdotfrac   :',Be_mdotfrac
          write(*,'(a,f6.2)') ' 6: start_mdot    :',start_mdot
          write(*,'(a,f9.5)') ' 7: fitm          :',fitm
          write(*,'(a,i2)') ' 8: ifitm         :',ifitm
          write(*,*) '------------------------------'
          write(*,*) 'Parameters to change (0 to skip or exit):'
          read(5,*) Change_params
          select case (Change_params)
          case (0)
            write(*,*) 'No more changes of SURFACE parameters'
          case (1)
            Temp_Var_Int = 99
            do while (Temp_Var_Int>=13)
              write(*,*) 'Possible values for IMLOSS'
              write(*,*) '------------------------------'
              write(*,*) ' 1: de Jager+ 1988'
              write(*,*) ' 2: mass loss in Msol/yr given by FMLOS'
              write(*,*) ' 3: Reimers formula with etaR given by FMLOS'
              write(*,*) ' 4: WR mass loss : as in papier V'
              write(*,*) ' 5: Kudritzki & Puls 2000'
              write(*,*) ' 6: Vink+ 2001'
              write(*,*) ' 9: Kudritzki 2002'
              write(*,*) '11: Vink+ 2001 modified by Markova & Puls 2008 + priv. comm. Puls (nov. 2010)'
              write(*,*) '12: Gormaz-Matamala+ 2022'
              write(*,*) '------------------------------'
              write(*,*) 'Enter the desired value (default 6):'
              read(5,*) Temp_Var_Int
            enddo
            imloss = Temp_Var_Int
          case (2)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*)'Enter the desired value for fmlos:'
              read(5,*) Temp_Var_real
            enddo
            fmlos = Temp_Var_real
          case (3)
            Temp_Var_Int = 99
            do while (Temp_Var_Int>=5)
              write(*,*) 'Possible values for RSG_MDOT'
              write(*,*) '------------------------------'
              write(*,*) ' 0: Sylvester (1998), van Loon 1999 (cf Crowther 2001)'
              write(*,*) ' 1: de Jager+ 1988'
              write(*,*) ' 2: Beasor & Davies 2020'
              write(*,*) ' 3: Kee+ 2021'
              write(*,*) ' 4: van Loon+ 2005'
              write(*,*) '------------------------------'
              write(*,*) 'Enter the desired value (default 0):'
              read(5,*) Temp_Var_Int
            enddo
            RSG_Mdot = Temp_Var_Int
          case (4)
            Temp_Var_char = ''
            do while (Temp_Var_char/='t' .and. Temp_Var_char/='f' &
                 .and. Temp_Var_char/='T' .and. Temp_Var_char/= 'F')
              write(*,*)'Enter the desired value for SupraEddMdot (T/F):'
              read(5,*) Temp_Var_char
            enddo
            if (Temp_Var_char=='t' .or. Temp_Var_char=='T') then
              SupraEddMdot = .true.
            elseif (Temp_Var_char=='f' .or. Temp_Var_char=='F') then
              SupraEddMdot = .false.
            endif
          case (5)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*) 'Enter the desired value for Be_Mdotfrac (recommended 0.1):'
              read(5,*) Temp_Var_real
            enddo
            Be_mdotfrac = Temp_Var_real
          case (6)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*) 'Enter the desired value for start_mdot (recommended 0.8):'
              read(5,*) Temp_Var_real
            enddo
            start_mdot = Temp_Var_real
          case (7)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*) 'Enter the desired value for fitm:'
              write(*,*) '     default: 0.980 for irot=0'
              write(*,*)'               0.999 for irot=1'
              read(5,*) Temp_Var_real
            enddo
            fitm = Temp_Var_real
            fitmi = fitm
          case (8)
            Temp_Var_Int = 99
            do while (Temp_Var_Int>=7)
              write(*,*) 'Possible values for IFITM'
              write(*,*) '------------------------------'
              write(*,*) ' 0: manual change'
              write(*,*) ' 1: automatic change, check for fitm being inside the outer CZ'
              write(*,*) ' 2: automatic change with polynomial fit (not if irot=1)'
              write(*,*) ' 3: automatic change, smoothly at each time step (default)'
              write(*,*) ' 4: automatic change, every 10 timesteps'
              write(*,*) ' 5: automatic change, only upwards'
              write(*,*) ' 6: automatic change, only downwards'
              write(*,*) '------------------------------'
              write(*,*) 'Enter the desired value:'
              read(5,*) Temp_Var_Int
            enddo
            ifitm = Temp_Var_Int
          case default
            write(*,*) 'Wrong number, should be an integer between 0 and 8'
          end select ! end SURFACE inputs selection
        enddo
      case (5) ! *** change of CONVECTION inputs
        Change_params = 99
        do while (Change_params /= 0)
          write(*,*) '------------------------------'
          write(*,*) '*** CONVECTION inputs ***'
          write(*,'(a,i2)') ' 1: iledou   :',iledou
          write(*,'(a,i2)') ' 2: iover    :',iover
          write(*,'(a,f7.3)') ' 3: dovhp    :',dovhp
          write(*,'(a,i2)') ' 4: iunder   :',iunder
          write(*,'(a,f7.3)') ' 5: dunder   :',dunder
          write(*,*) '------------------------------'
          write(*,*) 'Parameters to change (0 to skip or exit):'
          read(5,*) Change_params
          select case (Change_params)
          case (0)
            write(*,*) 'No more changes of CONVECTION parameters'
          case (1)
            Temp_Var_Int = 99
            do while (Temp_Var_Int/=0 .and. Temp_Var_Int/=1)
              write(*,*)'Enter the desired value for iledou (0,1):'
              read(5,*) Temp_Var_Int
            enddo
            iledou = Temp_Var_Int
          case (2)
            Temp_Var_Int = 99
            do while (Temp_Var_Int/=0 .and. Temp_Var_Int/=1 .and. Temp_Var_Int/=2)
              write(*,*) 'Possible values for IOVER'
              write(*,*) '------------------------------'
              write(*,*) ' 0: no overshoot'
              write(*,*) ' 1: fixed overshoot set by DOVHP'
              write(*,*) ' 2: variable overshoot from Baraffe+ 2023'
              write(*,*)'Enter the desired value for iover (0,1,2):'
              read(5,*) Temp_Var_Int
            enddo
            iover = Temp_Var_Int
          case (3)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*)'Enter the desired value for dovhp:'
              read(5,*) Temp_Var_real
            enddo
            dovhp = Temp_Var_real
          case (4)
            Temp_Var_Int = 99
            do while (Temp_Var_Int/=0 .and. Temp_Var_Int/=1)
              write(*,*)'Enter the desired value for iunder (0,1):'
              read(5,*) Temp_Var_Int
            enddo
            iunder = Temp_Var_Int
          case (5)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*)'Enter the desired value for dunder:'
              read(5,*) Temp_Var_real
            enddo
            dunder = Temp_Var_real
          case default
            write(*,*) 'Wrong number, should be an integer between 0 and 5'
          end select ! end CONVECTION inputs selection
        enddo
      case (6) ! *** change of CONVERGENCE inputs
        Change_params = 99
        do while (Change_params /= 0)
          write(*,*) '------------------------------'
          write(*,*) '*** CONVERGENCE inputs ***'
          write(*,*) 'Parameters to change (0 to skip or exit):'
          write(*,'(a,f7.3)') ' 1: gkorm   :',gkorm
          write(*,'(a,f7.3)') ' 2: alph    :',alph
          write(*,'(a,d11.5)') ' 3: faktor  :',faktor
          write(*,'(a,d11.5)') ' 4: dgrp    :',dgrp
          write(*,'(a,d11.5)') ' 5: dgrl    :',dgrl
          write(*,'(a,d11.5)') ' 6: dgry    :',dgry
          write(*,*) '------------------------------'
          read(5,*) Change_params
          select case (Change_params)
          case (0)
            write(*,*) 'No more changes of PHYSICS parameters'
          case (1)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*)'Enter the desired value for gkorm:'
              read(5,*) Temp_Var_real
            enddo
            gkorm = Temp_Var_real
          case (2)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*)'Enter the desired value for alph:'
              read(5,*) Temp_Var_real
            enddo
            alph = Temp_Var_real
          case (3)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*)'Enter the desired value for faktor:'
              read(5,*) Temp_Var_real
            enddo
            faktor = Temp_Var_real
          case (4)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*)'Enter the desired value for dgrp (default 0.010):'
              read(5,*) Temp_Var_real
            enddo
            dgrp = Temp_Var_real
          case (5)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*)'Enter the desired value for dgrl (default 0.010):'
              read(5,*) Temp_Var_real
            enddo
            dgrl = Temp_Var_real
          case (6)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*)'Enter the desired value for dgry (default 0.003):'
              read(5,*) Temp_Var_real
            enddo
            dgry = Temp_Var_real
          case default
            write(*,*) 'Wrong number, should be an integer between 0 and 6'
          end select ! end CONVERGENCE inputs selection
        enddo
      case (7) ! *** change of TIME CONTROLE inputs
        Change_params = 99
        do while (Change_params /= 0)
          write(*,*) '------------------------------'
          write(*,*) '*** TIME CONTROLE inputs ***'
          write(*,'(a,i2)') ' 1: icncst   :',icncst
          write(*,'(a,i2)') ' 2: tauH_fit :',tauH_fit
          write(*,*) '------------------------------'
          write(*,*) 'Parameters to change (0 to skip or exit):'
          read(5,*) Change_params
          select case (Change_params)
          case (0)
            write(*,*) 'No more changes of TIME CONTROLE parameters'
          case (1)
            Temp_Var_Int = 99
            do while (Temp_Var_Int/=0 .and. Temp_Var_Int/=1)
              write(*,*)'Enter the desired value for icncst (0,1):'
              read(5,*) Temp_Var_Int
            enddo
            icncst = Temp_Var_Int
          case (2)
            Temp_Var_Int = 99
            do while (Temp_Var_Int/=1 .and. Temp_Var_Int/=2)
              write(*,*) 'Possible values for TAUH_FIT'
              write(*,*) '------------------------------'
              write(*,*) ' 1: double fit with M_trans=10 Msol (default)'
              write(*,*) ' 2: single fit on the full mass domain'
              write(*,*) '------------------------------'
              write(*,*)'Enter the desired value:'
              read(5,*) Temp_Var_Int
            enddo
            tauH_fit = Temp_Var_Int
          case default
            write(*,*) 'Wrong number, should be 0,1, or 2'
          end select ! end TIME CONTROLE inputs selection
        enddo
      case (8) ! *** change of VARIOUS SETTINGS inputs
        Change_params = 99
        do while (Change_params /= 0)
          write(*,*) '------------------------------'
          write(*,*) '*** VARIOUS SETTINGS inputs ***'
          write(*,'(a,l2)') ' 1: display_plot :',display_plot
          write(*,'(a,l2)') ' 2: verbose      :',verbose
          write(*,'(a,i2)') ' 3: iauto        :',iauto
          write(*,'(a,i5)') ' 4: n_snap       :',n_snap
          write(*,'(a,i5)') ' 5: iprn         :',iprn
          write(*,'(a,l2)') ' 6: stop_deg     :',stop_deg
          write(*,'(a,l2)') ' 7: superv       :',superv
          write(*,'(a,l2)') ' 8: xyfiles      :',xyfiles
          write(*,'(a,i2)') ' 9: idebug       :',idebug
          write(*,'(a,i5)') '10: itests       :',itests
          write(*,*) '------------------------------'
          write(*,*) 'Parameters to change (0 to skip or exit):'
          read(5,*) Change_params
          select case (Change_params)
          case (0)
            write(*,*) 'No more changes of TIME CONTROLE parameters'
          case (1)
            Temp_Var_char = ''
            do while (Temp_Var_char/='t' .and. Temp_Var_char/='f' &
                 .and. Temp_Var_char/='T' .and. Temp_Var_char/= 'F')
              write(*,*)'Enter the desired value for display_plot (T/F):'
              read(5,*) Temp_Var_char
            enddo
            if (Temp_Var_char=='t' .or. Temp_Var_char=='T') then
              display_plot = .true.
            elseif (Temp_Var_char=='f' .or. Temp_Var_char=='F') then
              display_plot = .false.
            endif
          case (2)
            Temp_Var_char = ''
            do while (Temp_Var_char/='t' .and. Temp_Var_char/='f' &
                 .and. Temp_Var_char/='T' .and. Temp_Var_char/= 'F')
              write(*,*)'Enter the desired value for verbose (T/F):'
              read(5,*) Temp_Var_char
            enddo
            if (Temp_Var_char=='t' .or. Temp_Var_char=='T') then
              verbose = .true.
            elseif (Temp_Var_char=='f' .or. Temp_Var_char=='F') then
              verbose = .false.
            endif
          case (3)
            Temp_Var_Int = 99
            do while (Temp_Var_Int/=0 .and. Temp_Var_Int/=1 .and. Temp_Var_Int/=2)
              write(*,*) 'Enter the desired value for iauto (0,1,2):'
              read(5,*) Temp_Var_Int
            enddo
            iauto = Temp_Var_Int
          case (4)
            Temp_Var_Int = 99999
            write(*,*) 'Enter the desired value for n_snap (default 10):'
            read(5,*) Temp_Var_Int
            n_snap = Temp_Var_Int
          case (5)
            Temp_Var_Int = 99999
            write(*,*) 'Enter the desired value for iprn (default 10):'
            read(5,*) Temp_Var_Int
            iprn = Temp_Var_Int
          case (6)
            Temp_Var_char = ''
            do while (Temp_Var_char/='t' .and. Temp_Var_char/='f' &
                 .and. Temp_Var_char/='T' .and. Temp_Var_char/= 'F')
              write(*,*)'Enter the desired value for stop_deg (T/F):'
              read(5,*) Temp_Var_char
            enddo
            if (Temp_Var_char=='t' .or. Temp_Var_char=='T') then
              stop_deg = .true.
            elseif (Temp_Var_char=='f' .or. Temp_Var_char=='F') then
              stop_deg = .false.
            endif
          case (7)
            Temp_Var_char = ''
            do while (Temp_Var_char/='t' .and. Temp_Var_char/='f' &
                 .and. Temp_Var_char/='T' .and. Temp_Var_char/= 'F')
              write(*,*)'Enter the desired value for superv (T/F):'
              read(5,*) Temp_Var_char
            enddo
            if (Temp_Var_char=='t' .or. Temp_Var_char=='T') then
              superv = .true.
            elseif (Temp_Var_char=='f' .or. Temp_Var_char=='F') then
              superv = .false.
            endif
          case (8)
            Temp_Var_char = ''
            do while (Temp_Var_char/='t' .and. Temp_Var_char/='f' &
                 .and. Temp_Var_char/='T' .and. Temp_Var_char/= 'F')
              write(*,*)'Enter the desired value for xyfiles (T/F):'
              read(5,*) Temp_Var_char
            enddo
            if (Temp_Var_char=='t' .or. Temp_Var_char=='T') then
              xyfiles = .true.
            elseif (Temp_Var_char=='f' .or. Temp_Var_char=='F') then
              xyfiles = .false.
            endif
          case (9)
            Temp_Var_Int = 99
            do while (Temp_Var_Int > 4)
              write(*,*) 'Enter the desired value for idebug (between 0 and 4):'
              read(5,*) Temp_Var_Int
            enddo
            idebug = Temp_Var_Int
          case (10)
            Temp_Var_Int = 99
            write(*,*) 'Enter the desired value for itests:'
            read(5,*) Temp_Var_Int
            itests = Temp_Var_Int
          case default
            write(*,*)'Wrong number, should be an integer between 0 and 10'
          end select ! end VARIOUS SETTINGS inputs selection
        enddo
      case default
        write(*,*) 'Wrong number, should be an integer between 0 and 7'
      end select ! category selection
    enddo
  endif

end subroutine Ask_changes
!=======================================================================
end module inputparam
