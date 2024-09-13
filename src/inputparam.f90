module inputparam

  use io_definitions
  use evol,only: ldi,kindreal
  use caramodele,only: nwmd,xmini
  use storage, only: GenecStar

  implicit none

  interface Write_param
    module procedure Write_param_int
    module procedure Write_param_real
    module procedure Write_param_logical
  end interface Write_param

  integer,parameter:: &
          imagn_default=0,&
          inetwork_default=0,&
          ieos_default=0,&
          ianiso_default=0,&
          ipop3_default=0,&
          ibasnet_default=0,&
          iopac_default=3,&
          ikappa_default=5,&
          istati_default=0,&
          igamma_default=0,&
          n_M03_default=0,&
          nndr_default=1,&
          iledou_default=0,&
          idifcon_default=0,&
          iover_default=1,&
          iunder_default=0,&
          nbchx_default=200,&
          nrband_default=1,&
          icncst_default=0,&
          iprn_default=10,&
          iout_default=0,&
          itmin_default=5,&
          idebug_default=0,&
          itests_default=0,&
          tauH_fit_default=1,&
          n_mag_default=1,&
          nsmooth_default=5,&
          end_at_phase_default=4,&
          end_at_model_default=0,&
          iprezams_default=1,&
          n_snap_default=10
  real(kindreal),parameter:: &
          fenerg_default=1.0d0,&
          richac_default=1.0d0,&
          zsol_default=1.40d-2,&
          frein_default=0.0d0,&
          K_Kawaler_default=0.d0,&
          Omega_saturation_default=14.d0,&
          vwant_default=0.0d0,&
          A_M03_default=0.d0,&
          xfom_default=1.0d0,&
          Z_dep_default=0.85d0,&
          Xs_WR_default=0.3d0,&
          dunder_default=0.0d0,&
          dgro_default=0.010d0,&
          dgr20_default=0.010d0,&
          binm2_default=0.d0,&
          periodini_default=0.d0,&
          B_initial_default=0.d0,&
          add_diff_default=0.0d0,&
          Be_mdotfrac_default=0.0d0,&
          start_mdot_default=0.80d0,&
          alpha_F_default=6.d0,&
          end_at_time_default=4.418064d17,& ! 14 billion years
          D_clump_default = 10.d0
  logical,parameter:: &
          xyfiles_default=.false.,&
          bintide_default=.false.,&
          const_per_default=.true.,&
          var_rates_default=.false.,&
          verbose_default=.false.,&
          Add_Flux_default=.true.,&
          diff_only_default=.false.,&
          stop_deg_default=.true.,&
          SupraEddMdot_default=.true.,&
          qminsmooth_default=.true.,&
          dcirch_inclusion_default=.false.,&
          superv_default=.false.,&
          hardJump_default=.true.,&
          force_prescription_default=.false.,&
          print_winds_default=.false.,&
          winds_not_applied_default=.false.,&
          prezams_winds_not_applied_default=.false.

  ! if libgenec is set to .true., no input will be asked.
  logical,save:: &
          libgenec=.false.

! VARIABLES DE LECTURE
  integer,save:: &
          idern,&
          ichem,&
          itminc

! NAMELISTS VARIABLES
! **** Model characteristics
  integer,save:: &
          nwseq,&
          modanf,&
          nzmod,&
          end_at_phase=end_at_phase_default,&
          end_at_model=end_at_model_default
  character(256),save:: &
          starname
  real(kindreal),save:: &
          end_at_time=end_at_time_default
!-----------------------------------------------------------------------
  namelist /CharacteristicsParams/starname,nwseq,modanf,nzmod,end_at_phase,end_at_model
!-----------------------------------------------------------------------

! **** Physical inputs
  integer,save:: &
          irot,&
          isol,&
          imagn=imagn_default,&
          ieos = ieos_default,&
          inetwork = inetwork_default,&
          ialflu,&
          ianiso=ianiso_default,&
          ipop3=ipop3_default,&
          ibasnet=ibasnet_default,&
          phase,&
          iprezams=iprezams_default
  real(kindreal),save:: &
          binm2=binm2_default,&
          periodini=periodini_default
  logical,save:: &
          var_rates=var_rates_default,&
          bintide=bintide_default,&
          const_per=const_per_default
!-----------------------------------------------------------------------
  namelist /PhysicsParams/irot,isol,imagn,ieos,inetwork,ialflu,ianiso,ipop3,ibasnet,phase,var_rates,bintide,binm2,&
          periodini,const_per,iprezams
!-----------------------------------------------------------------------

! **** Chemical composition
  integer,save:: &
          iopac=iopac_default,&
          ikappa=ikappa_default
  real(kindreal),save:: &
          zinit,&
          zsol=zsol_default,&
          z
!-----------------------------------------------------------------------
  namelist /CompositionParams/zinit,zsol,z,iopac,ikappa
!-----------------------------------------------------------------------

! **** Rotation-linked parameters
  integer,save:: &
          idiff,&
          iadvec,&
          istati=istati_default,&
          icoeff,&
          igamma=igamma_default,&
          idialo,&
          idialu,&
          n_M03=n_M03_default,&
          n_mag=n_mag_default,&
          nsmooth=nsmooth_default
  real(kindreal),save:: &
          fenerg=fenerg_default,&
          richac=richac_default,&
          frein=frein_default,&
          K_Kawaler=K_Kawaler_default,&
          Omega_saturation=Omega_saturation_default,&
          rapcrilim,&
          vwant=vwant_default,&
          A_M03=A_M03_default,&
          xfom=xfom_default,&
          omega,&
          xdial,&
          B_initial=B_initial_default,&
          add_diff=add_diff_default,&
          alpha_F=alpha_F_default
  logical,save:: &
          Add_Flux=Add_Flux_default,&
          diff_only=diff_only_default,&
          qminsmooth=qminsmooth_default,&
          dcirch_inclusion=dcirch_inclusion_default
!-----------------------------------------------------------------------
  namelist /RotationParams/idiff,iadvec,istati,icoeff,fenerg,richac,igamma,frein,K_Kawaler,Omega_saturation,rapcrilim, &
          vwant,xfom,omega,xdial,idialo,idialu,Add_Flux,diff_only,B_initial,add_diff,&
          n_mag,alpha_F,nsmooth,qminsmooth,dcirch_inclusion,n_M03,A_M03
!-----------------------------------------------------------------------

! **** Winds parameters
  integer,save:: &
          imloss,&
          OB_Mdot,&
          RSG_Mdot,&
          WR_Mdot,&
          Fallback_Mdot
  real(kindreal),save:: &
          fmlos,&
          Be_mdotfrac=Be_mdotfrac_default,&
          start_mdot=start_mdot_default,&
          Z_dep=Z_dep_default,&
          Xs_WR=Xs_WR_default,&
          D_clump = D_clump_default
  logical,save:: &
          SupraEddMdot=SupraEddMdot_default,&
          hardJump=hardJump_default,&
          force_prescription=force_prescription_default,&
          print_winds=print_winds_default,&
          winds_not_applied=winds_not_applied_default,&
          prezams_winds_not_applied=prezams_winds_not_applied_default
!-----------------------------------------------------------------------
  namelist /WindsParams/fmlos,OB_Mdot,RSG_Mdot,WR_Mdot,Fallback_Mdot,Z_dep,Xs_WR, &
          SupraEddMdot,Be_mdotfrac,start_mdot,hardJump,force_prescription,print_winds,&
          D_clump,winds_not_applied,prezams_winds_not_applied 
!-----------------------------------------------------------------------

! **** Surface parameters
  integer,save:: &
          ifitm,&
          nndr=nndr_default
  real(kindreal),save:: &
          fitm,&
          fitmi,&
          fitmi_default,&
          deltal,&
          deltat
!-----------------------------------------------------------------------
  namelist /SurfaceParams/ifitm,fitm,fitmi,deltal,deltat,nndr
!-----------------------------------------------------------------------

! **** Convection-linked parameters
  integer,save:: &
          iledou=iledou_default,&
          idifcon=idifcon_default,&
          my,&
          iover=iover_default,&
          iunder=iunder_default
  real(kindreal),save:: &
          elph,&
          dovhp,&
          dunder=dunder_default
!-----------------------------------------------------------------------
  namelist /ConvectionParams/iledou,idifcon,iover,elph,my,dovhp,iunder,dunder
!-----------------------------------------------------------------------

! **** Convergence-linked parameters
  integer,save:: &
          nbchx=nbchx_default,&
          nrband=nrband_default
  real(kindreal),save:: &
          gkorm,&
          alph,&
          agdr,&
          faktor,&
          dgrp,&
          dgrl,&
          dgry,&
          dgrc,&
          dgro=dgro_default,&
          dgr20=dgr20_default
!-----------------------------------------------------------------------
  namelist /ConvergenceParams/gkorm,alph,agdr,faktor,dgrp,dgrl,dgry,dgrc,dgro,dgr20,nbchx,nrband
!-----------------------------------------------------------------------

! **** Timestep controle
  integer,save:: &
          islow,&
          icncst=icncst_default,&
          tauH_fit=tauH_fit_default
  real(kindreal),save:: &
          xcn
!-----------------------------------------------------------------------
  namelist /TimeControle/xcn,islow,icncst,tauH_fit
!-----------------------------------------------------------------------

! **** Other controles
  integer,save:: &
          iauto,&
          iprn=iprn_default,&
          iout=iout_default,&
          itmin=itmin_default,&
          idebug=idebug_default,&
          itests=itests_default,&
          n_snap=n_snap_default
  logical,save:: &
          display_plot,&
          xyfiles=xyfiles_default,&
          verbose=verbose_default,&
          stop_deg=stop_deg_default,&
          superv=superv_default
!-----------------------------------------------------------------------
  namelist /VariousSettings/display_plot,iauto,iprn,iout,itmin,xyfiles,idebug,&
      itests,verbose,stop_deg,n_snap,superv
!-----------------------------------------------------------------------

  integer:: &
          isugi=1
  real(kindreal),save:: &
          xtt,&
          agds,&
          agdp,&
          agdt

  public:: Write_param
  public
  private :: xtt
  private:: Write_param_int,Write_param_real,Write_param_logical
  public:: &
          imagn_default,&
          inetwork_default,&
          ieos_default,&
          ianiso_default,&
          ipop3_default,&
          ibasnet_default,&
          iopac_default,&
          ikappa_default,&
          istati_default,&
          igamma_default,&
          n_M03_default,&
          nndr_default,&
          iledou_default,&
          iunder_default,&
          nbchx_default,&
          nrband_default,&
          icncst_default,&
          iprn_default,&
          iout_default,&
          itmin_default,&
          fenerg_default,&
          richac_default,&
          zsol_default,&
          frein_default,&
          K_Kawaler_default,&
          Omega_saturation_default,&
          vwant_default,&
          A_M03_default,&
          xfom_default,&
          iover_default,&
          dunder_default,&
          dgr20_default,&
          xyfiles_default,&
          idebug_default,&
          bintide_default,&
          binm2_default,&
          periodini_default,&
          const_per_default,&
          tauH_fit_default,&
          var_rates_default,&
          verbose_default,&
          stop_deg_default,&
          n_mag_default,&
          alpha_F_default,&
          nsmooth_default,&
          SupraEddMdot_default,&
          Be_mdotfrac_default,&
          start_mdot_default,&
          iprezams_default,&
          n_snap_default,&
          superv_default,&
          Z_dep_default,&
          Xs_WR_default,&
          hardJump_default,&
          print_winds_default,&
          dcirch_inclusion_default,&
          winds_not_applied_default,&
          prezams_winds_not_applied_default

contains
!=======================================================================
subroutine Write_param_int(Unit,n_name,n_in,n_default)
!-----------------------------------------------------------------------
  implicit none

  integer,intent(in):: Unit,n_in,n_default
  character(*),intent(in):: n_name
!-----------------------------------------------------------------------
  if ((n_in /= n_default) .or. (modanf == 0)) then
    write(Unit,'(1x,a,i0)') trim(n_name),n_in
  endif

  return

end subroutine Write_param_int
!=======================================================================
subroutine Write_param_real(Unit,x_name,x_in,x_default)
!-----------------------------------------------------------------------
  implicit none

  integer,intent(in):: Unit
  real(kindreal),intent(in):: x_in,x_default
  character(*),intent(in):: x_name
!-----------------------------------------------------------------------
  if ((x_in /= x_default) .or. (modanf == 0)) then
    write(Unit,'(1x,a,d16.9)') trim(x_name),x_in
  endif

  return

end subroutine Write_param_real
!=======================================================================
subroutine Write_param_logical(Unit,a_name,a_in,a_default)
!-----------------------------------------------------------------------
  implicit none

  integer,intent(in):: Unit
  logical,intent(in):: a_in,a_default
  character(*),intent(in):: a_name
!-----------------------------------------------------------------------
  if ((a_in .neqv. a_default) .or. (modanf == 0)) then
    write(Unit,'(1x,a,l1)') trim(a_name),a_in
  endif

  return

end subroutine Write_param_logical
!=======================================================================
subroutine Write_namelist(Unit,nwseqnew,modanfnew,nzmodnew,xcnwant)
!-----------------------------------------------------------------------
  use const,only: um

  implicit none

  integer,intent(in):: Unit,nwseqnew,modanfnew,nzmodnew
  real(kindreal),intent(in):: xcnwant
!-----------------------------------------------------------------------
  if (.not. libgenec) then
    write(Unit,'(a)') "&CharacteristicsParams"
    write(Unit,'(1x,a,a)') "starname=","'"//trim(starname)//"'"
    write(Unit,'(1x,a,i0)') "nwseq=",nwseqnew
    write(Unit,'(1x,a,i0)') "modanf=",modanfnew
    write(Unit,'(1x,a,i0)') "nzmod=",nzmodnew
    call Write_param(Unit,"end_at_phase=",end_at_phase,end_at_phase_default)
    call Write_param(Unit,"end_at_model=",end_at_model,end_at_model_default)
    write(Unit,'("&END"/)')

    write(Unit,'(a)') "&PhysicsParams"
    write(Unit,'(1x,2(a,i0))') "irot=",irot,", isol=",isol
    call Write_param(Unit,"imagn=",imagn,imagn_default)
    call Write_param(Unit,"ieos=",ieos,ieos_default)
    call Write_param(Unit,"inetwork=",inetwork,inetwork_default)
    write(Unit,'(1x,a,i0)') "ialflu=",ialflu
    call Write_param(Unit,"ianiso=",ianiso,ianiso_default)
    if (modanf == 0) then
      if (abs(zinit) < epsilon(0.d0)) then
          ipop3 = 1
      else
          ipop3 = 0
      endif
    endif
    call Write_param(Unit,"ipop3=",ipop3,ipop3_default)
    call Write_param(Unit,"ibasnet=",ibasnet,ibasnet_default)
    write(Unit,'(1x,a,i0)') "phase=",phase
    if ((modanf == 0) .and. (irot > 0)) then
      iprezams = 1
    endif
    call Write_param(Unit,"iprezams=",iprezams,iprezams_default)
    call Write_param(Unit,"var_rates=",var_rates,var_rates_default)
    call Write_param(Unit,"bintide=",bintide,bintide_default)
    if (bintide .or. modanf == 0) then
      write(Unit,'(1x,a,es9.2)') "binM2=",binm2
      write(Unit,'(1x,a,es13.6)') "periodini=",periodini
      write(Unit,'(1x,a,l2)') "const_per=",const_per
    endif
    write(Unit,'("&END"/)')

    write(Unit,'(a)') "&CompositionParams"
    write(Unit,'(1x,a,es9.2)') "zinit=",zinit
    call Write_param(Unit,"zsol=",zsol,zsol_default)
    write(Unit,'(1x,a,d21.15)') "z=",z
    call Write_param(Unit,"iopac=",iopac,iopac_default)
    call Write_param(Unit,"ikappa=",ikappa,ikappa_default)
    write(Unit,'("&END"/)')

    write(Unit,'(a)') "&RotationParams"
    write(Unit,'(1x,2(a,i0))') "idiff=",idiff,", iadvec=",iadvec
    call Write_param(Unit,"istati=",istati,istati_default)
    write(Unit,'(1x,a,i0)') "icoeff=",icoeff
    call Write_param(Unit,"fenerg=",fenerg,fenerg_default)
    call Write_param(Unit,"richac=",richac,richac_default)
    call Write_param(Unit,"igamma=",igamma,igamma_default)
    call Write_param(Unit,"frein=",frein,frein_default)
    call Write_param(Unit,"K_Kawaler=",K_Kawaler,K_Kawaler_default)
    call Write_param(Unit,"Omega_saturation=",Omega_saturation,Omega_saturation_default)
    write(Unit,'(1x,a,f8.5)') "rapcrilim=",rapcrilim
    call Write_param(Unit,"vwant=",vwant,vwant_default)
    call Write_param(Unit,"A_M03=",A_M03,A_M03_default)
    call Write_param(Unit,"n_M03=",n_M03,n_M03_default)
    call Write_param(Unit,"xfom=",xfom,xfom_default)
    if (omega < 0.d0) then
      omega = 1.d-22
    endif
    write(Unit,'(1x,a,es21.15)') "omega=",omega
    write(Unit,'(1x,a,f6.3,2(a,i0))') "xdial=",xdial,", idialo=",idialo,", idialu=",idialu
    call Write_param(Unit,"Add_Flux=",Add_Flux,Add_Flux_default)
    call Write_param(Unit,"diff_only=",diff_only,diff_only_default)
    call Write_param(Unit,"B_initial=",B_initial,B_initial_default)
    call Write_param(Unit,"add_diff=",add_diff,add_diff_default)
    call Write_param(Unit,"n_mag=",n_mag,n_mag_default)
    call Write_param(Unit,"alpha_F=",alpha_F,alpha_F_default)
    call Write_param(Unit,"nsmooth=",nsmooth,nsmooth_default)
    call Write_param(Unit,"qminsmooth=",qminsmooth,qminsmooth_default)
    call Write_param(Unit,"dcirch_inclusion=",dcirch_inclusion,dcirch_inclusion_default)
    write(Unit,'("&END"/)')

    write(Unit,'(a)') "&WindsParams"
    write(Unit,'(1x,a,i0,a,i0)') "OB_Mdot=",OB_Mdot,", RSG_Mdot=",RSG_Mdot
    write(Unit,'(1x,a,i0,a,i0)') "WR_Mdot=",WR_Mdot,", Fallback_Mdot=",Fallback_Mdot
    write(Unit,'(1x,a,d10.3)') "fmlos=",fmlos
    call Write_param(Unit,"Z_dep=",Z_dep,Z_dep_default)
    call Write_param(Unit,"Xs_WR=",Xs_WR,Xs_WR_default)
    call Write_param(Unit,"D_clump=",D_clump,D_clump_default)
    call Write_param(Unit,"SupraEddMdot=",SupraEddMdot,SupraEddMdot_default)
    call Write_param(Unit,"Be_mdotfrac=",Be_mdotfrac,Be_mdotfrac_default)
    call Write_param(Unit,"start_mdot=",start_mdot,start_mdot_default)
    call Write_param(Unit,"hardJump=",hardJump,hardJump_default)
    call Write_param(Unit,"force_prescription=",force_prescription,&
                                            force_prescription_default)
    call Write_param(Unit,"print_winds=",print_winds,print_winds_default)
      call Write_param(Unit,"winds_not_applied=",winds_not_applied,winds_not_applied_default)
    call Write_param(Unit,"prezams_winds_not_applied=",prezams_winds_not_applied,prezams_winds_not_applied_default)
    write(Unit,'("&END"/)')

    if (irot > 0) then
      fitmi_default = 0.9990d0
    else
      fitmi_default = 0.980d0
    endif
    if ((fitmi == 0.0d0) .and. modanf == 0) then
      fitmi = fitmi_default
    endif
    write(Unit,'(a)') "&SurfaceParams"
    write(Unit,'(1x,a,i0,a,f12.9)') "ifitm=",ifitm,", fitm=",fitm
    call Write_param(Unit,"fitmi=",fitmi,fitmi_default)
    write(Unit,'(1x,2(a,f8.5))') "deltal=",deltal,", deltat=",deltat
    call Write_param(Unit,"nndr=",nndr,nndr_default)
    write(Unit,'("&END"/)')

    write(Unit,'(a)') "&ConvectionParams"
    call Write_param(Unit,"iledou=",iledou,iledou_default)
    write(Unit,'(1x,a,i0)') "idifcon=",idifcon
    write(Unit,'(1x,a,f0.3,a,i0)') "elph=",elph,", my=",my
    write(Unit,'(1x,a,i0,a,f6.3)') "iover=",iover,", dovhp=",dovhp
    call Write_param(Unit,"iunder=",iunder,iunder_default)
    call Write_param(Unit,"dunder=",dunder,dunder_default)
    write(Unit,'("&END"/)')

    write(Unit,'(a)') "&ConvergenceParams"
    write(Unit,'(1x,a,f0.3,a,f6.3)') "gkorm=",gkorm,", alph=",alph
    write(Unit,'(1x,a,d9.2,a,es9.2)') "agdr=",agdr,", faktor=",faktor
    if (modanf == 0) then
      write(Unit,'(1x,2(a,f7.4))') "dgrp=",dgrp,", dgrl=",dgrl
    else
      write(Unit,'(1x,2(a,f7.4))') "dgrp=",dgrp/um,", dgrl=",dgrl/um
    endif
    write(Unit,'(1x,3(a,f8.5))') "dgry=",dgry,", dgrc=",dgrc,", dgro=",dgro
    call Write_param(Unit,"dgr20=",dgr20,dgr20_default)
    call Write_param(Unit,"nbchx=",nbchx,nbchx_default)
    call Write_param(Unit,"nrband=",nrband,nrband_default)
    write(Unit,'("&END"/)')

    write(Unit,'(a)') "&TimeControle"
    write(Unit,'(1x,a,i0,a,f0.3)') "islow=",islow,", xcn=",xcnwant
    call Write_param(Unit,"icncst=",icncst,icncst_default)
    call Write_param(Unit,"tauH_fit=",tauH_fit,tauH_fit_default)
    write(Unit,'("&END"/)')

    write(Unit,'(a)') "&VariousSettings"
    write(Unit,'(1x,2(a,l2))') "display_plot=",display_plot
    write(Unit,'(1x,a,i2)') "iauto=",iauto
    call Write_param(Unit,"n_snap=",n_snap,n_snap_default)
    call Write_param(Unit,"iprn=",iprn,iprn_default)
    call Write_param(Unit,"iout=",iout,iout_default)
    call Write_param(Unit,"itmin=",itmin,itmin_default)
    call Write_param(Unit,"superv=",superv,superv_default)
    call Write_param(Unit,"xyfiles=",xyfiles,xyfiles_default)
    call Write_param(Unit,"idebug=",idebug,idebug_default)
    call Write_param(Unit,"itests=",itests,itests_default)
    call Write_param(Unit,"verbose=",verbose,verbose_default)
    call Write_param(Unit,"stop_deg=",stop_deg,stop_deg_default)
    write(Unit,'("&END")')
  endif !libgenec

  return

end subroutine Write_namelist
!=======================================================================
subroutine Read_namelist
!-----------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------
! * Parse the CharacteristicsParams namelist *
  read(*,nml=CharacteristicsParams)

! * Parse the PhysicsParams namelist *
  read(*,nml=PhysicsParams)

! * Parse the CompositionParams namelist *
  read(*,nml=CompositionParams)

! * Parse the RotationParams namelist *
  read(*,nml=RotationParams)

! * Parse the WindsParams namelist *
  read(*,nml=WindsParams)

! * Parse the SurfaceParams namelist *
  read(*,nml=SurfaceParams)
  if (irot > 0) then
    fitmi_default = 0.9990d0
  else
    fitmi_default = 0.980d0
  endif
  if (fitmi == 0.0d0) then
    fitmi = fitmi_default
  endif

! * Parse the ConvectionParams namelist *
  read(*,nml=ConvectionParams)

! * Parse the ConvergenceParams namelist *
  read(*,nml=ConvergenceParams)

! * Parse the TimeControle namelist *
  read(*,nml=TimeControle)

! * Parse the VariousSettings namelist *
  read(*,nml=VariousSettings)

  return

end subroutine Read_namelist
!=======================================================================
subroutine FITM_Change(teffvv,fitmIon,m,zensi,q,notFullyIonised,BaseZC)
!-----------------------------------------------------------------------
! Changement de FITM
!-----------------------------------------------------------------------
  use caramodele,only: teff,teffv

  implicit none

  integer,intent(in):: m
  real(kindreal),intent(in):: teffvv,BaseZC,fitmIon
  real(kindreal),dimension(ldi),intent(in):: q,zensi
  logical,intent(in):: notFullyIonised

  integer:: ijk,signf
  real(kindreal):: xteffprev,ffactor,xttfitm,fitmold,fitmf,FITMfactor
  logical:: ChangeTeff=.false.
!-----------------------------------------------------------------------
  xtt=log10(teff)
  xteffprev=log10(teffv)

  if (ifitm >= 1 .and. fitm <= 0.99990d0) then
    ffactor=0.d0
    xttfitm=3.8d0
    if ((((teff-teffv)<0.d0) .and. ((teffv-teffvv)<0.d0)) .or. (((teff-teffv)>=0.d0) .and. ((teffv-teffvv)>=0.d0))) then
      ChangeTeff = .true.
    else
      ChangeTeff = .false.
    endif
    if (xttfitm < 4.d0 .and. xtt < 4.2d0) then
      if (verbose) then
        write(*,'(a,f9.6,f7.4)')'fitm,xtt= ',fitm,xtt
      endif
      fitmold=fitm
      if (irot == 1) then
        fitmf=0.98d0
      else
        fitmf=0.97d0
      endif
      select case(ifitm)
! decrease of fitm following convective zone
        case(1)
          if (xtt < 4.d0) then
            write(io_logs,*)'fully ionised up to: ',fitmIon
            if (verbose) then
              write(*,*)'fully ionised up to: ',fitmIon
            endif
            do ijk=1,m
             if (zensi(ijk) < 0.d0) then
               fitm=(fitmi+(1.d0-exp(q(ijk)))*3.d0)/4.d0
               if (xtt <= xteffprev) then
                 if (abs(fitmold-fitm) > (fitmi-fitmold+1.d-4))then
                   signf=int((fitm-fitmold)/abs(fitm-fitmold))
                   fitm=fitmold+signf*(fitmi-fitmold+1.d-4)
                 endif
               else
                 if (abs(fitmold-fitm) > (fitmold-fitmf+1.d-4))then
                   signf=int((fitm-fitmold)/abs(fitm-fitmold))
                   fitm=fitmold+signf*(fitmold-fitmf+1.d-4)
                 endif
               endif
               if (verbose) then
                 print*,'xtt=',xtt,'xteffprev=',xteffprev
               endif
               exit
             endif
            enddo
            if (fitm  >  fitmi .and. fitm  <  fitmf) then
              write(*,'(3(a,f9.6),a)')'new fitm=',fitm,'(',fitmi,',',fitmf,')'
            endif
          endif
! polynomial decrease of fitm
        case(2)
          if (irot == 1) then
            write(io_input_changes,*)'Not a good choice of ifitm'
            stop 'this is not a good choice of ifitm'
          endif
          if (notFullyIonised) then
            ffactor=(fitmf-fitmi)/(xttfitm-4.d0)**2.d0
            if (xtt < 4.2d0 .and. xtt > xttfitm-0.2d0) then
              fitm=fitmi+(xtt-4.d0)**2.d0*ffactor
              if (verbose) then
                print*,'fitm= ',fitm,' xtt= ',xtt
              endif
            endif
            if (fitm /= fitmf .and. xtt < xttfitm) then
              fitm=fitmf
            endif
            if (fitm /= fitmi .and. xtt > 4.d0) then
              fitm=fitmi
            endif
          endif
! hand-like method
        case (3,4,5,6)
          if (ifitm == 3 .or. ifitm == 5) then
            FITMfactor = 1.d0
          else if (ifitm == 6) then
            FITMfactor = 0.2d0
          else
            FITMfactor = 10.d0
          endif
          if (verbose) then
            write(*,*) 'BASE ZC :', BaseZC
            write(*,*) 'FITMION : ', fitmIon
          endif
          if (xtt < 4.d0) then
            if ((irot==1 .and. ChangeTeff) .and. ((ifitm==3 .or. ifitm==5 .or. ifitm==6) &
                                             .or. (mod(nwmd,10)==0))) then
              if (xtt<xteffprev .and. fitm>fitmf .and. BaseZC>0.d0 .and. BaseZC<0.9995d0*fitm .and. ifitm /= 5) then
                if (fitm - 0.9988d0 > 1.d-9) then
                  fitm = fitm - FITMfactor*0.00001d0
                else if (fitm - 0.9986 > 1.d-9) then
                  fitm = fitm - FITMfactor*0.00002d0
                else if (fitm - 0.9980d0 > 1.d-9) then
                  fitm = fitm - FITMfactor*0.00003d0
                else if (fitm - 0.9880d0 > 1.d-9) then
                  fitm = fitm - FITMfactor*0.0001d0
                else if (fitm - 0.9860d0 > 1.d-9) then
                  fitm = fitm - FITMfactor*0.0002d0
                else if (fitm - 0.9800d0 > 1.d-9) then
                  fitm = fitm - FITMfactor*0.0003d0
                else
                  fitm = fitm - FITMfactor*0.0005d0
                endif
              else if (xtt >= xteffprev .and. fitm < fitmi .and. fitm < fitmIon .and. BaseZC >= 0.98d0*fitm .and. ifitm /=6) then
                if (fitm-0.9900d0 < -1.d-9) then
                  fitm = fitm + FITMfactor*0.00004d0
                else if (fitm-0.9988d0 < -1.d-9) then
                  fitm = fitm + FITMfactor*0.00002d0
                else
                  fitm = fitm + FITMfactor*0.00001d0
                endif
              endif
            endif
          endif
          if (fitm < fitmf) then
            fitm = fitmf
          endif
          if (fitm > fitmi) then
            fitm = fitmi
          endif
        case default
          stop 'Bad value for ifitm, should be 0, 1, 2, 3, 4, 5 or 6'
      end select

      if (fitmold < fitm .and. xtt <= xteffprev) then
        print*, 'set back fitm,fitmold=',fitmold,'fitm=',fitm
        fitm=fitmold
      else if (fitmold > fitm .and. xtt >= xteffprev+1.d-4) then
        print*, 'set back fitm,fitmold=',fitmold,'fitm=',fitm
        fitm=fitmold
      endif
      if (fitm < fitmf .and. xtt > xttfitm) then
        fitm=fitmf
      endif
      if (fitm > fitmi) then
        fitm=fitmi
      endif
      if (fitm < fitmf) then
        fitm=fitmf
      endif
      if (abs(fitm-fitmold) >= 1.d-9) then
        write(io_input_changes,'(i7.7,a7,f12.9)') nwmd+1,': FITM=',fitm
        write(io_logs,'(i7.7,a7,f12.9)') nwmd+1,': FITM=',fitm
        write(*,*)'NEW FITM: ',fitm
      else
        fitm=fitmold
      endif
    endif
  endif
! [/Modif FITM]

  return

end subroutine FITM_Change
!=======================================================================
subroutine IMLOSS_Change(Xc,supraEdd,vequat,logTeff)
!-----------------------------------------------------------------------
! Mdot prescription modifications
! For the massives: supra-Edd multiplication factor, or WR-type Mdot
! For the low-mass: Reimers Mdot on the giant branch
!-----------------------------------------------------------------------
  implicit none

  real(kindreal),intent(in):: Xc,vequat,logTeff
  logical,intent(in):: supraEdd

  real(kindreal),parameter:: fmlosrsg=3.0d0
!-----------------------------------------------------------------------
! FMLOS change after MS
  if (fmlos == 0.85d0 .and. logTeff < 4.d0 .and. vequat < 50.d0) then
    if (Xc < 1.d-5) then
      fmlos=1.d0
      write(io_input_changes,'(i7.7,a8,d10.3)') nwmd+1,': FMLOS=',fmlos
      print*,'FMLOS changed to 1'
    endif
  endif

! WR
  ! if (logTeff >= 4.d0 .and. Xsurf < 0.3d0) then
  !   if (Xsurf > 1.d-7 .and. imloss /= 8) then
  !     imloss=8
  !     write(io_input_changes,'(i7.7,a9,i2)') nwmd+1,': IMLOSS=',imloss
  !     write(io_input_changes,*)'X(surf)= ',Xsurf
  !     print*,'IMLOSS changed to ',imloss,',X(surf)= ',Xsurf
  !   else if (Xsurf <= 1.d-7 .and. imloss /= 7) then
  !     imloss=7
  !     write(io_input_changes,'(i7.7,a9,i2)') nwmd+1,': IMLOSS=',imloss
  !     write(io_input_changes,*)'X(surf)= ',Xsurf
  !     print*,'IMLOSS changed to ',imloss,',X(surf)= ',Xsurf
  !   endif
  ! endif

! SupraEdd
  if (xmini >= 20.d0 .and. supraEdd .and. SupraEddMdot .and. phase /= 1 .and. fmlos < fmlosrsg) then
    fmlos = fmlosrsg
    write(io_input_changes,'(i7.7,a,f5.1)') nwmd+1,':  SUPRA-EDD, fmlos= ',fmlos
    print*,'Supra-Edd: Mdot multiplied by ',fmlos
  endif
  if (xmini >= 20.d0 .and. fmlos == fmlosrsg .and. .not.supraEdd) then
    fmlos = 1.d0
    write(io_input_changes,'(i7.7,a,f5.1)') nwmd+1,': no more SUPRA-EDD, fmlos back to ',fmlos
    print*,'No more Supra-Edd: fmlos back to ',fmlos
  endif

! Red Giants
  ! if (Xc < 1.d-7 .and. logTeff < 3.8d0 .and. Llast > Lprev .and. imloss /= 3 .and. xmini < 8.5d0) then
  !   if (xmini < 5.5d0 .and. fmlos /= 0.5d0) then
  !     imloss = 3
  !     fmlos=0.5d0
  !     write(*,*) nwmd+1,': IMLOSS= 3, FMLOS= 0.500'
  !     write(io_input_changes,*) nwmd+1,': IMLOSS= 3, FMLOS= 0.500'
  !   else if (xmini >= 5.5d0 .and. fmlos /= 0.6d0) then
  !     imloss = 3
  !     fmlos=0.6d0
  !     write(*,*) nwmd+1,': IMLOSS= 3, FMLOS= 0.600'
  !     write(io_input_changes,*) nwmd+1,': IMLOSS= 3, FMLOS= 0.600'
  !   endif
  ! endif

  return

end subroutine IMLOSS_Change
!=======================================================================
subroutine INPUTS_Change(Xc,Yc,Cc,Nec,Oc,rapom2,m,nzmodini,nzmodnew)
!-----------------------------------------------------------------------
! Change of the input parameters at the end of a series.
!-----------------------------------------------------------------------
  use const,only: um

  implicit none

  integer,intent(in):: m,nzmodini
  integer,intent(inout):: nzmodnew
  real(kindreal),intent(in):: Xc,Yc,Cc,Nec,Oc,rapom2
!-----------------------------------------------------------------------
   nzmodnew=nzmodini

  select case (phase)
    case (1)
      if (iauto >= 2) then
!       increase gkorm/faktor when X_c <0.1 (near the end of the MS)
        if (Xc < 0.1d0 .and. Yc > 0.5d0) then
          if (gkorm < 0.2d0) then
            gkorm = 0.2d0
            write(io_input_changes,'(i7.7,a8,f5.2)') nwmd+1,': GKORM=',gkorm
            write(*,*) 'GKORM changed to 0.2'
          endif
          if (faktor < 5.d0) then
            faktor = 5.d0
            write(io_input_changes,'(i7.7,a9,1pd9.2)') nwmd+1,': FAKTOR=',faktor
            write(*,*) 'FAKTOR changed to 5'
          endif
        endif
      endif
      if (Xc < 5.d-3 .and. iadvec == 1) then
!       switch off the advection at the end of the MS
        iadvec = 0
        idialo = 0
        idialu= 0
        xdial = 0.d0
        write(io_input_changes,'(i7.7,a,i2)') nwmd+1,': IADVEC,IDIALO,IDIALU,XDIAL= ',iadvec
        write(*,*) 'IADVEC, IDIALO, IDIALU, XDIAL changed to 0'
      endif
      if (Xc < 1.d-8 .and. Yc > 0.5d0) then
!       end of MS: PHASE 1 --> 2 and usual changes for He-b
        if (xmini > 2.0d0) then
          phase = 2
        else
          phase = 10   ! quicker timestep for the red-giants branch climbing
        endif
        write(io_input_changes,*) "------------------------------------------------"
        write(io_input_changes,'(i7.7,a8,i2)') nwmd+1,': phase=',phase
        write(*,*) 'PHASE 1 --> 2'
        if (gkorm < 0.3d0) then
          gkorm = 0.3d0
          write(io_input_changes,'(7x,a)') '  gkorm = 0.3'
        endif
        if (faktor < 10.d0) then
          faktor = 10.d0
          write(io_input_changes,'(7x,a)') '  faktor = 10'
        endif
      endif
    case (2)
      if (Yc < 1.d-8 .and. (Cc+Oc) > 0.5d0) then
!       end of He-b: PHASE 2 --> 3 and usual changes for C-b
        phase = 3
        if (iover /= 0) iover=0
        if (dovhp /= 0.d0) dovhp=0.d0
        if (gkorm < 0.5d0) gkorm=0.5d0
        if (agdr > 1.d-6) agdr = 1.d-6
        if (faktor < 1.d4) faktor = 1.d4
        if (alph > 0.8d0) alph = 0.8d0
        write(io_input_changes,*) "------------------------------------------------"
        write(io_input_changes,'(i7.7,a2)') nwmd+1,': PHASE= 3 IOVER= 0 DOVHP= 0.00\n    AGDRSPT=  1.00E-06 FAKTOR=1.00E+04'
        write(*,*) 'PHASE 2 --> 3, IOVER --> 0 +fakt+agd...'
      endif
    case (3)
      if (Cc < 1.d-5) then
!       end of C-b: PHASE 3 --> 4 and usual changes for Ne-b
        phase=4
        faktor = faktor*10.d0
        write(io_input_changes,*) "------------------------------------------------"
        write(io_input_changes,*) nwmd+1,': PHASE= 4, FAKTOR*10: ',faktor
        write(*,*) nwmd+1,': PHASE= 4, FAKTOR*10: ',faktor
      endif
    case (4)
      if (Nec < 1.d-2) then
!       end of Ne-b: PHASE 4 --> 5 and usual changes for O-b
        phase=5
        faktor = faktor*10.d0
        write(io_input_changes,*) "------------------------------------------------"
        write(*,*) nwmd+1,': PHASE= 5, FAKTOR*10: ',faktor
        write(io_input_changes,*) nwmd+1,': PHASE= 5, FAKTOR*10: ',faktor
      endif
    case (5)
      if (nzmodnew <= 10 .and. mod(nwmd,20) == 0) then
        nzmodnew=20
        write(io_input_changes,*) nwmd+1,': NZMOD= ',nzmodnew
        write(*,*) 'NZMOD --> ',nzmodnew
      endif
      if (idifcon == 0) then
        idifcon=1
        if (idiff /= 1) idiff=1
        write(io_input_changes,*)nwmd+1,': IDIFF= ',idiff,' IDIFCON= ',idifcon
        write(*,*) 'IDIFCON (IDIFF) 0 --> 1'
      endif
      if (Oc < 0.03d0) then
!       PHASE changed to 6 to have the nuclear statistical equilibrium, even if O-b not finished
        phase=6
        faktor=faktor*10.d0
        write(io_input_changes,*) "------------------------------------------------"
        write(io_input_changes,*)nwmd+1,': PHASE= 6, FAKTOR*10:',faktor
        write(*,*) nwmd+1,': PHASE= 6, FAKTOR*10:',faktor
        alph = alph - 0.1d0
        write(io_input_changes,'(i7.7,a7,f5.2)') nwmd+1,': ALPH=',alph
        write(*,*) nwmd+1,': ALPH=',alph
      endif
    case (6)
      if (idifcon == 0) then
        idifcon=1
        if (idiff /= 1) idiff=1
        write(io_input_changes,*)nwmd+1,': IDIFF= ',idiff,' IDIFCON= ',idifcon
        write(*,*) 'IDIFCON (IDIFF) 0 --> 1'
      endif
    case (10)
      continue
    case default
      rewind(io_runfile)
      write(io_runfile,*) nwmd,": Problem with the phase number"
      stop "Problem with the phase number ==> STOP"
  end select

! DGRP-L-Y
  if (m > 1500) then
    if (dgrp < 0.1d0*um .and. dgrl < 0.1d0*um) then
      dgrp = dgrp + 0.01d0*um
      dgrl = dgrl + 0.01d0*um
      write(io_input_changes,'(i7.7,a12,f6.3)') nwmd+1,': DGRP,DGRL=',dgrp/um
    else
      if (dgry < 0.005d0) then
        dgry = dgry + 0.001d0
        write(io_input_changes,'(i7.7,a7,f6.3)') nwmd+1,': DGRY=',dgry
      else if (dgrp < 0.2d0*um .and. dgrl < 0.2d0*um) then
        dgrp = dgrp + 0.01d0*um
        dgrl = dgrl + 0.01d0*um
        write(io_input_changes,'(i7.7,a12,f6.3)') nwmd+1,': DGRP,DGRL=',dgrp/um
      endif
    endif
  endif

! pour ialflu=1, nrband jamais plus grand que 10 dans combustion de l'helium
  if (ialflu == 1 .and. phase /= 1 .and. Yc /= 0.d0) then
    nrband=min(10,nrband)
  endif

  if (irot == 1 .and. isol == 0) then
    if (rapcrilim /= 0.d0 .and. rapom2 > rapcrilim .and. gkorm < 0.3d0) then
      gkorm = 0.3d0
      alph = 0.8d0
    endif
    if (rapom2 < 0.9d0*rapcrilim .and. gkorm /= 0.1d0 .and. phase < 2) then
      if (iauto == 0) then
        gkorm = 0.1d0
        alph = 1.d0
      endif
    endif
  endif   ! irot

  if (nwmd == 10) then
    alph=1.d0
    gkorm=1.d0
    if (.not. libgenec) then
    write(io_input_changes,*) nwmd+1,': alph,gkorm=',alph,gkorm
    endif
  else if (modanf == 1) then
    gkorm=0.3d0
    if (.not. libgenec) then
    write(io_input_changes,*) nwmd+1,': gkorm=',gkorm
    endif
  else if (modanf == 5) then
    if (irot == 0) then
      gkorm = 0.10d0
      islow = 0
      if (.not. libgenec) then
      write(io_input_changes,*) nwmd+1,': gkorm,islow=',gkorm,islow
      endif
    endif
  endif

  return

end subroutine INPUTS_Change

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
      write(*,*) '  4: WINDS inputs'
      write(*,*) '  5: SURFACE inputs'
      write(*,*) '  6: CONVECTION inputs'
      write(*,*) '  7: CONVERGENCE inputs'
      write(*,*) '  8: TIME CONTROL inputs'
      write(*,*) '  9: VARIOUS SETTINGS inputs'
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
          write(*,'(a,f7.5)') ' 6: A_M03    :',A_M03
          write(*,'(a,i2)') ' 7: n_M03    :',n_M03
          write(*,'(a,d11.5)') ' 8: B_initial:',B_initial
          write(*,'(a,d11.5)') ' 9: add_diff :',add_diff
          write(*,'(a,i2)') '10: n_mag    :',n_mag
          write(*,'(a,f7.3)') ' 11: alpha_F  :',alpha_F
          write(*,'(a,i2)') '12: nsmooth  :',nsmooth
          write(*,'(a,l2)') '13: qminsmooth:',qminsmooth
          write(*,'(a,l2)') '14: dcirch_inclusion:',dcirch_inclusion
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
              write(*,*)'Enter the desired value for A_M03:'
              write(*,*)'     (old value of A_Dh=0.002 if A_M03=0.d0 and n_M03=0)'
              read(5,*) Temp_Var_real
            enddo
            A_M03 = Temp_Var_real
          case (7)
            Temp_Var_Int = -99
            do while (Temp_Var_Int < 0)
              write(*,*)'Enter the desired value for n_M03:'
              write(*,*)'     (old value of A_Dh=0.002 if A_M03=0.d0 and n_M03=0)'
              read(5,*) Temp_Var_Int
            enddo
            n_M03 = Temp_Var_Int
          case (8)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*)'Enter the desired value for B_initial (in Gauss):'
              read(5,*) Temp_Var_real
            enddo
            B_initial = Temp_Var_real
          case(9)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*)'Enter the desired value for add_diff:'
              read(5,*) Temp_Var_real
            enddo
            add_diff = Temp_Var_real
          case(10)
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
          case(11)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*)'Enter the desired value for alpha_F (default 1.0):'
              read(5,*) Temp_Var_real
            enddo
            alpha_F = Temp_Var_real
          case(12)
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
          case (13)
            Temp_Var_char = ''
            do while (Temp_Var_char/='t' .and. Temp_Var_char/='f' &
                 .and. Temp_Var_char/='T' .and. Temp_Var_char/= 'F')
              write(*,*)'Enter the desired value for qminsmooth (T/F, default T):'
              read(5,*) Temp_Var_char
            enddo
            if (Temp_Var_char=='t' .or. Temp_Var_char=='f') then
              qminsmooth = .true.
            elseif (Temp_Var_char=='0' .or. Temp_Var_char=='F') then
              qminsmooth = .false.
            endif
          case (14)
            Temp_Var_char = ''
            do while (Temp_Var_char/='t' .and. Temp_Var_char/='f' &
                 .and. Temp_Var_char/='T' .and. Temp_Var_char/= 'F')
              write(*,*)'Enter the desired value for dcirch_inclusion (T/F, default F):'
              read(5,*) Temp_Var_char
            enddo
            if (Temp_Var_char=='t' .or. Temp_Var_char=='f') then
              dcirch_inclusion = .true.
            elseif (Temp_Var_char=='0' .or. Temp_Var_char=='F') then
              dcirch_inclusion = .false.
            endif
          case default
            write(*,*) 'Wrong number, should be an integer between 0 and 14'
          end select ! end ROTATION inputs selection
        enddo
      case (4) ! *** change of WINDS inputs
        Change_params = 99
        do while (Change_params /= 0)
          write(*,*) '------------------------------'
          write(*,*) '*** WINDS inputs ***'
          write(*,'(a,i3)') ' 1: OB_Mdot       :',OB_Mdot
          write(*,'(a,i3)') ' 2: RSG_Mdot      :',RSG_Mdot
          write(*,'(a,i3)') ' 3: WR_Mdot       :',WR_Mdot
          write(*,'(a,i3)') ' 4: Fallback_Mdot :',Fallback_Mdot
          write(*,'(a,d12.5)') ' 5: fmlos         :',fmlos
          write(*,'(a,f6.2)') ' 6: Z_dep         :',Z_dep
          write(*,'(a,f6.2)') ' 7: Xs_WR         :',Xs_WR
          write(*,'(a,f6.2)') ' 8: D_clump         :',D_clump
          write(*,'(a,l2)') ' 9: SupraEddMdot  :',SupraEddMdot
          write(*,'(a,f6.2)') '10: Be_Mdotfrac   :',Be_mdotfrac
          write(*,'(a,f6.2)') '11: start_mdot    :',start_mdot
          write(*,'(a,l2)') '12: hardJump      :',hardJump
          write(*,'(a,l2)') '13: force_prescription:',force_prescription
          write(*,'(a,l2)') '14: print_winds   :',print_winds
          write(*,'(a,l2)') '15: prezams_winds_not_applied:',prezams_winds_not_applied
          write(*,*) '------------------------------'
          write(*,*) 'Parameters to change (0 to skip or exit):'
          read(5,*) Change_params
          select case (Change_params)
          case (0)
            write(*,*) 'No more changes of WINDS parameters'
          case (1)
            Temp_Var_Int = 99
            do while (Temp_Var_Int>=14)
              write(*,*) 'Possible values for OB_MDOT'
              write(*,*) '------------------------------'
              write(*,*) '  0: none'
              write(*,*) '  1: de Jager+ 1988'
              write(*,*) '  2: mass loss in Msol/yr given by FMLOS'
              write(*,*) '  3: de Jager+ 1988 (linear)'
              write(*,*) '  4: Vink+ 2001'
              write(*,*) '  5: Vink+ 2001 modified by Markova & Puls 2008'
              write(*,*) '  6: Kudritzki & Puls 2000'
              write(*,*) '  7: Kudritzki 2002'
              write(*,*) '  8: Bestenlehner+ 2020'
              write(*,*) '  9: Bjorklund+ 2023'
              write(*,*) ' 10: Gormaz-Matamala+ 2022'
              write(*,*) ' 11: Krticka+ 2021'
              write(*,*) ' 12: Sabhahit+ 2022'
              write(*,*) ' 13: Grafener 2021'
              write(*,*) '------------------------------'
              write(*,*) 'Enter the desired value:'
              read(5,*) Temp_Var_Int
            enddo
            OB_Mdot = Temp_Var_Int
          case (2)
            Temp_Var_Int = 99
            do while (Temp_Var_Int>=18)
              write(*,*) 'Possible values for RSG_MDOT'
              write(*,*) '------------------------------'
              write(*,*) '  0: none'
              write(*,*) '  1: de Jager+ 1988'
              write(*,*) '  2: mass loss in Msol/yr given by FMLOS'
              write(*,*) '  3: de Jager+ 1988 (linear)'
              write(*,*) '  4: Crowther 2001 (standard GENEC)'
              write(*,*) '  5: Beasor & Davies 2020'
              write(*,*) '  6: Kee+ 2021'
              write(*,*) '  7: Reimers 1975'
              write(*,*) '  8: van Loon+ 2005'
              write(*,*) '  9: Nieuwanhuijzen 1990'
              write(*,*) ' 10: Vanbeveren 1998'
              write(*,*) ' 11: Salasnich 1999'
              write(*,*) ' 12: Decin 2021'
              write(*,*) ' 13: Decin 2024'
              write(*,*) ' 14: Yang 2023'
              write(*,*) ' 15: Wachter 2002'
              write(*,*) ' 16: Schroder 2005'
              write(*,*) ' 17: Vink+ 2023'
              write(*,*) '------------------------------'
              write(*,*) 'Enter the desired value (default 4):'
              read(5,*) Temp_Var_Int
            enddo
            RSG_Mdot = Temp_Var_Int
          case (3)
            Temp_Var_Int = 99
            do while (Temp_Var_Int>=15)
              write(*,*) 'Possible values for WR_MDOT'
              write(*,*) '------------------------------'
              write(*,*) '  0: none'
              write(*,*) '  1: de Jager+ 1988'
              write(*,*) '  2: mass loss in Msol/yr given by FMLOS'
              write(*,*) '  3: de Jager+ 1988 (linear)'
              write(*,*) '  4: Graefener & Hammann 2008'
              write(*,*) '  5: Nugis & Lamers 2000'
              write(*,*) '  6: Schmutz 1997, except for WNL = Nugis+ 1998'
              write(*,*) '  7: Hainich 2015'
              write(*,*) '  8: Langer 1989'
              write(*,*) '  9: Yoon+ 2006'
              write(*,*) ' 10: Nugis & Lamers 2000, combined eq. for WN and WC'
              write(*,*) ' 11: Sander 2020'
              write(*,*) ' 12: Vink 2017'
              write(*,*) ' 13: Shenar 2019'
              write(*,*) ' 14: Tramper 2016'
              write(*,*) '------------------------------'
              write(*,*) 'Enter the desired value (default 4):'
              read(5,*) Temp_Var_Int
            enddo
            WR_Mdot = Temp_Var_Int
          case (4)
            Temp_Var_Int = 99
            do while (Temp_Var_Int>=4)
              write(*,*) 'Possible values for FALLBACK_MDOT'
              write(*,*) '------------------------------'
              write(*,*) ' 0: none'
              write(*,*) ' 1: de Jager+ 1988'
              write(*,*) ' 2: mass loss in Msol/yr given by FMLOS'
              write(*,*) ' 3: de Jager+ 1988 (linear)'
              write(*,*) '------------------------------'
              write(*,*)'Enter the desired value (recommended 1):'
              read(5,*) Temp_Var_Int
            enddo
            Fallback_Mdot = Temp_Var_Int
          case (5)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*)'Enter the desired value for fmlos:'
              read(5,*) Temp_Var_real
            enddo
            fmlos = Temp_Var_real
          case (6)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*)'Enter the desired value for Z_dep (recommended 0.85):'
              read(5,*) Temp_Var_real
            enddo
            Z_dep = Temp_Var_real
          case (7)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*) 'Enter the desired value for Xs_WR (recommended 0.3):'
              read(5,*) Temp_Var_real
            enddo
            Xs_WR = Temp_Var_real
          case (8)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*) 'Enter the desired value for D_clump (recommended 10.d0):'
              read(5,*) Temp_Var_real
            enddo
            D_clump = Temp_Var_real
          case (9)
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
          case (10)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*) 'Enter the desired value for Be_Mdotfrac (recommended 0.1):'
              read(5,*) Temp_Var_real
            enddo
            Be_mdotfrac = Temp_Var_real
          case (11)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*) 'Enter the desired value for start_mdot (recommended 0.8):'
              read(5,*) Temp_Var_real
            enddo
            start_mdot = Temp_Var_real
          case (12)
            Temp_Var_char = ''
            do while (Temp_Var_char/='t' .and. Temp_Var_char/='f' &
                 .and. Temp_Var_char/='T' .and. Temp_Var_char/= 'F')
              write(*,*)'Enter the desired value for hardJump (T=default/F):'
              read(5,*) Temp_Var_char
            enddo
            if (Temp_Var_char=='t' .or. Temp_Var_char=='T') then
              hardJump = .true.
            elseif (Temp_Var_char=='f' .or. Temp_Var_char=='F') then
              hardJump = .false.
            endif
          case (13)
            Temp_Var_char = ''
            do while (Temp_Var_char/='t' .and. Temp_Var_char/='f' &
                 .and. Temp_Var_char/='T' .and. Temp_Var_char/= 'F')
              write(*,*)'Enter the desired value for force_prescription (T/F=default):'
              read(5,*) Temp_Var_char
            enddo
            if (Temp_Var_char=='t' .or. Temp_Var_char=='T') then
              force_prescription = .true.
            elseif (Temp_Var_char=='f' .or. Temp_Var_char=='F') then
              force_prescription = .false.
            endif
          case (14)
            Temp_Var_char = ''
            do while (Temp_Var_char/='t' .and. Temp_Var_char/='f' &
                 .and. Temp_Var_char/='T' .and. Temp_Var_char/= 'F')
              write(*,*)'Enter the desired value for print_winds (T/F=default):'
              read(5,*) Temp_Var_char
            enddo
            if (Temp_Var_char=='t' .or. Temp_Var_char=='T') then
              print_winds = .true.
            elseif (Temp_Var_char=='f' .or. Temp_Var_char=='F') then
              print_winds = .false.
            endif
          case (15)
            Temp_Var_char = ''
            do while (Temp_Var_char/='t' .and. Temp_Var_char/='f' &
                 .and. Temp_Var_char/='T' .and. Temp_Var_char/= 'F')
              write(*,*)'Enter the desired value for prezams_winds_not_applied (T/F=default):'
              read(5,*) Temp_Var_char
            enddo
            if (Temp_Var_char=='t' .or. Temp_Var_char=='T') then
              prezams_winds_not_applied = .true.
            elseif (Temp_Var_char=='f' .or. Temp_Var_char=='F') then
              prezams_winds_not_applied = .false.
            endif
          case default
            write(*,*) 'Wrong number, should be an integer between 0 and 15'
          end select ! end WINDS inputs selection
        enddo
      case(5) ! *** change of SURFACE inputs
        Change_params = 99
        do while (Change_params /= 0)
          write(*,*) '------------------------------'
          write(*,*) '*** SURFACE inputs ***'
          write(*,'(a,f9.5)') ' 1: fitm          :',fitm
          write(*,'(a,i2)') ' 2: ifitm         :',ifitm
          write(*,*) '------------------------------'
          write(*,*) 'Parameters to change (0 to skip or exit):'
          read(5,*) Change_params
          select case (Change_params)
          case (0)
            write(*,*) 'No more changes of SURFACE parameters'
          case (1)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*) 'Enter the desired value for fitm:'
              write(*,*) '     default: 0.980 for irot=0'
              write(*,*)'               0.999 for irot=1'
              read(5,*) Temp_Var_real
            enddo
            fitm = Temp_Var_real
            fitmi = fitm
          case (2)
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
            write(*,*) 'Wrong number, should be an integer between 0 and 2'
          end select ! end SURFACE inputs selection
        enddo
      case (6) ! *** change of CONVECTION inputs
        Change_params = 99
        do while (Change_params /= 0)
          write(*,*) '------------------------------'
          write(*,*) '*** CONVECTION inputs ***'
          write(*,'(a,i2)') ' 1: iledou   :',iledou
          write(*,'(a,i2)') ' 2: iover    :',iover
          write(*,'(a,f7.3)') ' 3: dovhp    :',dovhp
          write(*,'(a,i2)') ' 4: iunder   :',iunder
          write(*,'(a,f7.3)') ' 5: dunder   :',dunder
          write(*,'(a,f7.3)') ' 6: elph     :',elph
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
              write(*,*)'Enter the desired value for overshooting parameter dovhp:'
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
              write(*,*)'Enter the desired value for undershooting parameter dunder:'
              read(5,*) Temp_Var_real
            enddo
            dunder = Temp_Var_real
          case (6)
            Temp_Var_real = -2.d0
            do while (Temp_Var_real < 0.d0)
              write(*,*)'Enter the desired value for MLT parameter elph (recommended 1.6):'
              read(5,*) Temp_Var_real
            enddo
            elph = Temp_Var_real
          case default
            write(*,*) 'Wrong number, should be an integer between 0 and 6'
          end select ! end CONVECTION inputs selection
        enddo
      case (7) ! *** change of CONVERGENCE inputs
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
      case (8) ! *** change of TIME CONTROLE inputs
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
      case (9) ! *** change of VARIOUS SETTINGS inputs
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
        write(*,*) 'Wrong number, should be an integer between 0 and 9'
      end select ! category selection
    enddo
  endif

end subroutine Ask_changes
!=======================================================================
end module inputparam
