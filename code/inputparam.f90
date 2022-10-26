module inputparam

  use io_definitions
  use evol,only: ldi,kindreal
  use caramodele,only: nwmd,xmini

  implicit none

  interface Write_param
    module procedure Write_param_int
    module procedure Write_param_real
    module procedure Write_param_logical
  end interface Write_param

  integer,parameter:: imagn_default=0,ianiso_default=0,ipop3_default=0,ibasnet_default=0,iopac_default=3,&
    ikappa_default=5,istati_default=0,igamma_default=0,nndr_default=1,iledou_default=0,idifcon_default=0,&
    iunder_default=0,nbchx_default=200,nrband_default=1,icncst_default=0,iprn_default=10,&
    iout_default=0,itmin_default=5,idebug_default=0,itests_default=0,tauH_fit_default=1,RSG_Mdot_default=0,&
    n_mag_default=1,nsmooth_default=1,end_at_phase_default=4,end_at_model_default=0,iprezams_default=0,&
    n_snap_default=10
  real(kindreal),parameter:: fenerg_default=1.0d0,richac_default=1.0d0,zsol_default=1.40d-2,frein_default=0.0d0,&
    K_Kawaler_default=0.d0,Omega_saturation_default=14.d0,vwant_default=0.0d0,xfom_default=1.0d0, &
    dunder_default=0.0d0,dgro_default=0.010d0,dgr20_default=0.010d0,binm2_default=0.d0,periodini_default=0.d0,&
    B_initial_default=0.d0,add_diff_default=0.0d0,Be_mdotfrac_default=0.0d0,start_mdot_default=0.80d0,&
    alpha_F_default=1.d0
  logical,parameter:: xyfiles_default=.false.,bintide_default=.false.,const_per_default=.true.,&
    var_rates_default=.false.,verbose_default=.false.,Add_Flux_default = .true.,&
    diff_only_default=.false.,stop_deg_default=.true.,noSupraEddMdot_default=.false.,&
    qminsmooth_default=.false.

! VARIABLES DE LECTURE
  integer,save:: idern,ichem,itminc

! NAMELISTS VARIABLES
! **** Model characteristics
  integer,save:: nwseq,modanf,nzmod,end_at_phase=end_at_phase_default, &
                 end_at_model=end_at_model_default
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
  namelist /PhysicsParams/irot,isol,imagn,ialflu,ianiso,ipop3,ibasnet,phase,var_rates,bintide,binm2,&
            periodini,const_per,iprezams
!-----------------------------------------------------------------------

! **** Chemical composition
  integer,save:: iopac=iopac_default,ikappa=ikappa_default
  real(kindreal),save:: zinit,zsol=zsol_default,z
!-----------------------------------------------------------------------
  namelist /CompositionParams/zinit,zsol,z,iopac,ikappa
!-----------------------------------------------------------------------

! **** Rotation-linked parameters
  integer,save:: idiff,iadvec,istati=istati_default,icoeff,igamma=igamma_default,idialo,idialu,n_mag=n_mag_default, &
       nsmooth=nsmooth_default
  real(kindreal),save:: fenerg=fenerg_default,richac=richac_default,frein=frein_default,K_Kawaler=K_Kawaler_default, &
                          Omega_saturation=Omega_saturation_default,rapcrilim,vwant=vwant_default,&
                          xfom=xfom_default,omega,xdial,B_initial=B_initial_default,add_diff=add_diff_default &
                          ,alpha_F=alpha_F_default
  logical,save:: Add_Flux=Add_Flux_default,diff_only=diff_only_default,qminsmooth=qminsmooth_default
!-----------------------------------------------------------------------
  namelist /RotationParams/idiff,iadvec,istati,icoeff,fenerg,richac,igamma,frein,K_Kawaler,Omega_saturation,rapcrilim, &
       vwant,xfom,omega,xdial,idialo,idialu,Add_Flux,diff_only,B_initial,add_diff,alpha_F,n_mag, &
       nsmooth,qminsmooth
!-----------------------------------------------------------------------

! **** Surface parameters
  integer,save:: imloss,ifitm,nndr=nndr_default,RSG_Mdot=RSG_Mdot_default
  real(kindreal),save:: fmlos,fitm,fitmi,fitmi_default,deltal,deltat,Be_mdotfrac=Be_mdotfrac_default,start_mdot=start_mdot_default
  logical,save:: noSupraEddMdot=noSupraEddMdot_default
!-----------------------------------------------------------------------
  namelist /SurfaceParams/imloss,fmlos,ifitm,fitm,fitmi,deltal,deltat,nndr,RSG_Mdot,noSupraEddMdot,Be_mdotfrac,start_mdot
!-----------------------------------------------------------------------

! **** Convection-linked parameters
  integer,save:: iledou=iledou_default,idifcon=idifcon_default,my,iover,iunder=iunder_default
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
  logical,save:: display_plot,xyfiles=xyfiles_default,verbose=verbose_default,stop_deg=stop_deg_default
!-----------------------------------------------------------------------
  namelist /VariousSettings/display_plot,iauto,iprn,iout,itmin,xyfiles,idebug,itests,verbose,stop_deg,n_snap
!-----------------------------------------------------------------------

  integer:: isugi=1
  real(kindreal),save:: xtt,agds,agdp,agdt

  public:: Write_param
  public
  private :: xtt
  private:: Write_param_int,Write_param_real,Write_param_logical
  private:: imagn_default,ianiso_default,ipop3_default,ibasnet_default,iopac_default,ikappa_default,istati_default,&
    igamma_default,nndr_default,iledou_default,iunder_default,nbchx_default,nrband_default, &
    icncst_default,iprn_default,iout_default,itmin_default,fenerg_default,richac_default,zsol_default, &
    frein_default,K_Kawaler_default,Omega_saturation_default,vwant_default,xfom_default,dunder_default,dgr20_default, &
    xyfiles_default,idebug_default,bintide_default,binm2_default,periodini_default,const_per_default,tauH_fit_default,&
    var_rates_default,verbose_default,stop_deg_default,n_mag_default,alpha_F_default,nsmooth_default,&
    RSG_Mdot_default,noSupraEddMdot_default,Be_mdotfrac_default,start_mdot_default,iprezams_default,n_snap_default

contains
!=======================================================================
subroutine Write_param_int(Unit,n_name,n_in,n_default)
!-----------------------------------------------------------------------
  implicit none

  integer,intent(in):: Unit,n_in,n_default
  character(*),intent(in):: n_name
!-----------------------------------------------------------------------
  if (n_in /= n_default) then
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
  if (x_in /= x_default) then
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
  if (a_in .neqv. a_default) then
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
  write(Unit,'(1x,a,i0)') "ialflu=",ialflu
  call Write_param(Unit,"ianiso=",ianiso,ianiso_default)
  call Write_param(Unit,"ipop3=",ipop3,ipop3_default)
  call Write_param(Unit,"ibasnet=",ibasnet,ibasnet_default)
  write(Unit,'(1x,a,i0)') "phase=",phase
  call Write_param(Unit,"iprezams=",iprezams,iprezams_default)
  call Write_param(Unit,"var_rates=",var_rates,var_rates_default)
  call Write_param(Unit,"bintide=",bintide,bintide_default)
  if (bintide) then
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
  write(Unit,'("&END"/)')

  if (irot > 0) then
    fitmi_default = 0.9990d0
  else
    fitmi_default = 0.980d0
  endif
  if (fitmi == 0.0d0) then
    fitmi = fitmi_default
  endif
  write(Unit,'(a)') "&SurfaceParams"
  write(Unit,'(1x,a,i0,a,d10.3)') "imloss=",imloss,", fmlos=",fmlos
  call Write_param(Unit,"RSG_Mdot=",RSG_Mdot,RSG_Mdot_default)
  call Write_param(Unit,"noSupraEddMdot=",noSupraEddMdot,noSupraEddMdot_default)
  call Write_param(Unit,"Be_mdotfrac=",Be_mdotfrac,Be_mdotfrac_default)
  call Write_param(Unit,"start_mdot=",start_mdot,start_mdot_default)
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
  write(Unit,'(1x,2(a,f7.4))') "dgrp=",dgrp/um,", dgrl=",dgrl/um
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
  call Write_param(Unit,"xyfiles=",xyfiles,xyfiles_default)
  call Write_param(Unit,"idebug=",idebug,idebug_default)
  call Write_param(Unit,"itests=",itests,itests_default)
  call Write_param(Unit,"verbose=",verbose,verbose_default)
  call Write_param(Unit,"stop_deg=",stop_deg,stop_deg_default)
  write(Unit,'("&END")')

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

! * Parse the SurfaceParams namelist *
  read(*,nml=SurfaceParams)

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
            if ((irot==1 .and. ChangeTeff) .and. ((ifitm==3 .or. ifitm==5 .or. ifitm==6) .or. (nwmd-nwseq)==9)) then
              if (xtt<xteffprev .and. fitm>fitmf .and. BaseZC>0.d0 .and. BaseZC<0.9995d0*fitm .and. ifitm /= 5) then
                if (fitm - 0.9998d0 > 1.d-9) then
                  fitm = fitm - FITMfactor*0.00001d0
                else if (fitm - 0.9996 > 1.d-9) then
                  fitm = fitm - FITMfactor*0.00002d0
                else if (fitm - 0.9990d0 > 1.d-9) then
                  fitm = fitm - FITMfactor*0.00003d0
                else if (fitm - 0.9980d0 > 1.d-9) then
                  fitm = fitm - FITMfactor*0.0001d0
                else if (fitm - 0.9960d0 > 1.d-9) then
                  fitm = fitm - FITMfactor*0.0002d0
                else if (fitm - 0.9900d0 > 1.d-9) then
                  fitm = fitm - FITMfactor*0.0003d0
                else
                  fitm = fitm - FITMfactor*0.0005d0
                endif
              else if (xtt >= xteffprev .and. fitm < fitmi .and. fitm < fitmIon .and. BaseZC >= 0.98d0*fitm .and. ifitm /=6) then
                if (fitm-0.9900d0 < -1.d-9) then
                  fitm = fitm + FITMfactor*0.00004d0
                else if (fitm-0.9998d0 < -1.d-9) then
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
subroutine IMLOSS_Change(Xc,Xsurf,Lprev,Llast,supraEdd,vequat,logTeff)
!-----------------------------------------------------------------------
! Mdot prescription modifications
! For the massives: supra-Edd multiplication factor, or WR-type Mdot
! For the low-mass: Reimers Mdot on the giant branch
!-----------------------------------------------------------------------
  implicit none

  real(kindreal),intent(in):: Xc,Xsurf,Lprev,Llast,vequat,logTeff
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
  if (logTeff >= 4.d0 .and. Xsurf < 0.3d0) then
    if (Xsurf > 1.d-7 .and. imloss /= 8) then
      imloss=8
      write(io_input_changes,'(i7.7,a9,i2)') nwmd+1,': IMLOSS=',imloss
      write(io_input_changes,*)'X(surf)= ',Xsurf
      print*,'IMLOSS changed to ',imloss,',X(surf)= ',Xsurf
    else if (Xsurf <= 1.d-7 .and. imloss /= 7) then
      imloss=7
      write(io_input_changes,'(i7.7,a9,i2)') nwmd+1,': IMLOSS=',imloss
      write(io_input_changes,*)'X(surf)= ',Xsurf
      print*,'IMLOSS changed to ',imloss,',X(surf)= ',Xsurf
    endif
  endif

! SupraEdd
  if (xmini >= 20.d0 .and. supraEdd .and. .not.noSupraEddMdot .and. phase /= 1 .and. fmlos < fmlosrsg) then
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
  if (Xc < 1.d-7 .and. logTeff < 3.8d0 .and. Llast > Lprev .and. imloss /= 3 .and. xmini < 8.5d0) then
    if (xmini < 5.5d0 .and. fmlos /= 0.5d0) then
      imloss = 3
      fmlos=0.5d0
      write(*,*) nwmd+1,': IMLOSS= 3, FMLOS= 0.500'
      write(io_input_changes,*) nwmd+1,': IMLOSS= 3, FMLOS= 0.500'
    else if (xmini >= 5.5d0 .and. fmlos /= 0.6d0) then
      imloss = 3
      fmlos=0.6d0
      write(*,*) nwmd+1,': IMLOSS= 3, FMLOS= 0.600'
      write(io_input_changes,*) nwmd+1,': IMLOSS= 3, FMLOS= 0.600'
    endif
  endif

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

  if (modanf == 0) then
    alph=1.d0
    gkorm=1.d0
  else if (modanf == 1) then
    gkorm=0.3d0
  else if (modanf == 5) then
    if (irot == 0) then
      gkorm = 0.10d0
      islow = 0
    endif
  endif

  return

end subroutine INPUTS_Change
!=======================================================================
end module inputparam
