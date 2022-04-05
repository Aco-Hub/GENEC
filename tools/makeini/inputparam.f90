module inputparam

  implicit none

! Definition de kindreal, pour eviter de devoir ajouter le module evol:
  integer, parameter:: kindreal = 8
  integer,parameter:: imagn_default=0,ianiso_default=0,ipop3_default=0,ibasnet_default=0,iopac_default=3,&
    ikappa_default=5,istati_default=0,igamma_default=0,nndr_default=1,iledou_default=0,idifcon_default=0,&
    iover_default=1,iunder_default=0,nbchx_default=200,nrband_default=1,icncst_default=0,iprn_default=99,&
    iout_default=0,itmin_default=5,idebug_default=0,itests_default=0,tauH_fit_default=1,RSG_Mdot_default=0
  real(kindreal),parameter:: fenerg_default=1.0d0,richac_default=1.0d0,zsol_default=1.40d-2,frein_default=0.0d0,&
    K_Kawaler_default=0.d0,Omega_saturation_default=14.d0,vwant_default=0.0d0,xfom_default=1.0d0, &
    dunder_default=0.0d0,dgro_default=0.010d0,dgr20_default=0.010d0,binm2_default=0.d0,periodini_default=0.d0,&
    B_initial_default=0.d0,add_diff_default=0.0d0,Be_mdotfrac_default=0.0d0,start_mdot_default=0.80d0
  logical,parameter:: xyfiles_default=.false.,bintide_default=.false.,const_per_default=.true.,&
    var_rates_default=.false.,verbose_default=.false.,Add_Flux_default = .true.,&
    diff_only_default=.false.,stop_deg_default=.true.,noSupraEddMdot_default=.false.
  logical,parameter:: amuseinterface_default=.false.

! NAMELISTS VARIABLES
! **** Model characteristics
  integer,save:: nwseq,modanf,nzmod
  logical,save:: amuseinterface=amuseinterface_default  
  character(256),save:: starname
!-----------------------------------------------------------------------
  namelist /CharacteristicsParams/starname,nwseq,modanf,nzmod
!-----------------------------------------------------------------------

! **** Physical inputs
  integer,save:: irot,isol,imagn=imagn_default,ialflu,ianiso=ianiso_default,ipop3=ipop3_default,&
      ibasnet=ibasnet_default,phase
  real(kindreal),save:: binm2=binm2_default,periodini=periodini_default
  logical,save:: var_rates=var_rates_default,bintide=bintide_default,const_per=const_per_default
!-----------------------------------------------------------------------
  namelist /PhysicsParams/irot,isol,imagn,ialflu,ianiso,ipop3,ibasnet,phase,var_rates,bintide,binm2,periodini,const_per
!-----------------------------------------------------------------------

! **** Chemical composition
  integer,save:: iopac=iopac_default,ikappa=ikappa_default
  real(kindreal),save:: zinit,zsol=zsol_default,z
!-----------------------------------------------------------------------
  namelist /CompositionParams/zinit,zsol,z,iopac,ikappa
!-----------------------------------------------------------------------

! **** Rotation-linked parameters
  integer,save:: idiff,iadvec,istati=istati_default,icoeff,igamma=igamma_default,idialo,idialu
  real(kindreal),save:: fenerg=fenerg_default,richac=richac_default,frein=frein_default,K_Kawaler=K_Kawaler_default, &
                          Omega_saturation=Omega_saturation_default,rapcrilim,vwant=vwant_default,&
                          xfom=xfom_default,omega,xdial,B_initial=B_initial_default,add_diff=add_diff_default
  logical,save:: Add_Flux=Add_Flux_default,diff_only=diff_only_default
!-----------------------------------------------------------------------
  namelist /RotationParams/idiff,iadvec,istati,icoeff,fenerg,richac,igamma,frein,K_Kawaler,Omega_saturation,rapcrilim, &
                           vwant,xfom,omega,xdial,idialo,idialu,Add_Flux,diff_only,B_initial,add_diff
!-----------------------------------------------------------------------

! **** Surface parameters
  integer,save:: imloss,ifitm,nndr=nndr_default,RSG_Mdot=RSG_Mdot_default
  real(kindreal),save:: fmlos,fitm,fitmi,deltal,deltat,fitmi_default,Be_mdotfrac=Be_mdotfrac_default,start_mdot=start_mdot_default
  logical,save:: noSupraEddMdot=noSupraEddMdot_default
!-----------------------------------------------------------------------
  namelist /SurfaceParams/imloss,fmlos,RSG_Mdot,noSupraEddMdot,ifitm,fitm,fitmi,deltal,deltat,nndr,Be_mdotfrac,start_mdot
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
      idebug=idebug_default,itests=itests_default
  logical,save:: display_plot,xyfiles=xyfiles_default,verbose=verbose_default,stop_deg=stop_deg_default
!-----------------------------------------------------------------------
  namelist /VariousSettings/display_plot,iauto,iprn,iout,itmin,xyfiles,idebug,itests,verbose,stop_deg
!-----------------------------------------------------------------------

  public
!  private:: Write_param_int,Write_param_real,Write_param_logical
  private:: imagn_default,ianiso_default,ipop3_default,ibasnet_default,iopac_default,ikappa_default,istati_default,&
    igamma_default,nndr_default,iledou_default,iover_default,iunder_default,nbchx_default,nrband_default, &
    icncst_default,iprn_default,iout_default,itmin_default,fenerg_default,richac_default,zsol_default, &
    frein_default,K_Kawaler_default,Omega_saturation_default,vwant_default,xfom_default,dunder_default,dgr20_default, &
    xyfiles_default,idebug_default,bintide_default,binm2_default,periodini_default,const_per_default, &
    var_rates_default,verbose_default,stop_deg_default,tauH_fit_default,noSupraEddMdot_default,RSG_Mdot_default,&
    Be_mdotfrac_default,start_mdot_default

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
  write(Unit,'("&END"/)')

  if (irot > 0) then
    fitmi_default = 0.9990d0
  else
    fitmi_default = 0.980d0
  endif
  write(Unit,'(a)') "&SurfaceParams"
  write(Unit,'(1x,a,i0,a,d10.3)') "imloss=",imloss,", fmlos=",fmlos
  write(Unit,'(1x,a,i0,1x,a,l2)') "RSG_Mdot=",RSG_Mdot,"noSupraEddMdot=",noSupraEddMdot
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
  write(Unit,'(1x,a,i0)') "iprn=",iprn
  write(Unit,'(1x,a,i0)') "iout=",iout
  write(Unit,'(1x,a,i0)') "itmin=",itmin
  write(Unit,'(1x,a,l1)') "xyfiles=",xyfiles
  write(Unit,'(1x,a,i0)') "idebug=",idebug
  write(Unit,'(1x,a,i0)') "itests=",itests
  write(Unit,'(1x,a,l1)') "verbose=",verbose
  write(Unit,'(1x,a,l1)') "stop_deg=",stop_deg
  write(Unit,'("&END")')

  return

end subroutine Write_namelist
!=======================================================================
end module inputparam
