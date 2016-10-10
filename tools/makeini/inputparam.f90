module inputparam

  implicit none

  interface Write_param
    module procedure Write_param_int
    module procedure Write_param_real
    module procedure Write_param_logical
  end interface Write_param

! Definition de kindreal, pour eviter de devoir ajouter le module evol:
  integer, parameter:: kindreal = 8

! **** Model characteristics
  namelist /CharacteristicsParams/starname,nwseq,modanf,nzmod
    character(256),save:: starname
    integer,save:: nwseq,modanf,nzmod

! **** Physical inputs
  namelist /PhysicsParams/irot,isol,imagn,ialflu,ianiso,ipop3,ibasnet,phase
    integer,save:: irot,isol,imagn=0,ialflu,ianiso,ipop3=0,ibasnet=0,phase

! **** Chemical composition
  namelist /CompositionParams/zinit,zsol,z,iopac,ikappa
    integer,save:: iopac=3,ikappa=5
    real(kindreal),save:: zinit,zsol=1.40d-2,z

! **** Rotation-linked parameters
  namelist /RotationParams/idiff,iadvec,istati,icoeff,fenerg,richac,igamma,frein,rapcrilim,vwant,xfom,omega,xdial,idialo,idialu
    integer,save:: idiff,iadvec,istati=0,icoeff,igamma=0,idialo,idialu
    real(kindreal),save:: fenerg=1.0d0,richac=1.0d0,frein=0.0d0,rapcrilim,vwant=0.0d0,xfom=1.0d0,omega,xdial

! **** Surface parameters
  namelist /SurfaceParams/imloss,fmlos,ifitm,fitm,deltal,deltat,nndr
    integer,save:: imloss,ifitm,nndr=1
    real(kindreal),save:: fmlos,fitm,deltal,deltat

! **** Convection-linked parameters
  namelist /ConvectionParams/iledou,idifcon,iover,elph,my,dovhp,iunder,dunder
    integer,save:: iledou=0,idifcon=0,my,iover=1,iunder=0
    real(kindreal),save:: elph,dovhp=0.30d0,dunder=0.0d0

! **** Convergence-linked parameters
  namelist /ConvergenceParams/gkorm,alph,agdr,faktor,dgrp,dgrl,dgry,dgrc,dgro,dgr20,nbchx,nrband
    integer,save:: nbchx=200,nrband=1
    real(kindreal),save:: gkorm,alph,agdr,faktor,dgrp,dgrl,dgry,dgrc,dgro=0.010d0,dgr20=0.010d0

! **** Timestep controle
  namelist /TimeControle/xcn,islow,icncst
    integer,save:: islow,icncst=0
    real(kindreal),save:: xcn

! **** Other controles
  namelist /VariousSettings/plot,refresh,iauto,iprn,ia,iout,itmin
    integer,save:: iauto,iprn=99,ia=1,iout=0,itmin=5
    logical,save:: plot,refresh

  integer,parameter:: imagn_default=0,ianiso_default=0,ipop3_default=0,ibasnet_default=0,iopac_default=3,ikappa_default=5, &
    istati_default=0,igamma_default=0,nndr_default=1,iledou_default=0,iunder_default=0,nbchx_default=200, &
    nrband_default=1,icncst_default=0,iprn_default=99,ia_default=1,iout_default=0,itmin_default=5
  real(kindreal),save:: fenerg_default=1.0d0,richac_default=1.0d0,zsol_default=1.40d-2,frein_default=0.0d0,vwant_default=0.0d0, &
    xfom_default=1.0d0,dunder_default=0.0d0,dgr20_default=0.010d0

  public:: Write_param
  public
  private:: Write_param_int,Write_param_real,Write_param_logical
  private:: imagn_default,ianiso_default,ipop3_default,ibasnet_default,iopac_default,ikappa_default,istati_default,&
    igamma_default,nndr_default,iledou_default,iunder_default,nbchx_default,nrband_default, &
    icncst_default,iprn_default,ia_default,iout_default,itmin_default,fenerg_default,richac_default,zsol_default,frein_default, &
    vwant_default,xfom_default,dunder_default,dgr20_default

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
    write(Unit,'(1x,a,d10.3)') trim(x_name),x_in
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
  call Write_param(Unit,"imagn=",imagn,imagn_default)
  write(Unit,'(1x,a,i0)') "ialflu=",ialflu
  call Write_param(Unit,"ianiso=",ianiso,ianiso_default)
  call Write_param(Unit,"ipop3=",ipop3,ipop3_default)
  call Write_param(Unit,"ibasnet=",ibasnet,ibasnet_default)
  write(Unit,'(1x,a,i0)') "phase=",phase
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
  write(Unit,'(1x,a,f7.5)') "rapcrilim=",rapcrilim
  call Write_param(Unit,"vwant=",vwant,vwant_default)
  call Write_param(Unit,"xfom=",xfom,xfom_default)
  write(Unit,'(1x,a,es21.15)') "omega=",omega
  write(Unit,'(1x,a,f5.3,2(a,i0))') "xdial=",xdial,", idialo=",idialo,", idialu=",idialu
  write(Unit,'("&END"/)')

  write(Unit,'(a)') "&SurfaceParams"
  write(Unit,'(1x,a,i0,a,d10.3)') "imloss=",imloss,", fmlos=",fmlos
  write(Unit,'(1x,a,i0,a,f11.9)') "ifitm=",ifitm,", fitm=",fitm
  write(Unit,'(1x,2(a,f7.5))') "deltal=",deltal,", deltat=",deltat
  call Write_param(Unit,"nndr=",nndr,nndr_default)
  write(Unit,'("&END"/)')

  write(Unit,'(a)') "&ConvectionParams"
  call Write_param(Unit,"iledou=",iledou,iledou_default)
  write(Unit,'(1x,a,i0)') "idifcon=",idifcon
  write(Unit,'(1x,a,f0.3,a,i0)') "elph=",elph,", my=",my
  write(Unit,'(1x,a,i0,a,f5.3)') "iover=",iover,", dovhp=",dovhp
  call Write_param(Unit,"iunder=",iunder,iunder_default)
  call Write_param(Unit,"dunder=",dunder,dunder_default)
  write(Unit,'("&END"/)')

  write(Unit,'(a)') "&ConvergenceParams"
  write(Unit,'(1x,2(a,f5.3))') "gkorm=",gkorm,", alph=",alph
  write(Unit,'(1x,a,d9.2,a,es9.2)') "agdr=",agdr,", faktor=",faktor
  write(Unit,'(1x,2(a,f5.3))') "dgrp=",dgrp,", dgrl=",dgrl
  write(Unit,'(1x,3(a,f5.3))') "dgry=",dgry,", dgrc=",dgrc,", dgro=",dgro
  call Write_param(Unit,"dgr20=",dgr20,dgr20_default)
  call Write_param(Unit,"nbchx=",nbchx,nbchx_default)
  call Write_param(Unit,"nrband=",nrband,nrband_default)
  write(Unit,'("&END"/)')

  write(Unit,'(a)') "&TimeControle"
  write(Unit,'(1x,a,i0,a,f0.3)') "islow=",islow,", xcn=",xcnwant
  call Write_param(Unit,"icncst=",icncst,icncst_default)
  write(Unit,'("&END"/)')

  write(Unit,'(a)') "&VariousSettings"
  write(Unit,'(1x,2(a,l2))') "plot=",plot,", refresh=",refresh
  write(Unit,'(1x,a,i2)') "iauto=",iauto
  write(Unit,'(1x,a)') "xyfiles=T"
  call Write_param(Unit,"iprn=",iprn,iprn_default)
  call Write_param(Unit,"ia=",ia,ia_default)
  call Write_param(Unit,"iout=",iout,iout_default)
  call Write_param(Unit,"itmin=",itmin,itmin_default)
  write(Unit,'("&END")')

  return

end subroutine Write_namelist
!=======================================================================
end module inputparam
