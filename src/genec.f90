!>   STELLAR EVOLUTION PROGRAM OF THE GENEVA GROUP
!!
!!  @author A. Maeder, G. Meynet, D. Schaerer, R. Hirschi, S. Ekstrom, C. Georgy
!!  @version  283
!!  @date     mars 2013
!!
!!  @brief Kippenhahn program modified for the effects of rotation, advanced phases, ...
! --------------------------------------------------------------------------
module genec

use io_definitions
use storage, only: GenecStar
use evol,only: kindreal,ldi,mmax,input_dir,npondcouche,npondcoucheAdv
use const,only: um,cst_a,lgLsol,cstlg_sigma,cstlg_G,lgMsol,cst_G,Msol,pi,lgRsol,Rsol,qapicg,xlsomo,year,day,Lsol,cstlg_K1, &
  cstlg_mH,cstlg_k,cst_sigma,Teffsol
use inputparam,only: modanf,nwseq,nzmod,iprn,iauto,ialflu,ianiso,imagn,ipop3,irot,isol,idiff,iadvec,icoeff, &
  igamma,ibasnet,istati,iledou,idifcon,iover,iunder,my,ikappa,iopac,ifitm,itmin,nndr,idialo,idialu,phase,isugi,nbchx, &
  nrband,iout,icncst,islow,ichem,zinit,zsol,z,frein,elph,dovhp,dunder,fmlos,fitm,rapcrilim,omega,xfom,vwant,gkorm,alph, &
  agdr,agds,agdp,agdt,faktor,deltal,deltat,dgrp,dgrl,dgry,dgrc,dgro,dgr20,xdial,fenerg,richac,xcn,idern,display_plot, &
  itminc,idebug,FITM_Change,IMLOSS_Change,INPUTS_Change,Write_namelist,Read_namelist,starname,xyfiles,idebug,&
  bintide,binm2,periodini,eccentricity_ini,verbose,Add_Flux,end_at_phase,end_at_model,iprezams,n_snap,libgenec,imloss, &
  winds_not_applied,prezams_winds_not_applied,ieos,init_synchronized,renorm_abund
use caramodele,only: xLtotbeg,dm_lost,inum,nwmd,xmini,firstmods,eddesc,hh6,glm,xLstarbefHen,hh1,xmdot,rhoc,tc,gls,teff, &
  glsv,teffv,ab,gms,zams_radius,Mdot_NotCorrected,xteffprev,xtefflast,xlprev,xllast,xrhoprev,xrholast,xcprev,xclast,xtcprev,&
  xtclast,modell,nwseqini,radius,xini,is_MS,is_OB,is_RSG,is_WR
use abundmod,only: x,y3,y,xc12,xc13,xc14,xn14,xn15,xo16,xo17,xo18,xf18,xf19,xne20,xne21,xne22,xna23,xmg24,xmg25,xmg26,xal26, &
  xal27,xsi28,xprot,xneut,xbid,xbid1,vx,vy3,vy,vxc12,vxc13,vxc14,vxn14,vxn15,vxo16,vxo17,vxo18,vxf18,vxf19,vxne20,vxne21,vxne22, &
  vxna23,vxmg24,vxmg25,vxmg26,vxal26g,vxal27,vxsi28,vxprot,vxneut,vxbid,vxbid1,ekrote,epote,ekine,erade,snube7,snub8, &
  nbelx,nbzel,nbael,zabelx,abels,abelx,vabelx,mbelx,maxCNO,abundCheck,lcnom,xmcno,scno
use equadiffmod,only: ccg1,ccg2,ccg3,ccz2,ccz3,gkorv,iprc,gkor,iter
use strucmod,only: m,q,p,t,r,s,vp,vt,vr,vs,e,rho,zensi,rprov,ccrad1,NPcoucheEff,id1,id2,drl,drte,dk,drp, &
  drt,drr,rlp,rlt,rlc,rrp,rrt,rrc,rtp,rtt,rtc,chem,ychem,neudr,fitmion,Nabla_mu,vna,vnr,k2_AMC
use rotmod,only: CorrOmega,dlelex,dlelexprev,suminenv,vsuminenv,vvsuminenv,omegi,vomegi,rapcri,xobla,rapom2,fMdot_rot,do1dr,bmomit,&
  btot,btotatm,Flux_remaining,BTotal_EndAdvect,BTotal_StartModel,dlelexsave,timestep_control,xldoex,ivcalc,rrro,ygmoye,vpsi
use timestep,only: alter,dzeitj,dzeit,dzeitv
use convection,only: bordn,jwint,xzc,ixzc,qbc,qmnc,CZdraw,BaseZC,iidraw,drawcon,r_core
use omegamod,only: vcritcalc,omescale,dlonew,omconv,momevo,omenex,om2old,momspe,xjspe1,xjspe2
use envelope,only: dreckf,dreck,notFullyIonised,supraEdd
use ionisation,only: abond,list,iatoms
use diffadvmod,only: tdiff,jdiff
use energy,only: enint,netinit,vmassen,rvect,t9n,pvect,epstot1,epsneut,dcoeff
use geomod, only: rpsi_min,initgeo,geomat,geomeang
use PGPlotModule, only: restart,InitPGplot,SavePlotData,EndPGplot,Chem_Species_Number,PlotEvol,Mass_Vector,Struc_Plotted
use SmallFunc,only: exphi
use LayersShift,only: fitmshift,schrit,mdotshift
use winds,only: aniso,xloss,xldote,corrwind,read_Mdot_prescriptions
use chemicals,only: netnew,chemeps,chemold
use diffusion,only: coedif,diffbr
use timestep,only: zeit,xcnwant,TimestepControle
use henyey_solver,only: henyey,nsugi,correction_message,henyey_last
use opacity,only: ioutable,rout,tout
use nablas,only: grapmui
use PrintAll, only: File_Unit,PrintCompleteStructure
use WriteSaveClose,only: OpenAll,CheckSchrit,write4,read4,SequenceClosing,&
  nzmodini,print_Snapshot,print_files,switch_outputfile,nzmodnew
use bintidemod,only: period, eccentricity, compute_k2_from_structure
use EOS,only: read_helm_table
use safestop, only: safe_stop

implicit none

real(kindreal):: bibib,bolm,dmneed,eddesm=0.0d0,fmain,glsvv,grav,h1,h2,hr, &
  opaesc,rapomm=0.0d0,raysl,teffvv=0.d0,teffel,teffpr,vcrit1=0.0d0,tzero,vcri2m=0.0d0, &
  vcri1m=0.0d0,vequat,vcrit2=0.0d0,vequam=0.0d0,xdilto,xdilex,xltof,xltod,xltot, &
  xmdotneed,xmdotwr,xo1,xtt,xtod2,zwi1,xdippp,zwi,rhocprev,Tcprev

integer:: i,ll,ii,iprnv,iterv,k,j,imlosssave

integer:: Iteration48,IterTriangle,ielemneg

real(kindreal):: summas

real(kindreal), dimension(5):: xnetalu
real(kindreal), dimension(npondcouche):: CorrZero
real(kindreal), dimension(Chem_Species_Number):: Species_PGplot
real(kindreal) :: Z_want, Z_current
character(*), parameter:: headx='                     mass                  radius             temperature                 &
  &density                pressure                  energy               eneutrino                  dcoeff                   &
  &zensi             x             y          xc12          xo16         xne20         xne22         xmg24         xsi28     &
  &    xni56', &
  heady='                     mass                  radius             temperature                 density                &
  &pressure                  energy               eneutrino                  dcoeff                   zensi             &
  &x            y3             y          xc12          xc13          xn14          xn15          xo16          xo17          &
  &xo18         xne20         xne22         xmg24         xmg25         xmg26         xsi28          xs32         xar36         &
  &xca40         xti44         xcr48         xfe52         xni56'

logical:: elemneg,checkVink=.true.,veryFirst,TriangleIteration,snap_printed

namelist/IniStruc/gms,alter,gls,teff,glsv,teffv,dzeitj,dzeit,dzeitv,summas,ab,m,q,p,t,r,s,vp,vt,vr,vs,x,y3,y,xc12,xc13,xn14,xn15,&
  xo16,xo17,xo18,xne20,xne22,xmg24,xmg25,xmg26,omegi

contains

subroutine initialise_genec
! --------------------------------------------------------------------------
  iprnv = 0
  snap_printed = .false.
  call getenv("GENEC_INPUT_DIR", input_dir)
  write(*,*) 'path to inputs directory:',trim(input_dir)
end subroutine initialise_genec

subroutine read_parameters
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reading the input parameters of the calculation.
! Choice of options.
  call Read_namelist
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  call OpenAll
end subroutine read_parameters

subroutine initialise_star
  if (idebug > 0) then
    verbose = .true.
  endif

  if (idebug > 1) then
    write(*,*) 'initialisations...'
  endif
  supraEdd = .false.
  ichem = 0

  rhocprev = 0.d0
  Tcprev = 0.d0

! [Modif CG]
! Initialisation of CorrOmega
  CorrOmega(:) = 0.d0
  xLtotbeg = 0.d0
  dlelex=0.d0
  dlelexprev = 0.d0
  Flux_remaining = 0.d0
  dlelexsave = 0.d0
  BTotal_EndAdvect = 0.d0
  timestep_control = 0.d0
  elemneg = .false.
  ielemneg = 0
  ioutable = 0
! We also initialise a table CorrZero, which consists of npondcouche lines with value zero.
! Intended for calls from momevo without correction.
  CorrZero = 0.d0
  Iteration48 = 1
  IterTriangle = 1
  TriangleIteration = .false.
! Initialisation of suminenv, the moment of inertia of the envelope
  suminenv = 0.d0
  veryFirst = .false.
! [/Modif]
  xteffprev=0.d0
  xtefflast=0.d0
  xlprev=0.d0
  xllast=0.d0
  xrhoprev=0.d0
  xrholast=0.d0
  xcprev=0.d0
  xclast=0.d0
  xcnwant=xcn

!***  IPRN=0   PRINTS ALL THE ITMIN ITERATIONS AND THE LAST ONE FOR EVERY MODEL.
!***  IPRN=1   PRINTS ONLY THE LAST ITERATION FOR EVERY MODEL ALTHOUGH THE ITMIN USUAL ITERATIONS ARE DONE.
!***  IPRN=2   PRINTS THE LAST ITERATION FOR THE 1ST MODEL, THE 3RD MODEL, THE 5TH MODEL AND SO ON.
!              THUS, IT SKIPS PRINTING EVERY 2ND. MODEL.

!***  IN GENERAL WHEN IPRN >= 1, THE PRINTED MODELS WILL BE THE 1ST ONE, THE (IPRN+1)TH ONE, THE (IPRN+2)TH ONE AND SO ON.

!***  IA=0   THE FINAL ATMOSPHERE AND ENVELOPE OF THE EQUILIBRIUM MODEL ARE NOT PRINTED.
!***  IA=1   PRINTS FINAL ATMOSPHERE AND ENVELOPE.C

!***  ALPH IS A FACTOR USED IN HENYEY TO DETERMINE THE CONVERGENCE FACTOR ALPH1 AND WE MUST MAKE SURE THAT ALPH <= 1.

!***  IDER IS USED IN CHEMISTRY FOR THE COMBUSTION OF C12
!***  XCN NOW USED IN ZEIT

  tzero=999999999.d0

  idern=0
  id1=0
  id2=5
  dm_lost=0.d0

  agdp = agdr    ! )
  agds = agdr    ! ) bounds on the corrections in henyey
  agdt = agdr    ! )

  if ((.not. libgenec) .or. (.not. GenecStar%initialised)) then
    dgrp = dgrp*um ! maximum allowed variation in Ln P
    dgrl = dgrl*um ! maximum allowed variation in Ln S
  endif

  write(*,*) 'RESTART AT NWSEQ',nwseq
  if (nwseq == 1) then
    if (idebug > 1) then
      write(*,*) 'initialisation of pgplot, very first run'
    endif
    restart = 0
  else
    if (idebug > 1) then
      write(*,*) 'initialisation of pgplot, continuing with nwseq=',nwseq
    endif
    restart = nwseq
  endif

  if (modanf == 0) then
    if (idebug > 1) then
      write(*,*) 'initial values check and corrections'
    endif
    if (faktor /= 1.d0) then
      faktor = 1.d0
      write(*,*)'First model: faktor put to 1'
    endif
    if (phase /= 1) then
      phase = 1
      write(*,*)'First model: phase put to 1'
    endif
    if (irot == 1 .and. isol /= 1) then
      isol = 1
      write(*,*)'First model: isol put to 1'
    endif
  endif

  if (modanf == 0) then
    if (.not. libgenec) then
      write(io_logs,'(a)') "==========   N E W   S E R I E S   =============="
      call Write_namelist(io_logs,nwseq,modanf,nzmod,xcn)
      write(io_logs,'(a)') "================================================="

      call Write_namelist(io_sfile,nwseq,modanf,nzmod,xcn)
      write(io_sfile,'(a)') "================================================="
    endif
  endif

  if (idebug > 1) then
    write(*,*) 'call netinit'
  endif
  call netinit(z)

  if (idebug > 1) then
    write(*,*) 'Reading of netalu'
  endif

  if (ialflu == 1) then
    if (.not. libgenec) then
    open(unit=io_network,file='netalu.dat')
    read(io_network,*)
    do i=1,5
     read(io_network,'(6x,d23.15)') xnetalu(i)
    enddo
    close(io_network)
    else ! libgenec
     xnetalu = GenecStar%xnetalu
    endif !.not. libgenec
    zabelx=zabelx-xnetalu(1)-xnetalu(2)-xnetalu(3)-xnetalu(4)
  endif

  if (.not. libgenec) then
  write(io_logs,*) z,' ?/= ',zabelx

  if (isugi >= 1 .and. nwseq  ==  1) then
    nsugi=mmax
  endif
  endif ! .not. libgenec

! Import Helmoltz (Rho,T) table in the case the Timmes EOS has been chosen
  if(ieos==1) THEN
    call read_helm_table
  endif

  inum=0
  modell = 1     ! comptage du modele dans la serie courante
  if (nzmod > 1) then
    modell = mod(nwseq,nzmod)     ! comptage du modele dans la serie courante
  else
    modell = 1
  endif
  nzmodini = nzmod
  if (.not. libgenec) then
  nwmd = nwseq   ! number of the first model of the new series
  endif
  nwseqini = nwseq
!=======================================================================
! modanf = 0 : 1st run : reading the structure in the ini_* file.
!        > 0 : Nth run : reading the structure in the .b file.
  if (modanf == 0) then
    inum=0
! security if initial file is missing the iprezams parameter
    if (vwant>epsilon(vwant) .and. iprezams==0) then
      write(*,*) 'VWANT/=0 --> IPREZAMS set to 1'
      iprezams=1
    endif
    if ((init_synchronized .and. bintide) .and. iprezams==0) then
      write(*,*) 'init_synchronized true --> IPREZAMS set to 1'
      iprezams=1
    endif
    if (idebug > 1) then
      write(*,*) 'Reading of initial structure'
    endif
    if (.not. libgenec) then
    read(*,nml=IniStruc)
    endif
    xmini=summas
    zams_radius = 0.d0
    xini = x(1)
    if (prezams_winds_not_applied) then
      winds_not_applied = .true.
    endif
    if (bintide) then
      period = periodini*day
      eccentricity = eccentricity_ini
    endif
    if (irot == 1 .and. isol>=1 .and. omega /= omegi(1)) then
      omegi(:) = omega
    endif
    if (alter == 0.d0) then
      firstmods = .true.
    else
      firstmods = .false.
    endif

    if (fitm /= exphi(q(1))) then
      q(1) = log10(1.d0 - fitm)
    endif

    if (ialflu == 1) then
      xf19(:)=xnetalu(1)
      xne21(:)=xnetalu(2)
      xna23(:)=xnetalu(3)
      xal26(:)=0.d0
      xal27(:)=xnetalu(4)
      xsi28(:)=xnetalu(5)
      xneut(:)=0.d0
      xprot(:)=0.d0
      xc14(:)=0.d0
      xf18(:)=0.d0
      xbid1(:)=0.d0

! Le 28Si est suivi dans abelx(1,ldi).
! Dans xsi28(ldi) on ne va suivre que ses modifications par le cycle Ne-Na-Mg-Al, d'ou la mise a zero suivante:
      xsi28(:)=0.d0
      bibib=1.d0-x(1)-y(1)-y3(1)-xc12(1)-xc13(1)-xn14(1)-xn15(1)-xo16(1)-xo17(1)-xo18(1)-xne22(1)-xmg24(1)-xmg25(1)-xmg26(1)- &
                 xne20(1)-xf19(1)-xne21(1)-xal27(1)-xsi28(1)-xna23(1)

!To get the correct metallicity one should mutliply all metals by Z_want / Z_current.
      do ii=1,nbelx
       bibib=bibib-abels(ii)
      enddo

      xbid(1:m)=bibib
    else
      xf19(:)=0.d0
      xne21(:)=0.d0
      xna23(:)=0.d0
      xal26(:)=0.d0
      xal27(:)=0.d0
      xsi28(:)=0.d0
      xneut(:)=0.d0
      xprot(:)=0.d0
      xc14(:)=0.d0
      xf18(:)=0.d0
      xbid1(:)=0.d0
    endif

    if (renorm_abund) then
      Z_want = 1.d0 - x(1) - y(1) -y3(1)

      Z_current = xc12(1) + xc13(1) +xn14(1)+xn15(1)+xo16(1)+xo17(1)+xo18(1)+xne22(1)+xmg24(1)+xmg25(1)+xmg26(1)+ &
      xne20(1)+xf19(1)+xne21(1)+xal27(1)+xsi28(1)+xna23(1)+ 1.e-16

      do ii=1,nbelx
        Z_current = Z_current + abels(ii)
      enddo
      write(io_logs,*) 'Z_want:',Z_want,' - Z_current:',Z_current
      write(io_logs,*) 'all metal abundances multiplied by ',Z_want/Z_current

      !Correct composition

      xc12(:) = Z_want/Z_current * xc12(:)
      xc13(:) = Z_want/Z_current * xc13(:)
      xn14(:) = Z_want/Z_current * xn14(:)
      xn15(:) = Z_want/Z_current * xn15(:)
      xo16(:) = Z_want/Z_current * xo16(:)
      xo17(:) = Z_want/Z_current * xo17(:)
      xo18(:) = Z_want/Z_current * xo18(:)
      xne22(:) = Z_want/Z_current * xne22(:)
      xmg24(:) = Z_want/Z_current * xmg24(:)
      xmg25(:) = Z_want/Z_current * xmg25(:)
      xmg26(:) = Z_want/Z_current * xmg26(:)
      xne20(:) = Z_want/Z_current * xne20(:)
      xf19(:)  = Z_want/Z_current * xf19(:)
      xne21(:) = Z_want/Z_current * xne21(:)
      xal27(:) = Z_want/Z_current * xal27(:)
      xsi28(:) = Z_want/Z_current * xsi28(:)
      xna23(:) = Z_want/Z_current * xna23(:)
      do ii=1,nbelx
        abels(ii) = Z_want/Z_current * abels(ii)
      enddo
    endif

! for each shell give same value
    zabelx=z
    do ii=1,nbelx
     abelx(ii,:)=abels(ii)
     zabelx=zabelx-abels(ii)
    enddo
    if (ialflu == 1) then
      zabelx=zabelx-xf19(1)-xne21(1)-xna23(1)-xal27(1)
    endif
    write(*,*) 'z,zabelx,m',z,zabelx,m

    if (isugi >= 1) then
      nsugi = m
    endif

! Initialisation du nombre de couche (correction du moment cinetique)
    NPcoucheEff = npondcouche

!  ----------------
!  Remplacer ecriture sur unite 4 par call suivant:
    call write4
!   -----------

    ab=ab*um
    q(:)=q(:)*um
    p(:)=p(:)*um
    t(:)=t(:)*um
    r(:)=r(:)*um
    s(:)=s(:)*um
    vp(:)=vp(:)*um
    vt(:)=vt(:)*um
    vr(:)=vr(:)*um
    vs(:)=vs(:)*um
    veryFirst = .true.

  else ! modanf > 0
!  -----------

    if (.not. libgenec) then
    if (idebug > 1) then
      write(*,*) 'reading .b file'
    endif

! Cas ou modanf > 0
!     Le modele initial est le dernier modele inscrit dans l'unite 'io_bfile_in' apres le run precedent.
!     On lit les parametres d'entree dans l'unite 'io_bfile_in', qui est utilisee pour stocker le dernier modele de chaque serie de
!     calculs.
    read(io_bfile_in) &
            gms,alter,gls,teff,glsv,teffv,dzeitj,dzeit,dzeitv,xmini,ab,dm_lost,m,&
            (q(i),p(i),t(i),r(i),s(i),x(i),y(i),xc12(i),&
            vp(i),vt(i),vr(i),vs(i),xo16(i),vx(i),vy(i),vxc12(i),vxo16(i),i=1,m),&
            drl,drte,dk,drp,drt,drr,rlp,rlt,rlc,rrp,rrt,&
            rrc,rtp,rtt,rtc,tdiff,vsuminenv,&
            (CorrOmega(i),i=1,npondcouche),&
            xLtotbeg,dlelexprev,zams_radius,xini



    read(io_bfile_in) &
            (y3(i),xc13(i),xn14(i),xn15(i),xo17(i),xo18(i),vy3(i),vxc13(i),&
            vxn14(i),vxn15(i),vxo17(i),vxo18(i),xne20(i),&
            xne22(i),xmg24(i),xmg25(i),xmg26(i),vxne20(i),&
            vxne22(i),vxmg24(i),vxmg25(i),vxmg26(i),omegi(i),vomegi(i),i=1,m)

    read(io_bfile_in) &
            (xf19(i),xne21(i),xna23(i),xal26(i),xal27(i),xsi28(i),vxf19(i),&
            vxne21(i),vxna23(i),vxal26g(i),vxal27(i),vxsi28(i),xneut(i),xprot(i),&
            xc14(i),xf18(i),xbid(i),xbid1(i),vxneut(i),vxprot(i),vxc14(i),vxf18(i),&
            vxbid(i),vxbid1(i),i=1,m)

    do ii=1,nbelx
     read(io_bfile_in) (abelx(ii,i),vabelx(ii,i),i=1,m)
    enddo


    read(io_bfile_in) xtefflast,xllast,xrholast,xclast,xtclast,inum,id1,imloss


    if (isugi >= 1) then
      read(io_bfile_in) nsugi
    endif

    if (bintide) then
      read(io_bfile_in) period,eccentricity,r_core,vna,vnr,k2_AMC
    endif




    write(3,*) 'A LA LECTURE: '
    write(3,*)'Corr(1), suminenv, xLtotbeg, dlelexprev: ',CorrOmega(1),vsuminenv,xLtotbeg,dlelexprev
    write(io_logs,*) 'A LA LECTURE: '
    write(io_logs,*)'Corr(1), suminenv, xLtotbeg, dlelexprev: ',CorrOmega(1),vsuminenv,xLtotbeg,dlelexprev
    endif ! .not. libgenec
    vvsuminenv = vsuminenv
    if (.not. libgenec) then
    if (bintide) then
      write(io_logs,*) 'Binary tides, initial and actual period:',periodini,period/day
    endif
    endif
    if (verbose) then
      write(*,*) 'A LA LECTURE: '
      write(*,*)'Corr(1), suminenv, xLtotbeg, dlelexprev: ',CorrOmega(1),vsuminenv,xLtotbeg,dlelexprev
      if (bintide) then
        write(*,*) 'Binary tides, initial and actual period:',periodini,period/day
      endif
    endif

! Initialisation du nombre de couche (correction du moment cinetique).
! Si phase  = 1, on corrige sur 200 couches. Sinon, seulement sur 20.
    if (phase > 1 .and. CorrOmega(npondcouche) < -100.d0) then
      NPcoucheEff = npondcoucheAdv
    else
      NPcoucheEff = npondcouche
    endif

    if (irot /= 0) then
! [ModifCG]
! Test pour savoir si FITM a change par rapport au modele precedent.
      if (idebug > 1) then
        write(*,*) 'call momevo'
      endif
      call momevo(r,vomegi,xltod,CorrOmega,.true.)
! Si FITM a change lors du passage de la sequence, il FAUT un appel supplementaire a momevo pour tenir compte de ce fait.
! Cet appel ANNULE l'appel suivant de momevo(que l'on ne peut pas completement deplacer ici car lors des modeles suivants,
! le programme ne passe plus ici). On calcule donc l'ancien FITM a l'aide de q(i) (qui n'a pas encore ete modifie pour
! prendre en compte le changement de FITM et on le compare avec le nouveau FITM, qui est lu dans l'unite 5.
! checkFITM est utilise pour ne pas faire un double appel de momevo plus loin.
    endif

    if (x(1) > 7.d-1) then
      if (x(m) > (xini - 5.d-3)) then
        firstmods = .true.
      else
        firstmods = .false.
      endif
    else
      firstmods = .false.
    endif
! [/Modif]

!> Pour augmenter progressivement le taux de rotation a la valeur voulue sur la ZAMS
    if (irot==1 .and. isol>=1 .and. (abs(vwant)>1.0d-5 .or. (init_synchronized .and. bintide))) then
      omegi(1:m)=sqrt(xfom)*omegi(1:m)
    endif

    if (.not. libgenec) then
    call write4

    if (idebug > 1) then
      write(*,*) 'call fitmshift'
    endif
    call fitmshift
    endif ! .not. libgenec

  endif ! modanf

  if (.not. libgenec) then
! PGplot initialisation
  call InitPGplot
  endif

! ftfp initialisation
  call initgeo

! reading Mdot_recipes.dat
  call read_Mdot_prescriptions

  if (ialflu==0 .and. xmini<=9.d0) then
    ichem = 1
  endif

! Ecriture du modele initial approximatif
  if (.not. libgenec) then
  write(io_logs,'(//1x,a,i6//1x,a,f8.4,9x,a,1pe13.5,4x,a,1pe9.2,3x,a,0pf8.0/45x,a,1pe8.2,3x,a,0pf7.0)') &
    'modele initial',nwseq-1,'gms=',gms,'alter=',alter,'GLS=',gls,'TEFF=',teff,'GLSV=',glsv,'TEFFV=',teffv

  if (verbose) then
    write(io_logs,'(/4x,"j",5x,"q",7x,"p",8x,"t",8x,"r",8x,"s",9x,"vp",7x,"vt",7x,"vr",7x,"vs",5x,"x",5x,"y3",6x,"y",5x,"xc12",5x,&
      &"xc13",4x,"xn14"/8x,"omega",31x,"xn15",5x,"xo16",6x,"xo17",5x,"xo18",6x,"xne20",11x,"xne22",11x,"xmg24",4x,"xmg25",5x,&
      &"xmg26"/)')
    write(io_logs,'(1x,i4,f8.4,4f9.4,1x,2f8.4,2f9.4,1x,f7.4,f9.6,f7.4,2e8.2,f9.6/8x,f11.8,24x,e8.1,2x,0p,e8.2,e8.2,1x,e8.2,1x,f9.6,&
      &5x,f9.6,6x,3f9.6)')(i,q(i)/um,p(i)/um,t(i)/um,r(i)/um,s(i)/um,vp(i)/um,vt(i)/um,vr(i)/um,vs(i)/um,x(i),y3(i),y(i),xc12(i), &
      xc13(i),xn14(i),omegi(i),xn15(i),xo16(i),xo17(i),xo18(i),xne20(i),xne22(i),xmg24(i),xmg25(i),xmg26(i),i=1,m)

    if (ialflu == 1) then
      write(io_logs,*)'  q,f19,ne21,na23,al26g,al27,si28,neu,pro,xc14,xf18,bid,bid1 - surf & centre:'
      write(io_logs,&
              '((1x,i4,1x,f9.4,12(1x,e9.3)))')1,q(1)/um,xf19(1),xne21(1),xna23(1),xal26(1),xal27(1),xsi28(1),xneut(1),xprot(1), &
              xc14(1),xf18(1),xbid(1),xbid1(1)
      write(io_logs,&
              '((1x,i4,1x,f9.4,12(1x,e9.3)))')m,q(m)/um,xf19(m),xne21(m),xna23(m),xal26(m),xal27(m),xsi28(m),xneut(m),xprot(m), &
              xc14(m),xf18(m),xbid(m),xbid1(m)
    endif

    write(io_logs,*)'    i,nbelx,abelxi - surf & centre:'
    write(io_logs,'(1x,i4,1x,i3,12(1x,e9.3))') 1,nbelx,(abelx(i,1),i=1,nbelx)
    write(io_logs,'(1x,i4,1x,i3,12(1x,e9.3))') m,nbelx,(abelx(i,m),i=1,nbelx)
  endif
  endif

end subroutine initialise_star

subroutine evolve
!******************* Model computation loop ************************
  do
   if (.not.TriangleIteration) then
     xmdot = 0.d0
! Age > 0
     if (.not.veryFirst) then
       if (mod(nwmd,10)==1) then
         call Mass_Vector
         if (display_plot) then
           call PlotEvol
         endif
       endif
       alter=alter+dzeitj   ! dzeitj : evolutionary timestep in years
       if (alter /= dzeitj) then
! To gradually increase the rotation rate
         if (irot==1 .and. isol==1 .and. (abs(vwant)>1.0d-5 .or. (init_synchronized .and. bintide))) then
           omegi(1:m)=sqrt(xfom)*omegi(1:m)
         endif
         if (idebug > 1) then
           write(*,*) 'call fitmshift'
         endif
         call fitmshift
         glsvv=glsv
         glsv=gls
! gls and teff of the new model are calculated by extrapolation from the glsv and teffv values of the previous model
         gls=exp((log(gls))+((log(gls))-log(glsvv))*dzeit/dzeitv)
         teffvv=teffv
         teffv=teff
         if (verbose) then
           write(io_logs,*) 'MAIN **** previous teff,teffvv,dzeit,dzeitv: ',log10(teff),log10(teffvv),dzeit,dzeitv
         endif
         teff= exp((log(teff))+(log(teff)-log(teffvv))*dzeit/dzeitv)
         if (verbose) then
           write(io_logs,*) 'extrapolated teff: ',log10(teff)
         endif
         if (log(teff)<0.d0) then
           write(*,*) 'teff<0 in main: teff,teffvv ',log(teff),log(teffvv)
           call safe_stop('teff<0 in main')
         endif

! computation of the Eddington factor (diffusion by free e-)
!    opaesc: opacity of free electrons diffusion cm^2/g
!    qapicg: 4pi c G
!    xlsomo: Lsol/Msol
         opaesc=0.2d0*(1.d0+x(1))

         eddesc=1.d0/qapicg*opaesc*gls/gms*xlsomo


         if (irot==1 .and. omegi(1)>1.d-15) then
           ivcalc = .true.
         else
           ivcalc = .false.
         endif

         if (idebug > 1) then
           write(*,*) 'call VcritCalc'
         endif
         call VcritCalc(ivcalc,vcrit1,vcrit2,vequat)

! values before convergence (written in the .g file)
! will be recomputed after convergence later
         vcri1m=vcrit1
         vcri2m=vcrit2
         eddesm=eddesc
         vequam=vequat
         rapomm=rapom2

       endif   !   alter /= dzeitj
     endif   !   not veryFirst
!---------------- autre entree pour prochain modele --------------------
!443 continue
     if (.not. libgenec) then
       write(io_logs,'(a)') "#################################################"
       write(io_logs,'("New timestep, model",i6)') nwmd
       write(io_logs,'(a)') "#################################################"
     endif

     if (.not.veryFirst) then
       if (irot /= 0) then
! [Modif CG]
! On ne souhaite appliquer la correction pour la conservation du moment cinetique que lorsque la diffusion (beaucoup plus robuste)
! est appliquee. Dans le cas ou le modele courant est un modele "advecte", alors on conserve comme moment cinetique initial
! du modele precedent.
! Ici on appelle momevo avec r et vomegi, car r est encore le rayon du pas de temps precedent
! alors que omegi a deja ete change, et celui du pas de temps precedent est vomegi
         if (idebug > 1) then
           write(*,*) 'call momevo with false'
         endif
         call momevo(r,vomegi,BTotal_StartModel,CorrOmega,.false.)

         if (iadvec == 0 .or. ((mod(nwmd,2) == 0 .or. xltotbeg < 1.d0).and. .not. elemneg)) then
           xltotbeg = BTotal_StartModel
         endif
         if (.not. libgenec) then
           write(io_logs,*) 'XLTOTBEG: ', xltotbeg
         endif
! [/Modif]
       endif
     endif

     write(*,*)'#################################################',nwmd
     write(*,*)'Model ',nwmd
     write(*,*)'#################################################',nwmd
     write(*,*)'    age=',alter,' gms= ',gms,' m= ',m
     write(*,'(a,f9.6,a,f9.6)') '    Teff = ',log10(teff),'     L = ',log10(gls)
     if (.not. libgenec) then
     write(io_logs,&
             '(a,f11.6,7x,a,1pe13.5,4x,a,0pf8.0,a,f8.0/46x,a,f8.0,a,f7.0//23x,a,1pe10.3,6x,a,e11.3/46x,a,1pe10.3)') ' gms=',gms, &
             'alter=',alter,'gls=',gls,'  teff=',teff,'glsv=',glsv,'  teffv=',teffv,'dzeitj=',dzeitj,'dzeit=',dzeit,'dzeitv=',dzeitv
     endif

! On initialise la densite centrale du precedent modele.
     if (.not.veryFirst) then
       rhocprev = rhoc
       Tcprev = Tc
! [Modif CG]
! Afin d'eviter les boucles infinies lorsque l'on sort du triangle lors de l'integration vers la surface
! (boucle 48 -> goto 48 -> 48, ...), on introduit un test supplementaire comptant ces iteration et en autorisant un
! maximum de 20. Le parametre est Iteration48. On l'initialise ici.
       Iteration48 = 1
! [/Modif]

! Impression d'un message si l'on est sorti des tables  d'opacite pendant le calcul du dernier modele.
       if (ioutable >= 1) then
         write(*,'(1x,a,i5,a,f6.2,a,f8.2)')'Sortie des tables ',ioutable,' fois avec: log(rho) = ',&
                  3.d0*log10(tout)+rout,' et logT = ',log10(tout)+6.d0
         ioutable = 0
       endif
!++----------------------------------------------------------------------
! Detection de spikes dans combustion de He:
       if (phase == 2) then
         xcprev=xclast
         xclast=y(m)
         if (nwmd > nwseq) then
           if (xclast > xcprev*1.05d0 .and. xclast > 0.08d0 .and. verbose) then
             print *,'Spike de ',100.d0*(xclast-xcprev)/xcprev,' % detecte !'
           endif
           if (xclast > xcprev .and. verbose) then
             print *,'Spike de ',100.d0*(xclast-xcprev)/xcprev,'%','         Yc=  ',xclast
           endif
         endif
       endif
!++-----------------------------------------------------------------------
     endif   !   not veryFirst
     henyey_last = .false.

     hh6=s(1)+log(faktor)
     do i=1,m
      zwi1=abs(s(i))-hh6
      zwi=0.d0
      if (zwi1 > -88.d0) then
        zwi=exp(zwi1)
      endif
 ! sign(A,B): value of A with the sign of B
      s(i)=log(1.d0+sign(zwi,s(i)))
      if (s(i) == 0.d0 .and. i /= m) then
        s(i) = 1.d-30
      endif
      zwi1=abs(vs(i))-hh6
      zwi=0.d0
      if (zwi1 > -88.d0) then
        zwi=exp(zwi1)
      endif
      vs(i)=log(1.d0+sign(zwi,vs(i)))
      if (vs(i) == 0.d0 .and. i /= m) then
        vs(i) = 1.d-30
      endif
     enddo

! [ModifCG]
! initialisation de dmneed et xmdotneed:
     dmneed=0.d0
     xmdotneed=0.d0
     Mdot_NotCorrected = 0.d0
! [/Modif]

     if (.not. winds_not_applied) then
       call xloss
       dm_lost=-xmdot*dzeit/year
       if (.not. libgenec) then
         write(io_logs,*) 'dm= ',dm_lost
       endif
       gms=gms+dm_lost
       write(*,*) 'GMS AFTER WINDS CALC:',gms
       if (.not. libgenec) then
         write(io_logs,*) 'GMS AFTER WINDS CALC:',gms
       endif
     else
       dm_lost = 0.d0
       write(*,*) 'NO WINDS, GMS:',gms
       if (.not. libgenec) then
         write(io_logs,*) 'NO WINDS, GMS:',gms
       endif
     endif

! BEFORE CALLING HENYEY, STORE PREVIOUS ABUNDANCES FOR APPLICATION OF THE IMPLICIT METHOD OF ITERATION ON ABUNDANCES IN SUB.
! NETWKI (NETWKI WILL BE CALLED WITHIN HENYEY).
     if (alter <= dzeitj) then
       vx(1:m)=x(1:m)
       vy3(1:m)=y3(1:m)
       vy(1:m)=y(1:m)
       vxc12(1:m)=xc12(1:m)
       vxc13(1:m)=xc13(1:m)
       vxn14(1:m)=xn14(1:m)
       vxn15(1:m)=xn15(1:m)
       vxo16(1:m)=xo16(1:m)
       vxo17(1:m)=xo17(1:m)
       vxo18(1:m)=xo18(1:m)
       vxne20(1:m)=xne20(1:m)
       vxne22(1:m)=xne22(1:m)
       vxmg24(1:m)=xmg24(1:m)
       vxmg25(1:m)=xmg25(1:m)
       vxmg26(1:m)=xmg26(1:m)
       if (ialflu == 1) then
         vxc14(1:m)=xc14(1:m)
         vxf18(1:m)=xf18(1:m)
         vxf19(1:m)=xf19(1:m)
         vxne21(1:m)=xne21(1:m)
         vxal26g(1:m)=xal26(1:m)
         vxal27(1:m)=xal27(1:m)
         vxsi28(1:m)=xsi28(1:m)
         vxna23(1:m)=xna23(1:m)
       endif
       vabelx(1:nbelx,1:m)=abelx(1:nbelx,1:m)
       vomegi(1:m)=omegi(1:m)
     endif

! The following subroutine account for the effects due to the stellar wind anisotropies, its output is the computation of Lexcess,
! the difference between the angular momentum lost supposing that mass is lost isotropically and the angular momentum lost
! when the anisotropies of the winds are accounted for.
     if (irot == 1) then
       if (idebug > 1) then
         write(*,*) 'call aniso'
       endif
       call aniso(1.d0/xobla,ygmoye,rrro)

! [Modif CG]
! On corrige ici (APRES le calcul de l'anisotropie) la vitesse de rotation de la derniere couche,
! si la diffusion est appliquee a ce modele.
       if (isol == 0) then
         if (iadvec == 0 .or. (mod(nwmd,2) == 0 .and. rapcrilim > 1.d-5)) then
           if (verbose) then
             write(*,'(3(a,d14.8))') &
               'APPLICATION DE LA CORRECTION: xltotbeg: ',xLtotbeg,'omegi(1): ',omegi(1),'CorrOmega(1): ',CorrOmega(1)
           endif
           do i=1,NPcoucheEff
            vomegi(i) = vomegi(i) + CorrOmega(i)
            omegi(i) = omegi(i) + CorrOmega(i)
            if (omegi(i) <= 0.d0) then
              omegi(i) = 1.d-20
            endif
            if (vomegi(i) <= 0.d0) then
              vomegi(i) = 1.d-20
            endif
           enddo
           if (phase > 1 .and. CorrOmega(npondcouche) > -100.d0) then
             NPcoucheEff = npondcoucheAdv
             CorrOmega = 0.d0
             CorrOmega(npondcouche) = -200.d0
             write(*,*) 'NPcoucheEff set to ', NPcoucheEff
           endif
         endif ! iadvec
       endif ! isol == 0
! [/Modif]
       xo1=omegi(1)

       if (idebug > 1) then
         write(*,*) 'call xldote'
       endif
       call xldote(dm_lost,dmneed)
! [ModifCG]
       if (iprezams >= 1) then
         dmneed = 0.d0
         xldoex = 0.d0
       endif
! [/Modif]

       do1dr=(xo1-omegi(2))/(exp(r(1))-exp(r(2)))
       if (verbose .or. dmneed /= 0.d0) then
         write(*,*) 'gms, dmneed:', gms, dmneed
       endif
       gms=gms+dmneed
     endif

! [ModifCG]
     if (abs(dm_lost)>epsilon(dm_lost) .or. abs(dmneed)>epsilon(dmneed)) then
       if (idebug > 1) then
         write(*,*) 'call MdotShift'
       endif
       call MdotShift(dmneed)
       if (xmdot > 0.d0) then
         xmdot=log10(xmdot)
       else
         xmdot = -30.d0
       endif
       if (.not. libgenec) then
         write(io_logs,'(//,2x,a,f13.8,2(1x,a,e14.7),1x,a,f8.3//)') 'gms=',gms,'dm=',dm_lost,'dmneed=',dmneed,'mdot=',xmdot
       endif
     endif
     if (irot == 1) then
       if (dmneed /= 0.d0) then
         xmdotneed=log10(-dmneed*year/dzeit)
       endif
       if (idebug > 1) then
         write(*,*) 'call momevo'
       endif
       call momevo(r,vomegi,xtod2,CorrZero,.true.)
! On additionne ici la contribution du modele precedent a la perte actuelle.
! dlelexprev est nul lorsque c'est necessaire. Pour l'impression fichier, on conserve une sauvegarde de dlelex effectif
! (du seul modele en cours).
       dlelexsave = dlelex
       dlelex = dlelex + dlelexprev
       if (.not. libgenec) then
         write(io_logs,*) 'dlelex, dlelexprev: ', dlelex,dlelexprev
       endif
! [/Modif]
     endif

! Constantes utilisees
! dans les equations aux differences finies G1,2,3 pour l'interieur,
! dans les conditions limites Z1,2,3 au centre.
     glm=log(gms)+log(Msol)                                                ! Ln(masse de l'etoile, en grammes)
     ccrad1=-glm+log(3.d0/(4.d0*qapicg*cst_a))                             ! -Ln(M*)+Ln[3/(16.pi.a.c.G)] pour le gradient radiatif.
     ccg1=glm
     ccg2=2.d0*glm+log(cst_G/(4.d0*pi))                                    ! Ln(GMM/4pi)
     ccg3=glm+log(1.d0/(4.d0*pi))                                          ! Ln(M/4pi)
     ccz2=(2.d0/3.d0)*glm+log(1.d0/2.d0*(4.d0*pi/3.d0)**(1.d0/3.d0)*cst_G) ! Ln[0.5(4.pi/3)E+0.33*G*(M)E+0.66]
     ccz3=glm+log(3.d0/(4.d0*pi))                                          ! Ln(3M/4pi)

     call CheckSchrit("avant")
     if (idebug > 1) then
       write(*,*) 'call schrit'
     endif
     call schrit
     call CheckSchrit("apres")

     do j=1,m-1
      e(j)=exphi(0.5d0*(q(j)+q(j+1)))
     enddo

! Extrapolation du modele initial approche a partir du modele precedent.
     if (alter > dzeitj) then
       fmain=dzeit/dzeitv
     else
       fmain=0.d0
     endif

     vp(1:m)=vp(1:m)*fmain
     vt(1:m)=vt(1:m)*fmain
     vs(1:m)=s(1:m)
     p(1:m)=p(1:m)+vp(1:m)
     t(1:m)=t(1:m)+vt(1:m)

     do j=1,m
      hr=r(j)+(r(j)-vr(j))*fmain
      vr(j)=r(j)
      r(j)=hr
     enddo
     s(m)=0.d0
     gkorv=0.d0
     iterv=0
     if (dm_lost /= 0.d0) then
       id1=0
     endif
!-----------------------------------------------------------------------
!  MODIFICATIONS MAI/JUIN 1990, D.SCHAERER:
!    IONISATION PARTIELLE POUR PLUSIEURS ELEMENTS...
!  * A cet endroit les abondances utilisees par la sous-routine IONPART et SAHA sont mises a jour,
!    c.a.d. a la valeur correspondante a la couche exterieure du noyau.
!  * En meme lieu on decide des elements que l'on veut traiter. Ceci se fait selon leurs abondances en surface.
!    Pour traiter un element (des 6 possibles) il suffit de mettre son num dans la liste.
!    Pour terminer, si on n'utilise pas tous les elements, on met '0' dans la liste.

!  LA CONSTANTE 'IATOMS' defini le nombre d'elements traites au maximum.
!  LA CONSTANTE 'IONSTATES' defini le nombre max. d'etats d'ionisation que l'on traite, donc normalement le NOMBRE ATOMIQUE

!  LES ABONDANCES EN H, He, C, O, Ne, Mg dans cet ordre:
     chem=x(1)
     ychem=y(1)
     abond(1)=x(1)
     abond(2)=y(1)
     abond(3)=xc12(1)
     abond(4)=xo16(1)
     abond(5)=xne20(1)
     abond(6)=xmg24(1)

     j=0
     do i=1,iatoms
      if (abond(i) > 2.d-2) then
        j=j+1
        list(j)=i
      endif
     enddo
     if (j /= iatoms) then
       list(j+1)=0
     endif
!***********************************************************************
!            Fin de l'initialisation. Calcul du modele.
!            ******************************************
! Calcul des sommets du triangle qui contiennent le point approximatif du diagramme HR.
     if (idebug > 1) then
       write(*,*) 'call dreck'
     endif
     call dreck(0)
! Apres dreck, on a :       neudr = 0 : conditions initiales inchangees,
!                           neudr = 1 : conditions limites a recalculer.

   endif   !   not TriangleIteration

   omega=omegi(1)
   if (TriangleIteration) then
     TriangleIteration = .false.
   endif
   if (neudr == 1 .or. id1 /= 1) then
     if (idebug > 1) then
       write(*,*) 'call dreckf'
     endif
     call dreckf
   endif
   if (nwmd == 1) then
     vvsuminenv = suminenv
   endif
! Calcul, par dreckf, des nouvelles conditions limites, c'est a dire determination de
! (alpha i, beta i, gamma i) (i=1,2,3) des equations (27,28,29).

   if (id1 /= 2) then
     id1=1
   endif
! iprn : Tous les iprn modeles, on imprime le modele complet.
! iprc : Lorsque iprc=1, le dernier appel de Henyey (itminc=1) se fait en 524. On imprime.
!        Lorsque iprc=0, le dernier appel de Henyey (itminc=1) se fait en 54. On n'imprime pas.

   if (iprn == 0) then
     iprc=1
   else
     iprc=0
   endif

! itmin : Nombre minimum d'iterations a effectuer dans Henyey lorsque toutes les corrections sont deja inferieures a gkorm,
!         afin de prevenir une divergence.
   itminc=itmin

! Initialisation du gradient de poids moleculaire.
   Nabla_mu(:) = 0.d0
!   if (iledou == 1 .or. irot == 1 .or. idifcon == 1) then
   if (idebug > 1) then
     write(*,*) 'call grapmui'
   endif
   call grapmui
!   endif
   if (irot == 1 .and. isol == 0) then
     if (idebug > 1) then
       write(*,*) 'call dlonew'
     endif

     call dlonew
   endif

! Before Henyey, we call once again momevo, to have the total angular momentum
! difference between the initial one and the present one.
   if (irot == 1) then
     if (idebug > 1) then
       write(*,*) 'call momevo'
     endif
! NB: xLstarbefHen est redetermine dans omenew plus tard. Celui-ci n'est pas utilise.
! Les deux valeurs sont identiques.
     call momevo (vr,vomegi,xlstarbefHen,CorrZero,.true.)
     xlstarbefHen = xlstarbefHen*1.d53
   endif

!  -----------
   if (idebug > 1) then
     write(*,*) 'call henyey 1'
   endif
   Flux_remaining = 0.d0
   call henyey
!  -----------
   if (gkor >= gkorm) then
     if (alter <= dzeitj) then
       write(*,*) 'worst gkor = ', gkor
       stop 'bad initial structure'
     endif

! gkorm : Valeur absolue de la correction maximale toleree.
! gkor  : Valeur absolue de la correction maximale calculee dans Henyey.
! Si gkor > gkorm : Il faut revenir en arriere.
!    S'il existe un modele precedent (modanf > 0), on en repart avec un pas de temps dzeit/2.
!    On extrapole un modele approche a (t+dzeit/2) que l'on corrige avec Henyey.
!    S'il n'existe pas de modele precedent (modanf=0), on cherche un nouveau triangle et on recalcule le modele initial.

! Si gkor > gkorm :
     write(*,*) '!!!   ELEMENT NEGATIF   !!!'
     write(*,*) correction_message
     if (iauto >= 2) then
       if (verbose) then
         write(*,*) 'OLD: ',gkorm,faktor,alph
       endif
       if (ielemneg > 2) then
         if ((gkorm<0.5d0 .and. phase<=2) .or. (gkorm<1.0d0 .and. phase>=3) .or. (gkorm<1.5d0 .and. phase>=5)) then
           gkorm = gkorm + 0.1d0
           write(io_input_changes,'(i7.7,a8,f5.2)') nwmd,': GKORM=',gkorm
         endif
         if (phase <= 3 .and. faktor < 10.d0**(0.5d0*real(phase))) then
           faktor = 2.d0*real(phase)*faktor
           write(io_input_changes,'(i7.7,a9,1pd9.2)') nwmd,': FAKTOR=',faktor
         else if (phase > 3 .and. faktor < 10.d0**(1.5d0 + 3.d0*(real(phase)-3.d0))) then
           faktor = 6.d0*faktor
           write(io_input_changes,'(i7.7,a9,1pd9.2)') nwmd,': FAKTOR=',faktor
         endif
         if ((phase >= 2 .and. alph > 0.8d0) .or. (phase >= 5 .and. alph > 0.5d0)) then
           alph = alph - 0.1d0
           write(io_input_changes,'(i7.7,a7,f5.2)') nwmd,': ALPH=',alph
         endif
       endif
       if (verbose) then
         print*,'NEW: ',gkorm,faktor,alph
       endif
     endif
     elemneg = .true.
     ielemneg = ielemneg + 1
     Iteration48 = 1
     IterTriangle = 1
     if (.not. libgenec) then
     write(io_logs,'(//////,10x,a,//////)')'GOING BACK : corrections too big'
     endif
     iprnv= iprnv - 1

     modell=modell-1
     nwmd=nwmd-1

     call read4

     dzeitj = dzeitj/2.d0
     if (phase < 3 .and. dzeitj < 1.d-4) then
       rewind(io_runfile)
       write(io_runfile,*) nwmd,': time step too small'
       stop 'time step too small'
     endif
     jdiff=2
     dzeit=dzeit/2.d0
     dm_lost = dm_lost/2.d0

     x(1:m)=(x(1:m)+vx(1:m))/2.d0
     y(1:m)=(y(1:m)+vy(1:m))/2.d0
     y3(1:m)=(y3(1:m)+vy3(1:m))/2.d0
     xc12(1:m)=(xc12(1:m)+vxc12(1:m))/2.d0
     xc13(1:m)=(xc13(1:m)+vxc13(1:m))/2.d0
     xn14(1:m)=(xn14(1:m)+vxn14(1:m))/2.d0
     xn15(1:m)=(xn15(1:m)+vxn15(1:m))/2.d0
     xo16(1:m)=(xo16(1:m)+vxo16(1:m))/2.d0
     xo17(1:m)=(xo17(1:m)+vxo17(1:m))/2.d0
     xo18(1:m)=(xo18(1:m)+vxo18(1:m))/2.d0
     xne20(1:m)=(xne20(1:m)+vxne20(1:m))/2.d0
     xne22(1:m)=(xne22(1:m)+vxne22(1:m))/2.d0
     xmg24(1:m)=(xmg24(1:m)+vxmg24(1:m))/2.d0
     xmg25(1:m)=(xmg25(1:m)+vxmg25(1:m))/2.d0
     xmg26(1:m)=(xmg26(1:m)+vxmg26(1:m))/2.d0
     if (ialflu == 1) then
       xc14(1:m)=(xc14(1:m)+vxc14(1:m))/2.d0
       xf18(1:m)=(xf18(1:m)+vxf18(1:m))/2.d0
       xf19(1:m)=(xf19(1:m)+vxf19(1:m))/2.d0
       xne21(1:m)=(xne21(1:m)+vxne21(1:m))/2.d0
       xal26(1:m)=(xal26(1:m)+vxal26g(1:m))/2.d0
       xal27(1:m)=(xal27(1:m)+vxal27(1:m))/2.d0
       xsi28(1:m)=(xsi28(1:m)+vxsi28(1:m))/2.d0
       xna23(1:m)=(xna23(1:m)+vxna23(1:m))/2.d0
       xneut(1:m)=(xneut(1:m)+vxneut(1:m))/2.d0
       xprot(1:m)=(xprot(1:m)+vxprot(1:m))/2.d0
       xbid(1:m)=(xbid(1:m)+vxbid(1:m))/2.d0
       xbid1(1:m)=(xbid1(1:m)+vxbid1(1:m))/2.d0
     endif
     abelx(1:nbelx,1:m)=(abelx(1:nbelx,1:m)+vabelx(1:nbelx,1:m))/2.d0
     omegi(1:m)=(omegi(1:m)+vomegi(1:m))/2.d0

   else   !   gkor < gkorm

     if (gkor > gkorv) then
       gkorv=gkor
     endif
     if (iter <= iterv) then
       iterv=iter
     endif
     gls=-exp(hh6-log(Lsol))*exphi(s(1))
     if (verbose) then
       if (.not. libgenec) then
       write(io_logs,*) 'After Henyey, teff untouched=',log10(teff)
       endif
     endif
     teff=exp(rtp*p(1)+rtt*t(1)+rtc)
     if (verbose) then
       if (.not. libgenec) then
       write(io_logs,*) '              teff new=',log10(teff)
       write(io_logs,*) '              rtp,p(1),rtt,t(1),rtc: ',rtp,p(1),rtt,t(1),rtc
       endif
     endif
     write(*,*) "TEFF ESTIMATION: ",log10(teff),log10(gls)
     if (isnan(log10(teff))) then
       rewind(io_runfile)
       write(io_runfile,*) 'teff undefined in main 1159: rtp,rtt,rtc,p(1),t(1) ',rtp,rtt,rtc,p(1),t(1)
       stop 'teff undefined in main 996'
     endif
     if (log10(teff)<3.d0) then
       write(io_logs,*) 'teff<3, set to teffv'
! Instead of crashing, we simply take the previous value of Teff
!       write(io_runfile,*) 'teff<3 in main 1159: rtp,rtt,rtc,p(1),t(1) ',rtp,rtt,rtc,p(1),t(1)
!       stop 'teff<3 in main 996'
       teff = teffv
     endif
     if (log10(teff)>6.5d0) then
!       rewind(io_runfile)
       write(io_logs,*) 'teff>6.5, set to teffv'
!       write(io_runfile,*) 'teff>6.5 in main 1159: rtp,rtt,rtc,p(1),t(1) ',rtp,rtt,rtc,p(1),t(1)
!       stop 'teff>6.5 in main 996'
       teff = teffv
     endif
     if (idebug > 1) then
       write(*,*) 'call dreck'
     endif
     call dreck(nndr)

! neudr : Initialisation dans dreck.
!         Si neudr=0 : conditions limites inchangees.
!         Si neudr=1 : conditions limites changees.
! nndr  : Parametre d'entree.
! Si neudr=0 ou nndr=0, on ne recalcule pas le triangle, et on a le modele le plus proche. Sinon, on revient en 48.
     if (neudr /= 0 .and. nndr /= 0) then

! [Modif CG]
! On teste ici Iteration48. Si superieur a  20: arret de l'execution et
! affichage d'un message d'erreur.
       if (IterTriangle > 12 .and. iauto == 2) then
         write(*,*) 'Convergence problems in the envelope... Triangle reinitialisation.'
         if (.not. libgenec) then
         write(io_logs,*) 'Convergence problems in the envelope... Triangle reinitialisation.'
         endif
         IterTriangle=0
         id1 = 2
       endif
       if (Iteration48 > 36) then
         write(*,*)
         write(*,*) '!*!*!*!*!*!*!*!*!'
         write(*,*) 'More than 36 iterations in model ',nwmd,':'
         write(*,*) 'convergence in the triangle not reached. Aborting...'
         write(*,*) '!*!*!*!*!*!*!*!*!'
         rewind(io_runfile)
         write(io_runfile,*) nwmd,': Problem with triangle convergence'
         call safe_stop('Problem with triangle convergence')
       endif
!-----------------------------------------------------------------------
       Iteration48 = Iteration48 + 1
       IterTriangle = IterTriangle + 1
! [/Modif]
       TriangleIteration = .true.
       cycle   !   ANCIEN GOTO 48
!-----------------------------------------------------------------------
     endif
     if (elemneg) then
       write(io_input_changes,'(i7.7,a,i2,a)') nwmd,': ',ielemneg,' times ELEM NEG'
       elemneg = .false.
       ielemneg = 0
     endif
     h1=log10(gls)   ! L*/Lsoleil
     h2=(rtp*p(1)+rtt*t(1)+rtc)/um
     radius=1.d0/2.d0*(lgLsol-log10(4.d0*pi)-cstlg_sigma)-lgRsol+0.5d0*h1-2.d0*h2
     grav=log10(4.d0*pi)+cstlg_sigma+cstlg_G+lgMsol-lgLsol+4.d0*h2-h1+log10(gms)
     bolm=4.77d0-2.5d0*h1

! computation of the Eddington coefficient (diffusion by free e-)
! opaesc: opacity for the diffusion by free electrons cm^2/g
! qapicg: 4pi c G
! xlsomo: Lsol/Msol
     opaesc=0.2d0*(1.d0+x(1))
     eddesc=1.d0/qapicg*opaesc*gls/gms*xlsomo

     if (irot == 1 .and. omegi(1)  >  1.d-15) then
       ivcalc = .true.
     else
       ivcalc = .false.
     endif

     write(*,*)'Modele converge'
     if (idebug > 1) then
       write(*,*) 'call VcritCalc'
     endif
     call VcritCalc(ivcalc,vcrit1,vcrit2,vequat)

     if (.not. libgenec) then
     write(io_logs,&
             '(/////a,f7.3,a,f7.4,a,f8.4,a,f6.3,a,f8.4/1x,a,f11.8,a,f12.8)') ' Equilibrium model for log l=',h1,'  logte=',h2, &
             '  log r=',radius,'  log g=',grav,' mbol=',bolm,' omega=',omega,' rapcri=',rapcri
     endif

!-------------- FINAL MODEL -----------
     itminc=1
     if (iprnv > 0) then     !  not printed
       iprc=0
       id2=1
     else                    !  will be printed
       iprc=1
       id2=6
     endif
     if (idebug > 1) then
       write(*,*) 'call dreckf'
     endif
     call dreckf
! Dans ce dreckf, on calcul le modele d'enveloppe apres convergence. On identifie vsuminenv
! a suminenv ici car on a enfin une valeur correcte.
     write(*,*) 'MAIN: vsuminenv=suminenv',vsuminenv,suminenv
     !if (abs(suminenv/vsuminenv -1.d0) < 0.02d0) then
     vsuminenv = suminenv
     !endif
     if (iprnv <= 0) then   ! final model will be printed
       if (irot == 1) then
         if (idebug > 1) then
           write(*,*) 'call momevo'
         endif
         call momevo(r,omegi,xltot,CorrZero,.true.)
       endif
     endif

!      -----------
     if (idebug > 0) then
       write(*,*) 'ini: ', BTotal_StartModel
       write(*,*) 'to lose: ', dlelexsave
       write(*,*) 'to reach: ', BTotal_StartModel-dlelexsave
       write(*,*) 'where we are :', BTotal_EndAdvect
       write(*,*) 'remaining: ', BTotal_StartModel-dlelexsave-BTotal_EndAdvect
     endif
     if (Add_Flux .and. rapcrilim>1.d-5) then
       Flux_remaining = (BTotal_StartModel-dlelexsave-BTotal_EndAdvect)/dzeit
     else
       Flux_remaining = 0.d0
     endif
     if (idebug > 1) then
       write(*,*) 'last call henyey',BTotal_EndAdvect,xltotbeg,BTotal_StartModel,Flux_remaining,dlelexsave
     elseif (idebug > 0) then
       write(*,*) 'last call henyey'
     endif
     henyey_last = .true.
     call henyey
     if (idebug > 1) then
       write(*,*) 'call VcritCalc'
     endif
     call VcritCalc(ivcalc,vcrit1,vcrit2,vequat)
     vvsuminenv = vsuminenv

!      -----------

     if (iprnv <= 0) then   ! iprnv <= 0
! Impression de la structure complete int+env+atm
       if (.not.libgenec) then
       call PrintCompleteStructure
       endif

! y-file similar to x-file but just for printed timesteps, but with the complete set of abundances (complete abelx)
       if (xyfiles) then
         write(io_yfile,'(i7,e23.16,e23.16,i5,4(1pe24.16))') nwmd,alter,gms,m,dzeit,dzeit/year,gls,teff
         write(io_yfile,'(a)') heady
         do i=1,m
          write(io_yfile,'(i5,9(1pe24.16),77(0pe14.7))') i,vmassen(i),rvect(i),t9n(i),exp(rho(i)),pvect(i),epstot1(i),epsneut(i), &
            dcoeff(i),zensi(i),x(i),y3(i),y(i),xc12(i),xc13(i),xn14(i),xn15(i),xo16(i),xo17(i),xo18(i),xne20(i),xne22(i), &
            xmg24(i),xmg25(i),xmg26(i),(abelx(ll,i),ll=1,nbelx)
         enddo
       endif

! iprnv, compteur de modeles imprimes, est reinitialise a iprn
       iprnv=iprn

       if (.not. libgenec) then
       write(io_logs,'(a,/,a)')'centre: m,x,y3,y,xc12,xc13,xn14,xn15,xo16,xo17,xo18','xne20,xne22,xmg24,xmg25,xmg26'
       write(io_logs,'(1x,i5,1p,10e11.3,/5e12.4)') m,x(m),y3(m),y(m),xc12(m),xc13(m),xn14(m),xn15(m),xo16(m),xo17(m),xo18(m), &
         xne20(m),xne22(m),xmg24(m),xmg25(m),xmg26(m)

       if (ialflu == 1) then
         write(io_logs,'(a)')'centre: xf19,xne21,xna23,xal26g,xal27,xsi28'
         write(io_logs,'(1x,1p,6e12.4)')xf19(m),xne21(m),xna23(m),xal26(m),xal27(m),xsi28(m)
       endif
       endif

     endif   ! iprnv

     iprnv=iprnv-1

! mod xfile
! ascii version
     if (xyfiles) then
       write(io_xfile,'(i7,e23.16,e23.16,i5,4(1pe24.16))') nwmd,alter,gms,m,dzeit,dzeit/year,gls,teff
       write(io_xfile,'(a)') headx
       do i=1,m
        write(io_xfile,'(i5,9(1pe24.16),77(0pe14.7))') i,vmassen(i),rvect(i),t9n(i),exp(rho(i)),pvect(i),epstot1(i),epsneut(i), &
          dcoeff(i),zensi(i),x(i),y(i),xc12(i), xo16(i),xne20(i),xne22(i),xmg24(i),abelx(1,i),abelx(8,i)
       enddo
     endif

     call maxCNO(m,q)

     if (irot == 1 .and. isol == 0) then
       if (idebug > 1) then
         write(*,*) 'call omescale'
       endif
       call omescale
     endif

! Abundance check but no measures taken if bad !
     if (idebug > 2) then
       call abundCheck(m,.true.)
     else
       call abundCheck(m,.false.)
     endif

! Determination of the convective zones for the .g file
     call CZdraw

     xtt=log10(teff)
     xteffprev=log10(teffv)

     if (phase > 1) then
       if (idebug > 1) then
         write(*,*) 'call FITM_Change'
       endif
       call FITM_Change(teffvv,fitmIon,m,zensi,q,notFullyIonised,BaseZC)
     endif

! [Modif IMLOSS]
! Passage si necessaire a une perte de masse WR ou supra-Edd pour les massives
! ou a la perte de masse geante rouge pour les petites
     if (imloss/=2 .and. iauto/=0) then
       if (idebug > 1) then
         write(*,*) 'call IMLOSS_Change'
       endif
       call IMLOSS_Change(x(m),supraEdd,vequat,xtt)
     endif
! [/Modif IMLOSS]

     if (alter > 0.d0) then
       dzeitv=dzeit
       hh1=abs(h2-log10(teffv))
       call zeit
     endif

     if (rhocprev /= 0.d0 .and. x(m) < 0.7d0) then
       if (abs(rhoc-rhocprev)/rhoc > 5.d-2 .or. abs(Tc-Tcprev)/Tc > 5.d-2) then
         write(*,*) 'Central density variation over the last time step too large: ',100.d0*abs(rhoc-rhocprev)/rhoc, '%'
         write(*,*) 'of central temperature variation too large: ',100.d0*abs(Tc-Tcprev)/Tc, '%'
         rewind(io_runfile)
         write(io_runfile,*)nwmd,': Variation of central conditions too large'
         stop
       endif
     endif

! Avant de sauvegarder definitivement les valeurs du modele converge dans les variables "v",
! on ajuste une derniere fois le profil du moment cinetique.
     vx(1:m)=x(1:m)
     vy3(1:m)=y3(1:m)
     vy(1:m)=y(1:m)
     vxc12(1:m)=xc12(1:m)
     vxc13(1:m)=xc13(1:m)
     vxn14(1:m)=xn14(1:m)
     vxn15(1:m)=xn15(1:m)
     vxo16(1:m)=xo16(1:m)
     vxo17(1:m)=xo17(1:m)
     vxo18(1:m)=xo18(1:m)
     vxne20(1:m)=xne20(1:m)
     vxne22(1:m)=xne22(1:m)
     vxmg24(1:m)=xmg24(1:m)
     vxmg25(1:m)=xmg25(1:m)
     vxmg26(1:m)=xmg26(1:m)

     if (ialflu == 1) then
       vxf19(1:m)=xf19(1:m)
       vxne21(1:m)=xne21(1:m)
       vxal26g(1:m)=xal26(1:m)
       vxal27(1:m)=xal27(1:m)
       vxsi28(1:m)=xsi28(1:m)
       vxna23(1:m)=xna23(1:m)
       vxneut(1:m)=xneut(1:m)
       vxprot(1:m)=xprot(1:m)
       vxc14(1:m)=xc14(1:m)
       vxf18(1:m)=xf18(1:m)
       vxbid(1:m)=xbid(1:m)
       vxbid1(1:m)=xbid1(1:m)
     endif

     vabelx(1:nbelx,1:m)=abelx(1:nbelx,1:m)

     vomegi(1:m)=omegi(1:m)

     idern=1
     jdiff = 2
     inum=inum+1

! Estimation de la composition chimique au pas temporel suivant
     if (idebug > 1) then
       write(*,*) 'call chemold'
     endif
     call chemold
     if (ichem == 1 .and. idifcon == 0) then
       if (idebug > 1) then
         write(*,*) 'call chemeps'
       endif
       call chemeps
     endif
     if (idebug > 1) then
       write(*,*) 'call netnew'
     endif
     call netnew

     if (irot==1 .and. idiff/=0 .and. isol==0 .or. idifcon==1) then
       if (idebug > 1) then
         write(*,*) 'call coediff'
       endif
       call coedif
       if (idebug > 1) then
         write(*,*) 'call diffbr'
       endif
       call diffbr
     endif

     if (irot == 1 .and. isol == 0) then
       if (idebug > 1) then
         write(*,*) 'call om2old'
       endif
       call om2old
       if (alter > dzeitj) then
         fmain=dzeit/dzeitv
       else
         fmain=0.d0
       endif

       rprov(1:m)=r(1:m)+(r(1:m)-vr(1:m))*fmain

       if (idebug > 1) then
         write(*,*) 'call omenex'
       endif
       call omenex

       if (idifcon == 0) then
         if (idebug > 1) then
           write(*,*) 'call omconv'
         endif
         call omconv
       endif

       if (idebug > 1) then
         write(*,*) 'call momevo'
       endif
       call momevo(r,vomegi,xltof,CorrOmega,.true.)

       xdilto=xltod-xtod2
       xdilex=xtod2-xltof
       xdippp=dlelex/1.d+53
       write(io_sfile,'(1x,i7,2(a,f10.7),2(/1x,a,e12.6),a,f5.2,/1x,a,e12.6)') nwmd,' MOM. ANG. DEB=',xltod,' FIN=',xltof, &
         ' difference due Mdot   ISO=',xdilto,' L exces                  =',xdippp,' XCN=',xcn,' difference due Mdot ANISO=',xdilex
       xltod=xltof
     endif   !   irot+isol

     jdiff=0
     idern=0
     if (.not. libgenec) then
     write(io_logs,'(/////,a,1p,e11.2,2(4x,a,e12.4)/46x,a,e12.4,42x/45x,4(1x,a,e12.4)/45x,4(1x,a,e12.4)/44x,4(1x,a,e12.4))') &
       ' CHANGEMENT DE LA CHIMIE    DZEIT=',dzeit,'x(m)=',x(m),'y(m)=',y(m),'y3(m)=',y3(m),'xc12(m)=',xc12(m),'xc13(m)=',xc13(m), &
       'xn14(m)=',xn14(m),'xn15(m)=',xn15(m),'xo16(m)=',xo16(m),'xo17(m)=',xo17(m),'xo18(m)=',xo18(m),'xne20(m)=',xne20(m), &
       'xne22(m)=',xne22(m),'xmg24(m)=',xmg24(m),'xmg25(m)=',xmg25(m),'xmg26(m)=',xmg26(m)

     if (ialflu == 1) then
       write(io_logs,&
               '(45x,4(a,e12.4)/45x,2(a,e12.4)/45x,4(a,e12.4)/45x,2(a,e12.4))') 'f19(m)=',xf19(m),'ne21(m)=',xne21(m),'na23(m)', &
               xna23(m),'al26g(m)=',xal26(m),'al27(m)=',xal27(m),'si28(m)=',xsi28(m),'neu(m)=',xneut(m),'prot(m)=',xprot(m), &
               ' c14(m)=',xc14(m),' f18(m) =',xf18(m),'bidon  =',xbid(m),'bidon1 =',xbid1(m)
     endif

     write(io_logs,'(10x,77("    (",i3,",",i3,") ",e11.3))')(nbzel(ii),nbael(ii),abelx(ii,m),ii=1,nbelx)
     endif

     dzeitj=dzeit/year

     do i=1,m
      zwi=exp(s(i))-1.d0
      s(i)=1.d0
      if (zwi /= 0.d0) then
        s(i)=hh6+log(abs(zwi))
      endif
      if (s(i) < 1.) then
        s(i)=1.d0
      endif
      s(i)=sign(s(i),zwi)
      zwi=exp(vs(i))-1.d0
      vs(i)=1.d0
      if (zwi /= 0.d0) then
        vs(i)=hh6+log(abs(zwi))
      endif
      if (vs(i) < 1.d0) then
        vs(i)=1.d0
      endif
      vs(i)=sign(vs(i),zwi)
     enddo

   endif   ! gkor

   call write4

!*******************************************************************************

   if (.not. libgenec) then
   write(io_logs,'(//a,1x,i7)') 'Result for model',nwmd
   endif

   if (idebug > 1) then
     write(*,*) 'call bordn'
   endif
   call bordn

! CORRECTIONS DE TEFF POUR LES ETOILES WR (CF. LANGER,1988)
   teffpr=0.d0
   if (is_WR > epsilon(is_WR)) then
!-----------------------------------------------------------------------
!   Modification D.Schaerer, juillet 1990:
!   Correction de TEFF pour vents stellaires tenant compte de diffusion par electrons libres (comme correction deja existante
!   ci-dessous) et des raies selon la theorie CAK (amelioree par Kudritzki, Pauldrac, Puls  A&A,219,205)
     raysl=10.d0**radius
     if (idebug > 1) then
       write(*,*) 'call corrwind'
     endif
     call corrwind(teffpr,teffel,xmdot,teff,raysl)
!-----------------Fin modification--------------------------------------
   endif

   xjspe1=0.d0
   xjspe2=0.d0
   if (irot /= 0) then
     if (idebug > 1) then
       write(*,*) 'call momevo'
     endif
     call momevo(r,vomegi,xltot,CorrOmega,.true.)

! pour calcul du moment specifique en xma1=3 et xma2=5
     if (gms >= 5.d0) then
       if (idebug > 1) then
         write(*,*) 'call momspe'
       endif
       call momspe(vomegi,xjspe1,xjspe2,gms)
     else if (gms >= 3.d0) then
       if (idebug > 1) then
         write(*,*) 'call momspe'
       endif
       call momspe(vomegi,xjspe1,xjspe2,gms)
       xjspe2=0.d0
     endif ! gms
   endif ! irot

   if (idebug > 1) then
     write(*,*) 'call enint'
   endif
   call enint

   if (.not. elemneg) then
! [Modif CG]
! Dans le cas "diffusion tout le temps (iadvec = 0), dlelex ne doit pas etre sauve pour le modele suivant.
! Dans le cas "diffusion-advection-diffusion-...", il doit etre sauve lors du passage "advection --> diffusion",
!                                                   et non sauvegarde lors du passage "diffusion --> advection".
     if (iadvec == 0 .or. mod(nwmd,2) == 1) then
       dlelexprev = 0.d0
     else
       dlelexprev = dlelex
     endif
     if (isol >= 1) then
       dlelexprev = 0.d0
     endif
! [/Modif]

     if (vxal26g(1)<1.d-75) then
       vxal26g(1)=0.d0
     endif
     if(snube7<1.d-75) then
       snube7 = 0.d0
     endif
     if(snub8<1.d-75) then
       snub8 = 0.d0
     endif
     do ii=iidraw,40
      drawcon(ii)=1.d0
     enddo
     
     call compute_k2_from_structure(k2_AMC)
     
     !vomegi(m-1) is printed as central rotation rate since vomegi(m) is not well computed when the core is radiative
     if (.not. libgenec) then
     write(io_buffer) &
       nwmd,alter,dzeitj,gms,gls,teff,teffpr,xmdot,rhoc,tc,jwint,(xzc(k),k=1,ixzc),qbc,qmnc,rapcri, &
       vomegi(1)+CorrOmega(1),vomegi(m-1),xobla,vequat,fMdot_rot,vcri1m,vcri2m,eddesm,vequam,rapomm,vcrit1,vcrit2,eddesc,rapom2, &
       dmneed,xmdotneed,dlelexsave,bmomit,btot,btotatm,xjspe1,xjspe2,ekrote,epote,ekine,erade, &
       vx(1),vy3(1),vy(1),vxc12(1),vxc13(1),vxc14(1),vxn14(1),vxn15(1),vxo16(1),vxo17(1),vxo18(1),vxf18(1),vxf19(1), &
       vxne20(1),vxne21(1),vxne22(1),vxna23(1),vxmg24(1),vxmg25(1),vxmg26(1),vxal26g(1),vxal27(1),vxsi28(1), &
       vx(m),vy3(m),vy(m),vxc12(m),vxc13(m),vxc14(m),vxn14(m),vxn15(m),vxo16(m),vxo17(m),vxo18(m),vxf18(m),vxf19(m), &
       vxne20(m),vxne21(m),vxne22(m),vxna23(m),vxmg24(m),vxmg25(m),vxmg26(m),vxal26g(m),vxal27(m),vxsi28(m), &
       vxneut(m),vxprot(m),vxbid(m),vxbid1(m),snube7,snub8,lcnom,xmcno,scno,(vabelx(ii,1),ii=1,nbelx),&
       (vabelx(ii,m),ii=1,nbelx),(drawcon(ii),ii=1,40),imloss,is_MS,is_OB,is_RSG,is_WR,k2_AMC
     endif

     xteffprev=xtefflast
     xtefflast=log10(teff)
     xlprev=xllast
     xllast=log10(gls)
     xrhoprev=xrholast
     xrholast=rhoc
     xcprev=xclast
     xtcprev=xtclast
     xtclast=tc

     select case (phase)
       case (1)
         xclast=vx(m)
       case (2,10)
         xclast=vy(m)
       case (3)
         xclast=vxc12(m)
       case (4)
         xclast=vxne20(m)
       case (5)
         xclast=vxo16(m)
       case (6)
         do ii=1,nbelx
          if (nbzel(ii) == 14 .and. nbael(ii) == 28) then
            xclast=vabelx(ii,m)   ! 28Si
          endif
         enddo
       case default
          rewind(io_runfile)
          write(io_runfile,*) nwmd,": Problem with the phase number"
          stop "Problem with the phase number ==> STOP"
     end select

! If pgplot is active, then call the needed routines.
     Species_PGplot(1) = vx(m)
     Species_PGplot(2) = vy(m)
     Species_PGplot(3) = vxc12(m)
     Species_PGplot(4) = vxn14(m)
     Species_PGplot(5) = vxo16(m)
     Species_PGplot(6) = vxne20(m)
     Species_PGplot(7) = vabelx(1,m)
     if (idebug > 1) then
       write(*,*) 'call SavePlotData'
     endif
     call SavePlotData(gms,gls,teff,nwmd,alter,tc,rhoc,Species_PGplot)
     if (mod(nwmd,10) == 0) then
       write(*,*) nwmd,mod(nwmd,10),': call TimestepControle'
       call TimestepControle
       xcn = xcnwant
       inum = 0
     else
       xcnwant = 1.d0
       xcn = 1.d0
     endif
     if (mod(nwmd,10) == 0) then
       if (iprezams == 2) then
         gkorm=0.10d0
         iprezams=0
       endif
     endif   ! nwmd % 10

! Computation of the ZAMS radius:
     if (x(m)<(x(1)-3.0d-3) .and. zams_radius <= 0.d0) then
       zams_radius = sqrt(gls*Lsol/(4.d0*pi*cst_sigma))/teff**2.d0
     endif

! Fin de la preZAMS automatique:
! Le programme boucle la serie et s'arrete
     if (abs(vwant) > 1.0d-5 .or. init_synchronized) then
       if (x(m)<(x(1)-3.0d-3)) then
         write(*,*) '***** End of preZAMS, usual changes of parameters *****'
         write(io_input_changes,*) nwmd,': ZAMS reached, usual changes of parameters'
         iprezams = 2
         vwant = 0.0d0
         init_synchronized = .false.
         xfom = 1.0d0
         islow = 0
         isol = 0
         if (istati == 1) then
           idiff=0
           iadvec=0
         else
           idiff = 1
         endif
         if (imagn == 0 .and. istati /=1 ) then
           iadvec = 1
           xdial = 1.0d0
           idialo = 1
           idialu = 1
         elseif (imagn > 0) then
           idialu = 1
         endif
         dgrp = 0.010d0*um
         dgrl = 0.010d0*um
         dgry = 0.0030d0
       endif ! x(m)<(x(1)-3.0d-3)

       if (iprezams==1 .and. (abs(vwant)>1.d-5 .or. (init_synchronized .and. bintide))) then
         if (idebug > 1) then
           write(*,*) 'calcul de xfom'
         endif
         ! In case of binaries only, possibility to initialize the spin angular velocity
         ! to the orbital angular velocity. In this case, the value of vwant is ignored
         ! (only need to give it a nonzero value otherwise a non-rotating star is created).
         if (init_synchronized .and. bintide) then
             write(*,*)'period',period
             xfom = min(2.d0*pi/(period*omegi(1)),1.2d0)
         else
           if (vwant > 1.0d0) then
             xfom = min(vwant/vequat,1.2d0)
           else if (vwant > 1.0d-5) then
             xfom =  min(vwant*vcrit1/vequat,1.2d0)
           else
             xfom = min(abs(vwant)/rapcri,1.2d0)
           endif
         endif
         write(*,*) 'xfom set to:',xfom
         write(io_input_changes,'(i6,a13,f9.5)') nwmd,': xfom set to',xfom
       endif ! iprezams==1
     endif ! abs(vwant) > 1.0d-5

     if (x(m)<(x(1)-3.0d-3) .and. prezams_winds_not_applied) then
       prezams_winds_not_applied = .false.
       winds_not_applied = .false.
     endif

     if (mod(nwmd,10)==0) then
       if (idebug > 1) then
         write(*,*) 'call INPUTS_Change'
       endif
       call INPUTS_Change(x(m),y(m),xc12(m),xne20(m),xo16(m),rapom2,m,nzmodini,nzmodnew)
     endif

     if (n_snap /= 0) then
       if (mod(nwmd,n_snap)==0) then
         if (idebug > 1) then
           write(*,*) 'call print_Snapshot, print_files, and switch_outputfile'
         endif
         call print_Snapshot
         snap_printed = .true.
         Struc_Plotted = .false.
         call print_files
         call switch_outputfile
       endif
     else ! n_snap == 0
       modanf = modanf + 1
     endif
!***********************************************************************
     if (&
             modell == nzmod&
             .or. phase==end_at_phase&
             .or. nwmd==end_at_model&
             ) then
       nzmodnew = nzmodini
       write(*,*) 'EXITING'
       exit   !   FIN DU BOUCLAGE DES MODELES, SERIE TERMINEE
     endif

!***********************************************************************
   endif ! ELEM NEG

   modell=modell+1
   nwmd=nwmd+1
   snap_printed = .false.
   write(*,*) 'Looping to new timestep, nwmd,modell:',nwmd, modell

! COUPURE QUAND LE MODELE FRAGMENTE LE PAS TEMPOREL INDEFINIMENT
   if (phase < 3) then
     if (dzeitj <= 1.0d-08) then
       rewind(io_runfile)
       write(io_runfile,*) nwmd,': time step too small'
       stop
     endif
   else
     if (dzeitj <= 1.0d-25) then
       rewind(io_runfile)
       write(io_runfile,*) nwmd,': time step too small'
       stop
     endif
   endif

   if (alter == 0.d0) then
     glsv=gls
     teffv=teff
     veryFirst=.false.
   endif

!******************* Fin boucle de calcul du modele ************************
  enddo
end subroutine evolve

subroutine finalise
  if (.not. snap_printed .and. n_snap /= 0) then
    call print_Snapshot
  endif
  if (idebug > 1) then
    write(*,*) 'call SequenceClosing'
  endif
  call SequenceClosing
end subroutine finalise

end module genec
