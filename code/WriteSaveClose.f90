module WriteSaveClose

use inputparam, only: writetofiles,readfromfiles
use evol,only: kindreal,ldi,npondcouche
use const,only: um
use inputparam,only: modanf,nwseq,nzmod,iprn,iauto,ialflu,ianiso,imagn,ipop3,irot,isol,idiff,iadvec,icoeff, &
  igamma,ibasnet,istati,iledou,idifcon,iover,iunder,my,ikappa,iopac,imloss,ifitm,itmin,nndr,idialo,idialu,phase,isugi,nbchx, &
  nrband,iout,icncst,islow,zinit,zsol,z,frein,dovhp,dunder,elph,fmlos,fitm,rapcrilim,omega,xfom,vwant,gkorm,alph,agdr, &
  agds,agdp,agdt,faktor,deltal,deltat,dgrp,dgrl,dgry,dgrc,dgro,dgr20,xdial,fenerg,richac,xcn,plot,refresh,starname, &
  Write_namelist,xyfiles,verbose
use caramodele,only: nwmd,glm,gms,gls,teff,glsv,teffv,ab,dm_lost,iwr,iprezams,xmini
use strucmod,only: m,q,p,t,r,s,vp,vt,vr,vs,drl,drte,drp,drt,drr,dk,rlp,rlt,rlc,rrp,rrt,rrc,rtp,rtt,rtc
use abundmod,only: x,y3,y,xc12,xc13,xc14,xn14,xn15,xo16,xo17,xo18,xf18,xf19,xne20,xne21,xne22,xna23,xmg24,xmg25,xmg26, &
                   xal26,xal27,xsi28,xprot,xneut,xbid,xbid1,ybe7,yb8,vx,vy3,vy,vxc12,vxc13,vxc14,vxn14,vxn15,vxo16,vxo17,vxo18, &
                   xal26,xal27,xsi28,xprot,xneut,xbid,xbid1,ybe7,yb8,vx,vy3,vy,vxc12,vxc13,vxc14,vxn14,vxn15,vxo16,vxo17,vxo18, &
                   vxf18,vxf19,vxne20,vxne21,vxne22,vxna23,vxmg24,vxmg25,vxmg26,vxal26g,vxal27,vxsi28,vxprot,vxneut,vxbid, &
                   vxbid1,nbael,nbzel,nbelx,abelx,vabelx,mbelx,lcnom,xmcno,scno
use rotmod,only: omegi,vomegi,CorrOmega
use convection,only: ixzc
use PGPlotModule,only: HRD_FileName,EndPGplot
use PrintAll,only: DataAll_FileName,File_Unit
use henyey_solver,only: nsugi
use diffadvmod,only: tdiff
use timestep,only: alter,dzeitj,dzeit,dzeitv

implicit none

logical,save:: quitafterclosing = .false.
integer,save:: nzmodini,nzmodnew,ichange,mold
real(kindreal),save:: xcprev,xclast,xteffprev,xtefflast,xrhoprev,xrholast,xageprev,xagelast,xtcprev,xtclast,xlprev,xllast
real(kindreal),save:: gmsold,alterold,glsold,teffold,glsvold,teffvold,dzeitjold,dzeitold,dzeitvold,abold,dmold,dkold,rlpold, &
                   rltold,rlcold,rrpold,rrtold,rrcold,rtpold,rttold,rtcold,tdiffold
real(kindreal),dimension(ldi),save:: qold,pold,told,rold,sold,xold,yold,xcold,vpold,vtold,vrold,vsold,xoold,vxold,vyold,vxc12old, &
                   vxo16old,y3old,xc13old,xn14old,xn15old,xo17old,xo18old,vy3old,vxc13old,vxn14old,vxn15old,vxo17old,vxo18old, &
                   x20old,x22old,x24old,x25old,x26old,vxne20old,vxne22old,vxmg24old,vxmg25old,vxmg26old,xneutold,vxneutold, &
                   xprotold,vxprotold,xc14old,vxc14old,xf18old,vxf18old,xbidold,vxbidold,xbid1old,vxbid1old,xf19old,xne21old, &
                   xal26gold,xal27old,xsi28old,xna23old,vxf19old,vxne21old,vxal26gold,vxal27old,vxsi28old,vxna23old,omold,vomold
real(kindreal),dimension(npondcouche),save::CorrOmegaOld
real(kindreal),dimension(3),save:: drlold,drteold,drpold,drtold,drrold
real(kindreal),dimension(mbelx,ldi),save:: abelxold,vabelxold

private
public:: OpenAll,SequenceClosing,CheckSchrit,CloseAll
public:: write4,read4
public:: nzmodnew,ichange,nzmodini
public:: xcprev,xclast,xteffprev
public:: quitafterclosing

contains
!=======================================================================
subroutine CheckSchrit(outString)
!-----------------------------------------------------------------------
implicit none

character(5),intent(in):: outString

integer:: ii,imb
real(kindreal):: xxmmbb,xxmmb1,xxmmbm,vxmmbb,vxmmb1,vxmmbm,xxmmhy,vxmmhy,xdm1
!-----------------------------------------------------------------------
! calcul du moment cinetique avant schrit
  xxmmbb=0.d0
  vxmmbb=0.d0
! pour controle conservation des abondances
  xxmmhy=0.d0
  vxmmhy=0.d0

  do imb=2,m-1
   xdm1=(exp(q(imb+1))-exp(q(imb-1)))/2.d0*exp(glm)
   xxmmb1=omegi(imb)*exp(2.d0*r(imb))*xdm1*2.d0/3.d0
   xxmmbb=xxmmbb+xxmmb1
   vxmmb1=vomegi(imb)*exp(2.d0*r(imb))*xdm1*2.d0/3.d0
   vxmmbb=vxmmbb+vxmmb1
   xxmmhy=xxmmhy+xo16(imb)*xdm1
   vxmmhy=vxmmhy+vxo16(imb)*xdm1
  enddo

  xxmmb1=2.d0/3.d0*exp(2.d0*r(1))*omegi(1)*(exp(q(2))-exp(q(1)))/2.d0*exp(glm)
  xxmmbm=2.d0/5.d0*exp(glm)*(1.d0-exp(q(m-1)))/2.d0*exp(r(m-1))/(2.d0)**(1.d0/3.d0)*exp(r(m-1))/(2.d0)**(1.d0/3.d0)*omegi(m)
  xxmmbb=xxmmbb+xxmmb1+xxmmbm

  vxmmb1=2.d0/3.d0*exp(2.d0*r(1))*vomegi(1)*(exp(q(2))-exp(q(1)))/2.d0*exp(glm)
  vxmmbm=2.d0/5.d0*exp(glm)*(1.d0-exp(q(m-1)))/2.d0*exp(r(m-1))/(2.d0)**(1.d0/3.d0)*exp(r(m-1))/(2.d0)**(1.d0/3.d0)*vomegi(m)
  vxmmbb=vxmmbb+vxmmb1+vxmmbm

  xxmmhy=xxmmhy+xo16(1)*(exp(q(2))-exp(q(1)))/2.d0*exp(glm)+xo16(m)*exp(glm)*(1.d0-exp(q(m-1)))/2.d0
  vxmmhy=vxmmhy+vxo16(1)*(exp(q(2))-exp(q(1)))/2.d0*exp(glm)+vxo16(m)*exp(glm)*(1.d0-exp(q(m-1)))/2.d0

  if (verbose) then
    write(*,'(1x,a5,1x,2(a,es24.17),a,i5)') outString,'SCHRIT, Ib=',xxmmbb,' VIb=',vxmmbb,' m=',m
    write(*,'(1x,a5,1x,2(a,es24.17),a,i5)') outString,'SCHRIT, O16=',xxmmhy,' O16=',vxmmhy,' m=',m
  endif

  do ii=1,nbelx
   xxmmhy=0.d0
   vxmmhy=0.d0
   do imb=2,m-1
    xdm1=(exp(q(imb+1))-exp(q(imb-1)))/2.d0*exp(glm)
    xxmmhy=xxmmhy+abelx(ii,imb)*xdm1
    vxmmhy=vxmmhy+vabelx(ii,imb)*xdm1
   enddo
   xxmmhy=xxmmhy+abelx(ii,1)*(exp(q(2))-exp(q(1)))/2.d0*exp(glm)+abelx(ii,m)*exp(glm)*(1.d0-exp(q(m-1)))/2.d0
   vxmmhy=vxmmhy+vabelx(ii,1)*(exp(q(2))-exp(q(1)))/2.d0*exp(glm)+vabelx(ii,m)*exp(glm)*(1.d0-exp(q(m-1)))/2.d0
   if (verbose) then
     write(*,'(1x,a5,1x,a,i3,2(a,es24.17))') outString,'SCHRIT,i= ',ii,' ab=',xxmmhy,' vab=',vxmmhy
   endif
  enddo

  return

end subroutine CheckSchrit
!=======================================================================
subroutine write4
!---------------------------------------------------------------------
  implicit none

  gmsold   =    gms
  alterold =    alter
  glsold   =    gls
  teffold  =    teff
  glsvold  =    glsv
  teffvold =    teffv
  dzeitjold=    dzeitj
  dzeitold =    dzeit
  dzeitvold=    dzeitv
  abold    =    ab
  dmold    =    dm_lost
  mold     =    m

  qold(1:m)    =     q(1:m)
  pold(1:m)    =     p(1:m)
  told(1:m)    =     t(1:m)
  rold(1:m)    =     r(1:m)
  sold(1:m)    =     s(1:m)
  xold(1:m)    =     x(1:m)
  yold(1:m)    =     y(1:m)
  xcold(1:m)   =     xc12(1:m)
  vpold(1:m)   =     vp(1:m)
  vtold(1:m)   =     vt(1:m)
  vrold(1:m)   =     vr(1:m)
  vsold(1:m)   =     vs(1:m)
  xoold(1:m)   =     xo16(1:m)
  vxold(1:m)   =     vx(1:m)
  vyold(1:m)   =     vy(1:m)
  vxc12old(1:m)=     vxc12(1:m)
  vxo16old(1:m)=     vxo16(1:m)
  y3old(1:m)   =     y3(1:m)
  xc13old(1:m) =     xc13(1:m)
  xn14old(1:m) =     xn14(1:m)
  xn15old(1:m) =     xn15(1:m)
  xo17old(1:m) =     xo17(1:m)
  xo18old(1:m) =     xo18(1:m)
  vy3old(1:m)  =     vy3(1:m)
  vxc13old(1:m)=     vxc13(1:m)
  vxn14old(1:m)=     vxn14(1:m)
  vxn15old(1:m)=     vxn15(1:m)
  vxo17old(1:m)=     vxo17(1:m)
  vxo18old(1:m)=     vxo18(1:m)
  x20old(1:m)  =     xne20(1:m)
  x22old(1:m)  =     xne22(1:m)
  x24old(1:m)  =     xmg24(1:m)
  x25old(1:m)  =     xmg25(1:m)
  x26old(1:m)  =     xmg26(1:m)
  vxne20old(1:m)=    vxne20(1:m)
  vxne22old(1:m)=    vxne22(1:m)
  vxmg24old(1:m)=    vxmg24(1:m)
  vxmg25old(1:m)=    vxmg25(1:m)
  vxmg26old(1:m)=    vxmg26(1:m)

  xneutold(1:m) =     xneut(1:m)
  vxneutold(1:m)=     vxneut(1:m)
  xprotold(1:m) =     xprot(1:m)
  vxprotold(1:m)=     vxprot(1:m)
  xc14old(1:m)  =     xc14(1:m)
  vxc14old(1:m) =     vxc14(1:m)
  xf18old(1:m)  =     xf18(1:m)
  vxf18old(1:m) =     vxf18(1:m)
  xbidold(1:m)  =     xbid(1:m)
  vxbidold(1:m) =     vxbid(1:m)
  xbid1old(1:m) =     xbid1(1:m)
  vxbid1old(1:m)=     vxbid1(1:m)
  xf19old(1:m)  =     xf19(1:m)
  xne21old(1:m) =     xne21(1:m)
  xal26gold(1:m)=     xal26(1:m)
  xal27old(1:m) =     xal27(1:m)
  xsi28old(1:m) =     xsi28(1:m)
  xna23old(1:m) =     xna23(1:m)
  vxf19old(1:m) =     vxf19(1:m)
  vxne21old(1:m)=     vxne21(1:m)
  vxal26gold(1:m)=    vxal26g(1:m)
  vxal27old(1:m)=     vxal27(1:m)
  vxsi28old(1:m)=     vxsi28(1:m)
  vxna23old(1:m)=     vxna23(1:m)

  omold(1:m)   =     omegi(1:m)
  vomold(1:m)  =     vomegi(1:m)

  abelxold(1:nbelx,1:m)=abelx(1:nbelx,1:m)
  vabelxold(1:nbelx,1:m)=vabelx(1:nbelx,1:m)

  CorrOmegaOld(1:npondcouche) = CorrOmega(1:npondcouche)

  drlold(1:3)   =    drl(1:3)
  drteold(1:3)  =    drte(1:3)
  drpold(1:3)   =    drp(1:3)
  drtold(1:3)   =    drt(1:3)
  drrold(1:3)   =    drr(1:3)

  dkold    =    dk
  rlpold   =    rlp
  rltold   =    rlt
  rlcold   =    rlc
  rrpold   =    rrp
  rrtold   =    rrt
  rrcold   =    rrc
  rtpold   =    rtp
  rttold   =    rtt
  rtcold   =    rtc
  tdiffold =    tdiff

  return

end subroutine write4
!=======================================================================
subroutine read4
!---------------------------------------------------------------------
  implicit none

  gms      =  gmsold
  alter    =  alterold
  gls      =  glsold
  teff     =  teffold
  glsv     =  glsvold
  teffv    =  teffvold
  dzeitj   =  dzeitjold
  dzeit    =  dzeitold
  dzeitv   =  dzeitvold
  ab       =  abold
  dm_lost  =  dmold
  m        =  mold

  q(1:m)     =   qold(1:m)
  p(1:m)     =   pold(1:m)
  t(1:m)     =   told(1:m)
  r(1:m)     =   rold(1:m)
  s(1:m)     =   sold(1:m)
  x(1:m)     =   xold(1:m)
  y(1:m)     =   yold(1:m)
  xc12(1:m)    =   xcold(1:m)
  vp(1:m)    =   vpold(1:m)
  vt(1:m)    =   vtold(1:m)
  vr(1:m)    =   vrold(1:m)
  vs(1:m)    =   vsold(1:m)
  xo16(1:m)    =   xoold(1:m)
  vx(1:m)    =   vxold(1:m)
  vy(1:m)    =   vyold(1:m)
  vxc12(1:m) =   vxc12old(1:m)
  vxo16(1:m) =   vxo16old(1:m)
  y3(1:m)    =   y3old(1:m)
  xc13(1:m)  =   xc13old(1:m)
  xn14(1:m)  =   xn14old(1:m)
  xn15(1:m)  =   xn15old(1:m)
  xo17(1:m)  =   xo17old(1:m)
  xo18(1:m)  =   xo18old(1:m)
  vy3(1:m)   =   vy3old(1:m)
  vxc13(1:m) =   vxc13old(1:m)
  vxn14(1:m) =   vxn14old(1:m)
  vxn15(1:m) =   vxn15old(1:m)
  vxo17(1:m) =   vxo17old(1:m)
  vxo18(1:m) =   vxo18old(1:m)
  xne20(1:m)   =   x20old(1:m)
  xne22(1:m)   =   x22old(1:m)
  xmg24(1:m)   =   x24old(1:m)
  xmg25(1:m)   =   x25old(1:m)
  xmg26(1:m)   =   x26old(1:m)
  vxne20(1:m)=   vxne20old(1:m)
  vxne22(1:m)=   vxne22old(1:m)
  vxmg24(1:m)=   vxmg24old(1:m)
  vxmg25(1:m)=   vxmg25old(1:m)
  vxmg26(1:m)=   vxmg26old(1:m)

  xneut(1:m) =   xneutold(1:m)
  vxneut(1:m)=   vxneutold(1:m)
  xprot(1:m) =   xprotold(1:m)
  vxprot(1:m)=   vxprotold(1:m)
  xc14(1:m)  =   xc14old(1:m)
  vxc14(1:m) =   vxc14old(1:m)
  xf18(1:m)  =   xf18old(1:m)
  vxf18(1:m) =   vxf18old(1:m)
  xbid(1:m)  =   xbidold(1:m)
  vxbid(1:m) =   vxbidold(1:m)
  xbid1(1:m) =   xbid1old(1:m)
  vxbid1(1:m)=   vxbid1old(1:m)
  xf19(1:m)  =   xf19old(1:m)
  xne21(1:m) =   xne21old(1:m)
  xal26(1:m)=   xal26gold(1:m)
  xal27(1:m) =   xal27old(1:m)
  xsi28(1:m) =   xsi28old(1:m)
  xna23(1:m) =   xna23old(1:m)
  vxf19(1:m) =   vxf19old(1:m)
  vxne21(1:m)=   vxne21old(1:m)
  vxal26g(1:m)=  vxal26gold(1:m)
  vxal27(1:m)=   vxal27old(1:m)
  vxsi28(1:m)=   vxsi28old(1:m)
  vxna23(1:m)=   vxna23old(1:m)

  omegi(1:m) =   omold(1:m)
  vomegi(1:m)=   vomold(1:m)

  abelx(1:nbelx,1:m)=abelxold(1:nbelx,1:m)
  vabelx(1:nbelx,1:m)=vabelxold(1:nbelx,1:m)

  CorrOmega(1:npondcouche) = CorrOmegaOld(1:npondcouche)

  drl(1:3)     =   drlold(1:3)
  drte(1:3)    =   drteold(1:3)
  drp(1:3)     =   drpold(1:3)
  drt(1:3)     =   drtold(1:3)
  drr(1:3)     =   drrold(1:3)

  dk         =   dkold
  rlp        =   rlpold
  rlt        =   rltold
  rlc        =   rlcold
  rrp        =   rrpold
  rrt        =   rrtold
  rrc        =   rrcold
  rtp        =   rtpold
  rtt        =   rttold
  rtc        =   rtcold
  tdiff      =   tdiffold


  return

end subroutine read4
!=======================================================================
subroutine SequenceClosing
!-----------------------------------------------------------------------
use const,only: cstlg_K1,cstlg_mh,cstlg_k
use inputparam,only: INPUTS_Change,stop_deg
use timestep,only: TimestepControle

implicit none

integer:: error9,nm,ii,k,kim,lcnom,jwint
real(kindreal):: age9,mass9,ll9,teff9,x1,ne201,y1,c121,c131,n141,ne221,o161,o171,o181,xmdot,rhoc,tc,xm,ne20m,ym,c12m,c13m,n14m, &
  ne22m,o16m,o17m,o18m,qbc,qmnc,teffpr,rapcri,rot1,rotm,xobla,vequat,alpro6,xmcno,scno,dzeitj,vcri1m,vcri2m,eddesm,vequam,rapomm, &
  vcrit1,vcrit2,eddesc,rapom2,dmneed,xmdotneed,dlelex,bmomit,btot,ekrote,epote,ekine,erade,xjspe1,xjspe2,f191,ne211,al261, &
  al271,si281,na231,f19m,ne21m,al26m,al27m,si28m,na23m,y31,n151,mg241,mg251,mg261,y3m,n15m,mg24m,mg25m,mg26m,neutm,protm,c14m, &
  f18m,bidm,bid1m,btotatm,snube7,snub8,fluxbe7,fluxb8
real(kindreal):: PrintVelocity,xl,xte,xtt,xcnwant,tcdeg
real(kindreal),dimension(ldi):: abel9
real(kindreal),dimension(40):: drawc
real(kindreal),dimension(ixzc):: xzc
!-----------------------------------------------------------------------
  rewind 9

! If pgplot is active, terminate the module
  if (plot .and. quitafterclosing) then
    call EndPGplot
  endif

  if (writetofiles) write(3,*) 'm, isugi, nsugi:',m,isugi,nsugi

  if (writetofiles) write(10,'(/2x,"NB",6x,"AGE",8x,"MASS",3x,"LOGL",2x,"LOGTE",5x,"X",8x,"Y",7x,"C12",6x,"C13",6x,"N14",6x,"O16",&
          &6x,"O17",6x,"O18",&
    &5x,"NE20",5x,"NE22"/10x,"QCC",8x,"MDOT",3x,"RHOC",2x,"LOGTC"/10x," O ",8x," Ve ",3x," Fc "/)')

  error9 = 0
  do while (error9 == 0)
   read(9,iostat=error9) nm,age9,dzeitj,mass9,ll9,teff9,teffpr,xmdot,rhoc,tc,jwint,(xzc(k),k=1,ixzc),qbc,qmnc,rapcri,rot1,rotm, &
     xobla,vequat,alpro6,vcri1m,vcri2m,eddesm,vequam,rapomm,vcrit1,vcrit2,eddesc,rapom2,dmneed,xmdotneed,dlelex,bmomit,btot, &
     btotatm,xjspe1,xjspe2,ekrote,epote,ekine,erade,x1,y31,y1,c121,c131,n141,n151,o161,o171,o181,ne201,ne221,mg241,mg251,mg261, &
     xm,y3m,ym,c12m,c13m,n14m,n15m,o16m,o17m,o18m,ne20m,ne22m,mg24m,mg25m,mg26m,f191,ne211,na231,al261,al271,si281,f19m,ne21m, &
     na23m,al26m,al27m,si28m,neutm,protm,c14m,f18m,bidm,bid1m,snube7,snub8,lcnom,xmcno,scno

   if (error9 == 0) then
     read(9) (abel9(ii),ii=1,2*nbelx)

     read(9) (drawc(ii),ii=1,40)

     if (irot == 1) then
       if (vcrit2 /= 0.d0) then
         PrintVelocity = vequat/min(vcrit1,vcrit2)
       else
         PrintVelocity = vequat/vcrit1
       endif
     else
       PrintVelocity = 0.d0
     endif
     xl=log10(ll9)
     xte=log10(teff9)
     xtt=xte
     if (x(1) < 0.30d0 .and. xte > 4.0d0 .and. teffpr /= 0.d0) then
       xtt=teffpr
     endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ECRITURE DU .S :
     if (writetofiles) write(10,'(i6,1pe14.7,0pf9.4,2(1x,f6.3),1x,f9.6,1x,f9.6,8(1x,1pe8.2)/5x,0pf7.4,1x,f6.3,2x,f7.3,2(1x,f6.3),&
             &1x,f9.6,1x,f9.6,&
       &8(1x,1pe8.2)/5x,0pf7.4,3x,1pe10.4,1x,0pf7.4,1x,0pf13.10,3x,i4,1x,f9.4,1x,f10.7/,5x,a,1x,e10.4,/,a,f8.2,1x,a,f8.2,1x,a,&
       &f8.2,1x,a,f8.2,1x,a,f9.6,/,a,f8.2,1x,a,f8.2,1x,a,f8.2,1x,a,f8.2,1x,a,f9.6,1x,a,f9.6,/1x,a,f10.3,1x,a,f10.3)') nm,age9, &
       mass9,xl,xtt,x1,y1,c121,c131,n141,o161,o171,o181,ne201,ne221,qmnc,xte,xmdot,rhoc,tc,xm,ym,c12m,c13m,n14m,o16m,o17m,o18m, &
       ne20m,ne22m,xobla,vequat,&
       rapcri,rot1,lcnom,xmcno,scno,'DELTA t=',dzeitj,'valeurs pour calcul Mdot: vcrit1=',vcri1m,'vcrit2=',vcri2m,'vequat=',&
       vequam,'omega/omegacrit=',rapomm,'EDDING. FAC=',eddesm,'valeurs bon modele      : vcrit1=',vcrit1,'vcrit2=',vcrit2,&
       'vequat=',vequat,'omega/omegacrit=',rapom2,'EDDING. FAC=',eddesc,'veq/vcrit=',PrintVelocity,'mom spe a 3Msol=',xjspe1,&
       'mom spe a 5Msol=',xjspe2

     if (ialflu == 1) then
       if (writetofiles) write(10,'(1x,6(a,e12.4)/1x,6(a,e12.4)/1x,6(a,e12.4))') 'f19(1)=',f191,'ne21(1)=',ne211,'na23(1)=',&
               na231,'al26g(1)=',al261,&
         'al27(1)=',al271,'si28(1)=',si281,'f19(m)=',f19m,'ne21(m)=',ne21m,'na23(m)=',na23m,'al26g(m)=',al26m,'al27(m)=',al27m, &
         'si28(m)=',si28m,&
         'neu(m)=',neutm,'pro(m)=',protm,'xc14(m)=',c14m,'xf18(m)=',f18m,'bidon(m)=',bidm,'bidon1=',bid1m
     endif

     if (writetofiles) write(10,'(77(1x,"(",i3,",",i3,")(1)= ",e11.4))') (nbzel(ii),nbael(ii),abel9(ii),ii=1,nbelx)
     if (writetofiles) write(10,'(77(1x,"(",i3,",",i3,")(m)= ",e11.4))') (nbzel(ii-nbelx),nbael(ii-nbelx),abel9(ii),&
             ii=nbelx+1,2*nbelx)
     if (writetofiles) write(10,*)
     if (snube7>=1.0d-60) then
       fluxbe7 = snube7/2.38d-10
     else
       fluxbe7 = 0.d0
       snube7 = 0.d0
     endif
     if (snub8>=1.0d-56) then
       fluxb8 = snub8/1.08d-06
     else
       fluxb8 = 0.d0
       snub8 = 0.d0
     endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! IMPRESSION COMPLEMENTAIRE POUR GRAPHIQUE (.G) :
     if (writetofiles) write(20,'(i6,1x,1pe22.15,0pf11.6,2(1x,f9.6),2(1x,e14.7),1p,9(1x,e14.7),1x,0pf7.4,3x,f9.6,1x,f7.3,&
             &2(1x,f9.6),2(1x,e14.7),1p,&
       &9(1x,e14.7),2(1x,e10.3),2(1x,e10.3),2(1x,e10.3),0pf12.8,6(1x,1pe10.3),1x,i4,1x,0pf9.4,1x,1pe9.2,2(1x,e10.4),0p,3x,&
       &3(1x,1pe8.2),0p,2(1x,f9.6),3(1x,1pe8.2),0p,2(1x,f9.6),9(1x,1pe14.7),0p,40f6.3,1x,1pe17.10)') nm,age9,mass9,xl,xtt,x1,y1, &
       y31,c121,c131,n141,o161,o171,o181,ne201,ne221,qmnc,xte,xmdot,rhoc,tc,xm,ym,y3m,c12m,c13m,n14m,o16m,o17m,o18m,ne20m,ne22m, &
       ybe7(m)*7.d0,yb8(m)*8.d0,fluxbe7,fluxb8,snube7,snub8,rapcri,rot1,rotm,xobla,al261,al26m,alpro6,lcnom,xmcno, &
       scno,xjspe1,xjspe2,vcri1m,vcri2m,vequam,rapomm,eddesm,vcrit1,vcrit2,vequat,rapom2,eddesc,dmneed,xmdotneed,dlelex/1.d53, &
       bmomit/1.d57,btot/1.d53,ekrote/1.d51,epote/1.d51,ekine/1.d51,erade/1.d51,(drawc(ii),ii=1,40),btotatm/1.d53

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! FICHIERS DES ABONDANCES (.A) :
     if (writetofiles) write(23,'(1x,i6,1x,1pe20.13,0pf9.4,64(1x,e12.6))') nm,age9,mass9,x1,y31,y1,c121,c131,n141,n151,o161,o171,&
             o181,ne201,ne221, &
       mg241,mg251,mg261,f191,ne211,na231,al261,al271,si281,(abel9(ii),ii=1,nbelx),xm,y3m,ym,c12m,c13m,n14m,n15m,o16m,o17m,o18m, &
       ne20m,ne22m,mg24m,mg25m,mg26m,f19m,ne21m,na23m,al26m,al27m,si28m,(abel9(ii),ii=nbelx+1,2*nbelx)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     xcprev=xclast

     select case (phase)
       case (1)
         xclast=xm
       case (2,10)
         xclast=ym
       case (3)
         xclast=c12m
       case (4)
         xclast=ne20m
       case (5)
         xclast=o16m
       case (6)
         do ii=nbelx+1,2*nbelx
          if (nbzel(ii-nbelx) == 14 .and. nbael(ii-nbelx) == 28) then
            xclast=abel9(ii)   ! 28Si
          endif
         enddo
       case default
         if (writetofiles) rewind(222)
         if (writetofiles) write(222,*) nwmd,": Problem with the phase number"
         stop "Problem with the phase number ==> STOP"
     end select

     xteffprev=xtefflast
     xtefflast=xte
     xrhoprev=xrholast
     xrholast=rhoc
     xageprev=xagelast
     xagelast=age9
     xtcprev=xtclast
     xtclast=tc
     xlprev=xllast
     xllast=xl

     if (iwr == 1)then
       if (writetofiles) write(10,'(10x,"LOG TEFF NON MODIFIEE  =", f6.3)') xte
     endif
     if (jwint == 0) then
       write(*,*) '  * ENTIEREMENT RADIATIVE'
     else
       do k=1,jwint
        kim=2*k-1
        if (jwint /= 1 .or. xzc(1) /= 10000.d0) then
          if (xzc(1) == 10000.d0) xzc(1)=0.d0
          if (writetofiles) write(10,'(3x,a,i3,1x,2(1x,a,f8.4))') 'ZONE',k,'MR/M INF=',xzc(kim),'SUP=',xzc(kim+1)
        endif
       enddo
     endif
   endif
  enddo   ! error9

  call INPUTS_Change(xm,ym,c12m,ne20m,O16m,rapom2,m,nzmodini,nzmodnew)

  call TimestepControle(xcprev,xclast,xteffprev,xtefflast,xlprev,xllast,xrhoprev,xrholast,nzmodini,xcnwant)

  write(*,'(25x,a,f14.10)') 'SURFACE H ABUNDANCE: ',x(1)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! ECRITURE DU .INPUT :
  if (writetofiles) call Write_namelist(31,nwseq+nzmodini,modanf+1,nzmodnew,xcnwant)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! [Modif CG]
! Arret de l'execution si le pas de temps devient trop petit sur la MS.
  if (ianiso /= 0 .and. phase == 1 .and. dzeitj <= 20.d0 .and. verbose) then
    write(*,*) '!*!*!*!*!*!*!*!*!'
    write(*,*) 'Time step less than 20 years, too short.'
    write(*,*) '*!*!*!*!*!*!*!*!*'
  endif
! [/Modif]

  if (phase == 3) then
! pour masses <= 7  M et > =1.5 M
    if (xmini < 7.d0 .and. xmini > 1.7d0 .and. stop_deg) then
! dans l'equation de Tdeg, on a pose mu=mu_e=2:
      tcdeg=(2.d0/3.d0)*xrholast+cstlg_K1+cstlg_mh-cstlg_k-(2.d0/3.d0)*log10(2.d0)
      if (xtclast < tcdeg) then
        write(*,*) 'Central T lower than Tdeg ==> STOP'
        if (writetofiles) rewind(222)
        if (writetofiles) write(222,*) nwmd,': Central T lower than Tdeg ==> STOP'
        call CloseAll
        stop 'Central T lower than Tdeg ==> STOP'
      endif
    else if (xmini < 1.7d0 .and. stop_deg) then
      if (xtclast >= 7.9d0) then
        if (writetofiles) rewind(222)
        if (writetofiles) write (222,*) nwmd,': Central T greater than 7.9 ==> STOP'
        call CloseAll
        stop 'Central T greater than 7.9 ==> STOP'
      endif
    endif
  endif
! file runfile written to continue calculation
  if (irot == 0 .or. iprezams /= 2) then
    if (writetofiles) rewind(222)
    if (writetofiles) write (222,*) 'running'
  endif

  call CloseAll

  if (nzmodini > 1) then
    write(*,*) 'Sequence ',nwseq,'-',nwseq+nzmodini-1
    if (quitafterclosing) then
        stop 'Sequence successfully computed ! '
    endif
  else
    write(*,*) 'Model ',nwseq
    if (quitafterclosing) then
        stop 'Model successfully computed ! '
    endif
  endif

  if (.not. quitafterclosing) then
      write(*,*) '** Old values: nwseq=',nwseq,' modanf=',modanf,' nzmod=',nzmod,' xcn=',xcn, ' nzmodini=',nzmodini
      nwseq = nwseq + nzmodini
      modanf = modanf + 1
      nzmod = nzmodnew
      nwmd = nwmd + 1
      xcn = xcnwant
      write(*,*) '** New values: nwseq=',nwseq,' modanf=',modanf,' nzmod=',nzmod,' xcn=',xcn,' nwmd=',nwmd
  endif
  return

end subroutine SequenceClosing
!=======================================================================
subroutine OpenAll
!-----------------------------------------------------------------------
use inputparam,only: const_per
implicit none

logical:: fexists=.true.
character(5):: fnamein,fnameout
character(7):: ffmodel
character(256):: fname3,fname10,fname20,fname23,fname29,fname31,fname39,fname51,fname52,fname997,fname81,fname998,fname999
!-----------------------------------------------------------------------
  write(ffmodel,'(i7.7)') nwseq
  write(fnamein,'(i5.5)') modanf
  write(fnameout,'(i5.5)') modanf+1

  fname3  =  trim(starname)//'.l'//ffmodel
  fname10 =  trim(starname)//'.s'//ffmodel
  fname20 =  trim(starname)//'.g'//ffmodel
  fname23 =  trim(starname)//'.a'//ffmodel
  fname29 =  trim(starname)//'.v'//ffmodel
  fname31 =  trim(starname)//'.input'
  fname39 =  trim(starname)//'.z'//ffmodel
  fname51 = trim(starname)//'.b'//fnamein
  fname52 = trim(starname)//'.b'//fnameout
  fname997 =  'input_changes.log'
  if (.not. const_per) then
    fname81 = trim(starname)//'.period_evol.dat'
  endif
  if (xyfiles) then
    fname998 =  trim(starname)//'.x'//ffmodel
    fname999 =  trim(starname)//'.y'//ffmodel
  endif
  HRD_FileName = ".PlotData_"//trim(starname)
  DataAll_FileName = trim(starname)//"_StrucData_"//ffmodel//".dat"

  if (readfromfiles) then
  if (fnamein /= '00000') then
    inquire(file=fname51,exist=fexists)
    if (.not. fexists) then
      inquire(file=fname51//'.gz',exist=fexists)
      if (fexists) then
        write(*,*) 'bfile ',fnamein,' is zipped, unzip before launching GENEC!'
        stop
      else
        write(*,*) 'bfile ',fnamein,' does not exist, check the input parameters!'
        stop
      endif
    endif
  endif
  endif !readfromfiles

  open(9, status='scratch',form='unformatted')
  open(997,file=fname997,status='unknown',form='formatted',access='append')
  if (readfromfiles) then
  open(51,file=fname51,status='unknown',form='unformatted')
  endif

  if (writetofiles) then
  open(3,file=fname3, status='unknown',form='formatted')
  open(10,file=fname10,status='unknown',form='formatted')
  open(20,file=fname20,status='unknown',form='formatted')
  open(23,file=fname23,status='unknown',form='formatted')
  open(29,file=fname29,status='unknown',form='formatted')
  open(31,file=fname31,status='unknown',form='formatted')
  open(39,file=fname39,status='unknown',form='formatted')
  open(52,file=fname52,status='unknown',form='unformatted')
  if (.not. const_per) then
    open(81,file=fname81,status='unknown',form='formatted',access='append')
  endif
  if (xyfiles) then
    open(998,file=fname998,status='unknown',form='formatted')
    open(999,file=fname999,status='unknown',form='formatted')
  endif
  open (222,file='runfile',status='unknown',form='formatted')
  open(unit=File_Unit,file=DataAll_FileName,status="unknown")
  endif !writetofiles

  xcprev=0.d0
  xclast=0.d0
  xteffprev=0.d0
  xtefflast=0.d0
  xllast=0.d0
  xrhoprev=0.d0
  xrholast=0.d0
  xtclast=0.d0
  xageprev=0.d0
  xagelast=0.d0

  return

end subroutine OpenAll
!=======================================================================
subroutine CloseAll
!-----------------------------------------------------------------------
use inputparam,only: const_per
implicit none
!-----------------------------------------------------------------------
  close(9)
  close(997)
  if (writetofiles) then
  close(222)
  close(3)
  close(10)
  close(20)
  close(23)
  close(29)
  close(31)
  close(39)
  close(51)
  close(52)
  if (.not. const_per) then
    close(81)
  endif
  close(File_Unit)

  endif !writetofiles
  return

end subroutine CloseAll
!=======================================================================
end module WriteSaveClose
