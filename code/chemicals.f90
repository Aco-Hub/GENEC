module chemicals

use io_definitions
use evol,only: ldi,kindreal
use inputparam,only: phase,irot,isol,idiff,idifcon,ialflu,nbchx,idern,nrband,ichem,ipop3,verbose,idebug
use caramodele,only: nwmd
use abundmod,only: x,y3,y,xc12,xc13,xc14,xn14,xn15,xo16,xo17,xo18,xf18,xf19,xne20,xne21,xne22,xna23,xmg24,xmg25,xmg26,xal26, &
  xal27,xsi28,xprot,xneut,xbid,xbid1,vx,vy3,vy,vxc12,vxc13,vxc14,vxn14,vxn15,vxo16,vxo17,vxo18,vxf18,vxf19,vxne20,vxne21,vxne22, &
  vxna23,vxmg24,vxmg25,vxmg26,vxal26g,vxal27,vxsi28,vxbid,vxbid1,vxneut,vxprot,vvx,vvy3,vvy,vvxc12,vvxc13,vvxc14,vvxn14,vvxn15, &
  vvxo16,vvxo17,vvxo18,vvxf18,vvxf19,vvxne20,vvxne21,vvxne22,vvxna23,vvxmg24,vvxmg25,vvxmg26,vvxal26g,vvxal27,vvxsi28,vvxprot, &
  vvxneut,vvxbid,vvxbid1,epsy,epsyy,epsyc,epsyo,epsc,b11,b33,b34,b112,b113,b114,b115a,b115g,b116,b117a,b117g,b118a,b118g,b119g, &
  b119a,b120,b121,b122,b123g,b123a,b124,b125g,b125m,b1mg26,b1al26,b127g,b127a,epcne,eps20,c144,c184,c224,c134,c12ago,o16agn, &
  epcna,e17an,e20ag,c224g,e12ng,e14np,e14ng,e14be,e15ag,e17ag,e18an,e18ng,e18be,e19ng,e21ag,e21ng,e21na,e22ng,e23ng,e24ag,e24ng, &
  e25an,e25ng,e26ng,e26be,e27ng,e28ng,a26ga,a26gp,ec14pg,ec14ag,ec14ng,ef18na,ef18np,e18pa,e19ap,e20ng,mbelx,nbelx,abelx,vabelx, &
  vvabelx,abelx,zabelx,fnucdif,eps
use equadiffmod,only: iter
use strucmod,only: m,q,zensi,t,rho
use timestep,only: alter,dzeit
use convection,only: jzint,ixzc,izc
use energy,only: netburning

implicit none

integer:: i,ii,iii,k,iint,nt,n,i1,i2

real(kindreal),parameter,private:: ey=5.804d+17,ec=1.728d+18,eo=1.142d+18,ecne=1.115d+18,ecna=5.408d+17,e23=9.980d+16, &
    e20=2.246d+18,e24=2.39d+18,w1628=2.296d+18,w1631=1.837d+18,w3128=4.588d+17,w28=3.734d+17
real(kindreal),private:: sumx,sumy,sumy3,sumxc12,sumc13,sumn14,sumn15,sumxo16,sumo17,sumo18,sumne20,sumne22,summg24,summg25, &
    summg26,sumdm,sumc14,sumf18,sumf19,sumne21,sumna23,sumal26g,sumal27,sumsi28,sumneut,sumprot,sumbid,sumbid1,dm,dmx,dmy3,dmy, &
    dmxc12,dmc13,dmn14,dmn15,dmxo16,dmo17,dmo18,dmxne20,dmxne22,dmxmg24,dmxmg25,dmxmg26,dmc14,dmf18,dmf19,dmne21,dmna23,dmal26g, &
    dmal27,dmsi28,dmneut,dmprot,dmbid,dmbid1,xm,ym,y3m,xc12m,xc13m,xn14m,xn15m,xo16m,xo17m,xo18m,xne20m,xne22m,xmg24m,xmg25m, &
    xmg26m,xc14m,xf18m,xf19m,xne21m,xna23m,xal27m,xsi28m,xneutm,xprotm,xbidm,xbid1m,d,dbis,sc,ep23,s23,s20,so,sumn21,sumn23, &
    suma26,suma27,sums28,sumpro,dmn21,dmn23,dma26,dma27,dms28,dmpro,xal26m,xal26gm
real(kindreal),dimension(mbelx),private:: dmabelx,mabelx,sumabelx

private
public :: netnew,chemeps,chemold

contains
!======================================================================
subroutine netnew
!-----------------------------------------------------------------------
! This routine computes the changes in chemical composition due to reaction
! rates.
!
! It solves a reaction network including pp-chains and cno-tricycle for
! abundances at t(n+1) by a fully implicit finite-difference method similar
! to that of Arnett+Truran (1969).
!
! If the NaNe-MgAl cycle is not followed, it puts to equilibrium the abundances
! of some elements above a given temperature.
!
! This routine is called within each henyey iteration
!
! It calls CHEMIE: routine that homogenises the convective zones. It is called
! only after the (itminc-1) first iterations and during the last iteration (itminc=1).
!
! The last call to NETWNEW, and hence to CHEMIE, is done to get an estimation
! of the chemical composition of the next model.
!
! nrband : number of intermediate time steps between model (n) and (n+1)
! (default=1)
!
! nbchx : number of iterations for the computation of chemical composition change
! for the combustion of hydrogene.
!     nbchx = 1 : implicit methode
!     For He-b : nbchx = 24

! Calcul du modele courant.
!---------------------------
! idern = 0.
! Le calcul dans NETNEW n'est effectue que lors des trois premieres
!   iterations (iter=1,2,3). CHEMIE homogeneise les zones convectives
!   et DIFFBE ou DIFFUSION traitent la diffusion.
!   On obtient alors les nouvelles abondances qui vont etre utilisees
!   pour calculer la structure interne du modele.
! La diffusion des elements "primordiaux" pour l'evolution (H,3He,4He)
!   doit etre traitee a ce moment, car ces elements interviennent
!   dans la structure interne de l'etoile.

! Approximation de la composition chimique du modele suivant.
!-------------------------------------------------------------
! idern = 1.
! L'appel a NETWKI (puis a CHEMIE, et a DIFFBE ou DIFFUSION)
!   est fait dans MAIN pour estimer la composition au pas temporel
!   suivant.
! La diffusion des elements "tests" (7Li,9Be) se fait ici, lorsque
!   l'on connait la structure interne du modele courant.

! Les taux de reactions nucleaires sont re-actualises.

! Derniere version : 22 janvier 1993
!-----------------------------------------------------------------------
  implicit none

  integer:: l,nbb,llim=0,ii,lw,lal26,ns,lflag,flag_girl
  real(kindreal),parameter:: tvieal=3.2786885d+13
  real(kindreal):: xsubd,ddeit,smev,smas,zs,dms,sm63,tbasec=0.d0

  real(kindreal),dimension(ldi):: d2
!-----------------------------------------------------------------------
  lflag = 0
  flag_girl = 0
  d2(:)=0.d0

! idern = 0 : Calcul du modele courant.
!             Appel de netwki dans henyey.
!       = 1 : Estimation de la composition chimique du modele suivant.
!             Appel de netwki dans main.
  if (idern /= 1) then
    if (ialflu == 1) then
      if (x(m) /= 0.d0) then
        nbb = nbchx ! nbchx default value = 200
      else
        nbb = 24
      endif
    else
      nbb = nbchx
    endif
    if (alter <= 0.d0 .or. iter >= nbb) then
      return
    endif
  endif ! idern

  xsubd=real(nrband)
  ddeit=dzeit/xsubd

  do l=m,1,-1
   if (zensi(l) < 0.d0) then
     llim=l-2
     exit
   endif
  enddo

! initialisation cf journal m40.j2 ceci ne doit etre fait
! que lorsque l'on diffuse les especes chimiques
  x(:)=vvx(:)
  y3(:)=vvy3(:)
  y(:)=vvy(:)
  xc12(:)=vvxc12(:)
  xc13(:)=vvxc13(:)
  xn14(:)=vvxn14(:)
  xn15(:)=vvxn15(:)
  xo16(:)=vvxo16(:)
  xo17(:)=vvxo17(:)
  xo18(:)=vvxo18(:)
  xne20(:)=vvxne20(:)
  xne22(:)=vvxne22(:)
  xmg24(:)=vvxmg24(:)
  xmg25(:)=vvxmg25(:)
  xmg26(:)=vvxmg26(:)
  if (ialflu==1) then
    xc14(:)=vvxc14(:)
    xf18(:)=vvxf18(:)
    xf19(:)=vvxf19(:)
    xne21(:)=vvxne21(:)
    xna23(:)=vvxna23(:)
    xal26(:)=vvxal26g(:)
    xal27(:)=vvxal27(:)
    xsi28(:)=vvxsi28(:)
    xneut(:)=vvxneut(:)
    xprot(:)=vvxprot(:)
    xbid(:)=vvxbid(:)
    xbid1(:)=vvxbid1(:)
  endif
  do ii=1,nbelx
   abelx(ii,:)=vvabelx(ii,:)
  enddo

  if (x(m) > 0.d0) then
    if (zensi(m-3) > 0.d0) then
      smev=0.d0
      smas=0.d0
      do lw=m-1,1,-1
       if (lw <= (llim+2)) then
         exit
       endif
       zs=0.5d0*(zensi(lw)+zensi(lw+1))
       dms=exp(q(lw+1))-exp(q(lw))
       smas=smas+dms
       smev=smev+zs*dms
      enddo
      if (smas /= 0.d0) then
        sm63=smev/smas
        write(io_logs,'(2x,a,1x,2(1x,f8.4))') 'ENERGIE PAR GR. TRANSF. X E-18 =',sm63,smas
      endif
    endif ! zensi
  endif ! x

  if (ialflu==1) then
  ! desintegration de l'aluminiun 26 dans les zones ou on ne passe
  ! pas dans neth_alu ou netflu
    do lal26=1,m
     if (xal26(lal26) == 0.d0 .or. y(lal26) == 0.d0) then
       cycle
     endif
     if (x(lal26) == 0.d0 .or. t(lal26) < log(4.d6)) then
       if (y(lal26) == 0.d0 .or. t(lal26) < 18.06398074d0) then
         xal26(lal26)=(1.d0/(1.d0+dzeit/tvieal))*vvxal26g(lal26)
         xmg26(lal26)=vvxmg26(lal26)+dzeit/(tvieal+dzeit)*vvxal26g(lal26)
       endif
     endif
    enddo
  endif

!-----------------------------------------------------------------------
! Loop on the interior layers.

  do ns=1,nrband
! loop from centre to surface:
   do l=m,1,-1

! case we use chemeps (ichem=1), homogeneisation of chemical composition
! else (ichem=0), go to 300 to skip homogeneisation
! NB: if ialflu=1 --> ichem=0 (beginning of main)
    if (ichem==1) then
      if (epsc(l) == 0.d0 .or. idifcon /= 1) then
        if (l == m) then
          tbasec=t(m)
          go to 300
        endif
        if (zensi(l) > 0.d0) then ! convective layer
          if (zensi(l+1) <= 0.d0) then
            tbasec=t(l)
            go to 300
          else
            if (tbasec > log(4.d6)) then
              if (x(l)<1.d-8 .and. y(l)<1.d-8 .and. (xc12(l)-xc12(l+1)>1.d-10 .or. xo16(l)-xo16(l+1)>1.d-10)) then
                write(io_logs,*)'better check,l= ',l, xc12(l)-xc12(l+1),xo16(l)-xo16(l+1)
              endif
              x(l)=x(l+1)
              y3(l)=y3(l+1)
              d2(l)=d2(l+1)
              y(l)=y(l+1)
              xc12(l)=xc12(l+1)
              xc13(l)=xc13(l+1)
              xn14(l)=xn14(l+1)
              xn15(l)=xn15(l+1)
              xo16(l)=xo16(l+1)
              xo17(l)=xo17(l+1)
              xo18(l)=xo18(l+1)
              xne20(l)=xne20(l+1)
              xne22(l)=xne22(l+1)
              xmg25(l)=xmg25(l+1)
              xmg26(l)=xmg26(l+1)
              xmg24(l)=xmg24(l+1)
              do ii=1,nbelx
               abelx(ii,l)=abelx(ii,l+1)
              enddo
              cycle
            endif ! tbasec
          endif ! inside convective layer
        endif ! zensi
      endif ! epsc or idifcon
    endif ! ichem

300 if (t(l) <= log(4.d6)) then
      cycle
    endif
    if (ipop3 == 0) then
      if (x(l) > 1.d-9.and.epsy(l) > 0.d0) then
        if (verbose) then
          print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          print*,'net',l,x(l),vvx(l),epsy(l),idern
          print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        endif
        write(io_logs,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(io_logs,'(a,i5,2(1x,f14.10),1x,d14.8,i3)') 'net ',l,x(l),vvx(l),epsy(l),idern
        write(io_logs,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      endif
    endif

    if (x(l) > 0.d0) then
      lflag = 0
      flag_girl = 0
      select case(ialflu)
      case (0)
        call neth(l,ns,llim,ddeit,lflag,flag_girl)
      case (1)
        call neth_alu(l,ns,llim,ddeit,lflag,flag_girl)
      case default
        stop 'Bad value for ialflu, should be 0 or 1'
      end select
      if (lflag /= 0) then
        exit
      endif

!-----------------------------------------------------------------------
    else  ! x(l)
! HE-BURNING
      if (epsy(l) > 0.d0) then
        if (y(l) <= 0.d0) then
          cycle
        endif
        flag_girl = 0
        select case (ialflu)
        case (0)
          call nethe(l,ns,ddeit,flag_girl)
        case (1)
          call nethe_alu(l,ns,ddeit,flag_girl)
        case default
          stop 'Bad value for ialflu, should be 0 or 1'
        end select

!-----------------------------------------------------------------------
      else   ! y(l)
! C-BURNING
        if (abs(epsc(l)) <= 0.d0) then
          cycle
        endif
!  assumes nrband=1
        call netc(l,ddeit)
      endif   ! y(l)
    endif   ! x(l)
   enddo ! l
  enddo ! ns

! Traitement du melange dans les zones convectives.
! Traitement de la diffusion.
  call chemie

  return

end subroutine netnew
!======================================================================
subroutine neth(l,ns,llim,ddeit,lflag,flag_girl)
!-----------------------------------------------------------------------
  use inputparam,only: ipop3
  use SmallFunc,only: girl

  implicit none

  integer,intent(in):: l,ns,llim
  real(kindreal),intent(in):: ddeit

  integer,intent(out):: lflag
  integer,intent(inout):: flag_girl

  integer,parameter:: idimneth=15

  integer:: iSE,jSE
  real(kindreal):: wy3_net,wy15_net,wy17_net,wy18_net,wy12_net,wy13_net, &
                   t11v1,t33v2,t33d,t34,t343,t112,t113,t114,t115a,t115g,t116,t117a,t117g,t118a,t118g,t44=0.d0, &
                   tc124=0.d0,t134=0.d0,t144=0.d0,t164=0.d0,t174=0.d0,t184=0.d0,t204=0.d0,t224n=0.d0,t224g=0.d0, &
                   av,av4,dweit
  real(kindreal),dimension(idimneth):: vyab
  real(kindreal),dimension(ldi):: d2
  real(kindreal),dimension(idimneth,idimneth+1):: a
  real(kindreal),dimension(idimneth,1):: c
!-----------------------------------------------------------------------
  lflag = 0
  flag_girl = 0
  vyab(:)=0.d0
  c(:,:)=0.d0
  wy3_net=1.d0
  wy15_net=1.d0
  wy17_net=1.d0
  wy18_net=1.d0
  wy12_net=1.d0
  wy13_net=1.d0

  if (ns == 1) then
    if (l < llim) then
      dweit=dzeit
    else
      dweit=ddeit
    endif
    vyab(1)=vvx(l)
    vyab(2)=vvy3(l)/3.d0
    d2(l)=0.d0
    vyab(3)=vvy(l)/4.d0
    vyab(4)=vvxc12(l)/12.d0
    vyab(5)=vvxc13(l)/13.d0
    vyab(6)=vvxn14(l)/14.d0
    vyab(7)=vvxn15(l)/15.d0
    vyab(8)=vvxo16(l)/16.d0
    vyab(9)=vvxo17(l)/17.d0
    vyab(10)=vvxo18(l)/18.d0
    if (ipop3 == 1)then
      vyab(11)=vvxne20(l)/20.d0
      vyab(12)=vvxne22(l)/22.d0
      vyab(13)=vvxmg24(l)/24.d0
      vyab(14)=vvxmg25(l)/25.d0
      vyab(15)=vvxmg26(l)/26.d0
    endif
    if (ns == nrband) then
      if (l == m) then
        write(io_logs,'(1x,a,i5,3(1x,f10.7),13(1x,e12.5))') 'BEFORE NETH',l,vvx(l),vvy3(l),vvy(l),vvxc12(l), &
                   vvxc13(l),vvxn14(l),vvxn15(l),vvxo16(l),vvxo17(l),vvxo18(l),vvxne20(l),vvxne22(l),vvxmg24(l), &
                   vvxmg25(l),vvxmg26(l),d2(l)
      endif
    endif

  else

    if (l < llim) then
      lflag = 1
      return
    endif
    dweit=ddeit
    vyab(1)=x(l)
    vyab(2)=y3(l)/3.d0
    vyab(3)=y(l)/4.d0
    vyab(4)=xc12(l)/12.d0
    vyab(5)=xc13(l)/13.d0
    vyab(6)=xn14(l)/14.d0
    vyab(7)=xn15(l)/15.d0
    vyab(8)=xo16(l)/16.d0
    vyab(9)=xo17(l)/17.d0
    vyab(10)=xo18(l)/18.d0
    if (ipop3 == 1) then
      vyab(11)=xne20(l)/20.d0
      vyab(12)=xne22(l)/22.d0
      vyab(13)=xmg24(l)/24.d0
      vyab(14)=xmg25(l)/25.d0
      vyab(15)=xmg26(l)/26.d0
    endif
    if (ns == nrband) then
      if (l == m) then
        write(io_logs,'(1x,a,i5,3(1x,f10.7),13(1x,e12.5))') 'BEFORE NETH',l,x(l),y3(l),y(l),xc12(l),xc13(l),xn14(l), &
                   xn15(l),xo16(l),xo17(l),xo18(l),xne20(l),xne22(l),xmg24(l),xmg25(l),xmg26(l),d2(l)
      endif
    endif
  endif

! bij(l) : [ij](l) : taux de la reaction (i,j).
! tij    : [ij](l)*pas de temps*Yi.
  t11v1=b11(l)*dweit*vyab(1)
  t33v2=b33(l)*dweit*vyab(2)
  t33d=b33(l)*dweit*(vyab(2)+d2(l))
  t34=b34(l)*dweit
  t343=t34*vyab(3)
  t112=b112(l)*dweit
  t113=b113(l)*dweit
  t114=b114(l)*dweit
  t115a=b115a(l)*dweit
  t115g=b115g(l)*dweit
  t116=b116(l)*dweit
  t117a=b117a(l)*dweit
  t117g=b117g(l)*dweit
  t118a=b118a(l)*dweit
  t118g=b118g(l)*dweit
  if (ipop3 == 1) then
    t44=epsyy(l)*dweit
    tc124=epsyc(l)*dweit
    t134=c134(l)*dweit
    t144=c144(l)*dweit
    t164=epsyo(l)*dweit
    t174=e17an(l)*dweit
    t184=c184(l)*dweit
    t204=e20ag(l)*dweit
    t224n=c224(l)*dweit
    t224g=c224g(l)*dweit
  endif

! Initialisation a zero de la matrice du reseau de reactions nucleaires
! pour chaque coquille. Diagonale initialisee a 1
  a=0.d0
  do i=1,idimneth
   a(i,i)=1.d0
  enddo

! Ecriture de la matrice a(y,z) du reseau de reactions nucleaires.
! z : numero de la colonne.
! z <= idimneth  : membre de droite (Yi, instant(n+1)) du reseau.
! z = idimneth+1 : membre de gauche (Yi, instant(n)) du reseau.

! Ligne 1 : Variation de l'hydrogene.
  a(1,1)=1.d0+3.d0*t11v1+t112*vyab(4)+t113*vyab(5)+t114*vyab(6)+(t115a+t115g)*vyab(7)+t116*vyab(8)+(t117a+t117g)*vyab(9)+ &
         (t118a+2.d0*t118g)*vyab(10)
  a(1,2)=-t33d-t33v2+t343
  a(1,3)=t34*vyab(2)
  a(1,4)=t112*vyab(1)
  a(1,5)=t113*vyab(1)
  a(1,6)=t114*vyab(1)
  a(1,7)=(t115a+t115g)*vyab(1)
  a(1,8)=t116*vyab(1)
  a(1,9)=(t117a+t117g)*vyab(1)
  a(1,10)=(t118a+2.d0*t118g)*vyab(1)
  av=0.d0
  do i=3,10
   av=av+a(1,i)*vyab(i)
  enddo
  a(1,idimneth+1)=(1.d0+1.5d0*t11v1)*vyab(1)-t33d*vyab(2)+av
! Ligne 2 : Variation de l'helium 3.
  a(2,1)=-t11v1
  a(2,2)=wy3_net+t33v2+t343+t33d
  a(2,3)=a(1,3)
  a(2,idimneth+1)=(wy3_net+t33d)*vyab(2)-0.5d0*t11v1*vyab(1)+a(1,3)*vyab(3)
! Ligne 3 : Variation de l'helium 4.
  a(3,1)=-t115a*vyab(7)-t117a*vyab(9)-(t118a+t118g)*vyab(10)
  a(3,2)=-0.5d0*t33d-0.5d0*t33v2-t343
  a(3,7)=-t115a*vyab(1)
  if (ipop3 == 1) then
    a(3,3)=1.d0-a(1,3)+1.5d0*t44*vyab(3)*vyab(3)+tc124*vyab(4)+t134*vyab(5)+t144*vyab(6)+t164*vyab(8)+t184*vyab(10)+ &
           (t224n+t224g)*vyab(12)+t174*vyab(9)+t204*vyab(11)
    a(3,4)=tc124*vyab(3)
    a(3,5)=t134*vyab(3)
    a(3,6)=t144*vyab(3)
    a(3,8)=t164*vyab(3)
    a(3,9)=-t117a*vyab(1)+t174*vyab(3)
    a(3,10)=-(t118a+t118g)*vyab(1)+t184*vyab(3)
    a(3,11)=t204*vyab(3)
    a(3,12)=(t224n+t224g)*vyab(3)
    av4=0.d0
    do i=4,12
     av4=av4+a(3,i)*vyab(i)
    enddo
    a(3,idimneth+1)=(1.d0+t44*vyab(3)*vyab(3)-t34*vyab(2))*vyab(3)-0.5d0*t33d*vyab(2)+av4
  else
    a(3,3)=1.d0-a(1,3)
    a(3,9)=-t117a*vyab(1)
    a(3,10)=-(t118a+t118g)*vyab(1)
    a(3,idimneth+1)=a(3,3)*vyab(3)-0.5d0*t33d*vyab(2)+a(3,7)*vyab(7)+a(3,9)*vyab(9)+a(3,10)*vyab(10)
  endif
! Ligne 4 : Variation du carbone 12.
  a(4,1)=t112*vyab(4)-t115a*vyab(7)
  a(4,7)=a(3,7)
  if (ipop3 == 1) then
    a(4,3)=-0.5d0*t44*vyab(3)*vyab(3)+tc124*vyab(4)
    a(4,4)=wy12_net+a(1,4)+tc124*vyab(3)
    a(4,idimneth+1)=a(4,4)*vyab(4)+a(4,7)*vyab(7)-(1.d0/3.d0)*t44*vyab(3)*vyab(3)*vyab(3)
  else
    a(4,4)=wy12_net+a(1,4)
    a(4,idimneth+1)=a(4,4)*vyab(4)+a(4,7)*vyab(7)
  endif
! Ligne 5 : Variation du carbone 13.
  a(5,1)=-t112*vyab(4)+t113*vyab(5)
  a(5,4)=-a(1,4)
  if (ipop3 == 1) then
    a(5,3)=t134*vyab(5)
    a(5,5)=wy13_net+a(1,5)+t134*vyab(3)
  else
    a(5,5)=wy13_net+a(1,5)
  endif
  a(5,idimneth+1)=a(5,5)*vyab(5)-a(1,4)*vyab(4)
! Ligne 6 : Variation de l'azote 14.
  a(6,1)=-t113*vyab(5)+t114*vyab(6)-t117a*vyab(9)
  a(6,5)=-a(1,5)
  a(6,9)=-t117a*vyab(1)
  if (ipop3 == 1) then
    a(6,3)=t144*vyab(6)
    a(6,6)=1.d0+a(1,6)+t144*vyab(3)
  else
    a(6,6)=1.d0+a(1,6)
  endif
  a(6,idimneth+1)=a(6,6)*vyab(6)-a(1,5)*vyab(5)+a(6,9)*vyab(9)
! Ligne 7 : Variation de l'azote 15.
  a(7,1)=-t114*vyab(6)+(t115a+t115g)*vyab(7)-t118a*vyab(10)
  a(7,6)=-a(1,6)
  a(7,7)=wy15_net+a(1,7)
  a(7,10)=-t118a*vyab(1)
  a(7,idimneth+1)=a(7,7)*vyab(7)-a(1,6)*vyab(6)+a(7,10)*vyab(10)
! Ligne 8 : Variation de l'oxygene 16.
  a(8,1)=-t115g*vyab(7)+t116*vyab(8)-t118g*vyab(10)
  a(8,7)=-t115g*vyab(1)
  a(8,10)=-t118g*vyab(1)
  if (ipop3 == 1) then
    a(8,3)=-tc124*vyab(4)-t134*vyab(5)+t164*vyab(8)
    a(8,4)=-tc124*vyab(3)
    a(8,5)=-t134*vyab(3)
    a(8,8)=1.d0+a(1,8)+t164*vyab(3)
    a(8,idimneth+1)=a(8,8)*vyab(8)+a(8,7)*vyab(7)+a(8,10)*vyab(10)+a(8,4)*vyab(4)+a(8,5)*vyab(5)
  else
    a(8,8)=1.d0+a(1,8)
    a(8,idimneth+1)=a(8,8)*vyab(8)+a(8,7)*vyab(7)+a(8,10)*vyab(10)
  endif
! Ligne 9 : Variation de l'oxygene 17.
  a(9,1)=-t116*vyab(8)+(t117a+t117g)*vyab(9)
  a(9,8)=-a(1,8)
  if (ipop3 == 1) then
    a(9,3)=t174*vyab(9)
    a(9,9)=wy17_net+a(1,9)+t174*vyab(3)
  else
    a(9,9)=wy17_net+a(1,9)
  endif
  a(9,idimneth+1)=a(9,9)*vyab(9)-a(1,8)*vyab(8)
! Ligne 10 : Variation de l'oxygene 18.
  a(10,1)=-t117g*vyab(9)+(t118a+t118g)*vyab(10)
  a(10,9)=-t117g*vyab(1)
  if (ipop3 == 1) then
    a(10,3)=-t144*vyab(6)+t184*vyab(10)
    a(10,6)=-t144*vyab(3)
    a(10,10)=wy18_net+t184*vyab(3)+(t118a+t118g)*vyab(1)
    a(10,idimneth+1)=a(10,10)*vyab(10)+a(10,9)*vyab(9)+a(10,6)*vyab(6)
  else
    a(10,10)=wy18_net-a(3,10)
    a(10,idimneth+1)=a(10,10)*vyab(10)+a(10,9)*vyab(9)
  endif

  if (ipop3 == 1) then
! Ligne 11 : Variation du neon 20.
    a(11,3)=-t164*vyab(8)+t204*vyab(11)-t174*vyab(9)
    a(11,8)=-t164*vyab(3)
    a(11,9)=-t174*vyab(3)
    a(11,11)=1.d0+t204*vyab(3)
    a(11,idimneth+1)=a(11,11)*vyab(11)+a(11,8)*vyab(8)+a(11,9)*vyab(9)
! Ligne 12 : Variation du neon 22.
    a(12,3)=-t184*vyab(10)+(t224n+t224g)*vyab(12)
    a(12,10)=-t184*vyab(3)
    a(12,12)=1.d0+(t224n+t224g)*vyab(3)
    a(12,idimneth+1)=a(12,12)*vyab(12)+a(12,10)*vyab(10)
! Ligne 13 : Variation du magnesium 24.
    a(13,3)=-t204*vyab(11)
    a(13,11)=-t204*vyab(3)
    a(13,13)=1.d0
    a(13,idimneth+1)=vyab(13)+a(13,11)*vyab(11)
! Ligne 14 : Variation du magnesium 25.
    a(14,3)=-t224n*vyab(12)
    a(14,12)=-t224n*vyab(3)
    a(14,14)=1.d0
    a(14,idimneth+1)=vyab(14)+a(14,12)*vyab(12)
! Ligne 15 : Variation du magnesium 26.
    a(15,3)=-t224g*vyab(12)
    a(15,12)=-t224g*vyab(3)
    a(15,15)=1.d0
    a(15,idimneth+1)=vyab(15)+a(15,12)*vyab(12)
  endif   ! ipop3

!---  INVERT MATRIX A (15X15) AND MULTIPLY BY THE R.H.S. (KNOWN TERMS)
!     A(idimneth,idimneth+1) . AFTER MULTIPLICATION, THE FIRST COLUMN OF VECTOR C WILL
!     CONTAIN THE NEW ABUNDANCES AT TIME T(N+1)
  flag_girl = 0
  call girl(a,c,idimneth,1,flag_girl)
  if (flag_girl /= 0) then
    if (idebug>0) then
      write(*,'("neth - matrix a(",i2","i2"),flag:",i1)') idimneth,idimneth+1,flag_girl
      do iSE=1,idimneth
       do jSE=1,idimneth+1
        write(*,'("a(",i2,",",i2,") :",d22.12)') iSE,jSE,a(iSE,jSE)
       enddo
      enddo
    endif
    rewind(io_runfile)
    write(io_runfile,*) nwmd,':girl crashes in neth with matrix a(15,16)'
    stop
  endif

! Nouvelles abondances dues a la combustion de l'hydrogene.
! Fractions de masse (Xi = Yi*Ai).
  x(l)=c(1,1)

  if (ns == nrband) then
    if (x(l) <= 1.d-09 .and. idern == 1) then
      x(l)=0.d0
    endif
  endif

  y3(l)=c(2,1)*3.d0
  d2(l)=(y3(l)/3.d0)-vyab(2)
  y(l)=c(3,1)*4.d0
  xc12(l)=c(4,1)*12.d0
  xc13(l)=c(5,1)*13.d0
  xn14(l)=c(6,1)*14.d0
  xn15(l)=c(7,1)*15.d0
  xo16(l)=c(8,1)*16.d0
  xo17(l)=c(9,1)*17.d0
  xo18(l)=c(10,1)*18.d0
  if (ipop3 == 1) then
    xne20(l)=c(11,1)*20.d0
    xne22(l)=c(12,1)*22.d0
    xmg24(l)=c(13,1)*24.d0
    xmg25(l)=c(14,1)*25.d0
    xmg26(l)=c(15,1)*26.d0
  endif

  if (ns == nrband) then
    if (l == m) then
      write(io_logs,'(1x,a,i5,15(1x,e12.5))') 'AFTER NETH ',l,x(l),y3(l)/3.d0,y(l)/4.d0,xc12(l)/12.d0,xc13(l)/13.d0, &
                 xn14(l)/14.d0,xn15(l)/15.d0,xo16(l)/16.d0,xo17(l)/17.d0,xo18(l)/18.d0,xne20(l)/20.d0,xne22(l)/22.d0, &
                 xmg24(l)/24.d0,xmg25(l)/25.d0,xmg26(l)/26.d0
    endif
  endif

  return
end subroutine neth
!======================================================================
subroutine neth_alu(l,ns,llim,ddeit,lflag,flag_girl)
!-----------------------------------------------------------------------
  use SmallFunc,only: girl

  implicit none

  integer,intent(in):: l,ns,llim
  real(8),intent(in):: ddeit
  integer,intent(out):: lflag
  integer,intent(inout):: flag_girl

  integer,parameter:: idimnetha=21

  integer:: i,ini,iSE,jSE
  real(kindreal):: wy3_net,wy15_net,wy17_net,wy18_net,wy12_net,wy13_net,t34,t343,t112,t113,t114,t115a,t115g,t116,t117a,t117g, &
    t118a,t118g,t44=0.d0,tc124=0.d0,t134=0.d0,t144=0.d0,t164=0.d0,t174=0.d0,t184=0.d0,t204=0.d0,t224n=0.d0,t224g=0.d0,av,av4, &
    t11v1,t33v2,t119g,t119a,t120,t121,t122,t123g,t123a,t124,t125g,t125m,t1mg26,t1al26,t127g,t127a,t15ag=0.d0,t174g=0.d0, &
    t18an=0.d0,t19ap=0.d0,t214=0.d0,t24ag=0.d0,t25an=0.d0,dweit

  real(kindreal), dimension(idimnetha):: vyab
  real(kindreal), dimension(ldi):: d2
  real(kindreal), dimension(idimnetha,1):: c
  real(kindreal), dimension(idimnetha,idimnetha+1):: b
!-----------------------------------------------------------------------
  lflag=0
  flag_girl=0
  wy3_net=1.d0
  wy15_net=1.d0
  wy17_net=1.d0
  wy18_net=1.d0
  wy12_net=1.d0
  wy13_net=1.d0

  if (ns == 1) then
    if (l < llim) then
      dweit=dzeit
    else
      dweit=ddeit
    endif
    vyab(1)=vvx(l)
    vyab(2)=vvy3(l)/3.d0
    d2(l)=0.d0
    vyab(3)=vvy(l)/4.d0
    vyab(4)=vvxc12(l)/12.d0
    vyab(5)=vvxc13(l)/13.d0
    vyab(6)=vvxn14(l)/14.d0
    vyab(7)=vvxn15(l)/15.d0
    vyab(8)=vvxo16(l)/16.d0
    vyab(9)=vvxo17(l)/17.d0
    vyab(10)=vvxo18(l)/18.d0
    vyab(11)=vvxf19(l)/19.d0
    vyab(12)=vvxne20(l)/20.d0
    vyab(13)=vvxne21(l)/21.d0
    vyab(14)=vvxne22(l)/22.d0
    vyab(15)=vvxna23(l)/23.d0
    vyab(16)=vvxmg24(l)/24.d0
    vyab(17)=vvxmg25(l)/25.d0
    vyab(18)=vvxmg26(l)/26.d0
    vyab(19)=vvxal26g(l)/26.d0
    vyab(20)=vvxal27(l)/27.d0
    vyab(21)=vvxsi28(l)/28.d0

    if (ns == nrband) then
      if (l == m) then
        write(io_logs,'(1x,a,i5,1x,f10.7,1x,e8.2,1x,f10.7,13(1x,e8.2))') 'AVANT NETWKI',l,vvx(l),vvy3(l),vvy(l), &
                   vvxc12(l),vvxc13(l),vvxn14(l),vvxn15(l),vvxo16(l),vvxo17(l),vvxo18(l),vvxne20(l),vvxne22(l), &
                   vvxmg24(l),vvxmg25(l),vvxmg26(l),d2(l)
      endif
    endif

  else
    if (l < llim) then
      lflag = 1
      return
    endif
    dweit=ddeit
    vyab(1)=x(l)
    vyab(2)=y3(l)/3.d0
    vyab(3)=y(l)/4.d0
    vyab(4)=xc12(l)/12.d0
    vyab(5)=xc13(l)/13.d0
    vyab(6)=xn14(l)/14.d0
    vyab(7)=xn15(l)/15.d0
    vyab(8)=xo16(l)/16.d0
    vyab(9)=xo17(l)/17.d0
    vyab(10)=xo18(l)/18.d0
    vyab(11)=xf19(l)/19.d0
    vyab(12)=xne20(l)/20.d0
    vyab(13)=xne21(l)/21.d0
    vyab(14)=xne22(l)/22.d0
    vyab(15)=xna23(l)/23.d0
    vyab(16)=xmg24(l)/24.d0
    vyab(17)=xmg25(l)/25.d0
    vyab(18)=xmg26(l)/26.d0
    vyab(19)=xal26(l)/26.d0
    vyab(20)=xal27(l)/27.d0
    vyab(21)=xsi28(l)/28.d0

    if (ns == nrband) then
      if (l == m) then
        write(io_logs,'(1x,a,i5,1x,f10.7,1x,e8.2,1x,f10.7,13(1x,e8.2))') 'AVANT NETH_ALU',l,x(l),y3(l),y(l),xc12(l), &
                   xc13(l),xn14(l),xn15(l),xo16(l),xo17(l),xo18(l),xne20(l),xne22(l),xmg24(l),xmg25(l),xmg26(l),d2(l)
      endif
    endif
  endif

!---  COMPUTATION OF PREVIOUS ABUNDANCES
!---  T=B*DT
  t11v1=b11(l)*dweit*vyab(1)
  t33v2=b33(l)*dweit*vyab(2)
  t34=b34(l)*dweit
  t343=T34*vyab(3)
  t112=b112(l)*dweit
  t113=b113(l)*dweit
  t114=b114(l)*dweit
  t115a=b115a(l)*dweit
  t115g=b115g(l)*dweit
  t116=b116(l)*dweit
  t117a=b117a(l)*dweit
  t117g=b117g(l)*dweit
  t118a=b118a(l)*dweit
  t118g=b118g(l)*dweit
  t119g=b119g(l)*dweit
  t119a=b119a(l)*dweit
  t120=b120(l)*dweit
  t121=b121(l)*dweit
  t122=b122(l)*dweit
  t123g=b123g(l)*dweit
  t123a=b123a(l)*dweit
  t124=b124(l)*dweit
  t125g=b125g(l)*dweit
  t125m=b125m(l)*dweit
  t1mg26=b1mg26(l)*dweit
  t1al26=b1al26(l)*dweit
  t127g=b127g(l)*dweit
  t127a=b127a(l)*dweit
  if (ipop3 == 1) then
    t44=epsyy(l)*dweit
    tc124=epsyc(l)*dweit
    t134=c134(l)*dweit
    t144=c144(l)*dweit
    t15ag=e15ag(l)*dweit
    t164=epsyo(l)*dweit
    t174=e17an(l)*dweit
    t174g=e17ag(l)*dweit
    t184=c184(l)*dweit
    t18an=e18an(l)*dweit
    t19ap=e19ap(l)*dweit
    t204=e20ag(l)*dweit
    t214=e21ag(l)*dweit
    t224n=c224(l)*dweit
    t224g=c224g(l)*dweit
    t24ag=e24ag(l)*dweit
    t25an=e25an(l)*dweit
  endif

!---  INITIALIZE MATRIX A TO ZERO FOR EVERY MASS SHELL.
  b=0.d0

!--- EQUATION FOR H1
  b(1,1) = 1.d0+3.d0*t11V1+t112*vyab(4)+t113*vyab(5)+t114*vyab(6)+(t115A+t115G)*vyab(7)+t116*vyab(8)+(t117A+t117G)*vyab(9)+ &
          (t118A+t118G)*vyab(10)+(t119A+t119G)*vyab(11)+t120*vyab(12)+t121*vyab(13)+t122*vyab(14)+(t123A+t123G)*vyab(15)+ &
          t124*vyab(16)+(t125M+t125G)*vyab(17)+t1mg26*vyab(18)+t1al26*vyab(19)+(t127A+t127G)*vyab(20)
  b(1,2)=-2.d0*t33v2+t343
  if (ipop3 == 1) then
    b(1,3)=t34*vyab(2)-t19ap*vyab(11)
    b(1,11)=(t119A+t119G)*vyab(1)-t19ap*vyab(3)
    ini=4
  else
    b(1,3)=t34*vyab(2)
    b(1,11)=(t119A+t119G)*vyab(1)
    ini=3
  endif
  b(1,4)=t112*vyab(1)
  b(1,5)=t113*vyab(1)
  b(1,6)=t114*vyab(1)
  b(1,7)=(t115A+t115G)*vyab(1)
  b(1,8)=t116*vyab(1)
  b(1,9)=(t117A+t117G)*vyab(1)
  b(1,10)=(t118A+t118G)*vyab(1)
  b(1,12)=t120*vyab(1)
  b(1,13)=t121*vyab(1)
  b(1,14)=t122*vyab(1)
  b(1,15)=(t123A+t123G)*vyab(1)
  b(1,16)=t124*vyab(1)
  b(1,17)=(t125M+t125G)*vyab(1)
  b(1,18)=t1MG26*vyab(1)
  b(1,19)=t1AL26*vyab(1)
  b(1,20)=(t127A+t127G)*vyab(1)
  av=0.d0
  do i=ini,20
    av=av+b(1,i)*vyab(i)
  enddo
  b(1,22)=(1.d0+1.5d0*t11v1)*vyab(1)-t33v2*vyab(2)+av
  if (ipop3 == 1) then
    b(1,22)=b(1,22)+t34*vyab(2)*vyab(3)
  endif
!---  EQUATION FOR HE3
  b(2,1)=-t11v1
! wy3
  b(2,2)=wy3_net+2.d0*t33v2+t343
  b(2,3)=t34*vyab(2)
  b(2,22)=(wy3_net+t33v2)*vyab(2)-0.5d0*t11v1*vyab(1)+t34*vyab(2)*vyab(3)
!---  EQUATION FOR HE4
  b(3,1)=-t115a*vyab(7)-t117a*vyab(9)-t118a*vyab(10)-t119a*vyab(11)-t123a*vyab(15)-t127a*vyab(20)
  b(3,2)=-t33v2-t343
  b(3,3)=1.d0-t34*vyab(2)
  if (ipop3 == 1) then
    b(3,3)=b(3,3)+1.5d0*t44*vyab(3)*vyab(3)+tc124*vyab(4)+t134*vyab(5)+t144*vyab(6)+t15ag*vyab(7)+t164*vyab(8)+ &
           (t184+t18an)*vyab(10)+t19ap*vyab(11)+t204*vyab(12)+t214*vyab(13)+(t224n+t224g)*vyab(14)+t24ag*vyab(16)+t25an*vyab(17)
    b(3,4)=tc124*vyab(3)
    b(3,5)=t134*vyab(3)
    b(3,6)=t144*vyab(3)
    b(3,7)=-t115a*vyab(1)+t15ag*vyab(3)
    b(3,8)=t164*vyab(3)
    b(3,9)=-t117a*vyab(1)+(t174+t174g)*vyab(3)
    b(3,10)=-t118a*vyab(1)+(t184+t18an)*vyab(3)
    b(3,11)=-t119a*vyab(1)+t19ap*vyab(3)
    b(3,12)=t204*vyab(3)
    b(3,13)=t214*vyab(3)
    b(3,14)=(t224n+t224g)*vyab(3)
    b(3,16)=t24ag*vyab(3)
    b(3,17)=t25an*vyab(3)
  else
    b(3,7)=-t115a*vyab(1)
    b(3,9)=-t117a*vyab(1)
    b(3,10)=-t118a*vyab(1)
    b(3,11)=-t119a*vyab(1)
  endif
  b(3,15)=-t123a*vyab(1)
  b(3,20)=-t127a*vyab(1)
  if (ipop3 == 1) then
    av4=0.d0
    do i=4,20
     av4=av4+b(3,i)*vyab(i)
    enddo
    b(3,22)=(1.d0+t44*vyab(3)*vyab(3)-t34*vyab(2))*vyab(3)-0.5d0*t33v2*vyab(2)+av4
  else
    b(3,22)=b(3,3)*vyab(3)-0.5d0*t33v2*vyab(2)+b(3,7)*vyab(7)+b(3,9)*vyab(9)+b(3,10)*vyab(10)+b(3,11)*vyab(11)+ &
             b(3,15)*vyab(15)+b(3,20)*vyab(20)
  endif
!---  EQUATION FOR C12
  b(4,1)=t112*vyab(4)-t115a*vyab(7)
  if (ipop3 == 1) then
    b(4,3)=-0.5d0*t44*vyab(3)*vyab(3)+tc124*vyab(4)
    b(4,4)=wy12_net+t112*vyab(1)+tc124*vyab(3)
  else
    b(4,4)=wy12_net+t112*vyab(1)
  endif
  b(4,7)=-t115a*vyab(1)
  if (ipop3 == 1) then
    b(4,22)=b(4,4)*vyab(4)+b(4,7)*vyab(7)-(1.d0/3.d0)*t44*vyab(3)*vyab(3)*vyab(3)
  else
    b(4,22)=b(4,4)*vyab(4)+b(4,7)*vyab(7)
  endif
!---  EQUATION FOR C13
! immediate N13(,e+nu)C13 is supposed:
  b(5,1)=-t112*vyab(4)+t113*vyab(5)
  if (ipop3 == 1) then
    b(5,3)=t134*vyab(5)
    b(5,5)=wy13_net+t113*vyab(1)+t134*vyab(3)
  else
    b(5,5)=wy13_net+t113*vyab(1)
  endif
  b(5,4)=-t112*vyab(1)
  b(5,22)=b(5,5)*vyab(5)+b(5,4)*vyab(4)
!---  EQUATION POUR N14
  b(6,1)=-t113*vyab(5)+t114*vyab(6)-t117a*vyab(9)
  if (ipop3 == 1) then
    b(6,3)=t144*vyab(6)
    b(6,6)=1.d0+t114*vyab(1)+t144*vyab(3)
  else
    b(6,6)=1.d0+t114*vyab(1)
  endif
  b(6,5)=-t113*vyab(1)
  b(6,9)=-t117a*vyab(1)
  b(6,22)=b(6,6)*vyab(6)+b(6,5)*vyab(5)+b(6,9)*vyab(9)
!---  EQUATION FOR N15
! immediate O15(,e+nu)N15 is supposed:
  b(7,1)=-t114*vyab(6)+(t115a+t115g)*vyab(7)-t118a*vyab(10)
  if (ipop3 == 1) then
    b(7,3)=t15ag*vyab(7)
    b(7,7)=wy15_net+(t115a+t115g)*vyab(1)+t15ag*vyab(3)
  else
    b(7,7)=wy15_net+(t115a+t115g)*vyab(1)
  endif
  b(7,6)=-t114*vyab(1)
! wy15
  b(7,10)=-t118a*vyab(1)
  b(7,22)=b(7,7)*vyab(7)+b(7,6)*vyab(6)+b(7,10)*vyab(10)
!---  EQUATION FOR O16
  b(8,1)=-t115g*vyab(7)+t116*vyab(8)-t119a*vyab(11)
  if (ipop3 == 1) then
    b(8,3)=t164*vyab(8)-tc124*vyab(4)-t134*vyab(5)
    b(8,4)=-tc124*vyab(3)
    b(8,5)=-t134*vyab(3)
    b(8,8)=1.d0+t116*vyab(1)+t164*vyab(3)
  else
    b(8,8)=1.d0+t116*vyab(1)
  endif
  b(8,7)=-t115g*vyab(1)
  b(8,11)=-t119a*vyab(1)
  if (ipop3 == 1) then
    b(8,22)=b(8,8)*vyab(8)+b(8,4)*vyab(4)+b(8,5)*vyab(5)+b(8,7)*vyab(7)+b(8,11)*vyab(11)
  else
    b(8,22)=b(8,8)*vyab(8)+b(8,7)*vyab(7)+b(8,11)*vyab(11)
  endif
!---  EQUATION FOR O17
! immediate F17(,e+nu)O17 is supposed:
  b(9,1)=-t116*vyab(8)+(t117a+t117g)*vyab(9)
  if (ipop3 == 1) then
    b(9,3)=(t174+t174g)*vyab(9)
    b(9,9)=wy17_net+(t117a+t117g)*vyab(1)+(t174+t174g)*vyab(3)
  else
    b(9,9)=wy17_net+(t117a+t117g)*vyab(1)
  endif
  b(9,8)=-t116*vyab(1)
! wy17
  b(9,22)=b(9,9)*vyab(9)+b(9,8)*vyab(8)
!---  EQUATION FOR O18
! immediate F18(,e+nu)O18 is supposed:
  b(10,1)=-t117g*vyab(9)+(t118a+t118g)*vyab(10)
  if (ipop3 == 1) then
    b(10,3)=(t184+t18an)*vyab(10)-t144*vyab(6)
    b(10,6)=-t144*vyab(3)
    b(10,10)=wy18_net+(t118a+t118g)*vyab(1)+(t184+t18an)*vyab(3)
  else
    b(10,10)=wy18_net+t118a*vyab(1)
  endif
  b(10,9)=-t117g*vyab(1)
! wy18
  b(10,22)=b(10,10)*vyab(10)+b(10,9)*vyab(9)
  if (ipop3 == 1) then
    b(10,22)=b(10,22)+b(10,6)*vyab(6)
  endif
!---  EQUATION FOR F19
  b(11,1)=(t119a+t119g)*vyab(11)-t118g*vyab(10)
  if (ipop3 == 1) then
    b(11,3)=t19ap*vyab(11)-t15ag*vyab(7)
    b(11,7)=-t15ag*vyab(3)
    b(11,11)=1.d0+(t119a+t119g)*vyab(1)+t19ap*vyab(3)
  else
    b(11,11)=1.d0+(t119a+t119g)*vyab(1)
  endif
  b(11,10)=-t118g*vyab(1)
  b(11,22)=b(11,11)*vyab(11)+b(11,10)*vyab(10)
  if (ipop3 == 1) then
    b(11,22)=b(11,22)+b(11,7)*vyab(7)
  endif
!---  EQUATION FOR NE20
  b(12,1)=t120*vyab(12)-t119g*vyab(11)-t123a*vyab(15)
  if (ipop3 == 1) then
    b(12,3)=t204*vyab(12)-t164*vyab(8)-t174*vyab(9)
    b(12,8)=-t164*vyab(3)
    b(12,9)=-t174*vyab(3)
    b(12,12)=1.d0+t120*vyab(1)+t204*vyab(3)
  else
    b(12,12)=1.d0+t120*vyab(1)
  endif
  b(12,11)=-t119g*vyab(1)
  b(12,15)=-t123a*vyab(1)
  if (ipop3 == 1) then
    b(12,22)=b(12,12)*vyab(12)+b(12,8)*vyab(8)+b(12,9)*vyab(9)+b(12,11)*vyab(11)+b(12,15)*vyab(15)
  else
    b(12,22)=b(12,12)*vyab(12)+b(12,11)*vyab(11)+b(12,15)*vyab(15)
  endif
!---  EQUATION FOR NE21
! immediate NA21(,e+nu)NE21 is supposed:
  b(13,1)=t121*vyab(13)-t120*vyab(12)
  if (ipop3 == 1) then
    b(13,3)=t214*vyab(13)-t18an*vyab(10)
    b(13,10)=-t18an*vyab(3)
    b(13,13)=1.d0+t121*vyab(1)+t214*vyab(3)
  else
    b(13,13)=1.d0+t121*vyab(1)
  endif
  b(13,12)=-t120*vyab(1)
  if (ipop3 == 1) then
    b(13,22)=b(13,13)*vyab(13)+b(13,10)*vyab(10)+b(13,12)*vyab(12)
  else
    b(13,22)=b(13,13)*vyab(13)+b(13,12)*vyab(12)
  endif
!---  EQUATION FOR NE22
! immediate NA22(,e+nu)NE22 is supposed:
  b(14,1)=t122*vyab(14)-t121*vyab(13)
  if (ipop3 == 1) then
    b(14,3)=(t224n+t224g)*vyab(14)-t184*vyab(10)-t19ap*vyab(11)
    b(14,10)=-t184*vyab(3)
    b(14,11)=-t19ap*vyab(3)
    b(14,14)=1.d0+t122*vyab(1)+(t224n+t224g)*vyab(3)
  else
    b(14,14)=1.d0+t122*vyab(1)
  endif
  b(14,13)=-t121*vyab(1)
  if (ipop3 == 1) then
    b(14,22)=b(14,14)*vyab(14)+b(14,10)*vyab(10)+b(14,11)*vyab(11)+b(14,13)*vyab(13)
  else
    b(14,22)=b(14,14)*vyab(14)+b(14,13)*vyab(13)
  endif
!---  EQUATION FOR NA23
  b(15,1)=(t123a+t123g)*vyab(15)-t122*vyab(14)
  b(15,14)=-t122*vyab(1)
  b(15,15)=1.d0+(t123a+t123g)*vyab(1)
  b(15,22)=b(15,15)*vyab(15)+b(15,14)*vyab(14)
!---  EQUATION FOR MG24
  b(16,1)=t124*vyab(16)-t123g*vyab(15)-t127a*vyab(20)
  if (ipop3 == 1) then
    b(16,3)=t24ag*vyab(16)-t204*vyab(12)
    b(16,12)=-t204*vyab(3)
    b(16,16)=1.d0+t124*vyab(1)+t24ag*vyab(3)
  else
    b(16,16)=1.d0+t124*vyab(1)
  endif
  b(16,15)=-t123g*vyab(1)
  b(16,20)=-t127a*vyab(1)
  if (ipop3 == 1) then
    b(16,22)=b(16,16)*vyab(16)+b(16,12)*vyab(12)+b(16,15)*vyab(15)+b(16,20)*vyab(20)
  else
    b(16,22)=b(16,16)*vyab(16)+b(16,15)*vyab(15)+b(16,20)*vyab(20)
  endif
!---  EQUATION FOR MG25
! immediate AL25(,e+nu)MG25 is supposed:
  b(17,1)=(t125g+t125m)*vyab(17)-t124*vyab(16)
  if (ipop3 == 1) then
    b(17,3)=t25an*vyab(17)-t224n*vyab(14)-t214*vyab(13)
    b(17,13)=-t214*vyab(3)
    b(17,14)=-t224n*vyab(3)
    b(17,17)=1.d0+(t125g+t125m)*vyab(1)+t25an*vyab(3)
  else
    b(17,17)=1.d0+(t125g+t125m)*vyab(1)
  endif
  b(17,16)=-t124*vyab(1)
  if (ipop3 == 1) then
    b(17,22)=b(17,17)*vyab(17)+b(17,13)*vyab(13)+b(17,14)*vyab(14)+b(17,16)*vyab(16)
  else
    b(17,22)=b(17,17)*vyab(17)+b(17,16)*vyab(16)
  endif
!---  EQUATION FOR MG26
! immediate AL26M(,e+nu)MG26 is supposed:
  b(18,1)=t1mg26*vyab(18)-t125m*vyab(17)
  if (ipop3 == 1) then
    b(18,3)=-t224g*vyab(14)
    b(18,14)=-t224g*vyab(3)
  endif
  b(18,17)=-t125m*vyab(1)
  b(18,18)=1.d0+t1mg26*vyab(1)
! B DECAY OF AL26G : 1/TAU=2.98E-14 [S-1]
  b(18,19)=-3.05d-14*dweit
  b(18,22)=b(18,18)*vyab(18)+b(18,17)*vyab(17)
  if (ipop3 == 1) then
    b(18,22)=b(18,22)+b(18,14)*vyab(14)
  endif
!---  EQUATION FOR AL26G
  b(19,1)=t1al26*vyab(19)-t125g*vyab(17)
  b(19,17)=-t125g*vyab(1)
  b(19,19)=1.d0+t1al26*vyab(1)+3.05d-14*dweit
  b(19,22)=(1.d0+t1al26*vyab(1))*vyab(19)+b(19,17)*vyab(17)
!---  EQUATION FOR AL27
! immediate SI27(,e+nu)AL27 is supposed:
  b(20,1)=(t127a+t127g)*vyab(20)-t1mg26*vyab(18)-t1al26*vyab(19)
  b(20,18)=-t1mg26*vyab(1)
  b(20,19)=-t1al26*vyab(1)
  b(20,20)=1.d0+(t127a+t127g)*vyab(1)
  b(20,22)=b(20,20)*vyab(20)+b(20,18)*vyab(18)+b(20,19)*vyab(19)
!---  EQUATION FOR SI28
  b(21,1)=-t127g*vyab(20)
  if (ipop3 == 1) then
    b(21,3)=-t24ag*vyab(16)-t25an*vyab(17)
    b(21,16)=-t24ag*vyab(3)
    b(21,17)=-t25an*vyab(3)
  endif
  b(21,20)=-t127g*vyab(1)
  b(21,21)=1.d0
  b(21,22)=vyab(21)+b(21,20)*vyab(20)
  if (ipop3 == 1) then
    b(21,22)=b(21,22)+b(21,16)*vyab(16)+b(21,17)*vyab(17)
  endif

!---  INVERT MATRIX A (21X21) AND MULTIPLY BY THE R.H.S. (KNOWN TERMS)
!     OF THE DIFFERENCE EQNS. WHICH ARE STORED IN b(1,22),b(2,22),...,
!     b(21,22) . AFTER MULTIPLICATION, THE FIRST COLUMN OF VECTOR C WILL
!     CONTAIN THE NEW ABUNDANCES AT TIME T(N+1)
  flag_girl = 0
  call girl(b,c,idimnetha,1,flag_girl)
  if (flag_girl /= 0) then
    if (idebug>0) then
      write(*,'("neth_alu - matrix b(",i2","i2"),flag:",i1)') idimnetha,idimnetha+1,flag_girl
      do iSE=1,idimnetha
       do jSE=1,idimnetha+1
        write(*,'("b(",i2,",",i2,") :",d22.12)') iSE,jSE,b(iSE,jSE)
      enddo
     enddo
    endif
    rewind(io_runfile)
    write(io_runfile,*) nwmd,':girl crash in neth_alu with matrix b(21,22)'
    stop
  endif

  x(l)=c(1,1)
  if (l < llim) then
    if (x(l) <= 1.d-09 .and. idern == 1) then
      x(l)=0.d0
    endif
  else
    if (ns == nrband) then
      if (x(l) <= 1.d-09 .and. idern == 1) then
        x(l)=0.d0
      endif
    endif
  endif

  Y3(l)=c(2,1)*3.d0
  Y(l)=c(3,1)*4.d0
  xc12(l)=c(4,1)*12.d0
  xc13(l)=c(5,1)*13.d0
  xn14(l)=c(6,1)*14.d0
  xn15(l)=c(7,1)*15.d0
  xo16(l)=c(8,1)*16.d0
  xo17(l)=c(9,1)*17.d0
  xo18(l)=c(10,1)*18.d0
  xf19(l)=c(11,1)*19.d0
  xne20(l)=c(12,1)*20.d0
  xne21(l)=c(13,1)*21.d0
  xne22(l)=c(14,1)*22.d0
  xna23(l)=c(15,1)*23.d0
  xmg24(l)=c(16,1)*24.d0
  xmg25(l)=c(17,1)*25.d0
  xmg26(l)=c(18,1)*26.d0
  xal26(l)=c(19,1)*26.d0
  xal27(l)=c(20,1)*27.d0
  xsi28(l)=c(21,1)*28.d0

  if (ns == nrband) then
    if (l == m) then
      write(io_logs,'(1x,a,i5,1x,f10.7,1x,e8.2,1x,f10.7,12(1x,e8.2))') 'APRES NETH_ALU',l,x(l),y3(l),y(l),xc12(l),xc13(l), &
                 xn14(l),xn15(l),xo16(l),xo17(l),xo18(l),xne20(l),xne22(l),xmg24(l),xmg25(l),xmg26(l)
    endif
  endif

  return

end subroutine neth_alu
!======================================================================
subroutine nethe(l,ns,ddeit,flag_girl)
!-----------------------------------------------------------------------
  use SmallFunc,only: girl

  implicit none

  integer,intent(in):: l,ns
  integer,intent(inout):: flag_girl
  real(kindreal),intent(in):: ddeit

  integer,parameter:: idimnethe=12

  integer:: iSE,jSE
  real(kindreal):: d44,d4411,d124,d1242,d1241,d134,d1341,d1343,d164,d1641,d1644,d144,d1441,d1446, &
                   d184,d1841,d1847,d224n,d224g,d22ng,dweit

! vyab : Yi = Xi/Ai
!        Xi : fraction de masse de l'element i
!        Ai : masse atomique de l'element i
  real(kindreal),dimension(idimnethe):: vyab
  real(kindreal),dimension(idimnethe,idimnethe+1):: b
  real(kindreal),dimension(idimnethe,1):: d
!-----------------------------------------------------------------------
  d(:,:)=0.0d0
  b(:,:)=0.0d0
  vyab(:)=0.d0

  dweit=ddeit
  if (ns == 1) then
    vyab(1)=vvy(l)/4.d0
    vyab(2)=vvxc12(l)/12.d0
    vyab(3)=vvxc13(l)/13.d0
    vyab(4)=vvxo16(l)/16.d0
    vyab(5)=vvxne20(l)/20.d0
    vyab(6)=vvxn14(l)/14.d0
    vyab(7)=vvxo18(l)/18.d0
    vyab(8)=vvxne22(l)/22.d0
    vyab(9)=vvxmg25(l)/25.d0
    vyab(10)=vvxmg26(l)/26.d0
    vyab(11)=vvxo17(l)/17.d0
    vyab(12)=vvxmg24(l)/24.d0
    if (ns == nrband) then
      if (l == m) then
        write(io_logs,'(1x,a,i4,10(1x,f10.7))') 'AVANT NETHE',l,vy(l),vxc12(l),vxc13(l),vxn14(l),vxo16(l),vxo17(l), &
                   vxo18(l),vxne20(l),vxne22(l),vxmg24(l)
      endif
    endif
  else
    vyab(1)=y(l)/4.d0
    vyab(2)=xc12(l)/12.d0
    vyab(3)=xc13(l)/13.d0
    vyab(4)=xo16(l)/16.d0
    vyab(5)=xne20(l)/20.d0
    vyab(6)=xn14(l)/14.d0
    vyab(7)=xo18(l)/18.d0
    vyab(8)=xne22(l)/22.d0
    vyab(9)=xmg25(l)/25.d0
    vyab(10)=xmg26(l)/26.d0
    vyab(11)=xo17(l)/17.d0
    vyab(12)=xmg24(l)/24.d0
    if (ns == nrband) then
      if (l == m) then
        write(io_logs,'(1x,a,i4,10(1x,f10.7))') 'AVANT NETHE',l,vy(l),vxc12(l),vxc13(l),vxn14(l),vxo16(l),vxo17(l), &
                   vxo18(l),vxne20(l),vxne22(l),vxmg24(l)
      endif
    endif
  endif
  d44=epsyy(l)*dweit
  d4411=vyab(1)*vyab(1)*d44
  d124=epsyc(l)*dweit
  d1242=d124*vyab(2)
  d1241=d124*vyab(1)
  d134=c134(l)*dweit
  d1341=d134*vyab(1)
  d1343=d134*vyab(3)
  d164=epsyo(l)*dweit
  d1641=d164*vyab(1)
  d1644=d164*vyab(4)
  d144=c144(l)*dweit
  d1441=d144*vyab(1)
  d1446=d144*vyab(6)
  d184=c184(l)*dweit
  d1841=d184*vyab(1)
  d1847=d184*vyab(7)
  d224n=c224(l)*dweit
  d224g=c224g(l)*dweit
  d22ng=d224n+d224g

! Initialisation a zero de la matrice du reseau de reactions nucleaires
! pour chaque coquille.
  b=0.d0

! Ecriture de la matrice b(y,z) du reseau de reactions nucleaires.
  b(3,3)=1.d0+d1341
  b(3,1)=d1343
  b(3,13)=vyab(3)*b(3,3)
  b(4,4)=1.d0+d1641
  b(4,1)=d1644-d1242-d1343
  b(4,2)=-d1241
  b(4,3)=-d1341
  b(4,13)=vyab(4)*b(4,4)-d1241*vyab(2)-d1341*vyab(3)
  b(5,5)=1.d0+vyab(1)*e20ag(l)*dweit
  b(5,4)=-d1641
  b(5,1)=-d1644-vyab(11)*dweit*e17an(l)+vyab(5)*e20ag(l)*dweit
  b(5,11)=-vyab(1)*dweit*e17an(l)
  b(5,13)=vyab(5)*b(5,5)-vyab(1)*vyab(11)*dweit*e17an(l)-d1641*vyab(4)
  b(6,6)=1.d0+d1441
  b(6,1)=d1446
  b(6,13)=vyab(6)*b(6,6)
  b(7,7)=1.d0+d1841
  b(7,1)=d1847-d1446
  b(7,6)=-d1441
  b(7,13)=vyab(7)*b(7,7)-d1441*vyab(6)
  b(8,8)=1.d0+vyab(1)*d22ng
  b(8,1)=vyab(8)*d22ng-d1847
  b(8,7)=-d1841
  b(8,13)=vyab(8)*b(8,8)-d1841*vyab(7)
  b(9,9)=1.d0
  b(9,8)=-vyab(1)*d224n
  b(9,1)=-vyab(8)*d224n
  b(9,13)=vyab(9)+vyab(8)*b(9,8)
  b(10,10)=1.d0
  b(10,8)=-vyab(1)*d224g
  b(10,1)=-vyab(8)*d224g
  b(10,13)=vyab(10)+vyab(8)*b(10,8)
  b(11,11)=1.d0+vyab(1)*dweit*e17an(l)
  b(11,1)=vyab(11)*dweit*e17an(l)
  b(11,13)=vyab(11)*(1.d0+vyab(1)*dweit*e17an(l))
  b(12,1)=-vyab(5)*e20ag(l)*dweit
  b(12,5)=-vyab(1)*e20ag(l)*dweit
  b(12,12)=1.d0
  b(12,13)=vyab(12)+vyab(5)*b(12,5)
  b(2,2)=1.d0+d1241
  b(2,1)=d1242-0.5d0*d4411
  b(2,13)=vyab(2)*b(2,2)-(1.d0/3.d0)*vyab(1)*d4411
  b(1,1)=1.d0+1.5d0*d4411+d1242+d1343+d1644+d1446+d1847+d22ng*vyab(8)+vyab(11)*dweit*e17an(l)+vyab(5)*e20ag(l)*dweit
  b(1,2)=d1241
  b(1,3)=d1341
  b(1,4)=d1641
  b(1,5)=vyab(1)*e20ag(l)*dweit
  b(1,6)=d1441
  b(1,7)=d1841
  b(1,8)=d22ng*vyab(1)
  b(1,11)=vyab(1)*dweit*e17an(l)
  b(1,13)=vyab(1)*(1.d0+d1242+d1343+d1644+d1446+d1847+d22ng*vyab(8))+d4411*vyab(1)+vyab(1)*vyab(11)*dweit*e17an(l)+ &
          vyab(1)*vyab(5)*e20ag(l)*dweit
! Inversion de la matrice b du reseau de reactions nucleaires.

! Multiplication par le R.H.S. des equations aux differences, stocke
! dans les colonnes (z > 11) de la matrice b.
! Apres multiplication, la premiere colonne du vecteur d contient les
! abondances au temps t(n+1).
  flag_girl = 0
  call girl(b,d,idimnethe,1,flag_girl)
  if (flag_girl /= 0) then
    if (idebug>0) then
      write(*,'("nethe - matrix b(",i2","i2"),flag:",i1)') idimnethe,idimnethe+1,flag_girl
      do iSE=1,idimnethe
       do jSE=1,idimnethe+1
        write(*,'("b(",i2,",",i2,") :",d22.12)') iSE,jSE,b(iSE,jSE)
       enddo
      enddo
    endif
    rewind(io_runfile)
    write(io_runfile,*) nwmd,':girl crashes in nethe with matrix b(12,13)'
    stop
  endif

! Nouvelles abondances
  y(l)=d(1,1)*4.d0
  xc12(l)=12.d0*d(2,1)
  xc13(l)=13.d0*d(3,1)
  xo16(l)=16.d0*d(4,1)
  xne20(l)=20.d0*d(5,1)
  xn14(l)=14.d0*d(6,1)
  xo18(l)=18.d0*d(7,1)
  xne22(l)=22.d0*d(8,1)
  xmg25(l)=25.d0*d(9,1)
  xmg26(l)=26.d0*d(10,1)
  xo17(l)=17.d0*d(11,1)
  xmg24(l)=24.d0*d(12,1)
  if (ns == nrband) then
    if (l == m) then
      write(io_logs,'(1x,a,i4,10(1x,f10.7))') 'APRES NETHE',l,y(l),xc12(l),xc13(l),xn14(l),xo16(l),xo17(l),xo18(l), &
                 xne20(l),xne22(l),xmg24(l)
    endif

    if (xc13(l) < 1.0d-75) then
      xc13(l)=0.d0
    endif
    if (xn14(l) < 1.0d-75) then
      xn14(l)=0.d0
    endif
    if (xo17(l) < 1.0d-75) then
      xo17(l)=0.d0
    endif
    if (xo18(l) < 1.0d-75) then
      xo18(l)=0.d0
    endif
  endif

  return

end subroutine nethe
!======================================================================
subroutine nethe_alu(l,ns,ddeit,flag_girl)
!-----------------------------------------------------------------------
  use SmallFunc,only: girl

  implicit none

  integer,intent(in):: l,ns
  integer,intent(inout):: flag_girl
  real(kindreal),intent(in):: ddeit

  integer,parameter:: idimnethea=24

  integer:: iSE,jSE
  real(kindreal):: dweit,d44,d4411,d124,d1242,d1241,d134,d1341,d1343,d164,d1641,d1644,d144,d1441,d1446,&
    d184,d1841,d1847,d224n,d224g,d22ng,uno,d14np1,d14np2,d14npl,d14ng1,d14ng2,d14ngl,d12ng1,d12ng2,d12ngl,&
    d19ng1,d19ng2,d19ngl,d12pg1,d12pg2,d12pgl,d13pg1,d13pg2,d13pgl,dn14p1,dn14p2,dn14pl,d15pg1,d15pg2,d15pgl,&
    d15pa1,d15pa2,d15pal,d16pg1,d16pg2,d16pgl,d17pa1,d17pa2,d17pal,d17pg1,d17pg2,d17pgl,d18pg1,d18pg2,d18pgl,&
    d14pg1,d14pg2,d14pgl,d14ag1,d14ag2,d14agl,d18an1,d18an2,d18anl,d18na1,d18na2,d18nal,d15ag1,d15ag2,d15agl,&
    d18np1,d18np2,d18npl,d18pa1,d18pa2,d18pal,d19ap1,d19ap2,d19apl,d20ng1,d20ng2,d20ngl,d21ng1,d21ng2,d21ngl,&
    d22ng1,d22ng2,d22ngl,d23ng1,d23ng2,d23ngl,d24ng1,d24ng2,d24ngl,d25ng1,d25ng2,d25ngl,d26ng1,d26ng2,d26ngl,&
    d18ng1,d18ng2,d18ngl,dc14n1,dc14n2,dc14nl,d24ag1,d24ag2,d24agl,d17ag1,d17ag2,d17agl,d21ag1,d21ag2,d21agl,&
    d21na1,d21na2,d21nal,d25an1,d25an2,d25anl,d27ng1,d27ng2,d27ngl,d28ng1,d28ng2,d28ngl,da26a1,da26a2,da26al,&
    da26g1,da26g2,da26gl,dc14be,df18be,da26be,b55,b66,b77,b88

! vyab : Yi = Xi/Ai
!        Xi : fraction de masse de l'element i
!        Ai : masse atomique de l'element i
  real(kindreal), dimension(idimnethea):: vyab
  real(kindreal), dimension(600):: bb
  real(kindreal), dimension(idimnethea,idimnethea+1):: b,b_before
  real(kindreal), dimension(idimnethea,1):: d

  equivalence(b,bb)
!-----------------------------------------------------------------------
  flag_girl = 0
  dweit=ddeit
  b(:,:) = 0.d0
  if (ns == 1) then
    vyab(1)=vvy(l)/4.d0
    vyab(2)=vvxc12(l)/12.d0
    vyab(3)=vvxc13(l)/13.d0
    vyab(4)=vvxo16(l)/16.d0
    vyab(5)=vvxne20(l)/20.d0
    vyab(6)=vvxn14(l)/14.d0
    vyab(7)=vvxo18(l)/18.d0
    vyab(8)=vvxne22(l)/22.d0
    vyab(9)=vvxmg25(l)/25.d0
    vyab(10)=vvxmg26(l)/26.d0
    vyab(11)=vvxo17(l)/17.d0
    vyab(12)=vvxmg24(l)/24.d0
    vyab(13)=vvxf18(l)/18.d0
    vyab(14)=vvxc14(l)/14.d0
    vyab(15)=vvxneut(l)
    vyab(16)=vvxprot(l)
    vyab(17)=vvxn15(l)/15.d0
    vyab(18)=vvxne21(l)/21.d0
    vyab(19)=vvxf19(l)/19.d0
    vyab(20)=vvxna23(l)/23.d0
    vyab(21)=vvxal27(l)/27.d0
    vyab(22)=vvxal26g(l)/26.d0
    vyab(23)=vvxbid(l)/40.d0
    vyab(24)=vvxbid1(l)/41.d0
    if (ns == nrband) then
      if (l == m) then
        write(io_logs,&
                '(1x,a,i4,1x,f10.7,11(1x,e8.2),/,11x,3(1x,e8.2),9(1x,f10.7),/,11x,1(1x,f10.7))') 'BEFORE NETHE_ALU',l,vvy(l), &
                vvxc12(l),vvxc13(l),vvxn14(l),vvxn15(l),vvxo16(l),vvxo17(l),vvxo18(l),vvxne20(l),vvxne22(l),vvxmg24(l), &
                vvxmg25(l),vvxmg26(l),vvxf18(l),vvxc14(l),vvxneut(l),vvxprot(l),vvxn15(l),vvxne21(l),vvxf19(l),vvxna23(l), &
                vvxal27(l),vvxal26g(l),vvxbid(l),vvxbid1(l)
      endif
    endif
  else
    vyab(1)=y(l)/4.d0
    vyab(2)=xc12(l)/12.d0
    vyab(3)=xc13(l)/13.d0
    vyab(4)=xo16(l)/16.d0
    vyab(5)=xne20(l)/20.d0
    vyab(6)=xn14(l)/14.d0
    vyab(7)=xo18(l)/18.d0
    vyab(8)=xne22(l)/22.d0
    vyab(9)=xmg25(l)/25.d0
    vyab(10)=xmg26(l)/26.d0
    vyab(11)=xo17(l)/17.d0
    vyab(12)=xmg24(l)/24.d0
    vyab(13)=xf18(l)/18.d0
    vyab(14)=xc14(l)/14.d0
    vyab(15)=xneut(l)
    vyab(16)=xprot(l)
    vyab(17)=xn15(l)/15.d0
    vyab(18)=xne21(l)/21.d0
    vyab(19)=xf19(l)/19.d0
    vyab(20)=xna23(l)/23.d0
    vyab(21)=xal27(l)/27.d0
    vyab(22)=xal26(l)/26.d0
    vyab(23)=xbid(l)/40.d0
    vyab(24)=xbid1(l)/41.d0
    if (ns == nrband) then
      if (l == m) then
        write(io_logs,'(1x,a,i4,1x,f10.7,11(1x,e8.2),/,11x,3(1x,e8.2),9(1x,f10.7),/,11x,1(1x,f10.7))') 'BEFORE NETHE_ALU',l,y(l), &
                   xc12(l),xc13(l),xn14(l),xn15(l),xo16(l),xo17(l),xo18(l),xne20(l),xne22(l),xmg24(l),xmg25(l),xmg26(l), &
                   xf18(l),xc14(l),xneut(l),xprot(l),xn15(l),xne21(l),xf19(l),xna23(l),xal27(l),xal26(l),xbid(l),xbid1(l)
      endif
    endif
  endif
  d44=epsyy(l)*dweit
  d4411=vyab(1)*vyab(1)*d44
  d124=epsyc(l)*dweit
  d1242=d124*vyab(2)
  d1241=d124*vyab(1)
  d134=c134(l)*dweit
  d1341=d134*vyab(1)
  d1343=d134*vyab(3)
  d164=epsyo(l)*dweit
  d1641=d164*vyab(1)
  d1644=d164*vyab(4)
  d144=c144(l)*dweit
  d1441=d144*vyab(1)
  d1446=d144*vyab(6)
  d184=c184(l)*dweit
  d1841=d184*vyab(1)
  d1847=d184*vyab(7)
  d224n=c224(l)*dweit
  d224g=c224g(l)*dweit
  d22ng=d224n+d224g

  uno=e14np(l)*dweit
  d14np1=vyab(15)*uno
  d14np2=vyab(6)*uno
  d14npl=vyab(15)*d14np2

  uno=e14ng(l)*dweit
  d14ng1=vyab(15)*uno
  d14ng2=vyab(6)*uno
  d14ngl=vyab(15)*d14ng2

  uno=e12ng(l)*dweit
  d12ng1=vyab(15)*uno
  d12ng2=vyab(2)*uno
  d12ngl=vyab(15)*d12ng2

  uno=e19ng(l)*dweit
  d19ng1=vyab(15)*uno
  d19ng2=vyab(19)*uno
  d19ngl=vyab(15)*d19ng2

  uno=b112(l)*dweit
  d12pg1=vyab(16)*uno
  d12pg2=vyab(2)*uno
  d12pgl=vyab(16)*d12pg2

  uno=b113(l)*dweit
  d13pg1=vyab(16)*uno
  d13pg2=vyab(3)*uno
  d13pgl=vyab(16)*d13pg2

  uno=b114(l)*dweit
  dn14p1=vyab(16)*uno
  dn14p2=vyab(6)*uno
  dn14pl=vyab(16)*dn14p2

  uno=b115g(l)*dweit
  d15pg1=vyab(16)*uno
  d15pg2=vyab(17)*uno
  d15pgl=vyab(16)*d15pg2

  uno=b115a(l)*dweit
  d15pa1=vyab(16)*uno
  d15pa2=vyab(17)*uno
  d15pal=vyab(16)*d15pa2

  uno=b116(l)*dweit
  d16pg1=vyab(16)*uno
  d16pg2=vyab(4)*uno
  d16pgl=vyab(16)*d16pg2

  uno=b117a(l)*dweit
  d17pa1=vyab(16)*uno
  d17pa2=vyab(11)*uno
  d17pal=vyab(16)*d17pa2

  uno=b117g(l)*dweit
  d17pg1=vyab(16)*uno
  d17pg2=vyab(11)*uno
  d17pgl=vyab(16)*d17pg2

  uno=b118g(l)*dweit
  d18pg1=vyab(16)*uno
  d18pg2=vyab(7)*uno
  d18pgl=vyab(16)*d18pg2

  uno=ec14pg(l)*dweit
  d14pg1=vyab(16)*uno
  d14pg2=vyab(14)*uno
  d14pgl=vyab(16)*d14pg2

  uno=ec14ag(l)*dweit
  d14ag1=vyab(1)*uno
  d14ag2=vyab(14)*uno
  d14agl=vyab(1)*d14ag2

  uno=e18an(l)*dweit
  d18an1=vyab(1)*uno
  d18an2=vyab(7)*uno
  d18anl=vyab(1)*d18an2

  uno=ef18na(l)*dweit
  d18na1=vyab(15)*uno
  d18na2=vyab(13)*uno
  d18nal=vyab(15)*d18na2

  uno=e15ag(l)*dweit
  d15ag1=vyab(1)*uno
  d15ag2=vyab(17)*uno
  d15agl=vyab(1)*d15ag2

  uno=ef18np(l)*dweit
  d18np1=vyab(15)*uno
  d18np2=vyab(13)*uno
  d18npl=vyab(15)*d18np2

  uno=e18pa(l)*dweit
  d18pa1=vyab(16)*uno
  d18pa2=vyab(7)*uno
  d18pal=vyab(16)*d18pa2

  uno=e19ap(l)*dweit
  d19ap1=vyab(1)*uno
  d19ap2=vyab(19)*uno
  d19apl=vyab(1)*d19ap2

  uno=e20ng(l)*dweit
  d20ng1=vyab(15)*uno
  d20ng2=vyab(5)*uno
  d20ngl=vyab(15)*d20ng2

  uno=e21ng(l)*dweit
  d21ng1=vyab(15)*uno
  d21ng2=vyab(18)*uno
  d21ngl=vyab(15)*d21ng2

  uno=e22ng(l)*dweit
  d22ng1=vyab(15)*uno
  d22ng2=vyab(8)*uno
  d22ngl=vyab(15)*d22ng2

  uno=e23ng(l)*dweit
  d23ng1=vyab(15)*uno
  d23ng2=vyab(20)*uno
  d23ngl=vyab(15)*d23ng2

  uno=e24ng(l)*dweit
  d24ng1=vyab(15)*uno
  d24ng2=vyab(12)*uno
  d24ngl=vyab(15)*d24ng2

  uno=e25ng(l)*dweit
  d25ng1=vyab(15)*uno
  d25ng2=vyab(9)*uno
  d25ngl=vyab(15)*d25ng2

  uno=e26ng(l)*dweit
  d26ng1=vyab(15)*uno
  d26ng2=vyab(10)*uno
  d26ngl=vyab(15)*d26ng2

  uno=e18ng(l)*dweit
  d18ng1=vyab(15)*uno
  d18ng2=vyab(7)*uno
  d18ngl=vyab(15)*d18ng2

  uno=ec14ng(l)*dweit
  dc14n1=vyab(15)*uno
  dc14n2=vyab(14)*uno
  dc14nl=vyab(15)*dc14n2

  uno=e24ag(l)*dweit
  d24ag1=vyab(1)*uno
  d24ag2=vyab(12)*uno
  d24agl=vyab(1)*d24ag2

  uno=e17ag(l)*dweit
  d17ag1=vyab(1)*uno
  d17ag2=vyab(11)*uno
  d17agl=vyab(1)*d17ag2

  uno=e21ag(l)*dweit
  d21ag1=vyab(1)*uno
  d21ag2=vyab(18)*uno
  d21agl=vyab(1)*d21ag2

  uno=e21na(l)*dweit
  d21na1=vyab(15)*uno
  d21na2=vyab(18)*uno
  d21nal=vyab(15)*d21na2

  uno=e25an(l)*dweit
  d25an1=vyab(1)*uno
  d25an2=vyab(9)*uno
  d25anl=vyab(1)*d25an2

  uno=e27ng(l)*dweit
  d27ng1=vyab(15)*uno
  d27ng2=vyab(21)*uno
  d27ngl=vyab(15)*d27ng2

  uno=e28ng(l)*dweit
  d28ng1=vyab(15)*uno
  d28ng2=vyab(23)*uno
  d28ngl=vyab(15)*d28ng2

  uno=a26ga(l)*dweit
  da26a1=vyab(15)*uno
  da26a2=vyab(22)*uno
  da26al=vyab(15)*da26a2

  uno=a26gp(l)*dweit
  da26g1=vyab(15)*uno
  da26g2=vyab(22)*uno
  da26gl=vyab(15)*da26g2
  dc14be=e14be(l)*dweit
  df18be=e18be(l)*dweit
  da26be=e26be(l)*dweit

! initialisation
  bb(:)=0.d0

  b(3,3)=1.d0+d1341+d13pg1
  b(3,2)=-d12pg1-d12ng1
  b(3,15)=-d12ng2
  b(3,16)=-d12pg2+d13pg2
  b(3,1)=d1343
  b(3,25)=vyab(3)*b(3,3)-d12pgl-d12ngl

  b(4,4)=1.d0+d1641+d16pg1
  b(4,1)=d1644-d1242-d1343
  b(4,2)=-d1241
  b(4,3)=-d1341
  b(4,16)=-d15pg2+d16pg2
  b(4,17)=-d15pg1
  b(4,25)=vyab(4)*b(4,4)-d1241*vyab(2)-d1341*vyab(3)-d15pgl

  b55=1.d0+vyab(1)*e20ag(l)*dweit
  b(5,5)=b55+d20ng1
  b(5,4)=-d1641
  b(5,1)=-d1644-vyab(11)*dweit*e17an(l)+vyab(5)*e20ag(l)*dweit
  b(5,11)=-vyab(1)*dweit*e17an(l)
  b(5,15)=+d20ng2-d19ng2
  b(5,19)=-d19ng1
  b(5,25)=vyab(5)*b55-vyab(1)*vyab(11)*dweit*e17an(l)-d1641*vyab(4)+d20ngl-d19ngl

  b66=1.d0+d1441
  b(6,6)=b66+d14np1+dn14p1+d14ng1
  b(6,1)=d1446
  b(6,3)=-d13pg1
  b(6,11)=-d17pa1
  b(6,14)=-dc14be
  b(6,15)=+d14np2+d14ng2
  b(6,16)=-d13pg2+dn14p2-d17pa2
  b(6,25)=vyab(6)*b66-d13pgl+d14npl+dn14pl-d17pal+d14ngl

  b77=1.d0+d1841
  b(7,7)=b77+d18an1+d18pa1+d18ng1+d18pg1
  b(7,1)=d1847-d14ag2+d18an2
  b(7,13)=-d18np1-df18be
  b(7,14)=-d14ag1
  b(7,15)=-d18np2-d21na2+d18ng2
  b(7,16)=+d18pa2+d18pg2
  b(7,18)=-d21na1
  b(7,25)=vyab(7)*b77-d14agl+d18anl-d18npl+d18pal-d21nal+d18ngl+d18pgl

  b88=1.d0+vyab(1)*d22ng
  b(8,8)=b88+d22ng1
  b(8,1)=vyab(8)*d22ng-d1847-d19ap2
  b(8,7)=-d1841
  b(8,15)=-d21ng2+d22ng2
  b(8,18)=-d21ng1
  b(8,19)=-d19ap1
  b(8,25)=vyab(8)*b88-d1841*vyab(7)-d19apl-d21ngl+d22ngl

  b(9,9)=1.d0+d25ng1+d25an1
  b(9,8)=-vyab(1)*d224n
  b(9,1)=-vyab(8)*d224n-d21ag2+d25an2
  b(9,12)=-d24ng1
  b(9,15)=-d24ng2+d25ng2
  b(9,18)=-d21ag1
  b(9,25)=vyab(9)+vyab(8)*b(9,8)-d24ngl+d25ngl-d21agl+d25anl

  b(10,10)=1.d0+d26ng1
  b(10,8)=-vyab(1)*d224g
  b(10,9)=-d25ng1
  b(10,1)=-vyab(8)*d224g
  b(10,15)=+d26ng2-d25ng2-da26g2
  b(10,22)=-da26g1-da26be
  b(10,25)=vyab(10)+vyab(8)*b(10,8)+d26ngl-d25ngl-da26gl

  b(11,11)=1.d0+vyab(1)*dweit*e17an(l)+d17ag1+d17pa1+d17pg1
  b(11,1)=vyab(11)*dweit*e17an(l)+d17ag2
  b(11,4)=-d16pg1
  b(11,16)=-d16pg2+d17pa2+d17pg2
  b(11,25)=vyab(11)*(1.d0+vyab(1)*dweit*e17an(l))+d17agl-d16pgl+d17pal+d17pgl

  b(12,1)=-vyab(5)*e20ag(l)*dweit+d24ag2
  b(12,5)=-vyab(1)*e20ag(l)*dweit
  b(12,12)=1.d0+d24ng1+d24ag1
  b(12,15)=-d23ng2+d24ng2
  b(12,20)=-d23ng1
  b(12,25)=vyab(12)+vyab(5)*b(12,5)-d23ngl+d24ngl+d24agl

  b(2,2)=1.d0+d1241+d12pg1+d12ng1
  b(2,15)=+d12ng2
  b(2,16)=+d12pg2-d15pa2
  b(2,17)=-d15pa1
  b(2,1)=d1242-0.5d0*d4411
  b(2,25)=vyab(2)*b(2,2)-(1.d0/3.d0)*vyab(1)*d4411-d15pal

  b(1,1)=1.d0+1.5d0*d4411+d1242+d1343+d1644+d1446+d1847+d22ng*vyab(8)+vyab(11)*dweit*e17an(l) &
         +vyab(5)*e20ag(l)*dweit+d14ag2+d18an2+d15ag2+d19ap2+d24ag2+d17ag2+d21ag2+d25an2
  b(1,2)=d1241
  b(1,3)=d1341
  b(1,4)=d1641
  b(1,5)=vyab(1)*e20ag(l)*dweit
  b(1,6)=d1441
  b(1,7)=d1841+d18an1-d18pa1
  b(1,8)=d22ng*vyab(1)
  b(1,9)=+d25an1
  b(1,11)=vyab(1)*dweit*e17an(l)+d17ag1-d17pa1
  b(1,12)=+d24ag1
  b(1,13)=-d18na1
  b(1,14)=+d14ag1
  b(1,15)=-d18na2-d21na2-da26a2
  b(1,16)=-d18pa2-d15pa2-d17pa2
  b(1,17)=+d15ag1-d15pa1
  b(1,18)=+d21ag1-d21na1
  b(1,19)=+d19ap1
  b(1,22)=-da26a1
  b(1,25)=vyab(1)*(1.d0+d1242+d1343+d1644+d1446+d1847+d22ng*vyab(8))+d4411*vyab(1)+ &
          vyab(1)*vyab(11)*dweit*e17an(l)+vyab(1)*vyab(5)*e20ag(l)*dweit+d14agl+d18anl- &
          d18nal-d18pal+d15agl+d19apl+d24agl+d17agl+d21agl-d21nal+d25anl-da26al-d15pal-d17pal

  b(13,1)=-d1446
  b(13,6)=-d1441
  b(13,11)=-d17pg1
  b(13,13)=+d18na1+d18np1+df18be
  b(13,15)=+d18na2+d18np2
  b(13,16)=-d17pg2
  b(13,25)=+d18nal+d18npl-vyab(1)*d1446-d17pgl

  b(14,1)=+d14ag2
  b(14,6)=-d14np1
  b(14,14)=+d14pg1+d14ag1+dc14n1+dc14be+1.d0
  b(14,15)=-d14np2+dc14n2
  b(14,16)=+d14pg2
  b(14,25)=-d14npl+d14pgl+d14agl+dc14nl+vyab(14)

  b(15,1)=-d18an2-d25an2-d1343-e17an(l)*dweit*vyab(11)-d224n*vyab(8)
  b(15,2)=+d12ng1
  b(15,3)=-d1341
  b(15,5)=+d20ng1
  b(15,6)=+d14np1+d14ng1
  b(15,7)=-d18an1+d18ng1
  b(15,8)=+d22ng1-d224n*vyab(1)
  b(15,9)=+d25ng1-d25an1
  b(15,10)=+d26ng1
  b(15,11)=-vyab(1)*e17an(l)*dweit
  b(15,12)=+d24ng1
  b(15,13)=+d18na1+d18np1
  b(15,14)=+dc14n1
  b(15,15)=+d14np2+d18na2+d18np2+d20ng2+d21ng2+d22ng2+d23ng2+d24ng2+d25ng2+d26ng2+ &
            d18ng2+dc14n2+d21na2+d27ng2+d28ng2+da26a2+da26g2+d12ng2+d14ng2+d19ng2
  b(15,18)=+d21ng1+d21na1
  b(15,19)=+d19ng1
  b(15,20)=+d23ng1
  b(15,21)=+d27ng1
  b(15,22)=+da26a1+da26g1
  b(15,23)=+d28ng1
  b(15,25)=+d14npl-d18anl+d18nal+d18npl+d20ngl+d21ngl+d22ngl+d23ngl+d24ngl+d25ngl+ &
            d26ngl+d18ngl+dc14nl+d21nal-d25anl+d27ngl+d28ngl+da26al+da26gl-vyab(1)*d1343-&
            e17an(l)*dweit*vyab(11)*vyab(1)-d224n*vyab(8)*vyab(1)+d12ngl+d14ngl+d19ngl

  b(16,1)=-d19ap2
  b(16,2)=+d12pg1
  b(16,3)=+d13pg1
  b(16,4)=+d16pg1
  b(16,6)=-d14np1+dn14p1
  b(16,7)=+d18pa1+d18pg1
  b(16,11)=+d17pa1+d17pg1
  b(16,13)=-d18np1
  b(16,14)=+d14pg1
  b(16,15)=-d14np2-d18np2-da26g2
  b(16,16)=+d14pg2+d18pa2+1.d0+d12pg2+d13pg2+dn14p2+d15pg2+d15pa2+d16pg2+d17pa2+d17pg2+d18pg2
  b(16,17)=+d15pg1+d15pa1
  b(16,19)=-d19ap1
  b(16,22)=-da26g1
  b(16,25)=-d14npl+d14pgl-d18npl+d18pal-d19apl-da26gl+vyab(16)+d12pgl+d13pgl+dn14pl+&
           d15pgl+d15pal+d16pgl+d17pal+d17pgl+d18pgl

  b(17,1)=+d15ag2
  b(17,6)=-dn14p1-d14ng1
  b(17,7)=-d18pa1
  b(17,13)=-d18na1
  b(17,14)=-d14pg1-dc14n1
  b(17,15)=-d18na2-dc14n2-d14ng2
  b(17,16)=-d14pg2-d18pa2-dn14p2+d15pg2+d15pa2
  b(17,17)=+d15ag1+1.d0+d15pg1+d15pa1
  b(17,25)=-d14pgl-d18nal-d18pal+d15agl-dc14nl+vyab(17)-dn14pl+d15pgl+d15pal-d14ngl

  b(18,1)=-d18an2-d17ag2+d21ag2
  b(18,5)=-d20ng1
  b(18,7)=-d18an1
  b(18,11)=-d17ag1
  b(18,15)=-d20ng2+d21ng2+d21na2
  b(18,18)=+d21ng1+d21ag1+d21na1+1.d0
  b(18,25)=-d18anl-d20ngl+d21ngl-d17agl+d21agl+d21nal+vyab(18)

  b(19,1)=-d15ag2+d19ap2
  b(19,7)=-d18ng1-d18pg1
  b(19,15)=-d18ng2+d19ng2
  b(19,16)=-d18pg2
  b(19,17)=-d15ag1
  b(19,19)=+d19ap1+1.d0+d19ng1
  b(19,25)=-d15agl+d19apl-d18ngl+vyab(19)-d18pgl+d19ngl

  b(20,8)=-d22ng1
  b(20,15)=-d22ng2+d23ng2-da26a2
  b(20,20)=+d23ng1+1.d0
  b(20,22)=-da26a1
  b(20,25)=-d22ngl+d23ngl-da26al+vyab(20)

  b(21,10)=-d26ng1
  b(21,15)=-d26ng2+d27ng2
  b(21,21)=+d27ng1+1.d0
  b(21,25)=-d26ngl+d27ngl+vyab(21)

  b(22,15)=+da26a2+da26g2
  b(22,22)=+da26a1+da26g1+da26be+1.d0
  b(22,25)=+da26al+da26gl+vyab(22)

  b(23,1)=-d24ag2-d25an2
  b(23,9)=-d25an1
  b(23,12)=-d24ag1
  b(23,15)=-d27ng2+d28ng2
  b(23,21)=-d27ng1
  b(23,23)=+d28ng1+1.d0
  b(23,25)=-d24agl-d25anl-d27ngl+d28ngl+vyab(23)

  b(24,15)=-d28ng2
  b(24,23)=-d28ng1
  b(24,24)=1.d0
  b(24,25)=-d28ngl+vyab(24)
  b_before(:,:) = b(:,:)

  flag_girl = 0
  call girl(b,d,idimnethea,1,flag_girl)
  if (flag_girl /= 0) then
    if (idebug>0) then
      write(*,'("nethe_alu, layer ",i4," - matrix b(",i2","i2"),flag:",i1)') &
            l,idimnethea,idimnethea+1,flag_girl
      write(*,*) 'x(l),t(l),n:',x(l),t(l),vyab(15)
      do iSE=1,idimnethea
       do jSE=1,idimnethea+1
        write(*,'("b(",i2,",",i2,") :",d22.12)') iSE,jSE,b_before(iSE,jSE)
      enddo
     enddo
    endif
    rewind(io_runfile)
    write(io_runfile,*) nwmd,':girl crash in nethe_alu with matrix b(24,25)'
    write(*,*) nwmd,':girl crash in nethe_alu with matrix b(24,25)'
    stop
  endif

  y(l)=d(1,1)*4.d0
  xc12(l)=12.d0*d(2,1)
  xc13(l)=13.d0*d(3,1)
  xo16(l)=16.d0*d(4,1)
  xne20(l)=20.d0*d(5,1)
  xn14(l)=14.d0*d(6,1)
  xo18(l)=18.d0*d(7,1)
  xne22(l)=22.d0*d(8,1)
  xmg25(l)=25.d0*d(9,1)
  xmg26(l)=26.d0*d(10,1)
  xo17(l)=17.d0*d(11,1)
  xmg24(l)=24.d0*d(12,1)
  xf18(l)=18.d0*d(13,1)
  xc14(l)=14.d0*d(14,1)
  xneut(l)=d(15,1)
  xprot(l)=d(16,1)
  xn15(l)=d(17,1)*15.d0
  xne21(l)=21.d0*d(18,1)
  xf19(l)=19.d0*d(19,1)
  xna23(l)=23.d0*d(20,1)
  xal27(l)=27.d0*d(21,1)
  xal26(l)=26.d0*d(22,1)
  xbid(l)=40.d0*d(23,1)
  xbid1(l)=41.d0*d(24,1)
  xbid(l)=xbid(l)+xbid1(l)
  xbid1(l)=0.d0

  if (ns == nrband) then
    if (l == m) then
      write(io_logs,'(1x,a,i4,1x,f10.7,11(1x,e8.2),/,11x,3(1x,e8.2),9(1x,f10.7),/,11X,1(1x,f10.7))') 'AFTER NETHE_ALU',l,y(l), &
                 xc12(l),xc13(l),xn14(l),xn15(l),xo16(l),xo17(l),xo18(l),xne20(l),xne22(l),xmg24(l),xmg25(l),xmg26(l), &
                 xf18(l),xc14(l),xneut(l),xprot(l),xn15(l),xne21(l),xf19(l),xna23(l),xal27(l),xal26(l),xbid(l),xbid1(l)
    endif
    if (xprot(l) < 0.0d0) then
      xprot(l)=0.d0
    endif
    if (xneut(l) < 0.0d0) then
      xneut(l)=0.d0
    endif
    if (xal26(l) < 1.0d-75) then
      xal26(l)=0.d0
    endif
  endif

  return

end subroutine nethe_alu
!======================================================================
subroutine netc(l,ddeit)
!-----------------------------------------------------------------------
  implicit none

  integer,intent(in):: l
  real(kindreal),intent(in):: ddeit

  integer,parameter:: idimnetc=15

  integer:: i,ii
  real(kindreal):: sumvxab,t9

  real(kindreal),dimension(idimnetc):: vxab

!-----------------------------------------------------------------------
  vxab(1)  = vvx(l)
  vxab(2)  = vvy3(l)
  vxab(3)  = vvy(l)
  vxab(4)  = vvxc12(l)
  vxab(5)  = vvxc13(l)
  vxab(6)  = vvxn14(l)
  vxab(7)  = vvxn15(l)
  vxab(8)  = vvxo16(l)
  vxab(9)  = vvxo17(l)
  vxab(10) = vvxo18(l)
  vxab(11) = vvxne20(l)
  vxab(12) = vvxne22(l)
  vxab(13) = vvxmg24(l)
  vxab(14) = vvxmg25(l)
  vxab(15) = vvxmg26(l)

  sumvxab=1.d0-zabelx
  do i=1,15
   sumvxab=sumvxab-vxab(i)
  enddo
  do ii=1,nbelx
   sumvxab=sumvxab-vvabelx(ii,l)
  enddo
  if (abs(sumvxab) > 1.d-2 .and. verbose) then
    print*, l,'sumvxab= ', sumvxab
  endif

  if (l >= m) then
    write(io_logs,'(1p,a,i4,77(1x,e17.10))') 'BEFORE NETBURN',l,(vxab(i),i=1,15),(vvabelx(ii,m),ii=1,nbelx)
    write(io_logs,'(i4,1p,e12.5)') l,t9
  endif

  t9=exp(t(l)-log(1.d9))
  fnucdif = 0.0d0
  if (phase >= 5.and.idifcon == 1) then
    fnucdif=0.5d0
  endif
  call netburning(l,t9,ddeit,vxab,1)

  if (l >= m) then
    write(io_logs,'(1x,a,i4,77(1x,e17.10))') 'AFTER NETBURN',l,(vxab(i),i=1,15),(abelx(ii,m),ii=1,nbelx)
  endif

  x(l)     = vxab(1)
  y3(l)    = vxab(2)
  y(l)     = vxab(3)
  xc12(l)  = vxab(4)
  xc13(l)  = vxab(5)
  xn14(l)  = vxab(6)
  xn15(l)  = vxab(7)
  xo16(l)  = vxab(8)
  xo17(l)  = vxab(9)
  xo18(l)  = vxab(10)
  xne20(l) = vxab(11)
  xne22(l) = vxab(12)
  xmg24(l) = vxab(13)
  xmg25(l) = vxab(14)
  xmg26(l) = vxab(15)
  if (  y3(l) < 1.0d-75) then
    y3(l)=0.d0
  endif
  if (xc13(l) < 1.0d-75) then
    xc13(l)=0.d0
  endif
  if (xn14(l) < 1.0d-75) then
    xn14(l)=0.d0
  endif
  if (xn15(l) < 1.0d-75) then
    xn15(l)=0.d0
  endif
  if (xo17(l) < 1.0d-75) then
    xo17(l)=0.d0
  endif
  if (xo18(l) < 1.0d-75) then
    xo18(l)=0.d0
  endif
  sumvxab=1.d0
  sumvxab=sumvxab-x(l)- y3(l)-y(l)-xc12(l)-xc13(l)-xn14(l)-xn15(l)-xo16(l)-xo17(l)- xo18(l)-xne20(l)-xne22(l)- &
                  xmg24(l)-xmg25(l)-xmg26(l)-zabelx
  do ii=1,nbelx
   sumvxab=sumvxab-abelx(ii,l)
  enddo

  if (abs(sumvxab) > 1.d-2 .and. verbose) then
    print*, l,'sumvxab= ', sumvxab
  endif

  if (l >= m) then
    write(io_logs,*) 'impl. calc',l
    write(io_logs,*) 'Dxne20: ', xne20(l) - vvxne20(l)
    write(io_logs,*) 'Dxmg24: ', xmg24(l) - vvxmg24(l)
    write(io_logs,*) 'Dxo16 : ', xo16(l) - vvxo16(l)
    write(io_logs,*) 'Dxc12 : ', xc12(l) - vvxc12(l)
    write(io_logs,*) 'Dy  : ', y(l) - vvy(l)
    write(io_logs,*) l,'sumvxab= ', sumvxab
  endif

end subroutine netc
!======================================================================
subroutine chemie
!-----------------------------------------------------------------------
! mars 2006 : adaptation du programme pour le calcul des phases avancees
!             et des etoiles de Population III
  implicit none
!-----------------------------------------------------------------------
  if (idifcon /= 1) then

! m est le nombre de couches de l'interieur stellaire
! Tant que k <= m, on va de l'exterieur vers le centre.
    do k=1,m

     if (eps(k) > 0.d0 .or. epsy(k) > 0.d0) then
       cycle
     endif

! Cas ou l'helium ne brule pas
     if (idern == 0 .or. epcne(k) <= 0.d0) then
       cycle
     endif

     if (xc12(k) < 1.d-04) then
       cycle
     endif
! Calcul de l'energie par gm et par sec divisee par Q en Mev
! *9.6485E+17/4, sauf pour S23 ou on a (Q*9.6485E+17)/23

! Photodesintegration du neon ?
     d=epcne(k)/ecne
     dbis=epcna(k)/ecna
     sc=c12ago(k)/ec
     ep23=1.0611d0*epcna(k)
     s23=ep23/e23
     s20=eps20(k)/e20
     so=o16agn(k)/eo

! Les nouvelles abondances sont calculees a partir d'equations
! dXi/dt=-epsilon/(E/quantite d'element i entrant dans la reaction
! en grammes)

     if (k == m) then
       write(io_logs,*) 'expl. calc'
       write(io_logs,*)'Dx20: ',(5.d0*d+0.87d0*s23-5.d0*s20+5.d0*so)*dzeit
       write(io_logs,*) 'Dx24: ', 6.d0*s20 *dzeit
       write(io_logs,*) 'Dxo : ', -4.d0*(so-sc)*dzeit
       write(io_logs,*) 'Dxc : ', -(6.d0*dbis+6.d0*d+3.d0*sc)*dzeit
       write(io_logs,'(a,5(1x,e12.7))')'CENTRE: Y, 12C, 16O, 20Ne, 24Mg:',y(k), xc12(k), xo16(k), xne20(k), xmg24(k)
     endif

    enddo

! Lorsque m < k, on va du centre vers l'exterieur.
    k=m-1
    iint=1
    jzint=0
    do nt=1,ixzc
     izc(nt)=0
    enddo

    do while (k > 0)
     n=0
     sumx=0.d0
     sumy=0.d0
     sumy3=0.d0
     sumxc12=0.d0
     sumc13=0.d0
     sumn14=0.d0
     sumn15=0.d0
     sumxo16=0.d0
     sumo17=0.d0
     sumo18=0.d0
     sumne20=0.d0
     sumne22=0.d0
     summg24=0.d0
     summg25=0.d0
     summg26=0.d0
     sumdm=0.d0
     if (ialflu == 1) then
       sumc14=0.d0
       sumf18=0.d0
       sumf19=0.d0
       sumn21=0.d0
       sumn23=0.d0
       suma26=0.d0
       suma27=0.d0
       sums28=0.d0
       sumbid=0.d0
       sumpro=0.d0
     endif
     do ii=1,nbelx
      sumabelx(ii)=0.d0
     enddo

     do while (zensi(k) > 0.d0)
! Si zensi(k) > 0 : Couche convective.
!             < 0 : Couche radiative.

! Homogeneisation de la zone convective.
      n=n+1
      if (n == 1 .and. k < m-1) then
        dm=(exp(q(k+2))-exp(q(k+1)))/2.d0
        dmx=dm*x(k+1)
        dmy3=dm*y3(k+1)
        dmy=dm*y(k+1)
        dmxc12=dm*xc12(k+1)
        dmc13=dm*xc13(k+1)
        dmn14=dm*xn14(k+1)
        dmn15=dm*xn15(k+1)
        dmxo16=dm*xo16(k+1)
        dmo17=dm*xo17(k+1)
        dmo18=dm*xo18(k+1)
        dmxne20=dm*xne20(k+1)
        dmxne22=dm*xne22(k+1)
        dmxmg24=dm*xmg24(k+1)
        dmxmg25=dm*xmg25(k+1)
        dmxmg26=dm*xmg26(k+1)
        sumx=sumx+dmx
        sumy3=sumy3+dmy3
        sumy=sumy+dmy
        sumxc12=sumxc12+dmxc12
        sumc13=sumc13+dmc13
        sumn14=sumn14+dmn14
        sumn15=sumn15+dmn15
        sumxo16=sumxo16+dmxo16
        sumo17=sumo17+dmo17
        sumo18=sumo18+dmo18
        sumne20=sumne20+dmxne20
        sumne22=sumne22+dmxne22
        summg24=summg24+dmxmg24
        summg25=summg25+dmxmg25
        summg26=summg26+dmxmg26
        sumdm=sumdm+dm
        if (ialflu == 1) then
          dmc14=dm*xc14(k+1)
          dmf18=dm*xf18(k+1)
          dmf19=dm*xf19(k+1)
          dmn21=dm*xne21(k+1)
          dmn23=dm*xna23(k+1)
          dma26=dm*xal26(k+1)
          dma27=dm*xal27(k+1)
          dms28=dm*xsi28(k+1)
          dmpro=dm*xprot(k+1)
          dmbid=dm*xbid(k+1)
          sumc14=sumc14+dmc14
          sumf18=sumf18+dmf18
          sumf19=sumf19+dmf19
          sumn21=sumn21+dmn21
          sumn23=sumn23+dmn23
          suma26=suma26+dma26
          suma27=suma27+dma27
          sums28=sums28+dms28
          sumpro=sumpro+dmpro
          sumbid=sumbid+dmbid
        endif
        do ii=1,nbelx
         dmabelx(ii)=dm*abelx(ii,k+1)
         sumabelx(ii)=sumabelx(ii)+dmabelx(ii)
        enddo
      endif     ! c.f. if(n==1 .and. k<m-1)

      dm=exp(q(k+1))-exp(q(k))
      dmx=dm*(x(k+1)+x(k))/2.d0
      dmy3=dm*(y3(k+1)+y3(k))/2.d0
      dmy=dm*(y(k+1)+y(k))/2.d0
      dmxc12=dm*(xc12(k+1)+xc12(k))/2.d0
      dmc13=dm*(xc13(k+1)+xc13(k))/2.d0
      dmn14=dm*(xn14(k+1)+xn14(k))/2.d0
      dmn15=dm*(xn15(k+1)+xn15(k))/2.d0
      dmxo16=dm*(xo16(k+1)+xo16(k))/2.d0
      dmo17=dm*(xo17(k+1)+xo17(k))/2.d0
      dmo18=dm*(xo18(k+1)+xo18(k))/2.d0
      dmxne20=dm*(xne20(k+1)+xne20(k))/2.d0
      dmxne22=dm*(xne22(k+1)+xne22(k))/2.d0
      dmxmg24=dm*(xmg24(k+1)+xmg24(k))/2.d0
      dmxmg25=dm*(xmg25(k+1)+xmg25(k))/2.d0
      dmxmg26=dm*(xmg26(k+1)+xmg26(k))/2.d0
      sumx=sumx+dmx
      sumy3=sumy3+dmy3
      sumy=sumy+dmy
      sumxc12=sumxc12+dmxc12
      sumc13=sumc13+dmc13
      sumn14=sumn14+dmn14
      sumn15=sumn15+dmn15
      sumxo16=sumxo16+dmxo16
      sumo17=sumo17+dmo17
      sumo18=sumo18+dmo18
      sumne20=sumne20+dmxne20
      sumne22=sumne22+dmxne22
      summg24=summg24+dmxmg24
      summg25=summg25+dmxmg25
      summg26=summg26+dmxmg26
      sumdm=sumdm+dm
      if (ialflu == 1) then
        dmc14=dm*(xc14(k+1)+xc14(k))/2.d0
        dmf18=dm*(xf18(k+1)+xf18(k))/2.d0
        dmf19=dm*(xf19(k+1)+xf19(k))/2.d0
        dmn21=dm*(xne21(k+1)+xne21(k))/2.d0
        dmn23=dm*(xna23(k+1)+xna23(k))/2.d0
        dma26=dm*(xal26(k+1)+xal26(k))/2.d0
        dma27=dm*(xal27(k+1)+xal27(k))/2.d0
        dms28=dm*(xsi28(k+1)+xsi28(k))/2.d0
        dmpro=dm*(xprot(k+1)+xprot(k))/2.d0
        dmbid=dm*(xbid(k+1)+xbid(k))/2.d0
        sumc14=sumc14+dmc14
        sumf18=sumf18+dmf18
        sumf19=sumf19+dmf19
        sumn21=sumn21+dmn21
        sumn23=sumn23+dmn23
        suma26=suma26+dma26
        suma27=suma27+dma27
        sums28=sums28+dms28
        sumpro=sumpro+dmpro
        sumbid=sumbid+dmbid
      endif
      do ii=1,nbelx
       dmabelx(ii)=dm*(abelx(ii,k+1)+abelx(ii,k))/2.d0
       sumabelx(ii)=sumabelx(ii)+dmabelx(ii)
      enddo

      if (k < 1) then
        return
      else if (k == 1) then
        k=0
! Limite superieure. Arrivee dans la region exterieure.
        dm=exp(q(1))
        sumx=dm*x(1)+sumx
        sumy=dm*y(1)+sumy
        sumy3=dm*y3(1)+sumy3
        sumxc12=dm*xc12(1)+sumxc12
        sumc13=dm*xc13(1)+sumc13
        sumn14=dm*xn14(1)+sumn14
        sumn15=dm*xn15(1)+sumn15
        sumxo16=dm*xo16(1)+sumxo16
        sumo17=dm*xo17(1)+sumo17
        sumo18=dm*xo18(1)+sumo18
        sumne20=dm*xne20(1)+sumne20
        sumne22=dm*xne22(1)+sumne22
        summg24=dm*xmg24(1)+summg24
        summg25=dm*xmg25(1)+summg25
        summg26=dm*xmg26(1)+summg26
        sumdm=dm+sumdm
        if (ialflu.eq.1) then
          sumc14=dm*xc14(1)+sumc14
          sumf18=dm*xf18(1)+sumf18
          sumf19=dm*xf19(1)+sumf19
          sumn21=dm*xne21(1)+sumn21
          sumn23=dm*xna23(1)+sumn23
          suma26=dm*xal26(1)+suma26
          suma27=dm*xal27(1)+suma27
          sums28=dm*xsi28(1)+sums28
          sumpro=dm*xprot(1)+sumpro
          sumbid=dm*xbid(1)+sumbid
        endif
        do ii=1,nbelx
         sumabelx(ii)=sumabelx(ii)+dm*abelx(ii,1)
        enddo

        exit
      endif
! Si k < 1 : On est a la surface. Arret.
! Si k = 1 : On est en la limite superieure. Region exterieure.
! Si k > 1 : On est toujours dans l'interieur.
!            On continue a "remonter" du centre vers l'exterieur.

      k=k-1
     enddo ! while (zensi(k))
! Passage a la couche immediatement plus exterieure

     if (n > 0) then
       if (k /= 0) then
! Si n <= 0 : zone radiative
! Calcul de la composition moyenne dans la region convective.
         dm=(exp(q(k+1))-exp(q(k)))/2.d0
         dmx=dm*x(k+1)
         dmy3=dm*y3(k+1)
         dmy=dm*y(k+1)
         dmxc12=dm*xc12(k+1)
         dmc13=dm*xc13(k+1)
         dmn14=dm*xn14(k+1)
         dmn15=dm*xn15(k+1)
         dmxo16=dm*xo16(k+1)
         dmo17=dm*xo17(k+1)
         dmo18=dm*xo18(k+1)
         dmxne20=dm*xne20(k+1)
         dmxne22=dm*xne22(k+1)
         dmxmg24=dm*xmg24(k+1)
         dmxmg25=dm*xmg25(k+1)
         dmxmg26=dm*xmg26(k+1)
         sumx=sumx+dmx
         sumy3=sumy3+dmy3
         sumy=sumy+dmy
         sumxc12=sumxc12+dmxc12
         sumc13=sumc13+dmc13
         sumn14=sumn14+dmn14
         sumn15=sumn15+dmn15
         sumxo16=sumxo16+dmxo16
         sumo17=sumo17+dmo17
         sumo18=sumo18+dmo18
         sumne20=sumne20+dmxne20
         sumne22=sumne22+dmxne22
         summg24=summg24+dmxmg24
         summg25=summg25+dmxmg25
         summg26=summg26+dmxmg26
         sumdm=sumdm+dm
         if (ialflu == 1) then
           dmc14=dm*xc14(k+1)
           dmf18=dm*xf18(k+1)
           dmf19=dm*xf19(k+1)
           dmn21=dm*xne21(k+1)
           dmn23=dm*xna23(k+1)
           dma26=dm*xal26(k+1)
           dma27=dm*xal27(k+1)
           dms28=dm*xsi28(k+1)
           dmpro=dm*xprot(k+1)
           dmbid=dm*xbid(k+1)
           sumc14=sumc14+dmc14
           sumf18=sumf18+dmf18
           sumf19=sumf19+dmf19
           sumn21=sumn21+dmn21
           sumn23=sumn23+dmn23
           suma26=suma26+dma26
           suma27=suma27+dma27
           sums28=sums28+dms28
           sumpro=sumpro+dmpro
           sumbid=sumbid+dmbid
         endif
         do ii=1,nbelx
          dmabelx(ii)=dm*abelx(ii,k+1)
          sumabelx(ii)=sumabelx(ii)+dmabelx(ii)
         enddo

       endif   ! if (k)
       xm=sumx/sumdm
       ym=sumy/sumdm
       y3m=sumy3/sumdm
       xc12m=sumxc12/sumdm
       xc13m=sumc13/sumdm
       xn14m=sumn14/sumdm
       xn15m=sumn15/sumdm
       xo16m=sumxo16/sumdm
       xo17m=sumo17/sumdm
       xo18m=sumo18/sumdm
       xne20m=sumne20/sumdm
       xne22m=sumne22/sumdm
       xmg24m=summg24/sumdm
       xmg25m=summg25/sumdm
       xmg26m=summg26/sumdm
       if (ialflu == 1) then
         xc14m=sumc14/sumdm
         xf18m=sumf18/sumdm
         xf19m=sumf19/sumdm
         xne21m=sumn21/sumdm
         xna23m=sumn23/sumdm
         xal26m=suma26/sumdm
         xal27m=suma27/sumdm
         xsi28m=sums28/sumdm
         xprotm=sumpro/sumdm
         xbidm=sumbid/sumdm
       endif
       do ii=1,nbelx
        mabelx(ii)=sumabelx(ii)/sumdm
       enddo

       i1=k+1
       i2=k+n+1
       izc(iint)=i2
       izc(iint+1)=i1
       iint=iint+2
       jzint=(iint-1)/2

! Composition moyenne dans la region convective.
! Dans la zone convective, l'abondance de l'element i est la meme
! en chaque couche, et egale a la valeur calculee precedemment.
       do i=i1,i2
        if (ipop3 == 1) then
! population III: fusion He: on melange sans mettre a zero
          x(i)=xm
        else
          if (x(i) /= 0.0d0 .or. idern /= 0) then
            x(i)=xm
            if (x(i) <= 1.0d-09.and.idern == 1) then
              x(i)=0.d0
            endif
          endif
        endif
        if (epsc(i) == 0.0d0) then
          y(i)=ym
        endif
        y3(i)=y3m
        xc12(i)=xc12m
        xc13(i)=xc13m
        xn14(i)=xn14m
        xn15(i)=xn15m
        xo16(i)=xo16m
        xo17(i)=xo17m
        xo18(i)=xo18m
        xne20(i)=xne20m
        xne22(i)=xne22m
        xmg24(i)=xmg24m
        xmg25(i)=xmg25m
        xmg26(i)=xmg26m
        if (ialflu == 1) then
          xf19(i)=xf19m
          xne21(i)=xne21m
          xna23(i)=xna23m
          xal26(i)=xal26m
          xal27(i)=xal27m
          xsi28(i)=xsi28m
          xbid(i)=xbidm
        endif
        do ii=1,nbelx
         abelx(ii,i)=mabelx(ii)
        enddo

        if (y3(i) < 1.0d-75) then
          y3(i)=0.d0
        endif
        if (xc13(i) < 1.0d-75) then
          xc13(i)=0.d0
        endif
        if (xn14(i) < 1.0d-75) then
          xn14(i)=0.d0
        endif
        if (xn15(i) < 1.0d-75) then
          xn15(i)=0.d0
        endif
        if (xo17(i) < 1.0d-75) then
          xo17(i)=0.d0
        endif
        if (xo18(i) < 1.0d-75) then
          xo18(i)=0.d0
        endif
        if (xprot(i) < 0.0d0) then
          xprot(i)=0.d0
        endif
        if (xneut(i) < 0.0d0) then
          xneut(i)=0.d0
        endif
        if (xal26(i) < 1.0d-75) then
          xal26(i)=0.d0
        endif
       enddo
     endif   ! if (n)
! Si k <= 1 : On continue a "remonter" vers la surface.
! Si k  > 1 : On est a la limite exterieure.
      k=k-1
    enddo
! Limite exterieure
  endif

  do i=1,m
   if (x(i) < 1.0d-75) then
     x(i)=0.d0
   endif
   if (y3(i) < 1.0d-75) then
     y3(i)=0.d0
   endif
   if (y(i) < 1.0d-75) then
     y(i)=0.d0
   endif
   if (xc12(i) < 1.0d-75) then
     xc12(i)=0.d0
   endif
   if (xc13(i) < 1.0d-75) then
     xc13(i)=0.d0
   endif
   if (xn14(i) < 1.0d-75) then
     xn14(i)=0.d0
   endif
   if (xn15(i) < 1.0d-75) then
     xn15(i)=0.d0
   endif
   if (xo16(i) < 1.0d-75) then
     xo16(i)=0.d0
   endif
   if (xo17(i) < 1.0d-75) then
     xo17(i)=0.d0
   endif
   if (xo18(i) < 1.0d-75) then
     xo18(i)=0.d0
   endif
   if (xne20(i) < 1.0d-75) then
     xne20(i)=0.d0
   endif
   if (xne22(i) < 1.0d-75) then
     xne22(i)=0.d0
   endif
   if (xmg24(i) < 1.0d-75) then
     xmg24(i)=0.d0
   endif
   if (xmg25(i) < 1.0d-75) then
     xmg25(i)=0.d0
   endif
   if (xmg26(i) < 1.0d-75) then
     xmg26(i)=0.d0
   endif
   do ii=1,nbelx
    if (abelx(ii,i) < 1.d-75) then
      abelx(ii,i)=0.d0
    endif
   enddo

  enddo

  return

end subroutine chemie
!======================================================================
subroutine chemeps
!-----------------------------------------------------------------------
!     VERSION MODIFIEE METTANT A L EQUILIBRE LES ABONDANCES DE
!     CERTAINS ELEMENTS AU DESSUS D UNE CERTAINE TEMPERATURE
!---  REVISED VERSION THAT SOLVES A REACTION NETWORK INCLUDING PP-CHAINS
!     AND CNO-TRICYCLE FOR ABUNDANCES AT T(N+1) BY A FULLY IMPLICIT
!     FINITE-DIFFERENCE METHOD SIMILAR TO THAT OF ARNETT+TRURAN (1969).
!     THIS ROUTINE WILL BE CALLED WITHIN EACH HENYEY ITERATION WITH THE
!     BEST CURRENT VALUES OF THE REACTION RATES.
!     CHEMIE (FOR DOING HOMOGENIZATION IN CONVECTIVE ZONES) WILL NOW BE
!     CALLED WITHIN NETWKI ONLY AFTER THE (ITMINC-1)ST. ITERATION AND IN
!     THE LAST ITERATION WHERE ITMINC=1.
!     THE LAST CALL TO NETWKI (AND THUS TO CHEMIE) WILL BE MADE TO
!     ESTIMATE THE COMPOSITION FOR THE NEXT TIME STEP.
!-----------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------
  call mixe(m,zensi,q,b11)
  call mixe(m,zensi,q,b33)
  call mixe(m,zensi,q,b34)
  call mixe(m,zensi,q,b112)
  call mixe(m,zensi,q,b113)
  call mixe(m,zensi,q,b114)
  call mixe(m,zensi,q,b115a)
  call mixe(m,zensi,q,b115g)
  call mixe(m,zensi,q,b116)
  call mixe(m,zensi,q,b117a)
  call mixe(m,zensi,q,b117g)
  call mixe(m,zensi,q,b118a)
  call mixe(m,zensi,q,b118g)
  call mixe(m,zensi,q,epsyy)
  call mixe(m,zensi,q,epsyc)
  call mixe(m,zensi,q,epsyo)
  call mixe(m,zensi,q,c144)
  call mixe(m,zensi,q,c184)
  call mixe(m,zensi,q,c224)
  call mixe(m,zensi,q,c134)
  call mixe(m,zensi,q,e17an)
  call mixe(m,zensi,q,c224g)
  call mixe(m,zensi,q,e20ag)

  return

end subroutine chemeps
!======================================================================
subroutine mixe(m,zensi,q,x)
!-----------------------------------------------------------------------
  implicit none

  real(kindreal),intent(in),dimension(ldi):: zensi,q
  integer,intent(in):: m
  real(kindreal),intent(out),dimension(ldi):: x
!-----------------------------------------------------------------------
  k=m-1
  do while (k  >  0)
   n=0
   sumx=0.d0
   sumdm=0.d0
   do while (zensi(k)  > 0.d0)
    n=n+1
    if (n == 1 .and. k < m-1) then
      dm=(exp(q(k+2))-exp(q(k+1)))/2.d0
      dmx=dm*x(k+1)
      sumx=sumx+dmx
      sumdm=sumdm+dm
    endif
    dm=exp(q(k+1))-exp(q(k))
    dmx=dm*(x(k+1)+x(k))/2.d0
    sumx=sumx+dmx
    sumdm=sumdm+dm
    if (k < 1) then
      return
    else if (k == 1) then
      k=0
      dm=exp(q(1))
      sumx=dm*x(1)+sumx
      sumdm=dm+sumdm
      exit
    endif
    k=k-1
   enddo   ! while (zensi(k))
   if (n > 0) then
     if (k /= 0) then
       dm=(exp(q(k+1))-exp(q(k)))/2.d0
       dmx=dm*x(k+1)
       sumx=sumx+dmx
       sumdm=sumdm+dm
     endif
     xm=(sumx)/sumdm
     i1=k+1
     i2=k+n+1
     do i=i1,i2
      x(i)=xm
     enddo
   endif   ! if (n)
   k=k-1
  enddo   ! while (k  >  0)

  return

end subroutine mixe
!======================================================================
subroutine chemold
!-----------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------
  do i=1,m
   if (vxmg24(i) < 0.d0) then
    write(91,'(1x,a,f12.7,a, i3)') ' chemold, vxmg24(i)=',vxmg24(i),' i=',i
   endif
  enddo

  do iii=1,m
   vvx(iii)=vx(iii)
   vvy3(iii)=vy3(iii)
   vvy(iii)=vy(iii)
   vvxc12(iii)=vxc12(iii)
   vvxc13(iii)=vxc13(iii)
   vvxn14(iii)=vxn14(iii)
   vvxn15(iii)=vxn15(iii)
   vvxo16(iii)=vxo16(iii)
   vvxo17(iii)=vxo17(iii)
   vvxo18(iii)=vxo18(iii)
   vvxne20(iii)=vxne20(iii)
   vvxne22(iii)=vxne22(iii)
   vvxmg24(iii)=vxmg24(iii)
   vvxmg25(iii)=vxmg25(iii)
   vvxmg26(iii)=vxmg26(iii)
   if (ialflu == 1) then
     vvxc14(iii)=vxc14(iii)
     vvxf18(iii)=vxf18(iii)
     vvxf19(iii)=vxf19(iii)
     vvxne21(iii)=vxne21(iii)
     vvxna23(iii)=vxna23(iii)
     vvxal26g(iii)=vxal26g(iii)
     vvxal27(iii)=vxal27(iii)
     vvxsi28(iii)=vxsi28(iii)
     vvxneut(iii)=vxneut(iii)
     vvxprot(iii)=vxprot(iii)
     vvxbid(iii)=vxbid(iii)
     vvxbid1(iii)=vxbid1(iii)
   endif

   do ii=1,nbelx
    vvabelx(ii,iii)=vabelx(ii,iii)
   enddo

  enddo

  if (idifcon == 1) then
    return
  endif
  k=m-1
  iint=1
  jzint=0
  do nt=1,ixzc
   izc(nt)=0
  enddo

  do while (k > 0)
   n=0
   sumx=0.d0
   sumy=0.d0
   sumy3=0.d0
   sumxc12=0.d0
   sumc13=0.d0
   sumn14=0.d0
   sumn15=0.d0
   sumxo16=0.d0
   sumo17=0.d0
   sumo18=0.d0
   sumne20=0.d0
   sumne22=0.d0
   summg24=0.d0
   summg25=0.d0
   summg26=0.d0
   sumdm=0.d0
   if (ialflu == 1) then
     sumc14=0.d0
     sumf18=0.d0
     sumf19=0.d0
     sumne21=0.d0
     sumna23=0.d0
     sumal26g=0.d0
     sumal27=0.d0
     sumsi28=0.d0
     sumneut=0.d0
     sumprot=0.d0
     sumbid=0.d0
     sumbid1=0.d0
   endif
   do ii=1,nbelx
    sumabelx(ii)=0.d0
   enddo

   do while (zensi(k) > 0.d0)
! Si zensi(k) > 0 : Couche convective.
!             < 0 : Couche radiative.

! Homogeneisation de la zone convective.
    n=n+1
    if (n == 1 .and. k < m-1) then
      dm=(exp(q(k+2))-exp(q(k+1)))/2.d0
      dmx=dm*vx(k+1)
      dmy3=dm*vy3(k+1)
      dmy=dm*vy(k+1)
      dmxc12=dm*vxc12(k+1)
      dmc13=dm*vxc13(k+1)
      dmn14=dm*vxn14(k+1)
      dmn15=dm*vxn15(k+1)
      dmxo16=dm*vxo16(k+1)
      dmo17=dm*vxo17(k+1)
      dmo18=dm*vxo18(k+1)
      dmxne20=dm*vxne20(k+1)
      dmxne22=dm*vxne22(k+1)
      dmxmg24=dm*vxmg24(k+1)
      dmxmg25=dm*vxmg25(k+1)
      dmxmg26=dm*vxmg26(k+1)
      sumx=sumx+dmx
      sumy3=sumy3+dmy3
      sumy=sumy+dmy
      sumxc12=sumxc12+dmxc12
      sumc13=sumc13+dmc13
      sumn14=sumn14+dmn14
      sumn15=sumn15+dmn15
      sumxo16=sumxo16+dmxo16
      sumo17=sumo17+dmo17
      sumo18=sumo18+dmo18
      sumne20=sumne20+dmxne20
      sumne22=sumne22+dmxne22
      summg24=summg24+dmxmg24
      summg25=summg25+dmxmg25
      summg26=summg26+dmxmg26
      sumdm=sumdm+dm
      if (ialflu == 1) then
        dmc14=dm*vxc14(k+1)
        dmf18=dm*vxf18(k+1)
        dmf19=dm*vxf19(k+1)
        dmne21=dm*vxne21(k+1)
        dmna23=dm*vxna23(k+1)
        dmal26g=dm*vxal26g(k+1)
        dmal27=dm*vxal27(k+1)
        dmsi28=dm*vxsi28(k+1)
        dmneut=dm*vxneut(k+1)
        dmprot=dm*vxprot(k+1)
        dmbid=dm*vxbid(k+1)
        dmbid1=dm*vxbid1(k+1)
        sumc14=sumc14+dmc14
        sumf18=sumf18+dmf18
        sumf19=sumf19+dmf19
        sumne21=sumne21+dmne21
        sumna23=sumna23+dmna23
        sumal26g=sumal26g+dmal26g
        sumal27=sumal27+dmal27
        sumsi28=sumsi28+dmsi28
        sumneut=sumneut+dmneut
        sumprot=sumprot+dmprot
        sumbid=sumbid+dmbid
        sumbid1=sumbid1+dmbid1
      endif
      do ii=1,nbelx
       dmabelx(ii)=dm*vabelx(ii,k+1)
       sumabelx(ii)=sumabelx(ii)+dmabelx(ii)
      enddo
    endif     !  if(n==1 .and. k<m-1)

    dm=exp(q(k+1))-exp(q(k))
    dmx=dm*(vx(k+1)+vx(k))/2.d0
    dmy3=dm*(vy3(k+1)+vy3(k))/2.d0
    dmy=dm*(vy(k+1)+vy(k))/2.d0
    dmxc12=dm*(vxc12(k+1)+vxc12(k))/2.d0
    dmc13=dm*(vxc13(k+1)+vxc13(k))/2.d0
    dmn14=dm*(vxn14(k+1)+vxn14(k))/2.d0
    dmn15=dm*(vxn15(k+1)+vxn15(k))/2.d0
    dmxo16=dm*(vxo16(k+1)+vxo16(k))/2.d0
    dmo17=dm*(vxo17(k+1)+vxo17(k))/2.d0
    dmo18=dm*(vxo18(k+1)+vxo18(k))/2.d0
    dmxne20=dm*(vxne20(k+1)+vxne20(k))/2.d0
    dmxne22=dm*(vxne22(k+1)+vxne22(k))/2.d0
    dmxmg24=dm*(vxmg24(k+1)+vxmg24(k))/2.d0
    dmxmg25=dm*(vxmg25(k+1)+vxmg25(k))/2.d0
    dmxmg26=dm*(vxmg26(k+1)+vxmg26(k))/2.d0
    sumx=sumx+dmx
    sumy3=sumy3+dmy3
    sumy=sumy+dmy
    sumxc12=sumxc12+dmxc12
    sumc13=sumc13+dmc13
    sumn14=sumn14+dmn14
    sumn15=sumn15+dmn15
    sumxo16=sumxo16+dmxo16
    sumo17=sumo17+dmo17
    sumo18=sumo18+dmo18
    sumne20=sumne20+dmxne20
    sumne22=sumne22+dmxne22
    summg24=summg24+dmxmg24
    summg25=summg25+dmxmg25
    summg26=summg26+dmxmg26
    sumdm=sumdm+dm
    if (ialflu == 1) then
      dmc14=dm*(vxc14(k+1)+vxc14(k))/2.d0
      dmf18=dm*(vxf18(k+1)+vxf18(k))/2.d0
      dmf19=dm*(vxf19(k+1)+vxf19(k))/2.d0
      dmne21=dm*(vxne21(k+1)+vxne21(k))/2.d0
      dmna23=dm*(vxna23(k+1)+vxna23(k))/2.d0
      dmal26g=dm*(vxal26g(k+1)+vxal26g(k))/2.d0
      dmal27=dm*(vxal27(k+1)+vxal27(k))/2.d0
      dmsi28=dm*(vxsi28(k+1)+vxsi28(k))/2.d0
      dmneut=dm*(vxneut(k+1)+vxneut(k))/2.d0
      dmprot=dm*(vxprot(k+1)+vxprot(k))/2.d0
      dmbid=dm*(vxbid(k+1)+vxbid(k))/2.d0
      dmbid1=dm*(vxbid1(k+1)+vxbid1(k))/2.d0
      sumc14=sumc14+dmc14
      sumf18=sumf18+dmf18
      sumf19=sumf19+dmf19
      sumne21=sumne21+dmne21
      sumna23=sumna23+dmna23
      sumal26g=sumal26g+dmal26g
      sumal27=sumal27+dmal27
      sumsi28=sumsi28+dmsi28
      sumneut=sumneut+dmneut
      sumprot=sumprot+dmprot
      sumbid=sumbid+dmbid
      sumbid1=sumbid1+dmbid1
    endif
    do ii=1,nbelx
     dmabelx(ii)=dm*(vabelx(ii,k+1)+vabelx(ii,k))/2.d0
     sumabelx(ii)=sumabelx(ii)+dmabelx(ii)
    enddo

    if (k < 1) then
      return
    else if (k == 1) then
      k=0
! Passage a la couche immediatement plus exterieure
! Limite superieure. Arrivee dans la region exterieure.
      dm=exp(q(1))
      sumx=dm*vx(1)+sumx
      sumy=dm*vy(1)+sumy
      sumy3=dm*vy3(1)+sumy3
      sumxc12=dm*vxc12(1)+sumxc12
      sumc13=dm*vxc13(1)+sumc13
      sumn14=dm*vxn14(1)+sumn14
      sumn15=dm*vxn15(1)+sumn15
      sumxo16=dm*vxo16(1)+sumxo16
      sumo17=dm*vxo17(1)+sumo17
      sumo18=dm*vxo18(1)+sumo18
      sumne20=dm*vxne20(1)+sumne20
      sumne22=dm*vxne22(1)+sumne22
      summg24=dm*vxmg24(1)+summg24
      summg25=dm*vxmg25(1)+summg25
      summg26=dm*vxmg26(1)+summg26
      sumdm=dm+sumdm
      if (ialflu == 1) then
        sumc14=dm*vxc14(1)+sumc14
        sumf18=dm*vxf18(1)+sumf18
        sumf19=dm*vxf19(1)+sumf19
        sumne21=dm*vxne21(1)+sumne21
        sumna23=dm*vxna23(1)+sumna23
        sumal26g=dm*vxal26g(1)+sumal26g
        sumal27=dm*vxal27(1)+sumal27
        sumsi28=dm*vxsi28(1)+sumsi28
        sumneut=dm*vxneut(1)+sumneut
        sumprot=dm*vxprot(1)+sumprot
        sumbid=dm*vxbid(1)+sumbid
        sumbid1=dm*vxbid1(1)+sumbid1
      endif
      do ii=1,nbelx
       sumabelx(ii)=sumabelx(ii)+dm*vabelx(ii,1)
      enddo

      exit
    endif

! Si k < 1 : On est a la surface. Arret.
! Si k = 1 : On est en la limite superieure. Region exterieure.
! Si k > 1 : On est toujours dans l'interieur.
!            On continue a "remonter" du centre vers l'exterieur.
    k=k-1
   enddo   ! while (zensi(k))

   if (n > 0) then
     if (k /= 0) then
! Si n > 0 : zone convective
! Calcul de la composition moyenne dans la region convective.
       dm=(exp(q(k+1))-exp(q(k)))/2.d0
       dmx=dm*vx(k+1)
       dmy3=dm*vy3(k+1)
       dmy=dm*vy(k+1)
       dmxc12=dm*vxc12(k+1)
       dmc13=dm*vxc13(k+1)
       dmn14=dm*vxn14(k+1)
       dmn15=dm*vxn15(k+1)
       dmxo16=dm*vxo16(k+1)
       dmo17=dm*vxo17(k+1)
       dmo18=dm*vxo18(k+1)
       dmxne20=dm*vxne20(k+1)
       dmxne22=dm*vxne22(k+1)
       dmxmg24=dm*vxmg24(k+1)
       dmxmg25=dm*vxmg25(k+1)
       dmxmg26=dm*vxmg26(k+1)
       sumx=sumx+dmx
       sumy3=sumy3+dmy3
       sumy=sumy+dmy
       sumxc12=sumxc12+dmxc12
       sumc13=sumc13+dmc13
       sumn14=sumn14+dmn14
       sumn15=sumn15+dmn15
       sumxo16=sumxo16+dmxo16
       sumo17=sumo17+dmo17
       sumo18=sumo18+dmo18
       sumne20=sumne20+dmxne20
       sumne22=sumne22+dmxne22
       summg24=summg24+dmxmg24
       summg25=summg25+dmxmg25
       summg26=summg26+dmxmg26
       sumdm=sumdm+dm
       if (ialflu == 1) then
         dmc14=dm*vxc14(k+1)
         dmf18=dm*vxf18(k+1)
         dmf19=dm*vxf19(k+1)
         dmne21=dm*vxne21(k+1)
         dmna23=dm*vxna23(k+1)
         dmal26g=dm*vxal26g(k+1)
         dmal27=dm*vxal27(k+1)
         dmsi28=dm*vxsi28(k+1)
         dmneut=dm*vxneut(k+1)
         dmprot=dm*vxprot(k+1)
         dmbid=dm*vxbid(k+1)
         dmbid1=dm*vxbid1(k+1)
         sumc14=sumc14+dmc14
         sumf18=sumf18+dmf18
         sumf19=sumf19+dmf19
         sumne21=sumne21+dmne21
         sumna23=sumna23+dmna23
         sumal26g=sumal26g+dmal26g
         sumal27=sumal27+dmal27
         sumsi28=sumsi28+dmsi28
         sumneut=sumneut+dmneut
         sumprot=sumprot+dmprot
         sumbid=sumbid+dmbid
         sumbid1=sumbid1+dmbid1
       endif
       do ii=1,nbelx
        dmabelx(ii)=dm*vabelx(ii,k+1)
        sumabelx(ii)=sumabelx(ii)+dmabelx(ii)
       enddo
     endif

     xm=(sumx)/sumdm
     ym=(sumy)/sumdm
     y3m=sumy3/sumdm
     xc12m=(sumxc12)/sumdm
     xc13m=sumc13/sumdm
     xn14m=sumn14/sumdm
     xn15m=sumn15/sumdm
     xo16m=sumxo16/sumdm
     xo17m=sumo17/sumdm
     xo18m=sumo18/sumdm
     xne20m=sumne20/sumdm
     xne22m=sumne22/sumdm
     xmg24m=summg24/sumdm
     xmg25m=summg25/sumdm
     xmg26m=summg26/sumdm
     if (ialflu == 1) then
       xc14m=sumc14/sumdm
       xf18m=sumf18/sumdm
       xf19m=sumf19/sumdm
       xne21m=sumne21/sumdm
       xna23m=sumna23/sumdm
       xal26gm=sumal26g/sumdm
       xal27m=sumal27/sumdm
       xsi28m=sumsi28/sumdm
       xneutm=sumneut/sumdm
       xprotm=sumprot/sumdm
       xbidm=sumbid/sumdm
       xbid1m=sumbid1/sumdm
     endif
     do ii=1,nbelx
      mabelx(ii)=sumabelx(ii)/sumdm
     enddo

     i1=k+1
     i2=k+n+1
     izc(iint)=i2
     izc(iint+1)=i1
     iint=iint+2
     jzint=(iint-1)/2

! Composition moyenne dans la region convective.

! Dans la zone convective, l'abondance de l'element i est la meme
! en chaque couche, et egale a la valeur calculee precedemment.
     do i=i1,i2
      if (xm <= 1.0d-09) then
        xm=0.d0
      endif

      if (xm > 1.0d-09 .and. vvx(i) < 1.d-9) then
        if (verbose) then
          write(*,*) 'limit crossed here',i, xm,vvx(i)
        endif
        write(io_logs,*) 'limit crossed here',i, xm,vvx(i)
        if (xm < 1.0d-07) then
          xm=0.d0
        endif
      endif
      vvx(i)=xm

      if (epsc(i) == 0.0d0) then
        vvy(i)=ym
      endif
      vvy3(i)=y3m
      vvxc12(i)=xc12m
      vvxc13(i)=xc13m
      vvxn14(i)=xn14m
      vvxn15(i)=xn15m
      vvxo16(i)=xo16m
      vvxo17(i)=xo17m
      vvxo18(i)=xo18m
      vvxne20(i)=xne20m
      vvxne22(i)=xne22m
      vvxmg24(i)=xmg24m
      vvxmg25(i)=xmg25m
      vvxmg26(i)=xmg26m
      if (ialflu == 1) then
        vvxf19(i)=xf19m
        vvxne21(i)=xne21m
        vvxna23(i)=xna23m
        vvxal26g(i)=xal26gm
        vvxal27(i)=xal27m
        vvxsi28(i)=xsi28m
        vvxbid(i)=xbidm
        vvxbid1(i)=xbid1m
      endif
      do ii=1,nbelx
       vvabelx(ii,i)=mabelx(ii)
      enddo

     enddo
   endif
! Si k <= 1 : On continue a "remonter" vers la surface.
! Si k  > 1 : On est a la limite exterieure.
   k=k-1
  enddo

  return

end subroutine chemold
!======================================================================
end module chemicals
