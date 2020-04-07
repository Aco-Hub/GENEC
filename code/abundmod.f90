! module contenant les abondances, taux de reaction, opacites
module abundmod

  use evol,only: ldi,kindreal
  use inputparam,only: ialflu,verbose

  implicit none

  integer, save:: nbelx,lcnom
  integer, parameter:: mbelx=10
  integer, dimension (mbelx), save:: nbzel,nbael

! ABONDANCES:
  real(kindreal), dimension(ldi), save:: x,y3,y,xc12,xc13,xc14,xn14,xn15,xo16,xo17,xo18,xf18,xf19,xne20,xne21,xne22,xna23,xmg24, &
    xmg25,xmg26,xal26,xal27,xsi28,xprot,xneut,xbid,xbid1
  real(kindreal), dimension(ldi), save:: vx,vy3,vy,vxc12,vxc13,vxc14,vxn14,vxn15,vxo16,vxo17,vxo18,vxf18,vxf19,vxne20,vxne21, &
    vxne22,vxna23,vxmg24,vxmg25,vxmg26,vxal26g,vxal27,vxsi28,vxprot,vxneut,vxbid,vxbid1
  real(kindreal), dimension(ldi), save:: vvx,vvy3,vvy,vvxc12,vvxc13,vvxc14,vvxn14,vvxn15,vvxo16,vvxo17,vvxo18,vvxf18,vvxf19, &
    vvxne20,vvxne21,vvxne22,vvxna23,vvxmg24,vvxmg25,vvxmg26,vvxal26g,vvxal27,vvxsi28,vvxprot,vvxneut,vvxbid,vvxbid1
! en numerotation inverse (centre: 1 --> surface: m)
  real(kindreal), dimension(ldi), save:: wx,wy3,wy,wxc12,wxc13,wxc14,wxn14,wxn15,wxo16,wxo17,wxo18,wxf18,wxf19,wxne20,wxne21, &
    wxne22,wxna23,wxmg24,wxmg25,wxmg26,wxal26g,wxal27,wxsi28,wxprot,wxneut,wxbid,wxbid1
! TAUX DE REACTIONS, ENERGIE:
  real(kindreal), dimension(ldi), save:: eps,epsy,epsyy,epsyc,epsyo,epsc,b11,b33,b34,b112,b113,b114,b115a,b115g,b116,b117a,b117g, &
    b118a,b118g,epcne,eps20,c144,c184,c224,c134,c12ago,o16agn,epcna,eps_c_adv,eps_ne_adv,eps_o_adv,eps_si_adv,eps_grav,eps_nu
  real(kindreal), dimension(ldi), save:: b119a,b119g,b120,b121,b122,b123g,b123a,b124,b125g,b125m,b1mg26,b1al26,b127g,b127a,e24ag, &
    e17ag,e21ag,e18an,e21na,e25an,e20ng,e21ng,e22ng,e23ng,e24ng,e25ng,e26ng,e27ng,e28ng,a26ga,a26gp,e14np,ec14pg,ec14ag,ef18na, &
    e15ag,ef18np,e18pa,ec14ng,e19ap,e14be,e18be,e26be,e18ng
  real(kindreal), dimension(ldi), save:: e17an,c224g,e20ag,e12ng,e14ng,e19ng
  real(kindreal), dimension(ldi), save:: enpp,encno,densityj
  real(kindreal), dimension(13,ldi), save:: en13
  real(kindreal), save:: eg,egp,egp1,egt,egt1,en,enue,enuet,enuet1,enuep,enuep1,epsn,epsn1,epsp,epsp1,epst,epst1
  real(kindreal), save:: xmcno,scno
! NEUTRINOS:
  real(kindreal), dimension(ldi), save:: tauxbe7pg,tauxbe7el,ybe7,yb8,xnube7,xnub8,xnbrbe7,xnbrb8,b34neu,xfluxn
  real(kindreal)::snube7,snub8,snu
! ELEMENTS DES PHASES AVANCEES
  real(kindreal), save:: zabelx,fnucdif,xlostneu
  real(kindreal), dimension (mbelx,ldi), save:: abelx,vabelx,vvabelx,wabelx,sabelx
  real(kindreal), dimension (mbelx), save:: abels
  real(kindreal), dimension (ldi), save:: abelt,vabelt,vvabelt,wabelt,sabelt

! ENERGIE:
  real(kindreal), save:: ekrote,epote,ekine,erade

public

contains
!**********************************************************************
subroutine abundCheck(m,DoPrint)
!-----------------------------------------------------------------------
  use inputparam,only: idebug

  implicit none

  integer,intent(in):: m
  logical,intent(in):: DoPrint

  real(kindreal),parameter:: XTolMax=1.0d-3
  integer:: ii,lqt,lqtm
  real(kindreal):: cabmal,cabmad,cabm
  real(kindreal),dimension(ldi):: cab,cabal,cabad
!-----------------------------------------------------------------------
    lqtm=1
    cabmal=0.d0
    cabmad=0.d0
    cabm=x(1)+y3(1)+y(1)+xc12(1)+xc13(1)+xn14(1)+xn15(1)+xo16(1)+xo17(1)+xo18(1)+xne20(1)+xne22(1)+xmg24(1)+xmg25(1)+xmg26(1)+ &
        (xmg25(1)-xlostneu)/25.d0

    if (ialflu == 1) then
      cabmal=xf19(1)+xne21(1)+xal26(1)+xal27(1)+xna23(1)+xneut(1)+xprot(1)+xc14(1)+xf18(1)+xbid(1)+xsi28(1)+xbid1(1)
    endif

    do ii=1,nbelx
     cabmad=cabmad+abelx(ii,1)
    enddo
    if (ialflu /= 1) cabmad=cabmad+zabelx

    cabm=cabm+cabmal+cabmad
    cabm=abs(1.d0-cabm)

    do lqt=2,m
     cabal(lqt)=0.d0
     cabad(lqt)=0.d0
     cab(lqt)=x(lqt)+y3(lqt)+y(lqt)+xc12(lqt)+xc13(lqt)+xn14(lqt)+xn15(lqt)+xo16(lqt)+xo17(lqt)+xo18(lqt)+xne20(lqt)+xne22(lqt)+ &
              xmg24(lqt)+xmg25(lqt)+xmg26(lqt)+(xmg25(lqt)-xlostneu)/25.d0

     if (ialflu == 1) then
       cabal(lqt)=xf19(lqt)+xne21(lqt)+xna23(lqt)+xal26(lqt)+xal27(lqt)+xsi28(lqt)+xneut(lqt)+xprot(lqt)+xc14(lqt)+xf18(lqt)+ &
                  xbid(lqt)+xbid1(lqt)
     endif

     do ii=1,nbelx
      cabad(lqt)=cabad(lqt)+abelx(ii,lqt)
     enddo
     if (ialflu /= 1) cabad(lqt)=cabad(lqt)+zabelx

     cab(lqt)=cab(lqt)+cabal(lqt)+cabad(lqt)
     cab(lqt)=abs(1.d0-cab(lqt))

     if (DoPrint) then
       if (lqt == m) then
         write(3,'(2x,a,i4,77(1x,e17.10))') 'SOMME DES XI. ECAR MAX.',lqt,cab(lqt),xc12(lqt),xo16(lqt),xne20(lqt), &
                                                                      xmg24(lqt),(abelx(ii,lqt),ii=1,nbelx)
       endif
     endif
     if (cabm < cab(lqt)) then
       lqtm=lqt
       cabm=cab(lqt)
     endif
    enddo

    if (idebug > 2) then
      write(*,*)  'SOMME DES XI. ECAR MAX.',lqtm,cabm
    endif
    if (DoPrint) then
      write(3,'(2x,a,i4,2x,e14.7)') 'SOMME DES XI. ECAR MAX.',lqtm,cabm
      write(10,'(2x,a,i4,2x,e14.7)') 'SOMME DES XI. ECAR MAX.',lqtm,cabm
      if (abs(cabm - 1.d0) > XTolMax) then
        write(3,'(1x,3(2x,f11.8),2(2x,e8.2),2(2x,f11.8),3(2x,e8.2),16(2x,f11.8))') x(lqtm),y3(lqtm),y(lqtm),xc12(lqtm),xc13(lqtm), &
          xn14(lqtm),xn15(lqtm),xo16(lqtm),xo17(lqtm),xo18(lqtm),xne20(lqtm),xne22(lqtm),xmg24(lqtm),xmg25(lqtm),xmg26(lqtm), &
          xf19(lqtm),xne21(lqtm),xna23(lqtm),xal26(lqtm),xal27(lqtm),xsi28(lqtm),xneut(lqtm),xprot(lqtm),xc14(lqtm),xf18(lqtm), &
          xbid(lqtm)
      endif
      write(3,'(1x/,i5,1p,77(e11.3))')lqtm,(abelx(ii,lqtm),ii=1,nbelx)
    endif

    return

end subroutine abundCheck
!**********************************************************************
subroutine maxCNO(m,q)
!----------------------------------------------------------------------
! calcul du max de la somme C12+C13+N14+N15+O16+O17+O18 ds zone ou H brule
!----------------------------------------------------------------------
  use caramodele,only: gms

  implicit none

  integer,intent(in):: m
  real(kindreal),dimension(ldi),intent(in):: q

  integer:: i
  real(kindreal):: xaxn14
  real(kindreal),dimension(ldi):: cno
!----------------------------------------------------------------------
    lcnom=1
    xmcno=gms
    xaxn14=xn14(1)
    scno= xc12(1)+xc13(1)+xn14(1)+xn15(1)+xo16(1)+xo17(1)+xo18(1)
    do i=2,m
     cno(i)=xc12(i)+xc13(i)+xn14(i)+xn15(i)+xo16(i)+xo17(i)+xo18(i)
     if (xaxn14 < xn14(i)) then
       lcnom=i
       xaxn14=xn14(i)
       scno=cno(i)
       xmcno=(1.d0-exp(q(i)))*gms
     endif
    enddo



end subroutine maxCNO
!**********************************************************************
end module abundmod
