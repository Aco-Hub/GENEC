!> Module computing the EOS
!! @author Sylvia Ekstrom & Cyril Georgy
!! @version 278
!! @date 26.3.2013
!! @brief EOS computation with a mixture of perfect gas and radiation
!======================================================================
module EOS

  use evol,only: kindreal
  use const,only: cst_c,rgaz,cst_me

  implicit none

  integer,save:: num

  real(kindreal),save:: psi,pl,toni,rhe
  real(kindreal),save:: rh,rh1,rhp,rhp1,rht,rht1
  real(kindreal),save:: uta,rhpsi,rhpsip,rhpsit
  real(kindreal),save:: chi,hpsi
  real(kindreal):: tk,pg,vermy,rhes,rhete,pe,pes,pete,ue,ues,uete

private
public :: dichte
public :: num,rh,rh1,rhp,rhp1,rht,rht1,rhe,psi,rhpsi,rhpsip,rhpsit,toni,pl,uta

contains
!======================================================================
subroutine degen
!-----------------------------------------------------------------------
  use const,only: pi,cst_h,cst_mh

  implicit none

  integer:: i
  real(kindreal):: hchi,psi,sr,srs,srt=0.d0,sp,sps,spt=0.d0,su=0.d0,sus=0.d0,sut=0.d0,tu,wz,asr,asp,asu,xsq,x,wx,gol,phi1,ph1, &
    phi2,ph2,phi1d,ph1d,phi2d,ph2d,phi3,ph3,phi3d,ph3d

 real(kindreal)::cst_mch3,ccd,ccr,ccp,ccu
! Les data suivants (dega,degu,degex) sont tires de Kippenhahn, Weigert & Hofmeister 1967
! Methods in Computational Physics vol. 7, p. 178:
 real(kindreal), dimension(12):: &
! a_i:
  dega=(/7.28905d-2,2.04648d-1,2.54685d-1,1.96668d-1,1.04294d-1,3.95716d-2,1.09024d-2,2.21504d-3,3.14295d-4,3.64498d-5,2.42533d-6, &
         1.49830d-7/), &
! u_i:
  degu=(/0.11896d0,0.47652d0,1.07476d0,1.91722d0,3.00898d0,4.35703d0,5.96967d0,7.86147d0,10.02706d0,12.57936d0,15.19876d0, &
         18.87814d0/), &
! alpha_i:
  degex=(/8.87844d-1,6.20940d-1,3.41379d-1,1.47015d-1,4.93418d-2,1.28164d-2,2.55507d-3,3.85307d-4,4.41879d-5,3.44234d-6, &
          2.50763d-7,6.32888d-9/)
!-----------------------------------------------------------------------
  cst_mch3 = (cst_me*cst_c/cst_h)**3.d0
  ccd = pi**2.d0 / 6.d0
  ccr = 8.d0*pi*cst_mh*cst_mch3
  ccp = (16.d0*pi/6.d0)*cst_me*cst_c**2.d0*cst_mch3
  ccu = cst_me*cst_c**2.d0/cst_mh
  if (hpsi >= 0.d0) then
    hchi=1.d0
    psi=hpsi
    xsq=2.d0*tk*psi+(tk*psi)**2.d0
    x= sqrt(xsq)
    wx=1.d0+tk*psi
    gol=log(wx-x)
    phi1=xsq*x/3.d0
    ph1=tk*tk*(1.d0+2.d0*xsq)/x
    if ((x-0.5d0) > 0.d0) then
      phi2=x*wx*(2.d0*xsq-3.d0)*0.125d0-0.375d0*gol
    else
      phi2=0.125d0*x*xsq*xsq*(8.d0/5.d0-xsq*(4.d0/7.d0-xsq*(1.d0/3.d0-xsq*5.d0/22.d0)))
    endif
    ph2=tk**2.d0*3.d0*x*wx
    phi1d=x*wx
    ph1d=(2.d0*xsq-1.d0)*tk*tk*wx/(xsq*x)
    phi2d=x*xsq
    ph2d=(6.d0*xsq+3.d0)*tk*tk/x
    sr=phi1+ccd*ph1
    srs=-(phi1d+ccd*ph1d)*tk/hchi
    sp=phi2+ccd*ph2
    sps=-(phi2d+ccd*ph2d)*tk/hchi
    if ((hpsi+1.d0) >= 0.d0) then
      if ((x-0.5d0) > 0.d0) then
        phi3=x*wx*(2.d0*xsq+1.d0)*0.125d0-phi1+0.125d0*gol
      else
        phi3=x*xsq*xsq*(0.1d0-xsq*(1.d0/56.d0-xsq*(1.d0/144.d0-xsq*5.d0/1408.d0)))
      endif
      ph3=tk*tk*(wx*(1.d0+3.d0*xsq)-1.d0-2.d0*xsq)/x
      phi3d=x*wx*(wx-1.d0)
      ph3d=tk*tk*(xsq*(6.d0*xsq+3.d0)-1.d0-wx*(2.d0*xsq-1.d0))/(xsq*x)
      srt=(phi1d+ccd*ph1d)*psi+2.d0*ccd/tk*ph1
      spt=(phi2d+ccd*ph2d)*psi+2.d0*ccd/tk*ph2
      su=phi3+ccd*ph3
      sus=-(phi3d+ccd*ph3d)*tk/hchi
      sut=(phi3d+ccd*ph3d)*psi+2.d0*ccd/tk*ph3
    endif
  else
    hchi=chi
    if ((hchi-1.d0/exp(7.d0)) >= 0.d0) then
      sr=0.d0
      srs=0.d0
      srt=0.d0
      sp=0.d0
      sps=0.d0
      spt=0.d0
      su=0.d0
      sus=0.d0
      sut=0.d0
! PATENAUDE 1974, p. 51 sqq:
      do i=1,12
       tu=tk*degu(i)
       wz=sqrt (tk*(2.d0+tu))
       asr=(1.d0+tu)*tk*wz*dega(i)/(hchi+degex(i))
       asp=tu*wz*wz*wz*dega(i)/(hchi+degex(i))
       sr=sr+asr
       srs=srs-asr/(hchi+degex(i))
       sp=sp+asp
       sps=sps-asp/(hchi+degex(i))
       if ((hpsi+1.d0) < 0.d0) cycle
       srt=srt+(3.d0+7.d0*tu+3.d0*tu*tu)*tk/wz*dega(i)/(hchi+degex(i))
       spt=spt+(5.d0+4.d0*tu)*tu*wz*dega(i)/(hchi+degex(i))
       asu=(1.d0+tu)*tu*tk*wz*dega(i)/(hchi+degex(i))
       su=su+asu
       sus=sus-asu/(hchi+degex(i))
       sut=sut+(5.d0+10.d0*tu+4.d0*tu*tu)*tu*tk/wz*dega(i)/(hchi+degex(i))
      enddo
    else
      psi=-log(hchi)
      xsq=2.d0*tk*psi+(tk*psi)**2.d0
      x= sqrt(xsq)
      wx=1.d0+tk*psi
      gol=log(wx-x)
      phi1=xsq*x/3.d0
      ph1=tk*tk*(1.d0+2.d0*xsq)/x
      if ((x-0.5d0) > 0.d0) then
       phi2=x*wx*(2.d0*xsq-3.d0)*0.125d0-0.375d0*gol
      else
        phi2=0.125d0*x*xsq*xsq*(8.d0/5.d0-xsq*(4.d0/7.d0-xsq*(1.d0/3.d0-xsq*5.d0/22.d0)))
      endif
      ph2=tk**2.d0*3.d0*x*wx
      phi1d=x*wx
      ph1d=(2.d0*xsq-1.d0)*tk*tk*wx/(xsq*x)
      phi2d=x*xsq
      ph2d=(6.d0*xsq+3.d0)*tk*tk/x
      sr=phi1+ccd*ph1
      srs=-(phi1d+ccd*ph1d)*tk/hchi
      sp=phi2+ccd*ph2
      sps=-(phi2d+ccd*ph2d)*tk/hchi
      if ((hpsi+1.d0) >= 0.d0) then
        if ((x-0.5d0) > 0.d0) then
          phi3=x*wx*(2.d0*xsq+1.d0)*0.125d0-phi1+0.125d0*gol
        else
          phi3=x*xsq*xsq*(0.1d0-xsq*(1.d0/56.d0-xsq*(1.d0/144.d0-xsq*5.d0/1408.d0)))
        endif
        ph3=tk*tk*(wx*(1.d0+3.d0*xsq)-1.d0-2.d0*xsq)/x
        phi3d=x*wx*(wx-1.d0)
        ph3d=tk*tk*(xsq*(6.d0*xsq+3.d0)-1.d0-wx*(2.d0*xsq-1.d0))/(xsq*x)
        srt=(phi1d+ccd*ph1d)*psi+2.d0*ccd/tk*ph1
        spt=(phi2d+ccd*ph2d)*psi+2.d0*ccd/tk*ph2
        su=phi3+ccd*ph3
        sus=-(phi3d+ccd*ph3d)*tk/hchi
        sut=(phi3d+ccd*ph3d)*psi+2.d0*ccd/tk*ph3
      endif
    endif
  endif
  rhe=ccr*sr
  rhes=srs/sr
  pe=ccp*sp
  pes=sps/sp
  if ((hpsi+1.d0) >= 0.d0) then
    rhete=tk*srt/sr
    pete=tk*spt/sp
    ue=ccu*su/sr
    ues=(sr*sus-su*srs)/(sr*su)
    uete=tk*(sr*sut-su*srt)/(sr*su)
  endif

  return

end subroutine degen
!======================================================================
subroutine Compute_sepp(krueck,sepp,seppl,dchi)
!-----------------------------------------------------------------------
  implicit none

  integer,intent(in):: krueck
  real(kindreal),intent(out):: sepp,seppl,dchi
!-----------------------------------------------------------------------
  call degen

  select case (krueck)
  case (1)
    sepp=(pe+rgaz*vermy*toni*rhe)/pg
    seppl=(pe/pg)*chi*pes+(rgaz*vermy*toni*rhe)/pg*chi*rhes
    dchi=((1.d0-sepp)/seppl)*chi
  case (2)
    sepp=(pe+rgaz*vermy*toni*rhe)/pg
    seppl=(pe/pg)*pes
    dchi =sepp*log(sepp)/seppl
  case default
    stop 'Bad value for krueck in dichte'
  end select

  return

end subroutine Compute_sepp
!======================================================================
subroutine dichte
!-----------------------------------------------------------------------
  use const,only: cst_a,um,cst_k
  use inputparam,only: ialflu,z
  use abundmod,only: x,y3,y,xc12,xc13,xc14,xn14,xn15,xo16,xo17,xo18,xf18,xf19,xne20,xne21,xne22,xna23,xmg24,xmg25,xmg26, &
                     xal26,xal27,xsi28,zabelx,nbelx,abelx,nbzel,nbael
  use strucmod,only: j1,vmy1,vmyo,vmye,j,beta1,t,p,adi1,vmol,vny,x_env,vmion,vna,vmionp,vmiont
  use ionisation,only: ionpart

  implicit none

  integer:: ii,mist,krueck
  real(kindreal), save:: vx3
  real(kindreal):: refad,tion,pion,xcmp,sepp,seppl,dchi,rhes1,rpsi,hifa,pas,pate
!-----------------------------------------------------------------------
  if (j1 == 1) then
    num=0
    vx3=0.d0
    chi=55.d0
    hpsi=-3.d0
  endif

  vmy1 = 2.d0*x(j1)+y3(j1)+3.d0/4.d0*y(j1)+7.d0/12.d0*xc12(j1)+7.d0/13.d0*xc13(j1)+8.d0/14.d0*xn14(j1)+8.d0/15.d0*xn15(j1)+ &
         9.d0/16.d0*xo16(j1)+9.d0/17.d0*xo17(j1)+9.d0/18.d0*xo18(j1)+11.d0/20.d0*xne20(j1)+11.d0/22.d0*xne22(j1)+ &
         13.d0/24.d0*xmg24(j1)+13.d0/25.d0*xmg25(j1)+13.d0/26.d0*xmg26(j1)

  vmyo = x(j1)+y3(j1)/3.d0+y(j1)/4.d0+xc12(j1)/12.d0+xc13(j1)/13.d0+xn14(j1)/14.d0+xn15(j1)/15.d0+xo16(j1)/16.d0+ &
         xo17(j1)/17.d0+xo18(j1)/18.d0+xne20(j1)/20.d0+xne22(j1)/22.d0+xmg24(j1)/24.d0+xmg25(j1)/25.d0+xmg26(j1)/26.d0

  vmye=0.5d0*(1.d0+x(j1))+y3(j1)/6.d0

  if (ialflu == 1) then
    vmy1 = vmy1+0.5d0*xc14(j1)+10.d0/18.d0*xf18(j1)+10.d0/19.d0*xf19(j1)+11.d0/21.d0*xne21(j1)+12.d0/23.d0*xna23(j1)+ &
           14.d0/26.d0*xal26(j1)+14.d0/27.d0*xal27(j1)+15.d0/28.d0*xsi28(j1)
    vmyo = vmyo+xc14(j1)/14.d0+xf18(j1)/18.d0+xf19(j1)/19.d0+xne21(j1)/21.d0+xna23(j1)/23.d0+xal26(j1)/26.d0+ &
           xal27(j1)/27.d0+xsi28(j1)/28.d0
    vmye = x(j1)+2.d0/3.d0*y3(j1)+0.5d0*(y(j1)+xc12(j1)+xn14(j1)+xo16(j1)+xne20(j1)+xmg24(j1)+xal26(j)+xsi28(j)+xf18(j))+ &
           6.d0/13.d0*xc13(j1)+7.d0/15.d0*xn15(j1)+8.d0/17.d0*xo17(j1)+8.d0/18.d0*xo18(j1)+10.d0/22.d0*xne22(j1)+ &
           12.d0/25.d0*xmg25(j1)+12.d0/26.d0*xmg26(j1)+6.d0/14.d0*xc14(j1)+9.d0/19.d0*xf19(j1)+10.d0/21.d0*xne21(j1)+ &
           11.d0/23.d0*xna23(j1)+13.d0/27.d0*xal27(j1)
  endif

!  z elements taken ~ as Ca56: A=56 but (nbzel(ii)+1.)/nbael(ii)-->0.5
  vmy1= vmy1+0.5d0*zabelx
  vmyo= vmyo+zabelx/56.d0
  if (ialflu == 1) then
    vmye = x(j1)+2.d0/3.d0*y3(j1)+0.5d0*(y(j1)+xc12(j1)+xn14(j1)+xo16(j1)+xne20(j1)+xmg24(j1)+xal26(j)+xsi28(j)+ &
           xf18(j))+6.d0/13.d0*xc13(j1)+6.d0/14.d0*xc14(j1)+7.d0/15.d0*xn15(j1)+8.d0/17.d0*xo17(j1)+8.d0/18.d0*xo18(j1)+ &
           9.d0/19.d0*xf19(j1)+10.d0/21.d0*xne21(j1)+10.d0/22.d0*xne22(j1)+11.d0/23.d0*xna23(j1)+12.d0/25.d0*xmg25(j1)+ &
           12.d0/26.d0*xmg26(j1)+13.d0/27.d0*xal27(j1)+0.5d0*zabelx
  else
    vmye = x(j1)+2.d0/3.d0*y3(j1)+0.5d0*(y(j1)+xc12(j1)+xn14(j1)+xo16(j1)+xne20(j1)+xmg24(j1))+6.d0/13.d0*xc13(j1)+ &
           7.d0/15.d0*xn15(j1)+8.d0/17.d0*xo17(j1)+8.d0/18.d0*xo18(j1)+10.d0/22.d0*xne22(j1)+12.d0/25.d0*xmg25(j1)+ &
           12.d0/26.d0*xmg26(j1)+0.5d0*zabelx
  endif

  do ii=1,nbelx
   vmy1= vmy1+abelx(ii,j1)*real(nbzel(ii)+1)/real(nbael(ii))
   vmyo= vmyo+abelx(ii,j1)/real(nbael(ii))
   vmye= vmye+abelx(ii,j1)*real(nbzel(ii))/real(nbael(ii))
  enddo

  vmy1=-log(vmy1)
  vmyo=1.d0/vmyo
  vmye=1.d0/vmye
  vermy=vmye/vmyo

  beta1=1.d0- exp(log(cst_a)-log(3.d0)+4.d0*t(j1)-p(j1))
  beta1=max(beta1,1.d-5)

  refad=p(j1)-2.5d0*t(j1)

  refad=refad+3.68d0

  if (refad <= 0.d0) then

    if (num == 0) then
!-----------------------------------------------------------------------
!  Partie modifiee Mai 1990, D.Schaerer:
!    IONISATION PARTIELLE DE PLUSIEURS ELEMENTS...

!  The partial ionisation is computed for the first layers of the interior.
!  num=0 when entering with the layer j=1.
!  ionpart is called. As long as x_env(3)=0 or x_env(3) changes from one layer to the next
!  (comparison with vx3 that memorises the previous value), num keeps its value 0.
!  When x_env(3) doesn't change anymore from one layer to the other, num=j
!  and the partial ionisation is not computed anymore.
      tion=t(j1)/um
      pion=p(j1)/um
      vmol=vmyo
      vny(1)=vmyo*x(j1)
      vny(2)=vmyo*y(j1)/4.d0
      vny(3)=vny(2)
      xcmp=x(j1)

      call ionpart(pion,tion)

!---------------- Modifications de Schaerer 1990 -----------------------
! On prend les valeurs calculees par IONPART pour vmion,vna,vmionp,vmiont
!      if (x_env(3)-vx3)33,33,32
      if (x_env(3)-vx3 > 0.d0) then
        vmy1=log(vmion)
        adi1=vna
        vx3=x_env(3)

! modification apportee 19 XI 97
! Lorsque FITM a une valeur tres elevee on peut imaginer
! qu'il y ait des situations ou x_env(3)=0 et donc ou lors des
! premiers passages dans DICHTE x_env(3)=vx3. Les modifications apportees
! ici ont pour but de faire en sorte que dans cette situation on
! calcule les effets de l'ionisation partielle
      else
        if (x_env(3) /= 0.d0) then
          num=j1
        else
          vmy1=log(vmion)
          adi1=vna
          vx3=x_env(3)
        endif
      endif
    endif   ! num

    rh1=-log(rgaz)+vmy1+log(beta1)+p(j1)-t(j1)
    rhp1=1.d0/beta1+vmionp
    rht1=3.d0-4.d0/beta1+vmiont
!---------------- Fin des modifications --------------------------------
    return
!------------------------------------------------------------------------
  else   ! refad > 0.
! Calcul pour degenerescence
! psi <= 7: degenerescence electronique partielle
! psi >  7: degenerescence electronique totale
    pl=exp(p(j1))
    num = -1000
    pg=beta1*pl
    toni=exp(t(j1))
    tk=toni*cst_k/(cst_me*cst_c**2.d0)
    mist = 0
    if (hpsi > 0.d0) then
      krueck = 2
      hpsi = -log(chi)
      psi = hpsi
      call Compute_sepp(krueck,sepp,seppl,dchi)
    else
      krueck = 1
      call Compute_sepp(krueck,sepp,seppl,dchi)
    endif

    do while (mist <= 45)
     mist = mist+1
     if (krueck == 1) then
       if (chi+dchi <= 0.d0) then
         chi = 0.1d0*chi
         if (chi - 1.d-18 <= 0.d0) then
           krueck = 2
           hpsi = -log(chi)
           psi = hpsi
           call Compute_sepp(krueck,sepp,seppl,dchi)
           cycle
         else
           krueck = 1
           call Compute_sepp(krueck,sepp,seppl,dchi)
           cycle
         endif
       else
         chi = chi+dchi
         if (chi - 1.d-18 <= 0.d0) then
           krueck = 2
           hpsi = -log(chi)
           psi = hpsi
           call Compute_sepp(krueck,sepp,seppl,dchi)
           cycle
         else
           if (abs(1.d0-sepp)-0.002d0 >= 0.d0) then
             krueck = 1
             call Compute_sepp(krueck,sepp,seppl,dchi)
             cycle
           else
             psi=-log(chi)
             chi=chi*0.99d0
             call degen
             rhes1=rhes
             chi=chi/0.99d0
             hpsi=-1.d0
             call degen
             hpsi=-3.d0
             rpsi=100.d0*(rhes1-rhes)-rhes
             rhpsi=-chi*rhes
             rh1=log(rhe*vmye)
             hifa=rgaz*vermy*rhe*toni/pl
             pas=hifa*rhes+pe/pl*pes
             pate=hifa*(1.d0+rhete)+pe*pete/pl+4.d0*(1.d0-beta1)
             uta=(ue/toni)*(uete-ues*pate/pas)/vmye
             rhp1=rhes/pas
             rht1=rhete-rhp1*pate
             rhpsip=rpsi/pas
             rhpsit=-pate*rhpsip
             return
           endif
         endif
       endif
     endif

     do while (hpsi+dchi <= 0.d0)
      dchi = 0.5d0*dchi
     enddo
     hpsi = hpsi+dchi
     chi = exp(-hpsi)
     psi = hpsi
     if (hpsi-35.d0 <= 0.d0) then
       krueck = 1
       hpsi = -3.d0
       call Compute_sepp(krueck,sepp,seppl,dchi)
       cycle
     else
       if (abs(1.d0-sepp)-0.002d0 >= 0.d0) then
         krueck = 2
         call Compute_sepp(krueck,sepp,seppl,dchi)
         cycle
       else
         hpsi=hpsi+0.1d0
         call degen
         rhes1=rhes
         hpsi=hpsi-0.1d0
         call degen
         rpsi=(rhes1-rhes)/(0.1d0*chi)
         rhpsi=-rhes
         rh1=log(rhe*vmye)
         hifa=rgaz*vermy*rhe*toni/pl
         pas=hifa*rhes+pe/pl*pes
         pate=hifa*(1.d0+rhete)+pe*pete/pl+4.d0*(1.d0-beta1)
         uta=(ue/toni)*(uete-ues*pate/pas)/vmye
         rhp1=rhes/pas
         rht1=rhete-rhp1*pate
         rhpsip=rpsi/pas
         rhpsit=-pate*rhpsip
         return
       endif
     endif

    enddo

    select case (krueck)
    case (1)
      psi=-log(chi)
      chi=chi*0.99d0
      call degen
      rhes1=rhes
      chi=chi/0.99d0
      hpsi=-1.d0
      call degen
      hpsi=-3.d0
      rpsi=100.d0*(rhes1-rhes)-rhes
      rhpsi=-chi*rhes
    case (2)
      hpsi=hpsi+0.1d0
      call degen
      rhes1=rhes
      hpsi=hpsi-0.1d0
      call degen
      rpsi=(rhes1-rhes)/(0.1d0*chi)
      rhpsi=-rhes
    case default
      stop 'Bad value for krueck in dichte'
    end select
    rh1=log(rhe*vmye)
    hifa=rgaz*vermy*rhe*toni/pl
    pas=hifa*rhes+pe/pl*pes
    pate=hifa*(1.d0+rhete)+pe*pete/pl+4.d0*(1.d0-beta1)
    uta=(ue/toni)*(uete-ues*pate/pas)/vmye
    rhp1=rhes/pas
    rht1=rhete-rhp1*pate
    rhpsip=rpsi/pas
    rhpsit=-pate*rhpsip
    return
  endif

end subroutine dichte
!======================================================================
end module EOS
