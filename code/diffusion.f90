module diffusion
  use evol, only: kindreal,ldi
  use const,only: pi
  use inputparam,only: verbose
  use magmod, only: qmin

  implicit none

  real(kindreal):: gmsu
  real(kindreal), dimension(ldi):: M_r,D_chim,D_moychim,D_Omega,D_moyOm
! CONDITION DE COURAN
  integer, save:: nwpas
  real(kindreal), dimension(ldi), save:: dra
  real(kindreal), dimension(0:ldi), save:: omega_extended

  private
  public:: coedif
  public:: diffbr,diffom

contains
!***********************************************************************
!> Computes the diffusion coefficients D_h, D_shear, D_conv and K_ther.
!! Verifies the stability toward the Richardson criterion and the dynamical shear.
!! Uses the ICOEFF input parameter for the choice of the rotation prescriptions:
!! ICOEFF = X1 : D_h Zahn (1992) A&A 265,115
!! ICOEFF = X2 : D_h Maeder (2003) A&A 399,263
!! ICOEFF = X3 : D_h Mathis et al. (2004) A&A 425,243
!! ICOEFF = 1X : D_shear Maeder (1997) A&A 321,134
!! ICOEFF = 2X : D_shear Talon & Zahn (1997) A&A 317,749
!! The coefficient for Omega, D_Omega, and the chemical species, D_chim, are then constructed
!! from these coefficients
!-----------------------------------------------------------------------
subroutine coedif
!-----------------------------------------------------------------------
  use const,only: Msol,cst_G,cst_a,cst_c,Lsol
  use inputparam,only: iout,rapcrilim,icoeff,igamma,iadvec,istati,iledou,irot,fenerg,itminc,&
                       richac,xcn,imagn,add_diff
  use caramodele,only: inum,gms,glm,gls,hh6,nwmd
  use equadiffmod,only: iter,jterma
  use strucmod,only: m,q,pb,rb,tb,sb,zensi,Nabla_rad,Nabla_ad,delt,opac,rho,Nabla_mu,r,gravi,H_P
  use rotmod,only: omegi,dlodlr,xldoex,condbe,thext1,do1dr,thext2,vcirc,xmeg,ur1,gtilde
  use magmod,only: D_mago,D_magx,D_circh,Mag_diff
  use convection,only: rechzco,nzcon,nxzcon,iconra
  use diffadvmod,only: D_h,ucicoe,vcicoe,ursmooth,mtu,npasr,D_shear,D_conv,D_eff,Richardson,K_ther
  use geomod, only: rpsi_min,rpsi_max,geocalc
  use timestep,only: dzeit
  use advection,only: gbar,gtilgm
  use nagmod,only: c02agf
  use SmallFunc,only: neg_root
  use inputparam,only: mri !Adam MRI modification, fmu=0.05 or 1, for comparaison
  use inputparam,only: fmu


  implicit none

  integer, parameter:: nnrimax=120
  integer:: i,n,nn,ifirst,inosem,isemin,iconv,ideb,ifin,ncnc,ifail,jpos,ndegre,nroot,nn3,nuci,nxi,nli,ncci,ncce, &
            ncozo,iDh,iDshear,nnri,numricha
  integer, dimension(nnrimax):: n1r,n2r

  real(kindreal), parameter:: xpgam=0.10d0,xconv=1.5d0
  real(kindreal):: zwi1,xpsi,xgpsi,bcbcbc,dedede,fgfgfg,ababab,dshde,xnadm,richa,deltaR,dddccc,xfconv,vconv,dconml, &
     delna,delmu,croch1,delsh,aa0,aa1,aa2,aa3,xgam,gampol,dshun,rhom,xmst,xlumi,ura,urb,adramu,urc,xura,vmerid,xalpha, &
     xjojo,Cm,xbeta,xnut1,xnut2,xnut3,dr1,dr3,dr2,dU1,dU2,urn=0.d0,vrn,dmaxsh,dmaxef, &
     dbletimestep, bnmu, bnte !Adam added bnmu, bnte 
  real(kindreal), dimension(0:2):: apol2
  real(kindreal), dimension(0:3):: apol3
  real(kindreal), dimension(6):: www2
  real(kindreal), dimension(8):: www3
  real(kindreal), dimension(nnrimax):: rricha,drricha,domricha
  real(kindreal), dimension(ldi):: dV_z,Urho,D_sheardyn,admu,Urho_slope,lum,N_ad,N_mu,N_om,A_bc,B_bc,C_bc, &
     delta_bc,D_bcp,D_bcm,lambdab,mag_resist,etask,D_mri !Adam added lambdab, mag_resist, etask,D_mri,qmin
  real(kindreal), dimension(2,2):: zero2
  real(kindreal), dimension(2,3):: zero3
  real(kindreal) :: qmin_loc

  logical, parameter:: scale=.true.
!-----------------------------------------------------------------------
  gmsu = gms*Msol
  M_r(1:m) = (1.0d0-exp(q(1:m))) * gmsu
  lum(1:m) = (exp(sb(1:m))-1.d0)*exp(hh6)
  zwi1 = 1.d0/(exp(sb(1))-1.0d0)

  if (rapcrilim > 1.d-5) then
    nnri=nnrimax
  else
    nnri=80
  endif

  if (inum == 0) then
    dbletimestep = 1.0d0+xcn
  else
    dbletimestep = 2.d0
  endif

! choix des prescriptions de rotation Dh et Dshear:
  write(3,*) 'DANS COEDIF'

  iDh = mod(icoeff,10)
  iDshear = (icoeff-iDh)/10
  if (inum == 0 .and. iter == 1) then
    write(3,*)'PRESCRIPTIONS DE ROTATION: icoeff:',icoeff
    select case(icoeff)
      case(11)
        write(3,*)'    Dh: Zahn 1992'
        write(3,*)'    Dshear: Maeder 1997'
      case(12)
        write(3,*)'    Dh: Maeder 2003'
        write(3,*)'    Dshear: Maeder 1997'
      case(13)
        write(3,*)'    Dh: Mathis 2004'
        write(3,*)'    Dshear: Maeder 1997'
      case(21)
        write(3,*)'    Dh: Zahn 1992'
        write(3,*)'    Dshear: Talon & Zahn 1997'
      case(22)
        write(3,*)'    Dh: Maeder 2003'
        write(3,*)'    Dshear: Talon & Zahn 1997'
      case(23)
        write(3,*)'    Dh: Mathis 2004'
        write(3,*)'    Dshear: Talon & Zahn 1997'
      case(31)
        write(3,*)'    Dh: Zahn 1992'
        write(3,*)'    Dshear: Maeder 2013'
      case(32)
        write(3,*)'    Dh: Maeder 2003'
        write(3,*)'    Dshear: Maeder 2013'
      case(33)
        write(3,*)'    Dh: Mathis 2004'
        write(3,*)'    Dshear: Maeder 2013'
      case default
       stop 'Bad ICOEFF choice ! Must be 11,12,13,21,22,23,31,32 or 33.'
    end select
    write(3,*)'iDh:',iDh,' iDshear:',iDshear
  endif

  if (igamma==0 .and. (iadvec==0.or.jterma==1.or.istati==1)) then
    call gbar
    call gtilgm
  endif

  if (iledou == 1) then
    admu(1:m)=Nabla_ad(1:m)+1.d0/delt(1:m)*Nabla_mu(1:m)
  else
    admu(1:m)=Nabla_ad(1:m)
  endif

! calcul de la gravite
  gravi(m)=0.0d0

  do n=1,m-1
   xpsi=exp(rb(n))/((cst_G*M_r(n))/(omegi(n)*omegi(n)))**(1.0d0/3.0d0)
   if (xpsi >= rpsi_min .and. irot == 1) then
     if (xpsi > rpsi_max) then
       xpsi=0.9999999999d0*rpsi_max
     endif
! ancien appel: call geograv(xpsi,xgpsi)
     call geocalc(xpsi,xgpsi,1)
     gravi(n)=cst_G*M_r(n)*omegi(n)*omegi(n)*omegi(n)*omegi(n)
     gravi(n)=gravi(n)**(1.0d0/3.0d0)*xgpsi
   else
     gravi(n)=cst_G*M_r(n)/(exp(rb(n))*exp(rb(n)))
   endif
  enddo

! Calcule dans toute l'etoile de
! 1) l'echelle de pression H_P
! 2) du coefficient de diffusion radiative K_ther
! 3) du gradient de vitesse verticale dV_z
  do n=1,m
   if (gravi(n) /= 0.d0) then
     H_P(n)=exp(pb(n))/(exp(rho(n))*gravi(n))
   endif
  enddo
  K_ther(1:m)=4.d0*cst_a*cst_c*exp(4.d0*tb(1:m))*Nabla_ad(1:m)/(3.0d0*opac(1:m)*exp(rho(1:m))*exp(pb(1:m))*delt(1:m))
! 9 pi/32 = 0.8836: resultat de int_0^pi (sin^4(theta)dtheta)
!                              /int_0^pi (sin^3(theta)dtheta)
  dV_z(1:m)=abs((9.d0*pi/32.d0)*omegi(1:m)*dlodlr(1:m))

  H_P(m)=H_P(m-1)
! calcul de (3/(8pi) 1/(alpha (9 pi/32)^2 K rho))^(1/3)*Nad^{2/3}/r^2
! en fitm (couche 1)
  bcbcbc=gravi(1)*delt(1)/H_P(1)
  dedede=Nabla_ad(1)+1.d0/delt(1)*Nabla_mu(1)-Nabla_rad(1)
  fgfgfg=bcbcbc*dedede
  ababab=3.d0/(8.d0*pi*(9.d0*pi/32.d0)**2.d0)*fgfgfg/(fenerg*K_ther(1)*exp(rho(1)))
  ababab=neg_root(ababab,1.0d0/3.0d0)
  condbe=ababab/(exp(rb(1))*exp(rb(1)))*neg_root(xldoex,1.d0/3.d0)
  thext1=2.0d0/3.0d0*(exp(rb(1))*exp(rb(1)))*omegi(1)/gravi(1)*do1dr
  thext2=2.0d0/3.0d0*(exp(rb(1))*exp(rb(1)))*omegi(1)/gravi(1)*condbe

! critere de Richardson
! ici on verifie ici que le systeme est stable
! envers les instabilites de shear dynamiques
  do n=1,m
   Richardson(n)=0.0d0
   if ((n<m-1) .and. (zensi(n) <= 0.0d0 .and. Nabla_rad(n) >= admu(n))) then
   ! The sign of zensi at the centre is determined by Nabla_ad-Nabla_rad one layer above.
   ! There might be a descrepancy, hence the control only above m-1.
     if (itminc /= 1) then
       rewind(222)
       write (222,'(i7,a,i5)') nwmd,': problem in coedif l.204 in layer ',n
       write(*,*) 'problem in coedif l.204',n
       stop
     endif
   endif
   if (Nabla_rad(n) < admu(n)) then
     dshde=gravi(n)*delt(n)/H_P(n)
     xnadm=Nabla_ad(n)+1.d0/delt(n)*Nabla_mu(n)-Nabla_rad(n)
     richa=dshde*xnadm-richac*dV_z(n)*dV_z(n)
     if (richa < 0.0d0) then
       Richardson(n)=dshde*xnadm/(dV_z(n)*dV_z(n))-richac
     endif
   endif
  enddo

  D_sheardyn(1:m)=0.0d0
  D_conv(1:m)=0.0d0
  D_shear(1:m)=0.0d0
  D_h(1:m)=0.0d0
  D_eff(1:m)=0.0d0
  D_chim(1:m)=0.0d0
  D_circh(1:m)=0.0d0
  D_bcp(1:m)=0.0d0
  D_bcm(1:m)=0.0d0

  do n=1,m
   if (iadvec==1 .and. jterma==0) then
     ucicoe(n)=ursmooth(n)/dbletimestep
     vcicoe(n)=vcirc(n)
   endif
  enddo

  numricha=0
  n1r(1:nnri)=0
  n2r(1:nnri)=0
  rricha(1:nnri)=0.0d0
  drricha(1:nnri)=0.0d0
  domricha(1:nnri)=0.0d0

  do n=m,1,-1
   if (Richardson(n) < 0.00d0) then
     if (n == m) then
       numricha=numricha+1
       if (numricha > nnri) then
         stop 'prob. coedif richa l.250'
       endif
       n1r(numricha)=n
     else
       if (Richardson(n+1) == 0.0d0) then
         numricha=numricha+1
         if (numricha > nnri) then
           stop 'prob. coedif richa l.250'
         endif
         n1r(numricha)=n
       endif
     endif
   else ! xricha(n) >= 0
     if (n < m) then
       if (Richardson(n+1) < 0.0d0) then
         n2r(numricha)=n+1
       endif
     endif
   endif
  enddo
  if (Richardson(1) < 0.00d0) then
    n2r(numricha)=1
  endif

  write(3,*) 'numricha=',numricha,' <? nnri=',nnri
  do n=1,numricha
   drricha(n)=-(exp(r(min(n1r(n)+1,m)))+exp(r(n1r(n))))/2.d0+(exp(r(max(n2r(n)-1,1)))+exp(r(n2r(n))))/2.d0
   rricha(n)=(exp(r(n1r(n)))+exp(r(n2r(n))))/2.d0
   domricha(n)=(omegi(min(n1r(n)+1,m))+omegi(n1r(n)))/2.d0-(omegi(max(n2r(n)-1,1))+omegi(n2r(n)))/2.d0
   if (n1r(n) == m .and. n2r(n) == m) then
     drricha(n)=0.d0
   endif
   do nn= n1r(n),n2r(n),-1
    D_sheardyn(nn)=abs(1.d0/3.d0*rricha(n)*domricha(n)*drricha(n))
   enddo
  enddo

  ifirst=0
  inosem=0
  isemin=0
! recherche des zones convectives
  call rechzco

! etendue spatiale de la zone convective
  do iconv=1,nzcon
   ideb=nxzcon(2*iconv-1)
   ifin=nxzcon(2*iconv)
   deltaR=exp(r(ifin))-exp(r(ideb))
   dddccc=10.d0*deltaR*deltaR/dzeit
   do ncnc=ifin,ideb
    D_conv(ncnc)=dddccc
   enddo
  enddo
!**********************************************
! calcul de Dmago et de Dmagx 28 janvier 2003
  if (imagn == 1) then
    call Mag_diff(m,zensi,H_P,gravi,Nabla_mu,delt,Nabla_rad,Nabla_ad,rb,omegi,tb,dlodlr,rho,K_ther)
  endif
!**********************************************
  do n=1,m
! ZONE CONVECTIVE
   if (zensi(n) > 0.0d0) then
     if (lum(n) < 0.d0) then
       write(3,*) 'Neg. luminosity in coediff, abs value token.'
     endif
     xfconv=abs(lum(n))/(4.d0*pi*exp(rb(n))*exp(rb(n)))
     vconv=(0.25d0*xfconv*gravi(n)*Nabla_ad(n)*xconv*H_P(n)/exp(pb(n)))**(1.0d0/3.0d0)
     dconml=1.0d0/3.0d0*xconv*H_P(n)*vconv
     D_conv(n)=min(dconml,D_conv(n))
     cycle
   else   ! zensi =< 0
! AUTRES ZONES
     if (igamma == 1) then
! cas Dshear de Maeder 1997 (Papier II, Eq. 5.30)
! avec calcul des Gammas (terme 4(Nabla'-Nabla))
       delna=Nabla_ad(n)-Nabla_rad(n)
       delmu=Nabla_mu(n)/delt(n)
       croch1=dV_z(n)*dV_z(n)
       if (gravi(n) /= 0.0d0) then
         delsh=fenerg*H_P(n)/(gravi(n)*delt(n))*croch1
       else
         delsh=0.0d0
       endif
       aa0=12.d0*delmu
       aa1=2.d0*(delna+delmu)-6.d0*delsh
       aa2=2.d0*(delna+delmu)+4.d0*delna-delsh
       aa3=-delsh
       ifail=0
       jpos=0

       if (aa0 /= 0.d0) then
         ndegre=3
         apol3(0)=aa0
         apol3(1)=aa1
         apol3(2)=aa2
         apol3(3)=aa3
         call c02agf(apol3,ndegre,scale,zero3,www3,ifail)
         nroot=1
         do while (nroot <= ndegre)
          if (zero3(2,nroot) == 0.d0) then
            if (zero3(1,nroot) > 0.d0) then
              xgam=zero3(1,nroot)
              D_shear(n)=2.d0*K_ther(n)*zero3(1,nroot)
              jpos=jpos+1
            endif
          endif
          nroot=nroot+1
         enddo
       else   ! if aa0==0
         if (aa1 /= 0) then
           ndegre=2
           apol2(0)=aa1
           apol2(1)=aa2
           apol2(2)=aa3
           call c02agf(apol2,ndegre,scale,zero2,www2,ifail)
           nroot=1
           do while (nroot <= ndegre)
            if (zero2(2,nroot) == 0.d0) then
              if (zero2(1,nroot) > 0.d0) then
                xgam=zero2(1,nroot)
                D_shear(n)=2.d0*K_ther(n)*zero2(1,nroot)
                jpos=jpos+1
              endif
            endif
            nroot=nroot+1
           enddo
         else   ! if aa1==0
           if (aa2 /= 0.d0) then
             gampol=delsh/aa2
             if (gampol >= 0.d0) then
               D_shear(n)=2.d0*K_ther(n)*gampol
               jpos=jpos+1
             endif
           else
             jpos=-1
           endif   ! aa2
         endif   ! aa1
       endif   ! aa0

       cycle

     else   ! if igamma = 0
! cas Dshear de Talon et Zahn 1997
! fraction de l'energie du schear utilisee pour le melange=fenerg
! formule pour Dschear corrigee de l'effet discute dans
! Talon et Zahn AA 317, 749
       if (iter == 1 .and. itminc /= 1 .and. iDshear /=3) then
! cas Dshear de Maeder 1997 (Papier II, Eq. 5.30)
! sans calcul des Gammas (terme 4(Nabla'-Nabla))
         dshun=K_ther(n)*H_P(n)/(gravi(n)*delt(n))
         dshde=gravi(n)*delt(n)/H_P(n)
         xnadm=Nabla_ad(n)+1.d0/delt(n)*Nabla_mu(n)-Nabla_rad(n)
         croch1=dV_z(n)*dV_z(n)
         if (xnadm /= 0.0d0) then
           D_shear(n)=dshun/xnadm*(fenerg*croch1)
         else
           D_shear(n)=0.0d0
         endif
       else
         if (istati == 0) then
           if (iadvec == 0 .or. jterma == 1) then
! calcul de U simplifie
! calcul de rhom
             if (n == m) then
               rhom= exp(rho(n))
             else
               rhom=exp(glm)*(1.0d0-exp(q(n)))/(4.d0*pi/3.d0*exp(3.d0*rb(n)))
             endif
             xmst=exp(glm)*(1.0d0-exp(q(n)))*(1.0d0-omegi(n)*omegi(n)/(2.d0*pi*cst_G*rhom))
! L/(M* X g)
             xlumi=(exp(sb(n))-1.0d0)*zwi1*gls*Lsol
             ura=xlumi/(xmst*xmeg(n))
! P/(rho Cp T)
             urb=Nabla_ad(n)/delt(n)
! 1./(Ad-nab)
             adramu=Nabla_ad(n)-Nabla_rad(n)+Nabla_mu(n)/delt(n)
             if (adramu > 0.0d0) then
               urc=1.d0/(adramu)
             else
               urc=0.0d0
             endif
             ur1(n)=ura*urb*urc

! calcul de Ur simplifie
             xura=2.d0*gtilde(n)*(1.0d0-omegi(n)*omegi(n)/(2.d0*pi*cst_G*exp(rho(n))))
             ucicoe(n)=ur1(n)*xura
             vmerid=ucicoe(n)
           else ! iadvec == 1 and jterma /= 1
             vmerid=vcirc(n)
           endif
           xalpha=1.0d0+0.5d0*dlodlr(n)
           xjojo=abs(2.d0*vmerid-xalpha*ucicoe(n))

           select case(iDh)
             case (1)
! Dh Zahn 92
               D_h(n)=exp(rb(n))*xjojo
             case (2)
! Dh de Maeder (2003) A&A 399, 263
               Cm=1.0d0
               D_h(n)=0.002d0*Cm*exp(rb(n))**(4.0d0/3.0d0)*omegi(n)**(1.0d0/3.0d0)*abs(vmerid)**(1.0d0/3.0d0) &
                         *xjojo**(1.0d0/3.0d0)
             case (3)
! Dh de Mathis et al. (2004) A&A 425, 243
               xbeta=1.5d-6
               xnut1=sqrt(xbeta/10.d0)
               xnut2=sqrt(exp(rb(n))*exp(rb(n))*omegi(n))
               xnut3=sqrt(exp(rb(n))*xjojo)
               D_h(n)=xnut1*xnut2*xnut3
             case default
               stop 'Bad ICOEFF choice ! Must be 11,12,13,21,22,23,31,32 or 33.'
           end select

           select case(iDshear)
             case (1)
! Dshear Maeder 1997
               dshun=K_ther(n)*H_P(n)/(gravi(n)*delt(n))
               xnadm=Nabla_ad(n)+1.d0/delt(n)*Nabla_mu(n)-Nabla_rad(n)
             case (2)
! Dshear Talon & Zahn 1997
               dshun=(K_ther(n)+D_h(n))*H_P(n)/(gravi(n)*delt(n))
               if (D_h(n) /= 0.0d0) then
                 xnadm=Nabla_ad(n)+1.d0/delt(n)*Nabla_mu(n)*(1.d0+K_ther(n)/D_h(n))-Nabla_rad(n)
               else
                 xnadm=0.0d0
               endif
             case (3)
! Dshear Maeder 2013
               if (H_P(n) /= 0.0d0) then
                 N_ad(n) = gravi(n)*delt(n)/H_P(n) * (Nabla_ad(n)-Nabla_rad(n))
                 N_mu(n) = gravi(n)/H_P(n) * Nabla_mu(n)
               endif
               N_om(n) = 1.767146d0 * omegi(n)*omegi(n) * (dlodlr(n) + 2.0d0)
               N_om(n) = N_om(n) - 0.8d0 * richac * omegi(n)*omegi(n) * dlodlr(n)*dlodlr(n)

               A_bc(n) = N_ad(n) + N_mu(n) + N_om(n)
               B_bc(n) = N_ad(n)*D_h(n) + N_mu(n)*(K_ther(n)+D_h(n)) + N_om(n)*(K_ther(n)+2.0d0*D_h(n))
               C_bc(n) = N_om(n)*(D_h(n)*K_ther(n) + D_h(n)*D_h(n))
               delta_bc(n) = B_bc(n)*B_bc(n) - 4.0d0*A_bc(n)*C_bc(n)
               if (delta_bc(n) > 0.0d0) then
                 D_bcm(n) = abs((-B_bc(n)-sqrt(delta_bc(n)))/(A_bc(n))) !racine -, qui ne semble pas physique
                 D_bcp(n) = abs((-B_bc(n)+sqrt(delta_bc(n)))/(A_bc(n))) !racine +
               endif
               !D_shear(n) = D_bcm(n) !racine moins, semble moins physique
               D_shear(n) = D_bcp(n) !semble plus realiste
             case default
               stop 'Bad ICOEFF choice ! Must be 11,12,13,21,22,23,31,32 or 33.'
           end select

           dshde=gravi(n)*delt(n)/H_P(n)
           croch1=dV_z(n)*dV_z(n)
           if (iDshear /= 3) then
             if (xnadm /= 0.0d0) then
               D_shear(n)=dshun/xnadm*(fenerg*croch1)
             else
               D_shear(n)=0.0d0
             endif
           endif
         else   ! istati = 1

! calcul de U simplifie
! calcul de rhom
           if (n == m) then
             rhom= exp(rho(n))
           else
             rhom=exp(glm)*(1.0d0-exp(q(n)))/(4.d0*pi/3.d0*exp(3.d0*rb(n)))
           endif
           xmst=exp(glm)*(1.0d0-exp(q(n)))*(1.0d0-omegi(n)*omegi(n)/(2.d0*pi*cst_G*rhom))
! L/(M* X g)
           xlumi=(exp(sb(n))-1.d0)*zwi1*gls*Lsol
           ura=xlumi/(xmst*xmeg(n))
! P/(rho Cp T)
           urb=Nabla_ad(n)/delt(n)
! 1./(Ad-nab)
           adramu=Nabla_ad(n)-Nabla_rad(n)+Nabla_mu(n)/delt(n)
           if (adramu > 0.0d0) then
             urc=1.d0/(adramu)
           else
             urc=0.0d0
           endif
           ur1(n)=ura*urb*urc

! calcul de Ur simplifie
           xura=2.d0*gtilde(n)*(1.0d0-omegi(n)*omegi(n)/(2.d0*pi*cst_G*exp(rho(n))))
           ucicoe(n)=ur1(n)*xura
           D_h(n)=abs(exp(rb(n))*ucicoe(n))
! calcul de dshear
           select case(iDshear)
             case (1)
! Dshear Maeder 1997
               dshun=K_ther(n)*H_P(n)/(gravi(n)*delt(n))
               xnadm=Nabla_ad(n)+1.d0/delt(n)*Nabla_mu(n)-Nabla_rad(n)
             case (2)
! Dshear Talon & Zahn 1997
               dshun=(K_ther(n)+D_h(n))*H_P(n)/(gravi(n)*delt(n))
               xnadm=Nabla_ad(n)+1.d0/delt(n)*Nabla_mu(n)*(1.d0+K_ther(n)/D_h(n))-Nabla_rad(n)
             case (3)
! Dshear Maeder 2013
               if (H_P(n) /= 0.0d0) then
                 N_ad(n) = gravi(n)*delt(n)/H_P(n) * (Nabla_ad(n)-Nabla_rad(n))
                 N_mu(n) = gravi(n)/H_P(n) * Nabla_mu(n)
               endif
               N_om(n) = 1.767146d0 * omegi(n)*omegi(n) * (dlodlr(n) + 2.0d0)
               N_om(n) = N_om(n) - 0.8d0 * richac * omegi(n)*omegi(n) * dlodlr(n)*dlodlr(n)

               A_bc(n) = N_ad(n) + N_mu(n) + N_om(n)
               B_bc(n) = N_ad(n)*D_h(n) + N_mu(n)*(K_ther(n)+D_h(n)) + N_om(n)*(K_ther(n)+2.0d0*D_h(n))
               C_bc(n) = N_om(n)*(D_h(n)*K_ther(n) + D_h(n)*D_h(n))
               delta_bc(n) = B_bc(n)*B_bc(n) - 4.0d0*A_bc(n)*C_bc(n)
               if (delta_bc(n) > 0.0d0) then
                 D_bcm(n) = abs((-B_bc(n)-sqrt(delta_bc(n)))/(A_bc(n))) !racine -, qui ne semble pas physique
                 D_bcp(n) = abs((-B_bc(n)+sqrt(delta_bc(n)))/(A_bc(n))) !racine +
               endif
               !D_shear(n) = D_bcm(n) !racine moins, semble moins physique
               D_shear(n) = D_bcp(n) !plus realiste
             case default
               stop 'Bad ICOEFF choice ! Must be 11,12,13,21,22,23,31,32 or 33.'
           end select

           dshde=gravi(n)*delt(n)/H_P(n)
           croch1=dV_z(n)*dV_z(n)
           if (iDshear /= 3) then
             if (xnadm /= 0.0d0) then
               D_shear(n)=dshun/xnadm*(fenerg*croch1)
             else
               D_shear(n)=0.0d0
             endif
           endif
           if (D_shear(n) < 0.d0) then
             D_shear(n)=abs(D_shear(n))
           endif
         endif   ! istati
       endif   ! iter & itminc
     endif   !if igamma
   
  
  !   !Adam Implementation of MRI and advection, turned off for MRI+TS implementation
  !    if ((mri==1 .and. imagn==0)) then !If one wants to compute the MRI, not if the instabilit is active at point n !!! 
  !       if (H_P(n) /= 0.0d0) then
  !         ! bnmu: N_mu^2 (Paper 1, Eq. 1)
  !        bnmu=gravi(n)*Nabla_mu(n)/H_P(n)
  !         ! bnte: N_T^2 (Paper 1, Eq. 2)
  !         bnte=gravi(n)*delt(n)/H_P(n)*abs(Nabla_rad(n)-Nabla_ad(n))
  !       else
  !         bnmu=0.0d0
  !         bnte=0.0d0
  !       endif
  !       ! lambdab : Ln(Lambda)=-12.7+ln(T)-0.5ln(rho) as in Paper by Wheeler et al. 2015 eq (5), mag_resist at rest
  !       lambdab(n)=-12.7d0+tb(n)-0.5d0*rho(n)
  !       mag_resist(n)=5.2d0*(10.d0**11.d0)*lambdab(n)*exp(-1.5*tb(n))
    

  !       ! etask: eta/K
  !       etask(n)=mag_resist(n)/K_ther(n)

  !       ! MRI diffusion ceof as in Paper by Wheeler et al. 2015 eq (13), note that DmagO=DmagX in this case
  !       D_mri(n)= 0.02d0*abs(dlodlr(n))*omegi(n)*exp(rb(n))*exp(rb(n))

  !       qmin_loc=abs(-(etask(n)*bnte+fmu*bnmu)/(2.0d0*omegi(n)*omegi(n))) !MRI minimum shear to activate 
  !       if (abs(dlodlr(n)) > qmin_loc .and. (abs(dlodlr(n))<4) ) then ! ATTENTION this condition does not ask if Omega>alven, for simplicity alven not computed and this cond is always verified
  !         D_shear(n)=D_shear(n)+MIN(D_mri(n),10.d0**12.d0)
  !         qmin(n) = 1.
  !       else
  !         qmin(n) = -1.
  !       endif !q>qmin
  !     endif !mri subrout
    endif    ! zensi
  enddo



! calcul de dcirch
  if (iadvec == 0 .and. imagn == 1) then
    do n=1,m
     if (dlodlr(n) == 0.0d0 .or. zensi(n) > 0.0d0) then
       D_circh(n)=0.0d0
     else
      D_circh(n)=  abs(exp(rb(n))*ucicoe(n))  
      ! D_circh(n)=0.0d0
     endif
    enddo
  endif

! calcul de Deff
  if (istati == 1) then
    Urho(1:m)=exp(rho(1:m))*ursmooth(1:m)

    do i=2,m-1
     if (D_shear(i) >  0.0d0) then
       dr1=exp(rb(i-1))-exp(rb(i))
       dr2=exp(rb(i+1))-exp(rb(i))
       dr3=exp(rb(i-1))-exp(rb(i+1))
       dU2=(Urho(i)-Urho(i-1))
       dU1=(Urho(i)-Urho(i+1))
       Urho_slope(i)=dU2*dr1/(dr2*dr3) - dU1*dr2/(dr1*dr3)
     else   ! dshear <= 0
       Urho_slope(i)=0.0d0
     endif   ! dshear
    enddo
    Urho_slope(1)=Urho_slope(2)
    Urho_slope(m)=Urho_slope(m-1)

    do n=1,m
     if (D_shear(n) > 0.0d0) then
       vmerid=ursmooth(n)/3.d0+exp(rb(n))/(6.d0*exp(rho(n)))*Urho_slope(n)
       D_eff(n)=exp(rb(n))*ursmooth(n)*ursmooth(n)/30.d0
       xalpha=1.d0+0.5d0*dlodlr(n)
       xjojo=abs(2.d0*vmerid-xalpha*ursmooth(n))
       if (xjojo == 0.d0) then
         D_eff(n)=D_eff(n)/ursmooth(n)
       else
         D_eff(n)=D_eff(n)/xjojo
       endif
     else
       D_eff(n)=0.0d0
     endif   ! dshear
    enddo
  else   ! istati = 0

    if (iadvec == 1 .and. jterma == 0) then
      if (verbose) then
        write(*,*) ' iter=',iter,' itminc=',itminc
      endif
      if (iter /= 1 .or. itminc == 1) then
        nn3=npasr-mtu
        if (nn3 >= 3) then
          do n=mtu,npasr
           xalpha=1.d0+0.5d0*dlodlr(n)
           if (inum == 0) then
             urn=ursmooth(n)/(1.0d0+xcn)
           else if (inum > 0) then
             urn=ursmooth(n)/2.d0
           endif
           vrn=vcirc(n)
! Deff=1/30 (ru)^2/Dh
           D_eff(n)=1.d0/30.d0*exp(rb(n))*exp(rb(n))*urn*urn
           if (D_h(n) /= 0.0d0) then
             D_eff(n)=D_eff(n)/D_h(n)
           else
             D_eff(n)=0.0d0
           endif
          enddo
        endif
      endif

! pour calcul de Deff losque IADVEC=0 et IMAGN=1
    else   ! if iadvec = 0
      do n=1,m
       if (imagn == 1 .and. zensi(n) <= 0.0d0) then
         D_eff(n)=1.d0/30.d0*exp(rb(n))*exp(rb(n))*ucicoe(n)*ucicoe(n)
         if (D_h(n) /= 0.d0) then
           D_eff(n)=D_eff(n)/D_h(n)
         else
           D_eff(n)=0.0d0
         endif
       endif
      enddo
    endif   ! iadvec
  endif   ! istati

  if (D_conv(m-1) /= 0.0d0) then
    D_conv(m)=D_conv(m-1)
  endif
  D_shear(1)=D_shear(2)

  do n=1,m
   if (D_shear(n) > 0.0d0) then
! adjonction de maniere a eviter de trop grands coefficients de diff
     dmaxsh=1.0d+15
     dmaxef=1.0d+13
     D_shear(n)=min(D_shear(n),dmaxsh)
     D_eff(n)=min(D_eff(n),dmaxef)

     if (Richardson(n) < 0.0d0 .and. Richardson(n)+richac > 0.0d0) then
       D_shear(n)=D_shear(n)+D_sheardyn(n)
     endif
! si idifcon= 1, on peut entrer dans coedif meme si irot= 0:
     if (irot == 0) then
       D_shear(n)=0.0d0
       D_eff(n)=0.0d0
     endif

     D_chim(n)=D_shear(n)+D_eff(n)  !WARNING adam remove Deff No circulation mixing..
   else
     if (D_conv(n) < 0.0d0 .or. D_conv(n) > 1.0d99) then
       D_conv(n)=0.0d0
     endif
     D_chim(n)=D_conv(n)
   endif
  enddo
! Au bord des zones convectives on construit sur quelques coquilles
! un coefficient de diffusion tel que le gradient de ce coefficient
! soit doux
! 1) reperage du bord du noyau et de la base de l'enveloppe convective
  call rechzco
! Y a-t-il un noyau convectif ?
  if (iconra == 3 .or. iconra == 5) then
    nuci=nxzcon(2)
    if (iout /= 0) then
      nxi=nuci-iout
      do i=nxi,nuci
       D_chim(i)=D_chim(nxi)
      enddo
    endif
  endif
! Y a-t-il une enveloppe convective ?
  if (iconra == 3 .or. iconra == 4) then
    nuci=nxzcon(2*nzcon-1)
    if (iout /= 0) then
      nli=nuci+iout
      do i=nuci,nli
       D_chim(i)=D_chim(nli)
      enddo
    endif
  endif
! Y a-t-il d'autres zones convectives importantes ?
  if (iconra == 2 .and. nzcon >= 1) then
    do n=1,nzcon
! numero des couches des bords de la zone convective
! intermediaire n
     ncci=nxzcon(2*n-1)
     ncce=nxzcon(2*n)
! est-ce une zone convective suffisamment importante
     ncozo=ncci-ncce
     if (ncozo >= 10) then
       if (iout /= 0) then
         nli=ncci+iout
         nxi=ncce-iout
         if (nli > m) then
           nli=m
         endif
         if (nxi < 1) then
           nxi=1
         endif
         do i=ncci,nli
          D_chim(i)=D_chim(nli)
         enddo
         do i=nxi,ncce
          D_chim(i)=D_chim(nxi)
         enddo
       endif
     endif
    enddo
  endif

  if (iconra == 3 .and. nzcon >= 3) then
    do n=1,nzcon-2
! numero des couches des bords de la zone convective
! intermediaire n
     ncci=nxzcon(2+2*n-1)
     ncce=nxzcon(2+2*n)
! est-ce une zone convective suffisamment importante
     ncozo=ncci-ncce
     if (ncozo >= 10) then
       if (iout /= 0) then
         nli=ncci+iout
         nxi=ncce-iout
         if (nli > m) then
           nli=m
         endif
         if (nxi < 1) then
           nxi=1
         endif
         do i=ncci,nli
          D_chim(i)=D_chim(nli)
         enddo
         do i=nxi,ncce
          D_chim(i)=D_chim(nxi)
         enddo
       endif
     endif
    enddo
  endif   ! iconra

  if (iconra == 4 .and. nzcon >= 2) then
    do n=1,nzcon-1
! numero des couches des bords de la zone convective intermediaire n
     ncci=nxzcon(2*n-1)
     ncce=nxzcon(2*n)
! est-ce une zone convective suffisamment importante
     ncozo=ncci-ncce+1
     if (ncozo >= 10) then
       if (iout /= 0) then
         nli=ncci+iout
         nxi=ncce-iout
         if (nli > m) then
           nli=m
         endif
         if (nxi < 1) then
           nxi=1
         endif
         do i=ncci,nli
          D_chim(i)=D_chim(nli)
         enddo
         do i=nxi,ncce
          D_chim(i)=D_chim(nxi)
         enddo
       endif
     endif
    enddo
  endif

  if (iconra == 5 .and. nzcon >= 2) then
    do n=1,nzcon-1
! numero des couches des bords de la zone convective
! intermediaire n
     ncci=nxzcon(2+2*n-1)
     ncce=nxzcon(2+2*n)
! est-ce une zone convective suffisamment importante
     ncozo=ncci-ncce
     if (ncozo >= 10) then
       if (iout /= 0) then
         nli=ncci+iout
         nxi=ncce-iout
         if (nli > m) then
           nli=m
         endif
         if (nxi < 1) then
           nxi=1
         endif
         do i=ncci,nli
          D_chim(i)=D_chim(nli)
         enddo
         do i=nxi,ncce
          D_chim(i)=D_chim(nxi)
         enddo
       endif
     endif
    enddo
  endif
  if (imagn == 1) then
    D_chim(1:m) = D_chim(1:m) + D_magx(1:m)
  endif

! Coefficient de diffusion pour le transport du moment cinetique
! On doit prendre un coefficient de diffusion egal a deux fois le coefficient
! nominal lorsque IADVEC=1 car alors on alterne une fois sur deux la
! diffusion du moment angulaire et l'advection.
! Sinon on doit utiliser un coefficient de diffusion non multiplie par deux
! lorsque imagn eq 1 on utilise delta t et non 2*delta t, car on n'utilise pas
! la partie advective de l equation de transport du moment cinetique
  if (iadvec == 1 .and. imagn == 0) then
    D_Omega(1:m)=dbletimestep*(D_shear(1:m)+D_conv(1:m)) + add_diff
  else ! IADVEC= 0
    if (imagn == 0) then
      D_Omega(1:m)=D_shear(1:m)+D_conv(1:m) + add_diff
    else
      D_Omega(1:m)=D_shear(1:m)+D_conv(1:m)+D_mago(1:m)+D_circh(1:m) !Adam change put Deff instead of Dcirch
    endif
  endif ! IADVEC

  return

end subroutine coedif

!=======================================================================
subroutine couran
!-----------------------------------------------------------------------
  use inputparam,only: ipop3,idiff
  use abundmod,only: wx,wy,wxc12,wxo16,wxne20,wxne22,wxmg24,wxmg25,wxmg26,wxc13,wxn14,wxn15,wxo17,wxo18
  use strucmod,only: m,amu
  use diffadvmod,only: tdiff
  use timestep,only: dzeit

  implicit none

  integer:: i,ij,jj,nm,ielemen=0,nbshell=0,ncouch,nint,iint
  real(kindreal):: rapmax,rp,atm,tm,sous,dd
  real(kindreal), dimension(16):: rap
!-----------------------------------------------------------------------
  tdiff=dzeit

  do i=1,m-1
   nm = m-i+1
   dd=2.d0*D_moychim(nm)
   if (amu(i+1) == amu(i) .or. dd == 0.d0) then
     cycle
   endif
   rapmax=0.d0

   if (wx(i) == 0.d0) then
     rap(1)=wy(i+1)/max(wy(i),1.d-50)
     rap(2)=wxc12(i+1)/max(wxc12(i),1.d-50)
     rap(3)=1.d0
     rap(4)=1.d0
     rap(5)=wxo16(i+1)/max(wxo16(i),1.d-50)
     rap(6)=1.d0
     rap(7)=wxne20(i+1)/max(wxne20(i),1.d-50)
     rap(8)=wxne22(i+1)/max(wxne22(i),1.d-50)
     rap(9)=wxmg24(i+1)/max(wxmg24(i),1.d-50)
     rap(10)=wxmg25(i+1)/max(wxmg25(i),1.d-50)
     rap(11)=wxmg26(i+1)/max(wxmg26(i),1.d-50)
     ij=11
   else
     rap(1)=wx(i+1)/max(wx(i),1.d-50)
     rap(2)=wy(i+1)/max(wy(i),1.d-50)
     rap(3)=1.d0
     rap(5)=1.d0
     rap(6)=wxc12(i+1)/max(wxc12(i),1.d-50)
     rap(7)=wxc13(i+1)/max(wxc13(i),1.d-50)
     rap(8)=wxn14(i+1)/max(wxn14(i),1.d-50)
     rap(9)=wxn15(i+1)/max(wxn15(i),1.d-50)
     rap(10)=wxo16(i+1)/max(wxo16(i),1.d-50)
     rap(11)=wxo17(i+1)/max(wxo17(i),1.d-50)
     rap(4)=wxo18(i+1)/max(wxo18(i),1.d-50)
     ij=11
     if (wx(1) <= 0.d0 .or. ipop3 /= 0) then
       rap(12)=wxne20(i+1)/max(wxne20(i),1.d-50)
       rap(13)=wxne22(i+1)/max(wxne22(i),1.d-50)
       rap(14)=wxmg24(i+1)/max(wxmg24(i),1.d-50)
       rap(15)=wxmg25(i+1)/max(wxmg25(i),1.d-50)
       rap(16)=wxmg26(i+1)/max(wxmg26(i),1.d-50)
       ij=16
     endif
   endif

   do jj=1,ij
    rp=abs((rap(jj)-1.d0)/(rap(jj)+1.d0))
    rapmax=max(rapmax,rp)
    if (rapmax == rp) then
      ielemen=jj
      nbshell=i
    endif
   enddo

   if (rapmax <= 0.d0) then
     cycle
   endif
   atm=2.d0*log(dra(i))-log(dd)-log(rapmax)
   tm=exp(atm)
   tdiff=min(tdiff,tm)
   if (tdiff == tm) then
     ncouch=i
   endif
  enddo

  sous=dzeit/tdiff

  if (abs(sous) > 1.d9 .or. sous < 0.d0) then
    if (verbose) then
      print*,'sous= ',sous
    endif
    sous=100.d0
  endif
  nint=int(sous)
  iint=nint

  if (nint >= 99 .or. nint <= 0) then
    nint=99
    write(3,'(1x,a12,2i10,e13.5,e13.5,a2,i4,a2,i4)') ' NINT LIMITE',iint,nint,sous,tdiff,'el',ielemen,'sh',nbshell
  endif

  tdiff=dzeit/(nint*idiff)
  nwpas=nint*idiff
  write(3,'(a21,e13.5,a9,e13.5,a1)') ' temps de diffusion =',tdiff,' (dzeit= ',dzeit,')'
  write(3,'(1x,2(a7,i2),/)') ' idiff=',idiff,' nwpas=',nwpas

  return

end subroutine couran

!=======================================================================
subroutine courom
!-----------------------------------------------------------------------
  use inputparam,only: idiff
  use strucmod,only: m,rb
  use rotmod,only: vvomeg
  use diffadvmod,only: tdiff
  use timestep,only: dzeit

  implicit none

  integer:: i,ncouch,nint,iint
  real(kindreal):: dd,rapmax,rap,rp,atm,tm,sous
!-----------------------------------------------------------------------
  tdiff=dzeit
  do i=1,m-1
   dra(i)=exp(rb(i))-exp(rb(i+1))
  enddo

  do i=1,m-1
   dd=2.d0*D_moyOm(i)
   if (vvomeg(i+1) /= vvomeg(i) .and. dd /= 0.d0) then
     rapmax=0.d0
     rap=vvomeg(i+1)/max(vvomeg(i),1.0d-50)
     rp=(rap-1.d0)/(rap+1.d0)
     rp=abs(rp)
     rapmax=max(rapmax,rp)
     if (rapmax /= 0.d0) then
       atm=2.d0*log(dra(i))-log(dd)-log(rapmax)
       tm=exp(atm)
       tdiff=min(tdiff,tm)
       if (tdiff == tm) then
         ncouch=i
       endif
     endif
   endif
  enddo

  sous=dzeit/tdiff
  if (abs(sous) > 1.d9 .or. sous < 0.d0) then
    if (verbose) then
      write(*,*) 'dzeit/tdiff= ',sous
    endif
    sous=100.d0
  endif
  nint=int(sous)
  iint=nint

  if (nint >= 99 .or. nint <= 0) then
    nint=99
    write(3,'(a,2i10)') ' NINT LIMITE',iint,nint
    write(3,*) 'courom_impli.f, l.67'
  endif

  tdiff=dzeit/(nint*idiff)
  nwpas=nint*idiff
  write(3,'(a,e13.5)') 'temps de diffusion =',tdiff
  write(3,'(a,i2,/)') 'fraction courant=',idiff

  return

end subroutine courom

!=======================================================================
subroutine diffbr
!-----------------------------------------------------------------------
  use const,only: Msol
  use inputparam,only: z,ipop3,phase,ibasnet,ialflu,imagn,idifcon
  use abundmod,only: x,y3,y,xc12,xc13,xc14,xn14,xn15,xo16,xo17,xo18,xf18,xf19,xne20,xne21,xne22,xna23,xmg24,xmg25,xmg26,xal26, &
    xal27,xsi28,xprot,xneut,xbid,xbid1,wx,wy3,wy,wxc12,wxc13,wxc14,wxn14,wxn15,wxo16,wxo17,wxo18,wxf18,wxf19,wxne20,wxne21, &
    wxne22,wxna23,wxmg24,wxmg25,wxmg26,wxal26g,wxal27,wxsi28,wxprot,wxneut,wxbid,wxbid1,vvx,vvy3,vvy,vvxc12,vvxc13,vvxc14,vvxn14, &
    vvxn15,vvxo16,vvxo17,vvxo18,vvxf18,vvxf19,vvxne20,vvxne21,vvxne22,vvxna23,vvxmg24,vvxmg25,vvxmg26,vvxal26g,vvxal27,vvxsi28, &
    vvxprot,vvxneut,vvxbid,vvxbid1,epsc,epsy,nbelx,wabelx, &
    vvabelx,zabelx,fnucdif,mbelx,abelx,nbzel,nbael
  use strucmod,only: m,rb,t,rho
  use diffadvmod,only: tdiff,jdiff
  use convection,only: jzint,izc
  use timestep,only: dzeit
  use SmallFunc,only: tridiago
  use energy,only: netburning

  implicit none

  integer::i,n,nm,ii,jk,jbid,i1,i2,nnn
  real(kindreal):: sumx,sumy,sumy3,sumc12,sumc13,sumn14,sumn15,sumo16,sumo17,sumo18,sumne20,sumne22,summg24,summg25,summg26, &
    sumc14,sumf18,sumf19,sumne21,sumna23,sumal26g,sumal27,sumsi28,sumneut,sumprot,sumbid,sumbid1,r32,rho32,dr32,cc,aa, &
    bb,whc,wyc,wy3c,wxc12c,wxc13c,wxn14c,wxn15c,wxo16c,wxo17c,wxo18c,wxne20c,wxne22c,wxmg24c,wxmg25c,wxmg26c,wxc14c,wxf18c, &
    wxf19c,wxne21c,wxna23c,wxal26gc,wxal27c,wxsi28c,wxneutc,wxprotc,wxbidc,wxbid1c,val,ecart,suma,sumb,sum,sumc,sumcen, &
    differ,t9
  integer,parameter:: idimnetc=22 !Temporay switch back to 15
  real(kindreal), dimension(ldi):: at,bt,ct,br,dm,dqv
  real(kindreal), dimension(ldi):: wwx,wwy3,wwy,wwxc12,wwxc13,wwxn14,wwxn15,wwxo16,wwxo17,wwxo18,wwxne20,wwxne22,wwxmg24,wwxmg25, &
                                   wwxmg26,wwx19,wwx21,wwx23,wwxag,wwx27,wwx28,wwxne,wwxpr,wwx14,wwx18,wwxbi,wwxb1,wwtridagx
  real(kindreal), dimension(mbelx):: wabelxc,sumabelx
  real(kindreal), dimension(idimnetc):: vxab2
  logical :: normalise
!-----------------------------------------------------------------------
  if (tdiff == 0.d0) then
    tdiff=dzeit/2.d0
  endif

  wx(1:m)=vvx(m:1:-1)
  wy3(1:m)=vvy3(m:1:-1)
  wy(1:m)=vvy(m:1:-1)
  wxc12(1:m)=vvxc12(m:1:-1)
  wxc13(1:m)=vvxc13(m:1:-1)
  wxn14(1:m)=vvxn14(m:1:-1)
  wxn15(1:m)=vvxn15(m:1:-1)
  wxo16(1:m)=vvxo16(m:1:-1)
  wxo17(1:m)=vvxo17(m:1:-1)
  wxo18(1:m)=vvxo18(m:1:-1)
  wxne20(1:m)=vvxne20(m:1:-1)
  wxne22(1:m)=vvxne22(m:1:-1)
  wxmg24(1:m)=vvxmg24(m:1:-1)
  wxmg25(1:m)=vvxmg25(m:1:-1)
  wxmg26(1:m)=vvxmg26(m:1:-1)
  if (ialflu == 1) then
    wxc14(1:m)=vvxc14(m:1:-1)
    wxf18(1:m)=vvxf18(m:1:-1)
    wxf19(1:m)=vvxf19(m:1:-1)
    wxne21(1:m)=vvxne21(m:1:-1)
    wxna23(1:m)=vvxna23(m:1:-1)
    wxal26g(1:m)=vvxal26g(m:1:-1)
    wxal27(1:m)=vvxal27(m:1:-1)
    wxsi28(1:m)=vvxsi28(m:1:-1)
    wxneut(1:m)=vvxneut(m:1:-1)
    wxprot(1:m)=vvxprot(m:1:-1)
    wxbid(1:m)=vvxbid(m:1:-1)
    wxbid1(1:m)=vvxbid1(m:1:-1)
  endif
  wabelx(1:nbelx,1:m)=vvabelx(1:nbelx,m:1:-1)

  if (jzint /= 0) then
    do jk=1,2*jzint
     izc(jk)=m+1-izc(jk)
    enddo
  endif

  do i=1,m-1
   dra(i)=exp(rb(m-i))-exp(rb(m-i+1))
   dqv(i)=M_r(m-i)-M_r(m-i+1)
  enddo

! calcul du Dmoyen entre deux coquilles
  do i=2,m
   if (D_chim(i) > 0.d0 .and. D_chim(i-1) > 0.d0) then
     D_moychim(i)=D_chim(i)*D_chim(i-1)/(0.5d0*D_chim(i-1)+0.5d0*D_chim(i))
   else
     D_moychim(i)=0.d0
   endif
  enddo

!*******************************************************************
!      BOUCLE DE RESOLUTION A CHAQUE PAS DE TEMPS
!*******************************************************************
  jbid=0
  call couran
  if (verbose) then
    write(*,*) ' tdiff=',tdiff
  endif
  if (tdiff <= 0.0d0) then
    write(10,'(1x,a,d12.5)') ' tdiff= ',tdiff
  endif
! On met dans l'ordre habituel les variables

  do i=1,m-1
   r32=(exp(rb(i))+exp(rb(i+1)))/2.d0
   rho32=(exp(rho(i))+exp(rho(i+1)))/2.d0
   dr32=exp(rb(i+1))-exp(rb(i))
   if (dr32 /= 0.d0) then
     br(i)=4.0d0*pi*r32*r32*rho32*D_moychim(i+1)*tdiff/dr32
   else
     br(i) = 0.d0
   endif
  enddo
! calcul de la masse des coquilles
  dm(1)=(M_r(1)-M_r(2))/2.d0
  do i=2,m-1
   dm(i)=(M_r(i-1)-M_r(i+1))/2.d0
  enddo
  dm(m)=(M_r(m-1)-M_r(m))/2.d0
! i=1
  at(1)=0.0d0
  bt(1)=1.d0-br(1)/dm(1)
  ct(1)=br(1)/dm(1)
! i=2,m-1
  do i=2,m-1
   at(i)=br(i-1)/dm(i)
   bt(i)=1.d0-br(i-1)/dm(i)-br(i)/dm(i)
   ct(i)=br(i)/dm(i)
  enddo
! i=m
  at(m)=br(m-1)/dm(m)
  bt(m)=1.d0-br(m-1)/dm(m)
  ct(m)=0.0d0

  do while (jbid  <  nwpas)

   wwx(1:m)=wx(m:1:-1)
   wwy3(1:m)=wy3(m:1:-1)
   wwy(1:m)=wy(m:1:-1)
   wwxc12(1:m)=wxc12(m:1:-1)
   wwxc13(1:m)=wxc13(m:1:-1)
   wwxn14(1:m)=wxn14(m:1:-1)
   wwxn15(1:m)=wxn15(m:1:-1)
   wwxo16(1:m)=wxo16(m:1:-1)
   wwxo17(1:m)=wxo17(m:1:-1)
   wwxo18(1:m)=wxo18(m:1:-1)
   wwxne20(1:m)=wxne20(m:1:-1)
   wwxne22(1:m)=wxne22(m:1:-1)
   wwxmg24(1:m)=wxmg24(m:1:-1)
   wwxmg25(1:m)=wxmg25(m:1:-1)
   wwxmg26(1:m)=wxmg26(m:1:-1)
   if (ialflu == 1) then
     wwx14(1:m)=wxc14(m:1:-1)
     wwx18(1:m)=wxf18(m:1:-1)
     wwx19(1:m)=wxf19(m:1:-1)
     wwx21(1:m)=wxne21(m:1:-1)
     wwx23(1:m)=wxna23(m:1:-1)
     wwxag(1:m)=wxal26g(m:1:-1)
     wwx27(1:m)=wxal27(m:1:-1)
     wwx28(1:m)=wxsi28(m:1:-1)
     wwxne(1:m)=wxneut(m:1:-1)
     wwxpr(1:m)=wxprot(m:1:-1)
     wwxbi(1:m)=wxbid(m:1:-1)
     wwxb1(1:m)=wxbid1(m:1:-1)
   endif
! Nouvelle methode pour la resolution de l'equation de diffusion
! On ne renverse pas la numerotation des coquilles
! calcul des elements de matrice
   call tridiago(at,bt,ct,wwx,m)
   call tridiago(at,bt,ct,wwy3,m)
   call tridiago(at,bt,ct,wwy,m)
   call tridiago(at,bt,ct,wwxc12,m)
   call tridiago(at,bt,ct,wwxc13,m)
   call tridiago(at,bt,ct,wwxn14,m)
   call tridiago(at,bt,ct,wwxn15,m)
   call tridiago(at,bt,ct,wwxo16,m)
   call tridiago(at,bt,ct,wwxo17,m)
   call tridiago(at,bt,ct,wwxo18,m)
   call tridiago(at,bt,ct,wwxne20,m)
   call tridiago(at,bt,ct,wwxne22,m)
   call tridiago(at,bt,ct,wwxmg24,m)
   call tridiago(at,bt,ct,wwxmg25,m)
   call tridiago(at,bt,ct,wwxmg26,m)
   if (ialflu == 1) then !NO F18 ADN C14 ??  !WARNING CAUSED BUG IN EALRY PHASES ASK WHY
     if ( ( phase .gt. 3 )  .and.  ( t(m) .gt. log(3e8) ) ) then
        call tridiago(at,bt,ct,wwx14,m)
        call tridiago(at,bt,ct,wwx18,m)
     endif
     call tridiago(at,bt,ct,wwx19,m)
     call tridiago(at,bt,ct,wwx21,m)
     call tridiago(at,bt,ct,wwx23,m)
     call tridiago(at,bt,ct,wwxag,m)
     call tridiago(at,bt,ct,wwx27,m)
     call tridiago(at,bt,ct,wwx28,m)
     call tridiago(at,bt,ct,wwxbi,m)
     call tridiago(at,bt,ct,wwxb1,m)
   endif

   do ii=1,nbelx
    wwtridagx(1:m)=wabelx(ii,m:1:-1)
    call tridiago(at,bt,ct,wwtridagx,m)
    wabelx(ii,1:m)=wwtridagx(m:1:-1)
   enddo

   wx(1:m)=wwx(m:1:-1)
   wy3(1:m)=wwy3(m:1:-1)
   wy(1:m)=wwy(m:1:-1)
   wxc12(1:m)=wwxc12(m:1:-1)
   wxc13(1:m)=wwxc13(m:1:-1)
   wxn14(1:m)=wwxn14(m:1:-1)
   wxn15(1:m)=wwxn15(m:1:-1)
   wxo16(1:m)=wwxo16(m:1:-1)
   wxo17(1:m)=wwxo17(m:1:-1)
   wxo18(1:m)=wwxo18(m:1:-1)
   wxne20(1:m)=wwxne20(m:1:-1)
   wxne22(1:m)=wwxne22(m:1:-1)
   wxmg24(1:m)=wwxmg24(m:1:-1)
   wxmg25(1:m)=wwxmg25(m:1:-1)
   wxmg26(1:m)=wwxmg26(m:1:-1)
   if (ialflu == 1) then
     wxc14(1:m)=wwx14(m:1:-1)
     wxf18(1:m)=wwx18(m:1:-1)
     wxf19(1:m)=wwx19(m:1:-1)
     wxne21(1:m)=wwx21(m:1:-1)
     wxna23(1:m)=wwx23(m:1:-1)
     wxal26g(1:m)=wwxag(m:1:-1)
     wxal27(1:m)=wwx27(m:1:-1)
     wxsi28(1:m)=wwx28(m:1:-1)
     wxneut(1:m)=wwxne(m:1:-1)
     wxprot(1:m)=wwxpr(m:1:-1)
     wxbid(1:m)=wwxbi(m:1:-1)
     wxbid1(1:m)=wwxb1(m:1:-1)
   endif

   if (jzint /= 0) then
     do jk=1,jzint
      i1=izc(2*jk-1)
      i2=izc(2*jk)
      sumx=0.d0
      sumy=0.d0
      sumy3=0.d0
      sumn14=0.d0
      sumn15=0.d0
      sumc12=0.d0
      sumc13=0.d0
      sumo16=0.d0
      sumo17=0.d0
      sumo18=0.d0
      sumne20=0.d0
      sumne22=0.d0
      summg24=0.d0
      summg25=0.d0
      summg26=0.d0
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
      sumabelx(1:nbelx)=0.d0

      do i=i1,i2-1
       cc=dqv(i)/2.d0
       sumx=cc*(wx(i+1)+wx(i))+sumx
       sumy=cc*(wy(i+1)+wy(i))+sumy
       sumy3=cc*(wy3(i+1)+wy3(i))+sumy3
       sumc12=cc*(wxc12(i+1)+wxc12(i))+sumc12
       sumc13=cc*(wxc13(i+1)+wxc13(i))+sumc13
       sumn14=cc*(wxn14(i+1)+wxn14(i))+sumn14
       sumn15=cc*(wxn15(i+1)+wxn15(i))+sumn15
       sumo16=cc*(wxo16(i+1)+wxo16(i))+sumo16
       sumo17=cc*(wxo17(i+1)+wxo17(i))+sumo17
       sumo18=cc*(wxo18(i+1)+wxo18(i))+sumo18
       sumne20=cc*(wxne20(i+1)+wxne20(i))+sumne20
       sumne22=cc*(wxne22(i+1)+wxne22(i))+sumne22
       summg24=cc*(wxmg24(i+1)+wxmg24(i))+summg24
       summg25=cc*(wxmg25(i+1)+wxmg25(i))+summg25
       summg26=cc*(wxmg26(i+1)+wxmg26(i))+summg26
       if (ialflu == 1) then
       if ( ( phase .gt. 3 )  .and.  ( t(i) .gt. log(3e8) ) ) then
         sumc14=cc*(wxc14(i+1)+wxc14(i))+sumc14
         sumf18=cc*(wxf18(i+1)+wxf18(i))+sumf18
       endif
         sumf19=cc*(wxf19(i+1)+wxf19(i))+sumf19
         sumne21=cc*(wxne21(i+1)+wxne21(i))+sumne21
         sumna23=cc*(wxna23(i+1)+wxna23(i))+sumna23
         sumal26g=cc*(wxal26g(i+1)+wxal26g(i))+sumal26g
         sumal27=cc*(wxal27(i+1)+wxal27(i))+sumal27
         sumsi28=cc*(wxsi28(i+1)+wxsi28(i))+sumsi28
         sumneut=cc*(wxneut(i+1)+wxneut(i))+sumneut
         sumprot=cc*(wxprot(i+1)+wxprot(i))+sumprot
         sumbid=cc*(wxbid(i+1)+wxbid(i))+sumbid
         sumbid1=cc*(wxbid1(i+1)+wxbid1(i))+sumbid1
       endif
       do ii=1,nbelx
        sumabelx(ii)=cc*(wabelx(ii,i+1)+wabelx(ii,i))+sumabelx(ii)
       enddo
      enddo
      aa=0.d0

      if (i2 == m) then
        aa=gmsu-M_r(1)
      endif
      bb=aa+M_r(m-i2+1)-M_r(m-i1+1)
      whc=(aa*wx(m)+sumx)/bb
      wyc=(aa*wy(m)+sumy)/bb
      wy3c=(aa*wy3(m)+sumy3)/bb
      wxc12c=(aa*wxc12(m)+sumc12)/bb
      wxc13c=(aa*wxc13(m)+sumc13)/bb
      wxn14c=(aa*wxn14(m)+sumn14)/bb
      wxn15c=(aa*wxn15(m)+sumn15)/bb
      wxo16c=(aa*wxo16(m)+sumo16)/bb
      wxo17c=(aa*wxo17(m)+sumo17)/bb
      wxo18c=(aa*wxo18(m)+sumo18)/bb
      wxne20c=(aa*wxne20(m)+sumne20)/bb
      wxne22c=(aa*wxne22(m)+sumne22)/bb
      wxmg24c=(aa*wxmg24(m)+summg24)/bb
      wxmg25c=(aa*wxmg25(m)+summg25)/bb
      wxmg26c=(aa*wxmg26(m)+summg26)/bb
      if (ialflu == 1) then
      if ( ( phase .gt. 3 )  .and.  ( t(m) .gt. log(3e8) ) ) then
        wxc14c=(aa*wxc14(m)+sumc14)/bb
        wxf18c=(aa*wxf18(m)+sumf18)/bb
      endif
        wxf19c=(aa*wxf19(m)+sumf19)/bb
        wxne21c=(aa*wxne21(m)+sumne21)/bb
        wxna23c=(aa*wxna23(m)+sumna23)/bb
        wxal26gc=(aa*wxal26g(m)+sumal26g)/bb
        wxal27c=(aa*wxal27(m)+sumal27)/bb
        wxsi28c=(aa*wxsi28(m)+sumsi28)/bb
        wxneutc=(aa*wxneut(m)+sumneut)/bb
        wxprotc=(aa*wxprot(m)+sumprot)/bb
        wxbidc=(aa*wxbid(m)+sumbid)/bb
        wxbid1c=(aa*wxbid1(m)+sumbid1)/bb
      endif
      do ii=1,nbelx
       wabelxc(ii)=(aa*wabelx(ii,m)+sumabelx(ii))/bb
      enddo

      wx(i1:i2)=whc
      wy3(i1:i2)=wy3c
      wy(i1:i2)=wyc
      wxn14(i1:i2)=wxn14c
      wxn15(i1:i2)=wxn15c
      wxc12(i1:i2)=wxc12c
      wxc13(i1:i2)=wxc13c
      wxo16(i1:i2)=wxo16c
      wxo17(i1:i2)=wxo17c
      wxo18(i1:i2)=wxo18c
      wxne20(i1:i2)=wxne20c
      wxne22(i1:i2)=wxne22c
      wxmg24(i1:i2)=wxmg24c
      wxmg25(i1:i2)=wxmg25c
      wxmg26(i1:i2)=wxmg26c
      if (ialflu == 1) then
        if ( ( ( phase .gt. 3 )  .and.  ( t(m) .gt. log(3e8) ) ) ) then
          wxc14(i1:i2)=wxc14c
          wxf18(i1:i2)=wxf18c !again c14 and f18 left out. Mabye decay =? 
        endif
        wxf19(i1:i2)=wxf19c
        wxne21(i1:i2)=wxne21c
        wxna23(i1:i2)=wxna23c
        wxal26g(i1:i2)=wxal26gc
        wxal27(i1:i2)=wxal27c
        wxsi28(i1:i2)=wxsi28c
        wxbid(i1:i2)=wxbidc
        wxbid1(i1:i2)=wxbid1c
      endif
      do ii=1,nbelx
       wabelx(ii,i1:i2)= wabelxc(ii)
      enddo
     enddo
   endif ! jzint /= 0
!        ****************************************************
!         REINITIALISATION ET STOCKAGE DES DONNEES
!    ********************************************************
   jbid=jbid+1
  enddo

  val=1.d0
  ecart=0.d0
  nnn=0
  suma=0.d0
  sumb=0.d0
  sum=0.d0
  do n=1,m
   nm=m-n+1
    if  (  t(nm) < log(3e8)) then !Dont mix protons in advance phase
        x(nm)=x(nm)-vvx(nm)+wx(n)
    endif

   if (x(nm) <= 1.d-09 .and. t(nm) < log(3e8)) then
     x(nm)=0.d0
   endif
   if (epsc(nm) == 0.0d0 ) then !Never mix alphas in adv phase.
     y(nm)=y(nm)-vvy(nm)+wy(n)
   endif

   y3(nm)=y3(nm)-vvy3(nm)+wy3(n)
   xc12(nm)=xc12(nm)-vvxc12(nm)+wxc12(n)
   xc13(nm)=xc13(nm)-vvxc13(nm)+wxc13(n)
   xn14(nm)=xn14(nm)-vvxn14(nm)+wxn14(n)
   xn15(nm)=xn15(nm)-vvxn15(nm)+wxn15(n)
   xo16(nm)=xo16(nm)-vvxo16(nm)+wxo16(n)
   xo17(nm)=xo17(nm)-vvxo17(nm)+wxo17(n)
   xo18(nm)=xo18(nm)-vvxo18(nm)+wxo18(n)
   xne20(nm)=xne20(nm)-vvxne20(nm)+wxne20(n)
   xne22(nm)=xne22(nm)-vvxne22(nm)+wxne22(n)
   xmg24(nm)=xmg24(nm)-vvxmg24(nm)+wxmg24(n)
   xmg25(nm)=xmg25(nm)-vvxmg25(nm)+wxmg25(n)
   xmg26(nm)=xmg26(nm)-vvxmg26(nm)+wxmg26(n)
   if (ialflu == 1) then
    if ( ( phase .gt. 3 )  .and.  ( t(nm) .gt. log(3e8) ) ) then
     xc14(nm)=xc14(nm)-vvxc14(nm)+wxc14(n)
     xf18(nm)=xf18(nm)-vvxf18(nm)+wxf18(n)
    endif
     xf19(nm)=xf19(nm)-vvxf19(nm)+wxf19(n)
     xne21(nm)=xne21(nm)-vvxne21(nm)+wxne21(n)
     xna23(nm)=xna23(nm)-vvxna23(nm)+wxna23(n)
     if (wxal26g(n) < 1.d-20) then
       xal26(nm)=0.d0
     else
       xal26(nm)=xal26(nm)-vvxal26g(nm)+wxal26g(n)
     endif
     xal27(nm)=xal27(nm)-vvxal27(nm)+wxal27(n)
     xsi28(nm)=xsi28(nm)-vvxsi28(nm)+wxsi28(n)
     xneut(nm)=xneut(nm)-vvxneut(nm)+wxneut(n)
     if (xneut(nm) <= 1.d-30) then
       xneut(nm)=0.d0
     endif
     xprot(nm)=xprot(nm)-vvxprot(nm)+wxprot(n)
     if (xprot(nm) <= 1.d-30) then
       xprot(nm)=0.d0
     endif
     xbid(nm)=xbid(nm)-vvxbid(nm)+wxbid(n)
     xbid1(nm)=xbid1(nm)-vvxbid1(nm)+wxbid1(n)
   endif
   do ii=1,nbelx
  !  !Stop neutron diff
    if (nbzel(ii) /= 0) then
      abelx(ii,nm)=abelx(ii,nm)-vvabelx(ii,nm)+wabelx(ii,n)
    endif
    if (abelx(ii,nm) < 0.d0.or.abelx(ii,nm) > 1.d0) then
      write(10,*) 'abelx',nm,n, ii,abelx(ii,nm),vvabelx(ii,nm),wabelx(ii,n), nbzel(ii),nbael(ii)
    endif
   enddo

   sumb=0.d0
   sumc=0.d0
   suma=x(nm)+y3(nm)+y(nm)+xc12(nm)+xc13(nm)+xn14(nm)+xn15(nm)+xo16(nm)+xo17(nm)+xo18(nm)+xne20(nm)+xne22(nm)+ &
        xmg24(nm)+xmg25(nm)+xmg26(nm)
   if (ialflu == 1) then
     sumb=xc14(nm)+xf18(nm)+xf19(nm)+xne21(nm)+xna23(nm)+xal26(nm)+xal27(nm)+xsi28(nm)+xneut(nm)+xprot(nm)+ &
           xbid(nm)+xbid1(nm)
   endif
   do ii=1,nbelx
    sumc=sumc+abelx(ii,nm)
   enddo
   if (ialflu /= 1) then
     sumc=sumc+zabelx
   endif

   sum=suma+sumb+sumc

   if (nm == m) then
     sumcen=sum
   endif
   differ=abs(1.d0-sum)
   if (differ > ecart) then
     ecart=differ
     nnn=nm
     val=sum
   endif
    normalise = .False.
    if (normalise) then
      x(nm) = x(nm) / sum
      y3(nm) = y3(nm) / sum
      y(nm) = y(nm) / sum
      xc12(nm) = xc12(nm) / sum
      xc13(nm) = xc13(nm) / sum
      xn14(nm) = xn14(nm) / sum
      xn15(nm) = xn15(nm) / sum
      xo16(nm) = xo16(nm) / sum
      xo17(nm) = xo17(nm) / sum
      xo18(nm) = xo18(nm) / sum
      xne20(nm) = xne20(nm) / sum
      xne22(nm) = xne22(nm) / sum
      xmg24(nm) = xmg24(nm) / sum
      xmg25(nm) = xmg25(nm) / sum
      xmg26(nm) = xmg26(nm) / sum
      if (ialflu == 1) then
        xc14(nm) = xc14(nm) / sum
        xf18(nm) = xf18(nm) / sum
        xf19(nm) = xf19(nm) / sum
        xne21(nm) = xne21(nm) / sum
        xna23(nm) = xna23(nm) / sum
        xal26(nm) = xal26(nm) / sum
        xal27(nm) = xal27(nm) / sum
        xsi28(nm) = xsi28(nm) / sum
        xneut(nm) = xneut(nm) / sum
        xprot(nm) = xprot(nm) / sum
        xbid(nm) = xbid(nm) / sum
        xbid1(nm) = xbid1(nm) / sum
      endif
      do ii=1,nbelx
        abelx(ii,nm) = abelx(ii,nm) / sum
      enddo
    endif

  enddo
  write(3,'(a,0pf13.9,3x,a,f13.9,3x,a,f13.9,3x,a,i4)') 'Tracking abundances suma',suma,'sumb',sumb,'sumc',sumc,'nbelx',nbelx

! 26 septembre 2000, Georges Meynet
! Le probleme: en raison de la diffusion, de l'hydrogene
! peut etre ramene dans le coeur ou il avait disparu. Le pgm
! repasse alors dans le reseau de combustion de l'H et oscille
! entre H et He. Il se peut meme que des coquilles ou
! on a combustion de l'He s'intercalent entre des coquilles
! de combustion de l'hydrogene. Il faut a mon avis adopter
! dans ces cas la une regle claire qui est la suivante:
! dans une coquille ou on est passe une fois dans le
! reseau de combustion de l'helium, l'hydrogene est
! immediatement transforme en helium. C'est ce que font
! les lignes suivantes.
  if (ipop3 == 0 ) then !adam flag
    do n=1,m
     if (epsy(n) /= 0.d0 .and. x(n) /= 0.d0) then
       y(n)=y(n)+x(n)
       x(n)=0.d0
     endif
    enddo
  endif
  if (phase >= 3 .and. idifcon == 1) then !adam flag
    do n=1,m
     nm=m-n+1
     if (epsc(nm) /= 0.0d0) then
       t9=exp(t(nm)-log(1.d9))
       vxab2(1) =   x(nm)
       vxab2(2) =   y3(nm)
       vxab2(3) =   y(nm)
       vxab2(4) =   xc12(nm)
       vxab2(5) =   xc13(nm)
       vxab2(6) =   xn14(nm)
       vxab2(7) =   xn15(nm)
       vxab2(8) =   xo16(nm)
       vxab2(9) =   xo17(nm)
       vxab2(10)=   xo18(nm)
       vxab2(11)=   xne20(nm)
       vxab2(12)=   xne22(nm)
       vxab2(13)=   xmg24(nm)
       vxab2(14)=   xmg25(nm)
       vxab2(15)=   xmg26(nm)
       if (ialflu == 1 ) then
          vxab2(16) = xc14(nm)
          vxab2(17) = xf18(nm)
          vxab2(18) = xf19(nm)
          vxab2(19) = xna23(nm)
          vxab2(20) = xne21(nm)
          vxab2(21) = xal26(nm)
          vxab2(22) = xal27(nm)  
       else
          vxab2(16:22) = 0.d0
       endif
       if (nm >= m-2) then
         write(3,'(a,i4,77(1x,e17.10))') 'BEF. NETBURN',nm,(vxab2(i),i=1,idimnetc),(abelx(ii,m),ii=1,nbelx)
         write(3,'(i4,3(1p,e12.5))') nm,t9,dzeit,fnucdif
       endif
       call netburning(nm,t9,dzeit,vxab2,2)
       if (nm >= m-2) then
         write(3,'(a,i4,77(1x,e17.10))') 'AFT. NETBURN',nm,(vxab2(i),i=1,idimnetc),(abelx(ii,m),ii=1,nbelx)
       endif

       x(nm)   = vxab2(1)
       y3(nm)  = vxab2(2)
       y(nm)   = vxab2(3)
       xc12(nm)  = vxab2(4)
       xc13(nm)= vxab2(5)
       xn14(nm)= vxab2(6)
       xn15(nm)= vxab2(7)
       xo16(nm)  = vxab2(8)
       xo17(nm)= vxab2(9)
       xo18(nm)= vxab2(10)
       xne20(nm) = vxab2(11)
       xne22(nm) = vxab2(12)
       xmg24(nm) = vxab2(13)
       xmg25(nm) = vxab2(14)
       xmg26(nm) = vxab2(15)
       if ( ialflu == 1 ) then
          xc14(nm) = vxab2(16)
          xf18(nm) = vxab2(17)
          xf19(nm) = vxab2(18)
          xna23(nm) = vxab2(19)
          xne21(nm) = vxab2(20)
          xal26(nm) = vxab2(21)
          xal27(nm) = vxab2(22)
       endif
     endif
    enddo
  endif

! mise a zero des abondances negatives
! effet des reactions nucleaires et de la diffusion
  do n=1,m
! Si cela arrive que fait-on ?
! Cas ou X negatif et Y plus grand que un
! Que vaut l'ecart ?
   if (x(n) < 0.0d0) then
     write(10,'(a,i4,a,f10.6)') 'ATTENTION X NEG. COUCHE ',n,' X=',x(n)
   endif
   if (x(n) < 1.0d-9 .and. (t(n) < log(3e8))) then
     x(n)=0.0d0
   endif
   if (y(n) > 1.0d0) then
     write(10,'(a,i4,a,f10.6)') 'ATTENTION Y SUP A 1 COUCHE ',n,' Y=',y(n)
   endif
! attention phases avancees: on a besoin de l'He !
   if (y(n) < 1.0d-25 .and. phase <= 2) then
     y(n)=0.0d0
   endif
   if (y3(n) < 1.0d-25) then
     y3(n)=0.0d0
   endif
   if (ipop3 == 0 .or. phase >= 2) then
     if (xc12(n) < 1.0d-25) then
       xc12(n)=0.0d0
     endif
     if (xc13(n) < 1.0d-25) then
       xc13(n)=0.0d0
     endif
     if (xn14(n) < 1.0d-25) then
       xn14(n)=0.0d0
     endif
     if (xn15(n) < 1.0d-25) then
       xn15(n)=0.0d0
     endif
     if (xo16(n) < 1.0d-25) then
       xo16(n)=0.0d0
     endif
     if (xo17(n) < 1.0d-25) then
       xo17(n)=0.0d0
     endif
     if (xo18(n) < 1.0d-25) then
       xo18(n)=0.0d0
     endif
     if (xne20(n) < 1.0d-25) then
       xne20(n)=0.0d0
     endif
     if (xne22(n) < 1.0d-25) then
       xne22(n)=0.0d0
     endif
     if (xmg24(n) < 1.0d-25) then
       xmg24(n)=0.0d0
     endif
     if (xmg25(n) < 1.0d-25) then
       xmg25(n)=0.0d0
     endif
     if (xmg26(n) < 1.0d-25) then
       xmg26(n)=0.0d0
     endif
   else
     if (xc12(n) < 1.0d-75) then
       xc12(n)=0.0d0
     endif
     if (xc13(n) < 1.0d-75) then
       xc13(n)=0.0d0
     endif
     if (xn14(n) < 1.0d-75) then
       xn14(n)=0.0d0
     endif
     if (xn15(n) < 1.0d-75) then
       xn15(n)=0.0d0
     endif
     if (xo16(n) < 1.0d-75) then
       xo16(n)=0.0d0
     endif
     if (xo17(n) < 1.0d-75) then
       xo17(n)=0.0d0
     endif
     if (xo18(n) < 1.0d-75) then
       xo18(n)=0.0d0
     endif
     if (xne20(n) < 1.0d-75) then
       xne20(n)=0.0d0
     endif
     if (xne22(n) < 1.0d-75) then
       xne22(n)=0.0d0
     endif
     if (xmg24(n) < 1.0d-75) then
       xmg24(n)=0.0d0
     endif
     if (xmg25(n) < 1.0d-75) then
       xmg25(n)=0.0d0
     endif
     if (xmg26(n) < 1.0d-75) then
       xmg26(n)=0.0d0
     endif
   endif
   do ii=1,nbelx
    if (abelx(ii,n) < 0.0d0) then
      write(10,*) 'ATTENTION COUCHE: ',n,'el ',ii,': ab.=',abelx(ii,n)
    endif
    if (abelx(ii,n) < 1.d-50) then
      abelx(ii,n)= 0.d0
    endif
   enddo
  enddo

  if (jdiff /= 0) then
    write(3,'(a,0pf13.9,3x,a,i5,3x,a,f13.9)') 'WORST SUM OF ABUNDANCES',val,'COUCHE',nnn,'CENTRAL SUM',sumcen
  endif
  write(3,'(1x,a,1x,1pe10.3,a,1x,e10.3)') 'x(m):',x(m),' y(m):',y(m)

  return

end subroutine diffbr

!=======================================================================
subroutine diffom
!-----------------------------------------------------------------------
  use const,only: Lsol,cst_sigma
  use inputparam,only: iadvec,imagn,xcn,phase
  use caramodele,only: nwmd,glm,gls,teff
  use strucmod,only: m,q,rb,rho
  use rotmod,only: vvomeg,omegd,vsuminenv,xldoex,Flux_remaining
  use diffadvmod,only: tdiff
  use SmallFunc,only: tridiago
  use timestep,only: dzeit

  implicit none

  integer:: i,n,jbid
  real(kindreal):: btota,xmocin,dbrmr,dbrmr1,dbrmrm,r32,rho32,dr32,Rstar,&
    M_env,rhomoy_env,rmoy_env,dbrmr2,dbrmr12,dbrmrm2,btota2,inert1
  real(kindreal), dimension(0:ldi):: at,bt,ct,br,dm
  real(kindreal), dimension(0:ldi):: xmr,r2
  logical:: orderedR = .false.
!-----------------------------------------------------------------------
  if (tdiff == 0.0d0) then
    tdiff=dzeit/2.d0
  endif

  xmr(0:m) = 0.d0
  r2(0:m) = 0.d0

  xmr(1:m)=exp(glm)*(1.d0-exp(q(1:m)))
  xmr(0) = exp(glm)
  Rstar = sqrt(gls*Lsol/(4.d0*pi*cst_sigma))/teff**2.d0
  M_env = xmr(0)-xmr(1)
  rhomoy_env = 3.d0/(4.d0*pi)*M_env/(Rstar**3.d0-exp(3.d0*rb(1)))
  if (isnan(vsuminenv)) then
    write(*,*) 'vsuminenv=NaN'
    stop
  endif
  rmoy_env = (3.d0/2.d0)*vsuminenv/M_env
  orderedR = sqrt(rmoy_env) > exp(rb(1))
  if (.not. orderedR) then
      write(*,*) 'BEWARE - Radius inversion in diffom, rmoy_env set to r(1)+1%'
      write(*,*) '         Rstar,rmoy_env,r(1),vsuminenv: ',Rstar,sqrt(rmoy_env),exp(rb(1)),vsuminenv
      rmoy_env = (1.01d0*exp(rb(1)))**2.d0
  endif

  omega_extended(1:m)=vvomeg(1:m)
  omega_extended(0) = omega_extended(1)
  r2(1:m)=exp(rb(1:m))*exp(rb(1:m))
  r2(0) = rmoy_env

! calcul du moment angulaire total
  btota=0.d0
  do n=2,m-1
   xmocin=r2(n)*omega_extended(n)*2.d0/3.d0
   dbrmr=xmocin*(xmr(n-1)-xmr(n+1))/2.d0
   btota=btota+dbrmr
  enddo
!  dbrmr1=vsuminenv*omega_extended(0)
  xmocin=2.d0/3.d0*r2(1)*omega_extended(1)
  dbrmr1=xmocin*(xmr(1)-xmr(2))/2.d0
  dbrmr1=dbrmr1 + vsuminenv*omega_extended(0)
  dbrmrm=2.d0/5.d0*xmr(m-1)/2.d0*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*omega_extended(m)
  btota=btota+dbrmr1+dbrmrm

! calcul du Dmoyen entre deux coquilles
  do i=2,m
   if (D_Omega(i) > 0.d0 .and. D_Omega(i-1) > 0.d0) then
     D_moyOm(i)=D_Omega(i)*D_Omega(i-1)/(0.5d0*D_Omega(i-1)+0.5d0*D_Omega(i))
   else
     D_moyOm(i)=0.d0
   endif
  enddo
! On fixe le coefficient de diffusion entre la premiere couche et l'enveloppe a une valeur
! tres elevee pour assurer un profil plat de Omega.
  D_moyOm(1) = 1.0d20

!*******************************************************************
!      BOUCLE DE RESOLUTION A CHAQUE PAS DE TEMPS
!*******************************************************************
  jbid=0
  call courom
  if (tdiff == 0.0d0) then
    write(10,*)' DIFFOM: tdiff=0.'
  endif

  r32=(sqrt(r2(0))+exp(rb(1)))/2.d0
  rho32=(rhomoy_env+exp(rho(1)))/2.d0
  dr32=exp(rb(1))-sqrt(r2(0))
  br(0) = 4.d0*pi*r32**4.d0*rho32*D_moyOm(1)*tdiff/dr32
  do i=1,m-1
   r32=(exp(rb(i))+exp(rb(i+1)))/2.d0
   rho32=(exp(rho(i))+exp(rho(i+1)))/2.d0
   dr32=exp(rb(i+1))-exp(rb(i))
   if (dr32 /= 0.d0) then
     br(i)=4.d0*pi*r32*r32*rho32*D_moyOm(i+1)*tdiff/dr32
   else
     br(i) = 0.d0
   endif
   br(i)=br(i)*r32*r32
  enddo
! calcul du moment d'inertie des coquilles
  dm(0)=(xmr(0)-xmr(1))*r2(0)
  dm(1) =(xmr(1)-xmr(2))/2.d0*r2(1)
  do i=2,m-1
   dm(i)=(xmr(i-1)-xmr(i+1))/2.d0*r2(i)
  enddo
! Au centre, le facteur numerique du moment d'inertie vaut 2/5
  dm(m)=(xmr(m-1)-xmr(m))/5.d0*r2(m-1)/4.d0
! i=1
  at(0)=0.0d0
  bt(0)=1.d0-br(0)/dm(0)
  ct(0)=br(0)/dm(0)
! i=2,m-1
  do i=1,m-1
   at(i)=br(i-1)/dm(i)
   bt(i)=1.d0-br(i-1)/dm(i)-br(i)/dm(i)
   ct(i)=br(i)/dm(i)
  enddo
! i=m
  at(m)=br(m-1)/dm(m)
  bt(m)=1.d0-br(m-1)/dm(m)
  ct(m)=0.0d0
  if (omega_extended(0) < 0.d0) then
    omega_extended(0) = 1.d-20
    omega_extended(1) = 1.d-20
  endif
!  write(*,*) 'Entree dans tridiago: omega,xldoex,tdiff,dm:',omega_extended(0),xldoex,tdiff,dm(0)
  do while (jbid < nwpas)

! Nouvelle methode pour la resolution de l'equation de diffusion
! On ne renverse pas la numerotation des coquilles
! calcul des elements de matrice
   omega_extended(0) = omega_extended(0) + 3.d0*(xldoex+Flux_remaining)*tdiff/(2.d0*dm(0))
   if (omega_extended(0) < 0.d0) then
     if (phase <= 5 .and. dzeit < 1.d0) then
       stop 'omega neg before tridiago'
     else
       omega_extended(0) = 0.d0
     endif
   endif
   call tridiago(at(0:m),bt(0:m),ct(0:m),omega_extended(0:m),m+1)
   jbid=jbid+1

  enddo

  omegd(2:m)=omega_extended(2:m)
  inert1 = r2(1)*2.d0/3.d0*(xmr(0)-xmr(2))/2.d0
! On redefinit Omega(1) comme un omega moyen tel que (L0 + L1 = Omega(1)*(I0 + I1))
  omegd(1) = (inert1*omega_extended(1) + vsuminenv*omega_extended(0)) &
             /(inert1+vsuminenv)

! calcul du moment angulaire total
  btota2=0.d0
!  do n=1,m-1
  do n=2,m-1
   xmocin=r2(n)*omegd(n)*2.d0/3.d0
   dbrmr2=xmocin*(xmr(n-1)-xmr(n+1))/2.d0
   btota2=btota2+dbrmr2
  enddo
  xmocin=2.d0/3.d0*r2(1)*omega_extended(1)
  dbrmr12=xmocin*(xmr(1)-xmr(2))/2.d0
  dbrmr12=dbrmr12+vsuminenv*omegd(1)
  dbrmrm2=2.d0/5.d0*xmr(m-1)/2.d0*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*exp(rb(m-1))/(2.d0)**(1.d0/3.d0)*omegd(m)
  btota2=btota2+dbrmr12+dbrmrm2
!  write(*,*) 'DIFFOM APRES TRIDIAGO: btota,xo(1),vsuminenv:',btota2,omegd(1),vsuminenv
!  write(*,*) 'APRES TRIDIAGO: btota2',btota2
!  write(*,*)'rapport,difference:',btota2/btota,btota2-btota
!  write(*,*)'xldoex*dzeit:',xldoex*10000.d0*dzeit


  return

end subroutine diffom

!***********************************************************************
end module diffusion
