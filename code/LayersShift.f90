module LayersShift

use evol,only: kindreal

implicit none

private
public :: fitmshift,schrit,mdotshift

contains
!======================================================================
subroutine fitmshift
! DECALAGE DE LA NUMEROTATION DES NIVEAUX LORSQU ON DECROIT FITM
!----------------------------------------------------------------------
  use evol,only: ldi,npondcouche
  use const,only: Lsol,Msol
  use caramodele,only: gms,nwmd
  use inputparam,only: fitm,ialflu,idifcon,verbose
  use strucmod,only: m,q,r,p,t,s,vp,vt,vr,vs,zensi,NPcoucheEff
  use abundmod,only: x,y3,y,xc12,xc13,xc14,xn14,xn15,xo16,xo17,xo18,xf18,xf19,xne20,xne21,xne22,xna23,xmg24,xmg25,xmg26, &
    xal26,xal27,xsi28,xprot,xneut,xbid,xbid1,vx,vy3,vy,vxc12,vxc13,vxc14,vxn14,vxn15,vxo16,vxo17,vxo18,vxf18,vxf19, &
    vxne20,vxne21,vxne22,vxna23,vxmg24,vxmg25,vxmg26,vxal26g,vxal27,vxsi28,vxprot,vxneut,vxbid,vxbid1,nbelx,abelx,vabelx
  use rotmod,only: omegi,vomegi,CorrOmega,xinert,vsuminenv
  use SmallFunc,only: exphi

  implicit none

  integer:: i,ii,i1,i2,jo,k,n,no
  real(kindreal):: frac_fitm,old_xmr1,xlini,xlfin,CorrOm,sumine,sumcin,dmine,dmin1,dmcine,dmcin1,omegm
  real(kindreal),dimension(ldi):: xmr,xl,omegacorr
  real(kindreal),dimension(npondcouche):: vitcorrige
!----------------------------------------------------------------------
  xlini=0.d0
  xlfin=0.d0
  old_xmr1 = gms*Msol*(1.d0-exp(q(1)))
  do i=1,NPcoucheEff
   vitcorrige(i)=vomegi(i)+CorrOmega(i)
  enddo

  no=0
  if (abs(fitm/exphi(q(1)) - 1.d0)  >  1.d-10 .and. verbose) then
    write(*,*) 'FITM CHANGED TO : ', fitm
  endif
  if ((fitm-exphi(q(1))) <= 1.d-10) then
    do i=1,m-1
     if ((exphi(q(i))-fitm) > 1.d-10) then
       no=no+1
     endif
    enddo
    if (no > 0) then
      write(*,'(a,i3,a)')'fitm changed: ',no,' layers removed'
      write(3,'(a,i3,a)')'fitm changed: ',no,' layers removed'
      do jo=1,no+3
       xmr(jo)=gms*Msol*(1.d0-exp(q(jo)))
      enddo

! fraction of the layer in which FITM is located that is technically outside FITM. This
! will be used later to add a fraction of the angular momentum of this layer to the envelope.
      frac_fitm = (1.d0-exp(q(no))-fitm)/(exp(q(no+1))-exp(q(no)))

! calcul du moment d'inertie de chaque coquille du modele courant
      xl(1)=2.d0/3.d0*exp(2.d0*r(1))*(xmr(1)-xmr(2))/2.d0*vitcorrige(1)
      xlini=xl(1)
      do jo=2,no+2
       xl(jo)=2.d0/3.d0*exp(2.d0*r(jo))*(xmr(jo-1)-xmr(jo+1))/2.d0*vitcorrige(jo)
       xlini=xlini+xl(jo)
      enddo
      xlini = xlini + vsuminenv*vitcorrige(1)

      vsuminenv = vsuminenv + 2.d0/3.d0*exp(2.d0*r(1))*(xmr(1)-xmr(2))/2.d0
      do jo = 2,no
        vsuminenv = vsuminenv + 2.d0/3.d0*exp(2.d0*r(jo))*(xmr(jo-1)-xmr(jo+1))/2.d0
      enddo
! We add here only a fraction (determined above) of the momentum of inertia of the first layer.
      vsuminenv = vsuminenv + 2.d0/3.d0*exp(2.d0*r(no+1))*(xmr(no)-xmr(no+2))/2.d0*frac_fitm

      m=m-no
! Rq : m n'est pas ici le nombre de couches du modele initial modifie
      do i=1,m
! Rq : Valeurs exprimees en Ln(base e)
       q(i)=q(i+no)
       xmr(i)=gms*Msol*(1.d0-exp(q(i)))
       p(i)=p(i+no)
       t(i)=t(i+no)
       r(i)=r(i+no)
       s(i)=s(i+no)
       vp(i)=vp(i+no)
       vt(i)=vt(i+no)
       vr(i)=vr(i+no)
       vs(i)=vs(i+no)

       omegi(i)=omegi(i+no)
       vomegi(i)=vomegi(i+no)

       if (i+no  <=  NPcoucheEff) then
        vitcorrige(i) = vitcorrige(i+no)
       endif

       x(i)=x(i+no)
       y(i)=y(i+no)
       y3(i)=y3(i+no)
       xc12(i)=xc12(i+no)
       xc13(i)=xc13(i+no)
       xn14(i)=xn14(i+no)
       xn15(i)=xn15(i+no)
       xo16(i)=xo16(i+no)
       xo17(i)=xo17(i+no)
       xo18(i)=xo18(i+no)
       xne20(i)=xne20(i+no)
       xne22(i)=xne22(i+no)
       xmg24(i)=xmg24(i+no)
       xmg25(i)=xmg25(i+no)
       xmg26(i)=xmg26(i+no)

       vx(i)=vx(i+no)
       vy(i)=vy(i+no)
       vy3(i)=vy3(i+no)
       vxc12(i)=vxc12(i+no)
       vxc13(i)=vxc13(i+no)
       vxn14(i)=vxn14(i+no)
       vxn15(i)=vxn15(i+no)
       vxo16(i)=vxo16(i+no)
       vxo17(i)=vxo17(i+no)
       vxo18(i)=vxo18(i+no)
       vxne20(i)=vxne20(i+no)
       vxne22(i)=vxne22(i+no)
       vxmg24(i)=vxmg24(i+no)
       vxmg25(i)=vxmg25(i+no)
       vxmg26(i)=vxmg26(i+no)

       if (ialflu == 1) then
         xc14(i)=xc14(i+no)
         xf18(i)=xf18(i+no)
         xf19(i)=xf19(i+no)
         xne21(i)=xne21(i+no)
         xna23(i)=xna23(i+no)
         xal26(i)=xal26(i+no)
         xal27(i)=xal27(i+no)
         xsi28(i)=xsi28(i+no)
         xneut(i)=xneut(i+no)
         xprot(i)=xprot(i+no)
         xbid(i)=xbid(i+no)
         xbid1(i)=xbid1(i+no)

         vxc14(i)=vxc14(i+no)
         vxf18(i)=vxf18(i+no)
         vxf19(i)=vxf19(i+no)
         vxne21(i)=vxne21(i+no)
         vxna23(i)=vxna23(i+no)
         vxal26g(i)=vxal26g(i+no)
         vxal27(i)=vxal27(i+no)
         vxsi28(i)=vxsi28(i+no)
         vxneut(i)=vxneut(i+no)
         vxprot(i)=vxprot(i+no)
         vxbid(i)=vxbid(i+no)
         vxbid1(i)=vxbid1(i+no)
       endif

       do ii=1,nbelx
        abelx(ii,i)=abelx(ii,i+no)
        vabelx(ii,i)=vabelx(ii,i+no)
       enddo

      enddo
    endif   ! no
  endif     ! fitm

! Re-initialisation de q(1) en fitm du nouveau modele :
  q(1)=log(1.d0-fitm)
  xmr(1)=gms*Msol*(1.d0-exp(q(1)))

  if (no  >  0) then
!    CorrOm = (3.d0*xlini/2.d0-exp(2.d0*r(2))*vitcorrige(2)*(xmr(1)-xmr(3))/2.d0) / ((gms*Msol*exp(q(1))+ &
!             (xmr(1)-xmr(2))/2.d0)*exp(2.d0*r(1)))
    CorrOm = (xlini - 2.d0/3.d0*exp(2.d0*r(2))*vitcorrige(2)*(xmr(1)-xmr(3))/2.d0) / &
             (2.d0/3.d0*(xmr(1)-xmr(2))/2.d0*exp(2.d0*r(1)) + vsuminenv)
! [Modif CG]
! Lors du changement de FITM, on decale aussi la valeur de la correction sur
! Omega. Cependant, il faut aussi mettre a zero la correction des couches qui
! ne sont techniquement pas touchee par la correction calculee au modele precedent.
    if (no  >  NPcoucheEff) then
      write(*,*) 'WARNING: more than ', NPcoucheEff, ' shells',' removed while changing fitm. Aborting...'
      rewind(222)
      write (222,*) nwmd,': too many shells removed in fitmshift'
      stop
    else
      do i=1,NPcoucheEff-no
       CorrOmega(i) = CorrOmega(i+no)
      enddo
      do i=NPcoucheEff-no+1,NPcoucheEff
       CorrOmega(i) = 0.d0
      enddo
    endif
    CorrOmega(1)=CorrOm-vomegi(1)
! [/Modif]

    xl(1)=2.d0/3.d0*exp(2.d0*r(1))*(xmr(1)-xmr(2))/2.d0*(vomegi(1)+CorrOmega(1))
    xl(2)=2.d0/3.d0*exp(2.d0*r(2))*(xmr(1)-xmr(3))/2.d0*(vomegi(2)+CorrOmega(2))
    xlfin = xl(1) + xl(2) + vsuminenv*(vomegi(1)+CorrOmega(1))
    if (abs(xlini/xlfin-1.d0) > 1.d-9) then
      rewind(222)
      write(222,*)nwmd,': Problem of ang.mom. conserv. in fitmshift'
      write(*,*) 'Change: ', xlini/xlfin-1.d0
      stop         'WARNING: problem with momentum conservation while changing fitm'
    endif

! Si il y a des zones convectives a la surface, il faut faire attention a la correction.
! Ici, on applique la moyenne sur les zones convectives afin de corriger la correction dans
! celles-ci.
    do jo=1,NPcoucheEff-no
     omegacorr(jo) = vomegi(jo) + CorrOmega(jo)
    enddo
    do jo = NPcoucheEff-no+1,m
     omegacorr(jo) = vomegi(jo)
    enddo

    write(3,*) 'Omega(1) avant melange convectif = ', vomegi(1)

! calcul du moment d'inertie de chaque coquille du modele courant
! xMoCinScale contient le moment cinetique de l'iteration precedente.
    do jo=2,m-1
     xinert(jo)=2.d0/3.d0*exp(2.d0*r(jo))*(xmr(jo-1)-xmr(jo+1))/2.d0
    enddo
! La correction prend en compte l'enveloppe. Il faut donc maintenant aussi en prendre compte.
    xinert(1)=2.d0/3.d0*exp(2.d0*r(1))*(xmr(1)-xmr(2))/2.d0 + vsuminenv
    xinert(m)=2.d0/5.d0*exp(2.d0*r(m-1))/2.d0**(2.d0/3.d0)*xmr(m-1)/2.d0

    if (idifcon /= 1) then
      k=m-1
      couche2 : do while (k  >=  0)
       n=0
       sumine=0.d0
       sumcin=0.d0

       do while (zensi(k)  >  0.d0)
        n=n+1
        if (n == 1.and.k < m-1) then
          dmine=xinert(k+1)
          dmcine=omegacorr(k+1)*dmine
          sumine=sumine+dmine
          sumcin=sumcin+dmcine
        endif
        if (k == m-1) then
          dmine=xinert(k+1)
          dmcine=omegacorr(k+1)*dmine
          sumine=sumine+dmine
          sumcin=sumcin+dmcine
        endif
        dmin1=xinert(k)
        dmcin1=omegacorr(k)*dmin1
        sumine=sumine+dmin1
        sumcin=sumcin+dmcin1


        if (k-1  <  0) exit couche2
        if (k-1  == 0) then
          k=0
          exit
        endif
        k=k-1
       enddo
       if (n > 0) then
         omegm=sumcin/sumine
         i1=k+1
         i2=k+n+1
         do i=i1,i2
          omegacorr(i)=omegm
         enddo
       endif
       if (k-1  <=  0) exit
       k=k-1
      enddo couche2
    endif ! idifcon

! Une fois les valeurs moyennees sur la zone convective, on calcule la correction finale.
    do jo = 1,NPcoucheEff-no
     CorrOmega(jo) = omegacorr(jo) - vomegi(jo)
    enddo
    do jo = NPcoucheEff-no+1,m
     vomegi(jo) = omegacorr(jo)
    enddo
  else
! In case fitm is increased, we remove the added momentum of inertia of the first layer from
! the envelope
    vsuminenv = vsuminenv + (exp(2.d0*r(1)) + exp(2.d0*r(2)))*(old_xmr1 - xmr(1))/3.d0
    if (vsuminenv <= 0.d0) then
      rewind(222)
      write(222,*)nwmd,': Problem of ang.mom. conserv. in fitmshift'
      stop         'WARNING: problem with momentum conservation while changing fitm'
    endif
  endif

  return

end subroutine fitmshift
!======================================================================
subroutine MdotShift (dmneed)
! DECALAGE DE LA NUMEROTATION DES NIVEAUX LORS DE PERTE DE MASSE
!----------------------------------------------------------------------
  use caramodele, only: gms,dm_lost
  use strucmod, only: m,q,r
  use rotmod,only: omegi,vomegi
  use const,only: Msol

  implicit none

  real(kindreal),intent(in):: dmneed

  integer:: i,n,i1
  real(kindreal):: fm,f,hq,zerod,I_e,L_e,I_1,L_1,omega_1_old
!----------------------------------------------------------------------
    fm=(gms-dm_lost-dmneed)/gms
    hq = log(1.d0-(1.d0-exp(q(1)))/fm)
    i=2
    omega_1_old = omegi(1)
    do while (hq > q(i))
     i=i+1
    enddo
    f = (hq-q(i-1))/(q(i)-q(i-1))
    hq=q(1)
    call interx(1,i-1,f)
    q(1) = hq
    if (i > 2) then
      n=2
      do while (i <= m)
       zerod=0.d0
       call interx(n,i,zerod)
       n = n+1
       i = i+1
      enddo
      m=n-1
    endif
    hq = log(1.d0-(1.d0-exp(q(m-1)))/fm)
    f = (hq-q(m-1))/(q(m)-q(m-1))
    hq = q(m-1)
    call interx(m-1,m-1,f)
    q(m-1) = hq
    i1=m-2
    do i = 2,i1
     q(i) = log(1.d0-fm*(1.d0-exp(q(i))))
    enddo
! In the process of mass loss, the envelope is not well accounted for (and cannot really be).
! In some case where a strong gradient of Omega is present close to the first layer,
! the interpolation can lead to high value of Omega(1), and then, very fast rotating surface
! (because the envelope is by construction rotating at the same velocity as the first layer).
! To improve a bit the situation, let's assume that:
! 1. The rescaling of the total mass does not affect the angular velocity of the envelope.
!    Due to the relocalisation of fitm in the model with the new mass, a fraction fitm
!    of DelatM is actually lost by the interior (by repositioning the first layer and
!    maybe removing some layers), and a fraction (1-fitm) of deltaM is lost by the envelope.
!    If we assume that the mass loss carries out its own angular momentum, then the angular
!    velocity  of the envelope Omega_old(1) should not change.
! 2. The angular velocity Omega(1) of the first layer of the star is interpolated by repositionning
!    q(1) at a position corresponding to fitm in the model with the new mass.
! 3. Omega_old(1) and Omega(1) have no reason to be equal. However, by construction, they should.
! 4. In order to equalise these two values, we redristribute the angular momentum content
!    of the envelope and the first layer as follows:
!    L_env = I_env * Omega_old(1) = 2/3 M_env * r(1)**2 * Omega_old(1)
!          = 2/3 M_tot exp(q(1)) * r(1)**2 * Omega_old(1)
!          (here we consider the new mass of the envelope, thus M_tot is the new mass of
!          the model.)
!    L_1 = I_1 * Omega(1) = 2/3 (M(1)-M(2))/2 * r(1)**2 * Omega(1)
!        = 2/3 M_tot * (exp(q(2))-exp(q(1)))/2 * r(1)**2 * Omega(1)
!    Omega_equal(1) = (L_env + L_1)/(I_env + I_1)
!                   = (exp(q(1)) * Omega_old(1) + (exp(q(2))-exp(q(1)))/2 * Omega(1)) /
!                     (exp(q(1)) + (exp(q(2))-exp(q(1)))/2)
    I_e = 2.d0*gms*Msol*exp(q(1))*exp(2.d0*r(1))/3.d0
    L_e = I_e * omega_1_old
    I_1 = 2.d0*gms*Msol*(exp(q(2))-exp(q(1)))/2.d0*exp(2.d0*r(1))/3.d0
    L_1 = I_1*omegi(1)
    omegi(1) = (L_e + L_1)/(I_e + I_1)
! We proceed the same way for vomegi. At this point, they should be equal (vomegi is
! set to omegi in main slightly before).
    L_1 = I_1*vomegi(1)
    vomegi(1) = (L_e + L_1)/(I_e + I_1)

  return

end subroutine mdotshift
!======================================================================
subroutine schrit
! NOUVELLE VERSION QUI RAJOUTE DES PAS AU BORD DU NOYAU CONVECTIF
!----------------------------------------------------------------------
  use evol,only: mmax
  use inputparam,only: phase,dgrp,dgrl,dgry,dgrc,dgro,dgr20,verbose
  use caramodele,only:nwmd
  use abundmod,only: x,y,xc12,xo16,xne20,nbelx,mbelx,abelx
  use strucmod,only: m,p,q,r,s
  use rotmod,only: omegi
  use SmallFunc,only: exphi

  implicit none

  integer:: i,ii,ik,k,k1,mk,jschr,ischr
  real(kindreal),parameter:: dgrmr=0.3d0,dklmr=0.15d0,dgrx=0.010d0,dklx=0.004d0,dgrom=0.1d0,dklom=0.04d0
  real(kindreal):: hhabm,hhabm1,dgrra,dklra,dklp,dkll,dkly,dklc,dklo,dkl20,hhp,hhom,hhl,hhm,hhx,hhy,hhc,hho,hh20, &
                   hho1,hh201,hhp1,hhom1,hhl1,hhm1,hhx1,hhy1,hhc1,hhmr,hhmr1,hhra,hhra1,dgrm,dklm
  real(kindreal), dimension(mbelx):: hhab,hhab1
!----------------------------------------------------------------------
  dgrm = 0.008d0
  dklm = 0.003d0

! On enleve des coquilles trop proches en rayon
  dgrra = 0.020d0
  dklra = 0.001d0
  jschr = 0
  dklp = dgrp/2.5d0
  dkll = dgrl/2.5d0
  dkly = dgry/2.5d0
  dklc = dgrc/2.5d0
  dklo = dgro/2.5d0
  dkl20 = dgr20/2.5d0

  ischr = 1
  do while (ischr > 0)
   ischr = 0
   i = 1
   do while (m > i+1)
    hhp=p(i+1)-p(i)
    hhom=abs(omegi(i+1)-omegi(i))/(abs(omegi(i))+1.d-5)
    hhl=abs(s(i+1)-s(i))
    hhm =exp(q(i+1))-exp(q(i))
    hhx=abs(x(i+1)-x(i))
    hhy=abs(y(i+1)-y(i))
    hhc=abs(xc12(i+1)-xc12(i))
    hho=abs(xo16(i+1)-xo16(i))
    hh20=abs(xne20(i+1)-xne20(i))

    hhabm=0.d0
    hhabm1=0.d0
    do ii=1,nbelx
     hhab(ii) =abs(abelx(ii,i+1)-abelx(ii,i))
     hhab1(ii)=abs(abelx(ii,i+2)-abelx(ii,i))
     if (hhab(ii)  > hhabm ) hhabm = hhab(ii)
     if (hhab1(ii) > hhabm1) hhabm1= hhab1(ii)
    enddo

    hho1=abs(xo16(i+2)-xo16(i))
    hh201=abs(xne20(i+2)-xne20(i))
    hhp1=p(i+2)-p(i)

    hhom1=abs(omegi(i+2)-omegi(i))/(abs(omegi(i))+1.d-5)
    hhl1=abs(s(i+2)-s(i))
    hhm1 =exp(q(i+2))-exp(q(i))
    hhx1=abs(x(i+2)-x(i))
    hhy1=abs(y(i+2)-y(i))
    hhc1=abs(xc12(i+2)-xc12(i))
    hhmr = hhm/exphi(q(i))
    hhmr1=hhmr*hhm1/hhm

    hhra=abs(r(i+1)-r(i))
    hhra1=abs(r(i+2)-r(i))

    if (phase>=5 .and. exphi(q(i))<0.2d0 .and. xc12(i)<1.d-4) then
      if (i >= m-1 .and. verbose) then
        print*,'do you still want Ne-steps, schritgva.f?'
      endif
      if (xne20(i) < 1.d-3) then
        dgrm=0.004d0
        dklm=0.001d0
      else if (xne20(i) < 5.d-2) then
        dgrm=0.002d0
        dklm=0.0005d0
      else
        dgrm=0.006d0
        dklm=0.001d0
      endif
    else
      dgrm=0.008d0
      dklm=0.003d0
    endif

! si l'ecart entre deux couche est trop grand, on ajoute une couche:
    if (hhmr>=dgrmr .or. hhm>=dgrm .or. hhp>=dgrp .or. hhl>=dgrl .or. hhy>=dgry .or. hhc>=dgrc .or. hhx>=dgrx .or. &
        hho>=dgro .or. hh20>=dgr20 .or. hhom>=dgrom .or. hhabm>=dgro) then
      if (hhm >= 1.d-7) then
        if (mmax > m) then
          k1=m-i
          do k=1,k1
           mk=m-k
           call interx(mk+2,mk+1,0.d0)
          enddo
          call interx(i+1,i,0.5d0)
          write(3,*) ' COUCHE AJOUTEE DANS SCHRITT, i= ',i
          ischr=1
          m=m+1
        else
          jschr=jschr+1
        endif
      endif
! si l'ecart entre deux couche est trop petit, on enleve une couche:
    else if (hhmr<dklmr .and. hhm<dklm .and. hhp<dklp .and. hhl<dkll .and. hhy<dkly .and. hhc<dklc .and. hhmr1<dgrmr .and. &
             hhm1<dgrm .and. hho<dklo .and. hh20<dkl20 .and. hho1<dgro .and. hh201<dgr20 .and. hhom<dklom .and. &
             hhom1<dgrom .and. hhp1<dgrp .and. hhl1<dgrl .and. hhabm1<dgro .and. hhy1<dgry .and. hhc1<dgrc .and. &
             i<(m-2) .and. hhx<dklx .and. hhx1<dgrx .and. hhabm<dklo) then
      k1=m-i-1
      do k=1,k1
       ik=i+k
       call interx(ik,ik+1,0.d0)
      enddo
      write(3,*) ' COUCHE ENLEVEE DANS SCHRITT, i= ',i+1
      ischr=1
      m=m-1
    endif
    i=i+1
   enddo
  enddo
  if (jschr > 0) write(3,'(/1x,a,i3,a,//)') 'Try ',jschr,'times to add a layer, but maximal number reached'
  if (m  ==  mmax) then
    rewind(222)
    write (222,*) nwmd,': Max number of shells attained'
    stop '!!!! Max number of shells attained !!!!'
  endif

  return

end subroutine schrit
!======================================================================
subroutine interx(il,ir,f)
!----------------------------------------------------------------------
  use caramodele,only: nwmd
  use abundmod,only: x,y3,y,xc12,xc13,xc14,xn14,xn15,xo16,xo17,xo18,xf18,xf19,xne20,xne21,xne22,xna23,xmg24,xmg25,xmg26, &
    xal26,xal27,xsi28,xprot,xneut,xbid,vx,vy3,vy,vxc12,vxc13,vxc14,vxn14,vxn15,vxo16,vxo17,vxo18,vxf18,vxf19,vxne20, &
    vxne21,vxne22,vxna23,vxmg24,vxmg25,vxmg26,vxal26g,vxal27,vxsi28,vxprot,vxneut,vxbid,abelx,vabelx,nbelx
  use strucmod, only: m,q,p,t,r,s,vp,vt,vr,vs
  use inputparam, only: irot,isol,ialflu
  use rotmod, only: omegi,vomegi
  use SmallFunc,only: ValInterp,OmInterp

  implicit none

  integer,intent(in)::il,ir
  real(kindreal),intent(in)::f

  integer::ii
  real(kindreal):: xloil,xvloil
!----------------------------------------------------------------------
  if (abs(f)  >=  1.0d-10) then
    q  (il) = q  (ir)+f*(q  (ir+1)-q  (ir))
    p  (il) = p  (ir)+f*(p  (ir+1)-p  (ir))
    t  (il) = t  (ir)+f*(t  (ir+1)-t  (ir))
    r  (il) = r  (ir)+f*(r  (ir+1)-r  (ir))
    s  (il) = s  (ir)+f*(s  (ir+1)-s  (ir))
    vp (il) = vp (ir)+f*(vp (ir+1)-vp (ir))
    vt (il) = vt (ir)+f*(vt (ir+1)-vt (ir))
    vr (il) = vr (ir)+f*(vr (ir+1)-vr (ir))
    vs (il) = vs (ir)+f*(vs (ir+1)-vs (ir))
![Modif 2010-11]
! Interpolation en conservant l'abondance totale
    if (ir  /=  il) then
      if (ir+2  >  m) then
        write(222,*) nwmd,'Indice greater than m in interx'
        stop 'Indice greater than m in interx'
      endif
      x(il)=ValInterp(x(ir),x(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      y3(il)=ValInterp(y3(ir),y3(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      y(il)=ValInterp(y(ir),y(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      xc12(il)=ValInterp(xc12(ir),xc12(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      xc13(il)=ValInterp(xc13(ir),xc13(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      xn14(il)=ValInterp(xn14(ir),xn14(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      xn15(il)=ValInterp(xn15(ir),xn15(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      xo16(il)=ValInterp(xo16(ir),xo16(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      xo17(il)=ValInterp(xo17(ir),xo17(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      xo18(il)=ValInterp(xo18(ir),xo18(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      xne20(il)=ValInterp(xne20(ir),xne20(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      xne22(il)=ValInterp(xne22(ir),xne22(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      xmg24(il)=ValInterp(xmg24(ir),xmg24(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      xmg25(il)=ValInterp(xmg25(ir),xmg25(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      xmg26(il)=ValInterp(xmg26(ir),xmg26(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))

      vx(il)=ValInterp(vx(ir),vx(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      vy3(il)=ValInterp(vy3(ir),vy3(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      vy(il)=ValInterp(vy(ir),vy(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      vxc12(il)=ValInterp(vxc12(ir),vxc12(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      vxc13(il)=ValInterp(vxc13(ir),vxc13(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      vxn14(il)=ValInterp(vxn14(ir),vxn14(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      vxn15(il)=ValInterp(vxn15(ir),vxn15(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      vxo16(il)=ValInterp(vxo16(ir),vxo16(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      vxo17(il)=ValInterp(vxo17(ir),vxo17(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      vxo18(il)=ValInterp(vxo18(ir),vxo18(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      vxne20(il)=ValInterp(vxne20(ir),vxne20(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      vxne22(il)=ValInterp(vxne22(ir),vxne22(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      vxmg24(il)=ValInterp(vxmg24(ir),vxmg24(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      vxmg25(il)=ValInterp(vxmg25(ir),vxmg25(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      vxmg26(il)=ValInterp(vxmg26(ir),vxmg26(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))

      if (ialflu == 1) then
        xc14(il)=ValInterp(xc14(ir),xc14(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        xf18(il)=ValInterp(xf18(ir),xf18(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        xf19(il)=ValInterp(xf19(ir),xf19(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        xne21(il)=ValInterp(xne21(ir),xne21(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        xna23(il)=ValInterp(xna23(ir),xna23(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        xal26(il)=ValInterp(xal26(ir),xal26(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        xal27(il)=ValInterp(xal27(ir),xal27(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        xsi28(il)=ValInterp(xsi28(ir),xsi28(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        xneut(il)=ValInterp(xneut(ir),xneut(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        xprot(il)=ValInterp(xprot(ir),xprot(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        xbid(il)=ValInterp(xbid(ir),xbid(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        vxc14(il)=ValInterp(vxc14(ir),vxc14(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        vxf18(il)=ValInterp(vxf18(ir),vxf18(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        vxf19(il)=ValInterp(vxf19(ir),vxf19(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        vxne21(il)=ValInterp(vxne21(ir),vxne21(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        vxna23(il)=ValInterp(vxna23(ir),vxna23(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        vxal26g(il)=ValInterp(vxal26g(ir),vxal26g(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        vxal27(il)=ValInterp(vxal27(ir),vxal27(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        vxsi28(il)=ValInterp(vxsi28(ir),vxsi28(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        vxneut(il)=ValInterp(vxneut(ir),vxneut(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        vxprot(il)=ValInterp(vxprot(ir),vxprot(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
        vxbid(il)=ValInterp(vxbid(ir),vxbid(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      endif

      do ii=1,nbelx
       abelx(ii,il)=ValInterp(abelx(ii,ir),abelx(ii,ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
       vabelx(ii,il)=ValInterp(vabelx(ii,ir),vabelx(ii,ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)))
      enddo

      if (irot /= 0 .and. isol == 0) then
 ! Interpolation en conservant le moment cinetique total
        omegi(il)=OmInterp(omegi(ir),omegi(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)),exp(2.d0*r(ir)),exp(2.d0*r(ir+2)), &
                           exp(2.d0*r(il)))
        vomegi(il)=OmInterp(vomegi(ir),vomegi(ir+2),exp(q(ir)),exp(q(ir+2)),exp(q(il)),exp(2.d0*r(ir)),exp(2.d0*r(ir+2)), &
                             exp(2.d0*r(il)))
      else
        omegi(il)=omegi(ir)
        vomegi(il)=vomegi(ir)
      endif
![/Modifs 2010-11]
    else
      x  (il) = x  (ir)+f*(x  (ir+1)-x  (ir))
      y3 (il) = y3 (ir)+f*(y3 (ir+1)-y3 (ir))
      y  (il) = y  (ir)+f*(y  (ir+1)-y  (ir))
      xc12(il) = xc12(ir)+f*(xc12(ir+1)-xc12(ir))
      xc13(il)=xc13(ir)+f*(xc13(ir+1)-xc13(ir))
      xn14(il)=xn14(ir)+f*(xn14(ir+1)-xn14(ir))
      xn15(il)=xn15(ir)+f*(xn15(ir+1)-xn15(ir))
      xo16(il)= xo16(ir)+f*(xo16(ir+1)-xo16(ir))
      xo17(il)=xo17(ir)+f*(xo17(ir+1)-xo17(ir))
      xo18(il)=xo18(ir)+f*(xo18(ir+1)-xo18(ir))
      xne20(il)=xne20(ir)+f*(xne20(ir+1)-xne20(ir))
      xne22(il)=xne22(ir)+f*(xne22(ir+1)-xne22(ir))
      xmg24(il)=xmg24(ir)+f*(xmg24(ir+1)-xmg24(ir))
      xmg25(il)=xmg25(ir)+f*(xmg25(ir+1)-xmg25(ir))
      xmg26(il)=xmg26(ir)+f*(xmg26(ir+1)-xmg26(ir))
      vx(il)=vx(ir)+f*(vx(ir+1)-vx(ir))
      vy3(il)=vy3(ir)+f*(vy3(ir+1)-vy3(ir))
      vy(il)=vy(ir)+f*(vy(ir+1)-vy(ir))
      vxc12(il)=vxc12(ir)+f*(vxc12(ir+1)-vxc12(ir))
      vxc13(il)=vxc13(ir)+f*(vxc13(ir+1)-vxc13(ir))
      vxn14(il)=vxn14(ir)+f*(vxn14(ir+1)-vxn14(ir))
      vxn15(il)=vxn15(ir)+f*(vxn15(ir+1)-vxn15(ir))
      vxo16(il)=vxo16(ir)+f*(vxo16(ir+1)-vxo16(ir))
      vxo17(il)=vxo17(ir)+f*(vxo17(ir+1)-vxo17(ir))
      vxo18(il)=vxo18(ir)+f*(vxo18(ir+1)-vxo18(ir))
      vxne20(il)=vxne20(ir)+f*(vxne20(ir+1)-vxne20(ir))
      vxne22(il)=vxne22(ir)+f*(vxne22(ir+1)-vxne22(ir))
      vxmg24(il)=vxmg24(ir)+f*(vxmg24(ir+1)-vxmg24(ir))
      vxmg25(il)=vxmg25(ir)+f*(vxmg25(ir+1)-vxmg25(ir))
      vxmg26(il)=vxmg26(ir)+f*(vxmg26(ir+1)-vxmg26(ir))

      if (ialflu == 1) then
        xf19(il)=xf19(ir)+f*(xf19(ir+1)-xf19(ir))
        xne21(il)=xne21(ir)+f*(xne21(ir+1)-xne21(ir))
        xna23(il)=xna23(ir)+f*(xna23(ir+1)-xna23(ir))
        xal26(il)=xal26(ir)+f*(xal26(ir+1)-xal26(ir))
        xal27(il)=xal27(ir)+f*(xal27(ir+1)-xal27(ir))
        xsi28(il)=xsi28(ir)+f*(xsi28(ir+1)-xsi28(ir))
        vxf19(il)=vxf19(ir)+f*(vxf19(ir+1)-vxf19(ir))
        vxne21(il)=vxne21(ir)+f*(vxne21(ir+1)-vxne21(ir))
        vxna23(il)=vxna23(ir)+f*(vxna23(ir+1)-vxna23(ir))
        vxal26g(il)=vxal26g(ir)+f*(vxal26g(ir+1)-vxal26g(ir))
        vxal27(il)=vxal27(ir)+f*(vxal27(ir+1)-vxal27(ir))
        vxsi28(il)=vxsi28(ir)+f*(vxsi28(ir+1)-vxsi28(ir))
        xneut(il)=xneut(ir)+f*(xneut(ir+1)-xneut(ir))
        vxneut(il)=vxneut(ir)+f*(vxneut(ir+1)-vxneut(ir))
        xprot(il)=xprot(ir)+f*(xprot(ir+1)-xprot(ir))
        vxprot(il)=vxprot(ir)+f*(vxprot(ir+1)-vxprot(ir))
        xc14(il)=xc14(ir)+f*(xc14(ir+1)-xc14(ir))
        vxc14(il)=vxc14(ir)+f*(vxc14(ir+1)-vxc14(ir))
        xf18(il)=xf18(ir)+f*(xf18(ir+1)-xf18(ir))
        vxf18(il)=vxf18(ir)+f*(vxf18(ir+1)-vxf18(ir))
        xbid(il)=xbid(ir)+f*(xbid(ir+1)-xbid(ir))
        vxbid(il)=vxbid(ir)+f*(vxbid(ir+1)-vxbid(ir))
      endif

      do ii=1,nbelx
       abelx(ii,il)=abelx(ii,ir)+f*(abelx(ii,ir+1)-abelx(ii,ir))
       vabelx(ii,il)=vabelx(ii,ir)+f*(vabelx(ii,ir+1)-vabelx(ii,ir))
      enddo

      if (irot /= 0 .and. isol == 0) then
        if (omegi(ir) > 0.0d0.and.omegi(ir+1) > 0.0d0) then
          xloil=log(omegi(ir))+f*(log(omegi(ir+1))-log(omegi(ir)))
          xvloil=log(vomegi(ir))+f*(log(vomegi(ir+1))-log(vomegi(ir)))
          omegi(il)=exp(xloil)
          vomegi(il)=exp(xvloil)
        else
          omegi(il)=0.0d0
          vomegi(il)=0.0d0
        endif
      else
        omegi(il)=omegi(ir)
        vomegi(il)=vomegi(ir)
      endif
    endif

  else

    q  (il) = q  (ir)
    p  (il) = p  (ir)
    t  (il) = t  (ir)
    r  (il) = r  (ir)
    s  (il) = s  (ir)
    vp (il) = vp (ir)
    vt (il) = vt (ir)
    vr (il) = vr (ir)
    vs (il) = vs (ir)
    x  (il) = x  (ir)
    y3 (il) = y3 (ir)
    y  (il) = y  (ir)
    xc12(il) = xc12(ir)
    xc13(il)=xc13(ir)
    xn14(il)=xn14(ir)
    xn15(il)=xn15(ir)
    xo16(il)=xo16(ir)
    xo17(il)=xo17(ir)
    xo18(il)=xo18(ir)

    omegi(il)=omegi(ir)
    vomegi(il)=vomegi(ir)

    vx(il)=vx(ir)
    vy3(il)=vy3(ir)
    vy(il)=vy(ir)
    vxc12(il)=vxc12(ir)
    vxc13(il)=vxc13(ir)
    vxn14(il)=vxn14(ir)
    vxn15(il)=vxn15(ir)
    vxo16(il)=vxo16(ir)
    vxo17(il)=vxo17(ir)
    vxo18(il)=vxo18(ir)
    vxne20(il)=vxne20(ir)
    xne20(il)=xne20(ir)
    xmg24(il)=xmg24(ir)
    xne22(il)=xne22(ir)
    xmg25(il)=xmg25(ir)
    xmg26(il)=xmg26(ir)
    vxne22(il)=vxne22(ir)
    vxmg24(il)=vxmg24(ir)
    vxmg25(il)=vxmg25(ir)
    vxmg26(il)=vxmg26(ir)

    if (ialflu == 1) then
      xf19(il)=xf19(ir)
      xne21(il)=xne21(ir)
      xna23(il)=xna23(ir)
      xal26(il)=xal26(ir)
      xal27(il)=xal27(ir)
      xsi28(il)=xsi28(ir)
      vxf19(il)=vxf19(ir)
      vxne21(il)=vxne21(ir)
      vxna23(il)=vxna23(ir)
      vxal26g(il)=vxal26g(ir)
      vxal27(il)=vxal27(ir)
      vxsi28(il)=vxsi28(ir)
      xneut(il)=xneut(ir)
      vxneut(il)=vxneut(ir)
      xprot(il)=xprot(ir)
      vxprot(il)=vxprot(ir)
      xc14(il)=xc14(ir)
      vxc14(il)=vxc14(ir)
      xf18(il)=xf18(ir)
      vxf18(il)=vxf18(ir)
      xbid(il)=xbid(ir)
      vxbid(il)=vxbid(ir)
    endif

    do ii=1,nbelx
     abelx(ii,il) = abelx(ii,ir)
     vabelx(ii,il)=vabelx(ii,ir)
    enddo

  endif

  return

end subroutine interx
!======================================================================
end module LayersShift
