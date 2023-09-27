module ionisation

  use evol,only: kindreal
  use const,only: cst_k,cst_e

!  ELEMENTS INCLUS:
!    CHI(1)  H
!        2   He
!        3   C
!        4   O
!        5   Ne
!        6   Mg

  implicit none

!  DECLARATION DES CONSTANTES: IATOMS=NOMBRES D'ATOMS INCLUS
!                              IONSTATES=NOMBRE D'ETATS ION. AU MAX.
  integer, parameter:: iatoms=6,ionstates=12

!  POID ATOMIQUE
  real(kindreal), dimension(iatoms), parameter:: a_ion=(/1.00794d0,4.002602d0,12.0107d0,15.9994d0,20.1797d0,24.305d0/)
!  NOMBRE ATOMIQUE
  integer, dimension(iatoms), parameter:: iz=(/1,2,6,8,10,12/)
  real(kindreal), dimension(iatoms*ionstates), parameter:: &
!  POTENTIEL D'IONISATION (EN eV)
! Reference: Handbook of Chemistry and Physics, 84th edition (2003-2004)
    chi1=(/13.59844d0,24.58741d0,11.26030d0,13.61806d0,21.5646d0,7.64624d0,0.d0,54.41778d0,24.38332d0,35.1173d0,40.96328d0,&
           15.03528d0,0.d0,0.d0,47.8878d0,54.9355d0,63.45d0,80.1437d0,0.d0,0.d0,64.4939d0,77.41353d0,97.12d0,109.2655d0,0.d0,&
           0.d0,392.0870d0,113.8990d0,126.21d0,141.27d0,0.d0,0.d0,489.99334d0,138.1197d0,157.93d0,186.76d0,0.d0,0.d0,0.d0,&
           739.29d0,207.2759d0,225.02d0,0.d0,0.d0,0.d0,871.4101d0,239.0989d0,265.96d0,0.d0,0.d0,0.d0,0.d0,1195.8286d0,328.06d0,&
           0.d0,0.d0,0.d0,0.d0,1362.1995d0,367.5d0,0.d0,0.d0,0.d0,0.d0,0.d0,1761.805d0,0.d0,0.d0,0.d0,0.d0,0.d0,1962.665d0/),&
!  LOG DU RAPPORT DES FONCTIONS DE PARTITION ENTRE LES ETATS SUCCESSIFS
    vlu1=(/0.d0,0.602059991d0,0.124938737d0,-.051152522d0,1.079181246d0,0.602059991d0,0.d0,0.d0,-.477121255d0,0.653212514d0,&
           0.477121255d0,0.d0,0.d0,0.d0,0.602059991d0,0.124938737d0,-.051152522d0,1.079181246d0,0.d0,0.d0,0.d0,-0.477121255d0,&
           0.653212514d0,0.477121255d0,0.d0,0.d0,0.602059991d0,0.602059991d0,0.124938737d0,-.051152522d0,0.d0,0.d0,0.d0,0.d0,&
           -0.477121255d0,0.653212514d0,0.d0,0.d0,0.d0,0.602059991d0,0.602059991d0,0.124938737d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
           -0.477121255d0,0.d0,0.d0,0.d0,0.d0,0.602059991d0,0.602059991d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
           0.d0,0.602059991d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)
  real(kindreal), dimension(iatoms,0:ionstates-1), parameter:: chi = reshape (chi1, (/iatoms,ionstates/)),&
       vlu = reshape (vlu1, (/iatoms,ionstates/))

  integer, save:: ionized
  integer, dimension(iatoms), save:: list
  real(kindreal), save:: vmyion
  real(kindreal), dimension(iatoms), save:: abond,vnu
  real(kindreal), dimension(iatoms,0:ionstates), save:: xion

private
public :: ionpart
public :: iatoms,ionized,abond,list,a_ion,vnu,iz,vmyion

contains
!======================================================================
subroutine ionpart(p,t)
!-----------------------------------------------------------------------
! IONPART calcule l'ionisation partielle pour plusieurs elements.
! Entrees : p : log pression
!           t : log temperature
!
! Donnees atomiques pour les elements
!     iatoms : nombre d'elements consideres
!            : 6 (H, He, C, O, Ne, Mg)
!     ionstates : nombre maximum d'etats d'ionisation consideres
!            : 12
! Pour l'element i et l'etat d'ionisation j :
!     chi(i,j) : potentiels d'ionisation, en eV
!              : initialises dans block data ionpot
!     vlu(i,j) : log de 2 x rapport des fonctions de partition entre 2
!                etats d'ionisation successifs
!     a_ion(i)   : poids atomique
!     iz(i)  : nombre atomique (Zi)
!     abond(i) : abondance (fraction de masse)
!     vnu(i) : nombre relatif d'atomes de l'espece i
!     list(i)  : nombre d'elements reellement consideres
!     xion(i,j) : fraction d'ionisation pour l'element et l'etat
!                 consideres
!     vmyion : poids moleculaire si la matiere est totalement ionisee
!     ionized : = 0 si la matiere n'est pas totalement ionisee
!               = 1 si la matiere est totalement ionisee
! Variables locales :
!     e : nombre d'electrons par atome ("atome" = atome + ion)
!         liberes par l'ionisation
!     d, gi, h3, h4, h5, vngp, vngpp : temporaires
!     y :
! Derniere version : 2 septembre 1992
!-----------------------------------------------------------------------
  use const,only: cstlg_a
  use strucmod,only: vmion,vmol,beta_env,chem,ychem,x_env,cp,vna,vmionp,vmiont

  implicit none

  real(kindreal),intent(in):: p,t

  integer:: l,i,n,j,k
  real(kindreal):: phi,e,d,h3,h4,vng,vngp,vngpp,temp,gi,h5
  real(kindreal), dimension(iatoms,0:ionstates-1):: y
  real(kindreal), dimension(0:ionstates-1):: sum
!------------------------------------------------------------------------
! Initialisation

! Debut du calcul avec le nombre d'electrons par atome, e, correspondant
! a He ou au prochain element totalement ionise.
!  Si on a deja calcule l'ionisation on prendra le VMION comme meilleur
!  depart pour e.
  if (vmion /= 0.d0 .and. vmion /= vmol) then
    e=vmol/vmion-1.d0
  else
    e=0.d0
    do i=1,iatoms
     if (list(i) >= 2) then
       e=vnu(list(i))*abond(list(i))/a_ion(list(i))
     endif
     if (e /= 0.d0) exit
    enddo
  endif
!-----------------------------------------------------------------------
! beta_env = 1 + (aT^4/3P) = Pgaz/Ptotale = Pgaz/(Pgaz+Prad) (eq. B-5)
!        Ptotale = Pgaz + Prad = Pelect + Pionique + Prad

  beta_env= 1.d0-10.d0**(cstlg_a-log10(3.d0)+4.d0*t-p)
  beta_env=max(beta_env,1.0d-4)

! Eq. (B-6) : ionisation totale de H et He :

!      Dans les regions exterieures, l'ionisation de H et He peut etre +-
!      complete. La Pelect, le poids moleculaire moyen de la matiere
!      ionisee, ainsi que les quantites thermodynamiques, dependent des
!      degres d'ionisation des elements, principalement de l'hydrogene
!      et de l'helium.
! Lorsque la matiere n'est pas totalement ionisee et que le gaz
! d'electrons n'est pas degenere, l'equation classique de SAHA traduisant
! le nombre relatif d'atomes de l'element i entre 2 etats successifs
! d'ionisation est valable.

! Si la matiere n'est pas totalement ionisee (ionized = 0) :
  if (ionized == 0) then
    call saha(p,t,beta_env,e)
  endif

! L'equation classique de Saha n'est plus valable dans l'interieur
! profond de l'etoile ou les elements chimiques les plus abondants sont
! totalement ionises a cause des interactions electrostatiques (pressure
! ionisation).
  if (chem+ychem < 0.95d0) then
    if (p-t+log10(beta_env)-log10(1.d0+e) > 5.3447d0 .or. ionized == 1) then
      e=0.d0
      do n=1,iatoms
       if (list(n) == 0) then
         exit
       endif
       i=list(n)
       do j=0,iz(i)-1
        xion(i,j)=0.d0
       enddo
       xion(i,iz(i))=1.d0
       e=e+iz(i)*vnu(i)
      enddo
    endif
  else
    if (t-0.123d0*p > 4.09d0 .or. ionized == 1) then
! Si log T - 0.123 P > 4.09, alors l'ionisation par pression est importante,
!                            et on admet l'ionisation totale de H et He
! Mais si log T - 0.123 P < 4.09, le calcul des degres d'ionisation est
!                            effectue par SAHA, comme si l'ionisation par
!                            pression n'avait pas lieu.
      e=0.d0
      do n=1,iatoms
       if (list(n) == 0) then
         exit
       endif
       i=list(n)
       do j=0,iz(i)-1
        xion(i,j)=0.d0
       enddo
       xion(i,iz(i))=1.d0
       e=e+iz(i)*vnu(i)
      enddo
    endif
  endif

! Mise a jour des valeurs pour H et He dans les anciennes variables
  x_env(1)=xion(1,1)
  x_env(2)=xion(2,1)
  x_env(3)=xion(2,2)
!---------------------------------------------------------------------
! Calcul des grandeurs physiques
!--------------------------------
! Constantes
  d=  1.d0+e
  h3= 4.d0*(1.d0-beta_env)/beta_env + 2.5d0
  h4= (1.d0-beta_env)*(4.d0+beta_env)*d/(beta_env*beta_env)

  vng=0.d0
  vngp=0.d0
  vngpp=0.d0
  if (ionized /= 1) then

! Si l'ionisation n'est pas complete, recherche du terme dominant
    do n=1,iatoms
     if (list(n) == 0) then
       exit
     endif
     i=list(n)
     do j=0,iz(i)-1
      sum(j)=xion(i,j)+xion(i,j+1)
     enddo
     temp=0.d0
     k=0
     do j=0,iz(i)-1
      if (sum(j) > temp) then
        temp=sum(j)
        k=j
      endif
     enddo

     temp=0.d0
     if (k == 0) then
       l=0
     else
       l=1
     endif
     do j=iz(i)-1,k-l,-1
      temp=temp+xion(i,j+1)
      y(i,j)=temp
     enddo

! Calcul pour le terme dominant K
     if (k == 0) then
       temp=y(i,0)*(1.d0-y(i,0))
       if (d*e + vnu(i)*temp /= 0.d0) then
         gi= d*e*temp / (d*e + vnu(i)*temp)
       else
         gi = 0.d0
       endif
! cst_k * 1.d-7 : k boltzmann en SI.
       phi= h3 + chi(i,0)/10.d0**t/(cst_k*1.d-7/cst_e)
     else
       temp=y(i,k)*(y(i,k-1)-y(i,k))
       phi= h3 + chi(i,k)/10.d0**t/(cst_k*1.d-7/cst_e)
       if (d*e + vnu(i)*temp /= 0.d0) then
         gi= d*e*temp / (d*e*y(i,k-1) + vnu(i)*temp)
       else
         gi = 0.d0
       endif
     endif

     h5 = vnu(i)*gi
     vng= vng + h5
     vngp= vngp + h5*phi
     vngpp= vngpp + h5*phi*phi
    enddo

  endif

! Cp dans le cas de la ionisation partielle, cf. Patenaude (B.28)
  cp=   2.5d0*d + 4.d0*h4 + vngpp
  vna=  (d + h4 + vngp/beta_env)/cp
  vmionp= vng / (beta_env*d)
  vmiont= -vngp / d
  vmion=  vmol / d

  return

end subroutine ionpart
!======================================================================
subroutine saha(p,t,beta,e)
!-------------------------------------------------------------------
! SAHA resoud l'equation de ionisation pour plusieurs elements.

! Dans les regions exterieures de l'etoile, l'ionisation des elements
! (H et He principalement) peut etre plus ou moins complete.
! La pression electronique, le poids moleculaire moyen de la matiere
! ionisee ainsi que les quantites thermodynamiques dependent du degre
! d'ionisation de H et He.
! Generalement le gaz d'electrons est non degenere dans les regions
! exterieures de l'etoile (M > 1Msoleil), et l'equation classique de
! Saha traduisant le nombre relatif d'atomes de l'element i entre 2
! etats successifs d'ionisation est valable.

! Les equations de Saha pour un element i de degres d'ionisation j
! permettent de calculer les rapports entre etats d'ionisation
! successifs :
!     x(i,j+1)/x(1,j) = K(i,j) * (1+E/E),
! Cette expression de l'equation de Saha ne tient pas compte des
! differents etats excites (non ionises) des differents atomes ou ions
! dans le melange et neglige egalement les effets d'interaction
! electrostatique (pressure ionisation).
! Dans cette approximation, les fonctions de partition sont remplacees
! par les poids statistiques des niveaux fondamentaux des atomes ou ions.

! Les equations de Saha (sous cette forme) ne sont pas valables dans
! l'interieur profond d'une etoile ou les elements les plus abondants
! sont totalement ionises a cause des interactions electrostatiques.

! Les equations de Saha sont couplees entre elles et la determination des
! degres d'ionisation necessite une procedure iterative. Voir le manuel de
! Patenaude pour la description de la methode iterative utilisee

! Entrees : p : log10 Pression
!           t : log10 Temperature
!           beta : beta habituel
!           e : premiere estimation du nombre d'electrons par atome

! Donnees atomiques et variables communes
!-----------------------------------------
!     iatoms : nombre d'elements consideres
!            : 6 (H, He, C, O, Ne, Mg)
!     ionstates : nombre maximum d'etats d'ionisation consideres
!            : 12
! Pour l'element i et l'etat d'ionisation j :
!     chi(i,j) : potentiels d'ionisation, en eV
!              : initialises dans block data ionpot
!     vlu(i,j) : log de 2 x rapport des fonctions de partition entre 2
!                etats d'ionisation successifs
!     a_ion(i)   : poids atomique
!     iz(i)  : nombre atomique (Zi)
!     abond(i) : abondance (fraction de masse)
!     vnu(i) : nombre relatif d'atomes de l'espece i
!     list(i)  : nombre d'elements reellement consideres
!     xion(i,j) : fraction d'ionisation pour l'element et l'etat
!                 consideres
!     vmyion : poids moleculaire si la matiere est totalement ionisee
!     ionized : = 0 si la matiere n'est pas totalement ionisee
!               = 1 si la matiere est totalement ionisee

! Variables locales
!-------------------
!           thet : (log10(e))/kT, en eV, avec e = 2.7182818
!           h, temp : variables temporaires
!           enew : nouvelle valeur de e
!           xi   : xion, fraction de ionisation pour l'element et
!                  l'etat consideres
!           vlk  : log10 k(1+e)/2, de l'equation de Saha
!-----------------------------------------------------------------------
  use const, only: um,lgpi,cstlg_me,cstlg_k,cstlg_h

  implicit none

  real(kindreal),intent(in):: p,t,beta

  integer:: n,i,j,k,it
  real(kindreal):: e,thet,h,diff,enew,temp
  real(kindreal),dimension(iatoms,0:ionstates):: xi,vlk

  logical::TestConv
!-----------------------------------------------------------------------
  TestConv = .false.
  thet=cst_e/(1.d-7*cst_k*um*10.d0**(t))

  it=0

  h=2.5d0*t-p-log10(beta)+ (3.d0/2.d0)*(log10(2.d0)+lgpi+cstlg_me)+(5.d0/2.d0)*cstlg_k - (3.d0*cstlg_h)

! Debut du calcul, point de depart pour l'iteration
!---------------------------------------------------
! Rq : Le calcul est arrete lorsque la somme des xion est stable
  do while (.not. TestConv)
  diff=0.d0
  enew=0.d0
  e=max(e,1.d-10)

  it=it+1

! Calcul pour les 6 elements consideres dans la liste
   do n=1,iatoms
     if (list(n) == 0) exit
     i=list(n)

! Calcul des rapports entre etats de ionisation successifs
    do j=0,iz(i)-1
     vlk(i,j)=h-chi(i,j)*thet+vlu(i,j) + log10((1.d0+e)/e)
    end do

! Calcul des nombres pour chaque etat de ionisation
    nbatom: do j=0,iz(i)
     xi(i,j)=1.d0
     temp=0.d0
     do k=j-1,0,-1
      temp=temp-vlk(i,k)
      if (temp > 12.d0) then
        xi(i,j)=0.d0
        diff=diff+abs(xion(i,j)-xi(i,j))
        xion(i,j)=xi(i,j)
        cycle nbatom
      else
        xi(i,j)=xi(i,j) + 10.d0**temp
      endif
     enddo
     temp=0.d0
     do k=j,iz(i)-1
      temp=temp+vlk(i,k)
      if (temp > 12.d0) then
        xi(i,j)=0.d0
        diff=diff+abs(xion(i,j)-xi(i,j))
        xion(i,j)=xi(i,j)
        cycle nbatom
      else
        xi(i,j)=xi(i,j) + 10.d0**temp
      endif
     enddo
     xi(i,j)=1.d0/xi(i,j)
     diff=diff+abs(xion(i,j)-xi(i,j))
     xion(i,j)=xi(i,j)
    enddo nbatom

! Calcul de la contribution de l'element i au nombre d'electrons
! libres par atome
    temp=0.d0
    do k=1,iz(i)
      temp=temp+k*xi(i,k)
    end do
    enew=enew+vnu(i)*temp
   enddo

! Valeur du nombre d'electrons libres par atome
   e=enew

! Test pour l'iteration
   if (it > 40) return
   if (diff <= 5d-3) TestConv = .true.
! Lorsque diff < 5e-3, on considere que les degres d'ionisation ont
! converge vers leurs valeurs definitives.
  enddo

  return

end subroutine saha
!======================================================================
end module ionisation
