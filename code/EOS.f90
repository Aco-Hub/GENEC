!> Module computing the EOS
!! @author Sylvia Ekstrom & Cyril Georgy
!! @version 278
!! @date 26.3.2013
!! @brief EOS computation with a mixture of perfect gas and radiation
!======================================================================
MODULE EOS

  USE evol,ONLY: kindreal
  USE const,ONLY: cst_c,rgaz,cst_me
  use strucmod, only: m,adi1
  use inputparam, only: idebug

  IMPLICIT NONE

  INTEGER,SAVE:: num

  REAL(kindreal),SAVE:: psi,pl,toni,rhe
  REAL(kindreal),SAVE:: rh,rh1,rhp,rhp1,rht,rht1
  REAL(kindreal),SAVE:: uta,rhpsi,rhpsip,rhpsit
  REAL(kindreal),SAVE:: chi,hpsi
  REAL(kindreal):: tk,pg,vermy,rhes,rhete,pe,pes,pete,ue,ues,uete
  REAL(kindreal),SAVE:: cp_nablar_timmes,gamma1_Timmes,gamma1_dichte,adi1_timmes,entropy_timmes
  PRIVATE
  PUBLIC :: dichte, invert_helm_pt, read_helm_table
  PUBLIC :: num,rh,rh1,rhp,rhp1,rht,rht1,rhe,psi,rhpsi,rhpsip,rhpsit,toni,pl,uta,&
            cp_nablar_timmes,gamma1_Timmes,gamma1_dichte,adi1_timmes,entropy_timmes

CONTAINS

  !======================================================================
                            ! TIMMES EOS !
                                !2019!

                               !EOS==1!

  ! New EOS (created by Frank Timmes) introduced into GENEC to take into
  ! account pair creation in the intent to study Pair Instabilty Supernovae.
  !
  ! This EOS use table in independant variables (Rho,T) to dodge some problem
  ! at the edges of the tables. This assures thermodynamical constistency.
  !
  ! The invert_helm_pt routine compute the Timmes (Helmoltz (free energy)) EOS
  ! with (P,T) and uses a Newton-Raphson method to obtain the corresponding Rho
  ! in the tables.
  !
  ! Input: Total Pressure P
  !        Temperature T
  !        Chemical composition (see xmass and aion/zion)
  !
  ! Output: Density Rho
  !         Derivatives
  !         Degeneracy
  !         Thermodynamic exponents
  !
  ! The subroutine "invert_helm_pt" calls the "helmeos" subroutine that resolves
  ! the EOS and uses a Newton-Raphson method to derive Rho.

!======================= TIMMES EOS EXPLICATIONS ===========================

! Given a temperature temp [K], density den [g/cm**3], and a composition
! characterized by abar (average weight) and zbar (average charge),
! this routine returns all the other thermodynamic quantities.

! Of interest is the pressure [erg/cm**3], specific thermal energy [erg/gr],
! the entropy [erg/g/K], with their derivatives with respect to temperature,
! density, abar, and zbar.

! Other quantites such the normalized chemical potential eta (plus its
! derivatives), number density of electrons and positron pair (along
! with their derivatives), adiabatic indices, specific heats, and
! relativistically correct sound speed are also returned.

! This routine assumes planckian photons, an ideal gas of ions,
! and an electron-positron gas with an arbitrary degree of relativity
! and degeneracy. The full fermi-dirac integrals and their derivatives
! with respect to eta and beta are computed to machine precision, and
! all other derivatives are analytic.

! References: Cox & Giuli (c&g) chapter 24,
!             Timmes & Arnett, apj supp. 125, 277, 1999
!             Timmes & Swesty, apj supp. 126, 501, 2000

! All the input and output variables are in the file vector_eos.dek.
! The vector name is the scaler name appended with an "_row",
! for example, temp_row(i), den_row(i), and so on.


!======================= INPUTS AND OUTPUTS ===========================


! input:
! temp     = temperature
! den      = density
! abar     = average number of nucleons per nuclei
! zbar     = average number of protons per nuclei


! output:

! pres     = total pressure
! dpresdd  = derivative of total pressure with respect to density
! dpresdt  = derivative of total pressure with respect to temperature
! dpresda  = derivative of total pressure with respect to abar
! dpresdz  = derivative of total pressure with2002 respect to zbar

! ener     = total internal energy
! denerdd  = derivative of total energy with respect to density
! denerdt  = derivative of total energy with respect to temperature
! denerda  = derivative of total energy with respect to abar
! denerdz  = derivative of total energy with respect to zbar

! entr     = total entropy
! dentrdd  = derivative of total entropy with respect to density
! dentrdt  = derivative of total entropy with respect to temperature
! dentrda  = derivative of total entropy with respect to abar
! dentrdz  = derivative of total entropy with respect to zbar



! prad     = radiation pressure
! dpraddd  = derivative of the radiation pressure with density
! dpraddt  = derivative of the radiation pressure with temperature
! dpradda  = derivative of the radiation pressure with abar
! dpraddz  = derivative of the radiation pressure with zbar

! erad     = radiation energy
! deraddd  = derivative of the radiation energy with density
! deraddt  = derivative of the radiation energy with temperature
! deradda  = derivative of the radiation energy with abar
! deraddz  = derivative of the radiation energy with zbar

! srad     = radiation entropy
! dsraddd  = derivative of the radiation entropy with density
! dsraddt  = derivative of the radiation entropy with temperature
! dsradda  = derivative of the radiation entropy with abar
! dsraddz  = derivative of the radiation entropy with zbar

! radmult  = radiation multiplier (useful for turning radiation off/on)




! xni      = number density of ions
! dxnidd   = derivative of the ion number density with density
! dxnidt   = derivative of the ion number density with temperature
! dxnida   = derivative of the ion number density with abar
! dxnidz   = derivative of the ion number density with zbar

! pion     = ion pressure
! dpiondd  = derivative of the ion pressure with density
! dpiondt  = derivative of the ion pressure with temperature
! dpionda  = derivative of the ion pressure with abar
! dpiondz  = derivative of the ion pressure with zbar

! eion     = ion energy
! deiondd  = derivative of the ion energy with density
! deiondt  = derivative of the ion energy with temperature
! deionda  = derivative of the ion energy with abar
! deiondz  = derivative of the ion energy with zbar

! sion     = ion entropy
! dsiondd  = derivative of the ion entropy with density
! dsiondt  = derivative of the ion entropy with temperature
! dsionda  = derivative of the ion entropy with abar
! dsiondz  = derivative of the ion entropy with zbar

! ionmult  = ion multiplier (useful for turning ions off/on)


! etaele   = electron chemical potential
! detadd   = derivative of the electron chem potential with density
! detadt   = derivative of the electron chem potential with temperature
! detada   = derivative of the electron chem potential with abar
! detadz   = derivative of the electron chem potential with zbar

! etapos   = positron degeneracy parameter

! xne       = number density of electrons
! dxnedd    = derivative of the electron number density with density
! dxnedt    = derivative of the electron number density with temperature
! dxneda    = derivative of the electron number density with abar
! dxnedz    = derivative of the electron number density with zbar

! xnefer    = fermi integral electron number density
! dxneferdd = derivative of the fermi electron number density with density
! dxneferdt = derivative of the fermi electron number density with temperature
! dxneferda = derivative of the fermi electron number density with abar
! dxneferdz = derivative of the fermi electron number density with zbar

! xnpfer    = fermi integral positron number density
! dxnpferdd = derivative of the fermi positron number density with density
! dxnpferdt = derivative of the fermi positron number density with temperature
! dxnpferda = derivative of the fermi positron number density with abar
! dxnpferdz = derivative of the fermi positron number density with zbar

! pele      = electron pressure
! dpeledd   = derivative of the electron pressure with density
! dpeledt   = derivative of the electron pressure with temperature
! dpeleda   = derivative of the electron pressure with abar
! dpeledz   = derivative of the electron pressure with zbar

! eele     = electron energy
! deeledd   = derivative of the electron energy with density
! deeledt   = derivative of the electron energy with temperature
! deeleda   = derivative of the electron energy with abar
! deeledz   = derivative of the electron energy with zbar

! sele     = electron entropy
! dseledd   = derivative of the electron entropy with density
! dseledt   = derivative of the electron entropy with temperature
! dseleda   = derivative of the electron entropy with abar
! dseledz   = derivative of the electron entropy with zbar


! ppos     = positron pressure
! dpposdd   = derivative of the positron pressure with density
! dpposdt   = derivative of the positron pressure with temperature
! dpposda   = derivative of the positron pressure with abar
! dpposdz   = derivative of the positron pressure with zbar

! epos     = electron energy
! deposdd   = derivative of the positron energy with density
! deposdt   = derivative of the positron energy with temperature
! deposda   = derivative of the positron energy with abar
! deposdz   = derivative of the positron energy with zbar

! spos     = electron entropy
! dsposdd   = derivative of the positron entropy with density
! dsposdt   = derivative of the positron entropy with temperature
! dsposda   = derivative of the positron entropy with abar
! dsposdz   = derivative of the positron entropy with zbar

! pep      = electron + positron pressure
! dpepdd   = derivative of the electron+positron pressure with density
! dpepdt   = derivative of the electron+positron pressure with temperature
! dpepda   = derivative of the electron+positron pressure with abar
! dpepdz   = derivative of the electron+positron pressure with zbar

! eep      = electron + positron energy
! deepdd   = derivative of the electron+positron energy with density
! deepdt   = derivative of the electron+positron energy with temperature
! deepda   = derivative of the electron+positron energy with abar
! deepdz   = derivative of the electron+positron energy with zbar

! sep      = electron + positron entropy
! dsepdd   = derivative of the electron+positron entropy with density
! dsepdt   = derivative of the electron+positron entropy with temperature
! dsepda   = derivative of the electron+positron entropy with abar
! dsepdz   = derivative of the electron+positron entropy with zbar

! elemult  = electron multiplier (useful for turning e-e+ off/on)


! eip      = ionization potential ennergy
! deipdd   = derivative of ionization energy with density
! deipdt   = derivative of ionization energy with temperature
! deipda   = derivative of ionization energy with abar
! deipdz   = derivative of ionization energy with zbar


! sip      = ionization potential ennergy
! dsipdd   = derivative of ionization energy with density
! dsipdt   = derivative of ionization energy with temperature
! dsipda   = derivative of ionization energy with abar
! dsipdz   = derivative of ionization energy with zbar

! potmult  = ionization energy multiplier (useful for turning off ionization additions)



! pcoul    = coulomb pressure correction
! coulmult = coulomb component multiplier
! dpcouldd = derivative of the coulomb pressure with density
! dpcouldt = derivative of the coulomb pressure with temperature
! dpcoulda = derivative of the coulomb pressure with abar
! dpcouldz = derivative of the coulomb pressure with zbar

! ecoul    = coulomb energy correction
! decouldd = derivative of the coulomb energy with density
! decouldt = derivative of the coulomb energy with temperature
! decoulda = derivative of the coulomb energy with abar
! decouldz = derivative of the coulomb energy with zbar

! scoul    = coulomb entropy correction
! dscouldd = derivative of the coulomb entropy with density
! dscouldt = derivative of the coulomb entropy with temperature
! dscoulda = derivative of the coulomb entropy with abar
! dscouldz = derivative of the coulomb entropy with zbar


! kt       = kerg * temperature
! beta     = dimensionless ratio of kerg*temp/me*c^2

! chit     = temperature exponent in the pressure equation of state
! chid     = density exponent in the pressure equation of state
! cv       = specific heat at constant volume
! cp       = specific heat at constant pressure
! gam1     = first adiabatic exponent
! gam2     = second adiabatic exponent
! gam3     = third adiabatic exponent
! nabad    = adiabatic gradient
! sound    = relativistically correct adiabatic sound speed
! plasg    = ratio of electrostatic to thermal energy


! dse      = thermodynamic consistency check de/dt = t*ds/dt
! dpe      = thermodynamic consistency check p = d**2 de/dd + t*dpdt
! dsp      = thermodynamic consistency check dp/dt = - d**2 ds/dd

!======================= INVERSION SUBROUTINE ===========================

  subroutine invert_helm_pt
    USE const,ONLY: cst_h,cst_mh,cst_k,cst_me,cst_c,cst_a
    USE strucmod,ONLY: j1,j,t,p,x_env,beta1,vmy1,vmyo,vmye,vmol,vna
    USE abundmod,ONLY: x,y3,y,xc12,xc13,xc14,xn14,xn15,xo16,xo17,xo18,xf18,xf19,xne20,xne21,xne22,xna23,xmg24,xmg25,xmg26, &
         xal26,xal27,xsi28,zabelx,nbelx,nbzel,nbael,abelx
    USE inputparam, ONLY: ialflu

! call Timmes_config
    ! paths1=trim(input_dir)//trim('Timmes_EOS/implno.dek')
    ! paths2=trim(input_dir)//trim('Timmes_EOS/const.dek')
    ! paths3=trim(input_dir)//trim('Timmes_EOS/vector_eos.dek')
     include 'Timmes_EOS/implno.dek'
     include 'Timmes_EOS/const.dek'
     include 'Timmes_EOS/vector_eos.dek'
      ! include '/home/seb/Recherche/GENEC/GENEC_Timmes_LESTA/GENEC_Timmes/code/Timmes_EOS/implno.dek'
      ! include '/home/seb/Recherche/GENEC/GENEC_Timmes_LESTA/GENEC_Timmes/code/Timmes_EOS/const.dek'
      ! include '/home/seb/Recherche/GENEC/GENEC_Timmes_LESTA/GENEC_Timmes/code/Timmes_EOS/vector_eos.dek'
! given the pressure, temperature, and composition
! find everything else

! it is assumed that ptot_row(j), temp_row(j), abar_row(j),
! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
! been set before calling this routine.

! on input den_row(j) contains a guess for the density,
! on output den_row(j) contains the converged density.

! To get the greatest speed advantage, the eos should be fed a
! large pipe of data to work on.


! local variables
      integer          i,ii,j_bis,jlo_save,jhi_save
      double precision den,f,df,dennew,eostol,fpmin
      parameter        (eostol = 1.0d-11, &
                        fpmin  = 1.0d-14) !déclarer dans inputparam -> converg params
! local variables for computing GENEC-friendly VARIABLESINTEGER:: iINTEGER:: iii
      real(kindreal) :: ccd,tk,psi,xsq,x_bis,wx,gol,phi1,ph1,phi2,ph2,phi1d,ph1d,phi2d,ph2d,sr,srs,sp,&
      sps,phi3,ph3,phi3d,ph3d,srt,spt,su,sus,sut,rhes,cst_mch3,ccr,ccp,hchi,hpsi,tu,wz,asr,gamma1_dichte,&
      gamma_gas
      real(kindreal), DIMENSION(12):: &
           ! a_i:
           dega=(/7.28905d-2,2.04648d-1,2.54685d-1,1.96668d-1,1.04294d-1,3.95716d-2,1.09024d-2,&
           2.21504d-3,3.14295d-4,3.64498d-5,2.42533d-6,1.49830d-7/), &
           ! u_i:
           degu=(/0.11896d0,0.47652d0,1.07476d0,1.91722d0,3.00898d0,4.35703d0,5.96967d0,7.86147d0,&
           10.02706d0,12.57936d0,15.19876d0,18.87814d0/), &
           ! alpha_i:
           degex=(/8.87844d-1,6.20940d-1,3.41379d-1,1.47015d-1,4.93418d-2,1.28164d-2,2.55507d-3,&
           3.85307d-4,4.41879d-5,3.44234d-6, 2.50763d-7,6.32888d-9/)

! Chemical composition
       integer,parameter:: ionmax=49 !x,y,xc12,xc13,xn14,xn15,xo17,xo18,xne20, &
               !  xne22,xmg24,xmg25,xmg26,xf19,xne21,xna23,xal26,xal27,xsi28 AND abelx
       double precision:: xmass(ionmax),aion(ionmax),zion(ionmax),abar,zbar
       double precision:: input_norm

! if (j1==m-1) then
! ! write(3,*)j,j1,'Timmes START: T=',exp(t(j1)),'P=',exp(p(j1)),'rho=',exp(rh),'rh1=',exp(rh1)
! endif
! Chemical composition in friendly terms to Timmes EOS

aion(1)  = 1.0d0 !H
zion(1)  = 1.0d0
!!!!!!!!!!!!!!!!!
aion(23)  = 3.0d0 !He3
zion(23)  = 2.0d0
aion(2)  = 4.0d0  !He4
zion(2)  = 2.0d0
!!!!!!!!!!!!!!!!!
aion(3)  = 12.0d0 !C
zion(3)  = 6.0d0
aion(4)  = 13.0d0 !C13
zion(4)  = 6.0d0
aion(21)  = 14.0d0 !C14
zion(21)  = 6.0d0
!!!!!!!!!!!!!!!!! :N
aion(5)  = 14.0d0
zion(5)  = 7.0d0
aion(19)  = 15.0d0
zion(19)  = 7.0d0
!!!!!!!!!!!!!!!!! !O
aion(6)  = 16.0d0
zion(6)  = 8.0d0
aion(7)  = 17.0d0
zion(7)  = 8.0d0
aion(20)  = 18.0d0
zion(20)  = 8.0d0
!!!!!!!!!!!!!!!!! !Ne
aion(8)  = 20.0d0
zion(8)  = 10.0d0
aion(9)  = 22.0d0
zion(9)  = 10.0d0
!!!!!!!!!!!!!!!!! !Mg
aion(10)  = 24.0d0
zion(10)  = 12.0d0
aion(11)  = 25.0d0
zion(11)  = 12.0d0
aion(12)  = 26.0d0
zion(12)  = 12.0d0
!!!!!!!!!!!!!!!!! !F
aion(22)  = 18.0d0
zion(22)  = 9.0d0
aion(13)  = 19.0d0
zion(13)  = 9.0d0
!!!!!!!!!!!!!!!!! !Ne
aion(14)  = 21.0d0
zion(14)  = 10.0d0
!!!!!!!!!!!!!!!!! !Na
aion(15)  = 23.0d0
zion(15)  = 11.0d0
!!!!!!!!!!!!!!!!! !Al
aion(16)  = 26.0d0
zion(16)  = 13.0d0
aion(17)  = 27.0d0
zion(17)  = 13.0d0
!!!!!!!!!!!!!!!!! !Si
aion(18)  = 28.0d0
zion(18)  = 14.0d0
!!!!!!!!!!!!!!!!!
aion(24:ionmax) = nbael(:)
zion(24:ionmax) = nbzel(:)




xmass(1) = x(j1)!!!!!!!!!!!!!!!!!!
xmass(23)= y3(j1)
xmass(2) = y(j1)
xmass(3) = xc12(j1)
xmass(4) = xc13(j1)
xmass(21)= xc14(j1)
xmass(5) = xn14(j1)
xmass(19) = xn15(j1)
xmass(6) = xo16(j1)
xmass(7) = xo17(j1)
xmass(20) = xo18(j1)
xmass(8) = xne20(j1)
xmass(9) = xne22(j1)
xmass(10) = xmg24(j1)
xmass(11) = xmg25(j1)
xmass(12) = xmg26(j1)
xmass(22) = xf18(j1)
xmass(13) = xf19(j1)
xmass(14) = xne21(j1)
xmass(15) = xna23(j1)
xmass(16) = xal26(j1)
xmass(17) = xal27(j1)
xmass(18) = xsi28(j1)
xmass(24:ionmax) = abelx(:,j1)

! average atomic weight and charge

!Renormalise the mass fraction to 1

      

      abar   = 1.0d0/sum(xmass(1:ionmax)/aion(1:ionmax))
      zbar   = abar * sum(xmass(1:ionmax) * zion(1:ionmax)/aion(1:ionmax))


      if (j == m-1) then
         write(3,*) "ADAM eos stuff", j1
         write(3,*) "nbelx" , nbelx
         write(3,*) "aion", aion
         write(3,*) "zion", zion
         write(3,*) "abar", abar
         write(3,*) "zbar", zbar
         write(3,*) "rho", exp(rh1)
         write(3,*) "temp", exp(t(j1))
         write(3,*) "ptot", exp(p(j1))
      endif
! set the input vector. pipeline is only 1 element long
      abar_row(1) = abar
      zbar_row(1) = zbar
      jlo_eos = 1
      jhi_eos = 1
      !write(*,*),j, 'T=',t(j),'P=',p(j),'rho=',rh,'rh1=',rh1

! set the Temperature and pressure in friendly term to Timmes EOS
      den_row(1)  = exp(rh)   !0.44 !!! Initialisation pour permettre au Newton-Raphson de converger.
      temp_row(1) = exp(t(j1))
      ptot_row(1) = exp(p(j1))
      !write(*,*)j,t(j),p(j),rh1
! initialize
      jlo_save = jlo_eos
      jhi_save = jhi_eos
      do j_bis=jlo_eos, jhi_eos
       eoswrk01(j_bis) = 0.0d0
       eoswrk02(j_bis) = 0.0d0
       eoswrk03(j_bis) = ptot_row(j_bis)
       eoswrk04(j_bis) = den_row(j_bis)
      end do


! do the first newton loop with all elements in the pipe
      call helmeos

      do j_bis = jlo_eos, jhi_eos

       f     = ptot_row(j_bis)/eoswrk03(j_bis) - 1.0d0
       df    = dpd_row(j_bis)/eoswrk03(j_bis)
       eoswrk02(j_bis) = f/df

! limit excursions to factor of two changes
       den    = den_row(j_bis)
       dennew = min(max(0.5d0*den,den - eoswrk02(j_bis)),2.0d0*den)

! compute the error
       eoswrk01(j_bis)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
       den_row(j_bis)  = min(1.0d14,max(dennew,1.0d-11))
      enddo


! now loop over each element of the pipe individually
      do j_bis = jlo_save, jhi_save

       do i=2,100


! check for convergence
   if (idebug .gt. 0) then
         if (eoswrk01(j_bis) .lt. eostol ) then
            write(3,*) 'converged in invert_helm_pt'
            write(3,*) 'j_bis =',j_bis,'  eoswrk01(j_bis) =',eoswrk01(j_bis)
            write(3,*) 'eostol =',eostol,'f = ',f,'df = ',df
            write(3,*) 'Other conditions',eoswrk02(j_bis),fpmin
         endif

         if (abs(eoswrk02(j_bis)) .le. fpmin) then
            write(3,*) 'converged in invert_helm_pt'
            write(3,*) 'j_bis =',j_bis,'  eoswrk02(j_bis) =',eoswrk02(j_bis)
            write(3,*) 'fpmin =',fpmin
         endif
   endif
      !   if (eoswrk01(j_bis) .lt. eostol .or. &
      !       abs(eoswrk02(j_bis)) .le. fpmin) goto 20
      ! if (eoswrk01(j_bis) .lt. eostol )  then
      !    write(3,*) 'converged in invert_helm_pT AFTER', i 
      ! ENDIF
      if (eoswrk01(j_bis) .lt. eostol ) goto 20

        !Relax tolerance if needed.
        if ( i .ge. 20 .and. i .le. 40 ) then
            write(3,*) "struggle to converg 1 ", i, eoswrk01(j_bis), 10.d0*eostol
           if (eoswrk01(j_bis) .lt. 10.d0 * eostol ) goto 20
        else if ( i .ge. 40 .and. i .le. 60 ) then
           write(3,*) "struggle to converge 2 ", i, eoswrk01(j_bis), 100.d0*eostol
           if (eoswrk01(j_bis) .lt. 100.d0 * eostol ) goto 20
        else if ( i .ge. 60 .and. i .le. 80 ) then 
             write(3,*) "struggle to converge 3 ", i, eoswrk01(j_bis), 1000.d0*eostol
            if (eoswrk01(j_bis) .lt. 1000.d0 * eostol ) goto 20
         else if ( i .ge. 80 ) then
            write(3,*) "struggle to converge 4", i, eoswrk01(j_bis), 10000.d0*eostol
            if (eoswrk01(j_bis) .lt. 10000.d0 * eostol ) goto 20
        endif

        jlo_eos = j_bis
        jhi_eos = j_bis

        call helmeos

        f     = ptot_row(j_bis)/eoswrk03(j_bis) - 1.0d0
        df    = dpd_row(j_bis)/eoswrk03(j_bis)
        eoswrk02(j_bis) = f/df

! limit excursions to factor of two changes
        den    = den_row(j_bis)
        dennew = min(max(0.5d0*den,den - eoswrk02(j_bis)),2.0d0*den)

! compute the error
        eoswrk01(j_bis)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
        den_row(j_bis)  = min(1.0d14,max(dennew,1.0d-11))

! end of netwon loop
       end do


! we did not converge if we land here
      write(6,*)
      write(6,*) 'newton-raphson failed in routine invert_helm_pt in shell',j
      write(6,*) 'pipeline element',j_bis
      write(6,01) 'pwant  =',eoswrk03(j_bis),' temp =',temp_row(j_bis)
 01   format(1x,5(a,1pe16.8))
      write(6,01) 'error =',eoswrk01(j_bis), &
                  '  eostol=',eostol,'  fpmin =',fpmin
      write(6,01) 'den   =',den_row(j_bis),'  denold=',eoswrk04(j_bis)
      write(6,01) 'f/df  =',eoswrk02(j_bis),' f   =',f,    ' df    =',df
      write(6,*)
      !write(*,*), 'T=',t(j),'P=',p(j),'rho=',rh,'rh1=',rh1

      stop 'could not find a density in routine invert_helm_pt'

!gkorm trop grand
!flag pour sortir de EOS puis henyey reconnu comme gkorm > gkorm

! land here if newton loop converged, back for another pipe element
 20    continue
      end do



! call eos one more time with the converged value of the density

      jlo_eos = jlo_save
      jhi_eos = jhi_save

      call helmeos

!===== Part dedicated to transform Timmes EOS output into GENEC variables ======

! rhe and rhes computing for degenerate case

      !!!!! Définition des constantes !!!!!!
      cst_mch3 = (cst_me*cst_c/cst_h)**3.d0
      ccd = pi**2.d0 / 6.d0
      ccr = 8.d0*pi*cst_mh*cst_mch3
      ccp = (16.d0*pi/6.d0)*cst_me*cst_c**2.d0*cst_mch3
      toni=temp_row(1) !       toni=EXP(t(j))  !the non-log temperature value
      !write(*,*)'TONI = ', toni
      tk=toni*cst_k/(cst_me*cst_c**2.d0) !T_k=kT/m_ele*c²
      pl=EXP(p(j1))

      !!!!!!!!! Psi computed by Timmes !!!!!
      psi=etaele_row(1)
      hchi=exp(-psi)
      hpsi = psi


      !!!!!! FULLY DEGENERATE !!!!
      IF ((hchi-1.d0/EXP(7.d0)) >= 0.d0) THEN  !PATENAUDE 1974 p.52

         sr=0.d0
         srs=0.d0
         ! PATENAUDE 1974, p. 51 sqq:
         DO i=1,12
            tu=tk*degu(i)
            wz=SQRT(tk*(2.d0+tu))
            asr=(1.d0+tu)*tk*wz*dega(i)/(hchi+degex(i))
            sr=sr+asr
            srs=srs-asr/(hchi+degex(i))
         ENDDO

      !!!!!! PARTIALLY DEGENERATE !!!!
      ELSE
        xsq=2.d0*tk*psi+(tk*psi)**2.d0
        x_bis= SQRT(abs(xsq))
        wx=1.d0+tk*psi
        gol=LOG(wx-x_bis)
        phi1=xsq*x_bis/3.d0
        ph1=tk*tk*(1.d0+2.d0*xsq)/x_bis
        IF ((x_bis-0.5d0) > 0.d0) THEN
          phi2=x_bis*wx*(2.d0*xsq-3.d0)*0.125d0-0.375d0*gol
        ELSE
          phi2=0.125d0*x_bis*xsq*xsq*(8.d0/5.d0-xsq*(4.d0/7.d0-xsq*(1.d0/3.d0-xsq*5.d0/22.d0)))
        ENDIF
        ph2=tk**2.d0*3.d0*x_bis*wx
        phi1d=x_bis*wx
        ph1d=(2.d0*xsq-1.d0)*tk*tk*wx/(xsq*x_bis)
        phi2d=x_bis*xsq
        ph2d=(6.d0*xsq+3.d0)*tk*tk/x_bis
        sr=phi1+ccd*ph1
        srs=-(phi1d+ccd*ph1d)*tk/hchi
        sp=phi2+ccd*ph2
        sps=-(phi2d+ccd*ph2d)*tk/hchi
        IF ((hpsi+1.d0) >= 0.d0) THEN
         IF ((x_bis-0.5d0) > 0.d0) THEN
            phi3=x_bis*wx*(2.d0*xsq+1.d0)*0.125d0-phi1+0.125d0*gol
         ELSE
            phi3=x_bis*xsq*xsq*(0.1d0-xsq*(1.d0/56.d0-xsq*(1.d0/144.d0-xsq*5.d0/1408.d0)))
         ENDIF
         ph3=tk*tk*(wx*(1.d0+3.d0*xsq)-1.d0-2.d0*xsq)/x_bis
         phi3d=x_bis*wx*(wx-1.d0)
         ph3d=tk*tk*(xsq*(6.d0*xsq+3.d0)-1.d0-wx*(2.d0*xsq-1.d0))/(xsq*x_bis)
         srt=(phi1d+ccd*ph1d)*psi+2.d0*ccd/tk*ph1
         spt=(phi2d+ccd*ph2d)*psi+2.d0*ccd/tk*ph2
         su=phi3+ccd*ph3
         sus=-(phi3d+ccd*ph3d)*tk/hchi
         sut=(phi3d+ccd*ph3d)*psi+2.d0*ccd/tk*ph3
        ENDIF
      ENDIF

      rhe=ccr*sr
      !write(*,*)'RHE = ',rhe,'ccr=',ccr,'sr=',sr
      rhes=srs/sr
      !write(*,*),j1,'TIMMES: ','rhe=',rhe,'rhes=',rhes

! MODIF TO COMPUTE num, it trace when we are in partial ionisation. It was computed in degen before.
          num = -1000
          IF (j1 == 1) THEN
             num=0
          ENDIf
          IF (x_env(3) /= 0.d0) THEN
             num=j1
          ENDIf

! MODIF TO COMPUTE mu, as it is done in dichte (vmyo,vmy1,vmyhelio,...).

!-----------------------------------------------------------------------
vmy1 = 2.d0*x(j1)+y3(j1)+3.d0/4.d0*y(j1)+7.d0/12.d0*xc12(j1)+7.d0/13.d0*xc13(j1)+8.d0/14.d0*xn14(j1)+8.d0/15.d0*xn15(j1)+ &
     9.d0/16.d0*xo16(j1)+9.d0/17.d0*xo17(j1)+9.d0/18.d0*xo18(j1)+11.d0/20.d0*xne20(j1)+11.d0/22.d0*xne22(j1)+ &
     13.d0/24.d0*xmg24(j1)+13.d0/25.d0*xmg25(j1)+13.d0/26.d0*xmg26(j1)

vmyo = x(j1)+y3(j1)/3.d0+y(j1)/4.d0+xc12(j1)/12.d0+xc13(j1)/13.d0+xn14(j1)/14.d0+xn15(j1)/15.d0+xo16(j1)/16.d0+ &
     xo17(j1)/17.d0+xo18(j1)/18.d0+xne20(j1)/20.d0+xne22(j1)/22.d0+xmg24(j1)/24.d0+xmg25(j1)/25.d0+xmg26(j1)/26.d0

vmye=0.5d0*(1.d0+x(j1))+y3(j1)/6.d0

IF (ialflu == 1) THEN
   vmy1 = vmy1+0.5d0*xc14(j1)+10.d0/18.d0*xf18(j1)+10.d0/19.d0*xf19(j1)+11.d0/21.d0*xne21(j1)+12.d0/23.d0*xna23(j1)+ &
        14.d0/26.d0*xal26(j1)+14.d0/27.d0*xal27(j1)+15.d0/28.d0*xsi28(j1)
   vmyo = vmyo+xc14(j1)/14.d0+xf18(j1)/18.d0+xf19(j1)/19.d0+xne21(j1)/21.d0+xna23(j1)/23.d0+xal26(j1)/26.d0+ &
        xal27(j1)/27.d0+xsi28(j1)/28.d0
   vmye = x(j1)+2.d0/3.d0*y3(j1)+0.5d0*(y(j1)+xc12(j1)+xn14(j1)+xo16(j1)+xne20(j1)+xmg24(j1)+xal26(j1)+xsi28(j1)+xf18(j1))+ &
        6.d0/13.d0*xc13(j1)+7.d0/15.d0*xn15(j1)+8.d0/17.d0*xo17(j1)+8.d0/18.d0*xo18(j1)+10.d0/22.d0*xne22(j1)+ &
        12.d0/25.d0*xmg25(j1)+12.d0/26.d0*xmg26(j1)+6.d0/14.d0*xc14(j1)+9.d0/19.d0*xf19(j1)+10.d0/21.d0*xne21(j1)+ &
        11.d0/23.d0*xna23(j1)+13.d0/27.d0*xal27(j1)
ENDIF

!  z elements taken ~ as Ca56: A=56 but (nbzel(ii)+1.)/nbael(ii)-->0.5
vmy1= vmy1+0.5d0*zabelx
vmyo= vmyo+zabelx/56.d0
IF (ialflu == 1) THEN
   vmye = x(j1)+2.d0/3.d0*y3(j1)+0.5d0*(y(j1)+xc12(j1)+xn14(j1)+xo16(j1)+xne20(j1)+xmg24(j1)+xal26(j)+xsi28(j1)+ &
        xf18(j1))+6.d0/13.d0*xc13(j1)+6.d0/14.d0*xc14(j1)+7.d0/15.d0*xn15(j1)+8.d0/17.d0*xo17(j1)+8.d0/18.d0*xo18(j1)+ &
        9.d0/19.d0*xf19(j1)+10.d0/21.d0*xne21(j1)+10.d0/22.d0*xne22(j1)+11.d0/23.d0*xna23(j1)+12.d0/25.d0*xmg25(j1)+ &
        12.d0/26.d0*xmg26(j1)+13.d0/27.d0*xal27(j1)+0.5d0*zabelx
ELSE
   vmye = x(j1)+2.d0/3.d0*y3(j1)+0.5d0*(y(j1)+xc12(j1)+xn14(j1)+xo16(j1)+xne20(j1)+xmg24(j1))+6.d0/13.d0*xc13(j1)+ &
        7.d0/15.d0*xn15(j1)+8.d0/17.d0*xo17(j1)+8.d0/18.d0*xo18(j1)+10.d0/22.d0*xne22(j1)+12.d0/25.d0*xmg25(j1)+ &
        12.d0/26.d0*xmg26(j1)+0.5d0*zabelx
ENDIF

DO ii=1,nbelx
   vmy1= vmy1+abelx(ii,j1)*REAL(nbzel(ii)+1)/REAL(nbael(ii))
   vmyo= vmyo+abelx(ii,j1)/REAL(nbael(ii))
   vmye= vmye+abelx(ii,j1)*REAL(nbzel(ii))/REAL(nbael(ii))
ENDDO

vmy1=-LOG(vmy1)
vmyo=1.d0/vmyo
vmye=1.d0/vmye
vermy=vmye/vmyo
vmol=vmyo


! attributing Timmes EOS output to GENEC VARIABLES

            rh1=log(den_row(1)) ! Density obtained from (P,T)
            rhp1=(ptot_row(1)/den_row(1))*(1./dpd_row(1)) !! dln(Rho)/dln(P) = P/rho * 1/(dP/dRho) (ptot_row(1)/den_row(1))*
            rht1= -(exp(t(j)) / den_row(1)) * (dpt_row(1)/dpd_row(1)) !! dln(Rho)/dln(T) = -T/Rho * (dP/dT)/(dP/dRho)
            rhpsi= -exp(-psi)*rhes
            rhpsit=detadt_row(1)*temp_row(1)*exp(-psi)*rhes !! rhpsit = dPsi/dT * T * exp(-Psi) * rhes
            rhpsip=-1.*detadd_row(1)*rhp1*rh1*exp(-psi)*rhes !! rhpsip = -dPsi/dRho * dln(Rho)/dln(P) * rho * exp(-Psi) * rhes
            !MISTAKE HERE WAS WRITTEN rhp not rhp1
            cp_nablar_timmes=cp_row(1)
            adi1_timmes = nabad_row(1)

            beta1=1.d0- EXP(LOG(cst_a)-LOG(3.d0)+4.d0*t(j1)-p(j1))
            !write(*,*)'beta1=',beta1,'cst_a',cst_a,'t(j1)',t(j1),'p(j1)',p(j1)
            beta1=MAX(beta1,1.d-5)
            gamma1_Timmes=gam1_row(1)
            entropy_timmes = stot_row(1)
            !write(*,*)'beta1 = ',beta1
            !write(*,*)j,'Timmes END: rh1=',rh1
            !write(*,*)'TONI end = ', toni
            ! write(*,*)'Gamma 1 = ', j, gam1_row(1)-4./3.
            !ga1
            gamma_gas=32-24*beta1-3*beta1**2/(3*beta1*(8-7*beta1))
            gamma1_dichte=beta1+((4-3*beta1)**2*(gamma_gas-1)/(beta1+12*(gamma_gas-1)*(1-beta1)))


            if (j==m-1) then
               write(3,*) "SECOND EOS STUFF", j1
               write(3,*) "rho", den_row(1),exp(rh1)
               write(3,*) "T", temp_row(1)
               write(3,*) "Ptot", ptot_row(1)
               write(3,*) "abar", abar_row(1)
               write(3,*) "zbar", zbar_row(1)
            endif
!!!!! INTRODUCTION S et U
            ! srad=(4*a*t1**4/rh1)*(1/t1)
            ! sion=((avo*kt/abar)+(1.5*avo*kt/(abar*rh1)))*(1/t1)*((k*avo)/(abar*abar*kt))
            ! sele=-df_t*(max(1.0d-16,zbar/abar))
            ! plasg=(zbar**2)*(qe**2)*(1/kt1)*(4/3.*pi*avo*rh1/abar)**(1/3.)
            ! x=plasg**0.25
            ! scoul=(-1/(abar*kt))*(3.0*b1*(zbar**2*qe**2/kt1)
!!!!!!!!!!!!!!!

            ! write(*,*)'Gamma 1: ', j, gam1_row(1)
            ! write(*,*)'Gamma 1 GENEC way: ', j, gamma1_dichte
            ! write(*,*)'dp/drho',j,(ptot_row(1)/den_row(1))*(1/rhp1)
            ! if (gam1_row(1) <= 4./3.) then
            !   write(*,*)'Gamma 1 <= 4/3: ', j, gam1_row(1)
            !   write(*,*)'Gamma 1 <= 4/3 GENEC way: ', j, gamma1_dichte
            !   write(*,*) 'Eta_pos = ',etapos_row(1)
            ! endif
            !   if (ppos_row(1)>0.0) then
            !   write(*,*)'Positron pressure = ', j, ppos_row(1)
            ! endif

            ! endif
      return
      end


!======================= TIMMES EOS SUBROUTINE ===========================


      subroutine helmeos
        ! include '/home/seb/Recherche/GENEC/GENEC_Timmes_LESTA/GENEC_Timmes/code/Timmes_EOS/implno.dek'
        ! include '/home/seb/Recherche/GENEC/GENEC_Timmes_LESTA/GENEC_Timmes/code/Timmes_EOS/const.dek'
        ! include '/home/seb/Recherche/GENEC/GENEC_Timmes_LESTA/GENEC_Timmes/code/Timmes_EOS/vector_eos.dek'
        ! include '/home/seb/Recherche/GENEC/GENEC_Timmes_LESTA/GENEC_Timmes/code/Timmes_EOS/helm_table_storage.dek'
        include 'Timmes_EOS/implno.dek'
        include 'Timmes_EOS/const.dek'
        include 'Timmes_EOS/vector_eos.dek'
        include 'Timmes_EOS/helm_table_storage.dek'


        ! given a temperature temp [K], density den [g/cm**3], and a composition
        ! characterized by abar and zbar, this routine returns most of the other
        ! thermodynamic quantities. of prime interest is the pressure [erg/cm**3],
        ! specific thermal energy [erg/gr], the entropy [erg/g/K], along with
        ! their derivatives with respect to temperature, density, abar, and zbar.
        ! other quantites such the normalized chemical potential eta (plus its
        ! derivatives), number density of electrons and positron pair (along
        ! with their derivatives), adiabatic indices, specific heats, and
        ! relativistically correct sound speed are also returned.
        !
        ! this routine assumes planckian photons, an ideal gas of ions,
        ! and an electron-positron gas with an arbitrary degree of relativity
        ! and degeneracy. interpolation in a table of the helmholtz free energy
        ! is used to return the electron-positron thermodynamic quantities.
        ! all other derivatives are analytic.

      ! declare
            integer          i,j
            double precision temp,den,abar,zbar,ytot1,ye, &
                             x,y,zz,zzi,deni,tempi,xni,dxnidd,dxnida, &
                             dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt, &
                             dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt, &
                             deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt, &
                             dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion, &
                             sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd, &
                             dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp, &
                             gam1,gam2,gam3,chit,chid,nabad,sound,etaele, &
                             detadt,detadd,xnefer,dxnedt,dxnedd,s

            double precision pgas,dpgasdd,dpgasdt,dpgasda,dpgasdz, &
                             egas,degasdd,degasdt,degasda,degasdz, &
                             sgas,dsgasdd,dsgasdt,dsgasda,dsgasdz, &
                             cv_gas,cp_gas,gam1_gas,gam2_gas,gam3_gas, &
                             chit_gas,chid_gas,nabad_gas,sound_gas


            double precision sioncon,forth,forpi,kergavo,ikavo,asoli3,light2
            parameter        (sioncon = (2.0d0 * pi * amu * kerg)/(h*h), &
                              forth   = 4.0d0/3.0d0, &
                              forpi   = 4.0d0 * pi, &
                              kergavo = kerg * avo, &
                              ikavo   = 1.0d0/kergavo, &
                              asoli3  = asol/3.0d0, &
                              light2  = clight * clight)

      ! for the abar derivatives
            double precision dpradda,deradda,dsradda, &
                             dpionda,deionda,dsionda, &
                             dpepda,deepda,dsepda, &
                             dpresda,denerda,dentrda, &
                             detada,dxneda

      ! for the zbar derivatives
            double precision dpraddz,deraddz,dsraddz, &
                             dpiondz,deiondz,dsiondz, &
                             dpepdz,deepdz,dsepdz, &
                             dpresdz,denerdz,dentrdz, &
                             detadz,dxnedz

      ! for the interpolations
            integer          iat,jat
            double precision free,df_d,df_t,df_dd,df_tt,df_dt
            double precision xt,xd,mxt,mxd, &
                             si0t,si1t,si2t,si0mt,si1mt,si2mt, &
                             si0d,si1d,si2d,si0md,si1md,si2md, &
                             dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
                             dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
                             ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt, &
                             ddsi0d,ddsi1d,ddsi2d,ddsi0md,ddsi1md,ddsi2md, &
                             z,psi0,dpsi0,ddpsi0,psi1,dpsi1,ddpsi1,psi2, &
                             dpsi2,ddpsi2,din,h5,fi(36), &
                             xpsi0,xdpsi0,xpsi1,xdpsi1,h3, &
                             w0t,w1t,w2t,w0mt,w1mt,w2mt, &
                             w0d,w1d,w2d,w0md,w1md,w2md


      ! for the uniform background coulomb correction
            double precision dsdd,dsda,lami,inv_lami,lamida,lamidd, &
                             plasg,plasgdd,plasgdt,plasgda,plasgdz, &
                             ecoul,decouldd,decouldt,decoulda,decouldz, &
                             pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz, &
                             scoul,dscouldd,dscouldt,dscoulda,dscouldz, &
                             a1,b1,c1,d1,e1,a2,b2,c2,third,esqu
            parameter        (a1    = -0.898004d0, &
                              b1    =  0.96786d0, &
                              c1    =  0.220703d0, &
                              d1    = -0.86097d0, &
                              e1    =  2.5269d0, &
                              a2    =  0.29561d0, &
                              b2    =  1.9885d0, &
                              c2    =  0.288675d0, &
                              third =  1.0d0/3.0d0, &
                              esqu  =  qe * qe)


      ! quintic hermite polynomial statement functions
      ! psi0 and its derivatives
            psi0(z)   = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
            dpsi0(z)  = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
            ddpsi0(z) = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)


      ! psi1 and its derivatives
            psi1(z)   = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
            dpsi1(z)  = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
            ddpsi1(z) = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)


      ! psi2  and its derivatives
            psi2(z)   = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
            dpsi2(z)  = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
            ddpsi2(z) = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)


      ! biquintic hermite polynomial statement function
            h5(i,j,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)= &
                   fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t &
                 + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt &
                 + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t &
                 + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt &
                 + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t &
                 + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt &
                 + fi(13) *w1d*w0t   + fi(14) *w1md*w0t &
                 + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt &
                 + fi(17) *w2d*w0t   + fi(18) *w2md*w0t &
                 + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt &
                 + fi(21) *w1d*w1t   + fi(22) *w1md*w1t &
                 + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt &
                 + fi(25) *w2d*w1t   + fi(26) *w2md*w1t &
                 + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt &
                 + fi(29) *w1d*w2t   + fi(30) *w1md*w2t &
                 + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt &
                 + fi(33) *w2d*w2t   + fi(34) *w2md*w2t &
                 + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt



      ! cubic hermite polynomial statement functions
      ! psi0 & derivatives
            xpsi0(z)  = z * z * (2.0d0*z - 3.0d0) + 1.0
            xdpsi0(z) = z * (6.0d0*z - 6.0d0)


      ! psi1 & derivatives
            xpsi1(z)  = z * ( z * (z - 2.0d0) + 1.0d0)
            xdpsi1(z) = z * (3.0d0*z - 4.0d0) + 1.0d0


      ! bicubic hermite polynomial statement function
            h3(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) = &
                   fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t &
                 + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt &
                 + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t &
                 + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt &
                 + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t &
                 + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt &
                 + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t &
                 + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt



      ! popular format statements
      01    format(1x,5(a,1pe11.3))
      02    format(1x,a,1p4e16.8)
      03    format(1x,4(a,1pe11.3))
      04    format(1x,4(a,i4))



      ! start of pipeline loop, normal execution starts here
            eosfail = .false.
            do j=jlo_eos,jhi_eos

      !       if (temp_row(j) .le. 0.0) stop 'temp less than 0 in helmeos'
      !       if (den_row(j)  .le. 0.0) stop 'den less than 0 in helmeos'

             temp  = temp_row(j)
             den   = den_row(j)
             abar  = abar_row(j)
             zbar  = zbar_row(j)
             ytot1 = 1.0d0/abar
             ye    = max(1.0d-16,ytot1 * zbar)



      ! initialize
             deni    = 1.0d0/den
             tempi   = 1.0d0/temp
             kt      = kerg * temp
             ktinv   = 1.0d0/kt


      ! radiation section:
             prad    = asoli3 * temp * temp * temp * temp
             dpraddd = 0.0d0
             dpraddt = 4.0d0 * prad*tempi
             dpradda = 0.0d0
             dpraddz = 0.0d0

             erad    = 3.0d0 * prad*deni
             deraddd = -erad*deni
             deraddt = 3.0d0 * dpraddt*deni
             deradda = 0.0d0
             deraddz = 0.0d0

             srad    = (prad*deni + erad)*tempi
             dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
             dsraddt = (dpraddt*deni + deraddt - srad)*tempi
             dsradda = 0.0d0
             dsraddz = 0.0d0


      ! ion section:
              xni     = avo * ytot1 * den
              dxnidd  = avo * ytot1
              dxnida  = -xni * ytot1

              pion    = xni * kt
              dpiondd = dxnidd * kt
              dpiondt = xni * kerg
              dpionda = dxnida * kt
              dpiondz = 0.0d0

              eion    = 1.5d0 * pion*deni
              deiondd = (1.5d0 * dpiondd - eion)*deni
              deiondt = 1.5d0 * dpiondt*deni
              deionda = 1.5d0 * dpionda*deni
              deiondz = 0.0d0


      ! sackur-tetrode equation for the ion entropy of
      ! a single ideal gas characterized by abar
              x       = abar*abar*sqrt(abar) * deni/avo
              s       = sioncon * temp
              z       = x * s * sqrt(s)
              y       = log(z)

      !        y       = 1.0d0/(abar*kt)
      !        yy      = y * sqrt(y)
      !        z       = xni * sifac * yy
      !        etaion  = log(z)


              sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
              dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi &
                         - kergavo * deni * ytot1
              dsiondt = (dpiondt*deni + deiondt)*tempi - &
                        (pion*deni + eion) * tempi*tempi &
                        + 1.5d0 * kergavo * tempi*ytot1
              x       = avo*kerg/abar
              dsionda = (dpionda*deni + deionda)*tempi &
                        + kergavo*ytot1*ytot1* (2.5d0 - y)
              dsiondz = 0.0d0



      ! electron-positron section:


      ! assume complete ionization
              xnem    = xni * zbar


      ! enter the table with ye*den
              din = ye*den


      ! bomb proof the input
              if (temp .gt. t(jmax)) then
               write(6,01) 'temp=',temp,' t(jmax)=',t(jmax)
               write(6,*) 'temp too hot, off grid'
               write(6,*) 'setting eosfail to true and returning'
               eosfail = .true.
               return
              end if
              if (temp .lt. t(1)) then
               write(6,01) 'temp=',temp,' t(1)=',t(1)
               write(6,*) 'temp too cold, off grid'
               write(6,*) 'setting eosfail to true and returning'
               eosfail = .true.
               return
              end if
              if (din  .gt. d(imax)) then
               write(6,01) 'den*ye=',din,' d(imax)=',d(imax)
               write(6,*) 'ye*den too big, off grid'
               write(6,*) 'setting eosfail to true and returning'
               eosfail = .true.
               return
              end if
              if (din  .lt. d(1)) then
               write(6,01) 'ye*den=',din,' d(1)=',d(1)
               write(6,*) 'ye*den too small, off grid'
               write(6,*) 'setting eosfail to true and returning'
               eosfail = .true.
               return
              end if

      ! hash locate this temperature and density
              jat = int((log10(temp) - tlo)*tstpi) + 1
              jat = max(1,min(jat,jmax-1))
              iat = int((log10(din) - dlo)*dstpi) + 1
              iat = max(1,min(iat,imax-1))


      ! access the table locations only once
              fi(1)  = f(iat,jat)
              fi(2)  = f(iat+1,jat)
              fi(3)  = f(iat,jat+1)
              fi(4)  = f(iat+1,jat+1)
              fi(5)  = ft(iat,jat)
              fi(6)  = ft(iat+1,jat)
              fi(7)  = ft(iat,jat+1)
              fi(8)  = ft(iat+1,jat+1)
              fi(9)  = ftt(iat,jat)
              fi(10) = ftt(iat+1,jat)
              fi(11) = ftt(iat,jat+1)
              fi(12) = ftt(iat+1,jat+1)
              fi(13) = fd(iat,jat)
              fi(14) = fd(iat+1,jat)
              fi(15) = fd(iat,jat+1)
              fi(16) = fd(iat+1,jat+1)
              fi(17) = fdd(iat,jat)
              fi(18) = fdd(iat+1,jat)
              fi(19) = fdd(iat,jat+1)
              fi(20) = fdd(iat+1,jat+1)
              fi(21) = fdt(iat,jat)
              fi(22) = fdt(iat+1,jat)
              fi(23) = fdt(iat,jat+1)
              fi(24) = fdt(iat+1,jat+1)
              fi(25) = fddt(iat,jat)
              fi(26) = fddt(iat+1,jat)
              fi(27) = fddt(iat,jat+1)
              fi(28) = fddt(iat+1,jat+1)
              fi(29) = fdtt(iat,jat)
              fi(30) = fdtt(iat+1,jat)
              fi(31) = fdtt(iat,jat+1)
              fi(32) = fdtt(iat+1,jat+1)
              fi(33) = fddtt(iat,jat)
              fi(34) = fddtt(iat+1,jat)
              fi(35) = fddtt(iat,jat+1)
              fi(36) = fddtt(iat+1,jat+1)


      ! various differences
              xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
              xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
              mxt = 1.0d0 - xt
              mxd = 1.0d0 - xd

      ! the six density and six temperature basis functions
              si0t =   psi0(xt)
              si1t =   psi1(xt)*dt_sav(jat)
              si2t =   psi2(xt)*dt2_sav(jat)

              si0mt =  psi0(mxt)
              si1mt = -psi1(mxt)*dt_sav(jat)
              si2mt =  psi2(mxt)*dt2_sav(jat)

              si0d =   psi0(xd)
              si1d =   psi1(xd)*dd_sav(iat)
              si2d =   psi2(xd)*dd2_sav(iat)

              si0md =  psi0(mxd)
              si1md = -psi1(mxd)*dd_sav(iat)
              si2md =  psi2(mxd)*dd2_sav(iat)

      ! derivatives of the weight functions
              dsi0t =   dpsi0(xt)*dti_sav(jat)
              dsi1t =   dpsi1(xt)
              dsi2t =   dpsi2(xt)*dt_sav(jat)

              dsi0mt = -dpsi0(mxt)*dti_sav(jat)
              dsi1mt =  dpsi1(mxt)
              dsi2mt = -dpsi2(mxt)*dt_sav(jat)

              dsi0d =   dpsi0(xd)*ddi_sav(iat)
              dsi1d =   dpsi1(xd)
              dsi2d =   dpsi2(xd)*dd_sav(iat)

              dsi0md = -dpsi0(mxd)*ddi_sav(iat)
              dsi1md =  dpsi1(mxd)
              dsi2md = -dpsi2(mxd)*dd_sav(iat)

      ! second derivatives of the weight functions
              ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
              ddsi1t =   ddpsi1(xt)*dti_sav(jat)
              ddsi2t =   ddpsi2(xt)

              ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
              ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
              ddsi2mt =  ddpsi2(mxt)

      !        ddsi0d =   ddpsi0(xd)*dd2i_sav(iat)
      !        ddsi1d =   ddpsi1(xd)*ddi_sav(iat)
      !        ddsi2d =   ddpsi2(xd)

      !        ddsi0md =  ddpsi0(mxd)*dd2i_sav(iat)
      !        ddsi1md = -ddpsi1(mxd)*ddi_sav(iat)
      !        ddsi2md =  ddpsi2(mxd)


      ! the free energy
              free  = h5(iat,jat, &
                      si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
                      si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

      ! derivative with respect to density
              df_d  = h5(iat,jat, &
                      si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
                      dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)


      ! derivative with respect to temperature
              df_t = h5(iat,jat, &
                      dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
                      si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

      ! derivative with respect to density**2
      !        df_dd = h5(iat,jat,
      !     1          si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
      !     2          ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

      ! derivative with respect to temperature**2
              df_tt = h5(iat,jat, &
                    ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
                      si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

      ! derivative with respect to temperature and density
              df_dt = h5(iat,jat, &
                      dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
                      dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)



      ! now get the pressure derivative with density, chemical potential, and
      ! electron positron number densities
      ! get the interpolation weight functions
              si0t   =  xpsi0(xt)
              si1t   =  xpsi1(xt)*dt_sav(jat)

              si0mt  =  xpsi0(mxt)
              si1mt  =  -xpsi1(mxt)*dt_sav(jat)

              si0d   =  xpsi0(xd)
              si1d   =  xpsi1(xd)*dd_sav(iat)

              si0md  =  xpsi0(mxd)
              si1md  =  -xpsi1(mxd)*dd_sav(iat)


      ! derivatives of weight functions
              dsi0t  = xdpsi0(xt)*dti_sav(jat)
              dsi1t  = xdpsi1(xt)

              dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
              dsi1mt = xdpsi1(mxt)

              dsi0d  = xdpsi0(xd)*ddi_sav(iat)
              dsi1d  = xdpsi1(xd)

              dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
              dsi1md = xdpsi1(mxd)


      ! look in the pressure derivative only once
              fi(1)  = dpdf(iat,jat)
              fi(2)  = dpdf(iat+1,jat)
              fi(3)  = dpdf(iat,jat+1)
              fi(4)  = dpdf(iat+1,jat+1)
              fi(5)  = dpdft(iat,jat)
              fi(6)  = dpdft(iat+1,jat)
              fi(7)  = dpdft(iat,jat+1)
              fi(8)  = dpdft(iat+1,jat+1)
              fi(9)  = dpdfd(iat,jat)
              fi(10) = dpdfd(iat+1,jat)
              fi(11) = dpdfd(iat,jat+1)
              fi(12) = dpdfd(iat+1,jat+1)
              fi(13) = dpdfdt(iat,jat)
              fi(14) = dpdfdt(iat+1,jat)
              fi(15) = dpdfdt(iat,jat+1)
              fi(16) = dpdfdt(iat+1,jat+1)

      ! pressure derivative with density
              dpepdd  = h3(iat,jat, &
                             si0t,   si1t,   si0mt,   si1mt, &
                             si0d,   si1d,   si0md,   si1md)
              dpepdd  = max(ye * dpepdd,1.0d-30)



      ! look in the electron chemical potential table only once
              fi(1)  = ef(iat,jat)
              fi(2)  = ef(iat+1,jat)
              fi(3)  = ef(iat,jat+1)
              fi(4)  = ef(iat+1,jat+1)
              fi(5)  = eft(iat,jat)
              fi(6)  = eft(iat+1,jat)
              fi(7)  = eft(iat,jat+1)
              fi(8)  = eft(iat+1,jat+1)
              fi(9)  = efd(iat,jat)
              fi(10) = efd(iat+1,jat)
              fi(11) = efd(iat,jat+1)
              fi(12) = efd(iat+1,jat+1)
              fi(13) = efdt(iat,jat)
              fi(14) = efdt(iat+1,jat)
              fi(15) = efdt(iat,jat+1)
              fi(16) = efdt(iat+1,jat+1)


      ! electron chemical potential etaele
              etaele  = h3(iat,jat, &
                           si0t,   si1t,   si0mt,   si1mt, &
                           si0d,   si1d,   si0md,   si1md)


      ! derivative with respect to density
              x       = h3(iat,jat, &
                           si0t,   si1t,   si0mt,   si1mt, &
                          dsi0d,  dsi1d,  dsi0md,  dsi1md)
              detadd  = ye * x

      ! derivative with respect to temperature
              detadt  = h3(iat,jat, &
                          dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
                           si0d,   si1d,   si0md,   si1md)

      ! derivative with respect to abar and zbar
             detada = -x * din * ytot1
             detadz =  x * den * ytot1



      ! look in the number density table only once
              fi(1)  = xf(iat,jat)
              fi(2)  = xf(iat+1,jat)
              fi(3)  = xf(iat,jat+1)
              fi(4)  = xf(iat+1,jat+1)
              fi(5)  = xft(iat,jat)
              fi(6)  = xft(iat+1,jat)
              fi(7)  = xft(iat,jat+1)
              fi(8)  = xft(iat+1,jat+1)
              fi(9)  = xfd(iat,jat)
              fi(10) = xfd(iat+1,jat)
              fi(11) = xfd(iat,jat+1)
              fi(12) = xfd(iat+1,jat+1)
              fi(13) = xfdt(iat,jat)
              fi(14) = xfdt(iat+1,jat)
              fi(15) = xfdt(iat,jat+1)
              fi(16) = xfdt(iat+1,jat+1)

      ! electron + positron number densities
             xnefer   = h3(iat,jat, &
                           si0t,   si1t,   si0mt,   si1mt, &
                           si0d,   si1d,   si0md,   si1md)

      ! derivative with respect to density
             x        = h3(iat,jat, &
                           si0t,   si1t,   si0mt,   si1mt, &
                          dsi0d,  dsi1d,  dsi0md,  dsi1md)
             x = max(x,1.0d-30)
             dxnedd   = ye * x

      ! derivative with respect to temperature
             dxnedt   = h3(iat,jat, &
                          dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
                           si0d,   si1d,   si0md,   si1md)

      ! derivative with respect to abar and zbar
             dxneda = -x * din * ytot1
             dxnedz =  x  * den * ytot1


      ! the desired electron-positron thermodynamic quantities

      ! dpepdd at high temperatures and low densities is below the
      ! floating point limit of the subtraction of two large terms.
      ! since dpresdd doesn't enter the maxwell relations at all, use the
      ! bicubic interpolation done above instead of the formally correct expression
              x       = din * din
              pele    = x * df_d
              dpepdt  = x * df_dt
      !        dpepdd  = ye * (x * df_dd + 2.0d0 * din * df_d)
              s       = dpepdd/ye - 2.0d0 * din * df_d
              dpepda  = -ytot1 * (2.0d0 * pele + s * din)
              dpepdz  = den*ytot1*(2.0d0 * din * df_d  +  s)


              x       = ye * ye
              sele    = -df_t * ye
              dsepdt  = -df_tt * ye
              dsepdd  = -df_dt * x
              dsepda  = ytot1 * (ye * df_dt * din - sele)
              dsepdz  = -ytot1 * (ye * df_dt * den  + df_t)


              eele    = ye*free + temp * sele
              deepdt  = temp * dsepdt
              deepdd  = x * df_d + temp * dsepdd
              deepda  = -ye * ytot1 * (free +  df_d * din) + temp * dsepda
              deepdz  = ytot1* (free + ye * df_d * den) + temp * dsepdz




      ! coulomb section:

      ! uniform background corrections only
      ! from yakovlev & shalybkov 1989
      ! lami is the average ion seperation
      ! plasg is the plasma coupling parameter

              z        = forth * pi
              s        = z * xni
              dsdd     = z * dxnidd
              dsda     = z * dxnida

              lami     = 1.0d0/s**third
              inv_lami = 1.0d0/lami
              z        = -third * lami
              lamidd   = z * dsdd/s
              lamida   = z * dsda/s

              plasg    = zbar*zbar*esqu*ktinv*inv_lami
              z        = -plasg * inv_lami
              plasgdd  = z * lamidd
              plasgda  = z * lamida
              plasgdt  = -plasg*ktinv * kerg
              plasgdz  = 2.0d0 * plasg/zbar


      ! yakovlev & shalybkov 1989 equations 82, 85, 86, 87
              if (plasg .ge. 1.0) then
               x        = plasg**(0.25d0)
               y        = avo * ytot1 * kerg
               ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
               pcoul    = third * den * ecoul
               scoul    = -y * (3.0d0*b1*x - 5.0d0*c1/x &
                          + d1 * (log(plasg) - 1.0d0) - e1)

               y        = avo*ytot1*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
               decouldd = y * plasgdd
               decouldt = y * plasgdt + ecoul/temp
               decoulda = y * plasgda - ecoul/abar
               decouldz = y * plasgdz

               y        = third * den
               dpcouldd = third * ecoul + y*decouldd
               dpcouldt = y * decouldt
               dpcoulda = y * decoulda
               dpcouldz = y * decouldz


               y        = -avo*kerg/(abar*plasg)*(0.75d0*b1*x+1.25d0*c1/x+d1)
               dscouldd = y * plasgdd
               dscouldt = y * plasgdt
               dscoulda = y * plasgda - scoul/abar
               dscouldz = y * plasgdz


      ! yakovlev & shalybkov 1989 equations 102, 103, 104
              else if (plasg .lt. 1.0) then
               x        = plasg*sqrt(plasg)
               y        = plasg**b2
               z        = c2 * x - third * a2 * y
               pcoul    = -pion * z
               ecoul    = 3.0d0 * pcoul/den
               scoul    = -avo/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)

               s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
               dpcouldd = -dpiondd*z - pion*s*plasgdd
               dpcouldt = -dpiondt*z - pion*s*plasgdt
               dpcoulda = -dpionda*z - pion*s*plasgda
               dpcouldz = -dpiondz*z - pion*s*plasgdz

               s        = 3.0d0/den
               decouldd = s * dpcouldd - ecoul/den
               decouldt = s * dpcouldt
               decoulda = s * dpcoulda
               decouldz = s * dpcouldz

               s        = -avo*kerg/(abar*plasg)*(1.5d0*c2*x-a2*(b2-1.0d0)*y)
               dscouldd = s * plasgdd
               dscouldt = s * plasgdt
               dscoulda = s * plasgda - scoul/abar
               dscouldz = s * plasgdz
              end if


      ! bomb proof
              x   = prad + pion + pele + pcoul
              y   = erad + eion + eele + ecoul
              z   = srad + sion + sele + scoul

      !        write(6,*) x,y,z
      !        if (x .le. 0.0 .or. y .le. 0.0 .or. z .le. 0.0) then
              if (x .le. 0.0 .or. y .le. 0.0) then
      !        if (x .le. 0.0) then

      !         write(6,*)
      !         write(6,*) 'coulomb corrections are causing a negative pressure'
      !         write(6,*) 'setting all coulomb corrections to zero'
      !         write(6,*)

               pcoul    = 0.0d0
               dpcouldd = 0.0d0
               dpcouldt = 0.0d0
               dpcoulda = 0.0d0
               dpcouldz = 0.0d0
               ecoul    = 0.0d0
               decouldd = 0.0d0
               decouldt = 0.0d0
               decoulda = 0.0d0
               decouldz = 0.0d0
               scoul    = 0.0d0
               dscouldd = 0.0d0
               dscouldt = 0.0d0
               dscoulda = 0.0d0
               dscouldz = 0.0d0
              end if


      ! sum all the gas components
             pgas    = pion + pele + pcoul
             egas    = eion + eele + ecoul
             sgas    = sion + sele + scoul

             dpgasdd = dpiondd + dpepdd + dpcouldd
             dpgasdt = dpiondt + dpepdt + dpcouldt
             dpgasda = dpionda + dpepda + dpcoulda
             dpgasdz = dpiondz + dpepdz + dpcouldz

             degasdd = deiondd + deepdd + decouldd
             degasdt = deiondt + deepdt + decouldt
             degasda = deionda + deepda + decoulda
             degasdz = deiondz + deepdz + decouldz

             dsgasdd = dsiondd + dsepdd + dscouldd
             dsgasdt = dsiondt + dsepdt + dscouldt
             dsgasda = dsionda + dsepda + dscoulda
             dsgasdz = dsiondz + dsepdz + dscouldz




      ! add in radiation to get the total
             pres    = prad + pgas
             ener    = erad + egas
             entr    = srad + sgas

             dpresdd = dpraddd + dpgasdd
             dpresdt = dpraddt + dpgasdt
             dpresda = dpradda + dpgasda
             dpresdz = dpraddz + dpgasdz

             denerdd = deraddd + degasdd
             denerdt = deraddt + degasdt
             denerda = deradda + degasda
             denerdz = deraddz + degasdz

             dentrdd = dsraddd + dsgasdd
             dentrdt = dsraddt + dsgasdt
             dentrda = dsradda + dsgasda
             dentrdz = dsraddz + dsgasdz


      ! for the gas
      ! the temperature and density exponents (c&g 9.81 9.82)
      ! the specific heat at constant volume (c&g 9.92)
      ! the third adiabatic exponent (c&g 9.93)
      ! the first adiabatic exponent (c&g 9.97)
      ! the second adiabatic exponent (c&g 9.105)
      ! the specific heat at constant pressure (c&g 9.98)
      ! and relativistic formula for the sound speed (c&g 14.29)

             zz        = pgas*deni
             zzi       = den/pgas
             chit_gas  = temp/pgas * dpgasdt
             chid_gas  = dpgasdd*zzi
             cv_gas    = degasdt
             x         = zz * chit_gas/(temp * cv_gas)
             gam3_gas  = x + 1.0d0
             gam1_gas  = chit_gas*x + chid_gas
             nabad_gas = x/gam1_gas
             gam2_gas  = 1.0d0/(1.0d0 - nabad_gas)
             cp_gas    = cv_gas * gam1_gas/chid_gas
             z         = 1.0d0 + (egas + light2)*zzi
             sound_gas = clight * sqrt(gam1_gas/z)



      ! for the totals
             zz    = pres*deni
             zzi   = den/pres
             chit  = temp/pres * dpresdt
             chid  = dpresdd*zzi
             cv    = denerdt
             x     = zz * chit/(temp * cv)
             gam3  = x + 1.0d0
             gam1  = chit*x + chid

             nabad = x/gam1
             gam2  = 1.0d0/(1.0d0 - nabad)
             cp    = cv * gam1/chid
             z     = 1.0d0 + (ener + light2)*zzi
             sound = clight * sqrt(gam1/z)



      ! maxwell relations; each is zero if the consistency is perfect
             x   = den * den

             dse = temp*dentrdt/denerdt - 1.0d0

             dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0

             dsp = -dentrdd*x/dpresdt - 1.0d0


      ! store this row
              ptot_row(j)   = pres
              dpt_row(j)    = dpresdt
              dpd_row(j)    = dpresdd
              dpa_row(j)    = dpresda
              dpz_row(j)    = dpresdz

              etot_row(j)   = ener
              det_row(j)    = denerdt
              ded_row(j)    = denerdd
              dea_row(j)    = denerda
              dez_row(j)    = denerdz

              stot_row(j)   = entr
              dst_row(j)    = dentrdt
              dsd_row(j)    = dentrdd
              dsa_row(j)    = dentrda
              dsz_row(j)    = dentrdz


              pgas_row(j)   = pgas
              dpgast_row(j) = dpgasdt
              dpgasd_row(j) = dpgasdd
              dpgasa_row(j) = dpgasda
              dpgasz_row(j) = dpgasdz

              egas_row(j)   = egas
              degast_row(j) = degasdt
              degasd_row(j) = degasdd
              degasa_row(j) = degasda
              degasz_row(j) = degasdz

              sgas_row(j)   = sgas
              dsgast_row(j) = dsgasdt
              dsgasd_row(j) = dsgasdd
              dsgasa_row(j) = dsgasda
              dsgasz_row(j) = dsgasdz


              prad_row(j)   = prad
              dpradt_row(j) = dpraddt
              dpradd_row(j) = dpraddd
              dprada_row(j) = dpradda
              dpradz_row(j) = dpraddz

              erad_row(j)   = erad
              deradt_row(j) = deraddt
              deradd_row(j) = deraddd
              derada_row(j) = deradda
              deradz_row(j) = deraddz

              srad_row(j)   = srad
              dsradt_row(j) = dsraddt
              dsradd_row(j) = dsraddd
              dsrada_row(j) = dsradda
              dsradz_row(j) = dsraddz


              pion_row(j)   = pion
              dpiont_row(j) = dpiondt
              dpiond_row(j) = dpiondd
              dpiona_row(j) = dpionda
              dpionz_row(j) = dpiondz

              eion_row(j)   = eion
              deiont_row(j) = deiondt
              deiond_row(j) = deiondd
              deiona_row(j) = deionda
              deionz_row(j) = deiondz

              sion_row(j)   = sion
              dsiont_row(j) = dsiondt
              dsiond_row(j) = dsiondd
              dsiona_row(j) = dsionda
              dsionz_row(j) = dsiondz

              xni_row(j)    = xni

              pele_row(j)   = pele
              ppos_row(j)   = 0.0d0
              dpept_row(j)  = dpepdt
              dpepd_row(j)  = dpepdd
              dpepa_row(j)  = dpepda
              dpepz_row(j)  = dpepdz

              eele_row(j)   = eele
              epos_row(j)   = 0.0d0
              deept_row(j)  = deepdt
              deepd_row(j)  = deepdd
              deepa_row(j)  = deepda
              deepz_row(j)  = deepdz

              sele_row(j)   = sele
              spos_row(j)   = 0.0d0
              dsept_row(j)  = dsepdt
              dsepd_row(j)  = dsepdd
              dsepa_row(j)  = dsepda
              dsepz_row(j)  = dsepdz

              xnem_row(j)   = xnem
              xne_row(j)    = xnefer
              dxnet_row(j)  = dxnedt
              dxned_row(j)  = dxnedd
              dxnea_row(j)  = dxneda
              dxnez_row(j)  = dxnedz
              xnp_row(j)    = 0.0d0
              zeff_row(j)   = zbar

              etaele_row(j) = etaele
              detat_row(j)  = detadt
              detad_row(j)  = detadd
              detaa_row(j)  = detada
              detaz_row(j)  = detadz
              etapos_row(j) = 0.0d0

              pcou_row(j)   = pcoul
              dpcout_row(j) = dpcouldt
              dpcoud_row(j) = dpcouldd
              dpcoua_row(j) = dpcoulda
              dpcouz_row(j) = dpcouldz

              ecou_row(j)   = ecoul
              decout_row(j) = decouldt
              decoud_row(j) = decouldd
              decoua_row(j) = decoulda
              decouz_row(j) = decouldz

              scou_row(j)   = scoul
              dscout_row(j) = dscouldt
              dscoud_row(j) = dscouldd
              dscoua_row(j) = dscoulda
              dscouz_row(j) = dscouldz

              plasg_row(j)  = plasg

              dse_row(j)    = dse
              dpe_row(j)    = dpe
              dsp_row(j)    = dsp

              cv_gas_row(j)    = cv_gas
              cp_gas_row(j)    = cp_gas
              gam1_gas_row(j)  = gam1_gas
              gam2_gas_row(j)  = gam2_gas
              gam3_gas_row(j)  = gam3_gas
              nabad_gas_row(j) = nabad_gas
              cs_gas_row(j)    = sound_gas

              cv_row(j)     = cv
              cp_row(j)     = cp
              gam1_row(j)   = gam1
              gam2_row(j)   = gam2
              gam3_row(j)   = gam3
              nabad_row(j)  = nabad
              cs_row(j)     = sound

      ! end of pipeline loop
            enddo
            return
            end

            ! routine read_helm_table reads an electron helm free energy table
            ! routine helmeos computes the pressure, energy and entropy via tables


  !======================= TIMMES EOS TABLE READING ===========================


                  subroutine read_helm_table
                  ! include '/home/seb/Recherche/GENEC/GENEC_Timmes_LESTA/GENEC_Timmes/code/Timmes_EOS/implno.dek'
                  ! include '/home/seb/Recherche/GENEC/GENEC_Timmes_LESTA/GENEC_Timmes/code/Timmes_EOS/helm_table_storage.dek'
                  use evol, only: input_dir
                  include 'Timmes_EOS/implno.dek'
                  include 'Timmes_EOS/helm_table_storage.dek'

            ! this routine reads the helmholtz eos file, and
            ! must be called once before the helmeos routine is invoked.

            ! declare local variables
                  integer          i,j
                  double precision tsav,dsav,dth,dt2,dti,dt2i,dt3i, &
                                   dd,dd2,ddi,dd2i,dd3i


            ! open the file (use softlinks to input the desired table)

                   open(unit=19,file=trim(input_dir)//trim('Timmes_EOS/helm_table.dat'),&
                   status='old')


            ! for standard table limits
                   tlo   = 3.0d0
                   thi   = 13.0d0
                   tstp  = (thi - tlo)/float(jmax-1)
                   tstpi = 1.0d0/tstp
                   dlo   = -12.0d0
                   dhi   = 15.0d0
                   dstp  = (dhi - dlo)/float(imax-1)
                   dstpi = 1.0d0/dstp

            ! read the helmholtz free energy and its derivatives
                   do j=1,jmax
                    tsav = tlo + (j-1)*tstp
                    t(j) = 10.0d0**(tsav)
                    do i=1,imax
                     dsav = dlo + (i-1)*dstp
                     d(i) = 10.0d0**(dsav)
                     read(19,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j), &
                              fddt(i,j),fdtt(i,j),fddtt(i,j)
                    enddo
                   enddo
            !       write(6,*) 'read main table'


            ! read the pressure derivative with density table
                   do j=1,jmax
                    do i=1,imax
                     read(19,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
                    enddo
                   enddo
            !       write(6,*) 'read dpdd table'

            ! read the electron chemical potential table
                   do j=1,jmax
                    do i=1,imax
                     read(19,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
                    enddo
                   enddo
            !       write(6,*) 'read eta table'

            ! read the number density table
                   do j=1,jmax
                    do i=1,imax
                     read(19,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
                    enddo
                   enddo
            !       write(6,*) 'read xne table'

            ! close the file
                  close(unit=19)


            ! construct the temperature and density deltas and their inverses
                   do j=1,jmax-1
                    dth          = t(j+1) - t(j)
                    dt2         = dth * dth
                    dti         = 1.0d0/dth
                    dt2i        = 1.0d0/dt2
                    dt3i        = dt2i*dti
                    dt_sav(j)   = dth
                    dt2_sav(j)  = dt2
                    dti_sav(j)  = dti
                    dt2i_sav(j) = dt2i
                    dt3i_sav(j) = dt3i
                   end do
                   do i=1,imax-1
                    dd          = d(i+1) - d(i)
                    dd2         = dd * dd
                    ddi         = 1.0d0/dd
                    dd2i        = 1.0d0/dd2
                    dd3i        = dd2i*ddi
                    dd_sav(i)   = dd
                    dd2_sav(i)  = dd2
                    ddi_sav(i)  = ddi
                    dd2i_sav(i) = dd2i
                    dd3i_sav(i) = dd3i
                   enddo



            !      write(6,*)
            !      write(6,*) 'finished reading eos table'
            !      write(6,04) 'imax=',imax,' jmax=',jmax
            !04    format(1x,4(a,i4))
            !      write(6,03) 'temp(1)   =',t(1),' temp(jmax)   =',t(jmax)
            !      write(6,03) 'ye*den(1) =',d(1),' ye*den(imax) =',d(imax)
            !03    format(1x,4(a,1pe11.3))
            !      write(6,*)

                  return
                  end

!========================== END EOS TIMMES ==============================





  !========================================================================
                            ! Dichte EOS !
                               ! 2013 !

                               !EOS==0!

  ! GENEC Default equation of state.

  ! Contains: degen (for degeneracy computations)
  !           Compute_sepp
  !           dichte (EOS computation)

  !======================================================================

  SUBROUTINE degen
    !-----------------------------------------------------------------------
    USE const,ONLY: pi,cst_h,cst_mh

    IMPLICIT NONE

    INTEGER:: i
    REAL(kindreal):: hchi,psi,sr,srs,srt=0.d0,sp,sps,spt=0.d0,su=0.d0,sus=0.d0,sut=0.d0,tu,wz,asr,asp,asu,xsq,x,wx,gol,phi1,ph1, &
         phi2,ph2,phi1d,ph1d,phi2d,ph2d,phi3,ph3,phi3d,ph3d

    REAL(kindreal)::cst_mch3,ccd,ccr,ccp,ccu
    ! Les data suivants (dega,degu,degex) sont tires de Kippenhahn, Weigert & Hofmeister 1967
    ! Methods in Computational Physics vol. 7, p. 178:
    REAL(kindreal), DIMENSION(12):: &
         ! a_i:
         dega=(/7.28905d-2,2.04648d-1,2.54685d-1,1.96668d-1,1.04294d-1,3.95716d-2,1.09024d-2,&
         2.21504d-3,3.14295d-4,3.64498d-5,2.42533d-6,1.49830d-7/), &
         ! u_i:
         degu=(/0.11896d0,0.47652d0,1.07476d0,1.91722d0,3.00898d0,4.35703d0,5.96967d0,7.86147d0,&
         10.02706d0,12.57936d0,15.19876d0,18.87814d0/), &
         ! alpha_i:
         degex=(/8.87844d-1,6.20940d-1,3.41379d-1,1.47015d-1,4.93418d-2,1.28164d-2,2.55507d-3,&
         3.85307d-4,4.41879d-5,3.44234d-6, 2.50763d-7,6.32888d-9/)

    !-----------------------------------------------------------------------
    cst_mch3 = (cst_me*cst_c/cst_h)**3.d0
    ccd = pi**2.d0 / 6.d0
    ccr = 8.d0*pi*cst_mh*cst_mch3
    ccp = (16.d0*pi/6.d0)*cst_me*cst_c**2.d0*cst_mch3
    ccu = cst_me*cst_c**2.d0/cst_mh

    IF (hpsi >= 0.d0) THEN
       hchi=1.d0
       psi=hpsi
       xsq=2.d0*tk*psi+(tk*psi)**2.d0
       x= SQRT(xsq)
       wx=1.d0+tk*psi
       gol=LOG(wx-x)
       phi1=xsq*x/3.d0
       ph1=tk*tk*(1.d0+2.d0*xsq)/x
       IF ((x-0.5d0) > 0.d0) THEN
          phi2=x*wx*(2.d0*xsq-3.d0)*0.125d0-0.375d0*gol
       ELSE
          phi2=0.125d0*x*xsq*xsq*(8.d0/5.d0-xsq*(4.d0/7.d0-xsq*(1.d0/3.d0-xsq*5.d0/22.d0)))
       ENDIF
       ph2=tk**2.d0*3.d0*x*wx
       phi1d=x*wx
       ph1d=(2.d0*xsq-1.d0)*tk*tk*wx/(xsq*x)
       phi2d=x*xsq
       ph2d=(6.d0*xsq+3.d0)*tk*tk/x
       sr=phi1+ccd*ph1
       srs=-(phi1d+ccd*ph1d)*tk/hchi
       sp=phi2+ccd*ph2
       sps=-(phi2d+ccd*ph2d)*tk/hchi
       IF ((hpsi+1.d0) >= 0.d0) THEN
          IF ((x-0.5d0) > 0.d0) THEN
             phi3=x*wx*(2.d0*xsq+1.d0)*0.125d0-phi1+0.125d0*gol
          ELSE
             phi3=x*xsq*xsq*(0.1d0-xsq*(1.d0/56.d0-xsq*(1.d0/144.d0-xsq*5.d0/1408.d0)))
          ENDIF
          ph3=tk*tk*(wx*(1.d0+3.d0*xsq)-1.d0-2.d0*xsq)/x
          phi3d=x*wx*(wx-1.d0)
          ph3d=tk*tk*(xsq*(6.d0*xsq+3.d0)-1.d0-wx*(2.d0*xsq-1.d0))/(xsq*x)
          srt=(phi1d+ccd*ph1d)*psi+2.d0*ccd/tk*ph1
          spt=(phi2d+ccd*ph2d)*psi+2.d0*ccd/tk*ph2
          su=phi3+ccd*ph3
          sus=-(phi3d+ccd*ph3d)*tk/hchi
          sut=(phi3d+ccd*ph3d)*psi+2.d0*ccd/tk*ph3
       ENDIF
    ELSE
       hchi=chi
       IF ((hchi-1.d0/EXP(7.d0)) >= 0.d0) THEN
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
          DO i=1,12
             tu=tk*degu(i)
             wz=SQRT (tk*(2.d0+tu))
             asr=(1.d0+tu)*tk*wz*dega(i)/(hchi+degex(i))
             asp=tu*wz*wz*wz*dega(i)/(hchi+degex(i))
             sr=sr+asr
             srs=srs-asr/(hchi+degex(i))
             sp=sp+asp
             sps=sps-asp/(hchi+degex(i))
             IF ((hpsi+1.d0) < 0.d0) CYCLE
             srt=srt+(3.d0+7.d0*tu+3.d0*tu*tu)*tk/wz*dega(i)/(hchi+degex(i))
             spt=spt+(5.d0+4.d0*tu)*tu*wz*dega(i)/(hchi+degex(i))
             asu=(1.d0+tu)*tu*tk*wz*dega(i)/(hchi+degex(i))
             su=su+asu
             sus=sus-asu/(hchi+degex(i))
             sut=sut+(5.d0+10.d0*tu+4.d0*tu*tu)*tu*tk/wz*dega(i)/(hchi+degex(i))
          ENDDO
       ELSE
          psi=-LOG(hchi)
          xsq=2.d0*tk*psi+(tk*psi)**2.d0
          x= SQRT(xsq)
          wx=1.d0+tk*psi
          gol=LOG(wx-x)
          phi1=xsq*x/3.d0
          ph1=tk*tk*(1.d0+2.d0*xsq)/x
          IF ((x-0.5d0) > 0.d0) THEN
             phi2=x*wx*(2.d0*xsq-3.d0)*0.125d0-0.375d0*gol
          ELSE
             phi2=0.125d0*x*xsq*xsq*(8.d0/5.d0-xsq*(4.d0/7.d0-xsq*(1.d0/3.d0-xsq*5.d0/22.d0)))
          ENDIF
          ph2=tk**2.d0*3.d0*x*wx
          phi1d=x*wx
          ph1d=(2.d0*xsq-1.d0)*tk*tk*wx/(xsq*x)
          phi2d=x*xsq
          ph2d=(6.d0*xsq+3.d0)*tk*tk/x
          sr=phi1+ccd*ph1
          srs=-(phi1d+ccd*ph1d)*tk/hchi
          sp=phi2+ccd*ph2
          sps=-(phi2d+ccd*ph2d)*tk/hchi
          IF ((hpsi+1.d0) >= 0.d0) THEN
             IF ((x-0.5d0) > 0.d0) THEN
                phi3=x*wx*(2.d0*xsq+1.d0)*0.125d0-phi1+0.125d0*gol
             ELSE
                phi3=x*xsq*xsq*(0.1d0-xsq*(1.d0/56.d0-xsq*(1.d0/144.d0-xsq*5.d0/1408.d0)))
             ENDIF
             ph3=tk*tk*(wx*(1.d0+3.d0*xsq)-1.d0-2.d0*xsq)/x
             phi3d=x*wx*(wx-1.d0)
             ph3d=tk*tk*(xsq*(6.d0*xsq+3.d0)-1.d0-wx*(2.d0*xsq-1.d0))/(xsq*x)
             srt=(phi1d+ccd*ph1d)*psi+2.d0*ccd/tk*ph1
             spt=(phi2d+ccd*ph2d)*psi+2.d0*ccd/tk*ph2
             su=phi3+ccd*ph3
             sus=-(phi3d+ccd*ph3d)*tk/hchi
             sut=(phi3d+ccd*ph3d)*psi+2.d0*ccd/tk*ph3
          ENDIF
       ENDIF
    ENDIF
    rhe=ccr*sr
    rhes=srs/sr
    !write(*,*),'DICHTE: ','rhe=',rhe,'rhes=',rhes
    pe=ccp*sp
    pes=sps/sp
    IF ((hpsi+1.d0) >= 0.d0) THEN
       rhete=tk*srt/sr
       pete=tk*spt/sp
       ue=ccu*su/sr
       ues=(sr*sus-su*srs)/(sr*su)
       uete=tk*(sr*sut-su*srt)/(sr*su)
    ENDIF

    RETURN

  END SUBROUTINE degen
  !======================================================================
  SUBROUTINE Compute_sepp(krueck,sepp,seppl,dchi)
    !-----------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER,INTENT(in):: krueck
    REAL(kindreal),INTENT(out):: sepp,seppl,dchi
    !-----------------------------------------------------------------------
    CALL degen

    SELECT CASE (krueck)
    CASE (1)
       sepp=(pe+rgaz*vermy*toni*rhe)/pg
       seppl=(pe/pg)*chi*pes+(rgaz*vermy*toni*rhe)/pg*chi*rhes
       dchi=((1.d0-sepp)/seppl)*chi
    CASE (2)
       sepp=(pe+rgaz*vermy*toni*rhe)/pg
       seppl=(pe/pg)*pes
       dchi =sepp*LOG(sepp)/seppl
    CASE default
       STOP 'Bad value for krueck in dichte'
    END SELECT

    RETURN

  END SUBROUTINE Compute_sepp
  !======================================================================
  SUBROUTINE dichte
    !-----------------------------------------------------------------------
    USE const,ONLY: cst_a,um,cst_k
    USE inputparam,ONLY: ialflu,z
    USE abundmod,ONLY: x,y3,y,xc12,xc13,xc14,xn14,xn15,xo16,xo17,xo18,xf18,xf19,xne20,xne21,xne22,xna23,xmg24,xmg25,xmg26, &
         xal26,xal27,xsi28,zabelx,nbelx,abelx,nbzel,nbael
    USE strucmod,ONLY: j1,vmy1,vmyo,vmye,j,beta1,t,p,adi1,vmol,vny,x_env,vmion,vna,vmionp,vmiont
    USE ionisation,ONLY: ionpart

    IMPLICIT NONE

    INTEGER:: ii,mist,krueck
    REAL(kindreal), SAVE:: vx3
    REAL(kindreal):: refad,tion,pion,xcmp,sepp,seppl,dchi,rhes1,rpsi,hifa,pas,pate,gamma_GENEC,gamma_gas
    !-----------------------------------------------------------------------
    !write(*,*)j,'dichte start: T=',t(j),'P=',p(j)

    !-----------------------------------------------------------------------
    IF (j1 == 1) THEN
       num=0
       vx3=0.d0
       chi=55.d0
       hpsi=-3.d0
    ENDIF

    vmy1 = 2.d0*x(j1)+y3(j1)+3.d0/4.d0*y(j1)+7.d0/12.d0*xc12(j1)+7.d0/13.d0*xc13(j1)+8.d0/14.d0*xn14(j1)+8.d0/15.d0*xn15(j1)+ &
         9.d0/16.d0*xo16(j1)+9.d0/17.d0*xo17(j1)+9.d0/18.d0*xo18(j1)+11.d0/20.d0*xne20(j1)+11.d0/22.d0*xne22(j1)+ &
         13.d0/24.d0*xmg24(j1)+13.d0/25.d0*xmg25(j1)+13.d0/26.d0*xmg26(j1)

    vmyo = x(j1)+y3(j1)/3.d0+y(j1)/4.d0+xc12(j1)/12.d0+xc13(j1)/13.d0+xn14(j1)/14.d0+xn15(j1)/15.d0+xo16(j1)/16.d0+ &
         xo17(j1)/17.d0+xo18(j1)/18.d0+xne20(j1)/20.d0+xne22(j1)/22.d0+xmg24(j1)/24.d0+xmg25(j1)/25.d0+xmg26(j1)/26.d0

    vmye=0.5d0*(1.d0+x(j1))+y3(j1)/6.d0

    IF (ialflu == 1) THEN
       vmy1 = vmy1+0.5d0*xc14(j1)+10.d0/18.d0*xf18(j1)+10.d0/19.d0*xf19(j1)+11.d0/21.d0*xne21(j1)+12.d0/23.d0*xna23(j1)+ &
            14.d0/26.d0*xal26(j1)+14.d0/27.d0*xal27(j1)+15.d0/28.d0*xsi28(j1)
       vmyo = vmyo+xc14(j1)/14.d0+xf18(j1)/18.d0+xf19(j1)/19.d0+xne21(j1)/21.d0+xna23(j1)/23.d0+xal26(j1)/26.d0+ &
            xal27(j1)/27.d0+xsi28(j1)/28.d0
       vmye = x(j1)+2.d0/3.d0*y3(j1)+0.5d0*(y(j1)+xc12(j1)+xn14(j1)+xo16(j1)+xne20(j1)+xmg24(j1)+xal26(j)+xsi28(j)+xf18(j))+ &
            6.d0/13.d0*xc13(j1)+7.d0/15.d0*xn15(j1)+8.d0/17.d0*xo17(j1)+8.d0/18.d0*xo18(j1)+10.d0/22.d0*xne22(j1)+ &
            12.d0/25.d0*xmg25(j1)+12.d0/26.d0*xmg26(j1)+6.d0/14.d0*xc14(j1)+9.d0/19.d0*xf19(j1)+10.d0/21.d0*xne21(j1)+ &
            11.d0/23.d0*xna23(j1)+13.d0/27.d0*xal27(j1)
    ENDIF

    !  z elements taken ~ as Ca56: A=56 but (nbzel(ii)+1.)/nbael(ii)-->0.5
    vmy1= vmy1+0.5d0*zabelx
    vmyo= vmyo+zabelx/56.d0
    IF (ialflu == 1) THEN
       vmye = x(j1)+2.d0/3.d0*y3(j1)+0.5d0*(y(j1)+xc12(j1)+xn14(j1)+xo16(j1)+xne20(j1)+xmg24(j1)+xal26(j)+xsi28(j)+ &
            xf18(j))+6.d0/13.d0*xc13(j1)+6.d0/14.d0*xc14(j1)+7.d0/15.d0*xn15(j1)+8.d0/17.d0*xo17(j1)+8.d0/18.d0*xo18(j1)+ &
            9.d0/19.d0*xf19(j1)+10.d0/21.d0*xne21(j1)+10.d0/22.d0*xne22(j1)+11.d0/23.d0*xna23(j1)+12.d0/25.d0*xmg25(j1)+ &
            12.d0/26.d0*xmg26(j1)+13.d0/27.d0*xal27(j1)+0.5d0*zabelx
    ELSE
       vmye = x(j1)+2.d0/3.d0*y3(j1)+0.5d0*(y(j1)+xc12(j1)+xn14(j1)+xo16(j1)+xne20(j1)+xmg24(j1))+6.d0/13.d0*xc13(j1)+ &
            7.d0/15.d0*xn15(j1)+8.d0/17.d0*xo17(j1)+8.d0/18.d0*xo18(j1)+10.d0/22.d0*xne22(j1)+12.d0/25.d0*xmg25(j1)+ &
            12.d0/26.d0*xmg26(j1)+0.5d0*zabelx
    ENDIF

    DO ii=1,nbelx
       vmy1= vmy1+abelx(ii,j1)*REAL(nbzel(ii)+1)/REAL(nbael(ii))
       vmyo= vmyo+abelx(ii,j1)/REAL(nbael(ii))
       vmye= vmye+abelx(ii,j1)*REAL(nbzel(ii))/REAL(nbael(ii))
    ENDDO

    vmy1=-LOG(vmy1)
    vmyo=1.d0/vmyo
    vmye=1.d0/vmye
    vermy=vmye/vmyo

    beta1=1.d0- EXP(LOG(cst_a)-LOG(3.d0)+4.d0*t(j1)-p(j1))
    !write(*,*)'beta1=',beta1,'cst_a',cst_a,'t(j1)',t(j1),'p(j1)',p(j1)
    beta1=MAX(beta1,1.d-5)

    gamma_gas=32-24*beta1-3*beta1**2/(3*beta1*(8-7*beta1))
    gamma1_dichte=beta1+((4-3*beta1)**2*(gamma_gas-1)/(beta1+12*(gamma_gas-1)*(1-beta1)))

    refad=p(j1)-2.5d0*t(j1)

    refad=refad+3.68d0
     


    IF (refad <= 0.d0) THEN

       IF (num == 0) THEN
          !-----------------------------------------------------------------------
          !  Partie modifiee Mai 1990, D.Schaerer:
          !    IONISATION PARTIELLE DE PLUSIEURS ELEMENTS...

          !  The partial ionisation is computed for the first layers of the interior.
          !  num=0 when entering with the layer j=1.
          !  ionpart is called. As long as x_env(3)=0 or x_env(3) changes from one layer to the next
          !  (comparison wvmyoith vx3 that memorises the previous value), num keeps its value 0.
          !  When x_env(3) doesn't change anymore from one layer to the other, num=j
          !  and the partial ionisation is not computed anymore.
          tion=t(j1)/um
          pion=p(j1)/um
          vmol=vmyo
          vny(1)=vmyo*x(j1)
          vny(2)=vmyo*y(j1)/4.d0
          vny(3)=vny(2)
          xcmp=x(j1)

          CALL ionpart(pion,tion)

          !---------------- Modifications de Schaerer 1990 -----------------------
          ! On prend les valeurs calculees par IONPART pour vmion,vna,vmionp,vmiont
          !      if (x_env(3)-vx3)33,33,32
          IF (x_env(3)-vx3 > 0.d0) THEN
             vmy1=LOG(vmion)
             adi1=vna
             vx3=x_env(3)

             ! modification apportee 19 XI 97
             ! Lorsque FITM a une valeur tres elevee on peut imaginer
             ! qu'il y ait des situations ou x_env(3)=0 et donc ou lors des
             ! premiers passages dans DICHTE x_env(3)=vx3. Les modifications apportees
             ! ici ont pour but de faire en sorte que dans cette situation on
             ! calcule les effets de l'ionisation partielle
          ELSE
             IF (x_env(3) /= 0.d0) THEN
                num=j1
             ELSE
                vmy1=LOG(vmion)
                adi1=vna
                vx3=x_env(3)
             ENDIF
          ENDIF
       ENDIF   ! num

       rh1=-LOG(rgaz)+vmy1+LOG(beta1)+p(j1)-t(j1)
       rhp1=1.d0/beta1+vmionp
       rht1=3.d0-4.d0/beta1+vmiont
       !---------------- Fin des modifications --------------------------------
      !  write(*,*)j,'dichte END: rh1=',rh1

       RETURN
       !------------------------------------------------------------------------
    ELSE   ! refad > 0.
       ! Calcul pour degenerescence
       ! psi <= 7: degenerescence electronique partielle
       ! psi >  7: degenerescence electronique totale
       pl=EXP(p(j1))
       num = -1000
       pg=beta1*pl
       toni=EXP(t(j1))
       !if (j .gt. 700) then
       !write(*,*)'Toni dichte=',toni
      !  else if (j == 751) then
      !  write(*,*)'Toni dichte=',toni
      ! else if (j1 == 749) then
      !  write(*,*)'Toni dichte=',toni
     !endif
       tk=toni*cst_k/(cst_me*cst_c**2.d0)
       mist = 0
       IF (hpsi > 0.d0) THEN
          krueck = 2
          hpsi = -LOG(chi)
          psi = hpsi
          CALL Compute_sepp(krueck,sepp,seppl,dchi)
       ELSE
          krueck = 1
          CALL Compute_sepp(krueck,sepp,seppl,dchi)
       ENDIF

       DO WHILE (mist <= 45)
          mist = mist+1
          IF (krueck == 1) THEN
             IF (chi+dchi <= 0.d0) THEN
                chi = 0.1d0*chi
                IF (chi - 1.d-18 <= 0.d0) THEN
                   krueck = 2
                   hpsi = -LOG(chi)
                   psi = hpsi
                   CALL Compute_sepp(krueck,sepp,seppl,dchi)
                   CYCLE
                ELSE
                   krueck = 1
                   CALL Compute_sepp(krueck,sepp,seppl,dchi)
                   CYCLE
                ENDIF
             ELSE
                chi = chi+dchi
                IF (chi - 1.d-18 <= 0.d0) THEN
                   krueck = 2
                   hpsi = -LOG(chi)
                   psi = hpsi
                   CALL Compute_sepp(krueck,sepp,seppl,dchi)
                   CYCLE
                ELSE
                   IF (ABS(1.d0-sepp)-0.002d0 >= 0.d0) THEN
                      krueck = 1
                      CALL Compute_sepp(krueck,sepp,seppl,dchi)
                      CYCLE
                   ELSE
                      psi=-LOG(chi)
                      chi=chi*0.99d0
                      CALL degen
                      rhes1=rhes
                      chi=chi/0.99d0
                      hpsi=-1.d0
                      CALL degen
                      hpsi=-3.d0
                      rpsi=100.d0*(rhes1-rhes)-rhes
                      rhpsi=-chi*rhes
                      rh1=LOG(rhe*vmye)
                      hifa=rgaz*vermy*rhe*toni/pl
                      pas=hifa*rhes+pe/pl*pes
                      pate=hifa*(1.d0+rhete)+pe*pete/pl+4.d0*(1.d0-beta1)
                      uta=(ue/toni)*(uete-ues*pate/pas)/vmye
                      rhp1=rhes/pas
                      rht1=rhete-rhp1*pate
                      rhpsip=rpsi/pas
                      rhpsit=-pate*rhpsip
                      !write(*,*)j,'dichte END: rh1=',rh1

                      RETURN
                   ENDIF
                ENDIF
             ENDIF
          ENDIF

          DO WHILE (hpsi+dchi <= 0.d0)
             dchi = 0.5d0*dchi
          ENDDO
          hpsi = hpsi+dchi
          chi = EXP(-hpsi)
          psi = hpsi
          IF (hpsi-35.d0 <= 0.d0) THEN
             krueck = 1
             hpsi = -3.d0
             CALL Compute_sepp(krueck,sepp,seppl,dchi)
             CYCLE
          ELSE
             IF (ABS(1.d0-sepp)-0.002d0 >= 0.d0) THEN
                krueck = 2
                CALL Compute_sepp(krueck,sepp,seppl,dchi)
                CYCLE
             ELSE
                hpsi=hpsi+0.1d0
                CALL degen
                rhes1=rhes
                hpsi=hpsi-0.1d0
                CALL degen
                rpsi=(rhes1-rhes)/(0.1d0*chi)
                rhpsi=-rhes
                rh1=LOG(rhe*vmye)
                hifa=rgaz*vermy*rhe*toni/pl
                pas=hifa*rhes+pe/pl*pes
                pate=hifa*(1.d0+rhete)+pe*pete/pl+4.d0*(1.d0-beta1)
                uta=(ue/toni)*(uete-ues*pate/pas)/vmye
                rhp1=rhes/pas
                rht1=rhete-rhp1*pate
                rhpsip=rpsi/pas
                rhpsit=-pate*rhpsip
                !write(*,*)j,'dichte END: rh1=',rh1

                RETURN
             ENDIF
          ENDIF

       ENDDO

       SELECT CASE (krueck)
       CASE (1)
          psi=-LOG(chi)
          chi=chi*0.99d0
          CALL degen
          rhes1=rhes
          chi=chi/0.99d0
          hpsi=-1.d0
          CALL degen
          hpsi=-3.d0
          rpsi=100.d0*(rhes1-rhes)-rhes
          rhpsi=-chi*rhes
       CASE (2)
          hpsi=hpsi+0.1d0
          CALL degen
          rhes1=rhes
          hpsi=hpsi-0.1d0
          CALL degen
          rpsi=(rhes1-rhes)/(0.1d0*chi)
          rhpsi=-rhes
       CASE default
          STOP 'Bad value for krueck in dichte'
       END SELECT
       rh1=LOG(rhe*vmye)
       hifa=rgaz*vermy*rhe*toni/pl
       pas=hifa*rhes+pe/pl*pes
       pate=hifa*(1.d0+rhete)+pe*pete/pl+4.d0*(1.d0-beta1)
       uta=(ue/toni)*(uete-ues*pate/pas)/vmye
       rhp1=rhes/pas
       rht1=rhete-rhp1*pate
       rhpsip=rpsi/pas
       rhpsit=-pate*rhpsip
      !  write(*,*)j,'dichte END: rh1=',rh1
       RETURN
    ENDIF

  END SUBROUTINE dichte
  !======================================================================


END MODULE EOS
