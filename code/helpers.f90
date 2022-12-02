module helpers
   use evol,only: kindreal,ldi,mmax,input_dir,npondcouche,npondcoucheAdv
   use const,only: um,cst_a,lgLsol,cstlg_sigma,cstlg_G,lgMsol,cst_G,Msol,pi,lgRsol,Rsol,qapicg,xlsomo,year,day,Lsol,cstlg_K1, &
      cstlg_mH,cstlg_k,cst_sigma
   use inputparam,only: modanf,nwseq,nzmod,iprn,iauto,ialflu,ianiso,imagn,ipop3,irot,isol,idiff,iadvec,icoeff, &
      igamma,ibasnet,istati,iledou,idifcon,iover,iunder,my,ikappa,iopac,imloss,ifitm,itmin,nndr,idialo,idialu,phase,isugi,nbchx, &
      nrband,iout,icncst,islow,ichem,zinit,zsol,z,frein,elph,dovhp,dunder,fmlos,fitm,rapcrilim,omega,xfom,vwant,gkorm,alph, &
      agdr,agds,agdp,agdt,faktor,deltal,deltat,dgrp,dgrl,dgry,dgrc,dgro,dgr20,xdial,fenerg,richac,xcn,idern,display_plot, &
      itminc,idebug,FITM_Change,IMLOSS_Change,Write_namelist,Read_namelist,starname,xyfiles,idebug,&
      bintide,binm2,periodini,verbose,Add_Flux
   use inputparam,only: add_diff,B_initial,const_per,diff_only,itests,K_Kawaler,Omega_saturation,&
      stop_deg,tauH_fit,var_rates,RSG_Mdot,SupraEddMdot,Be_mdotfrac,start_mdot
   use inputparam,only: ipoly,n_snap
   use caramodele,only: xLtotbeg,dm_lost,inum,nwmd,xmini,firstmods,eddesc,hh6,glm,xLstarbefHen,hh1,iwr,xmdot,rhoc,tc,gls,teff, &
      glsv,teffv,ab,gms,zams_radius,Mdot_NotCorrected
   use abundmod,only: x,y3,y,xc12,xc13,xc14,xn14,xn15,xo16,xo17,xo18,xf18,xf19,xne20,xne21,xne22,xna23,xmg24,xmg25,xmg26,xal26, &
      xal27,xsi28,xprot,xneut,xbid,xbid1,vx,vy3,vy,vxc12,vxc13,vxc14,vxn14,vxn15,vxo16,vxo17,vxo18,vxf18,vxf19,vxne20,vxne21, &
      vxne22,vxna23,vxmg24,vxmg25,vxmg26,vxal26g,vxal27,vxsi28,vxprot,vxneut,vxbid,vxbid1,ekrote,epote,ekine,erade,snube7,snub8, &
      nbelx,nbzel,nbael,zabelx,abels,abelx,vabelx,mbelx,maxCNO,abundCheck,lcnom,xmcno,scno
   use equadiffmod,only: ccg1,ccg2,ccg3,ccz2,ccz3,gkorv,iprc,gkor,iter
   use strucmod,only: m,q,p,t,r,s,vp,vt,vr,vs,e,rho,zensi,rprov,ccrad1,NPcoucheEff,id1,id2,drl,drte,dk,drp, &
      drt,drr,rlp,rlt,rlc,rrp,rrt,rrc,rtp,rtt,rtc,chem,ychem,neudr,fitmion,Nabla_mu,vna,vnr
   use rotmod,only: CorrOmega,dlelex,suminenv,vsuminenv,vvsuminenv,omegi,vomegi,rapcri,xobla,rapom2,alpro6,do1dr,bmomit,&
      btot,btotatm,Flux_remaining,BTotal_EndAdvect,BTotal_StartModel,dlelexsave,timestep_control,xldoex
   use timestep,only: alter,dzeitj,dzeit,dzeitv
   use convection,only: bordn,jwint,xzc,ixzc,qbc,qmnc,CZdraw,BaseZC,iidraw,drawcon,r_core
   use omegamod,only: vcritcalc,omescale,dlonew,omconv,momevo,omenex,om2old,momspe,xjspe1,xjspe2
   use envelope,only: dreckf,dreck,notFullyIonised,supraEdd
   use ionisation,only: abond,list,iatoms
   use diffadvmod,only: tdiff,jdiff
   use energy,only: enint,netinit,vmassen,rvect,t9n,pvect,epstot1,epsneut,dcoeff
   use geomod, only: rpsi_min,initgeo,geomat,geomeang
   use PGPlotModule, only: restart,InitPGplot,SavePlotData,EndPGplot,Chem_Species_Number
   use SmallFunc,only: exphi
   use LayersShift,only: fitmshift,schrit,mdotshift
   use winds,only: aniso,xloss,xldote,corrwind
   use chemicals,only: netnew,chemeps,chemold
   use diffusion,only: coedif,diffbr
   use timestep,only: zeit
   use henyey_solver,only: henyey,nsugi,correction_message,henyey_last
   use opacity,only: ioutable,rout,tout
   use nablas,only: grapmui
   use PrintAll, only: File_Unit,PrintCompleteStructure
   use WriteSaveClose,only: OpenAll,CheckSchrit,write4,read4,SequenceClosing,nzmodini,nzmodnew
   use bintidemod,only: period
   use genec, only: elemneg,checkVink,ivcalc,veryFirst,TriangleIteration
   use genec, only: allam,bibib,bolm,fffff,dlelexprev,dmneed,eddesm,fmain,glsvv,h1,h2,hr,opaesc, &
      rap2,rap1,radius,rapg,rapomm,raysl,teffeq,rrro,teffvv,teffel,teffpr,vcrit1,tzero,vcri2m, &
      vcri1m,vequat,vcrit2,vequam,vpsi,xdilto,xdilex,xft,xgmoym,xini,xltof,xltod,xltot,xmdotneed,xmdotwr,xo1, &
      xogtef,xpsi,xrequa,xtt,xtod2,zwi1,ygmoye,xdippp,ygequa,zwi,rhocprev,Tcprev
   use inichemmod, only: inichem,idefaut,mainnam,xx,zini,znew,elemZ,elemA
   use const, only: pi,lgpi,cst_G,Msol,Rsol,Lsol,lgLsol,year,cst_mh,cst_k,cstlg_sigma
   use interpolation, only: fipoi
   use modinimod, only: diminipetit,dimini,dimdat,&
      polytrop,writetable,teffdat,lumdat,massdat,&
      Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,P1,P2,P3,P4,P5,P6,P7,P8,&
      T1,T2,T3,T4,T5,T6,T7,T8,R1,R2,R3,R4,R5,R6,R7,R8,&
      S1,S2,S3,S4,S5,S6,S7,S8
   !use State, only: conditioned_stop, stopping_condition

!  use inputparam
   implicit none
   real(8) :: mstar

   integer, parameter::n_dim=10001
   real(8), parameter::musol=0.6074202636615116d0
   real(8):: index_poly,Lstar,xteff,rstar,alpha,rhomoy,rhocrho,ka,mu,normC,deltaq
   real(8):: ztest
   real(8), allocatable::xi(:),theta(:),dthetadxi(:),pression(:),temp(:),xr(:),xmr(:),grav(:),xlum(:),qq(:)
   real(8), dimension(50)::rh

   integer:: i,ll,ii,iprnv,iterv,k,nfseq,j,imlosssave,modell

   integer:: Iteration48,IterTriangle,ielemneg

   real(kindreal):: summas
   real(kindreal), dimension(5):: xnetalu
   real(kindreal), dimension(npondcouche):: CorrZero

   namelist/IniStruc/gms,alter,gls,teff,glsv,teffv,dzeitj,dzeit,dzeitv,summas,ab,m,q,p,t,r,s,vp,vt,vr,vs,x,y3,y,xc12,xc13,&
      xn14,xn15,xo16,xo17,xo18,xne20,xne22,xmg24,xmg25,xmg26,omegi

contains
   subroutine makeini
      implicit none

      integer::i,jmax,longueur

      allocate(xi(n_dim))
      allocate(theta(n_dim))
      allocate(dthetadxi(n_dim))
      allocate(pression(n_dim))
      allocate(temp(n_dim))
      allocate(xr(n_dim))
      allocate(xmr(n_dim))
      allocate(grav(n_dim))
      allocate(xlum(n_dim))
      allocate(qq(n_dim))

      rho=0.d0
      xr=0.d0
      xmr=0.d0
      mu=0.d0
      grav=0.d0
      qq=0.d0

      ! TODO: require starname, mstar and zini to be set beforehand by AMUSE
      ! write(*,*)'Enter the star name:'
      ! read(5,*) starname
      ! write(*,*)'Enter the desired mass and metallicity:'
      ! read(5,*) mstar,zini

      Lstar=fipoi(mstar,dimdat,massdat,lumdat)
      xteff=fipoi(mstar,dimdat,massdat,teffdat)
      rstar=(Lstar+lgLsol-4.d0*xteff-log10(4.d0)-lgpi-cstlg_sigma)/2.d0
      rstar=10.d0**rstar

      write(*,*)'Star ',trim(starname)
      write(*,'(a,4(1x,f9.3))')'Mass, L/Lo, log Teff, Rstar:',mstar,Lstar,xteff,rstar/Rsol

      if (mstar<=40.d0) then
         elph=1.60d0
         my=0
      else
         elph=1.00d0
         my=1
      endif
      if (mstar<7.d0) then
         imloss=0
      else
         imloss=6
      endif
      if (mstar<11.d0) then
         ialflu=0
      else
         ialflu=1
      endif
      if (mstar>=1.7d0) then
         dovhp=0.10d0
      else if (mstar>=1.25d0) then
         dovhp=0.05d0
      else
         dovhp=0.d0
      endif
      if (mstar<=3.d0) then
         dzeitj=1.0d4
      else if (mstar<=7.d0) then
         dzeitj=5.0d3
      else if (mstar<=25.d0) then
         dzeitj=1.0d3
      else
         dzeitj=1.0d2
      endif
      dzeit=dzeitj*year
      dzeitv=dzeit/2.d0

      !TODO set this via AMUSE
      !write(*,*) 'Which rotation velocity on the ZAMS?'
      !read(5,*) vwant
      !END TODO
      if (abs(vwant) > epsilon(0.d0)) then
         irot=1
         isol=1
         fitm=0.99990d0
         ifitm=3
         rapcrilim=0.99d0
         omega=1.d-5
         !write(inifilename,'(a4,a,a4)') 'ini_',trim(starname),'.rot'
      else
         irot=0
         isol=0
         fitm=0.980d0
         ifitm=0
         rapcrilim=0.d0
         omega=0.d0
         vwant = 0.0d0
         !write(inifilename,'(a4,a,a4)') 'ini_',trim(starname),'.com'
      endif

      call inichem

      !TODO
      !select case (idefaut)
      !  case (0)
      !    write(*,*)'Do you want to compute a polytropic structure or use ',&
      !              'a pre-computed structure? Polytrope:1 - structure:0'
      !    read(5,*) ipoly
      !  case (1)
      !    ipoly = 0
      !end select
      !END TODO
      select case (ipoly)
       case (0)
         if(mstar <= 2.d0) then
            longueur=diminipetit
            q(1:longueur)=q1
            p(1:longueur)=p1
            t(1:longueur)=t1
            r(1:longueur)=r1
            s(1:longueur)=s1
         elseif(mstar <= 5.d0) then
            longueur=dimini
            q(1:longueur)=q2
            p(1:longueur)=p2
            t(1:longueur)=t2
            r(1:longueur)=r2
            s(1:longueur)=s2
         elseif(mstar <= 12.d0) then
            longueur=dimini
            q(1:longueur)=q3
            p(1:longueur)=p3
            t(1:longueur)=t3
            r(1:longueur)=r3
            s(1:longueur)=s3
         elseif(mstar <= 20.d0) then
            longueur=dimini
            q(1:longueur)=q4
            p(1:longueur)=p4
            t(1:longueur)=t4
            r(1:longueur)=r4
            s(1:longueur)=s4
         elseif(mstar <= 32.d0) then
            longueur=dimini
            q(1:longueur)=q5
            p(1:longueur)=p5
            t(1:longueur)=t5
            r(1:longueur)=r5
            s(1:longueur)=s5
         elseif(mstar <= 40.d0) then
            longueur=dimini
            q(1:longueur)=q6
            p(1:longueur)=p6
            t(1:longueur)=t6
            r(1:longueur)=r6
            s(1:longueur)=s6
         elseif(mstar <= 60.d0) then
            longueur=dimini
            q(1:longueur)=q7
            p(1:longueur)=p7
            t(1:longueur)=t7
            r(1:longueur)=r7
            s(1:longueur)=s7
         else
            longueur=dimini
            q(1:longueur)=q8
            p(1:longueur)=p8
            t(1:longueur)=t8
            r(1:longueur)=r8
            s(1:longueur)=s8
         endif
         q(1) = log10(1.d0-fitm)
       case (1)
         !write(*,*)'Enter the polytropic index (recommended: 2.5):'
         !read(5,*) index_poly
         longueur=50

         call polytrop(index_poly,xi,theta,dthetadxi,n_dim,jmax)

         do i=1,20
            mu=mu + xx(i)*(1.d0+elemZ(i))/elemA(i)
         enddo
         mu=1.d0/mu
         mu = musol
         write(*,*) 'mu=',mu

         alpha=xi(jmax-1)/rstar

         rhomoy=3.d0*mstar*Msol/(4.d0*pi*rstar**3.d0)
         rhocrho=-xi(jmax-1)/(3.d0*dthetadxi(jmax-1))
         rhoc=rhocrho*rhomoy

         ka=4.d0*pi*cst_G*rhoc**(1.d0-1.d0/index_poly)/((index_poly+1.d0)*alpha**2.d0)

         do i=1,jmax-1
            rho(i)=rhoc*theta(i)**index_poly
            pression(i)=ka*rho(i)**(1.d0+1.d0/index_poly)
            if (i ==1) then
               temp(i)=pression(i)*mu*cst_mh/(cst_k*rho(i))
            else
               temp(i)=temp(1)*theta(i)
            endif
            xr(i)=xi(i)/alpha
            if(i==1) then
               xmr(i) = 0.d0
            else if (i==2) then
               xmr(i)=(4.d0/3.d0)*pi*xr(i)**3.d0*((rho(i)+rho(i-1))/2.d0)
            else
               xmr(i)=xmr(i-1)+(4.d0/3.d0)*pi*(xr(i)**3.d0-xr(i-1)**3.d0)*((rho(i)+rho(i-1))/2.d0)
            endif
         enddo

         do i=1,jmax-1
            xmr(i)=xmr(i)*mstar*Msol/xmr(jmax-1)
            if(i<jmax-1) then
               qq(i)=log10(1.d0-(xmr(i)/(mstar*Msol)))
            endif
            grav(i) = grav(i-1)+(temp(i)+temp(i-1))*(xmr(i)-xmr(i-1))/2.d0
         enddo
         xmr(jmax-1)=mstar*Msol
         qq(jmax-1)=qq(jmax-2)

         normC = 10.d0**Lstar*Lsol/grav(jmax-1)
         do i=1,jmax-1
            pression(i)=log10(pression(i))
            temp(i)=log10(temp(i))-0.1d0
            rho(i)=log10(rho(i))
            if (i == 1) then
               xr(i) = 0.d0
               xlum(i) = 0.d0
            else
               xr(i)=log10(xr(i))
               xlum(i) = log10(grav(i)*normC)
            endif
         enddo

         q(longueur) = 0.d0
         q(longueur-1) = -2.d-4
         q(longueur-2) = -5.d-4
         q(longueur-3) = -9.d-4
         q(longueur-4) = -2.d-3
         q(longueur-5) = -5.d-3
         q(33)= log10(1.d0-1.d0/5.d0)
         q(16)= log10(1.d0-1.d0/2.d0)
         q(1)= log10(1.d0-fitm)

         r(longueur) = 7.0d0
         p(longueur) = pression(1)
         t(longueur) = temp(1)
         deltaq= (q(16)-q(1))/15.d0
         do i=2,15
            q(i)=q(i-1)+deltaq
         enddo
         deltaq= (q(33)-q(16))/17.d0
         do i=17,32
            q(i)=q(i-1)+deltaq
         enddo
         deltaq= (q(longueur-5)-q(33))/12.d0
         do i=34,longueur-6
            q(i)=q(i-1)+deltaq
         enddo

         do i=1,longueur-1
            r(i)=fipoi(q(i),jmax-1,qq(1:jmax-1),xr(1:jmax-1))
            s(i)=fipoi(q(i),jmax-1,qq(1:jmax-1),xlum(1:jmax-1))
            p(i)=fipoi(q(i),jmax-1,qq(1:jmax-1),pression(1:jmax-1))
            t(i)=fipoi(q(i),jmax-1,qq(1:jmax-1),temp(1:jmax-1))
            rh(i)=fipoi(q(i),jmax-1,qq(1:jmax-1),rho(1:jmax-1))
         enddo
         s(longueur) = 0.95d0*s(longueur-1)
       case default
         !stopping_condition = 'Bad choice for structure type, must be 0 or 1'
         !call conditioned_stop()
         stop
         return
      end select

      ! open(21,file=inifilename,iostat=ierror,status='unknown')

      nwseq = 1
      modanf = 0
      ianiso = 0
      phase = 1
      zinit = zini
      z = znew
      idiff = 0
      iadvec = 0
      icoeff = 11
      xdial = 0.d0
      idialo = 0
      idialu = 0
      fmlos = 0.85d0
      deltal = 0.02d0
      deltat = 0.02d0
      gkorm = 9.d0
      alph = 0.3d0
      agdr = 1.d-5
      faktor = 1.d0
      dgrp = 0.01d0
      dgrl = 0.01d0
      dgry = 0.003d0
      dgrc = 0.01d0
      islow = 2
      xcn = 1.d0
      display_plot = .false.
      iauto = 1
      ! call Write_namelist(21,nwseq,modanf,10,xcn)
      nzmod = 10
      if (abs(zinit) > epsilon(0.d0)) then
         ipop3 = 0
      else
         ipop3 = 1
      endif

      ! write(21,'(a)') ' &IniStruc'

      ! write(21,'(a,f9.3,a,1pd8.2,a,0pf6.0,a)') ' GMS=',mstar,&
      !    'd0, ALTER=0.d0, GLS=',10.d0**Lstar,', TEFF=',10.d0**xteff,'d0,'
      gms = mstar
      alter = 0.d0
      gls = 10.d0**Lstar
      teff = 10.d0**xteff
      ! write(21,'(25x,a,1pd8.2,a,0pf6.0,a)') 'GLSV=',10.d0**Lstar,', TEFFV=',10.d0**xteff,'d0,'
      glsv = 10.d0**Lstar
      teffv = 10.d0**xteff
      ! write(21,'(a,1pd8.2,a,d10.4,a)') ' DZEITJ=',dzeitj,', DZEIT=',dzeit,','
      ! write(21,'(18x,a,1pd10.4,a)') 'DZEITV=',dzeitv,','
      ! write(21,'(a,f9.3,a,i2,a)') ' SUMMAS=',mstar,'d0, AB=0.d0, M=',longueur,','
      summas = mstar
      xmini = mstar
      ab=0.d0
      m=longueur
      ! write(21,'(a)') ' Q='
      ! call writetable(q,longueur)
      ! write(21,'(a)') ' P='
      ! call writetable(p,longueur)
      ! write(21,'(a)') ' T='
      ! call writetable(t,longueur)
      ! write(21,'(a)') ' R='
      ! call writetable(r,longueur)
      ! write(21,'(a)') ' S='
      ! call writetable(s,longueur)
      ! write(21,'(a)') ' VP='
      ! write(21,'(1x,i2,a)') longueur,'*0.0d0,'
      vp = 0.d0
      ! write(21,'(a)') ' VT='
      ! write(21,'(1x,i2,a)') longueur,'*0.0d0,'
      vt = 0.d0
      ! write(21,'(a)') ' VR='
      ! write(21,'(1x,i2,a)') longueur,'*0.0d0,'
      vr = 0.d0
      ! write(21,'(a)') ' VS='
      ! write(21,'(1x,i2,a)') longueur,'*0.0d0,'
      vs=0.d0

      ztest=1.d0

      do i=1,longueur
         x(i) = xx(1)
         y3(i) = xx(2)
         y(i) = xx(3)
         xc12(i) = xx(4)
         xc13(i) = xx(5)
         xn14(i) = xx(6)
         xn15(i) = xx(7)
         xo16(i) = xx(8)
         xo17(i) = xx(9)
         xo18(i) = xx(10)
         xne20(i) = xx(11)
         xne22(i) = xx(12)
         xmg24(i) = xx(13)
         xmg25(i) = xx(14)
         xmg26(i) = xx(15)
      enddo
      do i=1,15
         !write(21,'(a6,1x,a1,i2,a1,1pd21.15,a1)') mainnam(i),'=',longueur,'*',xx(i),','
         ztest=ztest-xx(i)
      enddo
      !write(21,'(a6,1x,a1,i2,a1,1pd21.15,a1)')' omegi','=',longueur,'*',omega,','
      omegi = longueur * omega
      !write(21,'(a)') ' &END'
      write(*,*)'Ztest=',ztest,'=? Znew=',znew
      !close(21)

      !write (*,*) 'file: ',trim(inifilename),' done.'
   end subroutine makeini

   subroutine set_defaults  ! replaces Read_namelist
      ! For now: set defaults based on an ini file for a 7MSun star with solar metallicity (0.014) and no rotation
      ! TODO: replace with tool that calculates these setting based on stellar mass, metallicity and rotation rate
      ! (and other configuration options given in make_inifile)
      implicit none
      ! * CharacteristicsParams namelist *
      starname = "AmuseDefaultStar"
      nwseq = 1
      modanf = 0
      nzmod = 1 ! number of steps calculated before stopping
      n_snap = 0 ! don't write snapshots

      ! * PhysicsParams namelist *
      irot = 0
      isol = 0
      imagn = 0
      ialflu = 0
      ianiso = 0
      ipop3 = 0
      ibasnet = 0
      phase = 1
      var_rates = .false.
      bintide = .false.
      binM2 = 0.0D+00
      periodini = 0.0D+00
      const_per = .true.

      ! * CompositionParams namelist *
      zinit = 1.40D-02
      zsol = 0.140D-01
      z = 0.261956006494962D-02
      iopac = 3
      ikappa = 5

      ! * RotationParams namelist *
      idiff = 0
      iadvec = 0
      istati = 0
      icoeff = 11
      fenerg = 0.100D+01
      richac = 0.100D+01
      igamma = 0
      frein = 0.000D+00
      K_Kawaler = 0.000D+00
      Omega_saturation = 0.140D+02
      rapcrilim = 0.00000
      vwant = 0.000D+00
      xfom = 0.100D+01
      omega = 0.000000000000000E+00
      xdial = 0.000
      idialo = 0
      idialu = 0
      Add_Flux = .true.
      diff_only = .false.
      B_initial = 0.000D+00
      add_diff = 0.000D+00

      ! * SurfaceParams namelist *
      imloss = 6
      fmlos = 0.850D+00
      RSG_Mdot = 0
      SupraEddMdot = .true.
      ifitm = 0
      fitm = 0.980000000
      !fitmi = 0.980000000
      deltal = 0.02000
      deltat = 0.02000
      nndr = 1
      Be_mdotfrac = 0.0d0
      start_mdot = 0.80d0

      ! * ConvectionParams namelist *
      iledou = 0
      idifcon = 0
      elph = 1.600
      my = 0
      iover = 1
      dovhp = 0.100
      iunder = 0
      dunder = 0.000

      ! * ConvergenceParams namelist *
      gkorm = 9.000
      alph = 0.300
      agdr = 0.10D-04
      faktor = 1.00E+00
      dgrp = 0.0100
      dgrl = 0.0100
      dgry = 0.00300
      dgrc = 0.01000
      dgro = 0.01000
      dgr20 = 0.100D-01
      nbchx = 200
      nrband = 1

      ! * TimeControle namelist *
      islow = 2
      xcn = 1.000
      icncst = 0
      tauH_fit = 1

      ! * VariousSettings namelist *
      display_plot = .false.
      iauto = 1
      iprn = 10
      iout = 0
      itmin = 5
      xyfiles = .false.
      idebug = 0
      itests = 0
      verbose = .false.
      stop_deg = .true.

   end subroutine set_defaults

end module helpers
