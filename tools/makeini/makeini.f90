program makeini

  use const, only: pi,lgpi,cst_G,Msol,Rsol,Lsol,lgLsol,year,cst_mh,&
                   cst_k,cstlg_sigma
  use inichemmod, only: inichem,idefaut,mainnam,xx,zini,znew,elemZ,elemA
  use interpolmod, only: fipoi
  use modinimod, only: diminipetit,dimini,dimdat,&
                       polytrop,writetable,teffdat,lumdat,massdat,&
                       Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,P1,P2,P3,P4,P5,P6,P7,P8,&
                       T1,T2,T3,T4,T5,T6,T7,T8,R1,R2,R3,R4,R5,R6,R7,R8,&
                       S1,S2,S3,S4,S5,S6,S7,S8

  use inputparam

  implicit none

  integer, parameter::n_dim=10001
  integer::i,jmax,ierror,ipoly,longueur

  real(kindreal), parameter::musol=0.6074202636615116d0
  real(kindreal):: mstar,dzeitj,dzeit,dzeitv,n,Lstar,xteff,rstar,alpha,&
                     rhomoy,rhocrho,rhoc,ka,mu,normC,deltaq
  real(kindreal):: ztest
  !real(8), dimension(n_dim)::xi,theta,dthetadxi,rho,pression,temp,xr,xmr,grav,xlum,qq
  real(kindreal), allocatable::xi(:),theta(:),dthetadxi(:),rho(:),pression(:),temp(:),xr(:)
  real(kindreal), allocatable::xmr(:),grav(:),xlum(:),qq(:)
  real(kindreal), dimension(50)::q,r,s,p,t,rh

  character(256)::inifilename

  allocate(xi(n_dim))
  allocate(theta(n_dim))
  allocate(dthetadxi(n_dim))
  allocate(rho(n_dim))
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

  write(*,*)'Enter the star name:'
  read(5,*) starname
  write(*,*)'Enter the desired mass and metallicity:'
  read(5,*) mstar,zini

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

  write(*,*) 'Which rotation velocity on the ZAMS?'
  read(5,*) vwant
  if (abs(vwant) > epsilon(0.d0)) then
    iprezams=1
    irot=1
    isol=1
    fitm=0.99990d0
    ifitm=3
    rapcrilim=0.99d0
    omega=1.d-5
    write(inifilename,'(a4,a,a4)') 'ini_',trim(starname),'.rot'
  else
    iprezams=0
    irot=0
    isol=0
    fitm=0.980d0
    ifitm=0
    rapcrilim=0.d0
    omega=0.d0
    vwant = 0.0d0
    write(inifilename,'(a4,a,a4)') 'ini_',trim(starname),'.com'
  endif

  call inichem

  select case (idefaut)
    case (0)
      write(*,*)'Do you want to compute a polytropic structure or use ',&
                'a pre-computed structure? Polytrope:1 - structure:0'
      read(5,*) ipoly
    case (1)
      ipoly = 0
  end select
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
      write(*,*)'Enter the polytropic index (recommended: 2.5):'
      read(5,*) n
      longueur=50

      call polytrop(n,xi,theta,dthetadxi,n_dim,jmax)

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

      ka=4.d0*pi*cst_G*rhoc**(1.d0-1.d0/n)/((n+1.d0)*alpha**2.d0)

      do i=1,jmax-1
       rho(i)=rhoc*theta(i)**n
       pression(i)=ka*rho(i)**(1.d0+1.d0/n)
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
      stop 'Bad choice for structure type, must be 0 or 1'
  end select

  open(21,file=inifilename,iostat=ierror,status='unknown')

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
  call Write_namelist(21,nwseq,modanf,10,xcn)

  write(21,'(a)') ' &IniStruc'

  write(21,'(a,1pd21.15,a,d21.15,a,d21.15)') ' GMS=',mstar,&
     ', ALTER=0.d0, GLS=',10.d0**Lstar,', TEFF=',10.d0**xteff
  write(21,'(25x,a,1pd21.15,a,d21.15)') 'GLSV=',10.d0**Lstar,', TEFFV=',10.d0**xteff
  write(21,'(a,1pd21.15,a,d21.15,a)') ' DZEITJ=',dzeitj,', DZEIT=',dzeit,','
  write(21,'(18x,a,1pd21.15,a)') 'DZEITV=',dzeitv,','
  write(21,'(a,d21.15,a,i2,a)') ' SUMMAS=',mstar,', AB=0.d0, M=',longueur,','
  write(21,'(a)') ' Q='
  call writetable(q,longueur)
  write(21,'(a)') ' P='
  call writetable(p,longueur)
  write(21,'(a)') ' T='
  call writetable(t,longueur)
  write(21,'(a)') ' R='
  call writetable(r,longueur)
  write(21,'(a)') ' S='
  call writetable(s,longueur)
  write(21,'(a)') ' VP='
  write(21,'(1x,i2,a)') longueur,'*0.0d0,'
  write(21,'(a)') ' VT='
  write(21,'(1x,i2,a)') longueur,'*0.0d0,'
  write(21,'(a)') ' VR='
  write(21,'(1x,i2,a)') longueur,'*0.0d0,'
  write(21,'(a)') ' VS='
  write(21,'(1x,i2,a)') longueur,'*0.0d0,'

  ztest=1.d0
  do i=1,15
   write(21,'(a6,1x,a1,i2,a1,1pd21.15,a1)') mainnam(i),'=',longueur,'*',xx(i),','
   ztest=ztest-xx(i)
  enddo
  write(21,'(a6,1x,a1,i2,a1,1pd21.15,a1)')' omegi','=',longueur,'*',omega,','
  write(21,'(a)') ' &END'
  write(*,*)'Ztest=',ztest,'=? Znew=',znew
  close(21)

  write (*,*) 'file: ',trim(inifilename),' done.'
end program makeini
