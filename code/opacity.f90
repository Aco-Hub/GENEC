!*********************************************************************
!*********************************************************************
module opacity
! Version modifiee par Georges Meynet, Octobre 2000
! pour utilisation avec fichier contenant opacites OPAL
! + Alexander & Ferguson a basse temperature.
! (Fichier solar.dat)

!     Dans cette version les tables sont etendues a basse T de
!     log T = 3.00 a 3.70 en pas de 0.05 dex (15 lignes supplem.).

!     !!! A utiliser pour melange solaire uniquement conjointement
!     !!! avec la routine kappa96__solar.f.

!     Les modifications sont marquees par 'cds'...

! Double precision + modules avec les constantes : CG 02.04.09

!-----VERSION of October 5, 1995-----------------------------------------
!-----------------------------------------------------------------------

!          This subroutine contains instructions for using the subroutine
!     OPACGN93( z, xh, t6, R ) and OPAC(izi,mzin,xh,t6,r).
!          The purpose of these subroutines is to perform 3 or 4 variable
!     interpolation on log10(kappa).  The opacity tables to be interpolated
!     are known to have somewhat random numerical errors of a few percent.
!     Consequently, adjusting the data prior to performing the interpolation is
!     justified at this level.  These codes are set-up to read the original
!     (unsmoothed) tabular data, this data is then passed through a smoothing
!     filter; using a set of routines developed by Mike Seaton (see M. J. Seaton
!     MNRAS 265,L25(1993)). It is the adjusted data that is actually used in
!     carrying out the interpolations in OPACGN93 and OPAC.  The interpolation
!     method, described below,  also produces smooth values for logK and its
!     first derivatives. The initial adjustment step helps improve the smoothness
!     of the OPACGn93 and OPAC  output,particularly at the smallest values of R.
!     The medium to large R output is only slightly effected by this step.  It
!     takes only a few seconds to carryout the initial data smoothing step. This
!     step can be skipped by setting the parameter ismdata in subroutine
!     readco =1.
!GMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGM
! ICI ON A MIS READCO=1
! LES SS-ROUTINES UTILISEES SONT opacgn93,opac,readco,t6rinterp,
!                                la fonction quad et le block data
! NOUS N'UTILISONS PAS LE LISSAGE sous routine opaltab et suivantes
! Nous avons verifie que les interpolations fonctionnent.
! Ces test ont permis de corriger une erreur due a la mise
! en double precision. Dans opacgn93 on a declarer reel les variables
! kapz, kapz1 et kapz2. Mais du fait de la double precision implicite
! a-h,o-z ces variables ne sont pas en double precision. Cela pose
! probleme lorsque kapz (qui est un tableau) est appele en argument
! de quad. J'ai change le nom de cette variable devenue xkapz et
! le pgm fonctionne
!GMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGMGM
!          The interpolation variables are :

!          z       The metallicity, Z
!          xh      The hydrogen mass fraction, X
!          t6      The temperature in millions of degrees Kelvin, T6
!          r       =rho(gm/cc)/T6**3, R

!          Additional input to OPAC is:

!          izi     Keeps or recalculates table indices. The value 0 causes
!                  the table indices to be recalculated.  A value other than 0
!                  causes the previous indices to be used.

!          mzin    The integer value of i of the Z value to use.  The

!                  choices are:
!                  1=0.0  2=0.0001 3=0.0003 4=0.001 5=0.002 6=0.004 7=0.01
!                  8=0.02 9=0.03  10=0.04  11=0.06 12=0.08 13=0.1


!          An interpolation between overlapping quadratics is used to obtain
!     smoothed results.  A 4x4 grid in logT6 and logR is used to interpolate
!     in four different 3x3 sub-grids. Linear interpolation between quadratic
!     fits in these different sub-grids gives smoothed results in both log T6
!     and Log R directions. Compared to ealier versions of this code, the
!     interpolation in Z is in Kappa vs. Z; not log Kappa vs. Z.
!     The overlapping quadratic procedure produces results that are  smooth,
!     similar to bicubic spline interpolation, but require storage of only local
!     information.
!
!          The code OPACGN93 performs interpolation in Z, X, T6, and R. It calls
!     the subroutine OPAC at each Z.  If you are working with a fixed Z, as
!     listed above, it is more efficient to call OPAC directly. In this case use
!     izi=0 and the appropiate value of mzin. The  opacity data will be read
!     from unit 2 in the subroutine readco.  you will need to have the file
!     GN93hz available on disk.
!         Each of the individual tables in the GN93hz file cover 70 temperatures
!     in the range logT=3.75[T6=0.0056341325]( referred to as temperature 1) to
!     logT=8.7[T6=500]. and 19 values of log R in the range -8 (referred to as 1)
!     to +1; at half-integer steps.  (NOTE: earlier tables were explicitly in
!     terms of T6. For convenience the present tables tabulate log Kappa vs logT. The
!     interpolation however still uses T6 for the temperature variable)
!     For specialized problems, if storage space is a problem, a reduced set of
!     data can be input .  This requires a recompilation with altered parameter
!     values.  In order to limit the range of R, set the parameter nrb= index of
!     lowest value of log R to use(count from log R=-8).  Then set the parameter
!     nre to the index of the largest value of log R to use.  (NOTE: nre-nrb must
!     be at least 5). To ignore the first ntb-1 values of T6 (starting from
!     temperature 1) set the parameter ntb to desired value.  No other parameters
!     should be modified.

!     ***CAUTION***
!         As a result of the mixing procedure used to calculate the data a few
!     X=0.0, low T-small R, table values fell outside the range of T and R accessible
!     from the X=0.35 data directly calculated for this purpose.  These T-R
!     locations are filled in with 9.99 (or for diagnostic purposes in some cases
!     larger values.  At the same locations the derivatives are set to 99.9.  When
!     T-R falls in this region a message is issued by the interpolation code to
!     inform the user of this situation.  Presumable very few users will have
!     applications that take them into this region.

!          Your routine that performs the call to OPAC should include the
!      module xztrinmod:

!         The variables have the following meanings:
!
!         OPACT        Is the Log of the Rosseland mean opacity: Log(kappa)
!         DOPACT      Is Dlog(kappa)/Dlog(T6)   at constant R
!         DOPACR      Is Dlog(kappa)/Dlog(R)    at constant T
!         DOPACTD     Is Dlog(kappa)/Dlog(T6)   at constant Rho
!*********************************************************************
  use evol,only: ldi,kindreal,input_dir
  use inputparam,only: verbose

  implicit none

! Public variables, communicating with other routines:
  integer,save,public:: error_xzt

  real(kindreal),save,public:: opact,dopact,dopacr,dopactd

! Private variables, invisible from elsewhere out this module:

  integer,parameter,private:: mx=10,mz=13,nrm=19,nrb=1,nre=19,nr=nrm+1-nrb,ntm=85,ntb=1, &
                                nt=ntm+1-ntb,ip=100,ipr=20,nt2=57,nr2=17,nx2=3,nz2=12,irmaxtot=17,&
                                irmaxint=13,jzmax=12,nrh=25,nsum=10111
  integer,save,private:: mzz,m,mf,m1,l1,l2,l3,l4,k1,k2,k3,k4,ip_op,iq_op,ntemp,nsm,nrlow,nrhigh, &
                           nrl,nset
  integer,dimension(mx),parameter,private:: n=(/13,13,13,13,13,13,13,13,10,12/)
  integer,dimension(nr),parameter,private:: nta=(/85,85,85,85,85,85,85,85,85,85,85,85,85,85,84, &
                                                     79,75,73,72/)
  integer,dimension(mx,mz),save,private:: itab



  real(kindreal),save,private:: zval,xxh,dkap,rle,rls,tmax
  real(kindreal),dimension(4),save,private:: qxz,hxz
  real(kindreal),parameter,private:: rmin=-7.d0,rmax=1.d0
  real(kindreal),dimension(nx2),parameter:: xk_kap=(/0.7d0,0.35d0,0.0d0/)
  real(kindreal),dimension(6),parameter,private::gam=(/+0.0073469388d0,-0.0293877551d0,-0.0416326531d0, &
                                                        +0.1175510204d0,+0.1665306122d0,+0.2359183673d0/)

  real(kindreal),dimension(11),parameter,private:: alp=(/-0.0844897959d0,-0.0048979592d0, &
                   +0.0073469388d0,+0.0012244898d0,0.3379591837d0,+0.0195918367d0,-0.0293877551d0, &
                   +0.4787755102d0,0.0277551020d0,-0.0416326531d0,-0.0069387755d0/), &
             bet=(/-0.0048979592d0,-0.0661224490d0,-0.0293877551d0, &
                   +0.0195918367d0,0.2644897959d0,+0.1175510204d0,-0.0783673469d0,+0.0277551020d0, &
                    0.3746938776d0,+0.1665306122d0,-0.1110204082d0/)

  real(kindreal),dimension(mx),save,private:: xa=(/0.0d0,0.1d0,0.2d0,0.35d0,0.5d0,0.7d0,0.8d0,0.9d0, &
                                                  0.95d0,0.0d0/)
  real(kindreal),dimension(mx),save,private:: dfsx,xx
  real(kindreal),dimension(mz),parameter,private:: za=(/0.0d0,1.0d-4,3.0d-4,1.0d-3,2.0d-3,4.0d-3, &
                                                       1.0d-2,2.0d-2,3.0d-2,4.0d-2,6.0d-2,8.0d-2, &
                                                       1.0d-1/)
  real(kindreal),dimension(mz),save,private:: dfsz,zza
  real(kindreal),dimension(nr),save,private:: alr,dfsr
  real(kindreal),dimension(nrm),save,private:: alrf
  real(kindreal),dimension(nt),save,private:: t6list,alt,dfs
  real(kindreal),dimension(ntm),save,private:: t6listf
  real(kindreal),dimension(ip),save,private:: t6arr
  real(kindreal),dimension(mx,mz),save,private:: x,y,zz
  real(kindreal),dimension(nt,nr),save,private:: opk,opk2,xzf
  real(kindreal),dimension(ntm,ipr),save,private:: f,fx,fy,fxy
  real(kindreal),dimension(ip,nr),save,private:: xzff
  real(kindreal),dimension(mx,nt,nr),save,private:: opl
  real(kindreal),dimension(mx,mz,nt,nr),save,private:: xz
  real(kindreal),dimension(nz2),parameter,private:: &
    zk_kap=(/0.00d0,0.0001d0,0.0003d0,0.001d0,0.002d0,0.004d0,0.01d0,0.02d0,0.03d0,0.04d0,0.5d0,1.0d0/)
  real(kindreal),dimension(nr2),parameter,private:: &
    rk_kap=(/-7.0d0,-6.5d0,-6.0d0,-5.5d0,-5.d0,-4.5d0,-4.0d0,-3.5d0,-3.0d0,-2.5d0,-2.0d0,-1.5d0,-1.0d0,-0.5d0,0.00d0,0.5d0,1.0d0/)
  real(kindreal),dimension(nt2),save,private:: &
    tk_kap=(/0.002d0,0.003d0,0.004d0,0.005d0,0.006d0,0.007d0,0.008d0,0.009d0,0.010d0,0.011d0,0.012d0,0.014d0,0.016d0,0.018d0,&
             0.020d0,0.025d0,0.030d0,0.035d0,0.040d0,0.045d0,0.050d0,0.055d0,0.060d0,0.070d0,0.080d0,0.090d0,0.100d0,0.120d0,&
             0.150d0,0.20d0,0.250d0,0.30d0,0.40d0,0.50d0,0.60d0,0.80d0,1.0d0,1.20d0,1.50d0,2.0d0,2.50d0,3.0d0,4.0d0,5.0d0,6.0d0,&
             8.0d0,10.0d0,15.0d0,20.0d0,30.d0,40.d0,60.d0,80.d0,100.d0,300.d0,500.d0,1000.d0/)

  public:: opacgn93,kappa
  private:: opac,t6rinterp,readco,opaltab,fity,fitx,interp,smooth
  private:: kappa_out,condTest,cond

contains

! ***********************************************************************
!=======================================================================
!> Interpolation of the data along Z
!!  OPACT- opacity obtained from a quadraric interpolation at fixed log T6 at three values of log R;
!!         followed by quadratic interpolation along log T6.
!!         Results smoothed by mixing overlapping quadratics.
!!  DOPACT- is Dlog(k)/Dlog(T6) smoothed by mixing quadratics at fixed R
!!  DOPACR- is  Dlog(k)/Dlog(R) smoothed by mixing quadratics.
!!  DOPACTD- is Dlog(k)/Dlog(T6) smoothed by mixing quadratics at fixed rho
!------------------------------------------------------------------------
  subroutine opacgn93(z,xh,t6,r)
!------------------------------------------------------------------------
  use interpolation,only: quad
  use caramodele,only: nwmd

  implicit none

  save

  real(kindreal),intent(in):: z,xh,t6,r

  integer:: i,iCG,ihi,ilo,imd,is,itime,iw,iz,izi,izz,m2,m3,m4,mfm

  real(kindreal):: dix,dkapz1,dkapz2,dkapz3,dkapz4,kapz1,kapz2,zzl,zzltest
  real(kindreal),dimension(mz):: xkapz,dkapdtr,dkapdrt
!------------------------------------------------------------------------
  error_xzt = 0
  zval=z
  zzl=z   ! use zzl=log10(.0001+z) for log interpolation
  if (itime /= 12345678) then
    do i=1,mz
      zza(i)=za(i) ! use zza=log10(0.0001+za(i)) for log interpolation
    enddo
    itime=12345678
  endif

  do i=1,mz
    if (abs(z-za(i)) < 1.d-7 ) then
      izz=i
      call opac (0,izz,xh,t6,r)
         if (opact > 9.d0) then
           if (verbose) then
             write(*,'(2(a,f7.5),a,f10.5,a,e12.4)') ' logK > 9.0, X=',xh,' Z=',z,' T6=',t6,' R=',r
           endif
         endif
      return
    endif
  enddo

  ilo=2
  ihi=mz
  do while (ihi-ilo > 1)
    imd=(ihi+ilo)/2
    if (z <= za(imd)+1.d-7) then
      ihi=imd
    else
      ilo=imd
    endif
  enddo
  i=ihi
  m1=i-2
  m2=i-1
  m3=i
  m4=i+1
  mfm=m4

! check whether Z is near a table limit
  if ((z <= za(2)+1.d-7) .or. (z >= za(mz-1))) then
    mfm=m3
  endif
!.....  Check if Z+X interpolation sums exceed unity at needed indices.
!       If so, backup to lower Z indices to perform interpolation.
!       This should work OK, due to density of Z-grid points and the slow Z variation(except at very small Z)
  if (xh+za(mfm) > 1.d0) then
    mfm=m3
  endif
  if (xh+za(mfm) > 1.d0) then
    if (m1 <= 1) then
      write(*,*) 'xh,za(mfm),mfm,m3,m2,m1',xh,za(mfm),mfm,m3,m2,m1
      rewind(222)
      write (222,*) nwmd,': X,Z location not covered by opacities'
      stop 'special case: X,Z location not covered by logic'
    endif
    m1=m1-1
    m2=m2-1
    m3=m3-1
  endif

  izi=0
  do iz=m1,mfm
    izz=iz
    call opac(izi,izz,xh,t6,r)
    if (opact > 9.d0 .and. verbose) then
      write(*,'(2(a,f7.5),a,f10.5,a,e12.4)') ' logK > 9.0, X=',xh,' Z=',z,' T6=',t6,' R=',r
    endif
    izi=1
    xkapz(iz)=10.d0**opact ! converts logK to K
    dkapdtr(iz)=dopact
    dkapdrt(iz)=dopacr
  enddo
  is=0
  iw=1
  kapz1=quad(is,iw,zzl,xkapz(m1),xkapz(m2),xkapz(m3),zza(m1),zza(m2),zza(m3),dkap)
  is=1
  dkapz1=quad(is,iw,zzl,dkapdtr(m1),dkapdtr(m2),dkapdtr(m3),zza(m1),zza(m2),zza(m3),dkap)
  dkapz3=quad(is,iw,zzl,dkapdrt(m1),dkapdrt(m2),dkapdrt(m3),zza(m1),zza(m2),zza(m3),dkap)
  if (mfm == m3) then
    if (kapz1 <= 0.d0) then
      zzltest=1.d90
      do iCG = m1,mfm
!mod2010/09/init
        if (abs(zzl-zza(iCG))<=zzltest) then
          zzltest=zzl-zza(iCG)
          kapz1=xkapz(iCG)
        endif
      enddo
    endif
    opact=log10(kapz1)   ! converts K to logK
    dopact=dkapz1
    dopacr=dkapz3
    dopactd=-3.d0*dopacr+dopact
    is=0
    return
  endif
  is=0
  iw=2
  kapz2=quad(is,iw,zzl,xkapz(m2),xkapz(m3),xkapz(m4),zza(m2),zza(m3),zza(m4),dkap)
  is=1
  dkapz2=quad(is,iw,zzl,dkapdtr(m2),dkapdtr(m3),dkapdtr(m4),zza(m2),zza(m3),zza(m4),dkap)
  dkapz4=quad(is,iw,zzl,dkapdrt(m2),dkapdrt(m3),dkapdrt(m4),zza(m2),zza(m3),zza(m4),dkap)
  dix=(zza(m3)-zzl)*dfsz(m3)
  opact=log10(kapz1*dix+kapz2*(1.d0-dix))   ! converts K to logK
  dopact=dkapz1*dix+dkapz2*(1.d0-dix)
  dopacr=dkapz3*dix+dkapz4*(1.d0-dix)
  dopactd=-3.d0*dopacr+dopact
  is=0

  return

  end subroutine opacgn93
!=======================================================================
  subroutine opac (izi,mzin,xh,t6,r)
!..... The purpose of this subroutine is to interpolate log kappa
!      in X, T6, R
!        izi=0 recalulate table indices to use; =1 keep previous
!        mzin=index of za(i) in block data. za(i) are metallicities
!        t6=T6=temperature in millions of degrees kelvin
!        r=R=density(g/cm**3)/T6**3

!..... OPACT- opacity obtained from a quadraric interpolation at
!      fixed log T6 at three values of log R; followed by quadratic
!      interpolation along log T6. Results smoothed bt mixing
!      overlapping quadratics.
!..... DOPACT- is Dlog(k)/Dlog(T6) smoothed by mixing quadratics
!              at fixed R
!..... DOPACR- is  Dlog(k)/Dlog(R) smoothed by mixing quadratics.
!..... DOPACTD- is Dlog(k)/Dlog(T6) smoothed by mixing quadratics
!               at fixed rho
!------------------------------------------------------------------------
  use interpolation,only: quad

  implicit none

  save

  integer, intent(in):: izi,mzin
  real(kindreal), intent(in):: xh,t6,r

  integer::i,iadvance,ihi,ilo,imd,iop,ir,is,it,itime,iw,k1in,k3s,kmin,l3s,mf2,mfin,mg,mh,mi,mxend, &
           ntlimit

  real(kindreal)::dfsxmx1,dixr,ri,slr,slt,t6i,xamx1,xxi,xxmx1,xxx,z
!------------------------------------------------------------------------
  iop=1
! provides smoothed interpolations; iop=0 gives no smoothing

  mzz=mzin
  z=za(mzz)

  if (nre < 6) then
    write(*,*) 'Too few R values; NRE+1-NRB < 6'
    return
  endif

  if ((izi == 0) .and. (z+xh-1.d-6 > 1.d0 )) then
    rewind(222)
    write(222,*)'STOP xztrin'
    stop 'Mass fractions exceed unity'
  endif
  if ((izi /= 0) .and. (zval+xh-1.d-6 > 1.d0 )) then
    rewind(222)
    write(222,*)'STOP xztrin'
    stop 'Mass fractions exceed unity'
  endif
  xxh=xh
  xxi=xh
  t6i=t6
  ri=r

  xxx=log10(.005d0+xh)
  slt=log10(t6)
  slr=log10(r)

  if (itime /= 12345678) then
    itime=12345678
    do  i=1,mx
      xx(i)=log10(.005d0+xa(i))
    enddo
!..... this is the first time throught. Calculate the decadic
!      log of the perimeter points shifted by Z. m refers to
!      xa(m); the hydrogen table value.

!..... read the data files
    call readco
    xamx1=xa(mx-1)
    xxmx1=xx(mx-1)
    dfsxmx1=dfsx(mx-1)
  endif
  mxend=mx
  xa(mx)=1.d0-z
  xa(mx-1)=xamx1
  xx(mx-1)=xxmx1
  dfsx(mx-1)=dfsxmx1
  if (xa(mx) < xa(mx-1)) then
    mxend=mx-1
    xa(mxend)=xa(mx)
  endif
  if (xh >= 0.8d0 ) then
    xx(mxend)=log10 (0.005d0+xa(mxend))
    dfsx(mxend)=1.d0/(xx(mxend)-xx(mxend-1))
  endif

!..... Determine log R and log T6 grid points to use in the
!      interpolation.
  if ((slt < alt(1)) .or. (slt > alt(nt))) then
    if(verbose) then
      write(*,*) ' T6/LogR outside of table range'
    endif
    error_xzt = 1
    return
  endif
  if ((slr < alr (1)) .or. (slr > alr(nre))) then
    if (verbose) then
      write(*,*) ' T6/LogR outside of table range'
    endif
    error_xzt = 1
    return
  endif

  if (izi == 0) then  ! freeze table indices
    ilo=2
    ihi=mx
    do while (ihi-ilo > 1)
      imd=(ihi+ilo)/2
      if (xh <= xa(imd)+1.d-7) then
        ihi=imd
      else
        ilo=imd
      endif
    enddo
    i=ihi
    mf=i-2
    mg=i-1
    mh=i
    mi=i+1
    mf2=mi
    if (xh < 1.d-6) then
      mh=1
      mg=1
      mi=2
      mf2=1
    endif
    if ((xh <= xa(2)+1.d-7) .or. (xh >= xa(mx-2)-1.d-7)) then
      mf2=mh
    endif

    ilo=2
    ihi=nre
    do while (ihi-ilo > 1)
      imd=(ihi+ilo)/2
      if (slr <= alr(imd)+1.d-7) then
        ihi=imd
      else
        ilo=imd
      endif
    enddo
    i=ihi
    l1=i-2
    l2=i-1
    l3=i
    l4=l3+1

    ilo=2
    ihi=nt
    do while (ihi-ilo > 1)
      imd=(ihi+ilo)/2
      if (t6 <= t6list(imd)+1.d-7) then
        ihi=imd
      else
        ilo=imd
      endif
    enddo
    i=ihi
    k1=i-2
    k2=i-1
    k3=i
    k4=k3+1
    l3s=l3+nrb-1
    k3s=k3+ntb-1
  endif

  kmin=0
  k1in=k1
  iadvance=0
  mfin=mf
  if ((mfin == 1) .and. (xz(1,mzz,k1,l1) > 9.d0)) then
! no data
    do i=1,6
      if (xz(1,mzz,i,l1) > 9.d0)  then
        if (xh < 0.1d0) then
          kmin=i+1
        else
          if (iadvance == 0) then
            iadvance=iadvance+1
            mf=mf+1
            mg=mg+1
            mh=mh+1
            mi=mi+1
            mf2=mf2+1
          endif
        endif
      endif
    enddo
    if ((iadvance == 0) .and. (k1 <= kmin) .and. (slt <= alt(kmin))) then
      k1=kmin
      if ((xz(1,mzz,kmin,l1+1) < 9.d0) .and. ((slr+0.01d0) > alr(l1+1))) then
        l1=l1+1
        kmin=0
        k1=k1in
        do i=1,6
          if (xz(1,mzz,i,l1) > 9.d0) then
            kmin=i+1
          endif
        enddo
        if ((kmin /= 0) .and. (k1in < kmin)) then
          k1=kmin
        endif
      endif
    endif
    if ((slt+0.001d0) < alt(k1)) then
      opact=30.d0
      dopact=99.d0
      dopacr=99.d0
      return
    endif
    l2=l1+1
    l3=l2+1
    l4=l3+1
    l3s=l3+nrb-1
    k2=k1+1
    k3=k2+1
    k4=k3+1
    k3s=k3+ntb-1
  endif

  do i=14,18   ! allows jagged edge at high T,rho
    if ((l3s > i) .and. (k3s > nta(i+1))) then
      if (verbose) then
        write(*,*) ' T6/LogR outside of table range'
      endif
      error_xzt = 1
      return
    endif
  enddo
  do m=mf,mf2
    ip_op=3
    iq_op=3
    ntlimit=nta(l3s)
    if ((k3 == ntlimit) .or. (iop == 0)) then
      ip_op=2
    endif
    if(t6 <= t6list(2)+1.d-7) then
      ip_op=2
    endif

    if ((l3 == nre) .or. (iop == 0)) then
      iq_op=2
    endif
    if ((l4 <= nr) .and. (xz(m,mzz,k3,l4) == 0.d0)) then
      iq_op=2
    endif
    if (slr <= alr(2)+1.d-7) then
      iq_op=2
    endif

    is=0

    do ir=l1,l1+iq_op
      do it=k1,k1+ip_op
        opl(m,it,ir)=xz(m,mzz,it,ir)
        is=1
      enddo
    enddo
  enddo
  if ((zz(mg,mzin) /= zz(mf,mzin)) .or. (zz(mh,mzin) /= zz(mf,mzin))) then
    rewind(222)
    write(222,*)'STOP xztrin'
    stop 'Z does not match Z in GN93hz files you are using'
  endif
  if (z /= zz(mf,mzin)) then
    rewind(222)
    write(222,*)'STOP xztrin'
    write(*,*) 'z, zz= ', z, zz
    stop 'Z does not match Z in codata* files you are using'
  endif

  is=0
  iw=1
  do ir=l1,l1+iq_op
    do it=k1,k1+ip_op
      if (mf2 == 1) then
        opk(it,ir)=opl(mf,it,ir)
      else
        opk(it,ir)=quad(is,iw,xxx,opl(mf,it,ir),opl(mg,it,ir),opl(mh,it,ir),xx(mf),xx(mg),xx(mh),dkap)
        is=1
      endif
    enddo
  enddo

  if (mi == mf2) then  ! interpolate between quadratics
    is=0
    iw=1
    dixr=(xx(mh)-xxx)*dfsx(mh)
    do ir=l1,l1+iq_op
      do it=k1,k1+ip_op
        opk2(it,ir)=quad(is,iw,xxx,opl(mg,it,ir),opl(mh,it,ir),opl(mi,it,ir),xx(mg),xx(mh),xx(mi),dkap)
        opk(it,ir)=opk(it,ir)*dixr+opk2(it,ir)*(1.d0-dixr)
        is=1
      enddo
    enddo
!     interpolate X between two overlapping quadratics
  endif

  is=0

!..... completed H,Z interpolation. Now interpolate T6 and log R on a
!      4x4 grid. (log(T6(i)),i=i1,i1+3),log(R(j)),j=j1,j1+3)).Procedure
!      mixes overlapping quadratics to obtain smoothed derivatives.

  call t6rinterp(slr,slt)

  return

  end subroutine opac
!=======================================================================
  subroutine t6rinterp(slr,slt)
!     The purpose of this subroutine is to interpolate in logT6 and logR
!------------------------------------------------------------------------
  use interpolation,only: quad

  implicit none

  save

  integer:: iu,is,iw,kx,lx

  real(kindreal), intent(in):: slr,slt

  real(kindreal)::dix,dix2,dkap1,dkap2,dkapq1,dkapq2,dkapq3,dopacrq,dopactq,opacr,opacrq,opact2,opactq, &
                opactq2
!------------------------------------------------------------------------
  iu=0
  is=0
  do kx=k1,k1+ip_op
    iw=1
    iu=iu+1
    hxz(iu)=quad(is,iw,slr,opk(kx,l1),opk(kx,l2),opk(kx,l3),alr(l1),alr(l2),alr(l3),dkap)
    if(iq_op == 3) then
      iw=2
      qxz(iu)=quad(is,iw,slr,opk(kx,l2),opk(kx,l3),opk(kx,l4),alr(l2),alr(l3),alr(l4),dkap)
    endif
    is=1
  enddo

  is=0
  iw=1
!..... k and Dlog(k)/dlog(T6) in lower-right 3x3(i=i1,i1+2 j=j1,j1+2)
  opact=quad(is,iw,slt,hxz(1),hxz(2),hxz(3),alt(k1),alt(k2),alt(k3),dkap)
  dopact=dkap
  dkap1=dkap
  if (iq_op == 3) then
!.....    k and Dlog(k)/Dlog(T6) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2)
    opactq=quad(is,iw,slt,qxz(1),qxz(2),qxz(3),alt(k1),alt(k2),alt(k3),dkap)
    dkapq1=dkap
  endif
  if (ip_op == 3) then
!.....    k and Dlog(k)/Dlog(T6) in lower-left 3x3.
    opact2=quad(is,iw,slt,hxz(2),hxz(3),hxz(4),alt(k2),alt(k3),alt(k4),dkap)
    dkap2=dkap
!.....    k and Dlog(k)/Dlog(T6) smoothed in left 3x4
    dix=(alt(k3)-slt)*dfs(k3)
    dopact=dkap1*dix+dkap2*(1.d0-dix)
    opact=opact*dix+opact2*(1.d0-dix)
  endif
  if (iq_op == 3) then
!.....    k and Dlog(k)/Dlog(T6) in upper-right 3x3.
    opactq2=quad(is,iw,slt,qxz(2),qxz(3),qxz(4),alt(k2),alt(k3),alt(k4),dkap)
    dkapq2=dkap
    dopactq=dkapq1*dix+dkapq2*(1.d0-dix)
    opactq=opactq*dix+opactq2*(1.d0-dix)
  endif

  iu=0
  do lx=l1,l1+iq_op
    iw=1
    iu=iu+1
    hxz(iu)=quad(is,iw,slt,opk(k1,lx),opk(k2,lx),opk(k3,lx),alt(k1),alt(k2),alt(k3),dkap)
    if (ip_op == 3) then
      iw=2
      qxz(iu)=quad(is,iw,slt,opk(k2,lx),opk(k3,lx),opk(k4,lx),alt(k2),alt(k3),alt(k4),dkap)
    endif
    is=1
  enddo

  is=0
  iw=1
!..... k and Dlog(k)/Dlog(R) in lower-left 3x3
  opacr=quad(is,iw,slr,hxz(1),hxz(2),hxz(3),alr(l1),alr(l2),alr(l3),dkap)
  dopacr=dkap
  if (ip_op == 3) then
    opacrq=quad(is,iw,slr,qxz(1),qxz(2),qxz(3),alr(l1),alr(l2),alr(l3),dkap)
!.....    k and Dlog(k)/Dlog(R) in upper-left 3x3.
    dkapq3=dkap
  endif
  if (iq_op == 3) then
!.....    k and Dlog(k)/Dlog(R) in lower-right 3x3.
    opact2=quad(is,iw,slr,hxz(2),hxz(3),hxz(4),alr(l2),alr(l3),alr(l4),dkap)
    dix2=(alr(l3)-slr)*dfsr(l3)
    dopacr=dopacr*dix2+dkap*(1.d0-dix2)
    if (ip_op == 3) then
!.....        k and Dlog(k)/Dlog(T6) smoothed in both log(T6) and log(R)
      dopact=dopact*dix2+dopactq*(1.d0-dix2)
      opact=opact*dix2+opactq*(1.d0-dix2)
    endif
  endif
  if (ip_op == 3) then
!.....    k and Dlog(k)/Dlog(R) in upper-right 3x3.
    opacrq=quad(is,iw,slr,qxz(2),qxz(3),qxz(4),alr(l2),alr(l3),alr(l4),dkap)
    if (iq_op == 3) then
!.....        Dlog(k)/Dlog(R) smoothed in both log(T6) and Log(R).
      dopacrq=dkapq3*dix2+dkap*(1.d0-dix2)
      dopacr=dopacr*dix+dopacrq*(1.d0-dix)
    endif
  endif
  dopactd=dopact-3.d0*dopacr
  if (opact > 1.d+15) then
    rewind(222)
    write(222,*)'STOP xztrin'
    stop 'Interpolation indices out of range'
  endif
  if (opact > 9.d0) then
    dopact=99.d0
    dopacr=99.d0
    dopactd=99.d0
  endif

  return

  end subroutine t6rinterp
!=======================================================================
  subroutine readco
! The purpose of this subroutine is to read the data tables
!------------------------------------------------------------------------
  use inputparam,only: iopac

  implicit none

  save

  integer,parameter::ismdata=1   ! modified
  integer:: error_readco,i,j,k,l,isett6,itimeco,kk,ll

  character(256):: opacfile
!------------------------------------------------------------------------
  if (itimeco /= 12345678) then
    do i=1,mx
      do j=1,mz
        do k=1,nt
          do l=1,nr
            xz(i,j,k,l)=1.d+35
          enddo
        enddo
      enddo
    enddo
    itimeco=12345678
  endif

  close (2)
!..... read  tables
!         open(2, FILE='GN93hz')
!         new file name + err option
  select case(iopac)
    case (1)
      opacfile=trim(input_dir)//'inputs/opaSol_GN93.dat'
    case (2)
      opacfile=trim(input_dir)//'inputs/opaAlph_GN93.dat'
    case (3)
      opacfile=trim(input_dir)//'inputs/opaSol_AspCun06.dat'
    case (4)
      opacfile=trim(input_dir)//'inputs/opaAlph_AspCun06.dat'
    case default
      stop 'Bad IOPAC choice ! Must be 1,2,3 or 4.'
  end select

  open(2, file=opacfile,iostat=error_readco)
! old goto 1234
  if (error_readco /= 0) then
    write(*,*) 'File ',trim(opacfile),' not found !!'
    rewind(222)
    write(222,*) 'File ',trim(opacfile),' not found !!'
    stop
  endif

  do i=1,240
    read (2,*)
  enddo

  do m=1,mx
    do i=1,n(m)
      read(2,*)
      read(2,'(7x,i3,26x,f6.4,3x,f6.4,3x,f6.4)') itab(m,i),x(m,i),y(m,i),zz(m,i)
      read(2,*)
      read(2,*)
      read(2,*)
      read(2,'(4x,f6.1,18f7.1)') (alrf(kk),kk=1,nrm)
      read(2,*)
      do k=1,ntm
        read(2,'(f4.2,19f7.3)') alt(k),(xzf(k,l), l=1,nrm)
        alt(k)=alt(k)-6.d0
        if (isett6 /= 1234567) then
          t6listf(k)=10.d0**alt(k)
          t6arr(k)=t6listf(k)
        endif
!        do ll=1,nrm   ! modified
          xzff(k,1:nrm)=xzf(k,1:nrm)
!        enddo
      enddo
      isett6=1234567

      if (ismdata == 0) then
        tmax=10.d0   ! modified
! nset does NOT need to be modified ! since nset counts the some number of lines with temperatures
! ABOVE 3.75
        nset=65
        rls=-8.d0
        nsm=1
        rle=1.d0
        nrlow=1
        nrhigh=2*int(rle-rls)+1
        call opaltab    !modified
      endif

      ll=1
      do kk=1,nre
        alr(ll)=alrf(kk)
        do k=1,nt
          t6list(k)=t6listf(k+ntb-1)
          if(ismdata == 0) then
! Following skip required because, due to missing data, the X=0  low T data cannot be smoothed
            if ((m  == 1) .and. (k <= 9)) then
              xz(m,i,k,ll)=xzf(k+ntb-1,kk)
            else
              xz(m,i,k,ll)=xzff(k+ntb-1,kk)
            endif
          else
            xz(m,i,k,ll)=xzf(k+ntb-1,kk)
          endif
        enddo
        ll=ll+1
      enddo
    enddo
  enddo

  do i=2,nt
    dfs(i)=1.d0/(alt(i)-alt(i-1))
  enddo
  do i=2,nr
    dfsr(i)=1.d0/(alr(i)-alr(i-1))
  enddo
  do i=2,mx-1
    dfsx(i)=1.d0/(xx(i)-xx(i-1))
  enddo
  do i=2,mz
    dfsz(i)=1.d0/(zza(i)-zza(i-1))
  enddo

  return

  end subroutine readco
!=======================================================================
  subroutine opaltab
!------------------------------------------------------------------------
!  CODE FOR FITTING AND SMOOTHING OPAL DATA. ADAPTED FROM A CODE
!     WRITTEN BY MIKE SEATON(obtained june 1993)

!     OPAL DATA.
!     ASSUMES FIRST T6=0.006, LAST T6=10.OR 0.04). Depending on position
!     in the table.
!     USES RECTANGULAR ARRAY FOR VARIABLES T6 AND LOG10(R)

!     (1) NSM=NUMBER OF PASSES THROUGH SMOOTHING FILTER.
!     USE OF NSM=1 OR 2 IS RECOMMENDED.
!     NO SMOOTHING WITH NSM=0
!     (2) RANGE FOR LOG10(R),
!     RLS=FIRST VALUE, RLE=LAST VALE
!     (RLS MUST BE FIRST VALUYE IN TABLE)

!  SUBROUTINE INTERP
!     AFTER PROCESSING, DATA ARE IN A FORM FOR USE OF
!               SUBROUTINE INTERP
!     WHICH GIVES LOG(ROSS) AND TWO FIRST DERIVATIVES FOR ANY
!     VALUES OF LOG(T) AND LOG(RHO). SEE BELOW FOR FURTHER
!     EXPLANATION.

!  OUTPUT FOR THE CASE OF NSM.GT.0.
!     INTERP IS USED TO OBTAIN SMOOTHED DATA INTERPOLATED
!     BACK TO THE ORIGINAL OPAL MESH. TWO FILES ARE WRITTEN.

!  THE SUBROUTINES SPLINE AND SPLINT ARE ADAPTED FROM THOSE GIVE BY
!  W.H. Press, S.A. Teulolsky, W.T. Vettering and B.P. Flannery,
!  "Numerical Recipes in FORTRAN", 2nd edn., 1992, C.U.P.
!  OTHER REFERENCES ARE MADE TO METHODS DESCRIBED IN THAT BOOK.
!------------------------------------------------------------------------
  use interpolation,only: spline,splint

  implicit none

  integer,parameter:: naf=15
  integer:: i,j,k,l,ns

  real(kindreal):: dgdrho,dgdt,flr,flrho,flt,g_int,t6
  real(kindreal),dimension(ip):: u,v,v2
  real(kindreal),dimension(ip,ipr):: rossl

  logical::ierr
!------------------------------------------------------------------------
  nrl=2*int(rle-rls)+1

!     STORE LOG10(T) IN U AND LOG10(ROSS) IN ROSSL
!     CHECK FIRST VALUE OF T6
  t6=t6arr(1)
  do j=1,nrl
    rossl(1,j)=xzff(1,j)
  enddo

  if (abs(t6-10.d0**3.75d0/1.d6) < 1.d-8) then
    u(1)=6.d0+log10(t6)
  endif
  i=1
  do
    i=i+1
    t6=t6arr(i)
    do j=1,nrl
      rossl(i,j)=xzff(i,j)
    enddo
    u(i)=6.d0+log10(t6)
    if(t6 >= tmax) then
      exit
    endif
  enddo
  ntemp=i
  if(ntemp > ip) then
    rewind(222)
    write(222,*)'STOP xztrin'
    write(*,*) ' REQUIRE PARAMETER IP OF AT LEAST ',ntemp
    stop
  endif

  do J=1,nrl
! FOR EACH LOG10(R), STORE LOG10(ROSS) IN V(I)
    do i=1,ntemp
      v(i)=rossl(i,j)
    enddo

! GET FIRST DERIVATIVES AT END POINTS

! GET SECOND DERIVATIVES FOR SPLINE FIT
    call spline(u,v,ntemp,v2)

! INTERPOLATE TO LOG10(T)=FLT, FLT=3.8(0.05)8.0
    do i=1,nset ! modified
      flt=3.75d0+0.05d0*i
      call splint(u,v,ntemp,v2,flt,f(i,j),fx(i,j))
    enddo
  enddo

!  OPTION FOR SMOOTHING
  if(nsm > 0) then
    do ns=1,nsm
      call smooth
    enddo
    call fitx
  endif

!  GET FY AND FXY
  call fity

!  THE ARRAYS F, FX, FY AND FXY ARE NOW STORED

! INTERPOLATE BACK TO OPAL POINTS
  if(nsm > 0) then
    do l=1,nrl
      xzff(1,l)=rossl(1,l)
    enddo

    do k=2,ntemp
      flt=u(k)
      do l=nrlow,nrhigh
        flr=rls+0.5d0*(l-1)
        flrho=flr-18.d0+3.d0*flt
        call interp(flt,flrho,g_int,dgdt,dgdrho,ierr)
        v(l)=g_int
      enddo
      t6=t6arr(k)
      do l=nrlow,nrhigh
        xzff(k,l)=v(l)
      enddo
    enddo
  endif

  return

  end subroutine opaltab
!=======================================================================
  subroutine fity
!  THIS ROUTINE MAKES SPLINE FITS FOR F AND FX, AND OBTAINS
!  FY AND FXY
!------------------------------------------------------------------------
  use interpolation,only: getd

  implicit none

  integer::i,j

  real(kindreal):: ap1,apn,bp1,bpn
  real(kindreal),dimension(ipr):: a_fit,b_fit,ad,bd
!------------------------------------------------------------------------
  do i=1,nset   ! modified
    do j=1,nrl
      a_fit(j)=f(i,j)
      b_fit(j)=fx(i,j)
    enddo

    call getd(a_fit,nrl,ad,ap1,apn)
    call getd(b_fit,nrl,bd,bp1,bpn)

    fy(i,1)=ap1
    fy(i,nrl)=apn
    fxy(i,1)=bp1
    fxy(i,nrl)=bpn
    do j=2,nrl-1
      fy(i,j)= -a_fit(j)+a_fit(j+1)-2.d0*ad(j)-ad(j+1)
      fxy(i,j)=-b_fit(j)+b_fit(j+1)-2.d0*bd(j)-bd(j+1)
    enddo
  enddo

  return

  end subroutine fity
!=======================================================================
  subroutine fitx
!  THIS ROUTINE IS USED ONLY AFTER SMOOTHING.
!  ITS FUNCTION IS TO RECOMPUTE FX USING SMOOTHED F.
!------------------------------------------------------------------------
  use interpolation,only: getd

  implicit none

  integer:: i,j

  real(kindreal):: ap1,apn
  real(kindreal),dimension(ntm):: a_fit,d_fit
!------------------------------------------------------------------------
  do j=1,nrl
    do i=1,nset ! modified
      a_fit(i)=f(i,j)
    enddo
    call getd(a_fit,nset,d_fit,ap1,apn)  ! modified
    fx(1,j)=ap1
    fx(nset,j)=apn   ! modified
    do i=2,nset-1  ! modified
      fx(i,j)=-a_fit(i)+a_fit(i+1)-2.d0*d_fit(i)-d_fit(i+1)
    enddo
  enddo

  return

  end subroutine fitx
!=======================================================================
  subroutine interp(flt,flrho,g_int,dgdt,dgdrho,ierr)
!  GIVEN F,FX,FY AND FXY ON THE GRID POINTS, THIS ROUTINE
!  DOES BI-CUBIC INTERPOLATIONS USING METHODS DESCRIBED IN
!  Numerical Recipes, PP. 118 TO 120
!------------------------------------------------------------------------
  implicit none

  integer:: i,j

  real(kindreal), intent(in):: flt,flrho
  real(kindreal), intent(out):: g_int,dgdt,dgdrho

  real(kindreal):: ff_interp,ffx_interp,ffy_interp,flr,u,v,x_int,y_int
  real(kindreal),dimension(16):: b_int

  logical, intent(out):: ierr
!------------------------------------------------------------------------
!  EXTREME LIMITS ALLOWED ARE:-
!     (3.800-0.0125) TO (8.000+0.0125) FOR LOG10(T)
!     (RLS-0.125) TO (RLE+0.1254) FOR LOG10(R)
!     (ALLOWING FOR SMALL EXTRAPOLATIONS BEYOND TABULAR VALUES)

  ierr=.false.
  i=0

  x_int=20.d0*(flt-3.800d0)+1.d0
  flr=flrho+18.d0-3.d0*flt
  y_int=2.d0*( flr - rls )+1.d0

  if (x_int < 2.d0) then
    if (x_int < 0.75d0) then
      ierr=.true.
    else
      i=1
    endif
  elseif (x_int > 84.d0) then
    if (x_int > 85.25d0) then
      ierr=.true.
    else
      i=84
    endif
  else
    i=int(x_int)
  endif
  u=x_int-i

  if (y_int < 2.d0) then
    if (y_int < 0.75d0) then
      ierr=.true.
    else
      j=1
    endif
  elseif (y_int > nrl-1.d0) then
    if (y_int > nrl+0.25d0) then
      ierr=.true.
    else
      j=nrl-1
    endif
  else
    j=int(y_int)
  endif
  v=y_int-j

  if (ierr) then
    g_int=9.999d0
    dgdt=9.999d0
    dgdrho=9.999d0
    return
  endif

!  GIVEN FUNCTIONS AND DERIVATIVES AT GRID POINTS, COMPUTE COEFFICIENTS.
  b_int(1)=f(i,j)
  b_int(2)=fy(i,j)
  b_int(3)=3.d0*(-f(i,j)+f(i,j+1))-2.d0*fy(i,j)-fy(i,j+1)
  b_int(4)=2.d0*(f(i,j)-f(i,j+1))+fy(i,j)+fy(i,j+1)

  b_int(5)=fx(i,j)
  b_int(6)=fxy(i,j)
  b_int(7)=3.d0*(-fx(i,j)+fx(i,j+1))-2.d0*fxy(i,j)-fxy(i,j+1)
  b_int(8)=2.d0*(fx(i,j)-fx(i,j+1))+fxy(i,j)+fxy(i,j+1)

  b_int(9)=3.d0*(-f(i,j)+f(i+1,j))-2.d0*fx(i,j)-fx(i+1,j)
  b_int(10)=3.d0*(-fy(i,j)+fy(i+1,j))-2.d0*fxy(i,j)-fxy(i+1,j)
  b_int(11)=9.d0*(f(i,j)-f(i+1,j)+f(i+1,j+1)-f(i,j+1))+6.d0*(fx(i,j)-fx(i,j+1)+fy(i,j)-fy(i+1,j)) + &
            4.d0*fxy(i,j)+3.d0*(fx(i+1,j)-fx(i+1,j+1)-fy(i+1,j+1)+fy(i,j+1)) + &
            2.d0*(fxy(i,j+1)+fxy(i+1,j))+fxy(i+1,j+1)
  b_int(12)=6.d0*(-f(i,j)+f(i+1,j)-f(i+1,j+1)+f(i,j+1))+4.d0*(-fx(i,j)+fx(i,j+1)) + &
            3.d0*(-fy(i,j)+fy(i+1,j)+fy(i+1,j+1)-fy(i,j+1)) + &
            2.d0*(-fx(i+1,j)+fx(i+1,j+1)-fxy(i,j)-fxy(i,j+1))-fxy(i+1,j)-fxy(i+1,j+1)
  b_int(13)=2.d0*(f(i,j)-f(i+1,j))+fx(i,j)+fx(i+1,j)
  b_int(14)=2.d0*(fy(i,j)-fy(i+1,j))+fxy(i,j)+fxy(i+1,j)
  b_int(15)=6.d0*(-f(i,j)+f(i+1,j)-f(i+1,j+1)+f(i,j+1))+4.d0*(-fy(i,j)+fy(i+1,j)) + &
            3.d0*(-fx(i,j)-fx(i+1,j)+fx(i+1,j+1)+fx(i,j+1)) + &
            2.d0*(fy(i+1,j+1)-fy(i,j+1)-fxy(i,j)-fxy(i+1,j))-fxy(i+1,j+1)-fxy(i,j+1)
  b_int(16)=4.d0*(f(i,j)-f(i+1,j)+f(i+1,j+1)-f(i,j+1))+2.d0*(fx(i,j)+fx(i+1,j)-fx(i+1,j+1)-fx(i,j+1) + &
            fy(i,j)-fy(i+1,j)-fy(i+1,j+1)+fy(i,j+1))+fxy(i,j)+fxy(i+1,j)+fxy(i+1,j+1)+fxy(i,j+1)

!  GET G=LOG10(ROSS), DGDT=d LOG10(ROSS)/d LOG10(T),
!      DGDRHO=d LOG10(ROSS)/d LOG10(RHO)

  ff_interp = b_int( 1)+v*(b_int( 2)+v*(b_int( 3)+v*b_int( 4))) + &
              u*( b_int( 5)+v*(b_int( 6)+v*(b_int( 7)+v*b_int( 8))) + &
              u*( b_int( 9)+v*(b_int(10)+v*(b_int(11)+v*b_int(12))) + &
              u*( b_int(13)+v*(b_int(14)+v*(b_int(15)+v*b_int(16))))))

  ffx_interp = b_int( 5)+v*(b_int( 6)+v*(b_int( 7)+v*b_int( 8))) + &
               u*(  2*(b_int( 9)+v*(b_int(10)+v*(b_int(11)+v*b_int(12)))) + &
               u*(  3*(b_int(13)+v*(b_int(14)+v*(b_int(15)+v*b_int(16))))))

  ffy_interp = b_int( 2)+u*(b_int( 6)+u*(b_int(10)+u*b_int(14))) + &
               v*(  2*(b_int( 3)+u*(b_int( 7)+u*(b_int(11)+u*b_int(15)))) + &
               v*(  3*(b_int( 4)+u*(b_int( 8)+u*(b_int(12)+u*b_int(16)))) ))


  g_int=ff_interp
  dgdt=20.d0*ffx_interp-6.d0*ffy_interp
  dgdrho=2.d0*ffy_interp


  return

  end subroutine interp
!=======================================================================
  subroutine smooth
!  THIS SUBROUTINE USES A 2-DIMENSIONAL GENERALISATION OF THE SMOOTHING
!  TECHNIQUES DESCRIBED ON PP. 644 TO 649 OF Numerical Recipes.

!  CONSIDER THE 25 POINTS DEFINED BY
!       I+n, n=-2,-1,0,1,2 AND J+m, m=-2,-1,0,1,2.
!  THE FUNCTION TO BE SMOOTHED IS FITTED TO A BI-CUBIC, INVOLVING
!  16 COEFFICIENTS, USING TECHNIQUES OF LEAST-SQUARES. THE SMOOTHED
!  FUNCTION (TEMPORARILY STORED IN FXY) IS GIVEN BY THE FITTED VALUE
!  AT THE POINT I AND J.

!  THE FITTING IS SHIFTED FOR POINTS CLOSE TO BOUNDARIES.
!------------------------------------------------------------------------
  implicit none

  integer:: i,j
!------------------------------------------------------------------------
  do i=3,nset-2
    j=1
    fxy(i,j)=alp(1)*( f(i-2,j  )+f(i+2,j  ) )+alp(2)*( f(i-2,j+1)+f(i+2,j+1)+f(i-2,j+3)+f(i+2,j+3) + &
             f(i-1,j+4)+f(i+1,j+4) )+alp(3)*( f(i-2,j+2)+f(i+2,j+2) )+alp(4)*( f(i-2,j+4)+f(i+2,j+4) ) + &
             alp(5)*( f(i-1,j  )+f(i+1,j  ) )+alp(6)*( f(i-1,j+1)+f(i+1,j+1)+f(i-1,j+3)+f(i+1,j+3) ) + &
             alp(7)*( f(i-1,j+2)+f(i+1,j+2) )+alp(8)*  f(i  ,j  )+alp(9)*( f(i  ,j+1)+f(i  ,j+3) ) + &
             alp(10)* f(i  ,j+2) +alp(11)*f(i  ,j+4)

    j=2
    fxy(i,j)=bet(1)*( f(i-2,j-1)+f(i+2,j-1)+f(i-2,j+3)+f(i+2,j+3) )+bet(2)*( f(i-2,j  )+f(i+2,j  ) ) + &
             bet(3)*( f(i-2,j+1)+f(i+2,j+1) )+bet(4)*( f(i-2,j+2)+f(i+2,j+2)+f(i-1,j-1)+f(i+1,j-1) + &
             f(i-1,j+3)+f(i+1,j+3) )+bet(5)*( f(i-1,j  )+f(i+1,j  ) )+bet(6)*( f(i-1,j+1)+f(i+1,j+1) ) + &
             bet(7)*( f(i-1,j+2)+f(i+1,j+2) )+bet(8)*( f(i  ,j-1)+f(i  ,j+3) ) + &
             bet(9)*f(i  ,j  ) +bet(10)*f(i  ,j+1) +bet(11)*f(i  ,j+2)

    do j=3,nrl-2
      fxy(i,j)=gam(1)*( f(i-2,j-2)+f(i-2,j+2)+f(i+2,j-2)+f(i+2,j+2) ) + &
               gam(2)*( f(i-2,j+1)+f(i-2,j-1)+f(i-1,j-2)+f(i-1,j+2) + &
               f(i+1,j-2)+f(i+1,j+2)+f(i+2,j-1)+f(i+2,j+1) ) + &
               gam(3)*( f(i-2,j  )+f(i+2,j  )+f(i  ,j-2)+f(i  ,j+2) ) + &
               gam(4)*( f(i-1,j-1)+f(i-1,j+1)+f(i+1,j-1)+f(i+1,j+1) ) + &
               gam(5)*( f(i-1,j  )+f(i  ,j-1)+f(i  ,j+1)+f(i+1,j  ) ) + &
               gam(6)*  f(i  ,j  )
    enddo

    j=nrl-1
    fxy(i,j)=bet(1)*( f(i-2,j+1)+f(i+2,j+1)+f(i-2,j-3)+f(i+2,j-3) )+bet(2)*( f(i-2,j  )+f(i+2,j  ) ) + &
             bet(3)*( f(i-2,j-1)+f(i+2,j-1) )+bet(4)*( f(i-2,j-2)+f(i+2,j-2)+f(i-1,j+1)+f(i+1,j+1) + &
             f(i-1,j-3)+f(i+1,j-3) )+bet(5)*( f(i-1,j  )+f(i+1,j  ) )+bet(6)*( f(i-1,j-1)+f(i+1,j-1) ) + &
             bet(7)*( f(i-1,j-2)+f(i+1,j-2) )+bet(8)*( f(i  ,j+1)+f(i  ,j-3) )+bet(9)*f(i  ,j  ) + &
             bet(10)*f(i  ,j-1) +bet(11)*f(i  ,j-2)

    j=nrl
    fxy(i,j)=alp(1)*( f(i-2,j  )+f(i+2,j  ) )+alp(2)*( f(i-2,j-1)+f(i+2,j-1)+f(i-2,j-3)+f(i+2,j-3) + &
             f(i-1,j-4)+f(i+1,j-4) )+alp(3)*( f(i-2,j-2)+f(i+2,j-2) )+alp(4)*( f(i-2,j-4)+f(i+2,j-4) ) + &
             alp(5)*( f(i-1,j  )+f(i+1,j  ) )+alp(6)*( f(i-1,j-1)+f(i+1,j-1)+f(i-1,j-3)+f(i+1,j-3) ) + &
             alp(7)*( f(i-1,j-2)+f(i+1,j-2) )+alp(8)*  f(i  ,j  )+alp(9)*( f(i  ,j-1)+f(i  ,j-3) ) + &
             alp(10)* f(i  ,j-2) +alp(11)*f(i  ,j-4)

  enddo

  do i=3,nset-2   ! modified
    do j=1,nrl
      f(i,j)=fxy(i,j)
    enddo
  enddo

  return

  end subroutine smooth

!=======================================================================
subroutine kappa_out(rh,t,rhp,rht,x_kap,y_kap,cap,capp,capt,jj1)
!------------------------------------------------------------------------
  use const,only: um
  use caramodele,only: nwmd
  use inputparam,only: ikappa,ibasnet,ioutable,rout,tout,ialflu
  use abundmod,only: abundCheck
  use interpolation, only: indic,flin,qua,quad_gg

  implicit none

  real(kindreal),intent(in):: rh,t,rhp,rht,x_kap,y_kap
  integer,intent(inout)::jj1
  real(kindreal),intent(inout):: cap,capp,capt

  integer,save:: lec=0
  integer:: i,k,minz,l,m=0,j,icase,icase2,icase3,icase4,jt,jr,jz,irmax,jtmin,ixmin,ixmax,izmax,iz,izz,&
            ixx,ir=0,it=0,ixxx
  real(kindreal):: z_kap,t6,r,captt,caprr,tmin,zkm,rkm,tkm,r1,r2,r3,t1,t2,t3,at,at1,c11,c12,c13,&
                   c21,c22,c23,c31,c32,c33,frt1,frt2,frt3,ftr1,ftr2,ftr3,cap10,t6_table

  real(kindreal),dimension(nr2,nt2,nx2,nz2),save:: opa=0.d0
  real(kindreal),dimension(3):: b11_kap,b12,b13,b21,b22,b23,b31,b32,b33_kap
  real(kindreal),dimension(3,3):: a11,a12,a13,a21,a22,a23,a31,a32,a33
!------------------------------------------------------------------------
! CAVEAT: pour le moment, l'extrapolation en cas de sortie des tables OPAL se fait sur des anciennes tables
!         avec abondances de Grevesse & Noels 1993...
!!! A MODIFIER DES QUE POSSIBLE !!!
! lecture des tables
  if (lec /= 1) then
    lec = 1
    ioutable = 0
    open(file=input_dir//'inputs/kappa93.dat',unit=22,status='old')
    do i=1,9
     read(22,*)
    enddo

! Pour les tables a z_kap= 0.5 , 1.0 on ne doit lire que pour X=0
    do k= 1,nx2
     if (k == 1 .or. k == 2) then
       minz = 2
     else
       minz = 0
     endif
     do l=1,nz2-minz
      do m=1,10
       read(22,*)
      enddo
      do j=1,nt2
       read(22,*) t6_table,(opa(i,j,k,l),i=1,nr2)
      enddo
     enddo
    enddo
  endif

  z_kap = 1.d0 - x_kap - y_kap
  icase  = 0
  icase2 = 0
  icase3 = 0
  icase4 = 0

! Conversion de rh,t en log10(R) = r et en t6
! -------------------------------------------
  t6 = exp(t)/1.0d6
  r  = (rh - 3.d0*t)/um + 18.d0

! On  n'autorise pas a sortir des tables pour R= -5.5,1.5 T6= 0.001, 1500
! Ceci est realise en posant R ou T6 egal a la valeur limite des tables
  if (r>1.5d0 .or. r<-7.5d0 .or. t6>1500.d0 .or. t6<0.001d0) then
    ioutable = ioutable + 1
    rout = r
    tout = t6
  endif
! On fixe la valeur de R et T6 a un peu plus ou moins que le maximum et minimum des tables
  r = min(r,1.5d0)
  r = max(r,-7.5d0)
  t6= min(t6,1500.d0)
  t6= max(t6,0.001d0)

! On cherche la position dans la table en T6
  if (t6 >= tk_kap(nt2)) then
    jt = nt2
  else if (t6 <= tk_kap(1)) then
    jt = 1
  else
    jt = indic(t6,tk_kap,nt2)
  endif

! On cherche la position dans la table en R
  if (r >= rk_kap(nr2)) then
    jr = nr2
  else if (r <= rk_kap(1)) then
    jr =  1
  else
    jr = indic(r,rk_kap,nr2)
  endif

! On cherche la position dans la table en Z
  if (z_kap >= zk_kap(nz2)) then
    jz = nz2
  else if (z_kap <= zk_kap(1)) then
    jz = 1
  else
    jz = indic(z_kap,zk_kap,nz2)
  endif

!  l'opacite demandee se trouve entre les r= rk_kap(jr) et rk_kap(jr+1)
!                                            tk_kap(jt)    tk_kap(jt+1)
!                                            xk_kap(jx)    tk_kap(jx+1)
!                                            zk_kap(jz)    tk_kap(jz+1)
  if (z_kap > 0.5d0 .and. x_kap > 0.0001d0) then
    write(3,*) 'La metallicite est en dehors des tables d''opacite'
    write(3,'(1x,"t6 =",f9.3," r = rho/t6^3 =",f8.3," ln(rho) =",f9.3,"ln(T) =",f9.3,/," X =",f9.3,"Y =",f9.3, &
             & " z_kap =",f9.3)') t6,r,rh,t,x_kap,y_kap,z_kap
  endif
  if (x_kap /= 0.0d0 .and. z_kap >= 0.750d0) then
    rewind(222)
    write(222,*) nwmd,": mixture not covered by the opacity table in kappa93.dat"
    stop "Mixture not covered by the opacity table in kappa93.dat"
  endif
!++++++++++++++++++++++++++++++++++++++++++++++++++
! impression lorsque z_kap est trop grand
  call abundCheck(m,.false.)

  if (z_kap >= zk_kap(nz2)) jz = nz2
  icase = 0
  irmax  =  irmaxtot
  jtmin = 1
  tmin  = 0.002d0
  ixmax = 3

  if (z_kap > 0.04d0 .and. x_kap < 0.00001d0) then
! les parametres sont fixes pour utiliser
! seulement les tables x=0 et interpoler entre 2 tables
! a z_kap=0.5 z_kap=1.0 ou 0.04
    ixmin = 3
    izmax = 2
    if (z_kap >= zk_kap(nz2)) then
      iz = -1
    else
      iz  = 0
    endif
    if (z_kap > 0.5d0) then
      icase = 1
! on devra poser tmin = 1.2  et jtmin = 38 pour izz=2
    endif
  else
    zkm = (zk_kap(jz+1)+zk_kap(jz))/2.d0
! Quand jz=jzmax-3=9 on a 0.04<= Z < 0.5
! On prefera donc d'utiliser les tables a Z=0.03, 0.04 et 0.5
! pour l'interpolation en Z que de prendre les tables
! Z=0.04, 0.5 et 1.0 --> Donc on met iz=-1
    if ((zkm > z_kap .and. jz > 1) .or. jz == jzmax-3) then
      iz = -1
    else
      iz= 0
    endif
    ixmin = 1
    izmax = 3
  endif

  if (r > rk_kap(irmaxtot)) then
     r1 = rk_kap(irmaxtot-1)
     r2 = rk_kap(irmaxtot)
     r3 = r
  else if (r < rk_kap(1)) then
     r1 = r
     r2 = rk_kap(1)
     r3 = rk_kap(2)
  endif

  if (t > tk_kap(nt2)) then
     t1 = tk_kap(nt2-1)
     t2 = tk_kap(nt2)
     t3 = t
  else if (t < tk_kap(jtmin)) then
     t1 = t
     t2 = tk_kap(jtmin)
     t3 = tk_kap(jtmin+1)
  endif

! --  Boucles pricipales
  do izz = 1,izmax
   ixx = ixmin
   do while (ixx <= ixmax)
    if (icase == 1 .and. izz == 2) then
      tmin  = 1.2d0
      jtmin = 38
    endif
    if (jr <= irmaxint-2 .and. r >= rmin)  then
      ixmax = 3
      rkm = (rk_kap(jr+1) + rk_kap(jr) )/2.d0
      if (rkm > r .and. rk_kap(jr) /= rmin) then
        ir = -1
      else
        ir= 0
      endif
      goto 110
    else if (jr == irmaxint-2) then
      ir = 0
      goto 110
    else if (jr == irmaxint-1) then
      ir =  -1
      goto 110
    endif

! ici viennent les cas ou on doit extrapoler pour reponse jr=8 negative
    if (r < rmin) then
      icase2 = 1
      ir =-1
! on doit extrapoler en rmin
      goto 110
    endif
    if (t6 > 0.04d0.and.t6 >= tmin) then
      goto 200
    endif
    if (t6 > tmin) then
      if (t6 > 0.004d0)then
        if (x_kap >= 0.35d0)then
          ixmax=2
          irmax=irmaxtot
        else
          ixmax=3
          if (ixx == 3) then
            irmax=irmaxint
          else
            icase3 = 1
            irmax  = irmaxtot
          endif
          if (z_kap > 0.5d0) then
            write(6,*)'on est trop en dehors des tables d''opacite'
            write(6,'(1x,"t6 =",f9.3," r = rho/t6^3 =",f8.3," ln(rho) =",f9.3,"ln(T) =",f9.3,/," X =",f9.3,"Y =",f9.3,&
                   & " z_kap =",f9.3)') t6,r,rh,t,x_kap,y_kap,z_kap
            call abundCheck(m,.false.)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          endif
!+++++++++++++++++++++++++++++++++++++++++++
        endif
        tkm = (tk_kap(jt) + tk_kap (jt+1) )/2.d0
        if (tkm > t6 .and. jt > jtmin) then
          it = -1
        else
          it =  0
        endif
      else
        jt = 1
        it = 0
        ixmax = 3
      endif

      if (r <= rk_kap(irmax)) then
        rkm = (rk_kap(jr+1) + rk_kap(jr) )/2.d0
        if (rkm > r .or. jr >= irmax-1) then
          ir = -1
          if (jr == irmax) then
            ir = -2
          endif
        else
          ir= 0
        endif
!  -- cas sans extrapolation mais seulement pour 2 valeurs de x: 0.7, 0.35
        goto 120
      else
! on va extrapoler en irmax = 17 ou 13
        ir = -1
        goto 210
      endif

    else
! on va extrapoler en tmin
      goto  300
    endif

110 continue

    if (t6 > tmin .and. jt <= 55) then
! cas sans extrapolation
      tkm = (tk_kap(jt) + tk_kap (jt+1) )/2.d0
      if (tkm > t6 .and. jt > jtmin) then
        it = -1
      else
        it =  0
      endif
      goto 120
    else
      if (t6 < tmin) then
! on  extrapole en tmin

! si l'on est dans le cas ou on utilise les tables
! a z_kap=0.5 et z_kap=1.0  on doit extrapoler en t6 = la
! valeur de jt pour pour la table a z_kap=0.5
! c'est a dire fabrique une table en 1.0 aux valeurs
! de jt de la table a z_kap = 0.5
        if (icase == 0) then
          it = -1
          if (icase2 == 0) then
!    x x x
!    . . .
!    . . .
!
            a21(ixx,izz) = opa(jr+ir  ,jt+it+1,ixx,jz+iz+izz-1)
            a22(ixx,izz) = opa(jr+ir+1,jt+it+1,ixx,jz+iz+izz-1)
            a23(ixx,izz) = opa(jr+ir+2,jt+it+1,ixx,jz+iz+izz-1)
            a31(ixx,izz) = opa(jr+ir  ,jt+it+2,ixx,jz+iz+izz-1)
            a32(ixx,izz) = opa(jr+ir+1,jt+it+2,ixx,jz+iz+izz-1)
            a33(ixx,izz) = opa(jr+ir+2,jt+it+2,ixx,jz+iz+izz-1)
            a11(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),a21(ixx,izz),a31(ixx,izz),t6)
            a12(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),a22(ixx,izz),a32(ixx,izz),t6)
            a13(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),a23(ixx,izz),a33(ixx,izz),t6)
            r1 = rk_kap(jr+ir)
            r2 = rk_kap(jr+ir+1)
            r3 = rk_kap(jr+ir+2)
            t1 = t6
            t2 = tk_kap(jt+it+1)
            t3 = tk_kap(jt+it+2)
          else
!                    on doit aussi extrapoler en rmin = -7
!    x x x
!    x . .
!    x . .
!
            a22(ixx,izz) = opa(jr+ir+1,jt+it+1,ixx,jz+iz+izz-1)
            a23(ixx,izz) = opa(jr+ir+2,jt+it+1,ixx,jz+iz+izz-1)
            a32(ixx,izz) = opa(jr+ir+1,jt+it+2,ixx,jz+iz+izz-1)
            a33(ixx,izz) = opa(jr+ir+2,jt+it+2,ixx,jz+iz+izz-1)
            a21(ixx,izz) = flin(rk_kap(jr+ir+1),rk_kap(jr+ir+2),a22(ixx,izz),a23(ixx,izz),r)
            a31(ixx,izz) = flin(rk_kap(jr+ir+1),rk_kap(jr+ir+2),a32(ixx,izz),a33(ixx,izz),r)
            a11(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),a21(ixx,izz),a31(ixx,izz),t6)
            a12(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),a22(ixx,izz),a32(ixx,izz),t6)
            a13(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),a23(ixx,izz),a33(ixx,izz),t6)
            r1 = r
            r2 = rk_kap(jr+ir+1)
            r3 = rk_kap(jr+ir+2)
            t1 = t6
            t2 = tk_kap(jt+it+1)
            t3 = tk_kap(jt+it+2)
          endif

        else
! on est dans le cas t6<tmin=1.2 et z>0.5
! on garde le meme it que pour le cas a z=0.5 (izz=1)
          if (icase2 == 0) then

!    x x x
!    x x x
!    x x x

!    . . .
!    . . .

            a11(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),opa(jr+ir,jtmin,3,nz2),opa(jr+ir ,jtmin+1,3,nz2),tk_kap(jt+it))
            a12(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),opa(jr+ir+1,jtmin,3,nz2),opa(jr+ir+1,jtmin+1,3,nz2),tk_kap(jt+it))
            a13(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),opa(jr+ir+2,jtmin,3,nz2),opa(jr+ir+2,jtmin+1,3,nz2),tk_kap(jt+it))
            a21(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),opa(jr+ir,jtmin,3,nz2),opa(jr+ir,jtmin+1,3,nz2),tk_kap(jt+it+1))
            a22(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),opa(jr+ir+1,jtmin,3,nz2),opa(jr+ir+1,jtmin+1,3,nz2),tk_kap(jt+it+1))
            a23(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),opa(jr+ir+2,jtmin,3,nz2),opa(jr+ir+2,jtmin+1,3,nz2),tk_kap(jt+it+1))
            a31(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),opa(jr+ir,jtmin,3,nz2),opa(jr+ir,jtmin+1,3,nz2),tk_kap(jt+it+2))
            a32(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),opa(jr+ir+1,jtmin,3,nz2),opa(jr+ir+1,jtmin+1,3,nz2),tk_kap(jt+it+2))
            a33(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),opa(jr+ir+2,jtmin,3,nz2),opa(jr+ir+2,jtmin+1,3,nz2),tk_kap(jt+it+2))

          else

!    x x x
!    x x x
!    x x x

!    x . .
!    x . .
! on ne modifie pas les r1,r2,r3,t1,t2,t3....
! que l'on a pris pour z=0.5 a izz=1
            at = flin(rk_kap(jr+ir+1),rk_kap(jr+ir+2),opa(jr+ir+1,jtmin,3,nz2),opa(jr+ir+2,jtmin,3,nz2),r)
            at1= flin(rk_kap(jr+ir+1),rk_kap(jr+ir+2),opa(jr+ir+1,jtmin+1,3,nz2),opa(jr+ir+2,jtmin+1,3,nz2),r)

            a11(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),at,at1,tk_kap(jt+it))
            a12(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),opa(jr+ir+1,jtmin,3,nz2),opa(jr+ir+1,jtmin+1,3,nz2),tk_kap(jt+it))
            a13(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),opa(jr+ir+2,jtmin,3,nz2),opa(jr+ir+2,jtmin+1,3,nz2),tk_kap(jt+it))
            a21(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),at,at1,tk_kap(jt+it+1))
            a22(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),opa(jr+ir+1,jtmin,3,nz2),opa(jr+ir+1,jtmin+1,3,nz2),tk_kap(jt+it+1))
            a23(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),opa(jr+ir+2,jtmin,3,nz2),opa(jr+ir+2,jtmin+1,3,nz2),tk_kap(jt+it+1))
            a31(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),at,at1,tk_kap(jt+it+2))
            a32(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),opa(jr+ir+1,jtmin,3,nz2),opa(jr+ir+1,jtmin+1,3,nz2),tk_kap(jt+it+2))
            a33(ixx,izz) = flin(tk_kap(jtmin),tk_kap(jtmin+1),opa(jr+ir+2,jtmin,3,nz2),opa(jr+ir+2,jtmin+1,3,nz2),tk_kap(jt+it+2))
          endif
        endif
        goto 500
! en 500 on fit en z s'il n'y a plus de boucle en x
      endif

! on passe ici si jt>= 56
      if (jt == 56) then
        it= -1
        goto 120
      else
!  -- on  extrapole en tmax, jt =57
        goto 130
      endif
    endif
120 continue
! cas sans extrapolation en t
    if (icase2 == 0) then
!   . . .
!   . . .
!   . . .
      a11(ixx,izz) = opa(jr+ir  ,jt+it  ,ixx,jz+iz+izz-1)
      a12(ixx,izz) = opa(jr+ir+1,jt+it  ,ixx,jz+iz+izz-1)
      a13(ixx,izz) = opa(jr+ir+2,jt+it  ,ixx,jz+iz+izz-1)
      a21(ixx,izz) = opa(jr+ir  ,jt+it+1,ixx,jz+iz+izz-1)
      a22(ixx,izz) = opa(jr+ir+1,jt+it+1,ixx,jz+iz+izz-1)
      a23(ixx,izz) = opa(jr+ir+2,jt+it+1,ixx,jz+iz+izz-1)
      a31(ixx,izz) = opa(jr+ir  ,jt+it+2,ixx,jz+iz+izz-1)
      a32(ixx,izz) = opa(jr+ir+1,jt+it+2,ixx,jz+iz+izz-1)
      a33(ixx,izz) = opa(jr+ir+2,jt+it+2,ixx,jz+iz+izz-1)
      r1 = rk_kap(jr+ir)
      r2 = rk_kap(jr+ir+1)
      r3 = rk_kap(jr+ir+2)
      t1 = tk_kap(jt+it)
      t2 = tk_kap(jt+it+1)
      t3 = tk_kap(jt+it+2)
    else
!   x . .
!   x . .
!   x . .
      a12(ixx,izz) = opa(jr+ir+1,jt+it  ,ixx,jz+iz+izz-1)
      a13(ixx,izz) = opa(jr+ir+2,jt+it  ,ixx,jz+iz+izz-1)
      a22(ixx,izz) = opa(jr+ir+1,jt+it+1,ixx,jz+iz+izz-1)
      a23(ixx,izz) = opa(jr+ir+2,jt+it+1,ixx,jz+iz+izz-1)
      a32(ixx,izz) = opa(jr+ir+1,jt+it+2,ixx,jz+iz+izz-1)
      a33(ixx,izz) = opa(jr+ir+2,jt+it+2,ixx,jz+iz+izz-1)
      a11(ixx,izz) = flin(rk_kap(jr+ir+1),rk_kap(jr+ir+2),a12(ixx,izz),a13(ixx,izz),r)
      a21(ixx,izz) = flin(rk_kap(jr+ir+1),rk_kap(jr+ir+2),a22(ixx,izz),a23(ixx,izz),r)
      a31(ixx,izz) = flin(rk_kap(jr+ir+1),rk_kap(jr+ir+2),a32(ixx,izz),a33(ixx,izz),r)
      r1 =  r
      r2 = rk_kap(jr+ir+1)
      r3 = rk_kap(jr+ir+2)
      t1 = tk_kap(jt+it)
      t2 = tk_kap(jt+it+1)
      t3 = tk_kap(jt+it+2)
    endif
    goto 500

130 continue
!  -- on  extrapole en tmax, jt =57
    it = -1
    if (icase2 == 0) then
!    . . .
!    . . .
!    x x x
      a11(ixx,izz) = opa(jr+ir  ,jt+it  ,ixx,jz+iz+izz-1)
      a12(ixx,izz) = opa(jr+ir+1,jt+it  ,ixx,jz+iz+izz-1)
      a13(ixx,izz) = opa(jr+ir+2,jt+it  ,ixx,jz+iz+izz-1)
      a21(ixx,izz) = opa(jr+ir  ,jt+it+1,ixx,jz+iz+izz-1)
      a22(ixx,izz) = opa(jr+ir+1,jt+it+1,ixx,jz+iz+izz-1)
      a23(ixx,izz) = opa(jr+ir+2,jt+it+1,ixx,jz+iz+izz-1)
      a31(ixx,izz) = flin(tk_kap(jt+it),tk_kap(jt+it+1),a11(ixx,izz),a21(ixx,izz),t6)
      a32(ixx,izz) = flin(tk_kap(jt+it),tk_kap(jt+it+1),a12(ixx,izz),a22(ixx,izz),t6)
      a33(ixx,izz) = flin(tk_kap(jt+it),tk_kap(jt+it+1),a13(ixx,izz),a23(ixx,izz),t6)
      r1 = rk_kap(jr+ir)
      r2 = rk_kap(jr+ir+1)
      r3 = rk_kap(jr+ir+2)
      t1 = tk_kap(jt+it)
      t2 = tk_kap(jt+it+1)
      t3 = t6
    else
!  on extrapole en rmin et tmax
!  x . .
!  x . .
!  x x x
      a12(ixx,izz) = opa(jr+ir+1,jt+it  ,ixx,jz+iz+izz-1)
      a13(ixx,izz) = opa(jr+ir+2,jt+it  ,ixx,jz+iz+izz-1)
      a22(ixx,izz) = opa(jr+ir+1,jt+it+1,ixx,jz+iz+izz-1)
      a23(ixx,izz) = opa(jr+ir+2,jt+it+1,ixx,jz+iz+izz-1)
      a11(ixx,izz) = flin(rk_kap(jr+ir+1),rk_kap(jr+ir+2),a12(ixx,izz),a13(ixx,izz),r)
      a21(ixx,izz) = flin(rk_kap(jr+ir+1),rk_kap(jr+ir+2),a22(ixx,izz),a23(ixx,izz),r)
      a31(ixx,izz) = flin(tk_kap(jt+it),tk_kap(jt+it+1),a11(ixx,izz),a21(ixx,izz),t6)
      a32(ixx,izz) = flin(tk_kap(jt+it),tk_kap(jt+it+1),a12(ixx,izz),a22(ixx,izz),t6)
      a33(ixx,izz) = flin(tk_kap(jt+it),tk_kap(jt+it+1),a13(ixx,izz),a23(ixx,izz),t6)
      r1 = r
      r2 = rk_kap(jr+ir+1)
      r3 = rk_kap(jr+ir+2)
      t1 = tk_kap(jt+it)
      t2 = tk_kap(jt+it+1)
      t3 = t6
    endif
    goto 500

140 continue
! on extrapole en irmaxint = 17 et jt = 57
!  . . x
!  . . x
!  x x x
    a11(ixx,izz) = opa(irmaxtot-1,56,ixx,jz+iz+izz-1)
    a12(ixx,izz) = opa(irmaxtot,56,ixx,jz+iz+izz-1)
    a13(ixx,izz) = flin(rk_kap(irmaxtot-1),rk_kap(irmaxtot),a11(ixx,izz),a12(ixx,izz),r)
    a21(ixx,izz) = opa(irmaxtot-1,57,ixx,jz+iz+izz-1)
    a22(ixx,izz) = opa(irmaxtot,57,ixx,jz+iz+izz-1)
    a23(ixx,izz) = flin(rk_kap(irmaxtot-1),rk_kap(irmaxtot),a21(ixx,izz), a22(ixx,izz),r)
    a31(ixx,izz) = flin(tk_kap(56),tk_kap(57),a11(ixx,izz),a21(ixx,izz),t6)
    a32(ixx,izz) = flin(tk_kap(56),tk_kap(57),a12(ixx,izz),a22(ixx,izz),t6)
    a33(ixx,izz) = flin(tk_kap(56),tk_kap(57),a13(ixx,izz),a23(ixx,izz),t6)
    r1 = rk_kap(irmaxtot-1)
    r2 = rk_kap(irmaxtot)
    r3 = r
    t1 = tk_kap(56)
    t2 = tk_kap(57)
    t3 = t6
    goto 500
!  -- cas ou t>0.040 jr >8
200 continue
    if (t6 > 300.d0) then
      if (t6 > 1000.d0) then
        if (r > rmax)  then
! on extrapole en irmax = 17 et jt = 57
          goto 140
        else
!  -- on extrapole en tmax=57
          rkm = ( rk_kap(jr+1) + rk_kap(jr) )/2.d0
          if (rkm > r .or. jr >= 12) then
            ir = -1
            if (jr == irmaxint) ir = -2
          else
            ir= 0
          endif
          goto 130
        endif
      else
        jt = 55
        it = 0
        if (r > rmax) goto 210
        rkm = ( rk_kap(jr+1) + rk_kap(jr) )/2.d0
        if (rkm > r .or. jr >= irmaxint-1) then
          ir = -1
          if (jr == irmaxint) ir = -2
        else
          ir= 0
        endif
        goto 120
      endif
    else
      irmax =  irmaxint
      if (jz+iz+izz-1 == nz2-1) then
        irmax = irmaxint
      endif
      ixmax =  3
      tkm = (tk_kap(jt) + tk_kap (jt+1) )/2.d0
      if (tkm > t6) then
        it = -1
      else
        it =  0
      endif
      goto 210
    endif

! on va extrapoler en irmax = 17 ou 13
210 continue
!  . . x
!  . . x
!  . . x
    if ((icase3 == 1 .and. ixx <= 2) .or. (x_kap >= 0.35d0 .and. t6 <= 0.04d0))then
      icase4 = 1
    endif
    if ((icase4 == 1 .and. ixx < 3) .or. (t6 > 0.040d0)) then
!  -- cas ou on extrapole en r pour les valeurs de ixx=1,2
! ou bien pour t6>0.040 et irmax=17 ou 13

      a11(ixx,izz) = opa(irmax-1,jt+it  ,ixx,jz+iz+izz-1)
      a12(ixx,izz) = opa(irmax  ,jt+it  ,ixx,jz+iz+izz-1)
      a21(ixx,izz) = opa(irmax-1,jt+it+1,ixx,jz+iz+izz-1)
      a22(ixx,izz) = opa(irmax  ,jt+it+1,ixx,jz+iz+izz-1)
      a31(ixx,izz) = opa(irmax-1,jt+it+2,ixx,jz+iz+izz-1)
      a32(ixx,izz) = opa(irmax  ,jt+it+2,ixx,jz+iz+izz-1)
      a13(ixx,izz) = flin(rk_kap(irmax-1),rk_kap(irmax),a11(ixx,izz),a12(ixx,izz),r)
      a23(ixx,izz) = flin(rk_kap(irmax-1),rk_kap(irmax),a21(ixx,izz),a22(ixx,izz),r)
      a33(ixx,izz) = flin(rk_kap(irmax-1),rk_kap(irmax),a31(ixx,izz),a32(ixx,izz),r)
      r1 = rk_kap(irmax-1)
      r2 = rk_kap(irmax)
      r3 = r
      t1 = tk_kap(jt+it)
      t2 = tk_kap(jt+it+1)
      t3 = tk_kap(jt+it+2)
    else if (icase4 == 0 .and. t6 <= 0.04d0) then
!  -- cas ou on extrapole en r seulement pour ixx= 3
! on doit donc extrapoler  pour les memes valeurs de jr et ir qu'a ixx=1,2
! . .    x x x
! . .    x x x
! . .    x x x
      a11(ixx,izz) = flin(rk_kap(irmaxint-5),rk_kap(irmaxint-4),opa(irmaxint-5,jt+it,3,jz+iz+izz-1),opa(irmaxint-4,&
                      jt+it,3,jz+iz+izz-1),rk_kap(jr+ir))
      a12(ixx,izz) = flin(rk_kap(irmaxint-5),rk_kap(irmaxint-4),opa(irmaxint-5,jt+it,3, jz+iz+izz-1),opa(irmaxint-4,&
                      jt+it,3,jz+iz+izz-1),rk_kap(jr+ir+1))
      a13(ixx,izz) = flin(rk_kap(irmaxint-5),rk_kap(irmaxint-4),opa(irmaxint-5,jt+it,3,jz+iz+izz-1),opa(irmaxint-4,&
                      jt+it,3,jz+iz+izz-1),rk_kap(jr+ir+2))
      a21(ixx,izz) = flin(rk_kap(irmaxint-5),rk_kap(irmaxint-4),opa(irmaxint-5,jt+it+1,3,jz+iz+izz-1),opa(irmaxint-4,&
                      jt+it+1,3,jz+iz+izz-1),rk_kap(jr+ir))
      a22(ixx,izz) = flin(rk_kap(irmaxint-5),rk_kap(irmaxint-4),opa(irmaxint-5,jt+it+1,3,jz+iz+izz-1),opa(irmaxint-4,&
                      jt+it+1,3,jz+iz+izz-1),rk_kap(jr+ir+1))
      a23(ixx,izz) = flin(rk_kap(irmaxint-5),rk_kap(irmaxint-4),opa(irmaxint-5,jt+it+1,3,jz+iz+izz-1),opa(irmaxint-4,&
                      jt+it+1,3,jz+iz+izz-1),rk_kap(jr+ir+2))
      a31(ixx,izz) = flin(rk_kap(irmaxint-5),rk_kap(irmaxint-4),opa(irmaxint-5,jt+it+2,3,jz+iz+izz-1),opa(irmaxint-4,&
                      jt+it+2,3,jz+iz+izz-1),rk_kap(jr+ir))
      a32(ixx,izz) = flin(rk_kap(irmaxint-5),rk_kap(irmaxint-4),opa(irmaxint-5,jt+it+2,3,jz+iz+izz-1),opa(irmaxint-4,&
                      jt+it+2,3,jz+iz+izz-1),rk_kap(jr+ir+1))
      a33(ixx,izz) = flin(rk_kap(irmaxint-5),rk_kap(irmaxint-4),opa(irmaxint-5,jt+it+2,3,jz+iz+izz-1),opa(irmaxint-4,&
                      jt+it+2,3,jz+iz+izz-1),rk_kap(jr+ir+2))

! on garde les memes valeurs de  r1,r2,r3,t1,t2,t3
    else if (icase4 == 1) then
!  -- cas ou on extrapole en r  pour ixx= 3
! on doit donc extrapoler pour les memes valeurs de
! jr et jr+1 et ir qu'a ixx=1,2
! mais pour jr+2 on doit prendre r
! . .    x x x
! . .    x x x
! . .    x x x
      a11(ixx,izz) = flin(rk_kap(irmaxint-5),rk_kap(irmaxint-4),opa(irmaxint-5,jt+it,3,jz+iz+izz-1),opa(irmaxint-4,&
                      jt+it,3,jz+iz+izz-1),rk_kap(jr+ir))
      a12(ixx,izz) = flin(rk_kap(irmaxint-5),rk_kap(irmaxint-4),opa(irmaxint-5,jt+it,3,jz+iz+izz-1),opa(irmaxint-4,&
                      jt+it,3,jz+iz+izz-1),rk_kap(jr+ir+1))
      a21(ixx,izz) = flin(rk_kap(irmaxint-5),rk_kap(irmaxint-4),opa(irmaxint-5,jt+it+1,3,jz+iz+izz-1),opa(irmaxint-4,&
                      jt+it+1,3,jz+iz+izz-1),rk_kap(jr+ir))
      a22(ixx,izz) = flin(rk_kap(irmaxint-5),rk_kap(irmaxint-4),opa(irmaxint-5,jt+it+1,3,jz+iz+izz-1),opa(irmaxint-4,&
                      jt+it+1,3,jz+iz+izz-1),rk_kap(jr+ir+1))
      a31(ixx,izz) = flin(rk_kap(irmaxint-5),rk_kap(irmaxint-4),opa(irmaxint-5,jt+it+2,3,jz+iz+izz-1),opa(irmaxint-4,&
                      jt+it+2,3,jz+iz+izz-1),rk_kap(jr+ir))
      a32(ixx,izz) = flin(rk_kap(irmaxint-5),rk_kap(irmaxint-4),opa(irmaxint-5,jt+it+2,3,jz+iz+izz-1),opa(irmaxint-4,&
                      jt+it+2,3,jz+iz+izz-1),rk_kap(jr+ir+1))
      a13(ixx,izz) = flin(rk_kap(irmax-1),rk_kap(irmax),a11(ixx,izz),a12(ixx,izz),r)
      a23(ixx,izz) = flin(rk_kap(irmax-1),rk_kap(irmax),a21(ixx,izz),a22(ixx,izz),r)
      a33(ixx,izz) = flin(rk_kap(irmax-1),rk_kap(irmax),a31(ixx,izz),a32(ixx,izz),r)
! on garde les memes valeurs de  r1,r2,r3,t1,t2,t3
    endif
    goto 500

300 continue

! on extrapole en tmin, dans ce cas irmax = 13 pour toutes les tables
    if (r > rk_kap(irmaxtot)) then
!  on extrapole en tmin et en rmax
! x x x
! . . x
! . . x
      a21(ixx,izz) = opa(irmaxtot-1,jtmin,ixx,jz+iz+izz-1)
      a22(ixx,izz) = opa(irmaxtot,jtmin,ixx,jz+iz+izz-1)
      a31(ixx,izz) = opa(irmaxtot-1,jtmin+1,ixx,jz+iz+izz-1)
      a32(ixx,izz) = opa(irmaxtot,jtmin+1,ixx,jz+iz+izz-1)
      a33(ixx,izz) = flin(rk_kap(irmaxtot-1),rk_kap(irmaxtot),a31(ixx,izz),a32(ixx,izz),r)
      a23(ixx,izz) = flin(rk_kap(irmaxtot-1),rk_kap(irmaxtot),a21(ixx,izz),a22(ixx,izz),r)
      a13(ixx,izz) = flin(tmin+1 ,tmin ,a33(ixx,izz),a23(ixx,izz),t6)
      a11(ixx,izz) = flin(tmin+1, tmin,a31(ixx,izz),a21(ixx,izz),t6)
      a12(ixx,izz) = flin(tmin+1 ,tmin,a32(ixx,izz),a22(ixx,izz),t6)
      r1 = rk_kap(irmaxtot-1)
      r2 = rk_kap(irmaxtot)
      r3 = r
      t1 = t6
      t2 = tk_kap(jtmin)
      t3 = tk_kap(jtmin+1)
    else
      rkm = (rk_kap(jr+1) + rk_kap(jr) )/2.d0
      if (rkm > r .or. jr >= irmaxtot-1) then
        ir = -1
        if (jr == irmaxtot) ir = -2
      else
        ir= 0
      endif
! x x x
! . . .
! . . .
      a21(ixx,izz) = opa(jr+ir  ,jtmin,ixx,jz+iz+izz-1)
      a22(ixx,izz) = opa(jr+ir+1,jtmin,ixx,jz+iz+izz-1)
      a23(ixx,izz) = opa(jr+ir+2,jtmin,ixx,jz+iz+izz-1)
      a31(ixx,izz) = opa(jr+ir  ,jtmin+1,ixx,jz+iz+izz-1)
      a32(ixx,izz) = opa(jr+ir+1,jtmin+1,ixx,jz+iz+izz-1)
      a33(ixx,izz) = opa(jr+ir+2,jtmin+1,ixx,jz+iz+izz-1)
      a11(ixx,izz) = flin(tmin+1 ,tmin ,a31(ixx,izz),a21(ixx,izz),t6)
      a12(ixx,izz) = flin(tmin+1 ,tmin ,a32(ixx,izz),a22(ixx,izz),t6)
      a13(ixx,izz) = flin(tmin+1 ,tmin  ,a33(ixx,izz),a23(ixx,izz),t6)
      r1 = rk_kap(jr+ir)
      r2 = rk_kap(jr+ir+1)
      r3 = rk_kap(jr+ir+2)
      t1 = t6
      t2 = tk_kap(jtmin)
      t3 = tk_kap(jtmin+1)
    endif

500 continue
    ixx = ixx +1
   enddo
  enddo

! fit en z
  if (izmax == 3) then
    do ixxx=1,ixmax
     b11_kap(ixxx) = qua(zk_kap(jz+iz),zk_kap(jz+iz+1),zk_kap(jz+iz+2),a11(ixxx,1),a11(ixxx,2),a11(ixxx,3),z_kap)
     b12(ixxx) = qua(zk_kap(jz+iz),zk_kap(jz+iz+1),zk_kap(jz+iz+2),a12(ixxx,1),a12(ixxx,2),a12(ixxx,3),z_kap)
     b13(ixxx) = qua(zk_kap(jz+iz),zk_kap(jz+iz+1),zk_kap(jz+iz+2),a13(ixxx,1),a13(ixxx,2),a13(ixxx,3),z_kap)
     b21(ixxx) = qua(zk_kap(jz+iz),zk_kap(jz+iz+1),zk_kap(jz+iz+2),a21(ixxx,1),a21(ixxx,2),a21(ixxx,3),z_kap)
     b22(ixxx) = qua(zk_kap(jz+iz),zk_kap(jz+iz+1),zk_kap(jz+iz+2),a22(ixxx,1),a22(ixxx,2),a22(ixxx,3),z_kap)
     b23(ixxx) = qua(zk_kap(jz+iz),zk_kap(jz+iz+1),zk_kap(jz+iz+2),a23(ixxx,1),a23(ixxx,2),a23(ixxx,3),z_kap)
     b31(ixxx) = qua(zk_kap(jz+iz),zk_kap(jz+iz+1),zk_kap(jz+iz+2),a31(ixxx,1),a31(ixxx,2),a31(ixxx,3),z_kap)
     b32(ixxx) = qua(zk_kap(jz+iz),zk_kap(jz+iz+1),zk_kap(jz+iz+2),a32(ixxx,1),a32(ixxx,2),a32(ixxx,3),z_kap)
     b33_kap(ixxx) = qua(zk_kap(jz+iz),zk_kap(jz+iz+1),zk_kap(jz+iz+2),a33(ixxx,1),a33(ixxx,2),a33(ixxx,3),z_kap)
    enddo

! on interpole en x
! si ixmax = 2 on interpole lineairement entre x= 0.7 et x= 0.35
! car on devrait faire une extrapolation pour x=0 dans ce domaine de r et t
    if (ixmax == 3) then
      c11 = qua(xk_kap(1),xk_kap(2),xk_kap(3),b11_kap(1),b11_kap(2),b11_kap(3),x_kap)
      c12 = qua(xk_kap(1),xk_kap(2),xk_kap(3),b12(1),b12(2),b12(3),x_kap)
      c13 = qua(xk_kap(1),xk_kap(2),xk_kap(3),b13(1),b13(2),b13(3),x_kap)
      c21 = qua(xk_kap(1),xk_kap(2),xk_kap(3),b21(1),b21(2),b21(3),x_kap)
      c22 = qua(xk_kap(1),xk_kap(2),xk_kap(3),b22(1),b22(2),b22(3),x_kap)
      c23 = qua(xk_kap(1),xk_kap(2),xk_kap(3),b23(1),b23(2),b23(3),x_kap)
      c31 = qua(xk_kap(1),xk_kap(2),xk_kap(3),b31(1),b31(2),b31(3),x_kap)
      c32 = qua(xk_kap(1),xk_kap(2),xk_kap(3),b32(1),b32(2),b32(3),x_kap)
      c33 = qua(xk_kap(1),xk_kap(2),xk_kap(3),b33_kap(1),b33_kap(2),b33_kap(3),x_kap)
    else
      c11 = flin(xk_kap(1),xk_kap(2),b11_kap(1),b11_kap(2),x_kap)
      c12 = flin(xk_kap(1),xk_kap(2),b12(1),b12(2),x_kap)
      c13 = flin(xk_kap(1),xk_kap(2),b13(1),b13(2),x_kap)
      c21 = flin(xk_kap(1),xk_kap(2),b21(1),b21(2),x_kap)
      c22 = flin(xk_kap(1),xk_kap(2),b22(1),b22(2),x_kap)
      c23 = flin(xk_kap(1),xk_kap(2),b23(1),b23(2),x_kap)
      c31 = flin(xk_kap(1),xk_kap(2),b31(1),b31(2),x_kap)
      c32 = flin(xk_kap(1),xk_kap(2),b32(1),b32(2),x_kap)
      c33 = flin(xk_kap(1),xk_kap(2),b33_kap(1),b33_kap(2),x_kap)
    endif
  else
! on a que 2 metallicites et x=0
    ixxx = 3
    c11 = flin(zk_kap(jz+iz),zk_kap(jz+iz+1),a11(ixxx,1),a11(ixxx,2),z_kap)
    c12 = flin(zk_kap(jz+iz),zk_kap(jz+iz+1),a12(ixxx,1),a12(ixxx,2),z_kap)
    c13 = flin(zk_kap(jz+iz),zk_kap(jz+iz+1),a13(ixxx,1),a13(ixxx,2),z_kap)
    c21 = flin(zk_kap(jz+iz),zk_kap(jz+iz+1),a21(ixxx,1),a21(ixxx,2),z_kap)
    c22 = flin(zk_kap(jz+iz),zk_kap(jz+iz+1),a22(ixxx,1),a22(ixxx,2),z_kap)
    c23 = flin(zk_kap(jz+iz),zk_kap(jz+iz+1),a23(ixxx,1),a23(ixxx,2),z_kap)
    c31 = flin(zk_kap(jz+iz),zk_kap(jz+iz+1),a31(ixxx,1),a31(ixxx,2),z_kap)
    c32 = flin(zk_kap(jz+iz),zk_kap(jz+iz+1),a32(ixxx,1),a32(ixxx,2),z_kap)
    c33 = flin(zk_kap(jz+iz),zk_kap(jz+iz+1),a33(ixxx,1),a33(ixxx,2),z_kap)
  endif

!  -- on interpole en r et t6
!     -----------------------

! on interpole d'abord en r pour t1 , t2 et t3
! --------------------------------------------
  frt1 = qua(r1,r2,r3,c11,c12,c13,r)
  frt2 = qua(r1,r2,r3,c21,c22,c23,r)
  frt3 = qua(r1,r2,r3,c31,c32,c33,r)

! on interpole en t6  pour r1 , r2 et r3
! --------------------------------------
  ftr1 = qua(t1,t2,t3,c11,c21,c31,t6)
  ftr2 = qua(t1,t2,t3,c12,c22,c32,t6)
  ftr3 = qua(t1,t2,t3,c13,c23,c33,t6)
  cap10 =  qua(t1,t2,t3,frt1,frt2,frt3,t6)
  if (isnan(cap10)) then
    rewind(222)
    write(222,*) nwmd,": NaN in kappa_out"
    stop "NaN in kappa_out"
  endif
! on transforme l'opacite de log10 en ln, car la routine doit retourner ln(kappa)
! et les opacites lues dans les tables sont en log10.
  cap   = cap10 * um

! ici on utilise les definitions suivantes
! captt   = [dln(kappa)/ dln(t) ]r=cte  =  [dlog10(kappa)/ dlog10(t) ]r=cte
! caprr   = [dln(kappa)/ dln(R) ]t=cte  =  [dlog10(kappa)/ dlog10(R) ]t=cte
! Avec les relations suivantes:
!   t = 10^6 * t6
!   log10(t) = 6. + log10(t6) = 6. + log(t6)/um
!   on evalue les derivees pour log10(t) puisque les opacites son donnees en log10
  captt = quad_gg(log(t1)/um + 6.d0,log(t2)/um + 6.d0,log(t3)/um + 6.d0,frt1,frt2,frt3,log(t6)/um + 6.d0)

! ri = log10(Ri)
  caprr = quad_gg(r1,r2,r3,ftr1,ftr2,ftr3,r)

! on transforme les derivees de r et t en p et t
! ----------------------------------------------
  capt =  caprr * (rht - 3.d0)  + captt
  capp =  caprr * rhp

! On calcule les opacités conductives et on les ajoute si
! kappa_cond > kappa_rad / 100
  call condTest(t,rh,rht,rhp,cap,x_kap,y_kap,jj1,capp,capt)

  return

end subroutine kappa_out

!=======================================================================
subroutine condTest(t,rh,rht,rhp,cap,x_kap,y_kap,jj1,capp,capt)
! Test s'il faut prendre en compte  l'opacite conductive
! ----------------------------------------------------------------------
 use const,only: um
 use abundmod,only: xc12,xn14,xo16,xne20,xne22,xmg24

 implicit none

 real(kindreal),intent(in):: t,rh,rht,rhp,x_kap,y_kap
 integer,intent(inout):: jj1
 real(kindreal),intent(inout):: cap,capp,capt

 real(kindreal):: tl,rho,caprd,ccon,conrt,contr,r1,captcond,cappcond
 real(kindreal),dimension(8):: xx
!------------------------------------------------------------------------
  tl  = t/um
  rho = rh/um
  if (rho < -4.d0 .or. tl < 4.3d0) return
  caprd=exp(cap)
  xx(1)=x_kap
  xx(2)=y_kap
  if (jj1 < 1) jj1=1
  xx(3)=xc12(jj1)
  xx(4)=xn14(jj1)
  xx(5)=xo16(jj1)
  xx(6)=xne20(jj1)
  xx(7)=xne22(jj1)
  xx(8)=xmg24(jj1)
  call cond(rho,tl,xx,caprd,ccon,conrt,contr)

  if (ccon/caprd > 100.d0) return
! ici on prend l'opacite moyenne cond + rad et pour les derivees on prend seulement celles conductives
! que l'on transforme de rho,temp en pression,temp de maniere a pouvoir les combines si necessaire aux derivees
! radiatives exprimees en fonction de p et t

  cap=dlog(caprd*ccon/(caprd+ccon))
  captcond = conrt * rht + contr
  cappcond = conrt * rhp

  if (caprd/ccon > 100.d0) then
! ici  on garde les derivees de opa_cond car opa_cond < opa_rad/100 est vraiment petit
    capt = captcond
    capp = cappcond
    return
  endif
! ici on est dans le cas intermediaire ou opa_rad/100 < opa_con < 100 * opa_rad
!  on  pondere aussi les derivees

  r1=ccon/(caprd+ccon)
  capt = r1*capt + (1.d0-r1)*captcond
  capp = r1*capp + (1.d0-r1)*cappcond

  return

end subroutine condTest

!=======================================================================
subroutine cond(rho,tl,xx,caprd,ccon,conrt,contr)
!---  CALCULATION OF CONDUCTIVE OPACITY (Iben 1975, ApJ 196, 525 Appendix A)
!------------------------------------------------------------------------
  use const,only: pi
  use SmallFunc,only: expf10

  implicit none

  real(kindreal),intent(in):: rho,tl,caprd
  real(kindreal),dimension(8),intent(in):: xx
  real(kindreal),intent(out):: ccon,conrt,contr

  integer::i,k

  real(kindreal):: vmye,vmyel,za,zb,xz2,zal,rolg,tlg,dlro,dlt,ro6lg,t6lg,dellg,del,eta0,eta2,a1,b1,a2,b2,rnedne,flg,penktl,blamr2,&
       alpha,thxlg,thylg,thclg,thx,thy,thc,zc,f,vkch=0.d0,vcond,ef,gam,efm,glg,vkcc
  real(kindreal),dimension(3):: vkc
  real(kindreal),dimension(8),parameter:: z=(/1.d0,2.d0,6.d0,7.d0,8.d0,10.d0,10.d0,12.d0/),&
                                              a=(/1.d0,4.d0,12.d0,14.d0,16.d0,20.d0,22.d0,24.d0/)
!------------------------------------------------------------------------
  vmye=2.d0/(1.d0+xx(1))
  vmyel=log10(vmye)
  za=0.d0
  zb=0.d0
  do i=1,8
   xz2=xx(i)*z(i)*z(i)
   za=za+xz2/(a(i)**(1.d0/3.d0))
   zb=zb+xz2/a(i)
  enddo
  zal=log10(za)

  do k=1,3
  rolg=rho
  tlg=tl
  select case (k)
    case (1)
      dlro=0.d0
      dlt=0.d0
    case (2)
      dlro=0.001d0
      dlt=0.d0
    case (3)
      dlro=0.d0
      dlt=0.0001d0
    case default
      stop 'problem with k in cond'
  end select
  rolg=rolg+dlro
  tlg=tlg+dlt
  ro6lg=rolg-6.d0
  t6lg=tlg-6.d0
  if (ro6lg <= 0.3d0) then
    dellg=rolg-1.5d0*t6lg-vmyel
    del=expf10(dellg)
    eta0=expf10(-0.52255d0+2.d0*dellg/3.d0)
    eta2=eta0*eta0
    a1=-3.29243d0+log10(del*(1.d0+0.02804d0*del))
    b1=-4.80946d0+log10(del*del*(1.d0+9.376d0/eta2))
    a2=log10(1.d0+0.021876d0*del)
    b2=log10(0.4d0*eta0*(1.d0+4.1124d0/eta2))
    if (del <= 40.d0) then
      rnedne=1.0d0-0.01d0*del*(2.8966d0-0.034838d0*del)
    else
      rnedne=(1.5d0/eta0)*(1.d0-0.8225d0/eta2)
    endif
    flg=a1
    if (dellg <= 0.645d0) flg=-3.2862d0+log10(del*(1.d0+0.024417d0*del))
    if (dellg >= 2.5d0) flg=b1
    if (dellg > 2.d0 .and. dellg < 2.5d0) flg=2.d0*a1*(2.5d0-dellg)+2.d0*b1*(dellg-2.d0)
    penktl=a2
    if (dellg>1.5d0 .and. dellg<2.d0) penktl=2.d0*a2*(2.d0-dellg)+2.d0*b2*(dellg-1.5d0)
    if (dellg >= 2.d0) penktl=b2
    blamr2=9.24735d-3*expf10(dellg-0.5d0*t6lg-penktl)*(vmye*zb+rnedne)
    alpha=0.5d0*log10(blamr2)
    thxlg=0.13d0-alpha*(0.745d0+0.105d0*alpha)
    thylg=0.24d0-alpha*(0.55d0+0.0689d0*alpha)
    if (alpha <= -3.d0) then
      thxlg=1.048d0-0.124d0*alpha
      thylg=0.937d0-0.111d0*alpha
    endif
    if (alpha > -1.d0) thxlg=0.185d0-0.558d0*alpha
    if (alpha > 0.d0) thylg=0.24d0-0.6d0*alpha
    thclg=0.727d0-alpha*(0.511d0+0.0778d0*alpha)
    if (alpha <= -2.5d0) thclg=1.27d0-0.1d0*alpha
    if (alpha > 0.5d0) thclg=0.843d0-0.785d0*alpha
    thx=expf10(thxlg)
    thy=expf10(thylg)
    thc=expf10(thclg)

!---  ZC IS AN EFFECTIVE ABUNDANCE OF THE ELEMENTS C,N,O,NE WEIGTHED
!     BY THE SQUARES OF THE CHARGES AND NORMALIZED WITH RESPECT TO C
!     SINCE ELECTRON CONDUCTION HAS BEEN CALCULATED BY HUBBARD+LAMPE
!     ONLY FOR PURE H,PURE HE AND PURE C RESP.
    zc=(3.d0*xx(3)+3.5d0*xx(4)+4.d0*xx(5)+5.d0*xx(6)+4.54d0*xx(7)+6.d0*xx(8))/3.d0
    vkch=(xx(1)*thx+xx(2)*thy+zc*thc)*expf10(-t6lg-flg)
  endif
  if (ro6lg <= 0.d0) then
    vcond=vkch
  else
    ef=sqrt(1.d0+expf10(2.d0*(ro6lg-vmyel)/3.d0)) - 1.d0
    gam=22.76d0*expf10(ro6lg/3.d0-t6lg)*za
    efm=min(1.d0,0.5d0+log10(ef))
    glg=(0.873d0-0.298d0*zal+(0.333d0-0.168d0*zal)*efm)*(1.d0-(1.d0+gam)**(-0.85d0))
    vkcc=6.753d-8*expf10(2.d0*t6lg-glg)*zb/(ef**2.5d0*(1.d0+ef))
    vcond=vkcc
    if (ro6lg <= 0.3d0) then
      f=0.5d0*(1.d0-cos(pi*ro6lg/0.3d0))
      vcond=expf10((1.d0-f)*log10(vkch)+f*log10(vkcc))
    endif
  endif
  if (k == 1) then
    ccon=vcond
    conrt=0.d0
    contr=0.d0
    if (ccon/caprd > 100.d0) return
  endif
  vkc(k)=log10(vcond)
  enddo
  conrt=1000.d0*(vkc(2)-vkc(1))
  contr=10000.d0*(vkc(3)-vkc(1))

  return

end subroutine cond

!=======================================================================
!> Computation of the opacity in a layer.
!! @brief The program takes rho, T, X and Y as input and returns kappa.
!!        The computation of kappa depends on the value of ikappa:
!!           IKAPPA= 9: electron scattering opacity
!!           IKAPPA= 5: tabulated opacities (OPAL + Alexander & Ferguson)
!!        The opacities are tabulated in terms of T_6 and R=rho/T_6^3. In the tables, they are in log10
!!        A quadratic interpolation of the form y= a + b X + c X**2 is performed in X, Z, R et T.
subroutine kappa(rh,t,rhp,rht,x_kap,y_kap,cap,capp,capt,jj1)
!------------------------------------------------------------------------
!    Program de gestion des tables d'opacites          version mars 1992
!    ----------------------------------------

!    Les tables se trouvent dans le fichier kappa.dat
!    Le programme  lit les tables au premier appel.
!    Les opacites lues dans les tables sont en log10(kappa)
!    L'opacite retournee cap:= ln(kappa)

!    Les tables sont donnees dans les 2 variables:
!          r = log10(R) = log10 ( dens/T6**3 )
!       et t6 (temperature en million de degres).

!    Il effectue une interpolation quadratique de la forme y= a + b X + c X**2
!    en X, Z, R et T. Lorsque que l'on se trouve en dehors des tables il fait
!    une extrapolation lineaire en donnant un message d'information.
!    On fait d'abord un fit en Z pour les X utilises, en general pour
!    les 3 valeurs de X et ceci pour les 9 points en r et t.
!    Les 9 variables a12(3,3), a12(3,3), ... ,a33(3,3) contiennent les valeurs
!    de l'opacite des 9 points representatifs (si l'on est dans un bord
!    certaines peuvent provenir d'une extrapolation lineaire) pour
!    les 3 valeurs de x et 3 valeurs de Z.
!    Les variables b11_kap(3), b12(3), ... , b33_kap(3) contiennent
!    les 9 points representatif  pour les 3 valeurs de X mais Z "fitter".
!    Les c11, c12, ..., c33 contiennnent les opacites fittees en x
!    a partir des bij.
!    On peut ainsi calculer l'opacite en r,t et les derivees.

!    definition des principales variables:
!    ------------------------------------
!    rh   := ln(densite)
!    rho  := log10(densite)
!    dens := densite (g/cm**3)
!    t    := ln(temperature)
!    t6   := temperature/10.**6
!    temp := log10(temperature)
!    R    := dens/T6**3
!    r    := log10 (R)
!    nt   := nombre de valeurs de t
!    nr   := nombre de valeurs de r
!    nx   := nombre de valeurs de x
!    nz   := nombre de valeurs de z
!    lec  := vaut 0 si les tables n'ont pas encore ete lues sinon vaut 1

!    jr,jt, jz: sont les positions dans les tables de la valeur de
!              respectivement r,t et z
!              si jr = 3  la valeur de r se trouve entre rk_kap(3) et rk_kap(4).
!    it,ir, iz: sont des variables de decalage en t,r et z.
!              On prend toujours les 3 points de la table defini par les
!              positions jr, jr + 1, jr + 2 mais ceci n'est valable que si
!              r est plus pres de jr+1 que de jr.
!              Pour ce faire on definit la variable
!    rkm, tkm, zkm: densite, temperature et metal. moyenne
!                   comme [rk_kap(j)+rk_kap(j+1)] /2.
!                   si rkm>r on pose ir=-1
!                            et on prend jr-ir, jr-ir + 1, jr-ir + 2.
!    izz,ixx variables comptant les trois tables utilisees en x ou en z.
!
!    ixmax nombre de tables en x utilisees:  si x>0.35 et r,t  sont en dehors de la table a x=0
!                                            mais tombent dans les tables a x=0.35,0.7 alors
!                                            ixmax = 2. De meme si x=0  et z>0.04 on
!                                            pose ixmin= 3 et ixmax = 3 ainsi on utilise seulement
!                                            les table pour x=0.
!    on utilise une representation symbolique des differents cas par des
!    commentaires. Par exemple on represente par des points "." les valeurs
!    lues dans les tables  et par une croix celles qui sont tirees d'une
!    extrapolation.
!    x x x  Dans ce cas on extrapole en tmin et en rmax
!    . . x
!    . . x

!    Les limites physiques des tables sont: logR = [-7;1] et  T6 = [0.002,1000]
!    Dans la nouvelle version on autorise une extrapolation lineaire
!    des tables jusque dans les limites:
!        logR = [-7.5;1.5] et  T6 = [0.001,1500]
!    Au dela de ces limites on prend kappa = constante = valeur du bord de la
!    table.

!    Si ces limites sont depassees lors de la derniere iteration du modele
!    en cours, un message est affiche a l'ecran apres celui du no du modele
!    et temps. Le nombre de couches concernees ainsi que la valeur extreme
!    de depassement de R et T6 est donnee.
!------------------------------------------------------------------------
  use const,only: um
  use inputparam,only: ikappa

  implicit none

  real(kindreal),intent(in):: rh,t,rhp,rht,x_kap,y_kap
  integer,intent(inout)::jj1
  real(kindreal),intent(out):: cap,capp,capt

  real(kindreal):: z_kap,t6,r,captt,caprr
!------------------------------------------------------------------------
! Calcul d'opacites de 2 manieres differentes:
!  1) ikappa=5:
!     On utilise version 95 d'OPAL + Alexander & Fergson 94 pour le domaine
!     (R,T6) a disposition mais UNIQUEMENT pour X>0. !
!     A l'exterieur du domaine on utilisera kappa93.dat,
!     qui est la version mars 1992 (Schaller et al. 1992) avec les correctes
!     opacites de Kurucz.

!     Comme Z=0.1 est la valeur maximale disponible de ces sources nous
!     utilisons les tables de Huebner (non modifiees contenues dans
!     kappa93.dat) pour les cas avec X=0. et donc Z>0.1.

!  2) ikappa=9:
!     Seulement la diffusion par electron libre est prise en compte
!-----------------------------------------------------------------------
  error_xzt = 0

  select case (ikappa)
  case (9)   ! electron-scattering opacity only
    cap=log(0.2d0*(1.d0+x_kap))
    capp=0.d0
    capt=0.d0
    return

  case (5)   ! tabulated opacities (OPAL + Alexander & Ferguson)
!  ATTENTION NOUS N'AVONS PAS DE VALEURS A FAIBLES T LORSQUE X > 0.80
    z_kap = 1.d0 - x_kap - y_kap

!    Conversion de rh,t en R = r et en t6
    t6 = exp(t)/10.d0**6.d0
    r  = (rh - 3.d0*t)/um + 18.d0

!  Utilisons 3) pour X=0 et z_kap > 0.1.
    if (x_kap <= 0.d0 .and. z_kap > 0.1d0) then
      call kappa_out(rh,t,rhp,rht,x_kap,y_kap,cap,capp,capt,jj1)
      return
    endif
! ATTENTION ON UTILISE ALORS MELANGE SOLAIRE
! POUR LE MOMENT ON NE PEUT FAIRE AUTREMENT

! Si on sort du domaine maximum de la table (log R= -8 .. 1 et
! log T= 3.00, 8.70) on utilisera directement 3)
! [note: dans ce cas une extrapolation sera faite...]
    if (r > 1.d0 .or. r < -8.d0 .or. t6 > 501.2d0 .or. t6 < 0.001d0) then
      call kappa_out(rh,t,rhp,rht,x_kap,y_kap,cap,capp,capt,jj1)
      return
    endif

! De meme pour le domaine a basse temperature (log T< 3.7) et log R < -7
! sans donnees de Alexander & Ferguson:
    if (r < -7.d0 .and. t6 < 0.005012d0) then
      call kappa_out(rh,t,rhp,rht,x_kap,y_kap,cap,capp,capt,jj1)
      return
    endif

    r = 10.d0**r

    call opacgn93 (z_kap,x_kap,t6,r)

! if outside of R-T6 plane (since not entirely rectangular) or
! log(kappa)>=9. (dummy value to signal incomplete table) use kappa93.dat

    if (error_xzt == 1 .or. opact >= 9.d0) then
      call kappa_out(rh,t,rhp,rht,x_kap,y_kap,cap,capp,capt,jj1)
      return
    endif

!   ici on utilise les definitions suivantes
!   captt   = [dln(kappa)/ dln(t) ]r=cte  =  [dlog10(kappa)/ dlog10(t) ]r=cte
!   caprr   = [dln(kappa)/ dln(R) ]t=cte  =  [dlog10(kappa)/ dlog10(R) ]t=cte
    cap = opact * um
    captt = dopact
    caprr = dopacr
!   on transforme les derivees de r et t en p et t
!   ----------------------------------------------
    capt = caprr * (rht - 3.d0) + captt
    capp = caprr * rhp

! On termine en testant s'il faut inclure l'opacite conductive...
    call condTest(t,rh,rht,rhp,cap,x_kap,y_kap,jj1,capp,capt)

    return

  case default

    stop 'kappa2009.f: Set ikappa=5 or 9 !!'

  end select

  call kappa_out(rh,t,rhp,rht,x_kap,y_kap,cap,capp,capt,jj1)

  return

end subroutine kappa
!=======================================================================

end module opacity
!*********************************************************************
