!> \file geomod.f95
!> \brief ftfp.dat interpolation module
!>
!> Used for interpolating in the ftfp table
module geomod

  use evol,only: kindreal
  use interpolation,only: fipoi,fipoi1,fipoi2

  implicit none

  integer,private,parameter:: n_geo=5021
  real(kindreal), dimension(n_geo),private,save::rpsi,oblat,fp,ft,rapt,ratp,sund,rsurf,gpsgmo,gmoym,xinerm, &
                                                 xgnorp,omegageo,GammaEddMax
  real(kindreal),save::rpsi_min,rpsi_max,oblat_min,oblat_max,fp_min,fp_max,ft_min,ft_max,rapt_min,rapt_max, &
                      ratp_min,ratp_max,sund_min,sund_max,rsurf_min,rsurf_max,gpsgmo_min,gpsgmo_max,gmoym_min, &
                      gmoym_max,xinerm_min,xinerm_max,xgnorp_min,xgnorp_max,omegageo_min,omegageo_max, &
                      GammaEddMax_min,GammaEddMax_max

  private
  public :: initgeo
  public :: geocalc
  public :: geom
  public :: geomat
  public :: geomeang
  public :: geomEdd
  public :: rpsi_min,rpsi_max,sund_min,sund_max,GammaEddMax_min

contains

!-----------------------------------------------------------------------
!> Read the ftfp file at the begginning of the execution
  subroutine initgeo
!  col  1 - rpsi : rayon moyen, given by volume_star = 4/3 pi rpsi^3
!  col  2 - oblat : rpol/req
!  col  3 - fp : f_p as in Eq. 4.21 in Maeder (2008) except for the non-geometrical terms
!  col  4 - ft : f_t as in Eq. 4.13 in Maeder (2008)
!  col  5 - rapt : f_p/f_t
!  col  6 - ratp : Paris subway company
!  col  7 - sund : surface of the isobar S_psi
!  col  8 - rsurf : ratio surface_isobar / 4 pi rpsi^2
!  col  9 - gpsgmo : effective gravity at the point where the radius of the isobar equals r_psi
!                   divided by the average gravity on the surface
!  col 10 - gmoym : surface averaged value of 1/g_eff
!  col 11 - xinerm : surface averaged value of (1/g_eff r_mean^2 sin^2(theta))
!                    divided by gmoym
!  col 12 - xgnorp : same as gpsgmo, but not divided by the average gravity on the surface
!  col 13 - omegageo : Omega/Omega_crit ratio
!  col 14 - GammaEdd : highest Eddington factor value possible for this Omega/Omega_crit

    use evol,only: input_dir

    implicit none
    integer:: i

    open(file=trim(input_dir)//'inputs/ftfp.dat',unit=30,status='old')
    do i=1,n_geo
     read (30,'(1x,14(2x,e20.10))') rpsi(i),oblat(i),fp(i),ft(i),rapt(i),ratp(i),sund(i),rsurf(i), &
                                    gpsgmo(i),gmoym(i),xinerm(i),xgnorp(i),omegageo(i),GammaEddMax(i)
    enddo
    close(30)

    rpsi_min = minval(rpsi)
    rpsi_max = maxval(rpsi)
    oblat_min = minval(oblat)
    oblat_max = maxval(oblat)
    fp_min = minval(fp)
    fp_max = maxval(fp)
    ft_min = minval(ft)
    ft_max = maxval(ft)
    rapt_min = minval(rapt)
    rapt_max = maxval(rapt)
    ratp_min = minval(ratp)
    ratp_max = maxval(ratp)
    sund_min = minval(sund)
    sund_max = maxval(sund)
    rsurf_min = minval(rsurf)
    rsurf_max = maxval(rsurf)
    gpsgmo_min = minval(gpsgmo)
    gpsgmo_max = maxval(gpsgmo)
    gmoym_min = minval(gmoym)
    gmoym_max = maxval(gmoym)
    xinerm_min = minval(xinerm)
    xinerm_max = maxval(xinerm)
    xgnorp_min = maxval(xgnorp)
    xgnorp_max = maxval(xgnorp)
    omegageo_min = minval(omegageo)
    omegageo_max = maxval(omegageo)
    GammaEddMax_min = minval(GammaEddMax)
    GammaEddMax_max = maxval(GammaEddMax)

  end subroutine initgeo

!-----------------------------------------------------------------------
!> Replaces the geocin (1), geograv (2), and geocrit (4) routines
!! Called by geomeang, geomat, and geomedd.
  subroutine geocalc(xx,yy,geotype)
! geotype: 1= ancien geograv:    xpsi --> xgpsi      fipoi1(xpsi,n_geo,rpsi,xgnorp,xgpsi)
!          2= ancien geomeang:   xpsi --> xrapg      fipoi1(xpsi,n_geo,rpsi,gpsgmo,xrapg)
!          3= ancien geocrit:    xpsi --> xobla      fipoi1(xpsi,n_geo,rpsi,oblat,xobla)
!          4= ancien geocrit:    spsi --> xobla      fipoi1(spsi,n_geo,sund,oblat,xobla)
!          5= ancien geomat (1): spsi --> xpsi       fipoi1(spsi,n_geo,sund,rpsi,xpsi)
!          6= ancien geomat (2): spsi --> rap1       fipoi1(spsi,n_geo,sund,rsurf,rap1)
!          7= ancien geomat (3): spsi --> xft        fipoi1(spsi,n_geo,sund,ft,xft)
!          8= ancien geomat (4): spsi --> rap2       fipoi1(spsi,n_geo,sund,gpsgmo,rap2)
!          9= ancien geomat (5): spsi --> xgmoym     fipoi1(spsi,n_geo,sund,gmoym,xgmoym)
!         10= ancien geomedd(1): eddesc --> xoblaMax fipoi1(eddesc,n_geo,GammaEddMax,oblat,xoblaMax)
!         11= ancien geomedd(2): eddesc --> OmegaMax fipoi1(eddesc,n_geo,GammaEddMax,omegageo,OmegaMax)

    implicit none

    integer, intent(in):: geotype
    real(kindreal), intent(inout):: xx
    real(kindreal), intent(out):: yy
    real(kindreal), dimension(n_geo):: interpx,interpy

    select case (geotype)
      case (1:3)
        if (xx < rpsi_min) then
          xx = rpsi_min
        endif
        interpx(:) = rpsi(:)
      case (4:9)
        if (xx < sund_min) then
          xx = sund_min
        endif
        interpx(:) = sund(:)
      case (10,11)
        if (xx < GammaEddmax_min) then
          xx = GammaEddMax_min
        endif
        interpx(:) = GammaEddMax(:)
      case default
        stop 'Bad choice of input variable for interpolation in ftfp'
    end select

    select case (geotype)
      case (1)
        interpy(:) = xgnorp(:)
      case (2,8)
        interpy(:) = gpsgmo(:)
      case (3,4,10)
        interpy(:) = oblat(:)
      case (5)
        interpy(:) = rpsi(:)
      case (6)
        interpy(:) = rsurf(:)
      case (7)
        interpy(:) = ft(:)
      case (9)
        interpy(:) = gmoym(:)
      case (11)
        interpy(:) = omegageo(:)
      case default
        stop 'Bad choice of output variable for interpolation in ftfp'
    end select

    call fipoi1(xx,n_geo,interpx,interpy,yy)

    return

  end subroutine geocalc

!-----------------------------------------------------------------------
  subroutine geomat(spsi,xpsi,rap1,xft,rap2,xgmoym)

    implicit none

    real(kindreal), intent(inout):: spsi
    real(kindreal), intent(out):: xpsi,rap1,xft,rap2,xgmoym

    call geocalc(spsi,xpsi,5)   ! xpsi   is the interpolation in rpsi
    call geocalc(spsi,rap1,6)   ! rap1   is the interpolation in rsurf
    call geocalc(spsi,xft,7)    ! xft    is the interpolation in ft
    call geocalc(spsi,rap2,8)   ! rap2   is the interpolation in gpsgmo
    call geocalc(spsi,xgmoym,9) ! xgmoym is the interpolation in gmoym

    return

  end subroutine geomat
!-----------------------------------------------------------------------
  subroutine geomeang(xpsi,xmeang)

    implicit none

    real(kindreal), intent(inout):: xpsi
    real(kindreal), intent(out):: xmeang
    real(kindreal):: xgpsi,xrapg

    call geocalc(xpsi,xgpsi,1)   ! xgpsi is the interpolation in xgnorp
    call geocalc(xpsi,xrapg,2)   ! xrapg is the interpolation in gpsgmo

    xmeang=xgpsi*(1.d0/xrapg)

    return

  end subroutine geomeang
!-----------------------------------------------------------------------
  subroutine geomEdd(eddesc,xoblaMax,OmegaMax)

    implicit none

    real(kindreal), intent(inout):: eddesc
    real(kindreal), intent(out):: xoblaMax,OmegaMax

    call geocalc(eddesc,xoblaMax,10)
    call geocalc(eddesc,OmegaMax,11)

    return

  end subroutine geomEdd

!-----------------------------------------------------------------------
  subroutine geom(xpsi,xfp,xratp,dxfp,dxratp)

    implicit none

    real(kindreal), intent(inout):: xpsi
    real(kindreal), intent(out):: xfp,xratp,dxfp,dxratp

    if (xpsi < rpsi_min) then
      xpsi = rpsi_min
    endif

    call fipoi2(xpsi,n_geo,rpsi,fp,xfp,dxfp)
    call fipoi2(xpsi,n_geo,rpsi,ratp,xratp,dxratp)

    return

  end subroutine geom

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module geomod
