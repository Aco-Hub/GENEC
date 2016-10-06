module const

  use evol, only: kindreal

  implicit none

  real(kindreal), parameter:: pi=3.141592653589793d0,cst_c=2.99792458d10,cst_G=6.67428d-8,&
     cst_h=6.62606896d-27,cst_k=1.3806504d-16,cst_a=7.56576738d-15,rgaz=8.314472d7,cst_sigma=5.67040d-5,&
     cst_me=9.1093826d-28,cst_avo=6.02214179d23,cst_u=1.660538782d-24,cst_mh=1.67372346d-24,&
     cst_thomson=6.652458558d-25,&
! charge de l'electron en coulomb.
     cst_e=1.602176487d-19,cst_ecgs=4.8032068d-10,&
! qapicg: 4*pi*c*G
     qapicg=2.514403597d4,&
! K1= (1/5)*(3/8pi)**(2/3)*h**2/(m_e*h**(-5/3))
     cst_K1=9.844840461d12

  real(kindreal), parameter:: lgpi=0.497149873d0,cstlg_c=10.4768207d0,cstlg_G=-7.175595578d0,&
     cstlg_h=-26.17874405d0, cstlg_k=-15.85991628d0,cstlg_a=-14.1210378d0,rgazlg=7.919834675d0,&
     cstlg_sigma=-4.246386304d0,cstlg_me=-27.04051106d0,cstlg_avo=23.77975098d0,cstlg_u=-23.77975098d0,&
     cstlg_mh=-23.7763163d0,cstlg_thomson=-24.17701782d0,cstlg_e=-18.79528965,lgqapicg=4.400434989,&
! K1= (1/5)*(3/8pi)**(2/3)*h**2/(m_e*h**(-5/3))
     cstlg_K1=1.29932d1

  real(kindreal), parameter:: Msol=1.9884d33,Rsol=6.9551d10,Lsol=3.8427d33,xlsomo=1.932558841d0,&
     uastr=1.49597870660d13,year=3.15569252d7,day=8.640d4,Omega_sol=2.9d-6

  real(kindreal), parameter:: lgMsol=33.29850375d0,lgRsol=10.84229713d0,lgLsol=33.58460257d0

  real(kindreal), parameter:: um=2.3025850929940457d0
      
  real(kindreal), parameter:: Q_H=26.229d0,Q_He=7.274d0,Q_C=4.617d0,convMeVerg=1.602176487d-6

end module const
