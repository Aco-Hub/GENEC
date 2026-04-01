"""Physical constants for GENEC -- ported from src/const.f90.

All values use float (Python float64) to match Fortran real(8).
"""

import math

# ---------------------------------------------------------------------------
# Fundamental constants (CODATA 2006)
# ---------------------------------------------------------------------------
pi = 3.141592653589793
cst_c = 2.99792458e10        # speed of light (cm/s)
cst_G = 6.67428e-8           # gravitational constant (cgs)
cst_h = 6.62606896e-27       # Planck constant (erg*s)
cst_k = 1.3806504e-16        # Boltzmann constant (erg/K)
cst_a = 7.56576738e-15       # radiation constant (erg/cm^3/K^4)
rgaz = 8.314472e7            # gas constant (erg/mol/K)
cst_sigma = 5.67040e-5       # Stefan-Boltzmann constant (cgs)
cst_me = 9.1093826e-28       # electron mass (g)
cst_avo = 6.02214179e23      # Avogadro number
cst_u = 1.660538782e-24      # atomic mass unit (g)
cst_mh = 1.67372346e-24      # proton mass (g)
cst_thomson = 6.652458558e-25  # Thomson cross-section (cm^2)
cst_e = 1.602176487e-19      # electron charge (Coulomb)
cst_ecgs = 4.8032068e-10     # electron charge (esu)
qapicg = 2.514403597e4       # 4*pi*c*G
cst_K1 = 9.844840461e12      # (1/5)*(3/8pi)^(2/3)*h^2/(m_e*h^(-5/3))

# ---------------------------------------------------------------------------
# Log10 of constants
# ---------------------------------------------------------------------------
lgpi = 0.497149873
cstlg_c = 10.4768207
cstlg_G = -7.175595578
cstlg_h = -26.17874405
cstlg_k = -15.85991628
cstlg_a = -14.1210378
rgazlg = 7.919834675
cstlg_sigma = -4.246386304
cstlg_me = -27.04051106
cstlg_avo = 23.77975098
cstlg_u = -23.77975098
cstlg_mh = -23.7763163
cstlg_thomson = -24.17701782
cstlg_e = -18.79528965
lgqapicg = 4.400434989
cstlg_K1 = 1.29932e1

# ---------------------------------------------------------------------------
# Solar parameters
# ---------------------------------------------------------------------------
Msol = 1.9884e33             # solar mass (g)
Rsol = 6.9551e10             # solar radius (cm)
Lsol = 3.8427e33             # solar luminosity (erg/s)
Teffsol = 5.777e3            # solar effective temperature (K)
xlsomo = 1.932558841
uastr = 1.49597870660e13     # astronomical unit (cm)
year = 3.1557600e7           # Julian year (s)
day = 8.640e4                # day (s)
Omega_sol = 2.9e-6           # solar angular velocity (rad/s)

lgMsol = 33.29850375
lgRsol = 10.84229713
lgLsol = 33.58460257

# ---------------------------------------------------------------------------
# Miscellaneous
# ---------------------------------------------------------------------------
um = 2.3025850929940457      # ln(10)

# Nuclear Q-values
Q_H = 26.229
Q_He = 7.274
Q_C = 4.617
convMeVerg = 1.602176487e-6
