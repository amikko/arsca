#here we have various methods for creating scattering cross-sections
#and scattering phase functions

import netCDF4
import scipy.interpolate
import numpy as np

def rayleigh_xsec(wl,extrapolate=False):
    wl_altius = np.array([300.0,  315.0,  351.0,  435.0,  442.0,  525.0,  600.0,  675.0,  943.0, 1020.0, 1700.0])
    rayleigh_xsec_altius = np.array([5.60283101e-26, 4.54991716e-26, 2.87884699e-26,
                            1.17740547e-26, 1.10228657e-26, 5.43800328e-27,
                            3.15614665e-27, 1.95738772e-27, 5.06993365e-28,
                            3.69241920e-28, 4.61162125e-29])
    rayleigh_interp_fun = scipy.interpolate.interp1d(wl_altius,rayleigh_xsec_altius,kind='linear')
    H = (wl_altius ** -4).reshape((wl_altius.size,1))
    # TODO: Add a weighting[  v--  ]matrix at this point to increase the accuracy on short wavelengths
    coeff = np.linalg.inv(H.T @ H) @ H.T @ (rayleigh_xsec_altius.reshape((rayleigh_xsec_altius.size,1)))
    LS_fun = lambda wl_ : coeff * wl_ ** -4

    m = 1.0001 #huom! wl-dependent!
    a = 1 # huom! tämä pitää korjata kans
    sigma_s = 128 * np.pi / 3 * np.abs((m ** 2 - 1)/(m ** 2 + 2)) * a ** 6 / wl ** 4
    # with some m, a
    if extrapolate:
        return LS_fun(wl)
    else:
        return rayleigh_interp_fun(wl)

def aerosol_xsec(wl,extrapolate=False):
    ds = netCDF4.Dataset('./files_used_in_examples/altius/altius_rtm_comparison_input_v3.nc')
    wl_altius = 1.0 * ds['wavelength'][:].data # the "1.0 *" ensures that the data is float
    aerosol_xsec_altius = ds['aerosol_scattering_cross_section'][:].data
    aerosol_interp_fun = scipy.interpolate.interp1d(wl_altius,aerosol_xsec_altius,kind='quadratic')
    angstrom_expo = 2.2 # this was selected by random checking
    H = (wl_altius ** -angstrom_expo).reshape((wl_altius.size,1))
    # TODO: Add a weighting[  v--  ]matrix at this point to increase the accuracy on short wavelengths
    coeff = np.linalg.inv(H.T @ H) @ H.T @ (aerosol_xsec_altius.reshape((aerosol_xsec_altius.size,1)))
    LS_fun = lambda wl_ : coeff * wl_ ** -4
    if extrapolate:
        return LS_fun(wl)
    else:
        return aerosol_interp_fun(wl)
