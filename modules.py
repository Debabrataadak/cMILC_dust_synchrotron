from __future__ import print_function
import healpy as hp
import numpy as np
import scipy.constants as constants
import scipy.integrate
import sys


def convert_units(unit1, unit2, nu):
    """Function to do unit conversions between Rayleigh-Jeans units, CMB
    units, and flux units.

    :param unit1: unit from which we are converting.
    :type unit1: str.
    :param unit2: unit to which we are converting.
    :type unit2: str.
    :param nu: frequency at which to calculate unit conversion.
    :type nu: float, np.ndarray.
    :returns: unit conversion coefficient - float, numpy.ndarray.
    """
    if "K_CMB" in unit1:
        #first deal with the unit conversion
        if "Jysr" in unit2:
            conversion_factor = K_CMB2Jysr(nu)
        elif "K_RJ" in unit2:
            conversion_factor = K_CMB2Jysr(nu) / K_RJ2Jysr(nu)
        elif "K_CMB" in unit2:
            conversion_factor = np.ones_like(nu)
        else:
            print("Incorrect format or unit.")

    elif "K_RJ" in unit1:
        if "Jysr" in unit2:
            conversion_factor = K_RJ2Jysr(nu)
        elif "K_CMB" in unit2:
            conversion_factor = K_RJ2Jysr(nu) / K_CMB2Jysr(nu)
        elif "K_RJ" in unit2:
            conversion_factor = np.ones_like(nu)
        else:
            print("Incorrect format or unit.")

    elif "Jysr" in unit1:
        if "Jysr" in unit2:
            conversion_factor = np.ones_like(nu)
        elif "K_RJ" in unit2:
            conversion_factor = 1. / K_RJ2Jysr(nu)
        elif "K_CMB" in unit2:
            conversion_factor = 1. / K_CMB2Jysr(nu)
        else:
            print("Incorrect format or unit.")

    # Now deal with the magnitude
    if "n" in unit1[0]:
        prefac = 1.e-9
    elif "u" in unit1[0]:
        prefac = 1.e-6
    elif "m" in unit1[0]:
        prefac = 1.e-3
    elif "k" in unit1[0]:
        prefac = 1.e3
    elif "M" in unit1[0]:
        prefac = 1.e6
    elif "G" in unit1[0]:
        prefac = 1.e9
    elif "K" in unit1[0]:
        prefac = 1.
    elif "J" in unit1[0]:
        prefac = 1.
    else:
        print("Invalid format for unit1 in convert_units")
        sys.exit(1)

    if "n" in unit2[0]:
        postfac = 1.e9
    elif "u" in unit2[0]:
        postfac = 1.e6
    elif "m" in unit2[0]:
        postfac = 1.e3
    elif "k" in unit2[0]:
        postfac = 1.e-3
    elif "M" in unit2[0]:
        postfac = 1.e-6
    elif "G" in unit2[0]:
        postfac = 1.e-9
    elif "K" in unit2[0]:
        postfac = 1.
    elif "J" in unit2[0]:
        postfac = 1.
    else:
        print("Invalid format for unit2 in convert_units")
        sys.exit(1)

    return np.array(conversion_factor * prefac * postfac)

def K_CMB2Jysr(nu):
    """Kelvin_CMB to Janskies per steradian. Nu is in GHz.

    :param nu: frequency in GHz at which to calculate unit conversion.
    :type nu: float, numpy.ndarray
    :return: unit conversion coefficient - float.
    
    """
    return dB(nu, 2.7255) * 1.e26

def K_RJ2Jysr(nu):
    """Kelvin_RJ to Janskies per steradian. Nu is in GHz.

    :param nu: frequency in GHz at which to calculate unit conversion.
    :type nu: float, numpy.ndarray.                                                                                                                      
    :return: unit conversion coefficient - float. 
    
    """
    return  2. * (nu * 1.e9 / constants.c) ** 2 * constants.k * 1.e26

def B(nu, T):
    """Planck function. 

    :param nu: frequency in GHz at which to evaluate planck function.
    :type nu: float.
    :param T: temperature of black body. 
    :type T: float.
    :return: float -- black body brightness.

    """
    x = constants.h * nu * 1.e9 / constants.k / T
    return 2. * constants.h * (nu * 1.e9) ** 3 / constants.c ** 2 / np.expm1(x)

def dB(nu, T):
    """Differential planck function. 
    
    :param nu: frequency in GHz at which to evaluate differential planck function.
    :type nu: float.
    :param T: temperature of black body. 
    :type T: float.
    :return: float -- differential black body function. 

    """
    x = constants.h * nu * 1.e9 / constants.k / T
    return B(nu, T) / T * x * np.exp(x) / np.expm1(x)




def smoother(bFWHM, map_array):
        """Function to smooth an array of N (T, Q, U) maps with N beams in
        units of arcmin.

        :param map_array:
        :type map_array:
        
        """
        b=bFWHM
        m=map_array
        smoothed_map_array = hp.smoothing(m, fwhm = np.pi / 180. * b / 60., verbose = False)
        return(smoothed_map_array)


def noiser(Nside,Sens_I, Sens_P,beam):
        """Calculate white noise maps for given sensitivities. Returns
          noise maps at the given nside in (T, Q, U). Input
        sensitivities are expected to be in uK_CMB amin.
        """
        
        npix = hp.nside2npix(Nside)

        
        # solid angle per pixel in amin2
        pix_amin2 = 4. * np.pi / float(hp.nside2npix(Nside)) * (180. * 60. / np.pi) ** 2
        """sigma_pix_I/P is std of noise per pixel. It is an array of length
        equal to the number of input maps."""
        sigma_pix_I = np.sqrt(Sens_I ** 2 / pix_amin2)
        sigma_pix_P = np.sqrt(Sens_P ** 2 / pix_amin2)
        #np.random.seed(seed = 1245)
        noise = np.zeros((3,npix))
        noise[0]= np.random.normal(loc=0.0, scale=sigma_pix_I, size=npix)
        noise[1]= np.random.normal(loc=0.0, scale=sigma_pix_P, size=npix)
        noise[2]= np.random.normal(loc=0.0, scale=sigma_pix_P, size=npix)
        #smoothed_noise = smoother(beam, noise)    
        #return smoothed_noise
        return noise


def cmb_map(cls,nside,beam):
        """Function for the calculation of lensed CMB maps directly from
        lensed Cls using healpix's synfast routine.
	parameters: CAMB cls, nside and beam in arcm
        """
        # get the spectra. These are in CAMB format, we discard the last
        # three corresponding to dd, dt, de, respectively.
        ell=cls[:,0]
        tt =cls[:,1]
        ee =cls[:,2]
        bb =cls[:,3]
        te =cls[:,4]
        lmax_cl = len(ell) + 1
        ell = np.arange(lmax_cl + 1)

        # in CAMB format so we must divide by the scaling factor
        factor = ell * (ell + 1.) / 2. / np.pi

        cl_teb = np.zeros((6, lmax_cl + 1))
        cl_teb[0, 2:] = tt / factor[2:]
        cl_teb[1, 2:] = ee / factor[2:]
        cl_teb[2, 2:] = bb / factor[2:]
        cl_teb[3, 2:] = te / factor[2:]
        cl_teb[4, 2:] = 0.
        cl_teb[5, 2:] = 0.

        #np.random.seed(1000)
        T, Q, U = hp.synfast(cl_teb, nside, pol=True, fwhm=np.radians(beam/60.), new=True, verbose=False)
        model=np.array([T,Q,U])

        return model

