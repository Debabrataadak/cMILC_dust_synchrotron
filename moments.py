import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/debabrata/softwares/PySM_public-master")
from modules import convert_units,noiser,smoother,cmb_map
import pysm
#from pysm  import models
from pysm.nominal import models
import scipy.constants as constants
import math

model="d1s1a2"#"TD"
nside=256
npix=hp.nside2npix(nside)

s = models("s1", nside)
sync=np.zeros((3,npix))
sync[0]=s[0]['A_I']
sync[1]=s[0]['A_Q']
sync[2]=s[0]['A_U']
bitas=s[0]['spectral_index']

s = models("d1", nside)
dust=np.zeros((3,npix))
dust[0]=s[0]['A_I']
dust[1]=s[0]['A_Q']
dust[2]=s[0]['A_U']
bitad=s[0]['spectral_index']
tempd =s[0]['temp']

sync=hp.read_map("./sim/wmap_planck_sync_uk_RJ_model%s_nu0030p00GHz_total_nside0%d.fits"%(model,nside),field=None)
dust=hp.read_map("./sim/wmap_planck_dust_uk_RJ_model%s_nu0353p00GHz_total_nside0%d.fits"%(model,nside),field=None)

frequency = np.array([23,33,41,61,94,30,44,70,100,143,217,353])
#frequency = np.array([30,44,70,100,143,217,353])
no_bands=np.size(frequency)


bita_d=1.53
bita_s = -3.11
Td=19.4
Tcmb=2.7255
xdd = constants.h * 353.0 * 1.e9 / constants.k / Td
###planck function ###
def B(nu):
    xd = constants.h * nu * 1.e9 / constants.k / Td
    return((nu/353.0)**(bita_d + 1.0)*(np.expm1(xdd)/np.expm1(xd)))
def planck(nu,T,bd):
    Temp_d=np.copy(T)
    spec_d=np.copy(bd)
    xd = constants.h * nu * 1.e9 / constants.k / Temp_d
    return((nu/353.0)**(spec_d + 1.0)*(np.expm1(xdd)/np.expm1(xd)))
############

sync_m = np.zeros((no_bands, 3, npix))
dust_m = np.zeros((no_bands, 3, npix))
sync_m_qu = np.zeros((no_bands, 3, npix))
dust_m_qu = np.zeros((no_bands, 3, npix))


###moments upto 2nd order in EB space ####
def synchrotron_moments(bita_s):
	for i in range(no_bands):
		delta_bita_s = bitas-bita_s
		c = (frequency[i]/30.0)**(-3.11)
		d = np.log(frequency[i]/30.0)
		sync_m[i] = np.array(sync)*(np.array(delta_bita_s) +np.array(delta_bita_s)*d*c + 0.5* np.array(delta_bita_s)**2*d**2*c) 
		sync_m_qu[i]=sync_m[i]
		alm=hp.map2alm(sync_m[i])
		for j in range(3):
			sync_m[i][j]=hp.alm2map(alm[j],nside,pol=False,verbose=False)
	return(sync_m,sync_m_qu)

def dust_moments(bita_d,Td):
	for i in range(no_bands):
		delta_bita_d = np.array(bitad-bita_d)
		delta_Td = np.array(tempd - Td)
		c=B(frequency[i])
		d=np.log(frequency[i]/353.0)
		xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
		e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
		f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
		dust_m[i] = np.array(dust)*(delta_bita_d*d*c + delta_Td*e*c + 0.5*delta_bita_d**2*d**2*c + 0.5*delta_Td**2*f*e*c + delta_bita_d*delta_Td*d*e*c)
		alm=hp.map2alm(dust_m[i])
		dust_m_qu[i]=dust_m[i]
		for j in range(3):
			dust_m[i][j]=hp.alm2map(alm[j],nside,pol=False,verbose=False)
		return(dust_m,dust_m_qu)

		
		
