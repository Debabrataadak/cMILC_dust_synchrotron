import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/debabrata/softwares/PySM_public-master")
from modules import convert_units,noiser,smoother,cmb_map
import pysm
from pysm.nominal import models
import scipy.constants as constants
import math
from dust_mixing_matrix_data_RJ import Mixingmatrix
#from dust_mixing_matrix_curve_sync_RJ import Mixingmatrix
import sys
import os
import moments 
from moments import dust_moments


model=str(sys.argv[3])
d_d=str(sys.argv[4])
s_d=str(sys.argv[5])
#model_ids(d_id,s_id,model)
moments.d_id=d_d
moments.s_id=s_d
moments.model_id=model
#x,y,z=moments.moment_ids()
#print(x,y,z)

model1=model

dir1="pdfs/wmap_planck_pysm/dust/%s"%(model)
dir2="maps/wmap_planck_pysm/dust/%s"%(model)
os.system( 'mkdir '+dir1)
os.system( 'mkdir '+dir2)


nside=256
npix=hp.nside2npix(nside)
#cmb = hp.read_map("./sim/wmap_planck_cmb_uk_RJ_nu0100p00GHz_total_nside0256.fits",field=None)
#cmb = hp.read_map("./sim/wmap_planck_cmb_uk_RJ_model%s_nu0100p00GHz_total_nside0256.fits"%(model),field=None)
sync=hp.read_map("./sim/wmap_planck_sync_uk_RJ_model%s_nu0030p00GHz_total_nside0256.fits"%(model),field=None)
dust=hp.read_map("./sim/wmap_planck_dust_uk_RJ_model%s_nu0353p00GHz_total_nside0256.fits"%(model),field=None)
###prepare mask ####
print("here")

msk = np.ones(hp.nside2npix(nside))
"""
msk = np.zeros(hp.nside2npix(nside))
th, ph = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
ph[np.where(ph > np.pi)[0]] -= 2 * np.pi
msk[np.where((th < 2.63) & (th > 1.86) &
             (ph > -np.pi / 4) & (ph < np.pi / 4))[0]] = 1.

##ns=256 mask
if(nside==256):
	msk=hp.read_map("planck_data/processed_gal_mask_fsky79_ns%d.fits"%nside)
#ns=64 mask
if(nside==64):
	msk=hp.read_map("planck_data/processed_gal_ps_mask_fsky73_ns%d.fits"%nside)
#hp.mollview(msk)
#plt.show()
"""
ind=np.where(msk==1)
npix1=np.size(ind)

##### rest of the performance will be over the unmasked regions ####

bita_d=1.53
Td=19.4
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
frequency = np.array([23,30,33,41,44,61,70,94,100,143,217,353])
#frequency = np.array([23,33,41,61,94,30,44,70,100,143,217,353])
#frequency = np.array([30,44,70,100,143,217,353])
no_bands=np.size(frequency)


##cmb in qu and EB space
maps=np.zeros((no_bands,3,npix),dtype=np.float64)
maps_qu=np.zeros((no_bands,3,npix),dtype=np.float64)
maps3=np.zeros((no_bands,3,npix),dtype=np.float64)
maps3_qu=np.zeros((no_bands,3,npix),dtype=np.float64)

"""
model="ffp10"
for i in range (no_bands):
    if(frequency[i]<100):
        maps3_qu[i]=hp.read_map("./sim/wmap_planck_total1_uk_RJ_model%s_nu00%dp00GHz_total_nside0256.fits"%(model,frequency[i]),field=None)
        maps3[i]=hp.read_map("./sim/wmap_planck_total1_uk_RJ_model%s_nu00%dp00GHz_total_nside0256.fits"%(model,frequency[i]),field=None)
    else:
        maps3_qu[i]=hp.read_map("./sim/wmap_planck_total1_uk_RJ_model%s_nu0%dp00GHz_total_nside0256.fits"%(model,frequency[i]),field=None)
        maps3[i]=hp.read_map("./sim/wmap_planck_total1_uk_RJ_model%s_nu0%dp00GHz_total_nside0256.fits"%(model,frequency[i]),field=None)

    alm=hp.map2alm(maps3[i])
    for j in range(3):
        maps3[i][j]=hp.alm2map(alm[j],nside,pol=False,verbose=False)

model=model1

"""
for i in range (no_bands):
    if(frequency[i]<100):
        maps3_qu[i]=hp.read_map("./sim/wmap_planck_cmb_uk_RJ_nu00%dp00GHz_total_nside0256.fits"%(frequency[i]),field=None)
        maps3[i]=hp.read_map("./sim/wmap_planck_cmb_uk_RJ_nu00%dp00GHz_total_nside0256.fits"%(frequency[i]),field=None)
    else:
        maps3_qu[i]=hp.read_map("./sim/wmap_planck_cmb_uk_RJ_nu0%dp00GHz_total_nside0256.fits"%(frequency[i]),field=None)
        maps3[i]=hp.read_map("./sim/wmap_planck_cmb_uk_RJ_nu0%dp00GHz_total_nside0256.fits"%(frequency[i]),field=None)
    alm=hp.map2alm(maps3[i])
    for j in range(3):
        maps3[i][j]=hp.alm2map(alm[j],nside,pol=False,verbose=False)
"""
for i in range (no_bands):
    if(frequency[i]<100):
        maps3_qu[i]=hp.read_map("./sim/wmap_planck_cmb_uk_RJ_model%s_nu00%dp00GHz_total_nside0256.fits"%(model,frequency[i]),field=None)
        maps3[i]=hp.read_map("./sim/wmap_planck_cmb_uk_RJ_model%s_nu00%dp00GHz_total_nside0256.fits"%(model,frequency[i]),field=None)
    else:
        maps3_qu[i]=hp.read_map("./sim/wmap_planck_cmb_uk_RJ_model%s_nu0%dp00GHz_total_nside0256.fits"%(model,frequency[i]),field=None)
        maps3[i]=hp.read_map("./sim/wmap_planck_cmb_uk_RJ_model%s_nu0%dp00GHz_total_nside0256.fits"%(model,frequency[i]),field=None)
    alm=hp.map2alm(maps[i])
    for j in range(3):
        maps3[i][j]=hp.alm2map(alm[j],nside,pol=False,verbose=False)
"""
###synchrotron in qu and EB space
maps1=np.zeros((no_bands,3,npix),dtype=np.float64)
maps1_qu=np.zeros((no_bands,3,npix),dtype=np.float64)
for i in range (no_bands):
    if(frequency[i]<100):
        maps1_qu[i]=hp.read_map("./sim/wmap_planck_sync_uk_RJ_model%s_nu00%dp00GHz_total_nside0256.fits"%(model,frequency[i]),field=None)
        maps1[i]=hp.read_map("./sim/wmap_planck_sync_uk_RJ_model%s_nu00%dp00GHz_total_nside0256.fits"%(model,frequency[i]),field=None)
    else:
        maps1_qu[i]=hp.read_map("./sim/wmap_planck_sync_uk_RJ_model%s_nu0%dp00GHz_total_nside0256.fits"%(model,frequency[i]),field=None)
        maps1[i]=hp.read_map("./sim/wmap_planck_sync_uk_RJ_model%s_nu0%dp00GHz_total_nside0256.fits"%(model,frequency[i]),field=None)
    
    alm=hp.map2alm(maps1[i])
    for j in range(3):
        maps1[i][j]=hp.alm2map(alm[j],nside,pol=False,verbose=False)

#######Dust in qu and EB space ####
maps2=np.zeros((no_bands,3,npix),dtype=np.float64)
maps2_qu=np.zeros((no_bands,3,npix),dtype=np.float64)
for i in range (no_bands):
    if(frequency[i]<100):
        maps2_qu[i]=hp.read_map("./sim/wmap_planck_dust_uk_RJ_model%s_nu00%dp00GHz_total_nside0256.fits"%(model,frequency[i]),field=None)
        maps2[i]=hp.read_map("./sim/wmap_planck_dust_uk_RJ_model%s_nu00%dp00GHz_total_nside0256.fits"%(model,frequency[i]),field=None)
    else:
        maps2_qu[i]=hp.read_map("./sim/wmap_planck_dust_uk_RJ_model%s_nu0%dp00GHz_total_nside0256.fits"%(model,frequency[i]),field=None)
        maps2[i]=hp.read_map("./sim/wmap_planck_dust_uk_RJ_model%s_nu0%dp00GHz_total_nside0256.fits"%(model,frequency[i]),field=None)

    alm=hp.map2alm(maps2[i])
    for j in range(3):
        maps2[i][j]=hp.alm2map(alm[j],nside,pol=False,verbose=False)

#######spin-Dust in qu and EB space ####
maps4=np.zeros((no_bands,3,npix),dtype=np.float64)
maps4_qu=np.zeros((no_bands,3,npix),dtype=np.float64)
for i in range (no_bands):
    if(frequency[i]<100):
        maps4_qu[i]=hp.read_map("./sim/wmap_planck_ame_uk_RJ_model%s_nu00%dp00GHz_total_nside0256.fits"%(model,frequency[i]),field=None)
        maps4[i]=hp.read_map("./sim/wmap_planck_ame_uk_RJ_model%s_nu00%dp00GHz_total_nside0256.fits"%(model,frequency[i]),field=None)
    else:
        maps4_qu[i]=hp.read_map("./sim/wmap_planck_ame_uk_RJ_model%s_nu0%dp00GHz_total_nside0256.fits"%(model,frequency[i]),field=None)
        maps4[i]=hp.read_map("./sim/wmap_planck_ame_uk_RJ_model%s_nu0%dp00GHz_total_nside0256.fits"%(model,frequency[i]),field=None)

    alm=hp.map2alm(maps4[i])
    for j in range(3):
        maps4[i][j]=hp.alm2map(alm[j],nside,pol=False,verbose=False)


##noise in qu and EB space
simdir    = "./sim"
noise=np.zeros((no_bands,3,npix),dtype=np.float64)
noise_qu=np.zeros((no_bands,3,npix),dtype=np.float64)
for i in range (no_bands):
    noise_qu[i]=hp.read_map("%s/wmap_planck_sim_noise_%dGHz_ns%d_uk_RJ.fits"%(simdir,int(frequency[i]),nside),field=None,verbose=False)
    noise[i]=hp.read_map("%s/wmap_planck_sim_noise_%dGHz_ns%d_uk_RJ.fits"%(simdir,int(frequency[i]),nside),field=None,verbose=False)
    alm=hp.map2alm(noise[i])
    for j in range(3):
        noise[i][j]=hp.alm2map(alm[j],nside,pol=False,verbose=False)

###total maps
for i in range (no_bands):
    maps[i] =  noise[i] + maps3[i] + maps1[i] + maps2[i] + maps4[i]
    maps_qu[i] =  noise_qu[i] + maps3_qu[i] + maps1_qu[i] + maps2_qu[i] + maps4_qu[i]


### covariance in EB/qu space 
cov_matrix=np.zeros((3,no_bands,no_bands),)
cov_matrix_qu=np.zeros((3,no_bands,no_bands))
for k in range(3):
    for i in range(no_bands):
        for j in range (no_bands):
            if (i<=j):
                ##eb
                cov=sum((maps[i,k][ind]-np.mean(maps[i,k][ind]))*(maps[j,k][ind]-np.mean(maps[j,k][ind])))/npix1
                cov_matrix[k,i,j] = cov
                ##qu
                cov=sum((maps_qu[i,k][ind]-np.mean(maps_qu[i,k][ind]))*(maps_qu[j,k][ind]-np.mean(maps_qu[j,k][ind])))/npix1
                cov_matrix_qu[k,i,j] = cov
            else:
                cov_matrix[k,i,j] = cov_matrix[k,j,i]
                cov_matrix_qu[k,i,j] = cov_matrix_qu[k,j,i]




inv_cov_matrix=np.zeros((3,no_bands,no_bands))
inv_cov_matrix_qu=np.zeros((3,no_bands,no_bands))
for k in range(3):
    inv_cov_matrix[k] = np.linalg.inv(cov_matrix[k])
    inv_cov_matrix_qu[k] = np.linalg.inv(cov_matrix_qu[k])

#### constrain matrix for CILC ####
ncrr=str(sys.argv[1])
ncr=int(sys.argv[2])
F=Mixingmatrix(ncrr) 
##################################
###
e=np.zeros(ncr)
e[0]=1
wc=np.zeros((3,no_bands),dtype=np.float64)
for i in range(3):
    x=np.zeros(ncr,dtype=np.float64)
    y=np.zeros((ncr,no_bands),dtype=np.float64)
    x=np.dot(F.T,np.dot(inv_cov_matrix[i],F))
    y=np.dot(F.T,inv_cov_matrix[i])
    #print(F.T)
    #print(np.linalg.det(x))
    #print(x)
    wc[i] = np.dot(e,np.dot(np.linalg.inv(x),y))

#####
a=np.array([B(frequency[i]) for i in range(no_bands)])
w=np.zeros((3,no_bands),dtype=np.float64)
w_qu=np.zeros((3,no_bands),dtype=np.float64)
for i in range(3):
    w[i]= np.dot(a,inv_cov_matrix[i])/np.dot(a,np.dot(inv_cov_matrix[i],a))
    w_qu[i]= np.dot(a,inv_cov_matrix_qu[i])/np.dot(a,np.dot(inv_cov_matrix_qu[i],a))


######PRILC method ###
cov_matrix_P=np.zeros((no_bands,no_bands),dtype=np.float64)
for i in range(no_bands):
        for j in range (no_bands):
            if (i<=j):
                ##P
                cov=sum(maps_qu[i,1][ind]*maps_qu[j,1][ind] + maps_qu[i,2][ind]*maps_qu[j,2][ind])/npix1
                cov_matrix_P[i,j] = cov
            else:
	        cov_matrix_P[i,j] = cov_matrix_P[j,i]
inv_cov_matrix_P=np.zeros((no_bands,no_bands),dtype=np.float64)
inv_cov_matrix_P = np.linalg.inv(cov_matrix_P)
w_P = np.dot(a,inv_cov_matrix_P)/np.dot(a,np.dot(inv_cov_matrix_P,a))	
w_P = np.ones((3,no_bands))*w_P
###PCILC ###
x=np.zeros(ncr,dtype=np.float64)
y=np.zeros((ncr,no_bands),dtype=np.float64)
x=np.dot(F.T,np.dot(inv_cov_matrix_P,F))
y=np.dot(F.T,inv_cov_matrix_P)
wpc = np.dot(e,np.dot(np.linalg.inv(x),y))
wpc = np.ones((3,no_bands))*wpc



####
"""
print("whether condition holds for ebILC")
print(np.dot(a,w[0]),np.dot(a,w[1]),np.dot(a,w[2]))
print("whether condition holds for quILC")
print(np.dot(a,w_qu[0]),np.dot(a,w_qu[1]),np.dot(a,w_qu[2]))
print("whether condition holds for pILC")
print(np.dot(a,w_P[0]),np.dot(a,w_P[1]),np.dot(a,w_P[2]))
print("whether condition holds for CILC")
print(np.dot(wc[0],F[:,0]),np.dot(wc[1],F[:,0]),np.dot(wc[2],F[:,0]))
print(np.dot(wc[0],F[:,1]),np.dot(wc[1],F[:,1]),np.dot(wc[2],F[:,1]))
print("PCILC")
print(np.dot(wpc[0],F[:,0]),np.dot(wpc[1],F[:,0]),np.dot(wpc[2],F[:,0]))
print(np.dot(wpc[0],F[:,1]),np.dot(wpc[1],F[:,1]),np.dot(wpc[2],F[:,1]))
"""
for i in range(int(ncr)):
        print(i)
        print(np.dot(wpc[0],F[:,i]),np.dot(wpc[1],F[:,i]),np.dot(wpc[2],F[:,i]))
##### recovered maps and residual maps ####

s_qu=np.zeros((3,npix),dtype=np.float64)
s=np.zeros((3,npix),dtype=np.float64)

sn_qu=np.zeros((3,npix),dtype=np.float64)
sn=np.zeros((3,npix),dtype=np.float64)

sf_qu=np.zeros((3,npix),dtype=np.float64)
sf=np.zeros((3,npix),dtype=np.float64)


s_P=np.zeros((3,npix),dtype=np.float64)
sn_P=np.zeros((3,npix),dtype=np.float64)
sf_P=np.zeros((3,npix),dtype=np.float64)

s_c=np.zeros((3,npix),dtype=np.float64)
sn_c=np.zeros((3,npix),dtype=np.float64)
sf_c=np.zeros((3,npix),dtype=np.float64)


s_pc=np.zeros((3,npix),dtype=np.float64)
sn_pc=np.zeros((3,npix),dtype=np.float64)
sf_pc=np.zeros((3,npix),dtype=np.float64)


map_moments = np.zeros((no_bands, 3, npix),dtype=np.float64)
map_moments_qu = np.zeros((no_bands, 3, npix),dtype=np.float64)
map_moments,map_moments_qu = dust_moments(bita_d,Td)
for k in range (3):
    for i in range(no_bands):
        ##recovered sync
        s[k]=s[k]+np.array(maps[i,k,:])*w[k,i]
        s_qu[k]=s_qu[k]+np.array(maps_qu[i,k,:])*w_qu[k,i]
        s_P[k]=s_P[k]+np.array(maps_qu[i,k,:])*w_P[k,i]
        s_c[k]=s_c[k]+np.array(maps[i,k,:])*wc[k,i]
        s_pc[k]=s_pc[k]+np.array(maps_qu[i,k,:])*wpc[k,i]
	##NOISE RESIDUAL
	sn[k]=sn[k]+np.array(noise[i,k,:])*w[k,i]
        sn_qu[k]=sn_qu[k]+np.array(noise_qu[i,k,:])*w_qu[k,i]
        sn_P[k]=sn_P[k]+np.array(noise_qu[i,k,:])*w_P[k,i]
        sn_c[k]=sn_c[k]+np.array(noise[i,k,:])*wc[k,i]
        sn_pc[k]=sn_pc[k]+np.array(noise_qu[i,k,:])*wpc[k,i]
        #foreground residual 
        sf[k]=sf[k]+np.array(maps[i,k,:] + maps1[i,k,:])*w[k,i]
        sf_qu[k]=sf_qu[k]+np.array(maps_qu[i,k,:] + maps1_qu[i,k,:])*w_qu[k,i]
        sf_P[k]=sf_P[k]+np.array(maps_qu[i,k,:] + maps1_qu[i,k,:])*w_P[k,i]
        sf_c[k]=sf_c[k]+np.array(maps[i,k,:] + maps1[i,k,:] + map_moments[i,k,:])*wc[k,i]
        sf_pc[k]=sf_pc[k]+np.array(maps3_qu[i,k,:] + map_moments_qu[i,k,:]+ maps1_qu[i,k,:] + maps4_qu[i,k,:])*wpc[k,i]

###input dust QU-EB transfer ##
dust_eb=np.zeros((3,hp.nside2npix(nside)))
alm=hp.map2alm(dust,pol=True)
for j in range(3):
	dust_eb[j]=hp.alm2map(alm[j],nside,pol=False,verbose=False)
for k in range(3):
	sf_c[k]=s_c[k] - dust_eb[k] - sn_c[k]
	sf_pc[k]=s_pc[k] - dust[k] - sn_pc[k]
#### qu -->eb transfromation of the recovered maps ####

"""
alm=hp.map2alm(s_qu)
almn=hp.map2alm(sn_qu)
almf=hp.map2alm(sf_qu)
almP=hp.map2alm(s_P)
almnP=hp.map2alm(sn_P)
almfP=hp.map2alm(sf_P)
almPc=hp.map2alm(s_pc)
almnPc=hp.map2alm(sn_pc)
almfPc=hp.map2alm(sf_pc)
almc=hp.map2alm(cmb)
alms=hp.map2alm(sync)
almd=hp.map2alm(dust)
cmb=np.zeros((3,npix))
sync=np.zeros((3,npix))
dust=np.zeros((3,npix))
for j in range(3):
	s_qu[j]=hp.alm2map(alm[j],nside,pol=False,verbose=False)
	sn_qu[j]=hp.alm2map(almn[j],nside,pol=False,verbose=False)
	sf_qu[j]=hp.alm2map(almf[j],nside,pol=False,verbose=False)
	cmb[j]=hp.alm2map(almc[j],nside,pol=False,verbose=False)
	sync[j]=hp.alm2map(alms[j],nside,pol=False,verbose=False)
	dust[j]=hp.alm2map(almd[j],nside,pol=False,verbose=False)
	s_P[j]=hp.alm2map(almP[j],nside,pol=False,verbose=False)
        sn_P[j]=hp.alm2map(almnP[j],nside,pol=False,verbose=False)
        sf_P[j]=hp.alm2map(almfP[j],nside,pol=False,verbose=False)	
        s_pc[j]=hp.alm2map(almPc[j],nside,pol=False,verbose=False)
        sn_pc[j]=hp.alm2map(almnPc[j],nside,pol=False,verbose=False)
        sf_pc[j]=hp.alm2map(almfPc[j],nside,pol=False,verbose=False) 
"""

##eb --->QU
alm=hp.map2alm(s,pol=False)
almn=hp.map2alm(sn,pol=False)
almf=hp.map2alm(sf,pol=False)
almc=hp.map2alm(s_c,pol=False)
almnc=hp.map2alm(sn_c,pol=False)
almfc=hp.map2alm(sf_c,pol=False)
#for j in range(3):
        
s=hp.alm2map(alm,nside,pol=True,verbose=False)
sn=hp.alm2map(almn,nside,pol=True,verbose=False)
sf=hp.alm2map(almf,nside,pol=True,verbose=False)
s_c=hp.alm2map(almc,nside,pol=True,verbose=False)
sn_c=hp.alm2map(almnc,nside,pol=True,verbose=False)
sf_c=hp.alm2map(almfc,nside,pol=True,verbose=False)

####plotting
fig = plt.figure(figsize = (12, 8))
hp.mollview(s[0]*msk,cmap=plt.cm.jet,sub = (131),title="ILC_T")
hp.mollview(s[1]*msk,cmap=plt.cm.jet,sub = (132),title="ILC_Q")
hp.mollview(s[2]*msk,cmap=plt.cm.jet,min=None,max=None,sub = (133),title="ILC_U")
plt.savefig("pdfs/wmap_planck_pysm/dust/%s/diff_dust_ebilc_ns%d.pdf"%(model,nside))
hp.write_map("maps/wmap_planck_pysm/dust/%s/dust_ebcilc_ns%d.fits"%(model,nside),s*msk)
"""
fig = plt.figure(figsize = (12, 8))
hp.mollview(s_qu[0]-cmb[0],cmap=plt.cm.jet,sub = (131),title="ILC_T")
hp.mollview(s_qu[1]-cmb[1],cmap=plt.cm.jet,sub = (132),title="ILC_Q")
hp.mollview(s_qu[2]-cmb[2],cmap=plt.cm.jet,sub = (133),title="ILC_U")
"""

fig = plt.figure(figsize = (12, 8))
hp.mollview(s_P[0]*msk,cmap=plt.cm.jet,sub = (131),title="ILC_T")
hp.mollview(s_P[1]*msk,cmap=plt.cm.jet,sub = (132),title="ILC_Q")
hp.mollview(s_P[2]*msk,cmap=plt.cm.jet,min=None,max=None,sub = (133),title="ILC_U")
plt.savefig("pdfs/wmap_planck_pysm/dust/%s/diff_dust_pilc_ns%d.pdf"%(model,nside))
hp.write_map("maps/wmap_planck_pysm/dust/%s/dust_pilc_ns%d.fits"%(model,nside),s_P*msk)

fig = plt.figure(figsize = (12, 8))
hp.mollview(s_c[0]*msk,cmap=plt.cm.jet,sub = (131),title="ILC_T")
hp.mollview(s_c[1]*msk,cmap=plt.cm.jet,sub = (132),title="ILC_Q")
hp.mollview(s_c[2]*msk,cmap=plt.cm.jet,min=None,max=None,sub = (133),title="ILC_U")
plt.savefig("pdfs/wmap_planck_pysm/dust/%s/diff_dust_cilc_nc%s_ns%d.pdf"%(model,ncrr,nside))
hp.write_map("maps/wmap_planck_pysm/dust/%s/dust_cilc_nc%s_ns%d.fits"%(model,ncrr,nside),s_c*msk)


fig = plt.figure(figsize = (12, 8))
hp.mollview(s_pc[0]*msk,cmap=plt.cm.jet,sub = (131),title="ILC_T")
hp.mollview(s_pc[1]*msk,cmap=plt.cm.jet,sub = (132),title="ILC_Q")
hp.mollview(s_pc[2]*msk,min=None,max=None,cmap=plt.cm.jet,sub = (133),title="ILC_U")
plt.savefig("pdfs/wmap_planck_pysm/dust/%s/diff_dust_pcilc_nc%s_ns%d.pdf"%(model,ncrr,nside))
hp.write_map("maps/wmap_planck_pysm/dust/%s/dust_pcilc_nc%s_ns%d.fits"%(model,ncrr,nside),s_pc*msk)

fig = plt.figure(figsize = (12, 8))
hp.mollview(dust[0]*msk,cmap=plt.cm.jet,sub = (131),title="ILC_T")
hp.mollview(dust[1]*msk,cmap=plt.cm.jet,sub = (132),title="ILC_Q")
hp.mollview(dust[2]*msk,min=None,max=None,cmap=plt.cm.jet,sub = (133),title="ILC_U")
plt.savefig("pdfs/wmap_planck_pysm/dust/%s/input_dust_ns%d.pdf"%(model,nside))
hp.write_map("maps/wmap_planck_pysm/dust/%s/input_dust_ns%d.fits"%(model,nside),s_pc*msk)

####fg residual ###
fig = plt.figure(figsize = (12, 8))
hp.mollview(sf_c[0]*msk,cmap=plt.cm.jet,sub = (131),title="ILC_T")
hp.mollview(sf_c[1]*msk,cmap=plt.cm.jet,sub = (132),title="ILC_Q",min=-0.02,max=0.02)
hp.mollview(sf_c[2]*msk,cmap=plt.cm.jet,min=-0.02,max=0.02,sub = (133),title="ILC_U")
plt.savefig("pdfs/wmap_planck_pysm/dust/%s/diff_dust_fgrs_cilc_nc%s_ns%d.pdf"%(model,ncrr,nside))
hp.write_map("maps/wmap_planck_pysm/dust/%s/dust_fgrs_cilc_nc%s_ns%d.fits"%(model,ncrr,nside),sf_c*msk)


fig = plt.figure(figsize = (12, 8))
hp.mollview(sf_pc[0]*msk,cmap=plt.cm.jet,sub = (131),title="ILC_T")
hp.mollview(sf_pc[1]*msk,cmap=plt.cm.jet,sub = (132),title="ILC_Q",min=-0.02,max=0.02)
hp.mollview(sf_pc[2]*msk,min=-0.02,max=0.02,cmap=plt.cm.jet,sub = (133),title="ILC_U")
plt.savefig("pdfs/wmap_planck_pysm/dust/%s/diff_dust_fgrs_pcilc_nc%s_ns%d.pdf"%(model,ncrr,nside))
hp.write_map("maps/wmap_planck_pysm/dust/%s/dust_fgrs_pcilc_nc%s_ns%d.fits"%(model,ncrr,nside),sf_pc*msk)

###noise residual ###
fig = plt.figure(figsize = (12, 8))
hp.mollview(sn_c[0]*msk,cmap=plt.cm.jet,sub = (131),title="ILC_T")
hp.mollview(sn_c[1]*msk,cmap=plt.cm.jet,min=-0.7,max=0.7,sub = (132),title="ILC_Q")
hp.mollview(sn_c[2]*msk,cmap=plt.cm.jet,min=-0.7,max=0.7,sub = (133),title="ILC_U")
plt.savefig("pdfs/wmap_planck_pysm/dust/%s/diff_dust_nrs_cilc_nc%s_ns%d.pdf"%(model,ncrr,nside))
hp.write_map("maps/wmap_planck_pysm/dust/%s/dust_nrs_cilc_nc%s_ns%d.fits"%(model,ncrr,nside),sn_c*msk)


fig = plt.figure(figsize = (12, 8))
hp.mollview(sn_pc[0]*msk,cmap=plt.cm.jet,sub = (131),title="ILC_T")
hp.mollview(sn_pc[1]*msk,cmap=plt.cm.jet,min=-0.7,max=0.7,sub = (132),title="ILC_Q")
hp.mollview(sn_pc[2]*msk,min=-0.7,max=0.7,cmap=plt.cm.jet,sub = (133),title="ILC_U")
plt.savefig("pdfs/wmap_planck_pysm/dust/%s/diff_dust_nrs_pcilc_nc%s_ns%d.pdf"%(model,ncrr,nside))
hp.write_map("maps/wmap_planck_pysm/dust/%s/dust_nrs_pcilc_nc%s_ns%d.fits"%(model,ncrr,nside),sn_pc*msk)

########## weights ####
fig = plt.figure(figsize = (12, 8))
plt.plot(frequency,wpc[1],color="b",label="Q")
plt.plot(frequency,wpc[2],color="g",label="U")
plt.legend()
plt.savefig("pdfs/wmap_planck_pysm/dust/%s/Pcilc_weights_QU_nc%s_ns%d.pdf"%(model,ncrr,nside))
np.savetxt("maps/wmap_planck_pysm/dust/%s/Pcilc_weights_QU_nc%s_ns%d.fits"%(model,ncrr,nside),wpc)


#plt.show()


###comapre with inputs ###

## for B-mode
cmb_sigma=np.std(dust[2])
nbins=100
plt.figure()
plt.title('B hist')
nx2, xbins2, ptchs = plt.hist(s[2][ind], bins=nbins, range=[-60,90])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_s=xbins2[0:nbins] + 0.5*width2
yd_s=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='EB space')

nx2, xbins2, ptchs = plt.hist(s_qu[2][ind], bins=nbins, range=[-60,90])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_qu=xbins2[0:nbins] + 0.5*width2
yd_qu=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='QU space')

nx2, xbins2, ptchs = plt.hist(s_P[2][ind], bins=nbins, range=[-60,90])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_p=xbins2[0:nbins] + 0.5*width2
yd_p=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='P space')

nx2, xbins2, ptchs = plt.hist(dust[2][ind], bins=nbins, range=[-60,90])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_in=xbins2[0:nbins] + 0.5*width2
yd_in=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='input')

nx2, xbins2, ptchs = plt.hist(s_c[2][ind], bins=nbins, range=[-60,90])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_c=xbins2[0:nbins] + 0.5*width2
yd_c=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='CILC')


nx2, xbins2, ptchs = plt.hist(s_pc[2][ind], bins=nbins, range=[-60,90])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_pc=xbins2[0:nbins] + 0.5*width2
yd_pc=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='pCILC')



plt.figure()
plt.title('B hist')
plt.plot(xd_s,yd_s,linewidth=2,label='EB space')
plt.plot(xd_qu,yd_qu,linewidth=2,label='QU space')
plt.plot(xd_p,yd_p,linewidth=2,label='P space')
plt.plot(xd_c,yd_c,linewidth=2,label='cILC')
plt.plot(xd_pc,yd_pc,linewidth=2,label='pcILC')
plt.plot(xd_in,yd_in,linewidth=2,label='input')
plt.legend()
plt.savefig("pdfs/wmap_planck_pysm/dust/%s/hist_dust_Bmode_nc%s_ns%d.pdf"%(model,ncrr,nside))
#plt.show()
np.savetxt("maps/wmap_planck_pysm/dust/%s/hist_dust_Bmode_nc%s_ns%d.txt"%(model,ncrr,nside),np.array([xd_s,yd_s,xd_qu,yd_qu,xd_p,yd_p,xd_c,yd_c,xd_pc,yd_pc,xd_in,yd_in]))
### B-mode residuals ###
##fg
cmb_sigma=1.0
nbins=100
plt.figure()
plt.title('B fg residual hist')
nx2, xbins2, ptchs = plt.hist(sf[2][ind]/cmb_sigma, bins=nbins, range=[-6,6])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_s=xbins2[0:nbins] + 0.5*width2
yd_s=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='EB space')

nx2, xbins2, ptchs = plt.hist(sf_qu[2][ind]/cmb_sigma, bins=nbins,range=[-6,6])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_qu=xbins2[0:nbins] + 0.5*width2
yd_qu=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='QU space')

nx2, xbins2, ptchs = plt.hist(sf_P[2][ind]/cmb_sigma, bins=nbins,range=[-6,6])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_p=xbins2[0:nbins] + 0.5*width2
yd_p=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='P space')

nx2, xbins2, ptchs = plt.hist(sf_pc[2][ind]/cmb_sigma, bins=nbins,range=[-0.01,0.01])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_pc=xbins2[0:nbins] + 0.5*width2
yd_pc=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='pCILC')

nx2, xbins2, ptchs = plt.hist(sf_c[2][ind]/cmb_sigma, bins=nbins,range=[-0.01,0.01])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_c=xbins2[0:nbins] + 0.5*width2
yd_c=nx2 /np.sum(nx2)


plt.figure()
plt.title('B fg residual hist')
plt.plot(xd_s,yd_s,linewidth=2,label='EB space')
plt.plot(xd_qu,yd_qu,linewidth=2,label='QU space')
plt.plot(xd_p,yd_p,linewidth=2,label='P space')
plt.plot(xd_c,yd_c,linewidth=2,label='cILC')
plt.plot(xd_pc,yd_pc,linewidth=2,label='pcILC')
plt.legend()
plt.savefig("pdfs/wmap_planck_pysm/dust/%s/hist_fg_Bmode_nc%s_ns%d.pdf"%(model,ncrr,nside))
#plt.show()
np.savetxt("maps/wmap_planck_pysm/dust/%s/hist_fg_Bmode_nc%s_ns%d.txt"%(model,ncrr,nside),np.array([xd_s,yd_s,xd_qu,yd_qu,xd_p,yd_p,xd_c,yd_c,xd_pc,yd_pc]))


###noise residual####
nbins=100
nx2, xbins2, ptchs = plt.hist(sn[2][ind]/cmb_sigma, bins=nbins,range=[-0.4,0.4])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_s=xbins2[0:nbins] + 0.5*width2
yd_s=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='EB space')

nx2, xbins2, ptchs = plt.hist(sn_qu[2][ind]/cmb_sigma, bins=nbins,range=[-0.4,0.4])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_qu=xbins2[0:nbins] + 0.5*width2
yd_qu=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='QU space')

nx2, xbins2, ptchs = plt.hist(sn_P[2][ind]/cmb_sigma, bins=nbins,range=[-0.4,0.4])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_p=xbins2[0:nbins] + 0.5*width2
yd_p=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='P space')

nx2, xbins2, ptchs = plt.hist(sn_pc[2][ind]/cmb_sigma, bins=nbins,range=[-0.7,0.7])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_pc=xbins2[0:nbins] + 0.5*width2
yd_pc=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='pCILC')

nx2, xbins2, ptchs = plt.hist(sn_c[2][ind]/cmb_sigma, bins=nbins,range=[-0.7,0.7])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_c=xbins2[0:nbins] + 0.5*width2
yd_c=nx2 /np.sum(nx2)


plt.figure()
plt.title('B noise residual hist')
plt.plot(xd_s,yd_s,linewidth=2,label='EB space')
plt.plot(xd_qu,yd_qu,linewidth=2,label='QU space')
plt.plot(xd_p,yd_p,linewidth=2,label='P space')
plt.plot(xd_c,yd_c,linewidth=2,label='cILC')
plt.plot(xd_pc,yd_pc,linewidth=2,label='pcILC')
plt.legend()
plt.savefig("pdfs/wmap_planck_pysm/dust/%s/hist_noise_Bmode_nc%s_ns%d.pdf"%(model,ncrr,nside))
#plt.show()
np.savetxt("maps/wmap_planck_pysm/dust/%s/hist_noise_Bmode_nc%s_ns%d.txt"%(model,ncrr,nside),np.array([xd_s,yd_s,xd_qu,yd_qu,xd_p,yd_p,xd_c,yd_c,xd_pc,yd_pc]))




#### E-mode ####
cmb_sigma=np.std(dust[1][ind])
nbins=100
plt.figure()
plt.title('E hist')
nx2, xbins2, ptchs = plt.hist(s[1][ind], bins=nbins, range=[-40,120])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_s=xbins2[0:nbins] + 0.5*width2
yd_s=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='EB space')

nx2, xbins2, ptchs = plt.hist(s_qu[1][ind], bins=nbins, range=[-40,120])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_qu=xbins2[0:nbins] + 0.5*width2
yd_qu=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='QU space')

nx2, xbins2, ptchs = plt.hist(s_P[1][ind], bins=nbins, range=[-40,120])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_p=xbins2[0:nbins] + 0.5*width2
yd_p=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='P space')

nx2, xbins2, ptchs = plt.hist(dust[1][ind], bins=nbins, range=[-40,120])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_in=xbins2[0:nbins] + 0.5*width2
yd_in=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='input')

nx2, xbins2, ptchs = plt.hist(s_c[1][ind], bins=nbins, range=[-40,120])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_c=xbins2[0:nbins] + 0.5*width2
yd_c=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='CILC')


nx2, xbins2, ptchs = plt.hist(s_pc[1][ind], bins=nbins, range=[-40,120])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_pc=xbins2[0:nbins] + 0.5*width2
yd_pc=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='pCILC')



plt.figure()
plt.title('E hist')
plt.plot(xd_s,yd_s,linewidth=2,label='EB space')
plt.plot(xd_qu,yd_qu,linewidth=2,label='QU space')
plt.plot(xd_p,yd_p,linewidth=2,label='P space')
plt.plot(xd_c,yd_c,linewidth=2,label='cILC')
plt.plot(xd_pc,yd_pc,linewidth=2,label='pcILC')
plt.plot(xd_in,yd_in,linewidth=2,label='input')
plt.legend()
plt.savefig("pdfs/wmap_planck_pysm/dust/%s/hist_dust_Emode_nc%s_ns%d.pdf"%(model,ncrr,nside))
#plt.show()
np.savetxt("maps/wmap_planck_pysm/dust/%s/hist_dust_Emode_nc%s_ns%d.txt"%(model,ncrr,nside),np.array([xd_s,yd_s,xd_qu,yd_qu,xd_p,yd_p,xd_c,yd_c,xd_pc,yd_pc,xd_in,yd_in]))


### B-mode residuals ###
##fg
cmb_sigma=1.0
nbins=100
plt.figure()
plt.title('E fg residual hist')
nx2, xbins2, ptchs = plt.hist(sf[1][ind]/cmb_sigma, bins=nbins, range=[-4,4])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_s=xbins2[0:nbins] + 0.5*width2
yd_s=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='EB space')

nx2, xbins2, ptchs = plt.hist(sf_qu[1][ind]/cmb_sigma, bins=nbins,range=[-4,4])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_qu=xbins2[0:nbins] + 0.5*width2
yd_qu=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='QU space')

nx2, xbins2, ptchs = plt.hist(sf_P[1][ind]/cmb_sigma, bins=nbins, range=[-4,4])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_p=xbins2[0:nbins] + 0.5*width2
yd_p=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='P space')

nx2, xbins2, ptchs = plt.hist(sf_pc[1][ind]/cmb_sigma, bins=nbins,range=[-0.01,0.01])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_pc=xbins2[0:nbins] + 0.5*width2
yd_pc=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='pCILC')

nx2, xbins2, ptchs = plt.hist(sf_c[1][ind]/cmb_sigma, bins=nbins,range=[-0.01,0.01])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_c=xbins2[0:nbins] + 0.5*width2
yd_c=nx2 /np.sum(nx2)


plt.figure()
plt.title('E fg residual hist')
plt.plot(xd_s,yd_s,linewidth=2,label='EB space')
plt.plot(xd_qu,yd_qu,linewidth=2,label='QU space')
plt.plot(xd_p,yd_p,linewidth=2,label='P space')
plt.plot(xd_c,yd_c,linewidth=2,label='cILC')
plt.plot(xd_pc,yd_pc,linewidth=2,label='pcILC')
plt.legend()
plt.savefig("pdfs/wmap_planck_pysm/dust/%s/hist_fg_Emode_nc%s_ns%d.pdf"%(model,ncrr,nside))
#plt.show()
np.savetxt("maps/wmap_planck_pysm/dust/%s/hist_fg_Emode_nc%s_ns%d.txt"%(model,ncrr,nside),np.array([xd_s,yd_s,xd_qu,yd_qu,xd_p,yd_p,xd_c,yd_c,xd_pc,yd_pc]))

###noise residual####
nbins=100
nx2, xbins2, ptchs = plt.hist(sn[1][ind]/cmb_sigma, bins=nbins,range=[-0.3,0.3])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_s=xbins2[0:nbins] + 0.5*width2
yd_s=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='EB space')

nx2, xbins2, ptchs = plt.hist(sn_qu[1][ind]/cmb_sigma, bins=nbins,range=[-0.3,0.3])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_qu=xbins2[0:nbins] + 0.5*width2
yd_qu=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='QU space')

nx2, xbins2, ptchs = plt.hist(sn_P[1][ind]/cmb_sigma, bins=nbins,range=[-0.3,0.3])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_p=xbins2[0:nbins] + 0.5*width2
yd_p=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='P space')

nx2, xbins2, ptchs = plt.hist(sn_pc[1][ind]/cmb_sigma, bins=nbins,range=[-0.7,0.7])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_pc=xbins2[0:nbins] + 0.5*width2
yd_pc=nx2 /np.sum(nx2)
#plt.plot(xd_in,yd_in,linewidth=2,label='pCILC')

nx2, xbins2, ptchs = plt.hist(sn_c[1][ind]/cmb_sigma, bins=nbins,range=[-0.7,0.7])
plt.clf()
width2 = xbins2[1] - xbins2[0]
xd_c=xbins2[0:nbins] + 0.5*width2
yd_c=nx2 /np.sum(nx2)


plt.figure()
plt.title('E noise residual hist')
plt.plot(xd_s,yd_s,linewidth=2,label='EB space')
plt.plot(xd_qu,yd_qu,linewidth=2,label='QU space')
plt.plot(xd_p,yd_p,linewidth=2,label='P space')
plt.plot(xd_c,yd_c,linewidth=2,label='cILC')
plt.plot(xd_pc,yd_pc,linewidth=2,label='pcILC')
plt.legend()
plt.savefig("pdfs/wmap_planck_pysm/dust/%s/hist_noise_Emode_nc%s_ns%d.pdf"%(model,ncrr,nside))
#plt.show()
np.savetxt("maps/wmap_planck_pysm/dust/%s/hist_noise_Emode_nc%s_ns%d.txt"%(model,ncrr,nside),np.array([xd_s,yd_s,xd_qu,yd_qu,xd_p,yd_p,xd_c,yd_c,xd_pc,yd_pc]))




