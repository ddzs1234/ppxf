from __future__ import division
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy import stats
from astropy.table import Table

np.seterr(divide='ignore', invalid='ignore')

def kk04(logn2o,log23,logo32):
    metal=[]
    for i in range(0,len(logn2o)):
        log_n2o=logn2o[i]
        log_23=log23[i]
        log_o32=logo32[i]
        if np.isinf(log_n2o):
            metal.append(np.nan)
        elif log_n2o>=-1.2:
            q=((32.81-1.153*log_o32**2+8.7*(-3.396-0.025*log_o32+0.1444*log_o32**2))/
               (4.603-0.3119*log_o32-0.163*log_o32**2+8.7*(-0.48-0.271*log_o32+0.02037*log_o32**2)))
            metal.append((9.72-0.777*log_23-0.951*log_23**2-0.072*log_23**3-0.811*log_23**4-
                np.log10(q)*(0.0737-0.0713*log_23-0.141*log_23**2+0.0373*log_23**3-0.058*log_23**4)))
        else:
            q=((32.81-1.153*log_o32**2+8.2*(-3.396-0.025*log_o32+0.1444*log_o32**2))/
               (4.603-0.3119*log_o32-0.163*log_o32**2+8.2*(-0.48-0.271*log_o32+0.02037*log_o32**2)))
            metal.append(9.4+4.65*log_23-3.17*log_23**2-np.log10(q)*(0.272+0.547*log_23-0.513*log_23**2))
    return np.array(metal)


def mask_m(a,b):
    return np.ma.array(a,mask=(b))

f=fits.open('/home/ashley/NewDisk/pair_galaxy/m-z-manga/result/stack/manga_mass>10/info_1.fits')
data=f[1].data
fmass=fits.open('/home/ashley/NewDisk/pair_galaxy/m-z-manga/data/mpl6_mass_sfr_morph.fits')
data_mass=fmass[1].data
plateifu_mass=data_mass.field('plateifu')
mass_s=data_mass.field('mass')
mass_l=data_mass.field('LMASS50_ALL')

oii_3726=data.field('OII_3726')
oii_3729=data.field('OII_3729')
oiii=data.field('OIII')
hbeta=data.field('Hbeta')
nii=data.field('NII')
halpha=data.field('Halpha')
plateifu=data.field('PLATEIFU')

mass=[]
for i in plateifu:
    index=list(plateifu_mass).index(i)
    if mass_s[index]>0:
        mass.append(mass_s[index])
    else:
        mass.append(mass_l[index])
mass=np.array(mass)
#########
logn2o=np.log10((nii)/(oii_3726+oii_3729))
logo32=np.log10(oiii/(oii_3726+oii_3729))
log23=np.log10((oii_3726+oii_3729+oiii*4/3)/hbeta)

kk04_1=np.array(kk04(logn2o,log23,logo32),dtype=float)
c=np.isnan(kk04_1)|(kk04_1<-10)|np.isinf(kk04_1)|np.isnan(mass)|(mass<0)

kk04_1=mask_m(kk04_1,c)
mass=mask_m(mass,c)
plateifu=mask_m(plateifu,c)

kk04_1=filter(None,kk04_1)
mass=filter(None,mass)
plateifu=filter(None,plateifu)

#fileresult=open('/home/ashley/NewDisk/pair_galaxy/m-z-manga/result/stack/KK04_mass10.txt','a+')

for i in range(0,len(kk04_1)):
    if kk04_1[i]>8.5 and kk04_1[i]<9.0:
        print(kk04_1[i],mass[i],plateifu[i])

#t=Table([plateifu,kk04_1,mass],names=['PLATEIFU','KK04_Z','MASS'],dtype=['str','f8','f8'])
#
#t.write('/home/ashley/NewDisk/pair_galaxy/m-z-manga/result/stack/KK04_plateifu_stack.fits',format='fits')

bin_means,bin_edges,binnumber=stats.binned_statistic(np.array(mass).astype('float'),np.array(kk04_1).astype('float'),'mean',bins=10)



plt.figure()
plt.scatter(mass,kk04_1,s=2,c='chocolate')
plt.hlines(bin_means,bin_edges[:-1],bin_edges[1:],color='b',lw=1.2,label='binned statistic of data')
plt.ylabel('$12+log(O/H) (KK04)$')
plt.xlabel('$log(M/M_{odot})$')
plt.ylim(6.0,10.0)
plt.legend()
#plt.savefig('/home/ashley/NewDisk/pair_galaxy/m-z-manga/result/ext_cor_2/m_z_KK04_v3.png',format='png')