#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 15:53:47 2018

@author: ashley
"""

from __future__ import division
import sympy as sy
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy import stats
from astropy.table import Table

np.seterr(divide='ignore', invalid='ignore')

def do2(xi):
    metal=[]
    for i in range(0,len(xi),1):
        x=xi[i]
        metal.append(9.12+0.73*x)
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

logn2=np.log10(nii/halpha)

do2_1=np.array(do2(logn2),dtype=float)
c=np.isnan(mass)|(mass<0)|np.isnan(do2_1)|np.isinf(do2_1)|(do2<-10)

do2_1=mask_m(do2_1,c)
mass=mask_m(mass,c)
plateifu=mask_m(plateifu,c)

do2_1=filter(None,do2_1)
mass=filter(None,mass)
plateifu=filter(None,plateifu)

#t=Table([plateifu,do2_1,mass],names=['PLATEIFU','DO2_1','MASS'],dtype=['str','f8','f8'])
#t.write('/home/ashley/NewDisk/pair_galaxy/m-z-manga/result/ext_cor_2/DO2_plateifu.fits',format='fits')
bin_means,bin_edges,binnumber=stats.binned_statistic(np.array(mass).astype('float'),np.array(do2_1).astype('float'),'mean',bins=10)





plt.figure()
plt.scatter(mass,do2_1,s=2,c='dodgerblue')
plt.hlines(bin_means,bin_edges[:-1],bin_edges[1:],color='b',lw=1.2,label='binned statistic of data')
plt.ylabel('$12+log(O/H) (D02)$')
plt.xlabel('$log(M/M_{odot})$')
plt.ylim(7,np.max(do2_1)*1.01)
plt.legend()
#plt.savefig('/home/ashley/NewDisk/pair_galaxy/m-z-manga/result/ext_cor_2/m_z_D02_v3.png',format='png')

