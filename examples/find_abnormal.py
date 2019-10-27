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




f=fits.open('/home/ashley/NewDisk/pair_galaxy/m-z-manga/result/ext_cor_2/KD02_plateifu.fits')
data=f[1].data
z_do2=data.field('KD02_Z')
mass=data.field('MASS')
plateifu=data.field('PLATEIFU')
f_do2=open('/home/ashley/NewDisk/pair_galaxy/m-z-manga/result/ext_cor_2/KD02_Z8.5.txt','a+' )
print('plateifu','mass','z_KD02','\t',file=f_do2)
for i in range(0,len(z_do2)):
    if z_do2[i]<8.5:
        print(plateifu[i],mass[i],z_do2[i],'\t',file=f_do2)
        
f_do2.close()






#plt.figure()
#plt.scatter(mass,do2_1,s=2,c='dodgerblue')
#plt.hlines(bin_means,bin_edges[:-1],bin_edges[1:],color='b',lw=1.2,label='binned statistic of data')
#plt.ylabel('$12+log(O/H) (D02)$')
#plt.xlabel('$log(M/M_{odot})$')
#plt.ylim(7.8,np.max(do2_1)*1.01)
#plt.legend()
#plt.savefig('/home/ashley/NewDisk/pair_galaxy/m-z-manga/result/ext_cor_2/m_z_D02_v3.png',format='png')