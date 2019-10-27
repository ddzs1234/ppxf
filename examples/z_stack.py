#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 31 10:11:39 2018

@author: ashley
"""

'''
1. read flux
2. correlate flux_redden
3. compare Z
'''


from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import os
import glob

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

def rred(wl,f_ha,f_hb):
    k_ha=2.659*(-1.857+1.040/0.6564)+4.88
    k_hb=2.659*(-2.156+1.509/0.4862-0.198/(0.4862**2)+0.011/(0.4862**3))+4.88
    e_bv=np.log10((f_ha/f_hb)/2.86)/(-0.4*(k_ha-k_hb))
    wl=wl/10000.0
    if 0.63<=wl<=2.2:
        k=2.659*(-1.857+1.040/wl)+4.88
    elif 0.12<=wl<=0.63:
        k=2.659*(-2.156+1.509/wl-0.198/(wl**2)+0.011/(wl**3))+4.88
    #print(k_ha,k_hb,e_bv,k,f_ha,f_hb)
    return 10**(0.4*k*e_bv)





f_info=open('/home/ashley/NewDisk/pair_galaxy/m-z-manga/result/stack/manga_mass>10/info.txt','a+')
print('PLATEIFU','NII','OIII','OII_3726','OII_3729','Hbeta','Halpha','NII_r','OIII_r','OII_3726_r','OII_3729_r','Hbeta_r','Halpha_r',file=f_info)
allfile=glob.glob('/home/ashley/NewDisk/pair_galaxy/m-z-manga/result/stack/manga_mass>10/0.5re/*.txt')
for i in range(0,len(allfile)):
    filename=allfile[i]
    if os.path.exists(filename):
        print('filename',filename)
        f=open(filename,'r')
        col0=filename[76:]
        line=f.readlines()
        info=[]
        for i in range(10,len(line)):
            info.append(line[i])
        
        flux1=[i.strip().split()[3] for i in info]
        flux1=map(float,flux1)
        flux1=list(flux1)
        name1=[i.strip().split()[2] for i in info]
        if 'Halpha' in name1 and 'Hbeta' in name1:
            if flux1[name1.index('Hbeta')]>0 and flux1[name1.index('[OII]3726')]>0 and flux1[name1.index('Halpha')]>0 and '[NII]6583_d' in name1 and '[OII]3726' in name1 and '[OIII]5007_d' in name1:
                nii=flux1[name1.index('[NII]6583_d')]
                oiii=flux1[name1.index('[OIII]5007_d')]
                oii_3726=flux1[name1.index('[OII]3726')]
                oii_3729=flux1[name1.index('[OII]3729')]
                hbeta=flux1[name1.index('Hbeta')]
                halpha=flux1[name1.index('Halpha')]
                if halpha/hbeta>2.86:
                    
                    nii_r=nii*rred(6583,halpha,hbeta)
                    oiii_r=oiii*rred(5007,halpha,hbeta)
                    oii_3726_r=oii_3726*rred(3726,halpha,hbeta)
                    oii_3729_r=oii_3729*rred(3729,halpha,hbeta)
                    hbeta_r=hbeta*rred(4861,halpha,hbeta)
                    halpha_r=halpha*rred(6563,halpha,hbeta)
                    
                    no=np.log10(nii/(oii_3726+oii_3729))
                    
                    r23=np.log10((oiii*4/3+oii_3726+oii_3729)/hbeta)
                    o32=np.log10(oiii/(oii_3726+oii_3729))
                    z=kk04([no],[r23],[o32])
                    print(col0,nii,oiii,oii_3726,oii_3729,hbeta,halpha,nii_r,oiii_r,oii_3726_r,oii_3729_r,hbeta_r,halpha_r,file=f_info)
                    
        else:
            print('loss Halpha or Hbeta flux')
f_info.close()
