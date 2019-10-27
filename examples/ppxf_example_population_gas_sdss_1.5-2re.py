#!/usr/bin/env python
##############################################################################
#
# Usage example for the procedure PPXF, which implements the
# Penalized Pixel-Fitting (pPXF) method originally described in
# Cappellari M., & Emsellem E., 2004, PASP, 116, 138
#     http://adsabs.harvard.edu/abs/2004PASP..116..138C
# and upgraded in Cappellari M., 2017, MNRAS, 466, 798
#     http://adsabs.harvard.edu/abs/2017MNRAS.466..798C
#
# This example shows how to study stellar population and include gas emission
# lines as templates instead of masking them using the GOODPIXELS keyword.
#

from time import perf_counter as clock
from os import path

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from astropy.io import fits

import ppxf as ppxf_package
from ppxf.ppxf import ppxf
import ppxf.ppxf_util as util
import ppxf.miles_util as lib
import glob
import os

##############################################################################

def ppxf_example_population_gas_sdss(tie_balmer, limit_doublets):
    
    

    ppxf_dir = path.dirname(path.realpath(ppxf_package.__file__))

    # Read SDSS DR8 galaxy spectrum taken from here http://www.sdss3.org/dr8/
    # The spectrum is *already* log rebinned by the SDSS DR8
    # pipeline and log_rebin should not be used in this case.
    #
    
    f_filename=fits.open('/home/zhaisai/zhaisai/data/z_float.fits')
    data_filename=f_filename[1].data
    z_filename=data_filename.field('col1')
    filename=data_filename.field('col0')
    f_sf=fits.open('/home/zhaisai/zhaisai/data/SF.fits')
    sf_name=f_sf[1].data.field('plateifu')
    for i in range(1,len(filename),1):
        dirname='/home/zhaisai/zhaisai/data/stack_file/'+str(filename[i])+'.fits'
        if os.path.exists(dirname) and str(filename[i]) in sf_name:
            print(dirname)
            f_sig=fits.open(dirname)
            t=f_sig[1].data
            if tie_balmer==False and limit_doublets==False:
                dirpic='/home/zhaisai/zhaisai/result/1.5-2re/'+str(filename[i])+'.eps'
            elif tie_balmer==True and limit_doublets==True:
                dirpic='/home/zhaisai/zhaisai/result/1.5-2re/'+str(filename[i])+'-Balmer.eps'
    
            dirfile='/home/zhaisai/zhaisai/result/1.5-2re/'+str(filename[i])+'.txt'
            
            if os.path.exists(dirfile):
                continue
            else:
                z =z_filename[i]# SDSS redshift estimate
            # Only use the wavelength range in common between galaxy and stellar library.
            #
                mask1 = (t['col0'] < 3540) | (t['col0'] > 7409)
                j=0
                for i in t.columns:
                    j+=1
                if j>=5:
                    flux=np.ma.array(t['col4'],mask=mask1)
                      # Normalize spectrum to avoid numerical issues
                    wave = np.ma.array(t['col0'],mask=mask1)
                    
                    flux=list(filter(None,flux))
                    flux=np.array(flux)
            #        print(flux)
                    wave=list(filter(None,wave))
                    wave=np.array(wave)
            #        print(flux.shape,wave.shape)
                    galaxy = flux/np.median(flux) 
            #        print(galaxy.shape)
    #                plt.plot(wave,flux)
                # The SDSS wavelengths are in vacuum, while the MILES ones are in air.
                # For a rigorous treatment, the SDSS vacuum wavelengths should be
                # converted into air wavelengths and the spectra should be resampled.
                # To avoid resampling, given that the wavelength dependence of the
                # correction is very weak, I approximate it with a constant factor.
                #
                    wave *= np.median(util.vac_to_air(wave)/wave)
            
            
                # The noise level is chosen to give Chi^2/DOF=1 without regularization (REGUL=0).
                # A constant noise is not a bad approximation in the fitted wavelength
                # range and reduces the noise in the fit.
                #
                    noise = np.full_like(galaxy, 0.01635)  # Assume constant noise per pixel here
            
                # The velocity step was already chosen by the SDSS pipeline
                # and we convert it below to km/s
                #
                    c = 299792.458  # speed of light in km/s
                    velscale = c*np.log(wave[1]/wave[0])  # eq.(8) of Cappellari (2017)
            #        print(c*np.log(3541.4/3540.5))
                    FWHM_gal = 2.76  # SDSS has an approximate instrumental resolution FWHM of 2.76A.
            
                #------------------- Setup templates -----------------------
            
                    pathname = ppxf_dir + '/miles_models/Mun1.30*.fits'
            
                # The templates are normalized to mean=1 within the FWHM of the V-band.
                # In this way the weights and mean values are light-weighted quantities
                    miles = lib.miles(pathname=pathname,velscale=velscale,FWHM_gal=FWHM_gal)
            
                # The stellar templates are reshaped below into a 2-dim array with each
                # spectrum as a column, however we save the original array dimensions,
                # which are needed to specify the regularization dimensions
                #
                    reg_dim = miles.templates.shape[1:]
                    stars_templates = miles.templates.reshape(miles.templates.shape[0], -1)
            
                # See the pPXF documentation for the keyword REGUL,
                    regul_err = 0.013  # Desired regularization error
            
                # Construct a set of Gaussian emission line templates.
                # Estimate the wavelength fitted range in the rest frame.
                #
                    lam_range_gal = np.array([np.min(wave), np.max(wave)])/(1 + z)
                    gas_templates, gas_names, line_wave = util.emission_lines(
                    miles.log_lam_temp, lam_range_gal, FWHM_gal,
                    tie_balmer=tie_balmer, limit_doublets=limit_doublets)
            
                # Combines the stellar and gaseous templates into a single array.
                # During the PPXF fit they will be assigned a different kinematic
                # COMPONENT value
                #
                    templates = np.column_stack([stars_templates, gas_templates])
            #        print(templates.shape)
            
                #-----------------------------------------------------------
            
                # The galaxy and the template spectra do not have the same starting wavelength.
                # For this reason an extra velocity shift DV has to be applied to the template
                # to fit the galaxy spectrum. We remove this artificial shift by using the
                # keyword VSYST in the call to PPXF below, so that all velocities are
                # measured with respect to DV. This assume the redshift is negligible.
                # In the case of a high-redshift galaxy one should de-redshift its
                # wavelength to the rest frame before using the line below as described
                # in PPXF_EXAMPLE_KINEMATICS_SAURON and Sec.2.4 of Cappellari (2017)
                #
                    c = 299792.458
                    dv = c*(miles.log_lam_temp[0] - np.log(wave[0]))  # eq.(8) of Cappellari (2017)
                    vel = c*np.log(1 + z)   # eq.(8) of Cappellari (2017)
                    start = [vel, 180.]  # (km/s), starting guess for [V, sigma]
            
                    n_temps = stars_templates.shape[1]
                    n_forbidden = np.sum(["[" in a for a in gas_names])  # forbidden lines contain "[*]"
                    n_balmer = len(gas_names) - n_forbidden
            
                # Assign component=0 to the stellar templates, component=1 to the Balmer
                # gas emission lines templates and component=2 to the forbidden lines.
                    component = [0]*n_temps + [1]*n_balmer + [2]*n_forbidden
                    gas_component = np.array(component) > 0  # gas_component=True for gas templates
            
                # Fit (V, sig, h3, h4) moments=4 for the stars
                # and (V, sig) moments=2 for the two gas kinematic components
                    moments = [4, 2, 2]
            
                # Adopt the same starting value for the stars and the two gas components
                    start = [start, start, start]
            
                # If the Balmer lines are tied one should allow for gas reddeining.
                # The gas_reddening can be different from the stellar one, if both are fitted.
                    gas_reddening = 0 if tie_balmer else None
            #        print(int(round(velscale/(c*np.log(3541.4/3540.5)))))
                # Here the actual fit starts.
                #
                # IMPORTANT: Ideally one would like not to use any polynomial in the fit
                # as the continuum shape contains important information on the population.
                # Unfortunately this is often not feasible, due to small calibration
                # uncertainties in the spectral shape. To avoid affecting the line strength of
                # the spectral features, we exclude additive polynomials (DEGREE=-1) and only use
                # multiplicative ones (MDEGREE=10). This is only recommended for population, not
                # for kinematic extraction, where additive polynomials are always recommended.
                #
                    t = clock()
                
                    pp = ppxf(dirfile,dirpic,templates, galaxy, noise,velscale, start,
                          plot=False, moments=moments, degree=-1, mdegree=10, vsyst=dv,
                          lam=wave, clean=False, regul=1./regul_err, reg_dim=reg_dim,
                          component=component, gas_component=gas_component,
                          gas_names=gas_names, gas_reddening=gas_reddening,
                          velscale_ratio=int(round(velscale/(c*np.log(3541.4/3540.5))))) 
            
                # When the two Delta Chi^2 below are the same, the solution
                # is the smoothest consistent with the observed spectrum.
                #print(pp.lam,pp.galaxy,pp.bestfit)
                #
            #    print('Desired Delta Chi^2: %.4g' % np.sqrt(2*galaxy.size))
            #    print('Current Delta Chi^2: %.4g' % ((pp.chi2 - 1)*galaxy.size))
            #    print('Elapsed time in PPXF: %.2f s' % (clock() - t))
            
                    weights = pp.weights[~gas_component]  # Exclude weights of the gas templates
                    weights = weights.reshape(reg_dim)/weights.sum()  # Normalized
            
                    miles.mean_age_metal(weights)
                    miles.mass_to_light(weights, band="r")
            
                # Plot fit results for stars and gas.
            #           plt.subplots_adjust(wspace =0, hspace =0.8)
                    plt.clf()
        #        plt.subplot(211)
                    pp.plot()
                    pp.gas_flux
                    pp.bestfit
                    pp.reddening
                    pp.reddening_func
 

    # Plot stellar population mass fraction distribution
#        plt.subplot(212)
#        miles.plot(weights,dirpic)
#        plt.tight_layout()
#        plt.pause(1)
#        plt.show()
##############################################################################

if __name__ == '__main__':

    print("\n===============================================\n" +
             " Fit with free Balmer lines and [SII] doublet: \n" +
             "===============================================\n")

    ppxf_example_population_gas_sdss(tie_balmer=False, limit_doublets=False)

#    print("\n=======================================================\n" +
#             " Fit with tied Balmer lines and limited [SII] doublet: \n" +
#             "=======================================================\n")
#
#    # Note tha the inclusion of a few extra faint Balmer lines is sufficient to
#    # decrease the chi2 of the fit, even though the Balmer decrement is fixed.
#    # In this case, the best-fitting gas reddening is at the E(B-V)=0 boundary.
#    ppxf_example_population_gas_sdss(tie_balmer=True, limit_doublets=True)
