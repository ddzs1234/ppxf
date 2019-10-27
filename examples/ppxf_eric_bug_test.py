#!/usr/bin/env python
#----------------------------------------------------------------------------
#
# Usage example for the procedure PPXF, which implements the
# Penalized Pixel-Fitting (pPXF) method originally described in
# Cappellari M., & Emsellem E., 2004, PASP, 116, 138
#     http://adsabs.harvard.edu/abs/2004PASP..116..138C
# and upgraded in Cappellari M., 2017, MNRAS, 466, 798
#     http://adsabs.harvard.edu/abs/2017MNRAS.466..798C
#
# This example it is useful to determine the desired value for
# the BIAS keyword of the pPXF procedure. This procedure generates
# a plot similar (but not identical) to Figure 6 in
# Cappellari & Emsellem (2004).
#
# A rough guideline to determine the BIAS value is the following: choose the *largest*
# value which make sure that in the range sigma>3*velScale and for (S/N)>30 the true values
# for the Gauss-Hermite parameters are well within the rms scatter of the measured values.
# See the documentation in the file ppxf.pro for a more accurate description.
#
# MODIFICATION HISTORY:
#   V1.0.0: By Michele Cappellari, Leiden, 28 March 2003
#   V1.1.0: Included in the standard PPXF distribution. After feedback
#       from Alejandro Garcia Bedregal. MC, Leiden, 13 April 2005
#   V1.1.1: Adjust GOODPIXELS according to the size of the convolution kernel.
#       MC, Oxford, 13 April 2010
#   V1.1.2: Use Coyote Graphics (http://www.idlcoyote.com/) by David W. Fanning.
#       The required routines are now included in NASA IDL Astronomy Library.
#       MC, Oxford, 29 July 2011
#   V2.0.0: Translated from IDL into Python. MC, Oxford, 9 December 2013
#   V2.0.1: Support both Python 2.6/2.7 and Python 3.x. MC, Oxford, 25 May 2014
#   V2.0.2: Support both Pyfits and Astropy to read FITS files.
#       MC, Oxford, 22 October 2015
#   V2.0.3: Use random input velocity to properly simulate situations with
#       undersampled LOSVD. MC, Oxford, 20 April 2016
#   V2.1.0: Replaced the Vazdekis-99 SSP models with the Vazdekis+10 ones.
#       Modified plot to emphasize undersampling effect. MC, Oxford, 3 May 2016
#   V2.1.1: Make files paths relative to this file, to run the example from
#       any directory. MC, Oxford, 18 January 2017
#   V2.1.2: Updated MILES file name. MC, Oxford, 29 November 2017
#   V2.1.3: Changed imports for pPXF as a package.
#       Make file paths relative to the pPXF package to be able to run the
#       example unchanged from any directory. MC, Oxford, 17 April 2018
#   V2.1.4: Dropped legacy Python 2.7 support. MC, Oxford, 10 May 2018
#   V4.0.3: Fixed clock DeprecationWarning in Python 3.7.
#       MC, Oxford, 27 September 2018
#
##############################################################################

from time import perf_counter as clock
from os import path

from astropy.io import fits
from scipy import ndimage, signal
import numpy as np
import matplotlib.pyplot as plt

import ppxf as ppxf_package
from ppxf.ppxf import ppxf, rebin
import ppxf.ppxf_util as util

#----------------------------------------------------------------------------

def ppxf_example_simulation():

    ppxf_dir = path.dirname(path.realpath(ppxf_package.__file__))

    hdu = fits.open(ppxf_dir + '/miles_models/Mun1.30Zp0.00T12.5893_iPp0.00_baseFe_linear_FWHM_2.51.fits')  # Solar metallicitly, Age=12.59 Gyr
    ssp = hdu[0].data
    h = hdu[0].header

    lamRange = h['CRVAL1'] + np.array([0.,h['CDELT1']*(h['NAXIS1']-1)])
    c = 299792.458 # speed of light in km/s
    velscale = c*h['CDELT1']/max(lamRange)   # Do not degrade original velocity sampling
    galaxy, logLam, velscale = util.log_rebin(lamRange, ssp, velscale=velscale)


    start = [1000., 100.]
    dx = 1

    pp = ppxf(galaxy, galaxy, galaxy*0+1, velscale, start, plot=1,
              goodpixels=np.arange(dx, galaxy.size - dx))

    print(pp.sol)

#----------------------------------------------------------------------------

if __name__ == '__main__':

    np.random.seed(123)  # For reprodcible results
    ppxf_example_simulation()
