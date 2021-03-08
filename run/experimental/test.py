import sys, os

sys.path.append(os.path.join(os.getcwd(), "python"))
os.environ["FF_DIR"] = os.getcwd()
os.environ["STELLARPOPMODELS_DIR"] = os.path.join(os.environ["FF_DIR"], "stellar_population_models")

import numpy as np
from astropy.io import fits
import astropy.cosmology as co
import firefly_setup as fs
import firefly_models as fm
import time

import matplotlib.pyplot as plt

hdul = fits.open(sys.argv[1])

cdelt1 = hdul[2].header['CDELT1']    #change in wavelength between pixels
crval1 = hdul[2].header['CRVAL1']   
crval2 = hdul[2].header['CRVAL2']       

wavelength = np.arange(crval1,crval2,cdelt1)  #create wavelength array
flux = hdul[2].data   
error = hdul[1].data

wavelength = wavelength.reshape((-1,))
flux       = flux.reshape((-1,))
error      = flux.reshape((-1,))

print(error.shape)
print(flux.shape)


redshift = hdul[0].header['redshift']
ra       = hdul[0].header['ra']
dec      = hdul[0].header['dec']
vdisp    = hdul[0].header['veldisp']


plt.plot(wavelength, flux)
plt.fill_between(wavelength, flux-error, flux+error, alpha = 0.5)
plt.show()

