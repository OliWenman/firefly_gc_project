import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

hdul = fits.open('MaStar_SSP_v0.2.fits.gz')

hdul.info()

#print(hdul[1].data)

print(hdul[3].header)
fluxgrid=hdul[3].data #fluxgrid matrix
flux=np.ndarray(len(fluxgrid[0,0,0,:]))

plt.plot(flux)
print(fluxgrid)
plt.show()

hdul.close()