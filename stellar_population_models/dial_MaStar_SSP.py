
# Python routine that plots MaStar SSP spectrum for given version number, library, age, metallicity, and IMF slope
# Call with: python dial_MaStar_SSP.py <version number> <library> <age> <metallicity> <IMF slope>
#
# Example MPL7: python dial_MaStar_SSP.py v0.2 th 10.5 0.1 1.55
#
###################################################################
# Example MPL9: python dial_MaStar_SSP.py v1.0 gold 10.5 0.1 2.35 #
###################################################################
#
# Daniel Thomas, 12/06/2020

import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import MaStar_SSP
import sys

ver=str(sys.argv[1])
lib_in=str(sys.argv[2])
tin=float(sys.argv[3])
Zin=float(sys.argv[4])
sin=float(sys.argv[5])

lib=lib_in
if (lib_in=='th' or lib_in=='Th' or lib_in=='TH' or lib_in=='theoretical'):
	lib='Th'
if (lib_in=='e' or lib_in=='E' or lib_in=='empirical'):
	lib='E'

t=MaStar_SSP.t(ver) #age array
print(t)
Z=MaStar_SSP.Z(ver) #metallicty array
s=MaStar_SSP.s(ver) #IMF slope array
wave=MaStar_SSP.wave(ver) #wavelength array
res=MaStar_SSP.res(ver) #resolution array
fluxgrid=MaStar_SSP.flux(ver,lib) #flux matrix

print()
print('**Plotting '+lib+'-MaStar SSP ('+ver+') spectrum')
print('**age='+sys.argv[3]+' Gyr, [Z/H]='+sys.argv[4]+' dex, IMF slope s='+sys.argv[5])
print()

flux=MaStar_SSP.inter(t,Z,s,fluxgrid,tin,Zin,sin)

plt.plot(wave, res)

#plt.plot(wave,flux)
#plt.xlabel('wavelength')
#plt.ylabel('flux')

plt.show()
