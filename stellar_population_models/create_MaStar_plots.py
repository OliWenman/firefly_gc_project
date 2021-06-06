import MaStar_SSP 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from astropy.io import fits
import math

plt.rcParams.update({'font.size': 26})

def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

ver    = "vMPL7"
lib_in = ["E", "Th"]
colour = []
tin    = [13]
Zin    = [-1.5, -1, -0.5]
sin    = 2.35 #2.35 = Salpeter, 1.3 = Kroupa

fig   = plt.figure(figsize=(11.69*2, 8.27*1))	

columns = 1
rows    = 1

ax1 = fig.add_subplot(1, 1, 1)

darken_amount = 0.5

print(tin)

for time in tin:
	for metal in Zin:
		for lib in lib_in:

			alpha = 1

			if lib == "E":
				colour = "lime"
			elif lib == "Th":
				colour = "royalblue"

			t        = MaStar_SSP.t(ver) #age array
			Z        = MaStar_SSP.Z(ver) #metallicty array
			s        = MaStar_SSP.s(ver) #IMF slope array
			wave     = MaStar_SSP.wave(ver) #wavelength array
			#res      = MaStar_SSP.res(ver) #resolution array
			fluxgrid = MaStar_SSP.flux(ver,lib) #flux matrix

			flux = MaStar_SSP.inter(t,Z,s,fluxgrid,time,metal,sin)

			wave = wave/10000

			ax1.plot(wave, flux, color = lighten_color(colour, darken_amount), alpha =alpha, label = lib + "-MaStar: [Z/H] = " + str(metal), lw = 3.5)

		darken_amount = darken_amount + 0.5


ax1.legend(loc = "upper right")
ax1.set_xlabel(r'Wavelength ($\mu$m)')
ax1.set_ylabel("Flux")
#ax1.set_title("Stellar Population Model MaStar with constant IMF = " + str(sin) + ", age = " + str(tin[0]) + "Gyr but as a function of [Z/H]")
plt.tight_layout()
plt.show()

fig.savefig("./plots/MaStar_metalicity.png")