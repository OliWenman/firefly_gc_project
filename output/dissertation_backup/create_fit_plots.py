import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import sys
import os
from astropy.io import fits
from math import floor, ceil

paths = [
	"output/dissertation/MASTAR_TH_VMPL7_KR",
	"output/dissertation/MASTAR_TH_VMPL9_KR",
	"output/dissertation/MASTAR_TH_VMPL11_KR",
	#"output/dissertation/MASTAR_TH_V0.2_KR/downgraded",
	"output/dissertation/MASTAR_E_VMPL7_KR",
	#"output/dissertation/MASTAR_E_V0.2_KR/downgraded",
	"output/dissertation/CONROY_E_KR/downgraded",
]

#SETTINGS
factor = 3
fontsize   = 6*factor
loc        = 'upper right' 
framealpha = 0.5

linewidth   = i_linewidth = linewidth_object = 4
d_linewidth = 0.7

alpha = 1

plt.rcParams.update({'font.size': fontsize})
mpl.rcParams['axes.linewidth'] = 2

rows    = 8
columns = 11

#rows    = 2
#columns = 4

#################################################
files_front = 0
files_end   = int((rows * columns) / 2)

continue_loop = True

file_index = 0

while continue_loop:

	fig = plt.figure(figsize=(11.69*factor, 8.27*factor))	
	file_index = file_index + 1 

	linewidth = i_linewidth

	for path in paths:
		
		linewidth = linewidth - d_linewidth


	 	#Get every file in the path of the folder and store in a list
		files = [os.path.join(path, f) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
		"""
		files = [
			os.path.join(path, "spFly-NGC0330_2015-11-02.fits"), 
			os.path.join(path, "spFly-NGC1846_2016-04-02.fits"),
			os.path.join(path, "spFly-NGC6528_2015-04-13.fits"),
			os.path.join(path, "spFly-NGC7099_2016-10-01.fits"),
		]
		"""

		if files_end > len(files):
			files_end = len(files) 
			continue_loop = False

		files = files[files_front:files_end]

		#index = 1

		#fig = plt.figure(figsize=(11.69, 8.27))	
		#spec = gridspec.GridSpec(ncols=columns, nrows=rows, figure=fig)

		"""
		for row in range(1, rows + 1):

			for column in range(1, columns + 1):

				print(row * column)

				ax = fig.add_subplot(rows, columns, (row*column))
				ax.set_yticklabels([])
				ax.set_xticklabels([])
				ax.tick_params(axis = "x", direction = "in")
				ax.tick_params(axis = "y", direction = "in")
		"""

		for file in files:
			print(file)

		for index in range(0, rows * columns):

			ax = fig.add_subplot(rows, columns, index + 1)
			#ax.set_yticklabels([])
			
			ax.tick_params(axis = "x", direction = "in")
			ax.tick_params(axis = "y", direction = "in")
			ax.set_xlim(0.36, 0.9)
			ax.grid(True)

			x = floor(index % columns)
			y = floor(index / columns)

			#print("Index ", index, x, y)

			#if y == rows - 1:
			#	y = y- 1
			#	current_row = current_row - 2
			#	print("here")


			if index > rows*(columns -columns -1): 
				ax.set_xlabel("Wavelength(\u03bcm)")
			else:
				ax.set_xticklabels([])


			if (y % 2) == 0:

				ax.set_ylim(-0.5, 3)

				if x == 0:
					ax.set_ylabel("Flux")
				else:
					ax.set_yticklabels([])

				y = int(y/2)

				new_index = x + (y * columns)

				#if index == new_index:
				#	print("the same")

				try:
					with fits.open(files[new_index], memmap=True) as hdul:

						wave       = hdul[1].data['wavelength'] / 10000
						flux       = hdul[1].data['original_data']
						model_flux = hdul[1].data['firefly_model']

						head, tail = os.path.split(files[new_index])
						name = tail[6:tail.find("_")]

						model = hdul[1].header['model']

						if "MPL7" in path:
							version = "(MPL7)"
						elif "VMPL9" in path:
							version = "(VMPL9)"
						else:
							version = ""

						model = model + version

						#Check the type of model used and assign a colour
						if "E-CONROY" in model.upper():
							color = "red"
						elif "TH-MASTAR" in model.upper():
							color = "royalblue"
							if "VMPL9" in model.upper():
								color = "navy"
							elif "VMPL11" in path.upper():
								color = "blueviolet"

						elif "E-MASTAR" in model.upper():
							color = "lime"	

						if path == paths[0]:
							ax.plot(wave, flux, label = name, linewidth = linewidth_object, alpha = alpha, c = "m")

						if index == 0:
							ax.plot(wave, model_flux, label = model, linewidth = linewidth, alpha = alpha, color = color)
						else:
							ax.plot(wave, model_flux, linewidth = linewidth, alpha = alpha, color = color)

						ax.legend(fontsize=fontsize, loc=loc, framealpha = framealpha)
						print("Index ", index, x, y)
						#print(new_index, "P")
				except(IndexError):
					print(new_index, x, y, "error - P")
					#pass

			else:

				ax.set_ylim(-0.75, 0.75)

				#if path == paths[0]:

				if x == 0:
					ax.set_ylabel("Residual Flux")
				else:
					ax.set_yticklabels([])

				y = int(ceil(y/2) -1)
				new_index = (x)+ (y * columns)
				
				try:
					with fits.open(files[new_index], memmap=True) as hdul:

						wave       = hdul[1].data['wavelength'] / 10000
						flux       = hdul[1].data['original_data']
						model_flux = hdul[1].data['firefly_model']

						head, tail = os.path.split(files[new_index])
						name = tail[6:tail.find("_")]

						model = hdul[1].header['model']

						#Check the type of model used and assign a colour
						if "E-CONROY" in model.upper():
							color = "red"
						elif "TH-MASTAR" in model.upper():
							color = "royalblue"
							if "VMPL9" in path.upper():
								color = "navy"
							elif "VMPL11" in path.upper():
								color = "blueviolet"

						elif "E-MASTAR" in model.upper():
							color = "lime"	

						ax.plot(wave, flux - model_flux, label = model, linewidth = linewidth, alpha = alpha, color = color)
						#ax.legend(fontsize=fontsize, loc=loc, framealpha = framealpha)

						print("Index ", index, x, y)
					
				except(IndexError):
					print(new_index, "error - R")
					#pass
				

			
			#if (y % 2) == 0:

				
	fig.tight_layout()
	fig.subplots_adjust(wspace=0, hspace=0)
	fig.savefig("output/dissertation/data/fits/" + "fit_plots" + str(file_index) + ".png")

	try:
		display_plot = sys.argv[1].lower() == 'true'
	except:
		display_plot = False

	if display_plot:
		plt.show()

	files_front = files_end

	files_end = files_end * 2





