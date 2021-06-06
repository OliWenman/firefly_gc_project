import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from astropy.io import fits
import math
from math import floor, ceil

plt.rcParams.update({'font.size': 16})

if __name__ == "__main__":

	#Get optinal argument to display plots.
	#If no argument given, default is True
	try:
		display_plot = sys.argv[1].lower() == 'true'
	except:
		display_plot = False

	#Paths (relative to firefly directory) to read in the data processed by firefly

	paths = [
		"output/dissertation/MASTAR_TH_VMPL7_KR",
		"output/dissertation/MASTAR_TH_VMPL9_KR",
		"output/dissertation/MASTAR_TH_VMPL11_KR",
		#"output/dissertation/MASTAR_TH_V0.2_KR/downgraded",
		"output/dissertation/MASTAR_E_VMPL7_KR",
		#"output/dissertation/MASTAR_E_V0.2_KR/downgraded",
		"output/dissertation/CONROY_E_KR/downgraded",
	]

	factor = 1

	fig   = plt.figure(figsize=(11.69*factor, 8.27*factor))	
	index = 0

	columns = len(paths)
	rows    = 2

	#Loop through the model data
	for path in paths:

		index = index + 1

		if path.split('/')[-1] == "downgraded" and not "CONROY" in path.split('/')[-2]:
			downgraded_to_conroy = "- downgraded"
		else:
			downgraded_to_conroy = ""
	 
	 	#Get every file in the path of the folder and store in a list
		files = [os.path.join(path, f) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

		age_lightW_array = []
		age_massW_array  = []

		metal_lightW_array = []
		metal_massW_array  = []

		#Check the type of model used and assign a colour
		if "CONROY_E" in path.upper():
			color = "red"
			model = "E-Conroy"

		elif "MASTAR_TH" in path.upper():
			color = "royalblue"
			model = "Th-MaStar"

			if "VMPL7" in path.upper():
				model = model + "(MPL7)"

			elif "VMPL9" in path.upper():
				model = model + "(MPL9)"
				color = "navy"

			elif "VMPL11" in path.upper():
				model = model + "(MPL11)"
				color = "blueviolet"

		elif "MASTAR_E" in path.upper():
			color = "lime"
			model = "E-MaStar"	

			if "VMPL7" in path.upper():
				model = model + "(MPL7)"

			elif "VMPL9" in path.upper():
				model = model + "(MPL9)"

			if "VMPL11" in path.upper():
				model = model + "(MPL11)"

		#Loop through files of processed firefly output 
		for file in files:

			path_, file_ = os.path.split(file)

			#Get just the name of the object observed
			object_name = file_[6:file_.find("_")]

			#Extract the model data
			hdul = fits.open(file)

			age_lightW     = float(hdul[1].header['age_lightW'])
			metal_lightW   = float(hdul[1].header['metallicity_lightW'])
			age_massW      = float(hdul[1].header['age_massW'])
			metal_massW    = float(hdul[1].header['metallicity_massW'])
			
			hdul.close()

			age_lightW_array.append(age_lightW)
			age_massW_array.append(age_massW)

			metal_lightW_array.append(metal_lightW)
			metal_massW_array.append(metal_massW)

		ax1 = fig.add_subplot(rows, columns, index)
		ax1.scatter(age_massW_array, age_lightW_array, label = model, color = color, facecolors='none')
		ax1.plot([0, 1], [0, 1], 'g--', transform=ax1.transAxes, color = color)
		ax1.set_xlim(-2, 1.5)
		ax1.set_ylim(-2, 1.5)
		ax1.set_aspect('equal')
		ax1.set_xlabel("Log Age (Gyr) - mass weighted")
		ax1.set_xticks(np.arange(-2, 2, 1))
		ax1.legend()

		if index == 1:
			ax1.set_ylabel("Log Age (Gyr) - light weighted")
			ax1.set_yticks(np.arange(-2, 2, 1))
		else:
			ax1.set_yticklabels([])

		x = floor(index % columns)
		y = floor(index / columns)

		y = y + 1

		new_index = x + (y * columns)

		ax2 = fig.add_subplot(rows, columns, new_index)
		ax2.scatter(metal_massW_array, metal_lightW_array, label = model, color = color, facecolors='none')
		ax2.plot([0, 1], [0, 1], 'g--', transform=ax2.transAxes, color = color)
		ax2.set_xlim(-2.5, 1)
		ax2.set_ylim(-2.5, 1)
		ax2.set_aspect('equal')
		ax2.set_xlabel("[Z/H] - mass weighted")

		if x == 1:
			ax2.set_ylabel("[Z/H] - light weighted")
			ax2.set_yticks(np.arange(-2, 2, 1.0))
		else:
			ax2.set_yticklabels([])

		ax2.set_xticks(np.arange(-2, 2, 1.0))
		ax2.legend()

	plt.tight_layout()
	plt.show()
	fig.savefig("output/dissertation/data/massW_vs_lightW/figure.png")

