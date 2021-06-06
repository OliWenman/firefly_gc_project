import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from astropy.io import fits
import math
from math import floor, ceil

plt.rcParams.update({'font.size': 18})

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

	fig   = plt.figure(figsize=(11.69*factor, 5*factor))	
	index = 0

	columns = len(paths)
	rows    = 1

	ax1 = fig.add_subplot(1, 1, 1)

	#Loop through the model data
	for path in paths:

		index = index + 1

		if path.split('/')[-1] == "downgraded" and not "CONROY" in path.split('/')[-2]:
			downgraded_to_conroy = "- downgraded"
		else:
			downgraded_to_conroy = ""
	 
	 	#Get every file in the path of the folder and store in a list
		files = [os.path.join(path, f) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

		EBV_array = []

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

			EBV  = float(hdul[1].header['EBV'])

			age_massW      = float(hdul[1].header['age_massW'])
			metal_massW    = float(hdul[1].header['metallicity_massW'])
			
			hdul.close()

			EBV_array.append(EBV)

			if EBV >= 0.15:
				to_print = model + " - " + object_name + " : E(B-V) = " + str(round(EBV, 2)) + ", age_massW = " + str(round(10**age_massW, 2)) + " Gyr, [Z/H] massW = " + str(round(metal_massW, 2))

				print(to_print)

		#ax1 = fig.add_subplot(rows, columns, index)
		ax1.hist(EBV_array, bins = 25, label = model, color = color, alpha = 0.7, histtype=u'step', linewidth = 4)
		#ax1.set_aspect('equal')
		#ax1.set_xlabel("Log Age (Gyr) - mass weighted")
		#ax1.set_xticks(np.arange(-1, 2, 1.0))
		ax1.legend()
		ax1.set_xlim(0, 0.35)
		ax1.set_ylim(0, 70)
		#ax1.set_title("Histogram of Firefly's derived E(B-V) values for 86 local galactic cluster data")

		ax1.set_xlabel("E(B-V)")
		ax1.set_ylabel("Frequency")
		"""
		if index == 1: 
			ax1.set_ylabel("Frequency")
		else:
			ax1.set_yticklabels([])
			ax1.tick_params(axis = "x", direction = "in")
		"""

	plt.tight_layout()
	plt.show()
	fig.savefig("output/dissertation/data/EBV/EBV_plot.png")