import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from astropy.io import fits
import math
from math import floor, ceil

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

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

	fig   = plt.figure(figsize=(8, 6))	
	index = 0

	columns = len(paths)
	rows    = 2

	std_values = {}

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
			model = "Conroy"

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

		min_bin     = -1.5
		max_bin     = 1.5
		bin_width    = 0.25

		#bins = np.arange(min_bin, max_bin + bin_width, bin_width)
		bins = np.linspace(start = min_bin, stop= max_bin, num=26)

		ax1 = fig.add_subplot(rows, columns, index)
		ax1.hist(np.array(age_massW_array) - np.array(age_lightW_array), label = model, color = color, bins = bins)
		std_values[model + "_age"] = np.std(np.array(age_massW_array) - np.array(age_lightW_array))
		
		ax1.legend(prop={'size': 12}, loc = 'upper right')
		ax1.grid()

		ax1.set_ylim(0, 80)
		ax1.set_aspect(1 / ax1.get_data_ratio())
		ax1.set_xlabel("Age (Log Gyr) MW - LW")
		if index == 1:
			ax1.set_ylabel("Frequency")		
		else:
			ax1.set_yticklabels([])

		x = floor(index % columns)
		y = floor(index / columns)

		y = y + 1

		new_index = x + (y * columns)

		ax2 = fig.add_subplot(rows, columns, new_index)
		ax2.hist(np.array(metal_massW_array) - np.array(metal_lightW_array), label = model, color = color, bins = bins)
		std_values[model + "_metal"] = np.std(np.array(metal_massW_array) - np.array(metal_lightW_array))

		#ax2.plot([0, 1], [0, 1], 'g--', transform=ax2.transAxes, color = color)
		#ax2.set_xlim(-2.5, 1)
		ax2.set_ylim(0, 55)
		ax2.set_aspect(1 / ax2.get_data_ratio())
		#ax2.set_xlabel("[Z/H] - MW")

		ax2.set_xlabel("[Z/H] MW - LW")
		if index == 1:
			ax2.set_ylabel("Frequency")		
		else:
			ax2.set_yticklabels([])

		#ax2.set_xticks(np.arange(-2, 2, 1.0))
		#ax2.legend(prop={'size': 12}, loc = 'lower right')
		ax2.legend(prop={'size': 12}, loc = 'upper right')
		ax2.grid()

	print(std_values)

	object_outside_std = []

	#Find associated spectra outside of the std_value range
	###########################################################################################################################

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
			model = "Conroy"

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
			object_name = file_[6:-5]

			#Extract the model data
			hdul = fits.open(file)

			age_lightW     = float(hdul[1].header['age_lightW'])
			metal_lightW   = float(hdul[1].header['metallicity_lightW'])
			age_massW      = float(hdul[1].header['age_massW'])
			metal_massW    = float(hdul[1].header['metallicity_massW'])
			
			hdul.close()
			"""
			if std_values[model + "_age"] <= (age_massW - age_lightW):
				print(object_name)

			if std_values[model + "_metal"] <= (metal_massW - metal_lightW):
				print(object_name)
			"""

			sigma = 2
			if (std_values[model + "_metal"]*sigma) <= (abs(metal_massW - metal_lightW)) or (std_values[model + "_age"]*sigma) <= (abs(age_massW - age_lightW)):
				object_outside_std.append(object_name)


	object_location = {
		"NGC0104":	"MW",
		"NGC0121":	"SMC",
		"NGC0330":	"SMC",
		"NGC0361":	"SMC",
		"NGC0362":	"MW",
		"NGC0416":	"SMC",
		"NGC0419":	"SMC",
		"NGC1261":	"MW",
		"NGC1786":	"LMC",
		"NGC1783":	"LMC",
		"NGC1846":	"LMC",
		"NGC1850":	"LMC",
		"NGC1856":	"LMC",
		"NGC1866":	"LMC",
		"NGC1851":	"MW",
		"NGC1868":	"LMC",
		"NGC1898":	"LMC",
		"NGC1916":	"LMC",
		"NGC1904":	"MW",
		"NGC1978":	"LMC",
		"NGC2004":	"LMC",
		"NGC2019":	"LMC",
		"NGC2100":	"LMC",
		"NGC2136":	"LMC",
		"NGC2808":	"MW",
		"NGC3201":	"MW",
		"NGC4147":	"MW",
		"NGC4590":	"MW",
		"NGC4833":	"MW",
		"NGC5024":	"MW",
		"NGC5139":	"MW",
		"NGC5272":	"MW",
		"NGC5286":	"MW",
		"NGC5634":	"MW",
		"NGC5694":	"MW",
		"NGC5824":	"MW",
		"NGC5904":	"MW",
		"NGC5927":	"MW",
		"NGC5986":	"MW",
		"NGC6093":	"MW",
		"NGC6121":	"MW",
		"NGC6139":	"MW",
		"NGC6171":	"MW",
		"NGC6218":	"MW",
		"NGC6254":	"MW",
		"NGC6266":	"MW",
		"NGC6273":	"MW",
		"NGC6284":	"MW",
		"NGC6293":	"MW",
		"NGC6304":	"MW",
		"NGC6316":	"MW",
		"NGC6333":	"MW",
		"NGC6342":	"MW",
		"NGC6356":	"MW",
		"NGC6352":	"MW",
		"NGC6362":	"MW",
		"NGC6388":	"MW",
		"NGC6397":	"MW",
		"NGC6440":	"MW",
		"NGC6441":	"MW",
		"NGC6522":	"MW",
		"NGC6528":	"MW",
		"NGC6541":	"MW",
		"NGC6553":	"MW",
		"NGC6569":	"MW",
		"NGC6584":	"MW",
		"NGC6624":	"MW",
		"NGC6637":	"MW",
		"NGC6652":	"MW",
		"NGC6656":	"MW",
		"NGC6681":	"MW",
		"NGC6715":	"MW",
		"NGC6717":	"MW",
		"NGC6723":	"MW",
		"NGC6752":	"MW",
		"NGC6809":	"MW",
		"NGC6838":	"MW",
		"NGC6864":	"MW",
		"NGC6934":	"MW",
		"NGC7006":	"MW",
		"NGC7078":	"MW",
		"NGC7089":	"MW",
		"NGC7099":	"MW",
	}

	#print("\n", object_outside_std)
	object_outside_std_counter = {i:object_outside_std.count(i) for i in object_outside_std}
	print("\n", object_outside_std_counter)

	final_data = []
	for obj in object_outside_std:

		if obj[:obj.find("_")] in object_location:
			final_data.append(object_location[obj[:obj.find("_")]] + " - " + obj) 

	final_data_count = {i:final_data.count(i) for i in final_data}
	print("\n")

	for i in final_data_count:
		print(i, ",", final_data_count[i])

	#for key in final_data_count:
	#	print

	plt.tight_layout()
	plt.show()
	fig.savefig("output/dissertation/data/massW_vs_lightW/LW_MW_hist.png")

