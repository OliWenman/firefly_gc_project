import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from astropy.io import fits
import math

plt.rcParams.update({'font.size': 18})

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

if __name__ == "__main__":

	#Get optinal argument to display plots.
	#If no argument given, default is True
	try:
		display_plot = sys.argv[1].lower() == 'true'
	except:
		display_plot = False

	try:
		absolute_value = sys.argv[2].lower() == 'true'
	except:
		absolute_value = False

	#light_weight = 0
	#mass_weight  = 1
	firefly_values_to_use = 1

	#Setup arrays to store calculated values
	model_array       = []
	lit_used_array    = []
	parameter_array   = []
	sample_size_array = []
	mean_array        = [[], []]
	medium_array      = [[], []]
	med_error_array   = [[], []]
	sigma_array       = [[], []]

	#Paths (relative to firefly directory) to read in the data processed by firefly

	paths = [
		"output/dissertation/MASTAR_TH_VMPL7_KR",
		"output/dissertation/MASTAR_TH_VMPL9_KR",
		"output/dissertation/MASTAR_TH_VMPL11_KR",
		"output/dissertation/MASTAR_E_VMPL7_KR",
		"output/dissertation/CONROY_E_KR/downgraded",
	]

	#Setup figure
	fig = plt.figure(figsize=(20,10))

	#Loop through the model data
	for path in paths:

		if path.split('/')[-1] == "downgraded" and not "CONROY" in path.split('/')[-2]:
			downgraded_to_conroy = "- downgraded"
		else:
			downgraded_to_conroy = ""
	 
	 	#Get every file in the path of the folder and store in a list
		files = [os.path.join(path, f) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

		#List of the literature data to read
		lit_files = [
			#"UsherGC.txt",
			"DeAngeli_GB.txt", 
			"DeAngeli_HST.txt",
		]

		#Initialise index for plotting graph
		index = 1

		"""
		if "MPL7" in path:
			version = "(MPL7)"
		elif "VMPL9" in path:
			version = "(VMPL9)"
		else:
			version = ""
		"""

		#Loop through literature data and plot/calculate stats
		for lit_file in lit_files:

			#Read in the literature data
			lit_table = pd.read_table(os.path.join(os.getcwd(), "output", "dissertation","literature_values", lit_file), 
									  delim_whitespace= True)

			#Initialise data of model and lit to be stored
			lit_age_array      = []
			lit_age_up_array   = []
			lit_age_low_array = []

			lit_metal_array     = []
			lit_metal_up_array  = []
			lit_metal_low_array = []

			#[0] is parameter_lightW
			#[1] is parameter_massW
			model_age_array     = [[], []]
			model_age_up_array  = [[], []]
			model_age_low_array = [[], []]
			
			model_metal_array     = [[], []]
			model_metal_up_array  = [[], []]
			model_metal_low_array = [[], []]

			#Loop through files of processed firefly output 
			for file in files:

				path_, file_ = os.path.split(file)

				#Get just the name of the object observed
				object_name = file_[6:file_.find("_")]
				
				#Search through literature table to see if it contains the spectra
				lit_values = lit_table.loc[lit_table['ID'] == object_name]

				#If the literature data doesn't contain the object, move to next file
				if lit_values.shape[0] == 0:
					continue

				#Check which file is being used and extract the data
				if lit_file == "UsherGC.txt":

					lit_age       = float(math.log10(lit_values['Age']))
					lit_age_error = 0

					lit_metal       = float(lit_values['[Fe/H]'])
					lit_metal_error = 0
				
				#DeAngeli_GB and HST have 2 values for metal and age. Have taken the average?
				elif lit_file == "DeAngeli_GB.txt":
					
					lit_age1      = float(lit_values['Age1'])
					lit_age2      = float(lit_values['Age2'])
					lit_age_error = abs(lit_age1 - lit_age2)
					
					lit_metal1      = float(lit_values['[Fe/H]zw'])
					lit_metal2      = float(lit_values['[Fe/H]cg'])
					lit_metal_error = abs(lit_metal1 - lit_metal2)

					lit_age   = (lit_age1 + lit_age2)/2
					lit_metal = (lit_metal1 + lit_metal2) / 2
					
				elif lit_file == "DeAngeli_HST.txt":
					
					lit_age1      = float(lit_values['Age1'])
					lit_age2      = float(lit_values['Age2'])
					lit_age_error = abs(lit_age1 - lit_age2)

					lit_metal1      = float(lit_values['[Fe/H]zw'])
					lit_metal2      = float(lit_values['[Fe/H]cg'])
					lit_metal_error = abs(lit_metal1 - lit_metal2)

					lit_age   = (lit_age1 + lit_age2)/2
					lit_metal = (lit_metal1 + lit_metal2) / 2

				#Extract the model data
				hdul = fits.open(file)

				model_age_lightW     = float(hdul[1].header['age_lightW'])
				model_age_lightW_up  = float(hdul[1].header['age_lightW_up_1sig'])
				model_age_lightW_low = float(hdul[1].header['age_lightW_low_1sig'])

				model_metal_lightW     = float(hdul[1].header['metallicity_lightW'])
				model_metal_lightW_up  = float(hdul[1].header['metallicity_lightW_up_1sig'])
				model_metal_lightW_low = float(hdul[1].header['metallicity_lightW_low_1sig'])

				model_age_massW     = float(hdul[1].header['age_massW'])
				model_age_massW_up  = float(hdul[1].header['age_massW_up_1sig'])
				model_age_massW_low = float(hdul[1].header['age_massW_low_1sig'])

				model_metal_massW     = float(hdul[1].header['metallicity_massW'])
				model_metal_massW_up  = float(hdul[1].header['metallicity_massW_up_1sig'])
				model_metal_massW_low = float(hdul[1].header['metallicity_massW_low_1sig'])

				#model = hdul[1].header["MODEL"] + version + downgraded_to_conroy
				
				hdul.close()

				#Check if the lit data is nan. If it isn't, save the lit and model data to array.
				if not np.isnan(lit_metal) and not np.isnan(lit_age):

					lit_age_array.append(lit_age)
					lit_age_up_array.append(lit_age_error)
					lit_age_low_array.append(lit_age_error)

					lit_metal_array.append(lit_metal)
					lit_metal_up_array.append(lit_metal_error)
					lit_metal_low_array.append(lit_metal_error)

					model_age_array[0].append(model_age_lightW)
					model_age_up_array[0].append(model_age_lightW_up)
					model_age_low_array[0].append(model_age_lightW_low)
					
					model_metal_array[0].append(model_metal_lightW)
					model_metal_up_array[0].append(model_metal_lightW_up)
					model_metal_low_array[0].append(model_metal_lightW_low)

					model_age_array[1].append(model_age_massW)
					model_age_up_array[1].append(model_age_massW_up)
					model_age_low_array[1].append(model_age_massW_low)
					
					model_metal_array[1].append(model_metal_massW)
					model_metal_up_array[1].append(model_metal_massW_up)
					model_metal_low_array[1].append(model_metal_massW_low)

				else:
					print("Missing:",object_name)

			sample_size = len(model_age_array[0])

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

			#Convert lists to numpy arrays
			model_age_array = np.array(model_age_array)
			lit_age_array   = np.array(lit_age_array)
			model_metal_array   = np.array(model_metal_array)
			lit_metal_array     = np.array(lit_metal_array)

			#Setup plotting options
			columns     = len(lit_files)
			rows        = 1
			capsize     = 3
			min_bin     = -3
			max_bin     = 3
			bin_width   = 0.25

			#Name figure
			#fig.suptitle("Model " + model + " compared with literature values.", fontweight='bold')		

			#Plot scatter plot data
			#Age
			ax3 = fig.add_subplot(rows, columns, index)
			ax3.set_title(lit_file[:-4] + "(Sample size = " + str(sample_size) + ")", fontweight='bold')

			markersize = 10

			ax3.errorbar(lit_age_array, 
						 model_age_array[0], 
						 yerr = (np.array(model_age_array[0]) - np.array(model_age_low_array[0]), np.array(model_age_up_array[0]) - np.array(model_age_array[0])),
						 xerr = (lit_age_low_array, lit_age_up_array),
						 color = color, 
						 fmt='o', 
						 ecolor=color, 
						 alpha = 0.75, 
						 markerfacecolor='none', 
						 capsize=capsize,
						 label = 'LW - ' + model,
						 linewidth = 1,
						 markersize = markersize)

			try:
				ax3.errorbar(lit_age_array, 
							 model_age_array[1], 
							 yerr = (np.array(model_age_array[1]) - np.array(model_age_low_array[1]), np.array(model_age_up_array[1]) - np.array(model_age_array[1])),
							 xerr = (lit_age_low_array, lit_age_up_array),
							 color = lighten_color(color, 1.2), 
							 fmt='o', 
							 ecolor=lighten_color(color, 1.2), 
							 alpha = 0.75, 
							 markerfacecolor='none', 
							 capsize=capsize,
							 marker = 'v',
							 label = 'MW - ' + model,
							 linewidth = 1,
							 markersize = markersize)
			except:
				ax3.errorbar(lit_age_array, 
						 model_age_array[1], 
						 yerr = (np.array(model_age_array[1]) - np.array(model_age_low_array[1]), np.array(model_age_up_array[1]) - np.array(model_age_array[1])),
						 xerr = (lit_age_low_array, lit_age_up_array),
						 color = "black", 
						 fmt='o', 
						 ecolor= "black", 
						 alpha = 0.75, 
						 markerfacecolor='none', 
						 capsize=capsize,
						 marker = 'v',
						 label = 'MW - ' + model,
						 linewidth = 1,
						 markersize = markersize)
			
			if index == 1:
				ax3.set_ylabel("Log Age (Gyr) - Model")
			else:
				ax3.set_yticklabels([])
				ax3.tick_params(axis = "y", direction = "in")

			ax3.set_xlim(-0.2, 1.15)
			ax3.set_ylim(-0.2, 1.15)

			ax3.set_xlabel("Log Age (Gyr) - Literature")
			ax3.legend(framealpha= 0.5)
			ax3.grid("black")

			ax3.plot([-5, 5], [-5, 5], 'g--', color = "black")

			ax3.axis([0.7, 1.15, -0.2, 1.15])

			index = index +1


	fig.tight_layout(rect=[0, 0.03, 1, 0.95])

	#Save plots
	fig.savefig("output/dissertation/data/comparison_to_lit/zoom/zoom_comparison_to_lit_" + (model if len(paths) == 1 else "(combined)") + ".png")

	if display_plot:
		plt.show()

