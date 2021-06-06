"""
Script for analysing processed firefly data and comparing the models used
to literature values. 
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from astropy.io import fits
import math

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

		#"output/dissertation/MASTAR_TH_V0.2_KR/downgraded",
		"output/dissertation/MASTAR_E_VMPL7_KR",
		#"output/dissertation/MASTAR_E_V0.2_KR/downgraded",
		"output/dissertation/CONROY_E_KR/downgraded",
	]

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
			"UsherGC.txt", 
			"DeAngeli_GB.txt", 
			"DeAngeli_HST.txt"
		]

		#Setup figure
		fig = plt.figure(figsize=(20,10))

		#Initialise index for plotting graph
		index = 1

		if "MPL7" in path:
			version = "(MPL7)"
		elif "VMPL9" in path:
			version = "(VMPL9)"
		elif "VMPL11" in path:
			version = "(VMPL11)"

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

				model = hdul[1].header["MODEL"] + version + downgraded_to_conroy
				
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
			columns     = 4
			rows        = len(lit_files)
			capsize     = 3
			min_bin     = -3
			max_bin     = 3
			bin_width    = 0.25

			#bins = np.arange(min_bin, max_bin + bin_width, bin_width)
			bins = np.linspace(start = min_bin, stop= max_bin, num=26)
			#Name figure
			fig.suptitle("Model " + model + " compared with literature values.", fontweight='bold')
			
			#Plot histogram data

			#Age
			ax1 = fig.add_subplot(rows, columns, index)

			age_difference = np.array(model_age_array[firefly_values_to_use]) - np.array(lit_age_array)
			if absolute_value:
				age_difference = np.absolute(age_difference)

			age_difference_lightW = np.array(model_age_array[0]) - np.array(lit_age_array)
			age_difference_massW  = np.array(model_age_array[1]) - np.array(lit_age_array)

			ax1.hist(age_difference_lightW, bins = bins, color = color, alpha = 0.5, label = "LW - " + model)
			ax1.hist(age_difference_lightW, bins = bins, color = color, alpha = 0.75, histtype=u'step', linewidth = 3)
			ax1.hist(age_difference_massW, bins = bins, color = "black", alpha = 0.5, label = "MW - " + model)
			ax1.hist(age_difference_massW, bins = bins, color = "black", alpha = 0.75, histtype=u'step', linewidth = 3)
			xabs_max = abs(max(ax1.get_xlim(), key=abs))

			ax1.set_xlim(xmin=min_bin, xmax=max_bin)
			ax1.set_ylabel("Frequency")
			ax1.vlines(0, linestyle = "dashed", ymin=0, ymax=100)

			if index == 1:
				ax1.legend(framealpha = 0.5)

			if "Usher" in lit_file:
				ax1.set_ylim(0, 40)
			elif "DeAngeli_GB" in lit_file:
				ax1.set_ylim(0, 15)
			elif "DeAngeli_HST" in lit_file:
				ax1.set_ylim(0, 20)

			if index > rows*columns -columns: 
				ax1.set_xlabel("Age (log Gyr) Model - Literature")
			else:
				ax1.set_xticklabels([])
				ax1.tick_params(axis = "x", direction = "in")

			ax1.annotate(lit_file[:-4] + "\n(Sample size = " + str(sample_size) + ")", 
						 xy=(0, 0.5), 
						 xytext=(-ax1.yaxis.labelpad - 75, 0), 
						 xycoords=ax1.yaxis.label, 
						 textcoords='offset points',
						 size='large', 
						 ha='center', 
						 va='center',
						 fontweight='bold')

			#Plot scatter plot data
			#Age
			index = index +1
			ax3 = fig.add_subplot(rows, columns, index)
			ax3.plot([0, 1], [0, 1], 'g--', transform=ax3.transAxes, color = "black")

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
						 linewidth = 1)

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
							 linewidth = 1)
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
						 linewidth = 1)
			
			ax3.set_ylabel("Log Age (Gyr) - Model")
			ax3.set_xlim(-2, 1.3)
			ax3.set_ylim(-2, 1.3)

			if index > rows*columns -columns: 
				ax3.set_xlabel("Log Age (Gyr) - Literature")
			else:
				ax3.set_xticklabels([])
				ax3.tick_params(axis = "x", direction = "in")

			if index == 2:
				ax3.legend(framealpha= 0.5)

			#Metal
			index = index +1

			ax2 = fig.add_subplot(rows, columns, index)

			metal_difference = np.array(model_metal_array[firefly_values_to_use]) - np.array(lit_metal_array)

			if absolute_value:
				metal_difference = np.absolute(metal_difference)

			metal_difference_lightW = np.array(model_metal_array[0]) - np.array(lit_metal_array)
			metal_difference_massW  = np.array(model_metal_array[1]) - np.array(lit_metal_array)

			ax2.hist(metal_difference_lightW, bins = bins, color = color, alpha = 0.5, label = "LW - " + model)
			ax2.hist(metal_difference_lightW, bins = bins, color = color, alpha = 0.75, histtype=u'step', linewidth = 3)
			ax2.hist(metal_difference_massW, bins = bins, color = "black", alpha = 0.5,label = "MW - " + model)
			ax2.hist(metal_difference_massW, bins = bins, color = "black", alpha = 0.75, histtype=u'step', linewidth = 3)
			xabs_max = abs(max(ax2.get_xlim(), key=abs))
			ax2.set_xlim(xmin=-3, xmax=3)
			ax2.vlines(0, linestyle = "dashed", ymin=0, ymax=100)

			if index == 3:
				ax2.legend(framealpha = 0.5)

			ax2.set_ylabel("Frequency")
			if index > rows*columns -columns: 
				ax2.set_xlabel("[Z/H] Model - Literature")
			else:
				ax2.set_xticklabels([])
				ax2.tick_params(axis = "x", direction = "in")

			if "Usher" in lit_file:
				ax2.set_ylim(0, 30)
			elif "DeAngeli_GB" in lit_file:
				ax2.set_ylim(0, 12)
			elif "DeAngeli_HST" in lit_file:
				ax2.set_ylim(0, 18)

			#Metal
			index = index +1
			ax4 = fig.add_subplot(rows, columns, index)
			#ax4.scatter(lit_metal_array, model_metal_array, color = color)
			ax4.plot([0, 1], [0, 1], 'g--', transform=ax4.transAxes, color = "black")

			ax4.errorbar(lit_metal_array, 
						 model_metal_array[0], 
						 yerr = (np.array(model_metal_array[0]) - np.array(model_metal_low_array[0]), np.array(model_metal_up_array[0]) - np.array(model_metal_array[0])),
						 xerr = (lit_metal_low_array, lit_metal_up_array),
						 color = color, 
						 fmt='o', 
						 ecolor=color, 
						 alpha = 0.75, 
						 markerfacecolor='none', 
						 capsize=capsize,
						 label = 'LW - ' + model)

			try:
				ax4.errorbar(lit_metal_array, 
							 model_metal_array[1], 
							 yerr = (np.array(model_metal_array[1]) - np.array(model_metal_low_array[1]), np.array(model_metal_up_array[1]) - np.array(model_metal_array[1])),
							 xerr = (lit_metal_low_array, lit_metal_up_array),
							 color = lighten_color(color, 1.3), 
							 fmt='o', 
							 ecolor= lighten_color(color, 1.3), 
							 alpha = 0.75, 
							 markerfacecolor='none', 
							 capsize=capsize,
							 marker = 'v',
							 label = 'MW - ' + model)
			except:
				ax4.errorbar(lit_metal_array, 
							 model_metal_array[1], 
							 yerr = (np.array(model_metal_array[1]) - np.array(model_metal_low_array[1]), np.array(model_metal_up_array[1]) - np.array(model_metal_array[1])),
							 xerr = (lit_metal_low_array, lit_metal_up_array),
							 color = "black", 
							 fmt='o', 
							 ecolor= "black", 
							 alpha = 0.75, 
							 markerfacecolor='none', 
							 capsize=capsize,
							 marker = 'v',
							 label = 'MW - ' + model)
			
			ax4.set_ylabel("[Z/H] - Model")
			ax4.set_xlim(-2.5, 0.5)
			ax4.set_ylim(-2.5, 0.5)

			if index == 4:
				ax4.legend(framealpha = 0.5)


			if index > rows*columns -columns: 
				ax4.set_xlabel("[Z/H] - Literature")
			else:
				ax4.set_xticklabels([])
				ax4.tick_params(axis = "x", direction = "in")

			index = index +1
			
			#Calculate stats
			#######################################################################

			#Light weight

			mean_age     = np.mean(age_difference_lightW)
			mean_metal   = np.mean(metal_difference_lightW)

			medium_age   = np.median(age_difference_lightW)
			medium_metal = np.median(metal_difference_lightW)

			#Medium error
			age_error_firefly      = np.maximum(model_age_up_array[0], model_age_low_array[0])
			age_error_literature   = np.maximum(lit_age_up_array, lit_age_low_array) 
			metal_error_firefly    = np.maximum(model_metal_up_array[0], model_metal_low_array[0])
			metal_error_literature = np.maximum(lit_metal_up_array, lit_metal_low_array) 
			med_error_age          = np.sqrt((np.median(age_error_firefly/model_age_array[0]))**2 + (np.median(age_error_literature/lit_age_array))**2) * mean_age
			med_error_metal        = np.sqrt((np.median(metal_error_firefly/model_metal_array[0]))**2 + (np.median(metal_error_literature/lit_metal_array))**2) * mean_metal

			std_age      = np.std(age_difference_lightW)
			std_metal    = np.std(metal_difference_lightW)

			#Save data to arrays
			sample_size_array.append(sample_size)
			sample_size_array.append(sample_size)

			model_array.append(model)
			model_array.append(model)

			lit_used_array.append(lit_file[:-4])
			lit_used_array.append(lit_file[:-4])

			parameter_array.append("Age (log Gyr)")
			parameter_array.append("[Z/H]")

			n_decimals= 2

			mean_array[0].append(round(mean_age, n_decimals))
			mean_array[0].append(round(mean_metal, n_decimals))			

			medium_array[0].append(round(medium_age, n_decimals))
			medium_array[0].append(round(medium_metal, n_decimals))

			med_error_array[0].append(abs(round(med_error_age, n_decimals)))
			med_error_array[0].append(abs(round(med_error_metal, n_decimals)))

			sigma_array[0].append(round(std_age, n_decimals))
			sigma_array[0].append(round(std_metal, n_decimals))

			#########################################################################################
			#Mass weight

			mean_age     = np.mean(age_difference_massW)
			mean_metal   = np.mean(metal_difference_massW)

			medium_age   = np.median(age_difference_massW)
			medium_metal = np.median(metal_difference_massW)

			age_error_firefly      = np.maximum(model_age_up_array[1], model_age_low_array[1])
			age_error_literature   = np.maximum(lit_age_up_array, lit_age_low_array) 
			metal_error_firefly    = np.maximum(model_metal_up_array[1], model_metal_low_array[1])
			metal_error_literature = np.maximum(lit_metal_up_array, lit_metal_low_array) 
			med_error_age          = np.sqrt((np.median(age_error_firefly/model_age_array[1]))**2 + (np.median(age_error_literature/lit_age_array))**2) * mean_age
			med_error_metal        = np.sqrt((np.median(metal_error_firefly/model_metal_array[1]))**2 + (np.median(metal_error_literature/lit_metal_array))**2) * mean_metal

			std_age      = np.std(age_difference_massW)
			std_metal    = np.std(metal_difference_massW)

			mean_array[1].append(round(mean_age, n_decimals))
			mean_array[1].append(round(mean_metal, n_decimals))			

			medium_array[1].append(round(medium_age, n_decimals))
			medium_array[1].append(round(medium_metal, n_decimals))

			med_error_array[1].append(abs(round(med_error_age, n_decimals)))
			med_error_array[1].append(abs(round(med_error_metal, n_decimals)))

			sigma_array[1].append(round(std_age, n_decimals))
			sigma_array[1].append(round(std_metal, n_decimals))

		fig.tight_layout(rect=[0, 0.03, 1, 0.95])

		#Save plots

		fig.savefig("output/dissertation/data/comparison_to_lit/" + ("absolute/" if absolute_value else "") + model + "_comparison_to_lit" + ("_absolute" if absolute_value else "") + ".png")

		if display_plot:
			plt.show()

	#Create dataframe and save the stats in a table
	data = {'Model': model_array, 'Literature data': lit_used_array, 'Parameter':parameter_array, 'Sample size': sample_size_array, 'Medium difference': medium_array[0], 'Medium error': med_error_array[0],'Standard deviation difference': sigma_array[0]}
	df = pd.DataFrame(data=data)
	df.to_csv("output/dissertation/data/comparison_to_lit/" + ("absolute/" if absolute_value else "") + "table_lightW" + ("_absolute" if absolute_value else "") +".csv", index = False, header=True)

	data = {'Model': model_array, 'Literature data': lit_used_array, 'Parameter':parameter_array, 'Sample size': sample_size_array, 'Medium difference': medium_array[1], 'Medium error': med_error_array[1],'Standard deviation difference': sigma_array[1]}
	df = pd.DataFrame(data=data)
	df.to_csv("output/dissertation/data/comparison_to_lit/" + ("absolute/" if absolute_value else "") + "table_massW" + ("_absolute" if absolute_value else "") +".csv", index = False, header=True)

	print(df.loc[(df['Literature data'] == 'UsherGC') & (df['Parameter'] == '[Z/H]')], "\n")
	print(df.loc[(df['Literature data'] == 'UsherGC') & (df['Parameter'] == 'Age (log Gyr)')], "\n")

	print(df.loc[(df['Literature data'] == 'DeAngeli_HST') & (df['Parameter'] == '[Z/H]')], "\n")
	print(df.loc[(df['Literature data'] == 'DeAngeli_HST') & (df['Parameter'] == 'Age (log Gyr)')], "\n")

	print(df.loc[(df['Literature data'] == 'DeAngeli_GB') & (df['Parameter'] == '[Z/H]')], "\n")
	print(df.loc[(df['Literature data'] == 'DeAngeli_GB') & (df['Parameter'] == 'Age (log Gyr)')], "\n")

	#print(df)

