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

if __name__ == "__main__":

	#Get optinal argument to display plots.
	#If no argument given, default is True
	try:
		display_plot = sys.argv[1].lower() == 'true'
	except:
		display_plot = False

	#Setup arrays to store calculated values
	model_array       = []
	lit_used_array    = []
	parameter_array   = []
	sample_size_array = []
	mean_array        = []
	medium_array      = []
	sigma_array       = []

	#Paths (relative to firefly directory) to read in the data processed by firefly

	paths = [
		"output/dissertation/MASTAR_TH_V0.2_KR",
		"output/dissertation/MASTAR_TH_V0.2_KR/downgraded",
		"output/dissertation/MASTAR_E_V0.2_KR",
		"output/dissertation/MASTAR_E_V0.2_KR/downgraded",
		"output/dissertation/CONROY_E_KR/downgraded",
	]

	#Loop through the model data
	for path in paths:

		if path.split('/')[-1] == "downgraded" and not "CONROY" in path.split('/')[-2]:
			downgraded_to_conroy = "(downgraded)"
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

		#Loop through literature data and plot/calculate stats
		for lit_file in lit_files:

			#Read in the literature data
			lit_table = pd.read_table(os.path.join(os.getcwd(), "output", "dissertation","literature_values", lit_file), 
									  delim_whitespace= True)

			#Initialise data of model and lit to be store
			lit_age_array   = []
			lit_metal_array = []

			model_age_array   = []
			model_metal_array = []

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
					lit_age   = float(lit_values['Age'])
					lit_metal = float(lit_values['[Fe/H]'])
				
				#DeAngeli_GB and HST have 2 values for metal and age. Have taken the average?
				elif lit_file == "DeAngeli_GB.txt":
					
					lit_age1   = float(10**lit_values['Age1'])
					lit_age2   = float(10**lit_values['Age2'])
					lit_metal1 = float(lit_values['[Fe/H]zw'])
					lit_metal2 = float(lit_values['[Fe/H]cg'])

					lit_age   = (lit_age1 + lit_age2)/2
					lit_metal = (lit_metal1 + lit_metal2) / 2
					
				elif lit_file == "DeAngeli_HST.txt":
					
					lit_age1   = float(10**lit_values['Age1'])
					lit_age2   = float(10**lit_values['Age2'])
					lit_metal1 = float(lit_values['[Fe/H]zw'])
					lit_metal2 = float(lit_values['[Fe/H]cg'])

					lit_age   = (lit_age1 + lit_age2)/2
					lit_metal = (lit_metal1 + lit_metal2) / 2

				#Extract the model data
				hdul        = fits.open(file)
				model_age   = float(10**hdul[1].header['age_lightW'])

				model_metal = float(hdul[1].header['metallicity_lightW'])
				model       = hdul[1].header["MODEL"] + downgraded_to_conroy
				
				hdul.close()

				#Check if the lit data is nan. If it isn't, save the lit and model data to array.
				if not np.isnan(lit_metal) and not np.isnan(lit_age):

					lit_age_array.append(lit_age)
					lit_metal_array.append(lit_metal)

					model_age_array.append(model_age)
					model_metal_array.append(model_metal)

				else:
					print("Missing:",object_name)

			#Check the type of model used and assign a colour
			if model.upper() == "E-CONROY":
				color = "red"
			elif model.upper() == "TH-MASTAR":
				color = "royalblue"
			elif model.upper() == "E-MASTAR":
				color = "navy"	

			#Convert age array to log
			lit_age_log_array   = [math.log10(i) for i in lit_age_array]
			model_age_log_array = [math.log10(i) for i in model_age_array]

			#Convert lists to numpy arrays
			model_age_log_array = np.array(model_age_log_array)
			lit_age_log_array   = np.array(lit_age_log_array)
			model_metal_array   = np.array(model_metal_array)
			lit_metal_array     = np.array(lit_metal_array)

			#Setup plotting options
			columns = 4
			rows    = len(lit_files)
			bins    = 15

			#Name figure
			fig.suptitle("Model " + model + " compared with literature values.", fontweight='bold')
			
			#Plot histogram data

			#Age
			ax1 = fig.add_subplot(rows, columns, index)
			ax1.hist(model_age_log_array - lit_age_log_array, bins = bins, color = color)
			xabs_max = abs(max(ax1.get_xlim(), key=abs))
			ax1.set_xlim(xmin=-xabs_max, xmax=xabs_max)
			ax1.set_ylabel("Frequency")

			if index > rows*columns -columns: 
				ax1.set_xlabel("Age (log Gyr) Model - Literature")

			ax1.annotate(lit_file[:-4], 
						 xy=(0, 0.5), 
						 xytext=(-ax1.yaxis.labelpad - 5, 0), 
						 xycoords=ax1.yaxis.label, 
						 textcoords='offset points',
						 size='large', 
						 ha='right', 
						 va='center',
						 fontweight='bold')

			#Metal
			index = index +1
			ax2 = fig.add_subplot(rows, columns, index)
			ax2.hist(model_metal_array - lit_metal_array, bins = bins, color = color)
			xabs_max = abs(max(ax2.get_xlim(), key=abs))
			ax2.set_xlim(xmin=-xabs_max, xmax=xabs_max)

			ax2.set_ylabel("Frequency")
			if index > rows*columns -columns: 
				ax2.set_xlabel("[Z/H] Model - Literature")

			#Plot scatter plot data
			#Age
			index = index +1
			ax3 = fig.add_subplot(rows, columns, index)
			ax3.scatter(lit_age_log_array, model_age_log_array, color = color)
			ax3.plot([0, 1], [0, 1], 'g--', transform=ax3.transAxes, color = color)
			ax3.set_ylabel("Log Age (Gyr) - Model")

			if index > rows*columns -columns: 
				ax3.set_xlabel("Log Age (Gyr) - Literature")

			#Metal
			index = index +1
			ax4 = fig.add_subplot(rows, columns, index)
			ax4.scatter(lit_metal_array, model_metal_array, color = color)
			ax4.plot([0, 1], [0, 1], 'g--', transform=ax4.transAxes, color = color)
			ax4.set_ylabel("[Z/H] - Model")

			if index > rows*columns -columns: 
				ax4.set_xlabel("[Z/H] - Literature")

			index = index +1
			
			#Calculate stats
			sample_size = len(model_age_array)
			print(len(model_age_array))
			print(len(model_metal_array))
			print(len(lit_age_array))
			print(len(lit_metal_array))

			mean_age     = np.mean(model_age_log_array - lit_age_log_array)
			mean_metal   = np.mean(model_metal_array - lit_metal_array)

			medium_age   = np.median(model_age_log_array - lit_age_log_array)
			medium_metal = np.median(model_metal_array - lit_metal_array)

			std_age      = np.std(model_age_log_array - lit_age_log_array)
			std_metal    = np.std(model_metal_array - lit_metal_array)

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

			mean_array.append(round(mean_age*100, n_decimals))
			mean_array.append(round(mean_metal*100, n_decimals))			

			medium_array.append(round(medium_age*100, n_decimals))
			medium_array.append(round(medium_metal*100, n_decimals))

			sigma_array.append(round(std_age*100, n_decimals))
			sigma_array.append(round(std_metal*100, n_decimals))

			#print("Medium Difference Age (log Gyr) =", medium_age)
			#print("Medium Difference [Z/H]         =", medium_metal)

			#print("STD Difference Age (log Gyr)    =", std_age)
			#print("STD Difference [Z/H]            =", std_metal)

		fig.tight_layout(rect=[0, 0.03, 1, 0.95])

		#Save plots
		fig.savefig("output/dissertation/data/" + model + "_comparison_to_lit.png")

		if display_plot:
			plt.show()

	#Create dataframe and save the stats in a table
	data = {'Model': model_array, 'Literature data': lit_used_array, 'Parameter':parameter_array, 'Sample size': sample_size_array,'Mean difference(%)': mean_array,'Medium difference(%)': medium_array, 'Standard deviation difference(%)': sigma_array}

	df = pd.DataFrame(data=data)

	df.to_csv("output/dissertation/data/table.csv", index = False, header=True)

	print(df.loc[(df['Literature data'] == 'UsherGC') & (df['Parameter'] == '[Z/H]')], "\n")
	print(df.loc[(df['Literature data'] == 'UsherGC') & (df['Parameter'] == 'Age (log Gyr)')], "\n")

	print(df.loc[(df['Literature data'] == 'DeAngeli_HST') & (df['Parameter'] == '[Z/H]')], "\n")
	print(df.loc[(df['Literature data'] == 'DeAngeli_HST') & (df['Parameter'] == 'Age (log Gyr)')], "\n")

	print(df.loc[(df['Literature data'] == 'DeAngeli_GB') & (df['Parameter'] == '[Z/H]')], "\n")
	print(df.loc[(df['Literature data'] == 'DeAngeli_GB') & (df['Parameter'] == 'Age (log Gyr)')], "\n")

	#print(df)

