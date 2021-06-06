import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from astropy.io import fits
from math import log10, floor

def round_sig(x, sig=2):
	if x == 0:
		return 0

	return round(x, sig-int(floor(log10(abs(x))))-1)

if __name__ == "__main__":

	#light_weight = 0
	#mass_weight  = 1
	firefly_values_to_use = 1

	#Setup arrays to store calculated values
	model_values_array = [[],[]]
	lit_values_array   = [[],[]]

	paths = [
		"output/dissertation/MASTAR_TH_VMPL7_KR",
		"output/dissertation/MASTAR_TH_VMPL9_KR",
		"output/dissertation/MASTAR_TH_VMPL11_KR",
		"output/dissertation/MASTAR_E_VMPL7_KR",
		"output/dissertation/CONROY_E_KR/downgraded",
	]

	#List of the literature data to read
	lit_files = [
		#"UsherGC.txt", 
		#"DeAngeli_GB.txt", 
		#"DeAngeli_HST.txt"
	]

	data = {}

	#Loop through the model data
	for path in paths:

		if path.split('/')[-1] == "downgraded" and not "CONROY" in path.split('/')[-2]:
			downgraded_to_conroy = "- downgraded"
		else:
			downgraded_to_conroy = ""
	 
	 	#Get every file in the path of the folder and store in a list
		files = [os.path.join(path, f) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

		model_values_array = [[],[]]
		objects_array      = []
		lit_array          = [[],[]]

		for file in files:

			#Extract the model data
			hdul = fits.open(file)

			path_, file_ = os.path.split(file)

			observed_object = file_[6:-5]

			observed_object = observed_object.replace("_", "(") + ")"

			sig_fig = 3

			model_age_lightW     = round_sig(10**float(hdul[1].header['age_lightW']), sig_fig)
			model_age_lightW_up  = round_sig(10**float(hdul[1].header['age_lightW_up_1sig']), sig_fig)
			model_age_lightW_low = round_sig(10**float(hdul[1].header['age_lightW_low_1sig']), sig_fig)

			model_metal_lightW     = round_sig(float(hdul[1].header['metallicity_lightW']), sig_fig)
			model_metal_lightW_up  = round_sig(float(hdul[1].header['metallicity_lightW_up_1sig']), sig_fig)
			model_metal_lightW_low = round_sig(float(hdul[1].header['metallicity_lightW_low_1sig']), sig_fig)

			model_age_massW     = round_sig(10**float(hdul[1].header['age_massW']), sig_fig)
			model_age_massW_up  = round_sig(10**float(hdul[1].header['age_massW_up_1sig']), sig_fig)
			model_age_massW_low = round_sig(10**float(hdul[1].header['age_massW_low_1sig']), sig_fig)

			model_metal_massW     = round_sig(float(hdul[1].header['metallicity_massW']), sig_fig)
			model_metal_massW_up  = round_sig(float(hdul[1].header['metallicity_massW_up_1sig']), sig_fig)
			model_metal_massW_low = round_sig(float(hdul[1].header['metallicity_massW_low_1sig']), sig_fig)

			#model = hdul[1].header["MODEL"] + version + downgraded_to_conroy
			
			hdul.close()

			error_age   =  round_sig(np.max([model_age_massW_up, model_age_massW_low]) - model_age_massW, sig_fig)
			error_metal =  round_sig(np.max([model_metal_massW_up, model_metal_massW_low]) - model_metal_massW, sig_fig)

			model_age = str(model_age_massW) + r' $\pm$ ' + str(error_age)
			model_metal = str(model_metal_massW) + r' $\pm$ ' + str(error_metal)

			model_values_array[0].append(model_age)
			model_values_array[1].append(model_metal)

			if paths[0] == path:
				objects_array.append(observed_object)


		if paths[0] == path:
			data[("Object", "", "")] = objects_array


		#Check the type of model used and assign a colour
		if "CONROY_E" in path.upper():
			model = "Conroy"
			version = ""

		elif "MASTAR_TH" in path.upper():
			model = "Th-MaStar"

			if "VMPL7" in path.upper():
				version = "MPL7"

			elif "VMPL9" in path.upper():
				version = "MPL9"

			elif "VMPL11" in path.upper():
				version = "MPL11"

		elif "MASTAR_E" in path.upper():
			model = "E-MaStar"
			version = ""	

		data[(model, version, "Age (Gyr)")] = model_values_array[0]
		data[(model, version,  "[Z/H]")]     = model_values_array[1]

		"""
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
		"""
		
	#print(data)	

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


		#Loop through files of processed firefly output 
		for file in files:

			path_, file_ = os.path.split(file)

			#Get just the name of the object observed
			object_name = file_[6:file_.find("_")]
			
			#Search through literature table to see if it contains the spectra
			lit_values = lit_table.loc[lit_table['ID'] == object_name]

			#If the literature data doesn't contain the object, move to next file
			if lit_values.shape[0] == 0:
				
				lit_age_array.append(None)
				lit_metal_array.append(None)
				continue

			#Check which file is being used and extract the data
			if lit_file == "UsherGC.txt":

				lit_age       = float(math.log10(lit_values['Age']))
				lit_age_error = 0

				lit_metal       = float(lit_values['[Fe/H]'])
				lit_metal_error = 0
			
			#DeAngeli_GB and HST have 2 values for metal and age. Have taken the average?
			elif lit_file == "DeAngeli_GB.txt":
				
				#print(lit_values['Age1'], type(lit_values['Age1']))
				lit_age1      = lit_values['Age1'].astype(float)
				lit_age2      = lit_values['Age2'].astype(float)
				lit_age_error = abs(lit_age1 - lit_age2)
				
				lit_metal1      = lit_values['[Fe/H]zw'].astype(float)
				lit_metal2      = lit_values['[Fe/H]cg'].astype(float)
				lit_metal_error = abs(lit_metal1 - lit_metal2)

				lit_age   = (lit_age1 + lit_age2) / 2
				lit_metal = (lit_metal1 + lit_metal2) / 2
				
			elif lit_file == "DeAngeli_HST.txt":
				
				lit_age1      = lit_values['Age1'].astype(float)
				lit_age2      = lit_values['Age2'].astype(float)
				lit_age_error = abs(lit_age1 - lit_age2)

				lit_metal1      = lit_values['[Fe/H]zw'].astype(float)
				lit_metal2      = lit_values['[Fe/H]cg'].astype(float)
				lit_metal_error = abs(lit_metal1 - lit_metal2)

				lit_age   = (lit_age1 + lit_age2) / 2
				lit_metal = (lit_metal1 + lit_metal2) / 2

			lit_age_array.append(round(float(10**lit_age), 2))
			lit_metal_array.append(round(float(lit_metal), 2))

		print(lit_file[:-4] + " Age (Gyr)", len(lit_age_array))

		data[(lit_file[:-4], "Age (Gyr)", "")] = lit_age_array
		data[(lit_file[:-4], "[Z/H]", "")]     = lit_metal_array

	#columns = pd.MultiIndex.from_product([['MaStar-Th','Literature'], ['Age(Gyr)','[Z/H]']])

	df = pd.DataFrame(data=data)

	print(df)

	df.to_csv("output/dissertation/data/comparison_to_lit/full_table_results.csv", index = False, header=True)

	"""		
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
	"""