import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from astropy.io import fits
import math

# take second element for sort
def takeFirst(elem):
    return elem[0]

def reorder_and_average(parameter_array, csp_light_array):

	parameter_array_sorted = sorted(parameter_array)

	csp_light_sorted = [x for _,x in sorted(zip(parameter_array,csp_light_array))]

	index_array = []

	parameter_array_sorted_ = []

	for index, param in enumerate(parameter_array_sorted):

		if param not in parameter_array_sorted_ :
			parameter_array_sorted_.append(param)

			if index == 0:
				index_array.append(index)
			else:
				index_array.append(index)

		if index == len(parameter_array_sorted) -1:
			index_array.append(index + 1)

	csp_light_average = []

	for i in range(len(parameter_array_sorted_)):

		average_array = csp_light_sorted[index_array[i] : index_array[i+1]]

		csp_light_average.append(np.mean(average_array))

	return parameter_array_sorted_, csp_light_average

paths = [
	"output/dissertation/MASTAR_TH_VMPL7_KR",
	"output/dissertation/MASTAR_TH_VMPL9_KR",
	"output/dissertation/MASTAR_TH_VMPL11_KR",
	"output/dissertation/MASTAR_E_VMPL7_KR",
	"output/dissertation/CONROY_E_KR/downgraded",
]

for path in paths:
	 
 	#Get every file in the path of the folder and store in a list
	files = [os.path.join(path, f) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

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

	fig = plt.figure(figsize=(20,10))	

	#age_bins = [[0, 2], [2, 4], [4, 6], [6, 8], [8, 10], [10, 1000]]
	#age_bins = [[0, 1], [1, 6], [6, 1000]]
	age_bins = [[0, 1], [1, 3], [3, 5], [5, 7], [7, 1000]]

	columns = len(age_bins)
	rows = 2
	index = 1

	object_checker_array = []
	object_full_array    = []

	for file in files:

		with fits.open(file, memmap=True) as hdul:

			age_value = 10**hdul[1].header["age_lightW"]
			head, tail = os.path.split(file)

			object_full_array.append((tail, age_value))

	for age_bin_index in range(len(age_bins)):

		age_array       = []
		metal_array     = []
		csp_light_array = []

		counter = 0
		other_counter = 0

		for file in files:

			with fits.open(file, memmap=True) as hdul:

				age_value = 10**hdul[1].header["age_lightW"]
				#print(age_value, other_counter)

				if age_bins[age_bin_index][0] <= age_value < age_bins[age_bin_index][1]: 

					csp_age   = np.ndarray(hdul[1].header['ssp_number'])
					csp_Z     = np.ndarray(hdul[1].header['ssp_number'])
					csp_light = np.ndarray(hdul[1].header['ssp_number'])
					csp_mass  = np.ndarray(hdul[1].header['ssp_number'])

					for i in range(len(csp_age)):
						csp_age[i]   = hdul[1].header['log_age_ssp_' + str(i)]
						csp_Z[i]     = hdul[1].header['metal_ssp_' + str(i)]
						csp_light[i] = hdul[1].header['weightLight_ssp_' + str(i)]
						csp_mass[i]  = hdul[1].header['weightMass_ssp_' + str(i)]

					for i, a in enumerate(csp_age):
						#print(csp_light[i])
						age_array.append(10**a)

					for i, m in enumerate(csp_Z):
						metal_array.append(m)

					for i in csp_light:
						csp_light_array.append(i)

					counter = counter + 1

					head, tail = os.path.split(file)

					#print(tail)

					object_checker_array.append((tail, age_value))

				other_counter = other_counter + 1

		fig.suptitle("Model " + model + " with associated distributions of Usher et al data")

		#print(counter, other_counter)

		alpha = 0.1

		ax1 = fig.add_subplot(rows, columns, index)
		ax2 = fig.add_subplot(rows, columns, index + columns)

		ax1.bar(age_array, csp_light_array, width=0.25, align='center', alpha = alpha, color = color, zorder = 0, label = 'Age distribution', edgecolor='orange')
		ax1.scatter(age_array, csp_light_array, color = 'orange', alpha = alpha, zorder = 1, label = 'Age distribution')
		age_array, csp_light_average = reorder_and_average(age_array, csp_light_array)
		ax1.plot(age_array, csp_light_average, color = 'black', linewidth = 2, zorder = 2, label = "Average age distribution")
		
		ax1.set_xlabel('Lookback time (Gyr)')
		if index == 1:
			ax1.set_ylabel('Frequency')
		else:
			ax1.set_yticklabels([])
		ax1.set_ylim(0, 1)
		
		title = str(age_bins[index -1][0]) + r'$\leq Age$'
		if index != len(age_bins):
			title = title + " < " + str(age_bins[index -1][1]) + " Gyr"
		else:
			title = title + " Gyr"  
		title = title + " (" + str(counter) + " spectra)"
		ax1.set_title(title)

		try:
			ax1.legend()
		except(IndexError):
			pass

		ax2.bar(metal_array, csp_light_array, width=0.1, align='center', alpha=alpha, color = color, zorder = 0, label = 'Metallicity distribution', edgecolor='orange')
		ax2.scatter(metal_array, csp_light_array, color = 'orange', alpha = alpha, zorder = 1, label = 'Metallicity distribution')
		metal_array, csp_light_average = reorder_and_average(metal_array, csp_light_array)
		ax2.plot(metal_array, csp_light_average, color = 'black', linewidth = 2, zorder = 2, label = "Average metallicity distribution")
		
		ax2.set_xlabel('[Z/H] (dex)')
		if index == 1 :
			ax2.set_ylabel('Frequency')
		else:
			ax2.set_yticklabels([])

		ax2.set_ylim(0, 1)

		index = index + 1

	plt.tight_layout(rect=[0, 0.03, 1, 0.95])
	
	try:
		display_plot = sys.argv[1].lower() == 'true'
	except:
		display_plot = False

	if display_plot:
		plt.show()

	fig.savefig("output/dissertation/data/look_backtime/" + model + "_lookback_time.png")

	print("Missing values in second list:", (set(object_full_array).difference(object_checker_array)))



