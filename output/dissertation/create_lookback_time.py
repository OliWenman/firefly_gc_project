import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from astropy.io import fits
import math

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

#age_bins = [[0, 2], [2, 4], [4, 6], [6, 8], [8, 10], [10, 1000]]
#age_bins = [[0, 1], [1, 6], [6, 1000]]
#age_bins = [[0, 0.5], [0.5, 3], [3, 5], [5, 10000]]
#age_bins  = [[0, 10000]]
log_true = True
#Log
#age_bins  = [[-2, -0.5], [-0.5, 0], [0, 0.5], [0.5, 1], [1, 1.5]]
#age_bins   = [[-2, 10000]]
age_bins   = [[-2, -1], [-1, 0], [0, 1], [1, 1000]]
metal_bins = [[-3, -2],   [-2, -1],    [-1, 0],    [0, 1]]
#metal_bins = [[-3, 2]] 

columns = len(age_bins)
rows = 1

paths = [
	"output/dissertation/MASTAR_TH_VMPL7_KR",
	"output/dissertation/MASTAR_TH_VMPL9_KR",
	"output/dissertation/MASTAR_TH_VMPL11_KR",
	"output/dissertation/MASTAR_E_VMPL7_KR",
	"output/dissertation/CONROY_E_KR/downgraded",
]


alpha = 0.5
width = 0.05


# take second element for sort
def takeFirst(elem):
    return elem[0]

def stack_bars(x_list, y_list):

	#print(1, x_list)
	#print(2, y_list, "\n")

	#Sort the lists in order according to the x_list
	y_list = [x for _,x in sorted(zip(x_list,y_list))]
	x_list = sorted(x_list)

	#print(3, x_list)
	#print(4, y_list, "\n")

	prev_value = None

	new_x = []
	new_y = []

	n_spectra = []

	counter = -1
	for i, x in enumerate(x_list):

		#If same x value, add on to y value 
		if x == prev_value:
			
			new_y[counter] = new_y[counter] + y_list[i]
			n_spectra[counter] = n_spectra[counter] + 1
		#If a new x value, store as new previous value and append new values to x and y lists 
		else:
			n_spectra.append(1)
			prev_value = x
			new_x.append(x)
			new_y.append(y_list[i])
			counter = counter + 1

	#print(5, x_list)
	#print(6, y_list, "\n")
	print(new_x, "\n")
	print(new_y, "\n")
	print(n_spectra, "\n")
	print(len(n_spectra), len(new_x))

	new_y_ = [b / m for b,m in zip(new_y, n_spectra)]

	return new_x, new_y_

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

def findMinDiff(arr):

	n = len(arr)
	# Initialize difference as infinite
	diff = 10**20

	# Find the min diff by comparing difference
	# of all possible pairs in given array
	for i in range(n-1):
		for j in range(i+1,n):
			if abs(arr[i]-arr[j]) < diff:
				diff = abs(arr[i] - arr[j])

	# Return min diff
	return diff

for path in paths:

	index = 1
	 
 	#Get every file in the path of the folder and store in a list
	files = [os.path.join(path, f) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
	#files = files[84:]

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

	print("=========================================================")
	print(model)

	fig = plt.figure(figsize=(8,6))	

	object_checker_array = []
	object_full_array    = []

	for file in files:

		with fits.open(file, memmap=True) as hdul:

			age_value = hdul[1].header["age_massW"]
			if not log_true:
				age_value = 10**age_value

			head, tail = os.path.split(file)

			object_full_array.append((tail, age_value))

	for age_bin_index in range(len(age_bins)):

		age_array           = []
		csp_light_age_array = []
		csp_mass_age_array  = []

		counter = 0

		for file in files:

			with fits.open(file, memmap=True) as hdul:

				age_value = hdul[1].header["age_massW"]

				if not log_true:
					age_value = 10**age_value

				if age_bins[age_bin_index][0] <= age_value < age_bins[age_bin_index][1]: 

					#print(age_bins[age_bin_index][0], "<=", age_value, "<", age_bins[age_bin_index][1])

					csp_age   = np.ndarray(hdul[1].header['ssp_number'])
					csp_light = np.ndarray(hdul[1].header['ssp_number'])
					csp_mass  = np.ndarray(hdul[1].header['ssp_number'])

					for i in range(len(csp_age)):
						csp_age[i]   = hdul[1].header['log_age_ssp_' + str(i)]
						csp_light[i] = hdul[1].header['weightLight_ssp_' + str(i)]
						csp_mass[i]  = hdul[1].header['weightMass_ssp_' + str(i)]
					
					for i, a in enumerate(csp_age):
						if not log_true:
							a = 10**a
						age_array.append(a)


					for i in csp_light:
						csp_light_age_array.append(i)

					for i in csp_mass:
						csp_mass_age_array.append(i)

					head, tail = os.path.split(file)

					object_checker_array.append((tail, age_value))

					counter = counter + 1

		#fig.suptitle("Model " + model + " with associated distributions of Usher et al data")

		print(index)

		ax1 = fig.add_subplot(rows, columns, index)

		age_array_, csp_light_age_array_ = stack_bars(age_array, csp_light_age_array)
		try:
			#width = abs(age_array[0] - age_array[-1])/(len(age_array)*2)
			width = findMinDiff(age_array_)
		except:
			width = 1

		if width >= 1.5:
			width = 1.5
		ax1.bar(age_array_, csp_light_age_array_, width = width, align='center', alpha = alpha, color = color, zorder = 0, label = 'Age LW components: ' + model, edgecolor=color)

	
		age_array_, csp_mass_age_array_ = stack_bars(age_array, csp_mass_age_array)
		try:
			#width = abs(age_array[0] - age_array[-1])/(len(age_array)*2)
			width = findMinDiff(age_array_)
		except:
			width = 1
		ax1.bar(age_array_, csp_mass_age_array_, width = width, align='center', alpha = alpha -0.2, color = "black", zorder = 0, label = 'Age MW components: ' + model, edgecolor="black")
		

		#ax1.set_ylim(bottom=0)
		ax1.set_ylim(0, 1)
		ax1.set_aspect(1 / ax1.get_data_ratio())

		ax1.set_xlabel("Age " + ('log(Gyr)' if log_true else 'Gyr'))
		if index == 1:
			ax1.set_ylabel('Weight')
			ax1.set_ylim(0, 1)
		else:
			ax1.set_yticklabels([])
		
		title = str(age_bins[index -1][0]) + r'$\leq Age$'
		if index != len(age_bins):
			title = title + " < " + str(age_bins[index -1][1]) + (' log(Gyr)' if log_true else ' Gyr')
		else:
			title = title + (' log(Gyr)' if log_true else ' Gyr')
		title = title + " (" + str(counter) + " spectra)"
		
		if len(age_bins) > 1:
			ax1.set_title(title)

		try:
			ax1.legend()
		except(IndexError):
			pass

		index = index + 1

	if rows > 1:
		print("starting metal")
		for metal_bin_index in range(len(metal_bins)):

			metal_array           = []
			csp_light_metal_array = []
			csp_mass_metal_array  = []

			counter = 0

			for file in files:

				with fits.open(file, memmap=True) as hdul:

					metal_mass_value  = hdul[1].header["metallicity_massW"]
					metal_light_value = hdul[1].header["metallicity_lightW"]

					#if not log_true:
					#	metal_value = 10**metal_value

					if metal_bins[metal_bin_index][0] <= metal_mass_value < metal_bins[metal_bin_index][1]: 

						csp_Z     = np.ndarray(hdul[1].header['ssp_number'])
						csp_light = np.ndarray(hdul[1].header['ssp_number'])
						csp_mass  = np.ndarray(hdul[1].header['ssp_number'])

						for i in range(len(csp_Z)):
							csp_Z[i]     = hdul[1].header['metal_ssp_' + str(i)]
							csp_light[i] = hdul[1].header['weightLight_ssp_' + str(i)]
							csp_mass[i]  = hdul[1].header['weightMass_ssp_' + str(i)]
						
						for i, m in enumerate(csp_Z):
							metal_array.append(m)

						for i in csp_light:
							csp_light_metal_array.append(i)

						for i in csp_mass:
							csp_mass_metal_array.append(i)

						head, tail = os.path.split(file)

						counter = counter + 1

			ax1 = fig.add_subplot(rows, columns, index)

			metal_array_, csp_light_metal_array_ = stack_bars(metal_array, csp_light_metal_array)
			try:
				#width = abs(metal_array[0] - metal_array[-1])/(len(metal_array)*2)
				width = findMinDiff(metal_array_)
			except:
				width = 1
			if width >= 1.5:
				width = 0.5
				ax1.set_xlim(-2, 2)

			ax1.bar(metal_array_, csp_light_metal_array_, width=width, align='center', alpha = alpha, color = color, zorder = 0, label = 'Metal LW components: ' + model, edgecolor=color)

			#ax1.set_ylim(bottom=0)
			ax1.set_ylim(0, 1)
			ax1.set_aspect(1 / ax1.get_data_ratio())

			metal_array_, csp_mass_metal_array_ = stack_bars(metal_array, csp_mass_metal_array)
			try:
				#width = abs(metal_array[0] - metal_array[-1])/(len(metal_array)*2)
				width = findMinDiff(metal_array_)
			except:
				width = 1

			if width >= 1.5:
				width = 0.5
				ax1.set_xlim(-2, 2)

			ax1.bar(metal_array_, csp_mass_metal_array_, width=width, align='center', alpha = alpha -0.2, color = "black", zorder = 0, label = 'Metal MW components: ' + model, edgecolor="black")

			ax1.set_xlabel('[Z/H] dex')
			if index == 1 + columns:
				ax1.set_ylabel('Weight')
			else:
				ax1.set_yticklabels([])
				
			title = str(metal_bins[index - columns -1][0]) + r'$\leq$ [Z/H]'
			if index  - columns != len(metal_bins):
				title = title + " < " + str(metal_bins[index - columns -1][1]) + " dex"
			else:
				title = title + " dex"
			title = title + " (" + str(counter) + " spectra)"
			
			if len(metal_bins) > 1:
				ax1.set_title(title)
			
			try:
				ax1.legend()
			except(IndexError):
				pass

			index = index + 1

	plt.tight_layout(rect=[0, 0.03, 1, 0.95])
	
	try:
		display_plot = sys.argv[1].lower() == 'true'
	except:
		display_plot = False

	if display_plot:
		plt.show()

	fig.savefig("output/dissertation/data/look_backtime/" + model + "_lookback_time" + ("_full" if len(age_bins) == 1 else "") +  ".png")

	print("Missing values in second list:", (set(object_full_array).difference(object_checker_array)))



