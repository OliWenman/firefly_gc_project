"""
PLEASE READ!
To test, do the following:
1) Download the conroy models and create a folder "SSP_CONROY" inside "stellar_population_models" -> stellar_population_models/SSP_CONROY
2) Create a folder within "SSP_CONROY" "E" and "Th" -> stellar_population_models/SSP_CONROY/E
													-> stellar_population_models/SSP_CONROY/Th
3) unzip file "vcj_ssp_v8.tar.gz" inside "stellar_population_models/SSP_CONROY/E"
4) unzip file "atlas_rfn_v3.tar.gz" inside "stellar_population_models/SSP_CONROY/Th"
5) Place this script inside dir "stellar_population_models" -> "stellar_population_models/read_stellar_populations.py"
6) Run using "python stellar_population_models/read_stellar_population.py" when inside firefly base directory
"""

import sys, os

sys.path.append(os.path.join(os.getcwd(), "python"))
os.environ["FF_DIR"] = os.getcwd()
os.environ["STELLARPOPMODELS_DIR"] = os.path.join(os.environ["FF_DIR"], "stellar_population_models")

import numpy as np
from astropy.io import fits
import astropy.cosmology as co
import pandas as pd
import copy
import glob

from firefly_estimations_3d import estimation
from firefly_dust import hpf, unred, determine_attenuation, dust_calzetti_py
from firefly_instrument import downgrade
from firefly_fitter import fitter
from firefly_library import airtovac, convert_chis_to_probs, light_weights_to_mass, calculate_averages_pdf, normalise_spec, match_data_models


#Defined some needed values for Firefly 
default_value = -9999
EPS = 10.E-10
ebv_mw = 0
dict_imfs = {'cha': 'Chabrier', 'ss': 'Salpeter', 'kr': 'Kroupa'}


#Inputs
age_limits       = [0,15]
Z_limits         = [-3., 3.]
models           = 'm11'
model_used       = 'MARCS'   #Have done empirical library only for time being.
imf_used         = 'kr'
data_wave_medium = "vacuum"
downgrade_models = True


def trylog10(value):
	if (value<EPS):
		logv = default_value
	else:
		logv = np.log10(value)
	return logv

# first the m11 case
if models == 'm11':
	first_file  = True
	model_files = []
	
	model_path 		= os.path.join(os.environ['STELLARPOPMODELS_DIR'],'SSP_M11_'+model_used ,'ssp_M11_' +model_used +'.' + imf_used)

	print("model_path", model_path)

	# Constructs the metallicity array of models :
	all_metal_files = sorted(glob.glob(model_path+'*'))
	print("all_metal_files", all_metal_files)
	
	metal_files 	= []
	metal 	    = [] #[-2.25, -1.35, -0.33, 0, 0.35]
	for z in range(len(all_metal_files)):
		zchar = all_metal_files[z][len(model_path):]
		if zchar == 'z001':
			#znum = -0.3
			znum = 10**(-0.33) #0.5
		elif zchar == 'z002':
			#znum = 0.0
			znum = 10**(0) #1.0
		elif zchar == 'z004':
			#znum = 0.3
			znum = 10**(0.35) #2.0
		elif zchar == 'z0001.bhb':
			#znum = -1.301
			znum = 10**(-1.35) #10**-1.301
		elif zchar == 'z0001.rhb':
			#znum = -1.302
			znum = 10**(-1.35) #10**-1.302
		elif zchar == 'z10m4.bhb':
			#znum = -2.301
			znum = 10**(-2.25) #10**-2.301
		elif zchar == 'z10m4.rhb':
			#znum = -2.302
			znum = 10**(-2.25) #10**-2.302
		elif zchar == 'z10m4':
			#znum = -2.300
			znum = 10**(-2.25) #10**-2.300
		elif zchar == 'z-0.6':
			znum = 10**-0.6
		elif zchar == 'z-0.9':
			znum = 10**-0.9
		elif zchar == 'z-1.2':
			znum = 10**-1.2
		elif zchar == 'z-1.6':
			znum = 10**-1.6
		elif zchar == 'z-1.9':
			znum = 10**-1.9
		else:
			raise NameError('Unrecognised metallicity! Check model file names.')

		if znum>10**(Z_limits[0]) and znum<10**(Z_limits[1]):
			metal_files.append(all_metal_files[z])
			metal.append(znum)
	
	print("metal", metal)
	
	#constructs the model array
	model_flux, age_model, metal_model = [],[],[]

	for zi,z in enumerate(metal_files):
		# print "Retrieving and downgrading models for "+z
		model_table = pd.read_table(z,converters={'Age':np.float64}, header=None ,usecols=[0,2,3], names=['Age','wavelength_model','flux_model'], delim_whitespace=True)
		age_data = np.unique(model_table['Age'].values.ravel())
		#print(age_data)

		for a in age_data:
			logyrs_a = trylog10(a)+9.0
			print("hello")
			## print "age model selection:", self.age_limits[0], logyrs_a, self.age_limits[1]
			if (((10**(logyrs_a-9)) < age_limits[0]) or ((10**(logyrs_a-9)) > age_limits[1])):
				continue
			else:
				spectrum = model_table.loc[model_table.Age == a, ['wavelength_model', 'flux_model'] ].values
				wavelength_int,flux = spectrum[:,0],spectrum[:,1]

				delta_lamda = (wavelength_int[-1] - wavelength_int[0])/20000

				print(delta_lamda)

				print()

				exit()

				# converts to air wavelength
				if data_wave_medium == 'vacuum':
					wavelength = airtovac(wavelength_int)
				else:
					wavelength = wavelength_int

				# downgrades the model
				if downgrade_models:
					#mf = downgrade(wavelength,flux,deltal,self.vdisp_round, wave_instrument, r_instrument)
					"""
					Can't do the downgrade part in this scipt as doesn't have input spec for vdisp, but should work in practice.
					"""
					mf = copy.copy(flux)
				else:
					mf = copy.copy(flux)

				# Reddens the models
				if ebv_mw != 0:
					attenuations = unred(wavelength,ebv=0.0-ebv_mw)
					model_flux.append(mf*attenuations)
				else:
					model_flux.append(mf)

				age_model.append(a)
				metal_model.append(metal[zi])
				first_model = False
	
	#self.model_wavelength, self.model_flux, self.age_model, self.metal_model = wavelength, model_flux, age_model, metal_model

	#return wavelength, model_flux, age_model, metal_model

#Conroy models
elif models == "CONROY":
	
	#Join path together to find needed files.
	model_path = os.path.join(os.environ['STELLARPOPMODELS_DIR'],'SPP_CONROY', model_used) #,'VCJ_v8_mcut0.08_')

	if model_used != "E": #and model_used != "Th":
		raise NameError("Unrecognised model_lib!")

	print("model_path", model_path)

	# Constructs the metallicity array of models :
	#all_metal_files = sorted(glob.glob(model_path + '*'))
	#print("all_metal_files", all_metal_files)

	#Used a different method, because all files will be used as it is not dependent on IMF as you select the IMF spectra to use inside the file.
	all_metal_files = [f for f in os.listdir(model_path) if os.path.isfile(os.path.join(model_path, f))]	
	print("all_metal_files", all_metal_files)
	
	metal_files = []
	metal 	    = []

	"""
	Get the metalicility from the file names, similiar to m11. 

	PLEASE CHECK THE ZNUM VALUES ARE CORRECT. 
	"""
	for z in range(len(all_metal_files)):
		
		zchar = all_metal_files[z][len("VCJ_v8_mcut0.08_") + 6: len(".ssp.imf_varydoublex.s100") +2]

		print("zchar =", zchar)

		if zchar == 'Zm0.5':
			znum = 10**(-0.5) 

		elif zchar == 'Zm1.0':
			znum = 10**(-1.0)

		elif zchar == 'Zm1.5':
			znum = 10**(-1.5) 

		elif zchar == 'Zp0.0':
			znum = 10**(0) 

		elif zchar == 'Zp0.2':
			znum = 10**(0.2) 

		else:
			raise NameError('Unrecognised metallicity! Check model file names.')

		if znum>10**(Z_limits[0]) and znum<10**(Z_limits[1]):
			metal_files.append(all_metal_files[z])
			metal.append(znum)

	print("metal", metal)
	print("metal_files", metal_files)

	age = []

	model_flux, age_model, metal_model = [],[],[]


	for zi,z in enumerate(metal_files):

		#IMF IS LOCATED INSIDE THE FILE NOT THE FILENAME.
		#Select the correct spectra columns to use dependent on imf selected (0 is the wavelength, the other 256 spectra are IMF dependent )		
		if imf_used == 'kr': 

			usecols = [0, 74]
			names    = ["wavelength", "flux1"] 

		elif imf_used == 'ss':

			usecols = [0,154] 
			names    = ["wavelength", "flux1"]
		else:
			#Open to add more IMF options for selected spectra if needed
			raise NameError("Unrecognised IMF!")

		#Read in the data using the selected columns
		model_table = pd.read_table(os.path.join(model_path, z), 
									usecols         = usecols,
									names           = names,
									header          = None, 
									delim_whitespace= True)
		
		#Get the number of spectra that are used
		n_spec = len(model_table.columns) -1

		#Get the age values from the filename, same as M11
		a = float(z[len("VCJ_v8_mcut0.08_t") : len(".ssp.imf_varydoublex.s100") - 6])


		#Loop through the number of spectra 
		for s in range(n_spec):

			logyrs_a = trylog10(a)+9.0

			if (((10**(logyrs_a-9)) < age_limits[0]) or ((10**(logyrs_a-9)) > age_limits[1])):
				continue
			else:
				
				#Get the wavelength values
				wavelength_int = model_table["wavelength"].values

				delta_lamda = (wavelength_int[-1] - wavelength_int[0])/1273.88

				print(delta_lamda)

				print()

				exit()
				
				#Get the selected flux in the loop
				flux = model_table["flux" + str(s + 1)].values

				# converts to air wavelength
				if data_wave_medium == 'vacuum':
					wavelength = airtovac(wavelength_int)
				else:
					wavelength = wavelength_int

				# downgrades the model
				if downgrade_models:
					"""
					Can't do the downgrade part in this scipt as doesn't have input spec for vdisp, but should work in practice.
					"""

					mf = copy.copy(flux)#downgrade(wavelength, flux,deltal, self.vdisp_round, wave_instrument, r_instrument)
				else:
					mf = copy.copy(flux)

				# Reddens the models
				if ebv_mw != 0:

					attenuations = unred(wavelength, ebv=0.0-ebv_mw)
					model_flux.append(mf*attenuations)
					#model_flux.append(mf)
				else:
					model_flux.append(mf)

				age_model.append(a)
				metal_model.append(metal[zi])
				first_model = False


print("complete!")

"""
print("age_model", age_model, "\n")
print("metal_model", metal_model, "\n")
print("model_flux", model_flux)

print(len(age_model), "\n")
print(len(metal_model), "\n")
print(len(model_flux))
"""

