"""
.. moduleauthor:: Daniel Thomas <daniel.thomas__at__port.ac.uk>
.. contributions:: Johan Comparat <johan.comparat__at__gmail.com>
.. contributions:: Violeta Gonzalez-Perez <violegp__at__gmail.com>

Firefly is initiated with this script. 
All input data and parmeters are now specified in this one file.

"""

import sys, os
import warnings
warnings.filterwarnings("ignore")

sys.path.append(os.path.join(os.getcwd(), "python"))
os.environ["FF_DIR"] = os.getcwd()
os.environ["STELLARPOPMODELS_DIR"] = os.path.join(os.environ["FF_DIR"], "stellar_population_models")

import numpy as np
import astropy.cosmology as co
import time
from astropy.io import fits
import matplotlib.pyplot as plt

from multiprocessing import Pool

import firefly_setup as fs
import firefly_models as fm
from firefly_instrument import downgrade, match_spectral_resolution

import inputs as inp

############################################################################################################################################################################
#Get the user values from the input file for Firefly to use

input_file = inp.input_file
if type(input_file) != list:
	input_file = [input_file] 

outputFolder          = inp.output_folder
r_instrument_value    = inp.r_instrument_value
model_key             = inp.model_key
model_lib             = inp.model_lib 
imfs                  = inp.imfs 
age_limits            = inp.age_limits
Z_limits              = inp.Z_limits
data_wave_medium      = inp.data_wave_medium 
flux_units            = inp.flux_units
N_angstrom_masked     = inp.N_angstrom_masked
emlines               = inp.emlines 
milky_way_reddening   = inp.milky_way_reddening
hpf_mode              = inp.hpf_mode
dust_law              = inp.dust_law 
max_ebv               = inp.max_ebv                   
num_dust_vals         = inp.num_dust_vals              
dust_smoothing_length = inp.dust_smoothing_length
max_iterations        = inp.max_iterations
pdf_sampling          = inp.pdf_sampling  
ask_to_override       = inp.ask_to_override
write_results         = inp.write_results

suffix = ""

cosmo = co.Planck15

downgrade_to_conroy = inp.downgrade_to_conroy

############################################################################################################################################################################
#Setup the input data for Firefly

check_inputs = True

#Display key inputs to the user before processing begins and ask if they want to continue.
#Helps to reduce mistakes before mass processing of files.
if check_inputs:

	print("")

	for file in input_file:
		
		if file == input_file[0]:
			print("input_file(s)       =", os.path.basename(file))
		else:
			print("                     ", os.path.basename(file))
	
	print("model_lib           =", model_lib)
	print("imfs                =", imfs)
	print("downgrade_to_conroy =", downgrade_to_conroy)
	print("data_wave_medium    =", data_wave_medium)
	print("flux_units          =", flux_units)

	answer = input("\nPlease check if correct. Continue?")
	
	if (answer=='N' or answer=='n'):
		sys.exit()

print('\nStarting firefly...')

#Create the output folders if they don't exist 
if os.path.isdir(outputFolder)==False:
	os.makedirs(outputFolder)

for file in input_file:

	#Construct the output file
	output_file = os.path.join( outputFolder , 'spFly-' + os.path.basename( file )[0:-5] ) + ".fits"

	#Ask the user if they want to override the file if setting is turned on
	if os.path.isfile(output_file) and ask_to_override:
		
		print()
		print('Warning: This object has already been processed, the file will be over-witten.')
		answer = input('** Do you want to continue? (Y/N)')
		
		if (answer=='N' or answer=='n'):
			sys.exit()
		
		os.remove(output_file)


	#Read in the data. Try reading using usher method first, if fails default back to the standard SDSS way.
	try:
		wavelength, restframe_wavelength, flux, error, redshift, ra, dec, vdisp = inp.read_usher(file)
	except:
		wavelength, restframe_wavelength, flux, error, redshift, ra, dec, vdisp = inp.read_sdss(file)

	#instrumental resolution
	r_instrument = np.full((len(wavelength), ), r_instrument_value)
	
	#If true, downgrade the input data (flux and r_instruemnt) to the same resolution as the CONROY models
	if downgrade_to_conroy:
		
		c = 3 * 10**5 #[km/h]
		conroy_sigma = 100 #[km/h] 
		r_conory = np.full( (len(wavelength), ), c/(2.35*conroy_sigma) )

		flux, r_instrument, sigma_offset, out_mask = match_spectral_resolution(wave          = wavelength,
																			   flux          = flux, 
																			   sres          = r_instrument, 
																			   new_sres_wave = wavelength, 
																			   new_sres      = r_conory)
	

	#################################################################################################################################################################
	#Core firefly, don't change

	t0 = time.time()

	age_min = age_limits[0]

	if type(age_limits[1]) == str:
		
		if age_limits[1] == 'AoU':
			age_max = cosmo.age(redshift).value
		
		elif age_limits[1] != 'AoU':
			print('Unrecognised maximum age limit. Try again.')
			sys.exit()
	else:
		age_max = age_limits[1]

	Z_min = Z_limits[0]
	Z_max = Z_limits[1]

	print()
	print( 'Output file: ', output_file                 )
	print()

	prihdr                          = fm.pyfits.Header()
	prihdr['FILE']                  = os.path.basename(output_file)
	prihdr['MODELS']	            = model_key
	prihdr['FITTER']	            = "FIREFLY"	
	prihdr['AGEMIN']	            = str(age_min)		
	prihdr['AGEMAX']	            = str(age_max)
	prihdr['ZMIN']	                = str(Z_min)
	prihdr['ZMAX']	                = str(Z_max)
	prihdr['redshift']	            = redshift
	prihdr['HIERARCH age_universe']	= np.round(cosmo.age(redshift).value,3)
	prihdu                          = fm.pyfits.PrimaryHDU(header=prihdr)
	tables                          = [prihdu]

	#define input object to pass data on to firefly modules and initiate run
	spec = fs.firefly_setup(file,
							milky_way_reddening = milky_way_reddening, 
							N_angstrom_masked   = N_angstrom_masked,
							hpf_mode            = hpf_mode)

	spec.openSingleSpectrum(wavelength, 
							flux, 
							error, 
							redshift, 
							ra, 
							dec, 
							vdisp, 
							emlines, 
							r_instrument)

	did_not_converge = 0.
	try :
		#prepare model templates
		model = fm.StellarPopulationModel(spec, 
										  output_file, 
										  cosmo, 
										  write_results         = write_results,
										  models                = model_key,
										  model_libs            = model_lib, 
										  imfs                  = imfs, 
										  age_limits            = [age_min,age_max], 
										  downgrade_models      = True, 
										  data_wave_medium      = data_wave_medium, 
										  Z_limits              = Z_limits, 
										  suffix                = suffix, 
										  use_downgraded_models = False,
										  dust_law              = dust_law, 
										  max_ebv               = max_ebv, 
										  num_dust_vals         = num_dust_vals, 
										  dust_smoothing_length = dust_smoothing_length,
										  max_iterations        = max_iterations, 
										  pdf_sampling          = pdf_sampling, 
										  flux_units            = flux_units)

		#initiate fit
		model.fit_models_to_data()
		tables.append( model.tbhdu )

	except (ValueError):

		tables.append( model.create_dummy_hdu() )
		did_not_converge +=1
		print('did not converge')

	if did_not_converge < 1 :
		complete_hdus = fm.pyfits.HDUList(tables)

		if os.path.isfile(output_file):
			os.remove(output_file)		
		
		complete_hdus.writeto(output_file)

	print()
	print ("Done... total time:", int(time.time()-t0) ,"seconds.")
	print()

############################################################################################################################################################################
