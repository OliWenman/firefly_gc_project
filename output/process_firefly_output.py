from astropy.io import fits
import matplotlib.pyplot as py
import matplotlib.gridspec as gridspec
import numpy as np
import sys
import os
import pandas as pd

sys.path.append(os.path.join(os.getcwd(), "../stellar_population_models"))
os.environ["FF_DIR"] = os.getcwd()
os.environ["STELLARPOPMODELS_DIR"] = os.path.join(os.environ["FF_DIR"], "stellar_population_models")

import MaStar_SSP


def create_firefly_plots(input_files, save = False, verbose = 0, preview = True):

	py.rcParams.update({'font.size': 23})

	if type(input_files) != list:
		input_files = [input_files]

	for file in input_files:

		if verbose > 0:
			print("Opening file", os.path.basename(file))

		with fits.open(file, memmap=True) as hdul:

			file = os.path.normpath(file)
			folder_list = file.split( os.sep )
			
			if verbose > 0:
				print("Opened, now loading in data")

			data       = hdul[1].data
			wave       = data['wavelength']
			flux       = data['original_data']
			model      = data['firefly_model']
			flux_error = data['flux_error']

			csp_age   = np.ndarray(hdul[1].header['ssp_number'])
			csp_Z     = np.ndarray(hdul[1].header['ssp_number'])
			csp_light = np.ndarray(hdul[1].header['ssp_number'])
			csp_mass  = np.ndarray(hdul[1].header['ssp_number'])

			for i in range(len(csp_age)):
				csp_age[i]   = hdul[1].header['log_age_ssp_'+str(i)]
				csp_Z[i]     = hdul[1].header['metal_ssp_'+str(i)]
				csp_light[i] = hdul[1].header['weightLight_ssp_'+str(i)]
				csp_mass[i]  = hdul[1].header['weightMass_ssp_'+str(i)]

			
			name = hdul[0].header['FILE']
			name = name[6:name.find("_")]

			#################################################################################################

			fig = py.figure(figsize = (12,8)) 
			gs = gridspec.GridSpec(nrows = 20, ncols =2, figure = fig)

			ax = fig.add_subplot(gs[0 : 8, :])

			ax.plot(wave,flux, label = name, c = "m", linewidth = 4)

			if hdul[1].header["MODEL"].upper() == "E-CONROY":
				color = "red"
			elif hdul[1].header["MODEL"].upper() == "TH-MASTAR":
				 color = "royalblue"
			elif hdul[1].header["MODEL"].upper() == "E-MASTAR":
				color = "navy"

			ax.plot(wave,model, label = hdul[1].header["MODEL"], c = color, linewidth = 3)
			ax.set_ylabel("Flux \n(Normalised)") #+ r"$(erg / s / A / cm^{2} )$")
			ax.set_ylim(np.amin(flux) - np.amax(flux) * 0.25, np.amax(flux) + np.amax(flux) * 0.1 )
			ax.legend(loc = 'lower right')

			ax.set_xlim(wave[0], wave[-1])
			
			ax.tick_params(labelbottom=False)
			fig.subplots_adjust(wspace=0, hspace=0)

			title = hdul[1].header["MODEL"] + "-" + hdul[1].header["IMF"] + " fitted to " + name

			if folder_list[-2] == "downgraded":
				title = title +"\n"+ r"(downgraded to $R \approx 1273$)"

			py.title(title)

			#################################################################################################

			ax = fig.add_subplot(gs[8 : 12, :])
			ax.plot(wave, flux -model, label = hdul[1].header["MODEL"], c = "m", linewidth = 3)
			ax.set_xlabel("Wavelength (Å)")
			ax.set_xlim(wave[0], wave[-1])
			ax.set_ylim([0.5, -0.5])
			ax.yaxis.set_major_locator(py.MaxNLocator(1))
			ax.yaxis.set_label_position("right")
			ax.set_xlim(wave[0], wave[-1])
			#py.rcParams['ytick.labelsize']=8
			ax.set_ylabel("Residual")
			#ax.set_yticks([-0.25, 0, 0.25])
			#ax.set_yticklabels(labels = [-0.25 , 0, 0.25], fontsize=20)
			#py.rcParams['ytick.labelsize']=24

			yabs_max = abs(max(ax.get_ylim(), key=abs))
			ax.set_ylim(ymin=-yabs_max, ymax=yabs_max)
			ax.tick_params(labelleft=False, labelright=True)
			ax.yaxis.tick_right()
			

			##################################################################################################

			ax = fig.add_subplot(gs[15 : 20, 0])
			ax.bar(10**(csp_age), csp_light, width=1, align='center', alpha=0.5)
			ax.scatter(10**(csp_age), csp_light,  color = "blue", linewidth = 5)
			ax.set_xlabel('Lookback time (Gyr)')
			ax.set_ylabel('Frequency')
			ax.set_ylim([0, 1])
			ax.set_xlim([0, 15])
			ax.xaxis.set_major_locator(py.MaxNLocator(4))
			ax.yaxis.set_major_locator(py.MaxNLocator(2))
			ax.grid()


			##################################################################################################

			ax = fig.add_subplot(gs[15 : 20, 1])
			ax.bar(csp_Z,csp_light, width=0.5, align='center', alpha=0.5, linewidth = 5)
			ax.scatter(csp_Z,csp_light, linewidth = 5, color = "blue")
			ax.set_xlabel('[Z/H] (dex)')
			ax.tick_params(labelleft = False, labelright = False)
			ax.set_xlim([-3, 3])
			ax.set_ylim([0, 1])
			ax.xaxis.set_major_locator(py.MaxNLocator(4))
			ax.yaxis.set_major_locator(py.MaxNLocator(2))
			ax.grid()

			age         = r'$Age = $'+ str(np.around(10**hdul[1].header['age_lightW'], decimals=2))  + ' Gyr'
			metallicity = r'$[Z/H] = $'+ str(np.around(hdul[1].header['metallicity_lightW'], decimals=2)) + ' dex'
			mass        = ""#r'$log (M/M_{sun}) = $'+ str(np.around(hdul[1].header['stellar_mass'], decimals=2)) 
			reddening   = r'$E(B-V) = $'+ str(np.around(hdul[1].header['EBV'], decimals=2)) + ' mag'
			
			text = age + "\n" + metallicity +  "\n" + reddening #+ "\n" + mass 

			box_prop = dict(boxstyle='round', facecolor='wheat', alpha=0.75)
		
			fig.text(0.35, 0.60, text, fontsize=22, bbox = box_prop)

			if preview:
				py.show()

			if save:
				path, file = os.path.split(file)

				file = file[:-5] + ".png"
				path = os.path.join(path, "analysis")

				if os.path.isdir(path)==False:
					os.makedirs(path)

				file = os.path.join(path, file)

				fig.savefig(file)


		if verbose > 0:
			print("Closed file")

def plot_models():

	py.rcParams.update({'font.size': 23})

	model_path = os.path.join(os.environ['STELLARPOPMODELS_DIR'],'SPP_CONROY', "E")

	tin      = [11.0, 11.0]
	Zin      = [-1.0, 0.2]
	#tin = [11.0]
	#Zin = [0.2]

	con_files = []

	for i in range(len(tin)):
		con_files.append("VCJ_v8_mcut0.08_t" + str(tin[i]) + "_Z" + (("p" + str(Zin[i])) if Zin[i] >= 0 else ("m" + str(Zin[i]*-1))) + ".ssp.imf_varydoublex.s100")


	#print(con_files)
	#sys.exit()

	con_wave    = []
	con_flux    = []

	mastar_wave = []
	mastar_flux = []

	fig = py.figure(figsize=(12,8))
	gs  = gridspec.GridSpec(nrows = 2, ncols = 1, figure = fig)

	for index in range(len(con_files)):

		usecols = [0, 74]
		names   = ["wavelength", "flux1"] 

		model_table = pd.read_table(os.path.join(con_files[index]), 
												usecols         = usecols,
												names           = names,
												header          = None, 
												delim_whitespace= True)

		con_wave.append(model_table["wavelength"])
		con_flux.append(model_table["flux1"]) #* 3.826*10**29)

		sin      = 1.3 
		ver      = "v0.2"
		lib      = "E"

		#MaStar
		t        = MaStar_SSP.t(ver) #age array
		Z        = MaStar_SSP.Z(ver) #metallicty array
		s        = MaStar_SSP.s(ver) #IMF slope array
		wave     = MaStar_SSP.wave(ver) #wavelength array
		res      = MaStar_SSP.res(ver) #resolution array
		fluxgrid = MaStar_SSP.flux(ver,lib) #flux matrix

		flux = MaStar_SSP.inter(t,Z,s,fluxgrid,tin[index],Zin[index],sin)
		
		if index == 0:
			flux = flux / np.max(flux*2.71)
		else:

			flux = flux /np.max(flux*5.3)

		mastar_wave.append(wave)
		mastar_flux.append(flux)
		"""
		ax  = fig.add_subplot(gs[index, 0])
		#ax.set_title("A comparison of Conroy and MaStar models in wavelength coverage\n Models are of age = " + str(tin[index]) + " Gy, [Fe/H] = " + str(Zin[index]) + " dex, IMF = " + str(sin))
		ax.plot(con_wave[index], con_flux[index], label = "E-CONROY", linewidth=2)
		ax.set_xlabel("Wavelength (Å)")
		ax.set_ylabel("Flux")
		ax.set_xlim(con_wave[0][0]-10, 10300)
		ax.set_yticks([])
		ax.legend()
		"""

		ax  = fig.add_subplot(gs[index, 0])
		ax.plot(mastar_wave[index], mastar_flux[index], label = "E-MaStar", linewidth=4, color = "navy")

		lib      = "Th"
		fluxgrid = MaStar_SSP.flux(ver,lib)
		flux     = MaStar_SSP.inter(t,Z,s,fluxgrid,tin[index],Zin[index],sin)

		if index == 0:
			flux = flux / np.max(flux*2.71)
		else:

			flux = flux /np.max(flux*5.3)

		ax.plot(wave, flux, label = "Th-MaStar", linewidth=3.5, color = "royalblue")

		#ax.set_title("A comparison of Conroy and MaStar models in wavelength coverage\n Models are of age = " + str(tin[index]) + " Gy, [Fe/H] = " + str(Zin[index]) + " dex, IMF = " + str(sin))
		ax.plot(con_wave[index], con_flux[index], label = "E-Conroy", linewidth=3, color = "red")
		ax.set_xlabel("Wavelength (Å)")
		ax.set_ylabel("Flux \n(Arbitrary Units)")
		ax.set_xlim(con_wave[0][0]-30, 10300)
		ax.set_yticks([])

		if index == 0:
			ax.legend(ncol=3)
		ax.xaxis.set_major_locator(py.MaxNLocator(4))


	fig.tight_layout()
	fig.savefig("Model_wavelength_comparison.png")
	py.show()

if __name__ == "__main__":


	input_files = [
		#"C:/Users/User/Documents/University/Fourth_Year/Project/firefly_gc_project/output/CONROY_E_KR/downgraded/spFly-NGC0330_2015-11-02.fits",
		#"C:/Users/User/Documents/University/Fourth_Year/Project/firefly_gc_project/output/MASTAR_E_KR/downgraded/spFly-NGC0330_2015-11-02.fits",
		"C:/Users/User/Documents/University/Fourth_Year/Project/firefly_gc_project/output/CONROY_E_KR/downgraded/spFly-NGC0121_2015-08-11.fits",
		"C:/Users/User/Documents/University/Fourth_Year/Project/firefly_gc_project/output/MASTAR_E_KR/downgraded/spFly-NGC0121_2015-08-11.fits",
		"C:/Users/User/Documents/University/Fourth_Year/Project/firefly_gc_project/output/MASTAR_Th_KR/downgraded/spFly-NGC0121_2015-08-11.fits",
	]

	#plot_models()

	
	create_firefly_plots(input_files, 
						 verbose = 0, 
						 save = True,
						 preview = True)
	
