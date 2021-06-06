from astropy.io import fits
import numpy as np
import os

dissertation_results = False

input_file      = [
    #os.path.join(os.getcwd(), "example_data/spec-0266-51602-0001.fits"),
    os.path.join(os.getcwd(), "usher_r2000/NGC0104_2015-01-30.fits"),
    #os.path.join(os.getcwd(), "usher_r2000/NGC0330_2015-11-02.fits"),
    #os.path.join(os.getcwd(), "usher_r2000/NGC0121_2015-08-11.fits"),                     
]

path = os.path.join(os.getcwd(), "usher_r2000")
#input_file = [os.path.join(path, f) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

ask_to_override     = False
write_results       = True

# choose model: 'm11', 'MaStar', 'CONROY'
model_key           = 'CONROY'
model_lib           = ['E-CONROY'] 
imfs                = ['kr']
#downgrade_to_conroy = False
MASTAR_VERSION      = "vMPL11"

r_instrument_value = 2000 

#Define limts
age_limits = [0,'AoU']
Z_limits   = [-3.,3.]

#Define flux variables
data_wave_medium  = 'vacuum' 

#Firefly assumes flux units of erg/s/A/cm^2. Choose factor in case flux is scaled (e.g. flux_units=10**(-17) for SDSS)
flux_units        = 10**(-17)

#Emission lines 
N_angstrom_masked = 0
emlines           = [

	'He-II',# 'He-II:  3202.15A, 4685.74'
	'Ne-V', #  'Ne-V:   3345.81, 3425.81'
    'O-II',#  'O-II:   3726.03, 3728.73'
    'Ne-III',# 'Ne-III: 3868.69, 3967.40'
    'H-ζ',#   'H-ζ:     3889.05'
    'H-ε', #  'H-ε:     3970.07'
    'H-δ',#   'H-δ:     4101.73'
    'H-γ',#   'H-γ:     4340.46'
    'O-III',# 'O-III:  4363.15, 4958.83, 5006.77'
    'Ar-IV',# 'Ar-IV:  4711.30, 4740.10'
    'H-β',#   'H-β:     4861.32'
    'N-I',#   'H-I:    5197.90, 5200.39'
    'He-I',#  'He-I:   5875.60'
    'O-I',#   'O-I:    6300.20, 6363.67'
    'N-II',#  'N-II:   6547.96, 6583.34'
    'H-α',#   'H-α:     6562.80'
    'S-II',#  'S-II:   6716.31, 6730.68'
    'Ar-III',#'Ar-III: 7135.67'

]

milky_way_reddening = True # set whether to correct for Milky Way reddening
hpf_mode            = 'on' # set parameters for dust determination: 'on', 'hpf_only' (i.e. E(B-V)=0)
dust_law            = 'calzetti' # 'calzetti', 'allen', 'prevot' 

# Only change the following parameters, if you know what you are doing.
max_ebv               = 1.5                   
num_dust_vals         = 200             
dust_smoothing_length = 200 
max_iterations        = 10
pdf_sampling          = 300  

############################################################################################################


if model_key.upper() == "CONROY":
    downgrade_to_conroy = True
else:
    downgrade_to_conroy = False

directory = model_key

for mod_lib in model_lib:

    if mod_lib == "E-MaStar":
        mod_lib = "E"
    elif mod_lib == "Th-MaStar":
        mod_lib = "Th"

    elif mod_lib == "E-CONROY":
        mod_lib = "E"
    elif mod_lib == "Th-CONROY":
        mod_lib = "Th"

if model_key == "MaStar":
    directory = directory + "_" + mod_lib + "_" + MASTAR_VERSION

for imf in imfs:

    directory = directory + "_" + imf

output_folder = os.path.join(os.getcwd(), "output")

if dissertation_results:
    output_folder = os.path.join(output_folder, "dissertation")
    
output_folder = os.path.join(output_folder, directory.upper())


if downgrade_to_conroy:
    output_folder = os.path.join(output_folder, "downgraded")

############################################################################################################

def read_usher(input_file):

    hdul = fits.open(input_file)    

    cdelt1 = hdul[2].header['CDELT1'] #change in wavelength between pixels
    crval1 = hdul[2].header['CRVAL1'] #Starting wavelength  
    crval2 = hdul[2].header['CRVAL2'] #End wavelength

    #create wavelength array
    wavelength = np.arange(crval1, crval2, cdelt1)  
    flux       = hdul[2].data   
    error      = hdul[1].data

    #Reformat it to correct shape to make it compatible  
    wavelength = wavelength.reshape((-1,))
    flux       = flux.reshape((-1,))
    error      = flux.reshape((-1,))

    error = flux - error

    redshift = hdul[0].header['redshift']
    ra       = hdul[0].header['ra']
    dec      = hdul[0].header['dec']
    vdisp    = hdul[0].header['veldisp']

    res = hdul[2].header['RES']
    name = hdul[2].header['OBJECT']

    hdul.close()

    restframe_wavelength = wavelength/(1+redshift)

    #Correct for the redshift if value is zero.
    #Causes issues when calculating mass.
    if redshift == 0:
        hdul = fits.open(os.path.join(os.getcwd(), "usher_r2000", "extra", "usher-gcs-z_v2-LH.fits"))
        object_names = hdul[1].data[hdul[1].columns[0].name].tolist()
        redshift_values = hdul[1].data[hdul[1].columns[1].name].tolist()

        index = object_names.index(name)
        redshift = redshift_values[index]

        restframe_wavelength = wavelength/(1+redshift)

        hdul.close()

    return wavelength, restframe_wavelength, flux, error, redshift, ra, dec, vdisp

def read_sdss(input_file):

    hdul = fits.open(input_file)   

    wavelength = 10**hdul[1].data['loglam']
    flux       = hdul[1].data['flux']
    error      = hdul[1].data['ivar']**(-0.5)

    ra       = hdul[0].header['RA']  
    dec      = hdul[0].header['DEC']
    redshift = hdul[2].data['Z'][0]
    vdisp    = hdul[2].data['VDISP'][0]

    hdul.close()

    restframe_wavelength = wavelength/(1+redshift)

    return wavelength, restframe_wavelength, flux, error, redshift, ra, dec, vdisp

"""
if __name__ == "__main__":

    input_file = os.path.join(os.getcwd(), "usher_r2000/NGC0104_2015-01-30.fits")
    hdul = fits.open(input_file)

    hdul.info()

    print(hdul[0].header)
    print(hdul[1].header)
    print(hdul[2].header)
    print(hdul[0].data)
    print(hdul[1].data)

"""