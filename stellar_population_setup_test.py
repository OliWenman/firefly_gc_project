import numpy as np
from astropy.io import fits
import glob
import pandas as pd
import os,sys

sys.path.append(os.path.join(os.getcwd(), "python"))
os.environ["FF_DIR"] = os.getcwd()
os.environ["STELLARPOPMODELS_DIR"] = os.path.join(os.environ["FF_DIR"], "stellar_population_models")

import firefly_models as fm

if __name__ == "__main__":

	print(os.getcwd())
	print(os.path.join(os.environ["FF_DIR"], "stellar_population_models"))

	#prepare model templates
	model = fm.StellarPopulationModel(spec, 
									  output_file, 
									  cosmo, 
									  models = ["Conroy"], 
									  model_libs = ["Emprical"])
