import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

path = "/Users/User/Documents/University/Fourth_Year/Project/firefly_gc_project/stellar_population_models/CONROY_SPP/E/"
age = "11.0"
metalicity = "p" + "0.0"

file = "VCJ_v8_mcut0.08_t" + age + "_Z" + metalicity + ".ssp.imf_varydoublex.s100"
input_file = path + file

#file = "VCJ_v8_mcut0.08_t01.0_Zm0.5.ssp.imf_varydoublex.s100"
#input_file = "/Users/User/Documents/University/Fourth_Year/Project/firefly_gc_project/stellar_population_models/SSP_M11_MILES/ssp_M11_MILES.krz-0.6"
#input_file = "/Users/User/Documents/University/Fourth_Year/Project/firefly_gc_project/stellar_population_models/CONROY_SPP/E/VCJ_v8_mcut0.08_t01.0_Zm0.5.ssp.imf_varydoublex.s100"

#input_file = "/Users/User/Documents/University/Fourth_Year/Project/firefly_gc_project/stellar_population_models/CONROY_SPP/E/VCJ_v8_mcut0.08_t13.5_Zp0.2.ssp.imf_varydoublex.s100"

model_table = pd.read_table(input_file, 
							delim_whitespace=True)

wavelength = model_table.iloc[:, 0]
flux       = model_table.iloc[:, -1]

#print(model_table)
print(wavelength)
print(model_table)

model_table.insert(0, "age", [age]*10565, True)

print(model_table) 

#plt.plot(wavelength[0: 6000], flux[0:6000])
plt.plot(wavelength, flux)
plt.show()