from astropy.io import fits
import matplotlib.pyplot as py
import matplotlib.gridspec as gridspec
import numpy as np
import sys
import os
import pandas as pd

def get_lit_values(object_name):

	files = [
		#"DeAngeli_GB.txt",
		#"DeAngeli_HST.txt",
		"UsherGC.txt"
	]

	lit_tables = []
	object_lit_values =[]

	for file in files:

		lit_table = pd.read_table(os.path.join(os.getcwd(), "output", "literature_values", file), 
									#usecols         = usecols,
									#names           = names,
									#header          = None, 
									delim_whitespace= True)
		lit_tables.append(lit_table)

		object_lit_values.append(lit_table.loc[lit_table['ID'] == object_name])

	return object_lit_values

if __name__ == "__main__":

	name = "NGC0330"

	lit_values = get_lit_values(name)

	print(lit_values)

	for value in lit_values:

		print(name, "Age =", float(value['Age']), "[Fe/H] =", float(value['[Fe/H]']), "Mass =", float(value['GCMass'])) 





