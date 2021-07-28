import os
import warnings

warnings.filterwarnings("ignore")

with open('protein_names.txt', "r") as f:
	for name in f.readlines():
		try:
			print(name[:-1])
			os.system('python one_protein.py ' + name[:-1] + '>nohup 2>errors') 
		except BaseException:
			continue

