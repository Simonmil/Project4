import numpy as np
import matplotlib.pyplot as plt
import os as os
from matplotlib import rc

files = os.listdir(".")
txtfiles = []
for names in files:
	if names.endswith(".txt"):
		txtfiles.append(names)


for file in txtfiles:
	name = file

	Temp = []
	Values = []

	file = open(file,'r')
	data = file.read()
	data = np.array(data.split()).astype(np.float)

	for i in range(len(data)):
		if(i%2==0):
			Temp.append(data[i])
		else:
			Values.append(data[i])

	file.closed
	plt.figure(txtfiles.index(name) + 1)
	plt.plot(Temp,Values,'-r')
	plt.xlabel("Temperature, kT/J")

	if name.startswith("Exp"):
		plt.ylabel("Expectationvalue of Energy, <E>")
	elif name.startswith("Abs"):
		plt.ylabel("Absolutevalue magnetic moment, <M>")
	elif name.startswith("Heat"):
		plt.ylabel("Heat capacity, Cv")
	elif name.startswith("Suc"):
		plt.ylabel("Susceptibility, X")
	
	plt.title(name[:-4])

plt.show()
