import numpy as np
import matplotlib.pyplot as plt
import os as os

print 0%2

files = os.listdir(".")
txtfiles = []
for names in files:
	if names.endswith(".txt"):
		txtfiles.append(names)


for file in txtfiles:
	name = file
	
	data = []
	y = []
	x = []
	file = open(file, 'r')
	data = file.read()
	data = np.array(data.split()).astype(np.float)
	for i in range(len(data)):
		if (i % 2 == 0):
			x.append(data[i])
		else:
			y.append(data[i])
	file.closed
	plt.figure(txtfiles.index(name)+1)
	if name.startswith('Pro'):
		x = np.linspace(-8,8,400)
		plt.bar(x,y,1./20,color="red")
		plt.xlabel("Energy, J")
		plt.ylabel("Probabilty")
	else:
		plt.plot(x,y,'-r')
		plt.xlabel("Temperature, kT/J")
		if name.startswith('Ab'):
			plt.ylabel("Magnetic moment, <M>")
		elif name.startswith('Ac'):
			plt.ylabel("Accepted Configurations")
		elif name.startswith('E'):
			plt.ylabel("Energy, <E>")
	plt.title(name[:-4])
	
plt.show()

""" 
plt.figure(1)
plt.plot(x,lines1,'r-', label='10 MonteCarlo-cycles')
plt.xlabel('Probability')
plt.ylabel('u(x)')
plt.title('10 MonteCarlo-cycles')
#plt.legend()


plt.figure(2)
plt.plot(x,lines2,'r-', label='100 MonteCarlo-cycles')
plt.xlabel('Probability')
plt.ylabel('u(x)')
plt.title('100 MonteCarlo-cycles')
#plt.legend()

plt.show()
"""