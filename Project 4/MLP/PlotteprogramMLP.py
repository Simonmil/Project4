import numpy as np
import matplotlib.pyplot as plt
import os as os


files = os.listdir(".")
txtfiles = []
for names in files:
	if names.endswith(".txt"):
		txtfiles.append(names)




for file in txtfiles:
	name = file
	name[1;3]+name[]
	file = open(file, 'r')
	data = file.read()
	data = np.array(data.split()).astype(np.float)
	
	plt.figure(txtfiles.index(name)+1)
	
	if name.startswith('Exp'):
		if (name[-7] == str(4) and name[-8] != str(1)):
			tempExp40 = []
			ExpValue40 = []
			for i in range(len(data)):
				if (i % 2 == 0):
					tempExp40.append(data[i])
				else:
					ExpValue40.append(data[i])	
			a += 1
			if a == 4:
				plt.plot(tempExp40,ExpValue40,'-r')
			
		elif (name[-7] == str(6)):
			tempExp60 = []
			ExpValue60 = []
			for i in range(len(data)):
				if (i % 2 == 0):
					tempExp60.append(data[i])
				else:
					ExpValue60.append(data[i])
			b += 1
			if b == 4:
				plt.plot(tempExp60,ExpValue60,'-r')

		elif (name[-7] == str(0) and name[-8] == str(1)):
			tempExp100 = []
			ExpValue100 = []
			for i in range(len(data)):
				if (i % 2 == 0):
					tempExp100.append(data[i])
				else:
					ExpValue100.append(data[i])	
			plt.plot(tempExp100,ExpValue100,'-r')

		elif (name[-7] == str(4) and name[-8] == str(1)):
			tempExp140 = []
			ExpValue140 = []
			for i in range(len(data)):
				if (i % 2 == 0):
					tempExp140.append(data[i])
				else:
					ExpValue140.append(data[i])	
			plt.plot(tempExp140,ExpValue140,'-r')

	if name.startswith('Abs'):
		if (name[-7] == str(4) and name[-8] != str(1)):
			tempAbs40 = []
			AbsValue40 = []
			for i in range(len(data)):
				if (i % 2 == 0):
					tempAbs40.append(data[i])
				else:
					AbsValue40.append(data[i])	
			plt.plot(tempAbs40,AbsValue40,'-r')

		elif (name[-7] == str(6)):
			tempAbs60 = []
			AbsValue60 = []
			for i in range(len(data)):
				if (i % 2 == 0):
					tempAbs60.append(data[i])
				else:
					AbsValue60.append(data[i])	
			plt.plot(tempAbs60,AbsValue60,'-r')

		elif (name[-7] == str(0) and name[-8] == str(1)):
			tempAbs100 = []
			AbsValue100 = []
			for i in range(len(data)):
				if (i % 2 == 0):
					tempAbs100.append(data[i])
				else:
					AbsValue100.append(data[i])	
			plt.plot(tempAbs100,AbsValue100,'-r')

		elif (name[-7] == str(4) and name[-8] == str(1)):
			tempAbs140 = []
			AbsValue140 = []
			for i in range(len(data)):
				if (i % 2 == 0):
					tempAbs140.append(data[i])
				else:
					AbsValue140.append(data[i])	
			plt.plot(tempAbs140,AbsValue140,'-r')

	if name.startswith('Heat'):
		if (name[-7] == str(4) and name[-8] != str(1)):
			tempHeat40 = []
			HeatValue40 = []
			for i in range(len(data)):
				if (i % 2 == 0):
					tempHeat40.append(data[i])
				else:
					HeatValue40.append(data[i])	
			plt.plot(tempHeat40,HeatValue40,'-r')

		elif (name[-7] == str(6)):
			tempHeat60 = []
			HeatValue60 = []
			for i in range(len(data)):
				if (i % 2 == 0):
					tempHeat60.append(data[i])
				else:
					HeatValue60.append(data[i])	
			plt.plot(tempHeat60,HeatValue60,'-r')

		elif( name[-7] == str(0) and name[-8] == str(1)):
			tempHeat100 = []
			HeatValue100 = []
			for i in range(len(data)):
				if (i % 2 == 0):
					tempHeat100.append(data[i])
				else:
					HeatValue100.append(data[i])
			plt.plot(tempHeat100,HeatValue100,'-r')	

		elif (name[-7] == str(4) and name[-8] == str(1)):
			tempHeat140 = []
			HeatValue140 = []
			for i in range(len(data)):
				if (i % 2 == 0):
					tempHeat140.append(data[i])
				else:
					HeatValue140.append(data[i])
			plt.plot(tempHeat140,HeatValue140,'-r')

	if name.startswith('Suc'):
		if (name[-7] ==str(4) and name[-8] != str(1)):
			tempSuc40 = []
			SucValue40 = []
			for i in range(len(data)):
				if (i % 2 == 0):
					tempSuc40.append(data[i])
				else:
					SucValue40.append(data[i])	
			plt.plot(tempSuc40,SucValue40,'-r')

		elif (name[-7] == str(6)):
			tempSuc60 = []
			SucValue60 = []
			for i in range(len(data)):
				if (i % 2 == 0):
					tempSuc60.append(data[i])
				else:
					SucValue60.append(data[i])	
			plt.plot(tempSuc60,SucValue60,'-r')

		elif (name[-7] == str(0) and name[-8] == str(1)):
			tempSuc100 = []
			SucValue100 = []
			for i in range(len(data)):
				if (i % 2 == 0):
					tempSuc100.append(data[i])
				else:
					SucValue100.append(data[i])	
			plt.plot(tempSuc100,SucValue100,'-r')

		elif (name[-7] == str(4) and name[-8] == str(1)):
			tempSuc140 = []
			SucValue140 = []
			for i in range(len(data)):
				if (i % 2 == 0):
					tempSuc140.append(data[i])
				else:
					SucValue140.append(data[i])	
			plt.plot(tempSuc140,SucValue140,'-r')

	
	#print tempExp40, ExpValue40
	file.closed
	plt.xlabel("Temperature, kT/J")
	if name.startswith('Exp'):
		plt.ylabel("Expectation value energy, <E>")
	elif name.startswith('Abs'):
		plt.ylabel("Expectation value magnetic moment, <M>")
	elif name.startswith('Heat'):
		plt.ylabel("Heat capacity, Cv")
	elif name.startswith('Suc'):
		plt.ylabel("Suceptibility")	
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