import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
	print "Incorrect number of arguments. Please use: "
	print "\tpython code.py [Min-He3+H3] [Min-H3] [Min-He3]"
	exit(-1)

PARs = []
for fi in sys.argv:
	if '.py' in fi: continue

	np_a = 0;
	np_b = 0;
	np_c = 0;
	np_d = 0;
	of_a = 0;
	with open(fi,"rb") as f:
		for line in f:
			if 'np_a' in line and "+/-" in line: np_a = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'np_b' in line and "+/-" in line: np_b = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'np_c' in line and "+/-" in line: np_c = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'np_d' in line and "+/-" in line: np_d = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'of_a' in line and "+/-" in line: of_a = float(line.strip().split("=")[-1].split("+/-")[0])
	PARs.append( [ np_a, np_b, np_c, np_d, of_a ] )

def np_phen(x,pars):
	return pars[0] + pars[1]*x + pars[2]*np.exp( pars[3]*(np.ones(len(x)) - x) )

x_axis = np.linspace(0,1,1000)
cols = ['orange','red','blue']
labs = ['He-3 and H-3','H-3','He-3']
for i in range(len(PARs)):
	plt.plot( x_axis , np_phen(x_axis,PARs[i]) , color=cols[i] , label=labs[i] ,linewidth=3)
og_pars = [-1.21721713,0.8622478,0.82047886,0.96399233]
plt.plot(x_axis,np_phen(x_axis,og_pars),label='Init F2n/F2p',linestyle='--',color='black')
plt.legend(numpoints=1,loc='best')

plt.ylabel('F2n/F2p')
plt.xlabel('x')
plt.ylim([0.2,1])
plt.grid(True)
plt.savefig('f2nf2p-off0.pdf',bbox_inches="tight")
plt.show()
