import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 3:
	print "Incorrect number of arguments. Please use: "
	print "\tpython code.py [dim of cov] [Min-He3+H3] [Min-H3] [Min-He3]"
	exit(-1)

COVs=[]
for fi in sys.argv[2:]:
	if '.py' in fi: continue
	if '.out' in fi: continue
	cov = []
	ctr1 = 0
	ctr2 = 0
	with open(fi,"rb") as f:
		arr = []
		for i in range(2): next(f)
		for line in f:
			if '*' in line: continue
			val = float(line.strip().split(" ")[0])
			arr.append(val)
			ctr2+=1
			if ctr2 == int(sys.argv[1]):
				ctr2 = 0
				ctr1 += 1
				cov.append(arr)
				arr = []

	COVs.append(cov)

PARs = []
for fi in sys.argv[2:]:
	if '.py' in fi: continue
	if '.txt' in fi: continue

	np_a = 0;
	np_b = 0;
	np_c = 0;
	with open(fi,"rb") as f:
		for line in f:
			if 'np_a' in line and "+/-" in line: np_a = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'np_b' in line and "+/-" in line: np_b = float(line.strip().split("=")[-1].split("+/-")[0])
			if 'np_c' in line and "+/-" in line: np_c = float(line.strip().split("=")[-1].split("+/-")[0])
	PARs.append( [ np_a, np_b, np_c ] )

def np_phen(x,pars):
	return pars[0] * np.power( np.ones(len(x)) - x , pars[1] ) + pars[2]
	#return pars[0] + pars[1]*x + pars[2]*np.exp( pars[3]*(np.ones(len(x)) - x) )
def np_err(x,pars,cov):
	a = pars[0] ; b = pars[1] ; c = pars[2] ;

	fa = np.power( (np.ones(len(x))-x) , b )						
	fb = a * np.power( (np.ones(len(x))-x) , b ) * np.log(np.ones(len(x))-x)		
	fc = 1.

	sig_a2 = cov[0][0];
	sig_b2 = cov[1][1];
	sig_c2 = cov[2][2];
	
	sig_ab = cov[0][1];
	sig_ba = cov[1][0];
	sig_ac = cov[0][2];
	sig_ca = cov[2][0];
	sig_bc = cov[1][2];
	sig_cb = cov[2][1];

	return np.sqrt( fa*fa*sig_a2 + fb*fb*sig_b2 + fc*fc*sig_c2 + \
			fa*fb*(sig_ab + sig_ba) + fa*fc*(sig_ac + sig_ca) + fb*fc*(sig_bc + sig_cb) )



x_axis = np.linspace(0.01,0.99,1000)
cols = ['orange','red','blue']
labs = ['He-3 and H-3','H-3','He-3']
for i in range(len(PARs)):
	plt.fill_between( x_axis, np_phen(x_axis,PARs[i]) - np_err(x_axis,PARs[i],COVs[i]) ,\
				np_phen(x_axis,PARs[i]) + np_err(x_axis,PARs[i],COVs[i]) , color=cols[i] , alpha = 0.15 )
for i in range(len(PARs)):
	#plt.plot( x_axis , np_phen(x_axis,PARs[i]) , color=cols[i] , label=labs[i] ,linewidth=3)
	plt.plot( x_axis , np_phen(x_axis,PARs[i]) - np_err(x_axis,PARs[i],COVs[i]) , color=cols[i] , linewidth=3 )
	plt.plot( x_axis , np_phen(x_axis,PARs[i]) + np_err(x_axis,PARs[i],COVs[i]) , color=cols[i] , linewidth=3 )
#og_pars = [-1.21721713,0.8622478,0.82047886,0.96399233]
og_pars = [0.46500111,2.68513156,0.47170629]
plt.plot(x_axis,np_phen(x_axis,og_pars),label='Init F2n/F2p',linestyle='--',color='black')
plt.legend(numpoints=1,loc='best')

plt.ylabel('F2n/F2p',fontsize=16)
plt.xlabel('x',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylim([0.2,1])
plt.grid(True)
plt.savefig('f2nf2p-constOff-LC.pdf',bbox_inches="tight")




plt.show()
