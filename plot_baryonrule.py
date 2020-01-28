import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
	print "Incorrect number of arguments. Please use: "
	print "\tpython code.py [sum rule output]"
	exit(-1)



Xs = []
He3s = []
H3s = []
for fi in sys.argv:
	if '.py' in fi: continue
	thisFi_x = []; thisFi_he3 = []; thisFi_h3 =[];
	with open(fi,"rb") as f:
		for line in f:
			parse = line.strip().split(" ")
			x = float(parse[0])
			he3 = float(parse[1])
			h3 = float(parse[2])

			thisFi_x.append(x)
			thisFi_he3.append(he3)
			thisFi_h3.append(h3)
	Xs.append( thisFi_x )
	He3s.append( thisFi_he3 )
	H3s.append( thisFi_h3 )

for i in range(len(Xs)):
	if i == 0 or i == 2:
		plt.plot(Xs[i],He3s[i],color='red',label='He-3',linestyle='-')
	if i == 0 or i == 1:
		plt.plot(Xs[i],H3s[i],color='blue',label='H-3',linestyle='-')
plt.legend(numpoints=1,loc='best')
plt.grid(True)
plt.axvline(x=0.825,color='black',linewidth=3)
plt.ylim([0,100])
plt.ylabel('Baryon Sum Rule [result / 3. * 100 ] (%)')
plt.xlabel('x')
plt.savefig('sumrule.pdf',bbox_inches='tight')
plt.show()
