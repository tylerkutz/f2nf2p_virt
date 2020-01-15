import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
	print "Incorrect number of arguments. Please use: "
	print "\tpython code.py [Integral-He3+H3] [Integral-H3] [Intgeral-He3]"
	exit(-1)

mar_x = []; mar_he3 = []; mar_he3e = []; mar_h3 = []; mar_h3e = [];
with open("exp_data/marathon_prelim_emc.txt","rb") as f:
        for line in f:
                if '#' in line: continue
                parse = line.strip().split("\t")
                x=float(parse[0])
                he3=float(parse[1])
                he3e=float(parse[2])
                h3=float(parse[3])
                h3e=float(parse[4])
                he3h3=float(parse[5])
                he3h3e=float(parse[6])

                mar_x.append(x)
                mar_he3.append(he3)
                mar_he3e.append(he3e)
                mar_h3.append(h3)
                mar_h3e.append(h3e)



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
			he3 = float(parse[2])
			h3 = float(parse[3])

			thisFi_x.append(x)
			thisFi_he3.append(he3)
			thisFi_h3.append(h3)
	Xs.append( thisFi_x )
	He3s.append( thisFi_he3 )
	H3s.append( thisFi_h3 )

labs = ['He-3 and H-3','H-3','He-3']
stys = ['-','--',':']
for i in range(len(Xs)):
	if i == 0 or i == 2:
		plt.plot(Xs[i],He3s[i],color='red',label=labs[i],linestyle=stys[i])
	if i == 0 or i == 1:
		plt.plot(Xs[i],H3s[i],color='blue',label=labs[i],linestyle=stys[i])
plt.legend(numpoints=1,loc='best')
plt.errorbar(mar_x,mar_he3,yerr=mar_he3e,fmt='ro')
plt.errorbar(mar_x,mar_h3,yerr=mar_h3e,fmt='bo')
plt.grid(True)
plt.savefig('emc_ratios.pdf',bbox_inches="tight")
plt.show()
