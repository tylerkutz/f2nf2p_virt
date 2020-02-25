import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
	print "Incorrect number of arguments. Please use: "
	print "\tpython code.py [Integral-He3+H3] "
	exit(-1)

mar_x = []; mar_he3 = []; mar_he3e = []; mar_h3 = []; mar_h3e = [];
with open("../exp_data/marathon_prelim_emc.txt","rb") as f:
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
x = []
he3_nom	= []
h3_nom 	= []
he3_iso	= []
h3_iso 	= []
he3_jnk	= []
h3_jnk 	= []
with open(sys.argv[1],"rb") as f:
	for line in f:
		parse = line.strip().split(" ")
		x.append(		float(parse[0])	)
		
		he3_nom.append(		float(parse[2]) )
		h3_nom .append(		float(parse[3]) )
		he3_iso.append(		float(parse[4]) )
		h3_iso .append(		float(parse[5]) )
		he3_jnk.append(		float(parse[6]) )
		h3_jnk .append(		float(parse[7]) )

plt.figure(1)
plt.title('EMC Ratios, SF Isospin-dep with constant-x offshell',fontsize=16)
plt.plot(x,he3_nom	,color='blue',linestyle='-',label='He-3, Min',linewidth=2)
plt.plot(x,h3_nom	,color='red',linestyle='-',label='H-3, Min',linewidth=2)

plt.plot(x,he3_iso	,color='blue',linestyle='--',label='He-3, Isoscalar',linewidth=3)
plt.plot(x,h3_iso	,color='red',linestyle='--',label='H-3, Isoscalar',linewidth=3)

#plt.plot(x,he3_jnk	,color='blue',linestyle=':',label='He-3, Junk',linewidth=3)
#plt.plot(x,h3_jnk	,color='red',linestyle=':',label='H-3, Junk',linewidth=3)

plt.legend(numpoints=1,loc='best')
plt.errorbar(mar_x,mar_he3,yerr=mar_he3e,fmt='bo')
plt.errorbar(mar_x,mar_h3,yerr=mar_h3e,fmt='ro')
plt.grid(True)
plt.ylim([0.7,1.5])
plt.xlabel('xB',fontsize=16)
plt.ylabel('A=3 EMC Ratios',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()

plt.savefig('emcratios-3point.pdf')
plt.show()
