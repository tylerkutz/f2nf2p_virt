
import numpy as np
import matplotlib.pyplot as plt
import sys
if len(sys.argv) != 3:
	print "Incorrect number of arguments. Please use: "
	print "\tpython code.py [Convolution-SF] [Convolution-LC]"
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



x_SF = []; 		x_LC  = [];
he3_SF = [];		h3_SF = [];
he3_LC = [];		h3_LC = [];
with open(sys.argv[1],"rb") as f:
	for line in f:
		parse = line.strip().split(" ")
		x_SF.append(		float(parse[0])	)
		he3_SF.append(	float(parse[2])	)
		h3_SF.append(		float(parse[3])	)
with open(sys.argv[2],"rb") as f:
	for line in f:
		parse = line.strip().split(" ")
		x_LC.append(		float(parse[0])	)
		he3_LC.append(	float(parse[2])	)
		h3_LC.append(		float(parse[3])	)


plt.figure(1)
plt.plot(x_SF,he3_SF,color='red',label='He-3 No Modification, SF')
plt.plot(x_SF,h3_SF,color='blue',label='H-3 No Modification,SF')
plt.plot(x_LC,he3_LC,color='red',label='He-3 No Modification, LC',linestyle='--')
plt.plot(x_LC,h3_LC,color='blue',label='H-3 No Modification,LC',linestyle='--')
#plt.plot(x,he3_onlySP,color='red',label='He-3 no off',linestyle='--')
#plt.plot(x,h3_onlySP,color='blue',label='H-3 no off',linestyle='--')
plt.legend(numpoints=1,loc='best')
plt.errorbar(mar_x,mar_he3,yerr=mar_he3e,fmt='ro')
plt.errorbar(mar_x,mar_h3,yerr=mar_h3e,fmt='bo')
plt.ylim([0.8,1.4])
plt.xlabel('xB',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel('A=3 EMC Ratio',fontsize=16)
plt.grid(True)
plt.savefig("emc-effect-without-off.pdf",bbox_inches="tight")
plt.show()
