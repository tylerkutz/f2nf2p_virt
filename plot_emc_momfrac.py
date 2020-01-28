import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2:
	print "Incorrect number of arguments. Please use: "
	print "\tpython code.py [Integral-He3+H3] "
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



thisFi_x = []; thisFi_he3 = []; thisFi_h3 =[];
thisFi_he3_2 = [];  thisFi_h3_2 = [];
thisFi_he3_4 = [];  thisFi_h3_4 = [];
thisFi_he3_6 = [];  thisFi_h3_6 = [];
thisFi_he3_8 = [];  thisFi_h3_8 = [];
with open(sys.argv[1],"rb") as f:
	for line in f:
		parse = line.strip().split(" ")
		x = float(parse[0])
		he3_0 = float(parse[2])
		h3_0  = float(parse[3])
		he3_2 = float(parse[4])
		h3_2  = float(parse[5])
		he3_4 = float(parse[6])
		h3_4  = float(parse[7])
		he3_6 = float(parse[8])
		h3_6  = float(parse[9])
		he3_8 = float(parse[10])
		h3_8  = float(parse[11])

		thisFi_x.append(x)
		thisFi_he3.append(he3_0)
		thisFi_h3.append(h3_0)
		thisFi_he3_2.append(he3_2)
		thisFi_he3_4.append(he3_4)
		thisFi_he3_6.append(he3_6)
		thisFi_he3_8.append(he3_8)
		thisFi_h3_2.append(h3_2)
		thisFi_h3_4.append(h3_4)
		thisFi_h3_6.append(h3_6)
		thisFi_h3_8.append(h3_8)


plt.figure(1)
plt.plot(thisFi_x,thisFi_he3,color='red',label='Mom > 0')
plt.plot(thisFi_x,thisFi_h3,color='blue',label='Mom > 0')

plt.legend(numpoints=1,loc='best')
plt.errorbar(mar_x,mar_he3,yerr=mar_he3e,fmt='ro')
plt.errorbar(mar_x,mar_h3,yerr=mar_h3e,fmt='bo')
plt.ylim([0.8,1.4])
plt.grid(True)
plt.savefig('emc-ratios.pdf',bbox_inches="tight")

plt.figure(2)
thisFi_he3=np.asarray(thisFi_he3)
thisFi_he3_2=np.asarray(thisFi_he3_2)
thisFi_he3_4=np.asarray(thisFi_he3_4)
thisFi_he3_6=np.asarray(thisFi_he3_6)
thisFi_he3_8=np.asarray(thisFi_he3_8)
thisFi_h3_2=np.asarray(thisFi_h3_2)
thisFi_h3_4=np.asarray(thisFi_h3_4)
thisFi_h3_6=np.asarray(thisFi_h3_6)
thisFi_h3_8=np.asarray(thisFi_h3_8)
plt.plot(thisFi_x,thisFi_he3_2/thisFi_he3,linestyle='solid',label='Mom > 50',linewidth=3,color='red')
plt.plot(thisFi_x,thisFi_he3_4/thisFi_he3,linestyle='dotted',label='Mom > 75',linewidth=3,color='red')
plt.plot(thisFi_x,thisFi_he3_6/thisFi_he3,linestyle='dashed',label='Mom > 100',linewidth=3,color='red')
plt.plot(thisFi_x,thisFi_he3_8/thisFi_he3,linestyle='dashdot',label='Mom > 200',linewidth=3,color='red')
plt.plot(thisFi_x,thisFi_h3_2/thisFi_h3,linestyle='solid',linewidth=3,color='blue')
plt.plot(thisFi_x,thisFi_h3_4/thisFi_h3,linestyle='dotted',linewidth=3,color='blue')
plt.plot(thisFi_x,thisFi_h3_6/thisFi_h3,linestyle='dashed',linewidth=3,color='blue')
plt.plot(thisFi_x,thisFi_h3_8/thisFi_h3,linestyle='dashdot',linewidth=3,color='blue')
plt.legend(numpoints=1,loc=3)
plt.grid(True)
plt.xlim([0.1,0.9])
plt.ylim([-0.2,1])
plt.savefig('emc-momfractions.pdf',bbox_inches="tight")

plt.show()
