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




x = []; 
he3_LC = [];		h3_LC = [];
he3_SF = [];		h3_SF = [];
he3_no_LC = [];		h3_no_LC = [];
he3_no_SF = [];		h3_no_SF = [];

he3_no_SF_50 = [];	h3_no_SF_50 = [];
he3_no_SF_100 = [];	h3_no_SF_100 = [];
he3_no_SF_200 = [];	h3_no_SF_200 = [];
he3_no_SF_400 = [];	h3_no_SF_400 = [];
he3_no_SF_500 = [];	h3_no_SF_500 = [];
he3_no_SF_750 = [];	h3_no_SF_750 = [];
he3_no_SF_1000 = [];	h3_no_SF_1000 = [];
he3_SF_50 = [];		h3_SF_50 = [];
he3_SF_100 = [];	h3_SF_100 = [];
he3_SF_200 = [];	h3_SF_200 = [];
he3_SF_400 = [];	h3_SF_400 = [];
he3_SF_500 = [];	h3_SF_500 = [];
he3_SF_750 = [];	h3_SF_750 = [];
he3_SF_1000 = [];	h3_SF_1000 = [];
with open(sys.argv[1],"rb") as f:
	for line in f:
		parse = line.strip().split(" ")
		x.append(		float(parse[0])	)
		
		he3_SF.append(		float(parse[2])	)
		h3_SF.append(		float(parse[3])	)
		he3_LC.append(		float(parse[4])	)
		h3_LC.append(		float(parse[5])	)
		he3_no_SF.append(		float(parse[6])	)
		h3_no_SF.append(		float(parse[7])	)
		he3_no_LC.append(		float(parse[8])	)
		h3_no_LC.append(		float(parse[9])	)

		he3_SF_50.append(		float(parse[10]) )	
		h3_SF_50.append(		float(parse[11]) )	
		he3_no_SF_50.append(		float(parse[12]) )	
		h3_no_SF_50.append(		float(parse[13]) )	

		he3_SF_100.append(		float(parse[14]) )	
		h3_SF_100.append(		float(parse[15]) )	
		he3_no_SF_100.append(		float(parse[16]) )	
		h3_no_SF_100.append(		float(parse[17]) )	

		he3_SF_200.append(		float(parse[18]) )	
		h3_SF_200.append(		float(parse[19]) )	
		he3_no_SF_200.append(		float(parse[20]) )	
		h3_no_SF_200.append(		float(parse[21]) )	

		he3_SF_400.append(		float(parse[22]) )	
		h3_SF_400.append(		float(parse[23]) )	
		he3_no_SF_400.append(		float(parse[24]) )	
		h3_no_SF_400.append(		float(parse[25]) )	

		he3_SF_500.append(		float(parse[26]) )	
		h3_SF_500.append(		float(parse[27]) )	
		he3_no_SF_500.append(		float(parse[28]) )	
		h3_no_SF_500.append(		float(parse[29]) )	

		he3_SF_750.append(		float(parse[30]) )	
		h3_SF_750.append(		float(parse[31]) )	
		he3_no_SF_750.append(		float(parse[32]) )	
		h3_no_SF_750.append(		float(parse[33]) )	

		he3_SF_1000.append(		float(parse[34]) )	
		h3_SF_1000.append(		float(parse[35]) )	
		he3_no_SF_1000.append(		float(parse[36]) )	
		h3_no_SF_1000.append(		float(parse[37]) )	

plt.figure(1)
plt.plot(x,he3_SF,color='blue',linestyle='-',label='He-3, SF, Constant Offshell')
plt.plot(x,h3_SF,color='red',linestyle='-',label='H-3, SF, Constant Offshell')
plt.plot(x,he3_LC,color='blue',linestyle='--',label='He-3, LC, Constant Offshell')
plt.plot(x,h3_LC,color='red',linestyle='--',label='H-3, LC, Constant Offshell')
plt.legend(numpoints=1,loc='best')
plt.errorbar(mar_x,mar_he3,yerr=mar_he3e,fmt='bo')
plt.errorbar(mar_x,mar_h3,yerr=mar_h3e,fmt='ro')
plt.grid(True)
plt.ylim([0.8,1.5])
plt.xlabel('xB',fontsize=16)
plt.ylabel('A=3 EMC Ratios',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

#plt.savefig('emc-ratios-linearOff-SF-LC.pdf',bbox_inches="tight")

plt.figure(2)
plt.plot(x,he3_no_SF,color='blue',linestyle='-',label='He-3, SF, No Offshell')
plt.plot(x, h3_no_SF,color='red',linestyle='-',label='H-3, SF, No Offshell')
plt.plot(x,he3_no_LC,color='blue',linestyle='--',label='He-3, LC, No Offshell')
plt.plot(x, h3_no_LC,color='red',linestyle='--',label='H-3, LC, No Offshell')
plt.legend(numpoints=1,loc='best')
plt.errorbar(mar_x,mar_he3,yerr=mar_he3e,fmt='bo')
plt.errorbar(mar_x,mar_h3,yerr=mar_h3e,fmt='ro')
plt.grid(True)
plt.ylim([0.8,1.5])
plt.xlabel('xB',fontsize=16)
plt.ylabel('A=3 EMC Ratios',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

#plt.savefig('emc-ratios-linearOff-noOff-SF-LC.pdf',bbox_inches="tight")

he3_SF = np.asarray(he3_SF)
he3_LC = np.asarray(he3_LC)
he3_no_SF = np.asarray(he3_no_SF)
he3_no_LC = np.asarray(he3_no_LC)
h3_SF = np.asarray(h3_SF)
h3_LC = np.asarray(h3_LC)
h3_no_SF = np.asarray(h3_no_SF)
h3_no_LC = np.asarray(h3_no_LC)
he3_no_SF_50  = np.asarray(he3_no_SF_50 )
he3_no_SF_100 = np.asarray(he3_no_SF_100)
he3_no_SF_200 = np.asarray(he3_no_SF_200)
he3_no_SF_400 = np.asarray(he3_no_SF_400)
he3_no_SF_500 = np.asarray(he3_no_SF_500)
he3_no_SF_750 = np.asarray(he3_no_SF_750)
he3_no_SF_1000 = np.asarray(he3_no_SF_1000)
he3_SF_50  = np.asarray(he3_SF_50 )
he3_SF_100 = np.asarray(he3_SF_100)
he3_SF_200 = np.asarray(he3_SF_200)
he3_SF_400 = np.asarray(he3_SF_400)
he3_SF_500 = np.asarray(he3_SF_500)
he3_SF_750 = np.asarray(he3_SF_750)
he3_SF_1000 = np.asarray(he3_SF_1000)
h3_SF_50  = np.asarray(h3_SF_50 )
h3_SF_100 = np.asarray(h3_SF_100)
h3_SF_200 = np.asarray(h3_SF_200)
h3_SF_400 = np.asarray(h3_SF_400)
h3_SF_500 = np.asarray(h3_SF_500)
h3_SF_750 = np.asarray(h3_SF_750)
h3_SF_1000 = np.asarray(h3_SF_1000)

#plt.figure(3)
#plt.plot(x,	he3_SF - he3_no_SF	,color='red',linestyle='-',label='He-3, SF, Offshell residual')
#plt.plot(x, 	 h3_SF -  h3_no_SF	,color='blue',linestyle='-',label='H-3, SF, Offshell residual')
#plt.legend(numpoints=1,loc='best')
#plt.grid(True)
#plt.xlabel('xB',fontsize=16)
#plt.ylabel('Full Convolution - No Offshell Convolution',fontsize=16)
#plt.xticks(fontsize=14)
#plt.yticks(fontsize=14)
#
#plt.figure(4)
#plt.plot(x,he3_SF-he3_no_SF	,label='All')
#plt.plot(x,he3_SF_50-he3_no_SF_50	,label='>50')	
#plt.plot(x,he3_SF_200-he3_no_SF_200	,label='>200')	
#plt.legend(numpoints=1,loc='best')
#plt.grid(True)
#plt.xlabel('xB',fontsize=16)
#plt.xticks(fontsize=14)
#plt.yticks(fontsize=14)

plt.figure(5)
plt.plot(x,(he3_SF_50-he3_no_SF_50)/(he3_SF-he3_no_SF)		,label='>50')
plt.plot(x,(he3_SF_100-he3_no_SF_100)/(he3_SF-he3_no_SF)	,label='>100')
plt.plot(x,(he3_SF_200-he3_no_SF_200)/(he3_SF-he3_no_SF)	,label='>200')
plt.plot(x,(he3_SF_400-he3_no_SF_400)/(he3_SF-he3_no_SF)	,label='>400')
plt.plot(x,(he3_SF_500-he3_no_SF_500)/(he3_SF-he3_no_SF)	,label='>500')
plt.plot(x,(he3_SF_750-he3_no_SF_750)/(he3_SF-he3_no_SF)	,label='>750')
plt.plot(x,(he3_SF_1000-he3_no_SF_1000)/(he3_SF-he3_no_SF)	,label='>1000')

plt.legend(numpoints=1,loc='best')
plt.grid(True)
plt.ylim([0,1])
plt.xlabel('xB',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)


#plt.savefig("emc-ratio-fraction-linearOff-SF-He3.pdf",bbox_inches="tight")

plt.show()
