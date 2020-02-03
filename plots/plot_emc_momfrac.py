import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2:
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
he3_full = [];		h3_full = [];
he3_onlySP = [];	h3_onlySP = [];
he3_full_1	= [];	he3_full_5	= [];
h3_full_1	= [];   h3_full_5	= [];
he3_onlySP_1	= [];   he3_onlySP_5	= [];
h3_onlySP_1	= [];   h3_onlySP_5	= [];
he3_full_2	= [];   he3_full_6	= [];
h3_full_2	= [];   h3_full_6	= [];
he3_onlySP_2	= [];   he3_onlySP_6	= [];
h3_onlySP_2	= [];   h3_onlySP_6	= [];
he3_full_3	= [];   he3_full_7	= [];
h3_full_3	= [];   h3_full_7	= [];
he3_onlySP_3	= [];   he3_onlySP_7	= [];
h3_onlySP_3	= [];   h3_onlySP_7	= [];
he3_full_4	= [];
h3_full_4	= [];
he3_onlySP_4	= [];
h3_onlySP_4	= [];
with open(sys.argv[1],"rb") as f:
	for line in f:
		parse = line.strip().split(" ")
		x.append(		float(parse[0])	)
		he3_full.append(	float(parse[2])	)
		h3_full.append(		float(parse[3])	)
		#he3_onlySP.append(	float(parse[4])	)
		#h3_onlySP.append(	float(parse[5])	)

		#he3_full_1.	append(	float(parse[6])	)
		#h3_full_1.	append(	float(parse[7])	)
		#he3_onlySP_1.	append(	float(parse[8])	)
		#h3_onlySP_1.	append(	float(parse[9])	)

		#he3_full_2.	append(	float(parse[10])	)
		#h3_full_2.	append(	float(parse[11])	)
		#he3_onlySP_2.	append(	float(parse[12])	)
		#h3_onlySP_2.	append(	float(parse[13])	)

		#he3_full_3.	append(	float(parse[14])	)
		#h3_full_3.	append(	float(parse[15])	)
		#he3_onlySP_3.	append(	float(parse[16])	)
		#h3_onlySP_3.	append(	float(parse[17])	)

		#he3_full_4.	append(	float(parse[18])	)
		#h3_full_4.	append(	float(parse[19])	)
		#he3_onlySP_4.	append(	float(parse[20])	)
		#h3_onlySP_4.	append(	float(parse[21])	)

		#he3_full_5.	append(	float(parse[22])	)
		#h3_full_5.	append(	float(parse[23])	)
		#he3_onlySP_5.	append(	float(parse[24])	)
		#h3_onlySP_5.	append(	float(parse[25])	)

		#he3_full_6.	append(	float(parse[26])	)
		#h3_full_6.	append(	float(parse[27])	)
		#he3_onlySP_6.	append(	float(parse[28])	)
		#h3_onlySP_6.	append(	float(parse[29])	)

		#he3_full_7.	append(	float(parse[30])	)
		#h3_full_7.	append(	float(parse[31])	)
		#he3_onlySP_7.	append(	float(parse[32])	)
		#h3_onlySP_7.	append(	float(parse[33])	)

he3_full = np.asarray(he3_full);	h3_full = np.asarray(h3_full);
he3_onlySP = np.asarray(he3_onlySP);	h3_onlySP = np.asarray(h3_onlySP);
he3_full_1 = np.asarray(he3_full_1);		h3_full_1 = np.asarray(h3_full_1);
he3_onlySP_1 = np.asarray(he3_onlySP_1);	h3_onlySP_1 = np.asarray(h3_onlySP_1);
he3_full_2 = np.asarray(he3_full_2);		h3_full_2 = np.asarray(h3_full_2);
he3_onlySP_2 = np.asarray(he3_onlySP_2);	h3_onlySP_2 = np.asarray(h3_onlySP_2);
he3_full_3 = np.asarray(he3_full_3);		h3_full_3 = np.asarray(h3_full_3);
he3_onlySP_3 = np.asarray(he3_onlySP_3);	h3_onlySP_3 = np.asarray(h3_onlySP_3);
he3_full_4 = np.asarray(he3_full_4);		h3_full_4 = np.asarray(h3_full_4);
he3_onlySP_4 = np.asarray(he3_onlySP_4);	h3_onlySP_4 = np.asarray(h3_onlySP_4);
he3_full_5 = np.asarray(he3_full_5);		h3_full_5 = np.asarray(h3_full_5);
he3_onlySP_5 = np.asarray(he3_onlySP_5);	h3_onlySP_5 = np.asarray(h3_onlySP_5);
he3_full_6 = np.asarray(he3_full_6);		h3_full_6 = np.asarray(h3_full_6);
he3_onlySP_6 = np.asarray(he3_onlySP_6);	h3_onlySP_6 = np.asarray(h3_onlySP_6);
he3_full_7 = np.asarray(he3_full_7);		h3_full_7 = np.asarray(h3_full_7);
he3_onlySP_7 = np.asarray(he3_onlySP_7);	h3_onlySP_7 = np.asarray(h3_onlySP_7);



plt.figure(1)
plt.plot(x,he3_full,color='red',label='He-3 No Modification')
plt.plot(x,h3_full,color='blue',label='H-3 No Modification')
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
'''
plt.figure(3)
plt.plot(x,	he3_full - he3_onlySP ,	label='All')
plt.plot(x,	he3_full_1 - he3_onlySP_1 ,	label='>50')
plt.plot(x,	he3_full_2 - he3_onlySP_2 ,	label='>100')
plt.plot(x,	he3_full_3 - he3_onlySP_3 ,	label='>200')
plt.plot(x,	he3_full_4 - he3_onlySP_4 ,	label='>400')
plt.grid(True)
plt.xlabel('xB')
plt.ylabel('(Integral with Offshell + Motion) - (Integral with out offshell)')
plt.legend(numpoints=1,loc='best')
plt.savefig("off-breakdown.pdf",bbox_inches="tight")

plt.figure(4)
plt.plot(x,	(he3_full_1 - he3_onlySP_1)/(he3_full - he3_onlySP) ,	label='>50')
plt.plot(x,	(he3_full_2 - he3_onlySP_2)/(he3_full - he3_onlySP) ,	label='>100')
plt.plot(x,	(he3_full_3 - he3_onlySP_3)/(he3_full - he3_onlySP) ,	label='>200')
plt.plot(x,	(he3_full_4 - he3_onlySP_4)/(he3_full - he3_onlySP) ,	label='>400')
plt.plot(x,	(he3_full_5 - he3_onlySP_5)/(he3_full - he3_onlySP) ,	label='>600')
plt.plot(x,	(he3_full_6 - he3_onlySP_6)/(he3_full - he3_onlySP) ,	label='>800')
plt.plot(x,	(he3_full_7 - he3_onlySP_7)/(he3_full - he3_onlySP) ,	label='>1000')
plt.ylim([0,1])
plt.xlim([0,1])
plt.legend(numpoints=1,loc='best')
plt.grid(True)
plt.xlabel('xB')
plt.ylabel('Fraction of offshell due to certain momenta')
plt.savefig("off-breakdown2.pdf",bbox_inches="tight")


plt.show()

'''
