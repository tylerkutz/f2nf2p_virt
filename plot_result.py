import numpy as np
import matplotlib.pyplot as plt


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

theo_x = []; theo_Q2 = []; theo_he3 = []; theo_h3 = [];
#with open("build/test-he3-h3-result.txt","rb") as f:
with open("build/test-no-h3-result.txt","rb") as f:
	for line in f:
		parse = line.strip().split(" ")
		theo_x.append(	float(parse[0])	)
		theo_Q2.append(	float(parse[1])	)
		theo_he3.append(	float(parse[2])	)
		theo_h3.append(		float(parse[3])	)


plt.errorbar(mar_x,mar_he3,yerr=mar_he3e,fmt='bo')
plt.plot(theo_x,theo_he3,color='blue',linestyle='--',linewidth=3)
plt.errorbar(mar_x,mar_h3,yerr=mar_h3e,fmt='ro')
plt.plot(theo_x,theo_h3,color='red',linestyle='--',linewidth=3)
plt.xlabel('x')
plt.ylabel('EMC Ratios')
plt.ylim([0.8,1.25])
plt.xlim([0.1,1])
plt.savefig("emc_he3.pdf",bbox_inches="tight")

plt.show()
