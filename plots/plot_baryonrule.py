import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
	print "Incorrect number of arguments. Please use: "
	print "\tpython code.py [sum rule output]"
	exit(-1)



bar = []
x = []
mom = []
with open(sys.argv[1],"rb") as f:
	for line in f:
		parse = line.strip().split(" ")
		x.append(	float(parse[0])	)
		bar.append(	float(parse[2])*100	)
		mom.append(	float(parse[3])*100	)



plt.plot(x,bar,label='Baryon Sum')
plt.plot(x,mom,label='Mom Sum')
plt.legend(numpoints=1,loc='best')
plt.grid(True)
plt.axvline(x=0.825,color='black',linewidth=3)
plt.ylim([50,100])
plt.ylabel('Sum Rule (%)')
plt.xlabel('x')
plt.savefig('sumrule.pdf',bbox_inches='tight')
plt.show()
