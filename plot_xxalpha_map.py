import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
	print "Incorrect number of arguments. Please use: "
	print "\tpython code.py [map output]"
	exit(-1)


x = []
xalpha = []
with open(sys.argv[1],"rb") as f:
	for line in f:
		parse = line.strip().split(" ")
		x.append(	float(parse[0])	)
		xalpha.append(	float(parse[1]) )



plt.plot(x,xalpha,linewidth=3)
xs=np.linspace(0,1,100)
plt.plot(xs,xs,linestyle='--',color='black')
plt.grid(True)
plt.ylim([0.1,1])
plt.xlim([0.1,1])
plt.xlabel('x')
plt.ylabel('Median x/alpha')
plt.savefig('xxalpha.pdf',bbox_inches='tight')
plt.show()
